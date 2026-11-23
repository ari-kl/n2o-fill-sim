# Vessels: the things that hold fluid and therefore carry state (bottles and tanks).
#
#   stratified  : a liquid zone, an ullage zone and a wall, each at its own temperature,
#                 exchanging mass at a finite rate. This is what lets a vent pulse drop the
#                 pressure below the liquid's vapour pressure and then rebound.
#   homogeneous : one well-mixed fluid at (m, U) plus the wall. The fallback when there is no
#                 ullage left (liquid-full) or the contents are close to / above the critical
#                 temperature. state_ρu covers compressed liquid and supercritical fluid, so an
#                 overfilled tank warming up simply shows its pressure climbing.
#
# The stratified model follows Zilliac & Karabeyoglu (AIAA 2005-3549) as written out in Zimmerman,
# Waxman, Cantwell & Zilliac (AIAA 2013-4045): liquid at saturation density for its temperature,
# real-gas ullage, a massless saturated surface at the N2O partial pressure, evaporation set by the
# heat arriving at that surface, natural-convection correlations from that review's Table 1.
# Additions for filling: ports draw liquid or ullage according to the liquid level, inflow below
# the level mixes into the liquid while inflow above it is flashed and split between the zones,
# and the ullage can hold air.

const G0 = 9.80665
const R_AIR = 287.05
const SIGMA = 5.670374e-8
const LEVEL_BAND = 0.01       # m of level over which a port goes from drawing ullage to liquid
# A sliver of liquid mass that keeps the liquid-zone equations finite as the liquid vanishes. It
# scales with the vessel (0.9 g for a 19 L tank) so that it stays negligible in a small one too.
m_reg(n) = 1e-4 * N2O_RHOC * n.V
const ULLAGE_MIN = 0.005      # ullage volume fraction below which a vessel is treated as liquid-full
const T_STRAT_MAX = SAT_TMAX - 1.0   # above this the two-zone picture is dropped (near critical)

"Ambient conditions and model coefficients shared by every volume, in SI."
mutable struct Env
    T_amb::Float64
    solar::Float64
    wind::Float64
    P_atm::Float64
    E::Float64
    τ_cond::Float64
    h_in_mult::Float64
    h_ext_mult::Float64
end

Base.@kwdef mutable struct Node
    name::String
    kind::Symbol                  # :vessel, :junction or :atm
    i::Int = 0                    # index of its first state; 0 for :junction and :atm
    V::Float64 = 0.0              # m³
    D::Float64 = 0.0              # m inside diameter
    H::Float64 = 0.0              # m height
    Cw::Float64 = 0.0             # wall heat capacity, J/K
    A_in::Float64 = 0.0           # m² inside wall area
    h_in::Float64 = 1.0           # this vessel's multiplier on wall-to-fluid heat transfer
    h_mixed::Float64 = 300.0      # W/m²K wall-to-fluid coefficient when the contents are one mixed fluid
    split::Bool = false           # the wall as two masses (thin tube, heavy ends) with a temperature each
    tube_frac::Float64 = 0.5      # share of the wall's heat capacity that is in the tube
    link::Float64 = 5.0           # W/K conducted between tube and ends
    G_ends::Float64 = 0.0         # W/K from the bulk of the end metal to its wetted face (0: no such limit)
    h_fill::Float64 = 0.0         # multiplier on wall heat transfer driven by the fill jet (0: off)
    h_boil::Float64 = 0.0         # multiplier on boiling against a wall hotter than saturation (0: off)
    A_inlet::Float64 = 7.85e-5    # m² area of the inlet the fill jet comes through
    absorb::Float64 = 0.0
    emiss::Float64 = 0.0
    shade::Float64 = 0.0
    heater_W::Float64 = 0.0
    heater_on::Bool = false
    strat::Bool = false
    ports::Vector{Float64} = Float64[]        # height of each connection, as a fraction of H
    # evaluated from the state each time the equations are evaluated
    P::Float64 = 0.0              # for a junction: solved so that its flows balance
    T::Float64 = 0.0              # liquid (stratified) or bulk temperature
    hl_sat::Float64 = 0.0         # saturated liquid / vapour enthalpy at P, for splitting inflow
    hv_sat::Float64 = 0.0
    hl::Float64 = 0.0
    hv::Float64 = 0.0
    port_h::Vector{Float64} = Float64[]       # specific enthalpy a flow leaving through each port carries
    port_fl::Vector{Float64} = Float64[]      # liquid mass fraction of that flow
    port_cover::Vector{Float64} = Float64[]   # 0..1: how far the liquid level covers the port
    # flow totals collected from the edges
    dm_l::Float64 = 0.0           # net mass rate into the liquid zone (or the whole volume)
    min_l::Float64 = 0.0          # inflow-only part of dm_l
    Hin_l::Float64 = 0.0          # enthalpy rate carried by that inflow
    dm_v::Float64 = 0.0           # net mass rate into the ullage
    dE_v::Float64 = 0.0           # net enthalpy rate into the ullage (or the whole volume)
    dair::Float64 = 0.0           # net air mass rate
    P_hold::Float64 = 0.0         # junction: its pressure at the last accepted step
end

Node(name, kind) = Node(; name, kind)

nstates(n::Node) = n.kind === :vessel ? (n.split ? 7 : 6) : 0

"Wall temperature: the one state, or with a split wall the mean of tube and ends weighted by heat capacity."
wall_T(n::Node, u) = n.split ? n.tube_frac * u[n.i+4] + (1 - n.tube_frac) * u[n.i+6] : u[n.i+4]

# Stratified state layout: [m_liquid, T_liquid, m_vapour, T_vapour, T_wall, m_air]
# Homogeneous vessel:      [m, U, -, -, T_wall, m_air]
# A junction has no state: `P` and `hl` (its mixture enthalpy) are found by `solve_junctions!`.

"Everything about a stratified vessel that follows from its state."
function strat(n::Node, u)
    i = n.i
    ml, mv, mair = max(u[i], 0.0), max(u[i+2], 1e-9), max(u[i+5], 0.0)
    Tl = clamp(u[i+1], N2O_TTRIP, SAT_TMAX)
    Tv = max(u[i+3], 150.0)
    z = clamp(z_of_T(Tl), _ZMIN, _ZMAX)
    dTdz = _spline(1, z)[2]
    ρl, dρl = _spline(3, z)
    ul, dul = _spline(5, z)
    dρl /= dTdz                              # along the saturation line, per kelvin
    cl = dul / dTdz
    Vl = min(ml / ρl, (1 - ULLAGE_MIN / 2) * n.V)
    Vv = n.V - Vl
    ρv = mv / Vv
    # The ullage may be slightly supersaturated before it condenses. Past the vapour spinodal the
    # EOS is meaningless, so evaluate no denser than 5 % above saturation and extend linearly.
    ρe = Tv < SAT_TMAX ? min(ρv, 1.05 * _spline(4, z_of_T(Tv))[1]) : ρv
    ev = eos(ρe, Tv)
    Pn = ev.P + ev.dPdρ * (ρv - ρe)          # N2O partial pressure
    P = Pn + mair * R_AIR * Tv / Vv
    level = Vl / (π * n.D^2 / 4)
    (; ml, Tl, mv, Tv, mair, ρl, dρl, ul, cl, Vl, Vv, ρv, ρe, ev, Pn, P, level, hl = ul + P / ρl, hv = ev.u + P / ρv)
end

"Fill in a vessel's (or the atmosphere's) pressure and what each port would deliver."
function evaluate!(n::Node, u, env::Env)
    n.dm_l = n.min_l = n.Hin_l = n.dm_v = n.dE_v = n.dair = 0.0
    n.kind === :junction && return
    if n.kind === :atm
        n.P, n.T = env.P_atm, env.T_amb
        n.hl = n.hv = vapour_PT(env.P_atm, env.T_amb).h
        fill!(n.port_h, n.hl)
        return
    end
    if n.strat
        s = strat(n, u)
        n.P, n.T, n.hl, n.hv = s.P, s.Tl, s.hl, s.hv
        for k in eachindex(n.ports)
            φ = clamp((s.level - n.ports[k] * n.H) / LEVEL_BAND, 0.0, 1.0)
            φ = φ^2 * (3 - 2φ)
            fl = φ * s.ρl / (φ * s.ρl + (1 - φ) * s.ρv)
            n.port_cover[k] = φ
            n.port_fl[k] = fl
            n.port_h[k] = fl * s.hl + (1 - fl) * s.hv
        end
    else
        f = homog(n, u)
        n.P, n.T, n.hl, n.hv = f.P + air_pressure(n, u, f), f.T, f.h, f.h
        fill!(n.port_h, f.h)
        fill!(n.port_fl, 1 - f.x)
    end
    st = sat_P(clamp(n.P, SAT_PMIN, SAT_PMAX))
    n.hl_sat, n.hv_sat = st.hl, st.hv
    nothing
end

homog(n::Node, u) = (m = max(u[n.i], 1e-12); state_ρu(m / n.V, u[n.i+1] / m))

# Air in a homogeneous vessel only matters while there is a vapour phase for it to sit in.
function air_pressure(n::Node, u, f::Fluid)
    mair = max(u[n.i+5], 0.0)
    mair == 0 && return 0.0
    Vv = f.twophase ? max(u[n.i], 0.0) * f.x / sat_T(f.T).ρv : f.x >= 1 ? n.V : 0.0
    Vv > ULLAGE_MIN * n.V ? mair * R_AIR * f.T / Vv : 0.0
end

"Record a flow `mdot` ≥ 0 leaving vessel `n` through its port `k`."
function outflow!(n::Node, k, mdot, u)
    n.kind === :vessel || return
    if n.strat
        fl = n.port_fl[k]
        n.dm_l -= fl * mdot
        n.dm_v -= (1 - fl) * mdot
        n.dE_v -= (1 - fl) * mdot * n.hv
        # Air is a massless passenger. It adds partial pressure and leaves with the
        # ullage gas, but is not counted in flow rates or energy. Fine for a tank's worth of air
        # (tens of grams); model it as a real species if N2 pressurisation is ever added.
        n.dair -= (1 - fl) * mdot * max(u[n.i+5], 0.0) / max(u[n.i+2], 1e-9)
    else
        n.dm_l -= mdot
        n.dE_v -= mdot * n.hl
        n.dair -= mdot * max(u[n.i+5], 0.0) / max(u[n.i], 1e-9)
    end
    nothing
end

"""
Record a flow `mdot` ≥ 0 with specific enthalpy `h` entering vessel `n` through its port `k`.

Through a port below the liquid level the stream joins the liquid zone with all its enthalpy: any
flash vapour bubbles up through the liquid and gives its heat to it, and the liquid then boils at
the surface if that leaves it superheated. Through a port in the ullage the stream is flashed at
tank pressure and split, vapour to the ullage and liquid falling into the liquid zone.
"""
function inflow!(n::Node, k, mdot, h)
    n.kind === :vessel || return
    if n.strat
        sub = n.port_cover[k] * mdot
        n.dm_l += sub
        n.min_l += sub
        n.Hin_l += sub * h
        mdot -= sub
        x = clamp((h - n.hl_sat) / (n.hv_sat - n.hl_sat), 0.0, 1.0)
        h_l = x > 0 ? n.hl_sat : h
        h_v = x < 1 ? n.hv_sat : h
        n.dm_l += (1 - x) * mdot
        n.min_l += (1 - x) * mdot
        n.Hin_l += (1 - x) * mdot * h_l
        n.dm_v += x * mdot
        n.dE_v += x * mdot * h_v
    else
        n.dm_l += mdot
        n.dE_v += mdot * h
    end
    nothing
end

# ---------------------------------------------------------------- heat transfer

# Natural convection, Nu = c·Ra^n (Zimmerman et al. 2013, Table 1, after Bergman et al.).
function _ra_group(ρ, cp, β, μ, k, ΔT)
    G0 * abs(β) * (abs(ΔT) + 1e-3) * ρ^2 * cp / (μ * k)      # Ra / L³
end
h_surface(ρ, cp, β, μ, k, ΔT) = 0.15 * k * cbrt(_ra_group(ρ, cp, β, μ, k, ΔT))
h_wall(ρ, cp, β, μ, k, ΔT, L) = 0.021 * k / L * (_ra_group(ρ, cp, β, μ, k, ΔT) * L^3)^0.4

"Film coefficient on the outside of a wall at `Tw`: natural convection and wind, W/m²K."
function h_outside(n::Node, Tw, env::Env)
    k, ν, α, Pr = 0.026, 1.6e-5, 2.2e-5, 0.71
    ΔT = env.T_amb - Tw
    L = n.H
    Ra = G0 / env.T_amb * (abs(ΔT) + 1e-3) * L^3 / (ν * α)
    h_nat = 0.59 * Ra^0.25 * k / L
    h_wind = 0.193 * (env.wind * n.D / ν + 1e-9)^0.618 * cbrt(Pr) * k / n.D
    env.h_ext_mult * cbrt(h_nat^3 + h_wind^3)
end
"""
Heat flow into a wall from outside [W]: convection, sunlight, infrared exchange and any heater.
Air properties are fixed at round room-temperature values. Convection is natural (0.59·Ra^¼,
vertical surface) combined with wind across a cylinder (Hilpert, 0.193·Re^0.618·Pr^⅓).
"""
function q_external(n::Node, Tw, env::Env)
    ΔT = env.T_amb - Tw
    A = π * n.D * n.H + π * n.D^2 / 2
    h_outside(n, Tw, env) * A * ΔT +
        n.absorb * (1 - n.shade) * env.solar * n.D * n.H +
        n.emiss * SIGMA * A * (env.T_amb^4 - Tw^4) +
        (n.heater_on ? n.heater_W : 0.0)
end
# The same in two pieces, for a split wall: what an outside area A at its own temperature takes by
# convection and radiation, and what lands on the tube alone (sun on its side, the heater).
q_outside(n::Node, Tw, env::Env, A) = h_outside(n, Tw, env) * A * (env.T_amb - Tw) + n.emiss * SIGMA * A * (env.T_amb^4 - Tw^4)
q_sun_heater(n::Node, env::Env) = n.absorb * (1 - n.shade) * env.solar * n.D * n.H + (n.heater_on ? n.heater_W : 0.0)

# A surface conductance `ua` in series with a conduction path `G` through the metal behind it
# (G = 0: no such path is modelled). A thick bulkhead gives up its heat only as fast as the metal
# conducts it to the wetted face, however hard the liquid boils there.
_through(ua, G) = G > 0 ? ua * G / (ua + G) : ua

"""
    boiling_flux(P, ΔT, st) -> W/m²

Heat flux from a wall `ΔT` kelvin hotter than saturation into liquid boiling on it at pressure
`P`. This is the pool-boiling data of Gorenflo and Kotthoff (2005) for carbon dioxide in the form
Zimmerman (2015, eqs 5.76-5.77) fitted it for tank models, taken over to nitrous oxide by reduced
pressure Pr = P/Pc (the two fluids are close analogues): the coefficient at 20 kW/m² is

    ln h20 [kW/m²K] = 0.3092 (ln Pr)³ + 1.649 (ln Pr)² + 3.641 ln Pr + 5.272

and the coefficient grows with flux as h ∝ q^n, n = 1 - 0.3 Pr^0.3, so that
q = (h20 ΔT)^(1/(1-n)) / q20^(n/(1-n)). It rises about as ΔT⁴: a wetted wall cannot stay more than
a kelvin or two above saturation. The data stop at Pr = 0.625, so warm tanks are an extrapolation.
The flux is capped at the critical heat flux (Zuber: 0.131 h_fg √ρv [σ g (ρl - ρv)]^¼), beyond
which a real surface dries out; that drying-out is not modelled. `st` is the saturation state at P.
"""
function boiling_flux(P, ΔT, st)
    ΔT > 0 || return 0.0
    Pr = clamp(P / N2O_PC, 0.05, 0.99)
    L = log(Pr)
    h20 = exp(((0.3092 * L + 1.649) * L + 3.641) * L + 5.272)       # kW/m²K at q = 20 kW/m²
    n = 1 - 0.3 * Pr^0.3
    q = 1e3 * (h20 * ΔT)^(1 / (1 - n)) / 20.0^(n / (1 - n))
    chf = 0.131 * (st.hv - st.hl) * sqrt(st.ρv) * sqrt(sqrt(surface_tension(st.T) * G0 * (st.ρl - st.ρv)))
    min(q, chf)
end

"""
Wall-to-liquid film coefficient for convection. Natural convection always. With `h_fill` on, the
stirring of the liquid by the incoming jet is added (Dittus-Boelter on the tank diameter, with the
jet's velocity decayed over one diameter as 6·d/D, combined as the cube root of the sum of cubes
so that it takes over smoothly). The multiplier is there to be fitted. Boiling is separate
(`q_boiling`): it is a heat flow set by the wall's superheat, not a coefficient on Tw - Tl.
"""
function h_wall_liquid(n::Node, s, cpl, βl, μl, kl, Tw, ss, lvl)
    h = h_wall(s.ρl, cpl, βl, μl, kl, Tw - s.Tl, max(lvl, 0.01))
    if n.h_fill > 0
        d = sqrt(4 * n.A_inlet / π)
        v = n.min_l / (s.ρl * n.A_inlet) * min(6 * d / n.D, 1.0)
        Re, Pr = s.ρl * v * n.D / μl, cpl * μl / kl
        h = cbrt(h^3 + (n.h_fill * 0.023 * Re^0.8 * Pr^0.4 * kl / n.D)^3)
    end
    h
end

# Two heat flows in series (both in the same direction, W): what gets through is less than either.
_series(a, b) = a + b > 0 ? a * b / (a + b) : 0.0

"""
Heat flow [W] by boiling from a wall at `Tw` into the liquid over wetted area `A`: zero at or
below saturation and rising smoothly from there (`boiling_flux`), scaled by the vessel's
`h_boil`. It depends on how far the wall is above saturation, not on the bulk liquid's
temperature. With `G` > 0 it is also limited by conduction through the metal behind the surface.
"""
function q_boiling(n::Node, ss, Tw, A, G = 0.0)
    q = n.h_boil * boiling_flux(ss.P, Tw - ss.T, ss) * A
    G > 0 ? _series(q, G * max(Tw - ss.T, 0.0)) : q
end
function derivative!(du, n::Node, u, env::Env)
    n.kind === :vessel || return 0.0
    i = n.i
    if !n.strat
        f = homog(n, u)
        Tw = u[i+4]
        # One fixed film coefficient while liquid-full / supercritical. Good enough to
        # get the warm-soak pressure rise; use the correlations if this regime ever matters more.
        if n.split
            Te, At, Ae = u[i+6], π * n.D * n.H, π * n.D^2 / 2
            hm = env.h_in_mult * n.h_in * n.h_mixed
            qt, qe, ql = hm * At * (Tw - f.T), _through(hm * Ae, n.G_ends) * (Te - f.T), n.link * (Te - Tw)
            qxt, qxe = q_outside(n, Tw, env, At) + q_sun_heater(n, env), q_outside(n, Te, env, Ae)
            du[i], du[i+1], du[i+2], du[i+3] = n.dm_l, n.dE_v + qt + qe, 0.0, 0.0
            du[i+4], du[i+5] = (qxt - qt + ql) / (n.Cw * n.tube_frac), n.dair
            du[i+6] = (qxe - qe - ql) / (n.Cw * (1 - n.tube_frac))
            return qxt + qxe
        end
        qw = env.h_in_mult * n.h_in * n.h_mixed * n.A_in * (Tw - f.T)
        qx = q_external(n, Tw, env)
        du[i], du[i+1], du[i+2], du[i+3] = n.dm_l, n.dE_v + qw, 0.0, 0.0
        du[i+4], du[i+5] = (qx - qw) / n.Cw, n.dair
        return qx
    end
    s = strat(n, u)
    Tw = u[i+4]
    wet = s.ml / (s.ml + m_reg(n))                 # fades liquid-side exchange as the liquid vanishes
    As = π * n.D^2 / 4
    # Surface: saturated at the N2O partial pressure.
    ss = sat_P(clamp(s.Pn, SAT_PMIN, SAT_PMAX))
    lp = sat_T(s.Tl)
    cpv = cp_of(s.ev, s.ρe, s.Tv)
    βv = s.ev.dPdT / (s.ρe * s.ev.dPdρ)
    μl, kl, μv, kv = visc_l(s.Tl), cond_l(s.Tl), visc_v(s.Tv), cond_v(s.Tv)
    βl = -s.dρl / s.ρl
    boil = s.Tl > ss.T ? env.E : 1.0           # heat reaches the surface far faster when the liquid boils
    q_ls = wet * boil * h_surface(s.ρl, lp.cpl, βl, μl, kl, s.Tl - ss.T) * As * (s.Tl - ss.T)
    q_sv = wet * h_surface(s.ρv, cpv, βv, μv, kv, ss.T - s.Tv) * As * (ss.T - s.Tv)
    m_evap = (q_ls - q_sv) / max(ss.hv - s.hl, 2e4)
    # Supersaturated ullage condenses in the bulk and rains into the liquid.
    m_cond, h_cond = 0.0, 0.0
    if s.Tv < SAT_TMAX
        sv = sat_T(s.Tv)
        m_cond = max(s.ρv - sv.ρv, 0.0) * s.Vv / env.τ_cond
        h_cond = sv.hl
    end
    # Wall.
    lvl = min(s.level, n.H)
    extra = n.h_fill > 0                           # the fill-jet term is off unless asked for
    dTw_tube = dTw_ends = 0.0
    if n.split
        # Tube: its wetted part to the liquid, its dry part to the ullage, sun and heater on its
        # outside. Ends: the bottom one to the liquid, the top one to the ullage. A link between.
        Te, c = u[i+6], env.h_in_mult * n.h_in
        hl_t = extra ? h_wall_liquid(n, s, lp.cpl, βl, μl, kl, Tw, ss, lvl) : h_wall(s.ρl, lp.cpl, βl, μl, kl, Tw - s.Tl, max(lvl, 0.01))
        hl_e = extra ? h_wall_liquid(n, s, lp.cpl, βl, μl, kl, Te, ss, lvl) : h_wall(s.ρl, lp.cpl, βl, μl, kl, Te - s.Tl, max(lvl, 0.01))
        hv_t = h_wall(s.ρv, cpv, βv, μv, kv, Tw - s.Tv, max(n.H - lvl, 0.01))
        hv_e = h_wall(s.ρv, cpv, βv, μv, kv, Te - s.Tv, max(n.H - lvl, 0.01))
        ql_t, ql_e = wet * c * hl_t * (π * n.D * lvl) * (Tw - s.Tl), _through(wet * c * hl_e * As, n.G_ends) * (Te - s.Tl)
        if n.h_boil > 0                            # boiling on the wetted tube and on the bottom end
            ql_t += wet * q_boiling(n, ss, Tw, π * n.D * lvl)
            ql_e += wet * q_boiling(n, ss, Te, As, n.G_ends)
        end
        qv_t, qv_e = c * hv_t * (π * n.D * (n.H - lvl)) * (Tw - s.Tv), _through(c * hv_e * As, n.G_ends) * (Te - s.Tv)
        link = n.link * (Te - Tw)
        qxt, qxe = q_outside(n, Tw, env, π * n.D * n.H) + q_sun_heater(n, env), q_outside(n, Te, env, 2As)
        q_wl, q_wv, qx = ql_t + ql_e, qv_t + qv_e, qxt + qxe
        dTw_tube = (qxt - ql_t - qv_t + link) / (n.Cw * n.tube_frac)
        dTw_ends = (qxe - ql_e - qv_e - link) / (n.Cw * (1 - n.tube_frac))
    else
        h_l = extra ? h_wall_liquid(n, s, lp.cpl, βl, μl, kl, Tw, ss, lvl) : h_wall(s.ρl, lp.cpl, βl, μl, kl, Tw - s.Tl, max(lvl, 0.01))
        q_wl = wet * env.h_in_mult * n.h_in * h_l * (π * n.D * lvl + As) * (Tw - s.Tl)
        n.h_boil > 0 && (q_wl += wet * q_boiling(n, ss, Tw, π * n.D * lvl + As))
        q_wv = env.h_in_mult * n.h_in * h_wall(s.ρv, cpv, βv, μv, kv, Tw - s.Tv, max(n.H - lvl, 0.01)) * (π * n.D * (n.H - lvl) + As) * (Tw - s.Tv)
        qx = q_external(n, Tw, env)
    end

    dml = n.dm_l - m_evap + m_cond
    Hl = n.Hin_l - n.min_l * s.hl - q_ls + q_wl + m_cond * (h_cond - s.hl)
    dTl = Hl / ((s.ml + m_reg(n)) * (s.cl - s.P * s.dρl / s.ρl^2))
    dVv = -(dml / s.ρl - s.ml * s.dρl * dTl / s.ρl^2)
    dmv = n.dm_v + m_evap - m_cond
    dUv = n.dE_v + m_evap * ss.hv - m_cond * h_cond + q_sv + q_wv - s.P * dVv
    dρv = (dmv - s.ρv * dVv) / s.Vv
    du_dρ = (s.ev.P - s.Tv * s.ev.dPdT) / s.ρe^2
    du[i] = dml
    du[i+1] = dTl
    du[i+2] = dmv
    du[i+3] = (dUv - s.ev.u * dmv - s.mv * du_dρ * dρv) / (s.mv * s.ev.cv)
    du[i+4] = n.split ? dTw_tube : (qx - q_wl - q_wv) / n.Cw
    du[i+5] = n.dair
    n.split && (du[i+6] = dTw_ends)
    qx
end

"N2O mass and total energy (fluid + wall) held by a node."
function inventory(n::Node, u)
    n.kind === :vessel || return 0.0, 0.0
    i = n.i
    if n.strat
        s = strat(n, u)
        return s.ml + s.mv, s.ml * s.ul + s.mv * s.ev.u + n.Cw * wall_T(n, u)
    end
    u[i], u[i+1] + n.Cw * wall_T(n, u)
end

"""
Switch a vessel between the stratified and homogeneous descriptions when it has to, rewriting its
states so mass and energy carry over. Returns true if it switched.
"""
function switch_mode!(n::Node, u)
    n.kind === :vessel || return false
    i = n.i
    if n.strat
        s = strat(n, u)
        (s.Vv < ULLAGE_MIN * n.V || s.Tl > T_STRAT_MAX || s.Tv > 2 * N2O_TC) || return false
        u[i], u[i+1] = s.ml + s.mv, s.ml * s.ul + s.mv * s.ev.u
        u[i+2] = u[i+3] = 0.0
        n.strat = false
        return true
    end
    f = homog(n, u)
    f.T < T_STRAT_MAX - 1.0 || return false
    m = u[i]
    if f.twophase
        st = sat_T(f.T)
        m * f.x / st.ρv > 2 * ULLAGE_MIN * n.V || return false
        u[i], u[i+1], u[i+2], u[i+3] = m * (1 - f.x), f.T, m * f.x, f.T
    elseif f.x >= 1
        u[i], u[i+1], u[i+2], u[i+3] = 0.0, f.T, m, f.T
    else
        return false                           # liquid-full: stay homogeneous
    end
    n.strat = true
    true
end
