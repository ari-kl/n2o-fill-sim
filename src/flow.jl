# Mass flux through a restriction, per unit effective area (Cd·A), for any upstream state.
#
#   liquid or two-phase upstream : Dyer et al. non-homogeneous non-equilibrium blend of the
#       incompressible (SPI) flux and the homogeneous-equilibrium (HEM) flux, in the form given as
#       eq. 7 of Zimmerman, Waxman, Cantwell & Zilliac, AIAA 2013-4045:
#           G = (κ·G_SPI + G_HEM) / (1 + κ),   κ = sqrt((P1 - P2) / (Psat(T1) - P2))
#       Subcooled liquid discharging above its vapour pressure cannot flash, so it is pure SPI.
#   vapour or supercritical upstream : HEM alone, which is real-gas isentropic nozzle flow
#       (including condensation in the jet) and reduces to the ideal-gas orifice formulas at low
#       pressure.
#
# HEM flux at a throat pressure P is ρ(P, s1)·sqrt(2·(h1 - h(P, s1))). It has a maximum over P; if
# the downstream pressure is below the pressure at that maximum, the flow is choked there.

const DP_REG = 300.0          # Pa. Scale of the smoothing of √ΔP around zero (see `massflux`).

# ρ·√(2Δh) after an isentropic expansion from `up` to pressure P. This is the innermost loop of
# the whole simulator, so the two-phase case reads only the six saturation columns it needs
# instead of building a full state.
# The same flux at a point inside the dome given by its place on the saturation table, or NaN if
# the isentrope is not inside the dome there. Six columns are read at one place.
function _g_dome(up::Fluid, z, inside = false)
    i, t = _locate(z)
    sl, sv = _val(9, i, t), _val(10, i, t)
    inside || sl <= up.s <= sv || return NaN
    x = clamp((up.s - sl) / (sv - sl), 0.0, 1.0)
    hl, hv = _val(7, i, t), _val(8, i, t)
    ρ = 1 / (x / _val(4, i, t) + (1 - x) / _val(3, i, t))
    ρ * sqrt(2 * max(up.h - (hl + x * (hv - hl)), 0.0))
end

# Where on the saturation table the isentrope through a single-phase state enters the dome as the
# pressure falls, or NaN if it never does above the triple-point pressure.
function _dome_entry(up::Fluid)
    s = up.s
    SAT[1, 9] <= s <= SAT[1, 10] && return _ZMIN                  # through the top, beside the critical point
    s < SAT[1, 9] && return s >= SAT[_NZ, 9] ? _invert(9, s) : NaN   # from the liquid side
    s <= SAT[_NZ, 10] ? _invert(10, s) : NaN                      # from the vapour side
end

# Has the expansion from `up` already passed the speed of sound when it reaches the dome at `ze`?
# The state there is saturated liquid or vapour, read from the table.
function _supersonic_at_entry(up::Fluid, ze)
    i, t = _locate(ze)
    liquid = up.s < SAT[1, 9]
    ρ, T, h = _val(liquid ? 3 : 4, i, t), _val(1, i, t), _val(liquid ? 7 : 8, i, t)
    e = eos(ρ, T)
    2 * (up.h - h) > e.dPdρ + T * e.dPdT^2 / (ρ^2 * e.cv)         # velocity² against (∂P/∂ρ) at constant entropy
end

function _g_hem(up::Fluid, P)
    if P < SAT_PMAX
        g = _g_dome(up, _invert(2, max(P, SAT_PMIN)))
        isnan(g) || return g
    end
    # Single phase after the expansion. Start Newton from where an ideal gas (γ ≈ 1.3) or an
    # incompressible liquid would end up, which is close enough to converge in two or three steps.
    f = if up.twophase
        state_Ps(P, up.s)
    elseif up.x >= 1
        r = P / up.P
        _solve_Py(P, up.s, up.ρ * r^0.77, up.T * r^0.23, true)
    else
        _solve_Py(P, up.s, up.ρ, up.T, true)
    end
    f.ρ * sqrt(2 * max(up.h - f.h, 0.0))
end

# Where a function with a single maximum on [a, b] peaks, and its value there: Brent's method
# (parabolic steps where the function is smooth, golden-section steps where it is not, such as at
# the kink where a subcooled liquid starts to flash). About a dozen evaluations.
function _argmax(f, a, b, tol)
    gold = 0.3819660112501051
    x = w = v = a + gold * (b - a)
    fx = fw = fv = -f(x)
    d = e = 0.0
    for _ in 1:60
        m = (a + b) / 2
        abs(x - m) <= 2tol - (b - a) / 2 && break
        parabolic = false
        if abs(e) > tol
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2 * (q - r)
            q > 0 && (p = -p)
            q = abs(q)
            last_e, e = e, d
            if abs(p) < abs(q * last_e / 2) && p > q * (a - x) && p < q * (b - x)
                d = p / q
                u = x + d
                (u - a < 2tol || b - u < 2tol) && (d = m > x ? tol : -tol)
                parabolic = true
            end
        end
        if !parabolic
            e = x >= m ? a - x : b - x
            d = gold * e
        end
        u = abs(d) >= tol ? x + d : x + (d > 0 ? tol : -tol)
        fu = -f(u)
        if fu <= fx
            u >= x ? (a = x) : (b = x)
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else
            u < x ? (a = u) : (b = u)
            if fu <= fw || w == x
                v, w = w, u
                fv, fw = fw, fu
            elseif fu <= fv || v == x || v == w
                v, fv = u, fu
            end
        end
    end
    x, -fx
end

"""
HEM flux from `up` to back pressure P2, choked if P2 is below the critical pressure.

The critical pressure and flux depend only on the upstream state, and the same upstream state is
asked about many times in a row (every sweep of a junction solve, every back-pressure
perturbation). `memo`, if given, is any mutable object with fields `m_Pc`, `m_G` and `m_Pu`; the
caller resets `m_Pc` to NaN and `m_Pu` to Inf when the upstream state changes. `m_Pc`/`m_G` keep the
choking point once it has been searched for; `m_Pu` is the lowest back pressure already shown
not to choke (the flux has a single maximum, so nothing above it chokes either).
"""
function hem_flux(up::Fluid, P2, memo = nothing)
    if memo !== nothing
        isnan(memo.m_Pc) || return P2 <= memo.m_Pc ? memo.m_G : _g_hem(up, P2)
        P2 >= memo.m_Pu && return _g_hem(up, P2)
    end
    g2 = _g_hem(up, P2)
    span = up.P - P2
    # Not choked if the flux is still rising as the throat pressure falls towards P2. The test
    # point sits close to P2 on purpose: with P2 just below the choking pressure this test still
    # says "rising", and the flux returned is short of the maximum by about (gap / pressure)².
    # A gap of 1e-3 of the span left steps of 1e-6 in the flux, which a junction solve working to
    # 1e-9 cannot get past.
    # Coming from single phase the flux can also peak where the expansion enters the dome (see
    # below), and a flux still rising at a P2 just beyond that peak is not the maximum.
    ze = up.twophase ? NaN : _dome_entry(up)
    Pe = isnan(ze) ? -Inf : _spline(2, ze)[1]
    ge = P2 < Pe < up.P ? _g_dome(up, ze, true) : 0.0
    if g2 >= ge && g2 > _g_hem(up, P2 + 3e-5 * span)
        memo === nothing || (memo.m_Pu = P2)
        return g2
    end
    TALLY.choke += 1
    if up.twophase && up.P < SAT_PMAX
        # From inside the dome the whole expansion stays inside it, so search along the saturation
        # table itself: no pressure has to be turned back into a place on the table at each trial.
        z1, z2 = _invert(2, up.P), _invert(2, P2)
        zc, Gc = _argmax(z -> _g_dome(up, z), z1, z2, 1e-7)       # 1e-7 in z is about 1e-6 of the pressure
        if !isnan(Gc)
            Pc = _spline(2, zc)[1]
            memo === nothing || (memo.m_Pc = Pc; memo.m_G = Gc)
            return Gc
        end
    end
    if !up.twophase
        # From single phase the expansion enters the dome at some pressure Pe (unless it is hot
        # gas that never condenses). Below Pe: search along the table as above. At Pe the speed of
        # sound drops sharply, so the flux often peaks exactly there: that point is evaluated,
        # not searched for. Above Pe the flux rises as the pressure falls for as long as the flow
        # is subsonic, and in single phase it normally still is at Pe: then nothing above Pe can
        # beat the value at Pe, and no equation-of-state search is needed at all.
        if P2 < Pe < up.P
            zc, Gc = _argmax(z -> _g_dome(up, z, true), ze, _invert(2, P2), 1e-7)
            Pc = _spline(2, zc)[1]
            ge >= Gc && ((Pc, Gc) = (Pe, ge))
            if ze == _ZMIN || _supersonic_at_entry(up, ze)
                Ps, Gs = _argmax(P -> _g_hem(up, P), Pe, up.P, 1e-6 * up.P)
                Gs > Gc && ((Pc, Gc) = (Ps, Gs))
            end
            memo === nothing || (memo.m_Pc = Pc; memo.m_G = Gc)
            return Gc
        end
    end
    Pc, Gc = _argmax(P -> _g_hem(up, P), P2, up.P, 1e-6 * up.P)   # the maximum is flat: this gives the flux to ~1e-11
    memo === nothing || (memo.m_Pc = Pc; memo.m_G = Gc)
    Gc
end

"""
    massflux(P1, h1, P2) -> G [kg/m²/s], ≥ 0

Flux per unit effective area from a reservoir at (P1, h1) to back pressure P2 < P1. Zero if
P2 ≥ P1. Close to ΔP = 0 the true flux goes like √ΔP, whose slope is infinite; that defeats
Newton's method around a stalled fill, so within about `DP_REG` of zero the flux is smoothed to
pass through zero with a finite slope.
"""
massflux(P1, h1, P2) = P1 > P2 ? massflux(state_Ph(P1, h1), P1, P2) : 0.0

const NEAR_CRIT = 2.0         # K below the top of the dome over which the non-equilibrium share fades out

"The same, given the upstream state already resolved (`up = state_Ph(P1, h1)`)."
function massflux(up::Fluid, P1, P2, memo = nothing)
    dP = P1 - P2
    dP <= 0 && return 0.0
    # Evaluate the model no closer than DP_REG, then reshape √ΔP into ΔP/(ΔP² + δ²)^¼, which is
    # smooth through zero and indistinguishable from √ΔP a few δ away (Modelica's regRoot).
    scale = dP / (sqrt(sqrt(dP^2 + DP_REG^2)) * sqrt(max(dP, DP_REG)))
    P2 = max(min(P2, P1 - DP_REG), SAT_PMIN)
    P2 < P1 || return 0.0                  # upstream below the triple-point pressure: nothing to model
    # Approaching the critical temperature liquid and vapour become one fluid, so the
    # non-equilibrium (metastable liquid) share fades out over the last NEAR_CRIT kelvin and what
    # is left is the equilibrium model, which is also what is used above it. Without the fade the
    # flux jumps by about 15 % where a supercritical supply crosses into the two-phase region, and
    # a junction solve straddling that jump cannot settle.
    fade = clamp((SAT_TMAX - up.T) / NEAR_CRIT, 0.0, 1.0)
    (up.x >= 1 || fade <= 0) && return scale * hem_flux(up, P2, memo)
    gspi = sqrt(2 * up.ρ * (P1 - P2))
    Pv = Psat(up.T)
    # Subcooled liquid that stays above its vapour pressure cannot flash: pure incompressible flow,
    # and no need to evaluate the equilibrium model at all.
    fade >= 1 && up.x <= 0 && P2 >= Pv && return scale * gspi
    ghem = hem_flux(up, P2, memo)
    w = fade * (P2 < Pv ? (κ = sqrt((P1 - P2) / (Pv - P2)); κ / (1 + κ)) : 1.0)
    gdyer = w * gspi + (1 - w) * ghem
    scale * ((1 - up.x) * gdyer + up.x * ghem)   # two-phase upstream: fade to HEM with quality
end

"""
Darcy friction factor, Churchill (1977), turbulent branch. Below Re = 4000 the value at 4000 is
used: laminar flow in these lines means grams per second, and the true laminar law (flow ∝ ΔP with
a very steep slope, giving way to √ΔP within a pascal) is a shape Newton's method cannot follow.
"""
function friction_factor(Re, relrough)
    Re = max(Re, 4000.0)
    A = _p16(-2.457 * log((7 / Re)^0.9 + 0.27 * relrough))
    B = _p16(37530 / Re)
    c = (8 / Re)^2
    c *= c                                             # (8/Re)^4
    8 * cbrt(sqrt(sqrt(c * c * c + 1 / ((A + B) * sqrt(A + B)))))
end
_p16(x) = (a = x * x; b = a * a; c = b * b; c * c)     # x^16 by squaring: this function is hot

"Viscosity of a fluid for a pipe Reynolds number (McAdams mixing rule when two-phase)."
visc_mix(up::Fluid) = 1 / (up.x / visc_v(up.T) + (1 - up.x) / visc_l(up.T))

const CV_AREA = 1.698e-5      # m² of Cd·A per unit of US valve flow coefficient Cv
