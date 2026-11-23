# Nitrous oxide properties. Everything is SI, per kg, absolute pressure.
#
# Equation of state: Lemmon & Span (2006), "Short fundamental equations of state for 20
# industrial fluids", J. Chem. Eng. Data 51, 785 (doi:10.1021/je050186n). Coefficients are copied
# from CoolProp's dev/fluids/NitrousOxide.json (the same EOS NIST REFPROP uses). test/runtests.jl
# checks this file against CoolProp itself.
#
# Why a native EOS and not CoolProp at runtime: the tank model needs states off the saturation line
# (superheated ullage, compressed liquid, liquid-full, supercritical soak), it has to be fast and
# never throw inside the ODE right-hand side, and the tool must run with no binaries to fetch.

const N2O_TC = 309.52                      # K
const N2O_RHOC = 10270 * 0.0440128         # kg/m³
const N2O_PC = 7.245e6                     # Pa
const N2O_R = 8.314472 / 0.0440128         # J/kg/K
const N2O_TTRIP = 182.33                   # K

const _N = (0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, -0.0036985,
            -0.23263, -0.00042859, -0.04281, -0.023038)
const _D = (1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2)
const _T = (0.25, 1.25, 1.5, 0.25, 0.875, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5)
const _L = (0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3)
const _A1, _A2 = -4.4262736272, 4.3120475243
const _PE_N = (2.1769, 1.6145, 0.48393)
const _PE_T = (2.839881106229, 7.663478935125356, 17.598216593435)

"""
Counts of the expensive operations since the last `reset_tally!`: equation-of-state evaluations,
flows through a part, choking-point searches, junction Jacobians built, right-hand sides and ODE
Jacobians. `test/bench.jl` prints them: unlike wall time they are the same on every run, so they
show what a change did.
"""
mutable struct Tally
    on::Bool                      # counting is off unless asked for: see below
    eos::Int; edge_flow::Int; choke::Int; newton::Int; rhs::Int; jac::Int
end
const TALLY = Tally(false, 0, 0, 0, 0, 0, 0)
# Counting writes to one shared place, which simulations running on several threads would all be
# writing to at once. So it is switched on only by the benchmark, which runs one at a time.
reset_tally!(on = true) = (TALLY.on = on; TALLY.eos = TALLY.edge_flow = TALLY.choke = TALLY.newton = TALLY.rhs = TALLY.jac = 0)

"Single-phase properties at (ρ, T) straight from the Helmholtz energy. Valid in metastable states too."
function eos(ρ, T)
    TALLY.on && (TALLY.eos += 1)
    δ, τ = ρ / N2O_RHOC, N2O_TC / T
    ar = ar_d = ar_t = ar_dd = ar_dt = ar_tt = 0.0
    # Every temperature exponent is a multiple of 1/8, so one eighth root and whole-number powers
    # replace twelve real powers, and the three exp(-δ^l) are worked out once each. This function
    # is a third of a run's time.
    τ8 = sqrt(sqrt(sqrt(τ)))
    a2 = τ8 * τ8; a4 = a2 * a2; a7 = a4 * a2 * τ8; a10 = τ * a2; a12 = τ * a4; a16 = τ * τ
    a17 = a16 * τ8; a19 = a17 * a2; a28 = a16 * a12; a38 = a28 * a10; a52 = a28 * a16 * τ; a100 = a52 * a38 * a10
    τp = (a2, a10, a12, a2, a7, a19, a16, a17, a28, a52, a38, a100)         # τ^_T[k]: _T in eighths is 2, 10, 12, 2, 7, 19, 16, 17, 28, 52, 38, 100
    δ2 = δ * δ; δ3 = δ2 * δ; δ4 = δ2 * δ2
    δp = (δ, δ, δ, δ3, δ4 * δ3, δ, δ2, δ4 * δ, δ, δ, δ4, δ2)                # δ^_D[k]
    δls = (0.0, δ, δ2, δ3)
    es = (1.0, exp(-δ), exp(-δ2), exp(-δ3))
    for k in 1:12
        n, d, t, l = _N[k], _D[k], _T[k], _L[k]
        δl = δls[l+1]
        term = n * δp[k] * τp[k] * es[l+1]
        q = d - l * δl                      # δ·∂/∂δ acting on the term
        ar += term
        ar_d += term * q
        ar_dd += term * (q * (q - 1) - l^2 * δl)
        ar_t += term * t
        ar_tt += term * (t * (t - 1))
        ar_dt += term * (q * t)
    end
    iδ, iτ = 1 / δ, 1 / τ                   # the divisions, once each instead of in every term
    ar_d *= iδ; ar_dd *= iδ * iδ; ar_t *= iτ; ar_tt *= iτ * iτ; ar_dt *= iδ * iτ
    a0 = log(δ) + _A1 + _A2 * τ + 2.5 * log(τ)
    a0_t = _A2 + 2.5 / τ
    a0_tt = -2.5 / τ^2
    for i in 1:3
        n, t = _PE_N[i], _PE_T[i]
        e = exp(t * τ)
        a0 += n * log(1 - 1 / e)
        a0_t += n * t / (e - 1)
        a0_tt -= n * t^2 * e / (e - 1)^2
    end
    R = N2O_R
    P = ρ * R * T * (1 + δ * ar_d)
    u = R * T * τ * (a0_t + ar_t)
    s = R * (τ * (a0_t + ar_t) - a0 - ar)
    cv = -R * τ^2 * (a0_tt + ar_tt)
    dPdρ = R * T * (1 + 2δ * ar_d + δ^2 * ar_dd)
    dPdT = ρ * R * (1 + δ * ar_d - δ * τ * ar_dt)
    (; P, u, h = u + P / ρ, s, cv, dPdρ, dPdT)
end

cp_of(e, ρ, T) = e.cv + T * e.dPdT^2 / (ρ^2 * e.dPdρ)

# ---------------------------------------------------------------- saturation line

# Ancillary equations (CoolProp JSON) — only used as starting guesses for the Maxwell solve.
function _anc_rhoL(T)
    θ = 1 - T / N2O_TC
    n = (2.248014962727216, 1354.8814407985476, -227009237.43183637, 312935961.0720787,
         -139155040.7743576, 295107831.0429585)
    t = (0.374, 6.796, 13.027, 13.205, 14.362, 18.638)
    N2O_RHOC * (1 + sum(n[i] * θ^t[i] for i in 1:6))
end
function _anc_rhoV(T)
    θ = 1 - T / N2O_TC
    n = (-3.439484106428665, -35.87348642351084, 98.43310989645504, -173.87569213455475,
         30857.724063434925, -31916.967231256218)
    t = (0.428, 2.405, 3.115, 4.433, 8.312, 8.395)
    N2O_RHOC * exp(N2O_TC / T * sum(n[i] * θ^t[i] for i in 1:6))
end

"Coexisting densities at T: equal pressure and equal Gibbs energy, by Newton."
function _maxwell(T, ρl = _anc_rhoL(T), ρv = _anc_rhoV(T))
    for _ in 1:60
        l, v = eos(ρl, T), eos(ρv, T)
        f1 = l.P - v.P
        f2 = (l.h - T * l.s) - (v.h - T * v.s)
        # ∂g/∂ρ at constant T is (∂P/∂ρ)/ρ
        a, b, c, d = l.dPdρ, -v.dPdρ, l.dPdρ / ρl, -v.dPdρ / ρv
        det = a * d - b * c
        dl, dv = (f1 * d - b * f2) / det, (a * f2 - c * f1) / det
        ρl -= clamp(dl, -0.1ρl, 0.1ρl)
        ρv -= clamp(dv, -0.1ρv, 0.1ρv)
        abs(dl) < 1e-11ρl && abs(dv) < 1e-11ρv && break
    end
    ρl > N2O_RHOC > ρv || error("saturation solve collapsed at T=$T")
    ρl, ρv
end

# The table is uniform in z = (1 - T/Tc)^(1/3). Near the critical point properties go like
# (1 - T/Tc)^0.35, so in z they are nearly straight lines and a cubic spline stays accurate.
const SAT_TMAX = 309.40                    # top of the tabulated dome; above this: single phase
const _ZMIN = cbrt(1 - SAT_TMAX / N2O_TC)
const _ZMAX = cbrt(1 - N2O_TTRIP / N2O_TC)
const _NZ = 2048
const _DZ = (_ZMAX - _ZMIN) / (_NZ - 1)
const _COLS = (:T, :P, :ρl, :ρv, :ul, :uv, :hl, :hv, :sl, :sv, :cpl, :cpv)
const SAT = let
    tab = zeros(_NZ, length(_COLS))
    for i in 1:_NZ
        T = N2O_TC * (1 - (_ZMIN + (i - 1) * _DZ)^3)
        ρl, ρv = _maxwell(T)
        l, v = eos(ρl, T), eos(ρv, T)
        tab[i, :] .= (T, l.P, ρl, ρv, l.u, v.u, l.h, v.h, l.s, v.s, cp_of(l, ρl, T), cp_of(v, ρv, T))
    end
    tab
end
# Slope of every column at every row (per grid step), by fourth-order differences. With these the
# cubic between rows is within 1e-11 of the equation of state, against 1e-7 for slopes taken from
# the two neighbours and 1e-4 in the end intervals. That matters where an expansion enters the
# dome: the table takes over from the equation of state there, and a mismatch is a step in the flux.
const SAT_SLOPE = let
    m = similar(SAT)
    n = _NZ
    for j in axes(SAT, 2)
        y = view(SAT, :, j)
        for i in 3:n-2
            m[i, j] = (-y[i+2] + 8y[i+1] - 8y[i-1] + y[i-2]) / 12
        end
        m[1, j] = (-25y[1] + 48y[2] - 36y[3] + 16y[4] - 3y[5]) / 12
        m[2, j] = (-3y[1] - 10y[2] + 18y[3] - 6y[4] + y[5]) / 12
        m[n, j] = (25y[n] - 48y[n-1] + 36y[n-2] - 16y[n-3] + 3y[n-4]) / 12
        m[n-1, j] = (3y[n] + 10y[n-1] - 18y[n-2] + 6y[n-3] - y[n-4]) / 12
    end
    m
end
const SAT_PMAX = SAT[1, 2]
const SAT_PMIN = SAT[end, 2]

# The table is read as a cubic Hermite on the uniform z grid. Reading several columns at one
# z is the common case (a flash needs six), so finding the interval is separate from evaluating.

"Interval `i` and position `t` in [0, 1] within it for `z`. Not finite: (1, NaN), so values come out NaN."
@inline function _locate(z)
    isfinite(z) || return 1, NaN           # never throw: the solver rejects a step that goes NaN
    s = clamp((z - _ZMIN) / _DZ, 0.0, _NZ - 1.0)
    i = min(floor(Int, s), _NZ - 2)
    i + 1, s - i
end

"The cubic of column j on interval i: y(t) = ((a·t + b)·t + m0)·t + y0."
@inline function _cubic(j, i)
    @inbounds y0, y1, m0, m1 = SAT[i, j], SAT[i+1, j], SAT_SLOPE[i, j], SAT_SLOPE[i+1, j]
    2y0 - 2y1 + m0 + m1, -3y0 + 3y1 - 2m0 - m1, m0, y0
end

"Value of column j at a place found by `_locate`."
@inline function _val(j, i, t)
    a, b, m0, y0 = _cubic(j, i)
    ((a * t + b) * t + m0) * t + y0
end

"Value and d/dz of column j at z."
function _spline(j, z)
    i, t = _locate(z)
    a, b, m0, y0 = _cubic(j, i)
    ((a * t + b) * t + m0) * t + y0, ((3a * t + 2b) * t + m0) / _DZ
end

z_of_T(T) = cbrt(max(1 - T / N2O_TC, 0.0))

"""
A place on the saturation line. Its properties (`T`, `P`, `ρl`, `ρv`, `ul`, `uv`, `hl`, `hv`, `sl`,
`sv`, `cpl`, `cpv`) are read from the table when asked for, so a caller that wants two of the
twelve pays for two.
"""
struct SatPoint
    z::Float64
    i::Int
    t::Float64
end
const _SAT_FIELD = (T = 1, P = 2, ρl = 3, ρv = 4, ul = 5, uv = 6, hl = 7, hv = 8, sl = 9, sv = 10, cpl = 11, cpv = 12)
@inline Base.getproperty(s::SatPoint, f::Symbol) =
    f === :z ? getfield(s, :z) : _val(getfield(_SAT_FIELD, f), getfield(s, :i), getfield(s, :t))
Base.propertynames(::SatPoint) = (:z, keys(_SAT_FIELD)...)

"All saturation properties at z (see `sat_T`, `sat_P`)."
function sat_z(z)
    z = clamp(z, _ZMIN, _ZMAX)
    i, t = _locate(z)
    SatPoint(z, i, t)
end
sat_T(T) = sat_z(z_of_T(T))

"z at which monotone column j equals y (clamped to the table)."
function _invert(j, y)
    inc = SAT[end, j] > SAT[1, j]
    lo, hi = 1, _NZ
    (inc ? y <= SAT[1, j] : y >= SAT[1, j]) && return _ZMIN
    (inc ? y >= SAT[end, j] : y <= SAT[end, j]) && return _ZMAX
    while hi - lo > 1
        mid = (lo + hi) >> 1
        (SAT[mid, j] < y) == inc ? (lo = mid) : (hi = mid)
    end
    # Newton on the one cubic that covers this interval, from the straight-line estimate.
    a, b, m0, y0 = _cubic(j, lo)
    t = (y - y0) / (SAT[hi, j] - y0)
    for _ in 1:4
        t = clamp(t - ((((a * t + b) * t + m0) * t + y0) - y) / ((3a * t + 2b) * t + m0), 0.0, 1.0)
    end
    _ZMIN + (lo - 1 + t) * _DZ
end

Psat(T) = _spline(2, clamp(z_of_T(T), _ZMIN, _ZMAX))[1]
Tsat(P) = _spline(1, _invert(2, P))[1]
sat_P(P) = sat_z(_invert(2, P))
"dPsat/dT at T."
function dPsat_dT(T)
    z = clamp(z_of_T(T), _ZMIN, _ZMAX)
    _spline(2, z)[2] / (-3 * N2O_TC * z^2)
end

# ---------------------------------------------------------------- state from two properties

"A resolved fluid state. `x` is vapour mass fraction; single-phase states have x = 0 or 1."
struct Fluid
    T::Float64
    P::Float64
    ρ::Float64
    u::Float64
    h::Float64
    s::Float64
    x::Float64
    twophase::Bool
end

function _single(ρ, T)
    e = eos(ρ, T)
    Fluid(T, e.P, ρ, e.u, e.h, e.s, ρ > N2O_RHOC ? 0.0 : 1.0, false)
end

function _mix(st, x)
    v = x / st.ρv + (1 - x) / st.ρl
    u = st.ul + x * (st.uv - st.ul)
    Fluid(st.T, st.P, 1 / v, u, st.hl + x * (st.hv - st.hl), st.sl + x * (st.sv - st.sl), x, true)
end

"State at saturation pressure P with quality x."
state_Px(P, x) = _mix(sat_P(P), x)

# Internal energy of a saturated mixture of specific volume v at z, its z-derivative, and quality.
function _umix(z, v)
    (ρl, dρl), (ρv, dρv), (ul, dul), (uv, duv) = _spline(3, z), _spline(4, z), _spline(5, z), _spline(6, z)
    vl, vv = 1 / ρl, 1 / ρv
    x = (v - vl) / (vv - vl)
    dx = (dρl / ρl^2 * (1 - x) + x * dρv / ρv^2) / (vv - vl)
    ul + x * (uv - ul), dul + x * (duv - dul) + dx * (uv - ul), x
end

"""
State of a homogeneous volume from density and specific internal energy. Inside the dome it is an
equilibrium liquid/vapour mixture; outside, the single-phase EOS. Liquid-full and supercritical
states fall out of the same call. Below the triple point the temperature is clamped (no solid).
"""
function state_ρu(ρ, u)
    (ρ > 0 && isfinite(u)) || return Fluid(NaN, NaN, NaN, NaN, NaN, NaN, NaN, false)
    v = 1 / ρ
    # zb: where this density crosses the dome boundary (colder than that, it is two-phase)
    zb = ρ >= SAT[1, 3] ? _invert(3, ρ) : ρ <= SAT[1, 4] ? _invert(4, ρ) : _ZMIN
    inside = SAT[end, 4] < ρ < SAT[end, 3]
    if inside && u < _umix(zb, v)[1]
        lo, hi = zb, _ZMAX                 # u of the mixture falls as z rises (colder)
        _umix(hi, v)[1] >= u && return _mix(sat_z(hi), clamp(_umix(hi, v)[3], 0.0, 1.0))
        z = (lo + hi) / 2
        for _ in 1:60
            f, df, _ = _umix(z, v)
            f -= u
            zn = z - f / df
            # converged: take the last Newton step so the result is good to round-off (the stiff
            # solver differentiates this function numerically)
            abs(f) < 1e-9 * (abs(u) + 1e3) && (z = clamp(zn, zb, _ZMAX); break)
            f > 0 ? (lo = z) : (hi = z)
            z = lo < zn < hi ? zn : (lo + hi) / 2
        end
        return _mix(sat_z(z), clamp(_umix(z, v)[3], 0.0, 1.0))
    end
    T = inside ? _spline(1, zb)[1] : 300.0
    for _ in 1:60
        e = eos(ρ, T)
        dT = (e.u - u) / e.cv
        T = clamp(T - clamp(dT, -40.0, 40.0), 120.0, 2000.0)
        abs(dT) < 1e-10 * T && break
    end
    _single(ρ, T)
end

# Newton for (ρ, T) such that P(ρ,T) = P and h or s (ρ,T) = y.
function _solve_Py(P, y, ρ, T, use_s)
    e = eos(ρ, T)
    for _ in 1:80
        f1, f2 = e.P - P, (use_s ? e.s : e.h) - y
        a, b = e.dPdρ, e.dPdT
        c = use_s ? -e.dPdT / ρ^2 : (e.P - T * e.dPdT) / ρ^2 + e.dPdρ / ρ - e.P / ρ^2
        d = use_s ? e.cv / T : e.cv + e.dPdT / ρ
        det = a * d - b * c
        dρ, dT = (f1 * d - b * f2) / det, (a * f2 - c * f1) / det
        ρ = max(ρ - clamp(dρ, -0.3ρ, 0.3ρ), 1e-6)
        T = clamp(T - clamp(dT, -25.0, 25.0), 120.0, 2000.0)
        e = eos(ρ, T)
        # Newton doubles the number of correct digits with every step, so after a step this small
        # the state just evaluated is exact to round-off and needs no further step to confirm it.
        abs(dρ) < 1e-8ρ && abs(dT) < 1e-8T && break
    end
    Fluid(T, e.P, ρ, e.u, e.h, e.s, ρ > N2O_RHOC ? 0.0 : 1.0, false)
end

# `near`, if given, is a state already solved close by (the same part a moment ago). When it is
# single phase on the same side, Newton starts from it and needs one or two steps instead of four.
function _state_Py(P, y, use_s, near = nothing)
    warm(liquid) = near !== nothing && !near.twophase && 0.8 < near.P / P < 1.25 && (near.ρ > N2O_RHOC) == liquid
    if P < SAT_PMAX
        st = sat_P(max(P, SAT_PMIN))
        yl, yv = use_s ? (st.sl, st.sv) : (st.hl, st.hv)
        yl <= y <= yv && return _mix(st, (y - yl) / (yv - yl))
        warm(y < yl) && return _solve_Py(P, y, near.ρ, near.T, use_s)
        return y < yl ? _solve_Py(P, y, st.ρl, st.T, use_s) : _solve_Py(P, y, st.ρv, st.T, use_s)
    end
    st = sat_z(_ZMIN)
    ymid = use_s ? (st.sl + st.sv) / 2 : (st.hl + st.hv) / 2
    warm(y < ymid) && return _solve_Py(P, y, near.ρ, near.T, use_s)
    y < ymid ? _solve_Py(P, y, 1.3st.ρl, st.T - 2, use_s) : _solve_Py(P, y, 0.7st.ρv, st.T + 10, use_s)
end

"State from pressure and specific enthalpy (an isenthalpic flash)."
state_Ph(P, h, near = nothing) = _state_Py(P, h, false, near)
"State from pressure and specific entropy (an isentropic expansion)."
state_Ps(P, s) = _state_Py(P, s, true)

"Density of single-phase vapour at (P, T) below saturation pressure; saturated vapour if P ≥ Psat."
function vapour_PT(P, T)
    ρ = T < SAT_TMAX ? min(P / (N2O_R * T), sat_T(T).ρv) : P / (N2O_R * T)
    T < SAT_TMAX && P >= Psat(T) && return _single(sat_T(T).ρv, T)
    for _ in 1:60
        e = eos(ρ, T)
        dρ = (e.P - P) / e.dPdρ
        ρ = max(ρ - clamp(dρ, -0.3ρ, 0.3ρ), 1e-9)
        abs(dρ) < 1e-11ρ && break
    end
    _single(ρ, T)
end

# ---------------------------------------------------------------- transport (saturation line)
#
# ESDU 91022 (1991), Eqs 4.9-4.12 and Table 7.1; they use Tc = 309.57 K. ESDU says not to
# extrapolate: viscosity is given for -90..30 °C and conductivity for -90..10 °C, so outside those
# ranges the nearest in-range value is used. Most of the underlying table values are flagged by
# ESDU as not backed by measurements. Each heat-transfer use has a multiplier for that reason.
"Surface tension, N/m: σ = 0.07087 (1 - T/Tc)^1.204 (Mulero, Cachadiña and Parra 2012, as carried by CoolProp)."
surface_tension(T) = 0.07087 * max(1 - T / N2O_TC, 0.0)^1.204

const _ESDU_TC = 309.57
"Saturated-liquid viscosity, Pa·s."
function visc_l(T)
    θ = (_ESDU_TC - 5.24) / (clamp(T, 182.33, 303.15) - 5.24)
    0.0293423e-3 * exp(1.6089 * cbrt(θ - 1) + 2.0439 * (θ - 1)^(4 / 3))
end
"Saturated-vapour viscosity, Pa·s."
function visc_v(T)
    a = _ESDU_TC / clamp(T, 182.33, 303.15) - 1
    1e-6 * exp(3.3281 - 1.18237 * cbrt(a) - 0.055155 * a^(4 / 3))
end
"Saturated-liquid thermal conductivity, W/m/K."
function cond_l(T)
    a = 1 - clamp(T, 182.33, 283.15) / _ESDU_TC
    72.35e-3 * (1 + 1.5 * cbrt(a) - 3.5 * a^(2 / 3) + 4.5 * a)
end
"Saturated-vapour thermal conductivity, W/m/K."
function cond_v(T)
    a = 1 - clamp(T, 182.33, 283.15) / _ESDU_TC
    1e-3 * exp(-7.08870 - 0.276962 * a^(-2 / 3) + 2.88672 * a^(-1 / 3) + 16.6116 * cbrt(a) - 11.8221 * a^(2 / 3))
end
