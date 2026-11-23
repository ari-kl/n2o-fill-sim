# Control policies: a controller that reads the instruments at a fixed tick and commands valves and
# heaters while a study's run goes on, in place of a schedule fixed in advance. See docs/STUDIES.md,
# "Control policies". Nothing here is used by a study without a [controller] table.
#
# The controller works only through what the network already offers: it moves the run from tick to
# tick with `advance!` (which lets the integrator run ahead and interpolates back to the tick, so a
# tick at which nothing is commanded costs no restart), reads the recorded history, and puts its
# commands into the pending schedule with `schedule!`, where `advance!` applies them like any other
# action.

"One command a controller sent: when, which policy, to which part, what."
struct Command
    t::Float64
    policy::String
    part::String
    action::String
    value::Float64
end

"""
A policy reads some instruments at each tick and may command some parts. Each kind is a struct
with at least `name` and `enabled`, and three methods: `_reads(p)` (the readings it looks at),
`_actuators(p)` (the parts it commands) and `_tick!(p, c, sim, t)`, which sends commands through
`_command!`. A setting that may differ from step to step is a `PerStep`.
"""
abstract type Policy end

# A setting by step: "" holds the value for every step that has none of its own.
const PerStep = Dict{String,Float64}

mutable struct Hysteresis <: Policy
    name::String
    enabled::Bool
    reading::String
    actuator::String
    open_above::PerStep
    close_below::PerStep
    min_open_s::PerStep
    min_shut_s::PerStep
    open::Bool                                # what it last commanded (it starts believing the valve shut)
    t_last::Float64                           # when it last commanded anything
end

mutable struct PulseWait <: Policy
    name::String
    enabled::Bool
    reading::String
    actuator::String
    above::PerStep
    pulse_s::PerStep
    wait_s::PerStep
    phase::Symbol                             # :idle, :pulse (valve open) or :wait (shut, not to reopen yet)
    t_phase::Float64
end

mutable struct HeaterHold <: Policy
    name::String
    enabled::Bool
    heater::String
    on::PerStep                               # 1 = switched on in that step, 0 = off
    commanded::Int                            # -1 before its first command, else 0 or 1
end

"A policy written in Julia (`kind = \"julia\"`); see `_julia_policy`."
mutable struct JuliaPolicy <: Policy
    name::String
    enabled::Bool
    fn::Any
    reads::Vector{String}
    actuators::Vector{String}
    settings::Dict{String,Any}
    state::Dict{String,Any}
end

_reads(p::Union{Hysteresis,PulseWait}) = [p.reading]
_reads(p::HeaterHold) = String[]
_reads(p::JuliaPolicy) = p.reads
_actuators(p::Union{Hysteresis,PulseWait}) = [p.actuator]
_actuators(p::HeaterHold) = [p.heater]
_actuators(p::JuliaPolicy) = p.actuators

"""
A controller: its tick, the delays, what it may read and how, its policies in tick order, and what
it has done. Built afresh for every run from the study's `[controller]` table.
"""
mutable struct Controller
    dt::Float64
    latency::Float64                          # a reading is what the instrument showed this long ago
    delay::Float64                            # a command reaches its part this long after it is sent
    source::Dict{String,String}               # reading => the recorded signal it is taken from
    stand_ins::Vector{String}                 # true signals a policy is allowed to read as if instruments
    noise::Dict{String,Float64}               # reading => standard deviation of the noise added to it (SI)
    rng::Dict{String,_Lcg}                    # one stream of noise per noisy reading
    policies::Vector{Policy}
    step::String                              # the step now running ("" without a sequence)
    obs::Dict{String,Float64}                 # this tick's readings
    log::Vector{Command}
    escape_hatch::Bool                        # a Julia policy took part
    signals::Vector{String}                   # the ctl.* signals recorded beside the run's history
    values::Vector{Float64}                   # their values since the last tick
end

const _CONTROLLER_KEYS = ("dt", "latency_s", "actuation_delay_s", "stand_ins", "noise", "policy")
const _POLICY_KEYS = Dict(
    "hysteresis" => ("name", "kind", "enabled", "reading", "actuator", "open_above", "close_below", "min_open_s", "min_shut_s"),
    "pulse_wait" => ("name", "kind", "enabled", "reading", "actuator", "above", "pulse_s", "wait_s"),
    "heater_hold" => ("name", "kind", "enabled", "heater", "on"),
    "julia" => ("name", "kind", "enabled", "file", "reads", "actuators", "settings"))

# Each named setting may be given plainly in SI, or by its name with _psig (a gauge pressure) or _C
# (a temperature), as metric arguments are. Only those that are levels take the suffixes.
const _LEVELS = ("open_above", "close_below", "above")

# A setting as a PerStep: one value for every step, or a table of values by step name.
function _setting(p, key, p_atm, where; default = nothing)
    keys_ = key in _LEVELS ? (key, key * "_psig", key * "_C") : (key,)
    given = [k for k in keys_ if haskey(p, k)]
    length(given) > 1 && error("$where: give one of $(join(given, " and "))")
    if isempty(given)
        default === nothing && error("$where needs `$key`" * (length(keys_) > 1 ? " (or $(key)_psig, $(key)_C)" : ""))
        return PerStep("" => default)
    end
    k = given[1]
    conv(v) = (v isa Real || error("$where: $k must be a number or a table of numbers by step"); x = Float64(v);
               endswith(k, "_psig") ? x * PSI + p_atm : endswith(k, "_C") ? x + 273.15 : x)
    v = p[k]
    v isa AbstractDict ? PerStep(String(s) => conv(x) for (s, x) in v) : PerStep("" => conv(v))
end

# A setting's value in a step. Every step a policy can run in has been checked to have one.
_at(v::PerStep, step) = get(v, step, get(v, "", NaN))

"""
    _controller(spec, scene, sim; seed = "", allow_julia = false, fns) -> Controller

Build a run's controller from its `[controller]` table (overrides already written into it) for a
run of `scene`, whose `Sim` is `sim`. Everything a policy reads is checked here: an instrument's
reading (a Gauge's `P_read` or `T_read`, a LoadCell's `m_read`), or a true signal listed in
`stand_ins`. `seed` is what the noise is seeded from. `fns` holds the functions of Julia policies,
loaded once when the study is read.
"""
function _controller(spec, scene::Scene, sim::Sim; seed = "", allow_julia = false, fns = Dict{String,Any}(), where = "[controller]")
    for k in keys(spec)
        k in _CONTROLLER_KEYS || error("$where: unknown key \"$k\". Keys: $(join(_CONTROLLER_KEYS, ", "))")
    end
    p_atm = scene.settings["P_atm_kPa"] * 1e3
    dt = Float64(get(spec, "dt", 0.2))
    dt > 0 || error("$where: dt must be positive")
    latency, delay = Float64(get(spec, "latency_s", 0.0)), Float64(get(spec, "actuation_delay_s", 0.0))
    latency >= 0 && delay >= 0 || error("$where: latency_s and actuation_delay_s cannot be negative")
    stand_ins = String.(get(spec, "stand_ins", String[]))
    for s in stand_ins
        haskey(sim.hist, s) || error("$where: stand_ins names $s, which this scene does not record. Signals: $(join(sort!(collect(keys(sim.hist))), ", "))")
    end
    policies = Policy[]
    for (i, p) in enumerate(get(spec, "policy", []))
        p isa AbstractDict || error("$where: each [[controller.policy]] must be a table")
        name = String(get(p, "name", ""))
        w = "$where: policy $(isempty(name) ? i : name)"
        (isempty(name) || occursin('.', name)) && error("$w needs a `name` without a '.'")
        kind = String(get(p, "kind", ""))
        haskey(_POLICY_KEYS, kind) || error("$w: kind must be one of $(join(sort!(collect(keys(_POLICY_KEYS))), ", "))")
        allowed = _POLICY_KEYS[kind]
        for k in keys(p)
            base = replace(k, r"_(psig|C)$" => "")
            (k in allowed || (base in _LEVELS && base in allowed)) || error("$w: unknown key \"$k\". Keys: $(join(allowed, ", "))")
        end
        enabled = Bool(get(p, "enabled", true))
        s(key; default = nothing) = _setting(p, key, p_atm, w; default)
        if kind == "hysteresis"
            q = Hysteresis(name, enabled, String(get(p, "reading", "")), String(get(p, "actuator", "")), s("open_above"), s("close_below"),
                           s("min_open_s"; default = 0.0), s("min_shut_s"; default = 0.0), false, -Inf)
        elseif kind == "pulse_wait"
            q = PulseWait(name, enabled, String(get(p, "reading", "")), String(get(p, "actuator", "")), s("above"), s("pulse_s"), s("wait_s"), :idle, -Inf)
        elseif kind == "heater_hold"
            q = HeaterHold(name, enabled, String(get(p, "heater", "")), s("on"; default = 1.0), -1)
        else
            allow_julia || error("$w: a Julia policy runs code from a file; it is only run with --allow-julia-policy")
            haskey(fns, p["file"]) || error("$w: $(p["file"]) was not loaded")
            q = JuliaPolicy(name, enabled, fns[p["file"]], String.(get(p, "reads", String[])), String.(get(p, "actuators", String[])),
                            Dict{String,Any}(get(p, "settings", Dict())), Dict{String,Any}())
        end
        any(o -> o.name == name, policies) && error("$w: two policies are called $name")
        push!(policies, q)
    end
    source = Dict{String,String}()
    for q in policies, r in _reads(q)
        source[r] = _reading_source(r, scene, sim, stand_ins, "$where: policy $(q.name)")
    end
    for q in policies, a in _actuators(q)
        _check_actuator(q, a, scene, "$where: policy $(q.name)")
    end
    for s in stand_ins
        s in values(source) || error("$where: stand_ins names $s, which nothing reads")
    end
    noise = Dict{String,Float64}()
    for (r, v) in get(spec, "noise", Dict())
        haskey(source, r) || error("$where: noise is given for $r, which no policy reads")
        v isa AbstractDict && length(v) == 1 || error("$where: noise for $r must be a table like { sigma_psi = 2 } or { sigma = 1000 }")
        k, x = first(v)
        k in ("sigma", "sigma_psi", "sigma_K") || error("$where: noise for $r: give sigma (SI), sigma_psi or sigma_K, not $k")
        x >= 0 || error("$where: noise for $r cannot be negative")
        noise[r] = Float64(x) * (k == "sigma_psi" ? PSI : 1.0)
    end
    rng = Dict(r => _Lcg(_fnv1a(seed * "|" * r)) for r in keys(noise))
    signals = String[]
    for q in policies
        append!(signals, ["ctl.$(q.name).enabled", "ctl.$(q.name).cmd"])
        isempty(_reads(q)) || push!(signals, "ctl.$(q.name).obs")
    end
    Controller(dt, latency, delay, source, stand_ins, noise, rng, policies, "", Dict{String,Float64}(), Command[],
               any(q -> q isa JuliaPolicy, policies), signals, fill(NaN, length(signals)))
end

# Where a policy's reading comes from. A gauge with no sensor parameter set records no reading;
# its P_read and T_read are then its P and T, exactly as in scoring (`_scored_signal`).
function _reading_source(r, scene, sim, stand_ins, where)
    r in stand_ins && return r
    m = match(r"^(.+)\.(P_read|T_read|m_read)$", r)
    i = m === nothing ? nothing : find_part(scene, m[1])
    kind = i === nothing ? "" : scene.parts[i].kind
    if kind == "Gauge" && m[2] != "m_read"
        haskey(sim.hist, r) && return r
        get(sim.net.gauges, m[1], 0) != 0 && return m[1] * "." * m[2][1:1]
        error("$where: $(m[1]) is not wired to anything, so it has no reading")
    elseif kind == "LoadCell" && m[2] == "m_read" && haskey(sim.hist, r)
        return r
    end
    error("$where reads $r, which is not an instrument's reading (a Gauge's P_read or T_read, a LoadCell's m_read). " *
          "To let it read a true signal, list it in [controller] stand_ins")
end

function _check_actuator(q, a, scene, where)
    i = find_part(scene, a)
    i === nothing && error("$where commands $a, which is not in the scene")
    kind = scene.parts[i].kind
    if q isa HeaterHold
        kind in VESSEL_KINDS || error("$where: $a is a $kind, not a Bottle or Tank with a heater")
    elseif q isa Union{Hysteresis,PulseWait}
        kind == "Valve" || error("$where: $a is a $kind, not a Valve")
        get(scene.parts[i].p, "manual", 0.0) > 0 && error("$where: $a is a hand valve, which takes no commands")
    end
end

# FNV-1a, 64 bits: a hash that is the same on every machine and every Julia version (`hash` is not).
function _fnv1a(s::AbstractString)
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ b) * 0x00000100000001b3
    end
    h
end

# What a run's noise is seeded from: its parameters other than the controller's, and its case. Two
# runs that differ only in controller settings therefore see the same noise.
_noise_seed(overrides, case) = join(["$k=$(string(reinterpret(UInt64, Float64(v)); base = 16))" for (k, v) in sort(collect(overrides); by = first)], ";") * "|" * case

# A draw from the standard normal distribution (Box-Muller; 1 - u is in (0, 1], so its log is finite).
_gauss!(g::_Lcg) = (u1 = _next!(g); u2 = _next!(g); sqrt(-2 * log(1 - u1)) * cos(2π * u2))

"""
    _observe!(c, sim, t)

Take every declared reading for the tick at `t`: the recorded signal at `t - latency`, plus noise.
Each noisy reading's stream is drawn once a tick whatever the policies do with it.
"""
function _observe!(c::Controller, sim::Sim, t)
    for (r, src) in c.source
        v = signal_at(sim, src, max(t - c.latency, 0.0))
        haskey(c.noise, r) && (v += c.noise[r] * _gauss!(c.rng[r]))
        c.obs[r] = v
    end
end

"Send a command: logged now, applied to the part `actuation_delay_s` later."
function _command!(c::Controller, sim::Sim, p::Policy, t, part, action, value = 0.0)
    push!(c.log, Command(t, p.name, part, action, value))
    schedule!(sim, Action(t + c.delay, part, action, value))
end

function _tick!(p::Hysteresis, c::Controller, sim::Sim, t)
    x = c.obs[p.reading]
    if !p.open && x > _at(p.open_above, c.step) && t - p.t_last >= _at(p.min_shut_s, c.step) - 1e-9
        _command!(c, sim, p, t, p.actuator, "open"); p.open, p.t_last = true, t
    elseif p.open && x < _at(p.close_below, c.step) && t - p.t_last >= _at(p.min_open_s, c.step) - 1e-9
        _command!(c, sim, p, t, p.actuator, "close"); p.open, p.t_last = false, t
    end
end

function _tick!(p::PulseWait, c::Controller, sim::Sim, t)
    if p.phase === :pulse && t - p.t_phase >= _at(p.pulse_s, c.step) - 1e-9
        _command!(c, sim, p, t, p.actuator, "close"); p.phase, p.t_phase = :wait, t
    elseif p.phase === :wait && t - p.t_phase >= _at(p.wait_s, c.step) - 1e-9
        p.phase = :idle
    end
    if p.phase === :idle && c.obs[p.reading] > _at(p.above, c.step)
        _command!(c, sim, p, t, p.actuator, "open"); p.phase, p.t_phase = :pulse, t
    end
end

# A heater is commanded only when what it should be changes: switching one restarts the integrator.
function _tick!(p::HeaterHold, c::Controller, sim::Sim, t)
    want = _at(p.on, c.step) > 0.5 ? 1 : 0
    want == p.commanded && return
    _command!(c, sim, p, t, p.heater, want == 1 ? "on" : "off")
    p.commanded = want
end

function _tick!(p::JuliaPolicy, c::Controller, sim::Sim, t)
    obs = (; t, step = c.step, readings = Dict(r => c.obs[r] for r in p.reads), settings = p.settings)
    for cmd in Base.invokelatest(p.fn, obs, p.state)
        part, action = String(cmd[1]), String(cmd[2])
        part in p.actuators || error("Julia policy $(p.name) commanded $part, which is not in its `actuators`")
        _command!(c, sim, p, t, part, action, length(cmd) >= 3 ? Float64(cmd[3]) : 0.0)
    end
end

"""
Load the function a Julia policy's file defines: the file's last expression, a function of
`(obs, state)` returning a list of commands `(part, action)` or `(part, action, value)`. `obs`
has `t`, `step`, `readings` (only those listed in its `reads`) and `settings`; `state` is a Dict of
its own, kept from tick to tick.
"""
function _julia_policy(path)
    fn = Base.include(Module(:JuliaPolicy), path)
    Base.invokelatest(applicable, fn, (; t = 0.0), Dict{String,Any}()) || error("$path must end with a function of (obs, state)")
    fn
end

# The last command a policy sent, as a number: 1 for open or on, 0 for close or off, a set value as itself.
function _last_cmd(c::Controller, p::Policy)
    k = findlast(x -> x.policy == p.name, c.log)
    k === nothing && return 0.0
    a = c.log[k]
    a.action in ("open", "on") ? 1.0 : a.action in ("close", "off", "deenergise") ? 0.0 : a.value
end

# The ctl.* signals have a value per recorded sample, held from one tick to the next. Samples
# recorded before `t` get the values as they stood; `upto = Inf` pads every sample there is.
function _pad!(c::Controller, sim::Sim, upto = Inf)
    n = searchsortedlast(sim.t_hist, upto - 1e-9)
    for (k, name) in enumerate(c.signals)
        h = get!(sim.hist, name, Float64[])
        while length(h) < n
            push!(h, c.values[k])
        end
    end
end

function _update_signals!(c::Controller)
    k = 0
    for p in c.policies
        c.values[k+1] = p.enabled ? 1.0 : 0.0
        c.values[k+2] = _last_cmd(c, p)
        k += 2
        isempty(_reads(p)) && continue
        c.values[k+1] = c.obs[_reads(p)[1]]
        k += 1
    end
end

# One tick: read, then each enabled policy in order.
function _tick!(c::Controller, sim::Sim, t)
    _pad!(c, sim, t)
    _observe!(c, sim, t)
    for p in c.policies
        p.enabled && _tick!(p, c, sim, t)
    end
    _update_signals!(c)
    _pad!(c, sim)
end

"""
    run_controlled!(sim, c, until; wall_limit, work_limit) -> sim

Run to `until` with the controller ticking every `dt` seconds from t = 0, or until the solver
fails or a limit is reached (as `run_until!`). At each tick the run is brought to the tick by
`advance!`, a sample is recorded there if the solver's own samples did not fall on it (so that a
reading is never older than the tick), and the policies act.
"""
function run_controlled!(sim::Sim, c::Controller, until; wall_limit = Inf, work_limit = typemax(Int))
    start = Base.time()
    reach(t) = (run_until!(sim, t; wall_limit = wall_limit - (Base.time() - start), work_limit); simtime(sim) >= t - 1e-9)
    try
        k = 0
        while k * c.dt <= until + 1e-9
            t = k * c.dt
            reach(t) || return sim
            sim.last_sample < simtime(sim) - 1e-9 && sample!(sim)
            _tick!(c, sim, t)
            k += 1
        end
        reach(until)
    finally
        _pad!(c, sim)                          # whatever stopped the run, the ctl.* signals cover it
    end
    sim
end
