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
    commanded::Int                            # 1 on, 0 off: what it last commanded (it starts believing the heater off, as a run starts)
end

"Marks the moment a reading has stayed inside a band for `hold_s`, once."
mutable struct BandWatch <: Policy
    name::String
    enabled::Bool
    reading::String
    lo::PerStep
    hi::PerStep
    hold_s::PerStep
    mark::String
    since::Float64                            # when the reading last came into the band (NaN while outside)
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

_reads(p::Union{Hysteresis,PulseWait,BandWatch}) = [p.reading]
_reads(p::HeaterHold) = String[]
_reads(p::JuliaPolicy) = p.reads
_actuators(p::Union{Hysteresis,PulseWait}) = [p.actuator]
_actuators(p::HeaterHold) = [p.heater]
_actuators(p::BandWatch) = String[]
_actuators(p::JuliaPolicy) = p.actuators

# Conditions: a step's `until`, compiled when the controller is built. Comparisons and rates are on
# declared readings, `elapsed_ge` is time in the step, `held` asks for a condition to have stayed
# true for a while, `mark` for a mark to have happened. Rates and holds keep what they have seen
# since their step began.
abstract type Cond end
struct Compare <: Cond
    reading::String
    op::Symbol                                # :gt, :ge, :lt or :le
    v::Float64
end
struct Elapsed <: Cond
    s::Float64
end
mutable struct Rate <: Cond
    reading::String
    window::Float64
    op::Symbol
    v::Float64
    seen::Vector{Tuple{Float64,Float64}}      # (t, reading) at each tick within the window
end
mutable struct Held <: Cond
    c::Cond
    s::Float64
    since::Float64
end
struct MarkCond <: Cond
    mark::String
    plus::Float64
end
struct AnyOf <: Cond
    cs::Vector{Cond}
end
struct AllOf <: Cond
    cs::Vector{Cond}
end

"One step of a sequence: what it does on entering and leaving, and what ends it."
struct Step
    name::String
    on_enter::Vector{Tuple{String,String,Float64}}    # (part, action, value)
    on_exit::Vector{Tuple{String,String,Float64}}
    until::Union{Nothing,Cond}
    duration::Float64                         # Inf for none
    timeout::Float64
    operator_delay::Float64                   # from `until` being met to the step ending
    enable::Vector{String}                    # policies switched on and off on entering it
    disable::Vector{String}
end

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
    step::String                              # the step now running, or the last one ("" without a sequence)
    steps::Vector{Step}
    k::Int                                    # which step is running: 0 before the first tick, length(steps) + 1 once all have ended
    t_start::Dict{String,Float64}
    t_end::Dict{String,Float64}
    timed_out::Dict{String,Bool}
    due_end::Float64                          # when the running step ends, once its `until` has been met
    marks::Dict{String,Float64}               # mark => when it happened
    mark_names::Vector{String}
    end_on::Tuple{String,String}              # (step, "start" or "end"): the run ends at that edge; ("", "") for none
    ended::Bool
    obs::Dict{String,Float64}                 # this tick's readings
    log::Vector{Command}
    escape_hatch::Bool                        # a Julia policy took part
    signals::Vector{String}                   # the ctl.* signals recorded beside the run's history
    values::Vector{Float64}                   # their values since the last tick
end

const _CONTROLLER_KEYS = ("dt", "latency_s", "actuation_delay_s", "stand_ins", "noise", "policy", "step")
const _STEP_KEYS = ("name", "on_enter", "on_exit", "until", "duration_s", "timeout_s", "operator_delay_s", "enable", "disable")
const _POLICY_KEYS = Dict(
    "hysteresis" => ("name", "kind", "enabled", "reading", "actuator", "open_above", "close_below", "min_open_s", "min_shut_s"),
    "pulse_wait" => ("name", "kind", "enabled", "reading", "actuator", "above", "pulse_s", "wait_s"),
    "heater_hold" => ("name", "kind", "enabled", "heater", "on"),
    "band_watch" => ("name", "kind", "enabled", "reading", "lo", "hi", "hold_s", "mark"),
    "julia" => ("name", "kind", "enabled", "file", "reads", "actuators", "settings"))

# Each named setting may be given plainly in SI, or by its name with _psig (a gauge pressure) or _C
# (a temperature), as metric arguments are. Only those that are levels take the suffixes.
const _LEVELS = ("open_above", "close_below", "above", "lo", "hi")

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
    v isa AbstractDict ? PerStep((s == "default" ? "" : String(s)) => conv(x) for (s, x) in v) : PerStep("" => conv(v))
end

# A level in a condition: `ge = 3.0e6`, `ge_psig = 420`, `lt_C = 5`.
_level(k, v, p_atm) = (x = Float64(v); endswith(k, "_psig") ? x * PSI + p_atm : endswith(k, "_C") ? x + 273.15 : x)

const _OPS = ("gt", "ge", "lt", "le")

# Compile a condition table; every reading it names is added to `reads`.
function _cond(d, p_atm, where, reads, marks)
    d isa AbstractDict || error("$where: a condition must be a table, like { reading = \"PT1.P_read\", ge_psig = 500 }")
    has(k) = haskey(d, k)
    if has("any") || has("all")
        length(d) == 1 || error("$where: `any` or `all` takes nothing beside its list")
        cs = Cond[_cond(x, p_atm, where, reads, marks) for x in (has("any") ? d["any"] : d["all"])]
        isempty(cs) && error("$where: an empty `any` or `all`")
        return has("any") ? AnyOf(cs) : AllOf(cs)
    elseif has("held")
        Set(keys(d)) == Set(["held", "for_s"]) || error("$where: `held` takes a condition and `for_s`")
        return Held(_cond(d["held"], p_atm, where, reads, marks), Float64(d["for_s"]), NaN)
    elseif has("elapsed_ge")
        length(d) == 1 || error("$where: `elapsed_ge` takes nothing beside it")
        return Elapsed(Float64(d["elapsed_ge"]))
    elseif has("mark")
        issubset(keys(d), ["mark", "plus_s"]) || error("$where: `mark` takes only `plus_s` beside it")
        d["mark"] in marks || error("$where: there is no mark called $(d["mark"]). Marks: $(join(marks, ", "))")
        return MarkCond(String(d["mark"]), Float64(get(d, "plus_s", 0.0)))
    end
    rate = has("rate")
    r = String(rate ? d["rate"] : get(d, "reading", ""))
    isempty(r) && error("$where: a condition needs `reading` (or `rate`, `elapsed_ge`, `held`, `mark`, `any`, `all`)")
    # a rate is a difference per second: plain SI, or _psi for psi per second (never a gauge zero)
    ops = [(k, op) for op in _OPS for k in (rate ? (op, op * "_psi") : (op, op * "_psig", op * "_C")) if has(k)]
    length(ops) == 1 || error("$where: a condition on $r needs exactly one of $(join(_OPS, ", "))" * (rate ? " (or with _psi)" : " (or with _psig, _C)"))
    (k, op) = ops[1]
    for key in keys(d)
        key in (rate ? ("rate", "window_s", k) : ("reading", k)) || error("$where: unknown key \"$key\" in a condition on $r")
    end
    push!(reads, r)
    rate || return Compare(r, Symbol(op), _level(k, d[k], p_atm))
    haskey(d, "window_s") || error("$where: a rate needs `window_s`")
    Rate(r, Float64(d["window_s"]), Symbol(op), Float64(d[k]) * (endswith(k, "_psi") ? PSI : 1.0), Tuple{Float64,Float64}[])
end

_compare(x, op, v) = op === :gt ? x > v : op === :ge ? x >= v : op === :lt ? x < v : x <= v

# Whether a condition holds at this tick. Every part of it is looked at, so that rates and holds
# keep their record whatever their neighbours say.
_holds(q::Compare, c, t) = _compare(c.obs[q.reading], q.op, q.v)
_holds(q::Elapsed, c, t) = t - c.t_start[c.step] >= q.s - 1e-9
_holds(q::MarkCond, c, t) = haskey(c.marks, q.mark) && t >= c.marks[q.mark] + q.plus - 1e-9
_holds(q::AnyOf, c, t) = (v = [_holds(x, c, t) for x in q.cs]; any(v))
_holds(q::AllOf, c, t) = (v = [_holds(x, c, t) for x in q.cs]; all(v))
function _holds(q::Held, c, t)
    _holds(q.c, c, t) ? (isnan(q.since) && (q.since = t)) : (q.since = NaN)
    !isnan(q.since) && t - q.since >= q.s - 1e-9
end
# The rate over the last `window_s`: from the earliest tick inside the window to now.
function _holds(q::Rate, c, t)
    x = c.obs[q.reading]
    push!(q.seen, (t, x))
    filter!(p -> p[1] >= t - q.window - 1e-9, q.seen)
    t0, x0 = q.seen[1]
    _compare((x - x0) / max(t - t0, 1e-6), q.op, q.v)
end
_reset!(q::Rate) = empty!(q.seen)
_reset!(q::Held) = (q.since = NaN; _reset!(q.c))
_reset!(q::Union{AnyOf,AllOf}) = foreach(_reset!, q.cs)
_reset!(q::Cond) = nothing

function _actions(list, scene, where)
    map(list) do a
        a isa AbstractDict && issubset(keys(a), ["part", "action", "value"]) && haskey(a, "part") && haskey(a, "action") ||
            error("$where: an action is a table like { part = \"AV1\", action = \"open\" } (value for set)")
        find_part(scene, a["part"]) === nothing && error("$where: $(a["part"]) is not in the scene")
        (String(a["part"]), String(a["action"]), Float64(get(a, "value", 0.0)))
    end
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
function _controller(spec, scene::Scene, sim::Sim; seed = "", allow_julia = false, fns = Dict{String,Any}(), where = "[controller]", end_on = ("", ""))
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
            q = HeaterHold(name, enabled, String(get(p, "heater", "")), s("on"; default = 1.0), 0)
        elseif kind == "band_watch"
            q = BandWatch(name, enabled, String(get(p, "reading", "")), s("lo"), s("hi"), s("hold_s"; default = 0.0), String(get(p, "mark", name)), NaN)
        else
            allow_julia || error("$w: a Julia policy runs code from a file; it is only run with --allow-julia-policy")
            haskey(fns, p["file"]) || error("$w: $(p["file"]) was not loaded")
            q = JuliaPolicy(name, enabled, fns[p["file"]], String.(get(p, "reads", String[])), String.(get(p, "actuators", String[])),
                            Dict{String,Any}(get(p, "settings", Dict())), Dict{String,Any}())
        end
        any(o -> o.name == name, policies) && error("$w: two policies are called $name")
        push!(policies, q)
    end
    mark_names = [q.mark for q in policies if q isa BandWatch]
    allunique(mark_names) || error("$where: two band_watch policies make the same mark")
    # the sequence
    steps, cond_reads = Step[], String[]
    for (i, d) in enumerate(get(spec, "step", []))
        d isa AbstractDict || error("$where: each [[controller.step]] must be a table")
        name = String(get(d, "name", ""))
        w = "$where: step $(isempty(name) ? i : name)"
        (isempty(name) || occursin('.', name) || name == "default") && error("$w needs a `name` without a '.' (and not \"default\")")
        any(x -> x.name == name, steps) && error("$w: two steps are called $name")
        for k in keys(d)
            k in _STEP_KEYS || error("$w: unknown key \"$k\". Keys: $(join(_STEP_KEYS, ", "))")
        end
        on_off = [String.(get(d, k, String[])) for k in ("enable", "disable")]
        for n in vcat(on_off...)
            any(q -> q.name == n, policies) || error("$w enables or disables $n, which is not a policy")
        end
        push!(steps, Step(name, _actions(get(d, "on_enter", []), scene, w), _actions(get(d, "on_exit", []), scene, w),
                          haskey(d, "until") ? _cond(d["until"], p_atm, w, cond_reads, mark_names) : nothing,
                          Float64(get(d, "duration_s", Inf)), Float64(get(d, "timeout_s", Inf)), Float64(get(d, "operator_delay_s", 0.0)), on_off...))
    end
    step_names = [x.name for x in steps]
    # a setting given step by step must name steps that exist, and have a value for each (or a default)
    for q in policies, f in fieldnames(typeof(q))
        v = getfield(q, f)
        v isa PerStep && length(v) > 1 || (v isa PerStep && !haskey(v, "") && !isempty(v)) || continue
        isempty(steps) && error("$where: policy $(q.name).$f is given step by step, but there is no [[controller.step]]")
        for k in keys(v)
            k == "" || k in step_names || error("$where: policy $(q.name).$f names a step $k that does not exist. Steps: $(join(step_names, ", "))")
        end
        haskey(v, "") || all(n -> haskey(v, n), step_names) ||
            error("$where: policy $(q.name).$f has no value for step $(first(n for n in step_names if !haskey(v, n))) (give it one, or a default)")
    end
    end_on[1] == "" || end_on[1] in step_names || error("$where: the run is to end at step $(end_on[1]), which does not exist")
    source = Dict{String,String}()
    for q in policies, r in _reads(q)
        source[r] = _reading_source(r, scene, sim, stand_ins, "$where: policy $(q.name)")
    end
    for r in cond_reads
        source[r] = _reading_source(r, scene, sim, stand_ins, "$where: a step's condition")
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
    signals = isempty(steps) ? String[] : ["ctl.step"]
    append!(signals, ["ctl.mark." * m for m in mark_names])
    for q in policies
        append!(signals, ["ctl.$(q.name).enabled", "ctl.$(q.name).cmd"])
        isempty(_reads(q)) || push!(signals, "ctl.$(q.name).obs")
    end
    Controller(dt, latency, delay, source, stand_ins, noise, rng, policies, "", steps, 0, Dict{String,Float64}(), Dict{String,Float64}(),
               Dict{String,Bool}(), NaN, Dict{String,Float64}(), mark_names, end_on, false, Dict{String,Float64}(), Command[],
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
_command!(c::Controller, sim::Sim, p::Policy, t, part, action, value = 0.0) = _command!(c, sim, p.name, t, part, action, value)
function _command!(c::Controller, sim::Sim, who::String, t, part, action, value = 0.0)
    push!(c.log, Command(t, who, part, action, value))
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

function _tick!(p::BandWatch, c::Controller, sim::Sim, t)
    haskey(c.marks, p.mark) && return
    x = c.obs[p.reading]
    if _at(p.lo, c.step) <= x <= _at(p.hi, c.step)
        isnan(p.since) && (p.since = t)
        t - p.since >= _at(p.hold_s, c.step) - 1e-9 && (c.marks[p.mark] = t)
    else
        p.since = NaN
    end
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

# Enter step k at t (k past the last: the sequence is over, and its last step's settings stay).
function _enter!(c::Controller, sim::Sim, k, t)
    c.k = k
    k > length(c.steps) && return
    st = c.steps[k]
    c.step, c.t_start[st.name], c.due_end = st.name, t, NaN
    st.until === nothing || _reset!(st.until)
    foreach(a -> _command!(c, sim, "sequence", t, a...), st.on_enter)
    for p in c.policies
        p.name in st.enable && (p.enabled = true)
        p.name in st.disable && (p.enabled = false)
    end
end

# The sequence at a tick: start it, or end the running step if its time has come and enter the
# next. One step ends at a tick at most.
function _sequence!(c::Controller, sim::Sim, t)
    isempty(c.steps) && return
    c.k == 0 && return _enter!(c, sim, 1, t)
    c.k > length(c.steps) && return
    st = c.steps[c.k]
    el = t - c.t_start[st.name]
    st.until !== nothing && isnan(c.due_end) && _holds(st.until, c, t) && (c.due_end = t + st.operator_delay)
    ending = el >= st.duration - 1e-9 || (!isnan(c.due_end) && t >= c.due_end - 1e-9)
    timed = !ending && el >= st.timeout - 1e-9
    ending | timed || return
    foreach(a -> _command!(c, sim, "sequence", t, a...), st.on_exit)
    c.t_end[st.name], c.timed_out[st.name] = t, timed
    _enter!(c, sim, c.k + 1, t)
end

# Whether the edge the run is to end at has come.
_at_end(c::Controller) = c.end_on[1] != "" && haskey(c.end_on[2] == "start" ? c.t_start : c.t_end, c.end_on[1])

function _update_signals!(c::Controller)
    k = 0
    if !isempty(c.steps)
        c.values[1] = c.k <= length(c.steps) ? c.k : 0
        k = 1
    end
    for m in c.mark_names
        c.values[k+1] = haskey(c.marks, m) ? 1.0 : 0.0
        k += 1
    end
    for p in c.policies
        c.values[k+1] = p.enabled ? 1.0 : 0.0
        c.values[k+2] = _last_cmd(c, p)
        k += 2
        isempty(_reads(p)) && continue
        c.values[k+1] = c.obs[_reads(p)[1]]
        k += 1
    end
end

# One tick: read, move the sequence on, then each enabled policy in order (so a policy switched on
# by a step that starts at this tick acts at this tick, with that step's settings).
function _tick!(c::Controller, sim::Sim, t)
    _pad!(c, sim, t)
    _observe!(c, sim, t)
    _sequence!(c, sim, t)
    for p in c.policies
        p.enabled ? _tick!(p, c, sim, t) : p isa BandWatch && (p.since = NaN)     # a watch counts only while it is on
    end
    _update_signals!(c)
    _pad!(c, sim)
    c.ended = _at_end(c)
end

"""
    run_controlled!(sim, c, until; wall_limit, work_limit) -> sim

Run to `until` with the controller ticking every `dt` seconds from t = 0, or until the step edge
the controller is to end at, or until the solver fails or a limit is reached (as `run_until!`). At each tick the run is brought to the tick by
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
            c.ended && return sim
            k += 1
        end
        reach(until)
    finally
        _pad!(c, sim)                          # whatever stopped the run, the ctl.* signals cover it
    end
    sim
end
