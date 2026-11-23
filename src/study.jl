# Studies: the same scene run many times with parameters varied, each run scored against a log
# and reduced to a few numbers. (Called a study because `sweep!` is already the network solver's
# pass over the flows.) A study is described by a TOML file; see docs/STUDIES.md.

"One run of a study: which parameters it changes and which extra actions it adds."
struct StudyRun
    id::String                                # readable and unique: the overrides and the case name
    case::String
    overrides::Vector{Pair{String,Float64}}
    actions::Vector{Action}
end

"""
One scored number: a channel over a phase, or over one time bin of a phase. Which log samples it
covers is worked out once, from the log, when the study is read; a run only fills in the value.
"""
struct Scored
    name::String                              # "PT202_psig.vent1", or "PT202_psig.hold.3" for the third bin
    channel::Int
    phase::Int
    idx::Vector{Int}                          # the log samples it covers
    sigma::Float64                            # how far a difference is trusted here, in the log's units:
                                              #   the reading's sigma and the phase's stated model error together
    binned::Bool
end

"A study file, read and checked."
struct Study
    dir::String                               # where results go
    scene::Scene
    schedule::Vector{Action}                  # the replay's schedule (empty without a log)
    t_log::Vector{Float64}                    # log time, seconds after t0 (empty without a log)
    cols::Dict{String,Vector{Float64}}
    columns::Vector{LogColumn}
    until::Float64
    sample_dt::Float64
    reltol::Float64
    wall_limit::Float64
    work_limit::Int                           # give up on a run after this many part-flow evaluations
    runs::Vector{StudyRun}
    varied::Vector{String}
    spans::Vector{Tuple{Float64,Float64,Bool}}   # each varied parameter's lowest and highest value, and whether it is on a ratio scale
    mode::String
    channels::Vector{NamedTuple}              # (log column, model signal, sigma in the log's units, weight, y: the log's readings,
                                              #   offset: "none", or its zero is unknown: "free" (one for the run) or "per_phase")
    phases::Vector{Pair{String,Tuple{Float64,Float64}}}
    scored::Vector{Scored}                    # every scored number, in the order of the residuals
    targets::Vector{NamedTuple}               # features scored beside them (see `_target`)
    metrics::Vector{Pair{String,Dict{String,Any}}}
    bands::NamedTuple                         # (on, signals, percentiles, grid_s): the spread over the runs to write (see `study`)
    controller::Any                           # nothing, or (; spec, fns, allow_julia, stand_ins, escape_hatch, signals): see control.jl
end

# A small generator of our own, so that a seeded study draws the same samples on every machine and
# every Julia version (the standard library's streams are allowed to change between versions).
mutable struct _Lcg
    s::UInt64
end
function _next!(g::_Lcg)                       # SplitMix64, uniform in [0, 1)
    g.s += 0x9e3779b97f4a7c15
    z = g.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    ((z ⊻ (z >> 31)) >> 11) / 2.0^53
end

# The values one varied parameter takes: a list, or { from, to, n } (with log = true for even ratios).
function _values(spec)
    spec isa AbstractVector && return Float64.(spec)
    spec isa Real && return [Float64(spec)]
    a, b, n = Float64(spec["from"]), Float64(spec["to"]), Int(get(spec, "n", 5))
    n == 1 && return [a]
    get(spec, "log", false) ? [a * (b / a)^(k / (n - 1)) for k in 0:n-1] : [a + (b - a) * k / (n - 1) for k in 0:n-1]
end
_span(spec) = (v = _values(spec); (minimum(v), maximum(v), spec isa AbstractDict && get(spec, "log", false)))

const _STUDY_KEYS = ("scene", "replay", "manifest", "t0", "until", "sample_dt", "reltol", "wall_limit", "work_limit", "mode", "n", "seed",
                     "vary", "case", "score", "targets", "metrics", "set", "sets", "bands", "controller")
const _BANDS_KEYS = ("signals", "percentiles", "grid_s")
const _SCORE_KEYS = ("channels", "phases", "bin_s")
const _CHANNEL_SCORE_KEYS = ("log", "model", "sigma", "weight", "offset")
const _PHASE_KEYS = ("from", "to", "bin_s", "discrepancy")
const _METRIC_FNS = ("final", "max", "min", "mean_over", "slope", "time_in_band", "first_crossing")
const _TARGET_KEYS = ("fn", "log", "model", "tol", "from", "to", "phase", "lo", "hi", "value")

"""
    read_study(path) -> Study

Read a study file: the scene, the log to replay and score against (optional), what to vary and
how, the cases, and the scores and metrics to record. Unknown keys, parameters and signals named
wrongly are errors here, before anything runs.
"""
read_study(path; allow_julia_policy = false) =
    _study(TOML.parsefile(path), dirname(abspath(path)), path, splitext(abspath(path))[1] * "_results", _STUDY_KEYS; allow_julia = allow_julia_policy)

# Overrides split into the scene's and the controller's (`controller.<policy>.<setting>[.<step>]`).
_split_overrides(ov) = ([p for p in ov if !startswith(first(p), "controller.")], [p for p in ov if startswith(first(p), "controller.")])

# Write a controller override into a copy of the [controller] table: `controller.dt`, or a
# policy's setting, or its value for one step (`controller.vent.open_above_psig.hold`).
function _controller_override!(spec, key, v)
    f = split(key, '.')[2:end]
    length(f) == 1 && f[1] in ("dt", "latency_s", "actuation_delay_s") && return (spec[f[1]] = v)
    2 <= length(f) <= 3 || error("cannot read \"$key\": expected controller.<policy>.<setting> or controller.<policy>.<setting>.<step>")
    k = findfirst(p -> get(p, "name", "") == f[1], get(spec, "policy", []))
    k === nothing && error("\"$key\": there is no policy named $(f[1])")
    p = spec["policy"][k]
    if length(f) == 2
        get(p, f[2], 0.0) isa AbstractDict && error("\"$key\": $(f[1]).$(f[2]) is given step by step; name the step, as in $key.<step>")
        p[f[2]] = v
    else
        haskey(p, f[2]) && !(p[f[2]] isa AbstractDict) && error("\"$key\": $(f[1]).$(f[2]) is one value for every step; remove the step from the name")
        get!(p, f[2], Dict{String,Any}())[f[3]] = v
    end
end

# A run's controller: the study's table with the run's controller overrides written in.
function _run_controller(st, ctl_ov, scene, sim; seed = "")
    spec = deepcopy(st.controller.spec)
    for (k, v) in ctl_ov
        _controller_override!(spec, k, v)
    end
    _controller(spec, scene, sim; seed, allow_julia = st.controller.allow_julia, fns = st.controller.fns)
end

# A study from its parsed file. `path` is only for messages; a fit builds one of these per dataset.
function _study(d, base, path, dir, allowed; allow_julia = false)
    for k in keys(d)
        k in allowed || error("$path: unknown key \"$k\". Keys: $(join(allowed, ", "))")
    end
    haskey(d, "scene") || error("$path: no `scene`")
    rel(f) = isabspath(f) ? f : joinpath(base, f)
    scene = load_scene(rel(d["scene"]))
    apply_overrides!(scene, ["$k=$v" for (k, v) in get(d, "set", Dict())])
    schedule, t_log, cols, columns, phases = Action[], Float64[], Dict{String,Vector{Float64}}(), LogColumn[], Pair{String,Tuple{Float64,Float64}}[]
    if haskey(d, "manifest")
        log = load_log(rel(d["manifest"]))
        t0 = Float64(get(d, "t0", log.t0))
        cols, columns = log.cols, log.columns
        t_log = cols[log.names[1]] .- t0
        schedule = log_schedule(cols, columns; t0)
        phases = [k => v .- (t0 - log.t0) for (k, v) in log.phases]
    elseif haskey(d, "replay")
        names, cols = read_csv(rel(d["replay"]); min_dt = 0.05)
        columns = guess_log_columns(names, scene)
        t0 = Float64(get(d, "t0", cols[names[1]][1]))
        t_log = cols[names[1]] .- t0
        schedule = log_schedule(cols, columns; t0)
    end
    until = Float64(get(d, "until", isempty(t_log) ? (isempty(scene.schedule) ? 0.0 : maximum(a.t for a in scene.schedule)) + 60 : t_log[end]))
    # what is varied
    vary = collect(get(d, "vary", Dict{String,Any}()))
    sort!(vary; by = first)                    # a TOML table has no order; fix one
    names_v = String[first(v) for v in vary]
    apply_overrides!(deepcopy(scene), ["$k=$(_values(s)[1])" for (k, s) in vary if !startswith(k, "controller.")])      # checks every name now
    spans = [_span(s) for (_, s) in vary]
    mode = get(d, "mode", "grid")
    haskey(d, "sets") && mode != "sets" && error("$path: `sets` is for mode = \"sets\"")
    sets = Vector{Pair{String,Float64}}[]
    prefix = String[]                          # put before a run's id where the values alone would not be unique
    if mode == "grid"
        vals = [_values(s) for (_, s) in vary]
        for combo in Iterators.product(vals...)
            push!(sets, [names_v[k] => combo[k] for k in eachindex(names_v)])
        end
    elseif mode == "oat"                       # one at a time around the scene's own values
        push!(sets, Pair{String,Float64}[])
        for (k, s) in vary, v in _values(s)
            push!(sets, [k => v])
        end
    elseif mode == "lhs"                       # Latin hypercube: each parameter's range cut into n strata, each used once
        n, g = Int(get(d, "n", 20)), _Lcg(UInt64(get(d, "seed", 1)))
        perms = [(p = collect(1:n); for i in n:-1:2; j = 1 + floor(Int, _next!(g) * i); p[i], p[j] = p[j], p[i]; end; p) for _ in vary]
        for i in 1:n
            push!(sets, [begin
                lo, hi, logscale = _span(vary[k][2])
                u = (perms[k][i] - 1 + _next!(g)) / n
                names_v[k] => (logscale ? lo * (hi / lo)^u : lo + (hi - lo) * u)
            end for k in eachindex(vary)])
        end
    elseif mode == "morris"
        # Screening by elementary effects (Morris, 1991). Each of n trajectories starts at a random
        # point on a grid of four levels per parameter (0, 1/3, 2/3 and 1 of its range) and then moves
        # one parameter at a time, in random order, by two thirds of its range: up from the lower two
        # levels, down from the upper two. That is k + 1 runs for k effects, one per parameter, each
        # taken somewhere different in the space.
        n, g = Int(get(d, "n", 20)), _Lcg(UInt64(get(d, "seed", 1)))
        all(v -> (s = _span(last(v)); s[2] > s[1]), vary) || error("$path: a morris study needs a range for every varied parameter, like { from = 0.5, to = 1.0 }")
        at(u) = [(s = _span(vary[k][2]); names_v[k] => (s[3] ? s[1] * (s[2] / s[1])^u[k] : s[1] + (s[2] - s[1]) * u[k])) for k in eachindex(vary)]
        for traj in 1:n
            u = [floor(4 * _next!(g)) / 3 for _ in vary]
            order = collect(eachindex(vary))
            for i in length(order):-1:2
                j = 1 + floor(Int, _next!(g) * i)
                order[i], order[j] = order[j], order[i]
            end
            push!(sets, at(u)); push!(prefix, "m$traj.0;")
            for (step, k) in enumerate(order)
                u[k] += u[k] < 0.5 ? 2 / 3 : -2 / 3
                push!(sets, at(u)); push!(prefix, "m$traj.$step;")
            end
        end
    elseif mode == "sets"
        # The parameter sets come from a file, one row each: those an earlier study could not rule
        # out (see `study_keep`), to be run here on what may be another scene and schedule.
        haskey(d, "sets") || error("$path: mode \"sets\" needs `sets`, a CSV file with a column for each parameter and a row for each set")
        isempty(vary) || error("$path: a sets study takes its parameters from the file; remove [vary]")
        columns_, rows = _table(rel(d["sets"]))
        append!(names_v, columns_)
        missing_ = [k for k in names_v if !startswith(k, "controller.") && try apply_overrides!(deepcopy(scene), ["$k=0"]); false catch; true end]
        isempty(missing_) || error("$path: $(d["sets"]) has columns that are not parameters of this scene: $(join(missing_, ", "))")
        isempty(rows) && error("$path: $(d["sets"]) has no sets")
        for (i, row) in enumerate(rows)
            length(row) == length(names_v) || error("$path: $(d["sets"]): set $i has $(length(row)) values for $(length(names_v)) columns")
            # an empty cell leaves that parameter as the scene has it
            push!(sets, [k => something(tryparse(Float64, cell), NaN) for (k, cell) in zip(names_v, row) if !isempty(cell)])
            all(isfinite ∘ last, sets[end]) || error("$path: $(d["sets"]): set $i has a value that is not a number")
            push!(prefix, "s$i;")              # two rows may hold the same set: the file is taken as it is
        end
        spans = [(v = [last(p) for set in sets for p in set if first(p) == k]; isempty(v) ? (NaN, NaN, false) : (minimum(v), maximum(v), false)) for k in names_v]
    else
        error("$path: mode must be \"grid\", \"oat\", \"lhs\", \"morris\" or \"sets\", not \"$mode\"")
    end
    isempty(prefix) && (prefix = fill("", length(sets)))
    cases = [(String(c["name"]), [Action(Float64(a["t"]), String(a["part"]), String(a["action"]), Float64(get(a, "value", 0.0))) for a in get(c, "actions", [])])
             for c in get(d, "case", [])]
    isempty(cases) && push!(cases, ("", Action[]))
    mode == "morris" && length(cases) > 1 && error("$path: a morris study takes one case at most (screen each case in its own study)")
    for (name, actions) in cases, a in actions
        find_part(scene, a.part) === nothing && error("$path: case $name acts on $(a.part), which is not in the scene")
    end
    label(v) = string(round(v, sigdigits = 6))
    runs = [StudyRun(pre * join(["$k=$(label(v))" for (k, v) in set], ";") * (isempty(name) ? "" : "|" * name), name, set, actions) for (pre, set) in zip(prefix, sets) for (name, actions) in cases]
    allunique(r.id for r in runs) || error("$path: two runs would be identical (a value listed twice, or two cases with one name)")
    controller = nothing
    if haskey(d, "controller")
        spec = d["controller"]
        spec isa AbstractDict || error("$path: [controller] must be a table")
        fns = Dict{String,Any}()
        for p in get(spec, "policy", [])
            p isa AbstractDict && get(p, "kind", "") == "julia" && allow_julia && haskey(p, "file") && (fns[p["file"]] = _julia_policy(rel(p["file"])))
        end
        controller = (; spec, fns, allow_julia, stand_ins = String[], escape_hatch = false, signals = String[])
        # Every run's controller is checked now, with its overrides: a name or a reading wrong in
        # one run should stop the study before any of it runs.
        sim0 = Sim(deepcopy(scene))
        c = nothing
        for ov in unique([last(_split_overrides(r.overrides)) for r in runs])
            c = try
                _run_controller((; controller), ov, scene, sim0)
            catch err
                err isa ErrorException ? error("$path: $(err.msg)") : rethrow()
            end
        end
        controller = (; spec, fns, allow_julia, stand_ins = c.stand_ins, escape_hatch = c.escape_hatch, signals = c.signals)
    end
    # what is scored: listed channels, or every log reading that a manifest maps to a model signal
    score = get(d, "score", Dict{String,Any}())
    for c in get(score, "channels", [])
        haskey(cols, c["log"]) || error("$path: the log has no column \"$(c["log"])\". Columns: $(join(sort!(collect(keys(cols))), ", "))")
    end
    # How far each reading is trusted: said here, else by the manifest, else the default for its kind.
    column(name) = columns[findfirst(c -> c.name == name, columns)]
    trust(name, given) = (s = Float64(given); isnan(s) && (s = column(name).sigma); isnan(s) ? default_sigma(column(name).role) : s)
    for k in keys(score)
        k in _SCORE_KEYS || error("$path: unknown key \"$k\" in [score]. Keys: $(join(_SCORE_KEYS, ", "))")
    end
    channels = NamedTuple[]
    for c in get(score, "channels", [])
        for k in keys(c)
            k in _CHANNEL_SCORE_KEYS || error("$path: unknown key \"$k\" in a [score] channel. Keys: $(join(_CHANNEL_SCORE_KEYS, ", "))")
        end
        log = String(c["log"])
        offset = String(get(c, "offset", "none"))
        offset in ("none", "free", "per_phase") || error("$path: channel $log: offset must be \"free\", \"per_phase\" or \"none\"")
        push!(channels, (; log, model = String(c["model"]), sigma = trust(log, get(c, "sigma", NaN)), weight = Float64(get(c, "weight", 1.0)), y = cols[log], offset))
    end
    isempty(channels) && !haskey(score, "channels") &&
        append!(channels, [(; log = c.name, model = c.target, sigma = trust(c.name, NaN), weight = 1.0, y = cols[c.name], offset = "none") for c in columns if !isempty(c.target)])
    all(c -> c.sigma > 0, channels) || error("$path: a channel's sigma must be positive")
    bin_s = Dict{String,Float64}()             # phase => bin length, where a phase gives its own
    stated = Dict{String,Any}()                # phase => stated model error: one number, or a table of them by log column
    sched = sort(vcat(schedule, scene.schedule); by = a -> a.t)          # the log's actions and the scene's own
    if haskey(score, "phases")
        phases = Pair{String,Tuple{Float64,Float64}}[]
        for (k, p) in score["phases"]
            if p isa AbstractVector
                push!(phases, String(k) => (Float64(p[1]), Float64(p[2])))
                continue
            end
            for key in keys(p)
                key in _PHASE_KEYS || error("$path: unknown key \"$key\" in phase $k. Keys: $(join(_PHASE_KEYS, ", "))")
            end
            f = _phase_time(p["from"], sched)
            push!(phases, String(k) => (f, _phase_time(p["to"], sched; after = f)))
            haskey(p, "bin_s") && (bin_s[String(k)] = Float64(p["bin_s"]))
            haskey(p, "discrepancy") && (stated[String(k)] = p["discrepancy"])
        end
    end
    isempty(phases) && push!(phases, "all" => (0.0, until))
    sort!(phases; by = first)
    # The scored numbers: for each channel and phase, the log samples that have a reading, as one
    # group or cut into bins. An empty group or bin is left out here, so that every run of the
    # study (and every evaluation of a fit) is scored on exactly the same numbers.
    for (phase, e) in stated, k in (e isa AbstractDict ? keys(e) : ())
        any(c -> c.log == k, channels) || error("$path: phase $phase states a discrepancy for \"$k\", which is not a scored log column")
    end
    scored = Scored[]
    for (ci, c) in enumerate(channels), (pi, (phase, (a, b))) in enumerate(phases)
        e = get(stated, phase, 0.0)
        e = Float64(e isa AbstractDict ? get(e, c.log, 0.0) : e)
        e >= 0 || error("$path: phase $phase: a discrepancy cannot be negative")
        sigma = sqrt(c.sigma^2 + e^2)          # independent errors add in quadrature
        idx = [i for (i, t) in enumerate(t_log) if a <= t <= b && t <= until && isfinite(c.y[i])]
        width = get(bin_s, phase, Float64(get(score, "bin_s", 0.0)))
        width >= 0 || error("$path: bin_s must be positive (or 0 for no bins)")
        if width == 0
            isempty(idx) || push!(scored, Scored("$(c.log).$phase", ci, pi, idx, sigma, false))
            continue
        end
        nb = max(1, round(Int, (b - a) / width))             # equal bins, as near to bin_s long as fit the phase
        bin(i) = min(nb, 1 + floor(Int, (t_log[i] - a) / (b - a) * nb))
        for k in 1:nb
            mine = [i for i in idx if bin(i) == k]
            isempty(mine) || push!(scored, Scored("$(c.log).$phase.$k", ci, pi, mine, sigma, true))
        end
    end
    targets = [_target(String(k), v, path, cols, columns, t_log, until, phases, sched) for (k, v) in get(d, "targets", Dict())]
    sort!(targets; by = g -> g.name)
    for (k, v) in get(d, "metrics", Dict())
        v isa AbstractDict || error("$path: metric $k must be a table like { fn = \"final\", signal = \"TANK.m\" }")
    end
    metrics = Pair{String,Dict{String,Any}}[String(k) => Dict{String,Any}(v) for (k, v) in get(d, "metrics", Dict())]
    sort!(metrics; by = first)
    for (name, m) in metrics
        get(m, "fn", "") in _METRIC_FNS || error("$path: metric $name needs fn, one of $(join(_METRIC_FNS, ", "))")
        haskey(m, "signal") || error("$path: metric $name needs a signal")
        m["fn"] == "slope" && haskey(m, "unit") && error("$path: metric $name: a slope is in the signal's SI units per second; `unit` would shift it, not scale it. Remove `unit`")
    end
    b = get(d, "bands", Dict{String,Any}())
    for k in keys(b)
        k in _BANDS_KEYS || error("$path: unknown key \"$k\" in [bands]. Keys: $(join(_BANDS_KEYS, ", "))")
    end
    bands = (; on = haskey(d, "bands"), signals = String.(get(b, "signals", String[])), percentiles = Float64.(get(b, "percentiles", [5, 50, 95])), grid_s = Float64(get(b, "grid_s", 1.0)))
    bands.grid_s > 0 || error("$path: [bands] grid_s must be positive")
    (!isempty(bands.percentiles) && issorted(bands.percentiles) && allunique(bands.percentiles) && all(p -> 0 <= p <= 100, bands.percentiles)) ||
        error("$path: [bands] percentiles must be numbers from 0 to 100 in rising order, like [5, 50, 95]")
    if !isempty(bands.signals)                 # the scene's signals are those a run of it records at its start
        known = vcat(collect(keys(Sim(deepcopy(scene)).hist)), controller === nothing ? String[] : controller.signals)
        unknown = [k for k in bands.signals if !(k in known)]
        isempty(unknown) || error("$path: [bands] names signals this scene does not record: $(join(unknown, ", ")). Signals: $(join(sort!(collect(known)), ", "))")
        allunique(bands.signals) || error("$path: [bands] names a signal twice")
    end
    Study(dir, scene, schedule, t_log, cols, columns, until, Float64(get(d, "sample_dt", 0.5)),
          Float64(get(d, "reltol", 1e-5)), Float64(get(d, "wall_limit", 600.0)), Int(get(d, "work_limit", 5_000_000)), runs, names_v, spans, mode, channels, phases, scored, targets, metrics, bands, controller)
end

# A target: one feature of a log channel (its minimum in a vent, its slope in a hold), to be
# compared with the same feature of the model. The log's side is worked out here, once.
function _target(name, spec, path, cols, columns, t_log, until, phases, sched)
    spec isa AbstractDict || error("$path: target $name must be a table like { fn = \"min\", log = \"PT202_psig\", model = \"TANK.P\", from = 236, to = 310, tol = 5 }")
    haskey(spec, "unit") && error("$path: target $name: a target is always in its log column's units (psi gauge, °C, kg); remove `unit`")
    for k in keys(spec)
        k in _TARGET_KEYS || error("$path: unknown key \"$k\" in target $name. Keys: $(join(_TARGET_KEYS, ", "))")
    end
    for k in ("fn", "log", "model", "tol")
        haskey(spec, k) || error("$path: target $name needs `$k`")
    end
    m = Dict{String,Any}(spec)
    m["fn"] in _METRIC_FNS || error("$path: target $name: fn must be one of $(join(_METRIC_FNS, ", "))")
    haskey(cols, m["log"]) || error("$path: target $name: the log has no column \"$(m["log"])\"")
    tol = Float64(m["tol"])
    tol > 0 || error("$path: target $name: tol must be positive")
    if haskey(m, "phase")                      # the window of a named phase, or from/to as for a phase
        k = findfirst(p -> first(p) == m["phase"], phases)
        k === nothing && error("$path: target $name: no phase called \"$(m["phase"])\"")
        m["from"], m["to"] = last(phases[k])
    else
        haskey(m, "from") && (m["from"] = _phase_time(m["from"], sched))
        haskey(m, "to") && (m["to"] = _phase_time(m["to"], sched; after = get(m, "from", -Inf)))
    end
    y = cols[m["log"]]
    keep = [i for i in eachindex(t_log) if t_log[i] <= until && isfinite(y[i])]     # times with no reading are skipped
    length(keep) >= 2 || error("$path: target $name: the log column \"$(m["log"])\" has no readings")
    value = _metric(t_log[keep], y[keep], m, 0.0)
    isfinite(value) || error("$path: target $name cannot be worked out from the log (no readings in its window, or the value is never crossed)")
    (; name, log = String(m["log"]), model = String(m["model"]), tol, m, t = t_log[keep], value,
       role = columns[findfirst(c -> c.name == m["log"], columns)].role)
end

# The recorded signal a channel or target is scored against. A gauge with no sensor parameter set
# reads the truth and records no reading, so its `P_read` and `T_read` are its `P` and `T`: a
# fit can then start a sensor parameter from its default.
function _scored_signal(sim::Sim, key)
    haskey(sim.hist, key) && return key
    m = match(r"^(.+)\.([PT])_read$", key)
    m !== nothing && get(sim.net.gauges, m[1], 0) != 0 && return m[1] * "." * m[2]
    error("no recorded signal named $key. Signals: $(join(sort!(collect(keys(sim.hist))), ", "))")
end

"A recorded signal at time `t`, by straight-line interpolation of the history."
signal_at(sim::Sim, key, t) = _interp(sim.t_hist, sim.hist[key], t)
function _interp(th, y, t)
    k = clamp(searchsortedlast(th, t), 1, length(th) - 1)
    length(th) == 1 && return y[1]
    y[k] + (y[k+1] - y[k]) * (clamp(t, th[k], th[k+1]) - th[k]) / (th[k+1] - th[k])
end

# A model signal in the units its log column is in (psi gauge, °C, kg).
_as_logged(role, v, p_atm) = role === :pressure ? (v - p_atm) / PSI : role === :temperature ? v - 273.15 : v

"""
    score_run(study, sim) -> (; rms, residuals, z, offsets, targets)

Compare a run with the study's log. The model is interpolated onto the log's times, and times with
no log reading are skipped.

`rms` has one entry per scored channel and phase, in the log's own units (psi, °C, kg): the number
to read.

`z` has one entry per scored number (see `Scored`), with no units: a difference divided by how far
that reading is trusted (`sigma`, widened by any model error stated for the phase), so that a
pressure and a temperature can be weighed against each other. For a channel over a whole phase it
is weight · rms / sigma. For a time bin it is weight · (mean of model − log over the bin) / sigma,
with its sign. The one-number `score` is the root mean square of `z`.

A channel with `offset = "free"` has an unknown zero (a load cell that was not tared). The mean of
model − log over all of that channel's scored samples in this run is taken out before anything is
scored, so only changes count, and is returned in `offsets` (what to add to the log to lay it on
the model). With `offset = "per_phase"` that is done for each phase by itself, for a reading whose
zero also jumps between phases: then only changes inside a phase count.

Each target adds one scored number, (model − log) / tol, where both sides are the same feature
(`fn`) taken over the same window, the model's from its signal at the log's own sample times and
in the log's units. `targets` has the two values. A feature the model does not have at all (a
level it never crosses) counts as 100 tolerances out.

`residuals` is what a fitter minimises. Its sum of squares is the sum of squares of `z`: a bin is
one entry, and a whole phase is one entry per sample, each divided by the square root of the number
of samples, so that a long quiet hold does not outvote a ten-second vent.
"""
function score_run(st::Study, sim::Sim)
    p_atm = st.scene.settings["P_atm_kPa"] * 1e3
    t_end = simtime(sim)
    e2, n = zeros(length(st.channels), length(st.phases)), zeros(Int, length(st.channels), length(st.phases))
    z, residuals, offsets = Pair{String,Float64}[], Float64[], Pair{String,Float64}[]
    for (ci, c) in enumerate(st.channels)
        model = _scored_signal(sim, c.model)
        role = st.columns[findfirst(col -> col.name == c.log, st.columns)].role
        d = fill(NaN, length(st.t_log))                       # model − log at each scored sample the run reached
        for g in st.scored, i in g.idx
            g.channel == ci && st.t_log[i] <= t_end && (d[i] = _as_logged(role, signal_at(sim, model, st.t_log[i]), p_atm) - c.y[i])
        end
        # the mean of d over the samples of these phases (all of them, or one), or nothing to take out
        mean_d(phases) = (v = [d[i] for g in st.scored if g.channel == ci && g.phase in phases for i in g.idx if isfinite(d[i])]; isempty(v) ? NaN : sum(v) / length(v))
        off = c.offset == "free" ? fill(mean_d(eachindex(st.phases)), length(st.phases)) :
              c.offset == "per_phase" ? [mean_d(k) for k in eachindex(st.phases)] : zeros(length(st.phases))
        c.offset == "none" || append!(offsets, last.(_offset_names(st, ci)) .=> (c.offset == "free" ? off[1:1] : off))
        for g in st.scored
            g.channel == ci || continue
            ds = [d[i] - off[g.phase] for i in g.idx if isfinite(d[i])]
            if isempty(ds)                                     # the run stopped before this one
                push!(z, "res." * g.name => NaN)
                continue
            end
            e2[ci, g.phase] += sum(abs2, ds)
            n[ci, g.phase] += length(ds)
            if g.binned
                push!(residuals, c.weight * sum(ds) / length(ds) / g.sigma)
                push!(z, "res." * g.name => residuals[end])
            else
                append!(residuals, c.weight .* ds ./ (g.sigma * sqrt(length(ds))))
                push!(z, "res." * g.name => c.weight * sqrt(sum(abs2, ds) / length(ds)) / g.sigma)
            end
        end
    end
    rms = ["rms.$(c.log).$phase" => (n[ci, pi] == 0 ? NaN : sqrt(e2[ci, pi] / n[ci, pi])) for (ci, c) in enumerate(st.channels) for (pi, (phase, _)) in enumerate(st.phases)]
    targets = Pair{String,Float64}[]
    for g in st.targets
        model = _scored_signal(sim, g.model)
        v = r = NaN
        if t_end >= st.until - 1e-6                            # a run that stopped early has no features to compare
            v = _metric(g.t, [_as_logged(g.role, signal_at(sim, model, t), p_atm) for t in g.t], g.m, 0.0)
            r = isnan(v) ? 100.0 : (v - g.value) / g.tol
            push!(residuals, r)
        end
        push!(z, "res.target.$(g.name)" => r)
        append!(targets, ["target.$(g.name).model" => v, "target.$(g.name).log" => g.value])
    end
    (; rms, residuals, z, offsets, targets)
end

"""
The one-number score of a run: the root mean square of its scored numbers `z`, each a difference
in units of that reading's `sigma`. About 1 means the model is, on average, as far from the log as
the readings themselves can be trusted; 3 means three times that.
"""
function study_score(s)
    v = [last(p) for p in s.z if !isnan(last(p))]            # a run that stopped early is scored on what it reached
    isempty(v) ? NaN : sqrt(sum(abs2, v) / length(v))
end

"How many numbers a complete run is scored on, and how long its residual vector is."
n_scored(st::Study) = length(st.scored) + length(st.targets)
n_residuals(st::Study) = sum(g -> g.binned ? 1 : length(g.idx), st.scored; init = 0) + length(st.targets)

# A metric's argument in SI: `lo_psig = 600` is a gauge pressure, `value_C = 10` a temperature.
function _si(m, key, p_atm)
    haskey(m, key) && return Float64(m[key])
    haskey(m, key * "_psig") && return Float64(m[key*"_psig"]) * PSI + p_atm
    haskey(m, key * "_C") && return Float64(m[key*"_C"]) + 273.15
    error("metric needs \"$key\" (or $(key)_psig, $(key)_C)")
end

"""
One number from a run: the last value of a signal, its extremes, a mean, a slope (per second, by
least squares), the time spent in a band, or when it first crosses a value. With `from` and `to`
(seconds) only that window of the run is looked at. Signals are in SI; `unit = "psig"` or `"C"`
converts a result that is a pressure or a temperature.
"""
function metric(sim::Sim, m, p_atm)
    haskey(sim.hist, m["signal"]) || error("no recorded signal named $(m["signal"])")
    v = _metric(sim.t_hist, sim.hist[m["signal"]], m, p_atm)
    unit = get(m, "unit", "")
    unit == "psig" ? (v - p_atm) / PSI : unit == "C" ? v - 273.15 : v
end
# The same for any series of times and values: a run's history, or a log channel.
# A run is recorded where the solver's steps fall, which can be seconds apart in a quiet stretch.
# So everything here treats a series as the straight lines between samples, not as the samples alone.
function _metric(t, y, m, p_atm)
    fn = m["fn"]
    if haskey(m, "from") || haskey(m, "to")    # only this window, with a point put at each end of it
        a, b = max(Float64(get(m, "from", t[1])), t[1]), min(Float64(get(m, "to", t[end])), t[end])
        b > a || return NaN
        inside = findall(x -> a < x < b, t)
        t, y = vcat(a, t[inside], b), vcat(_interp(t, y, a), y[inside], _interp(t, y, b))
    end
    fn == "final" && return y[end]
    fn == "max" && return maximum(y)
    fn == "min" && return minimum(y)
    length(t) >= 2 || return NaN
    if fn == "mean_over"
        return sum((t[i] - t[i-1]) * (y[i] + y[i-1]) / 2 for i in 2:length(t)) / (t[end] - t[1])
    elseif fn == "slope"
        # Least squares, each sample weighted by the time it stands for, so that a stretch where
        # samples are crowded together counts for its length and not for their number.
        w = [(t[min(i + 1, end)] - t[max(i - 1, 1)]) / 2 for i in eachindex(t)]
        tm, ym = sum(w .* t) / sum(w), sum(w .* y) / sum(w)
        return sum(w .* (t .- tm) .* (y .- ym)) / sum(w .* (t .- tm) .^ 2)
    elseif fn == "time_in_band"
        lo, hi = _si(m, "lo", p_atm), _si(m, "hi", p_atm)
        total = 0.0
        for i in 2:length(t)
            y0, y1 = y[i-1], y[i]
            if y0 == y1
                lo <= y0 <= hi && (total += t[i] - t[i-1])
            else                                # the part of this line that lies between lo and hi
                f1, f2 = minmax((lo - y0) / (y1 - y0), (hi - y0) / (y1 - y0))
                total += (t[i] - t[i-1]) * max(min(f2, 1.0) - max(f1, 0.0), 0.0)
            end
        end
        return total
    else                                       # first_crossing
        v = _si(m, "value", p_atm)
        i = findfirst(i -> (y[i-1] - v) * (y[i] - v) <= 0 && y[i-1] != y[i], 2:length(y))
        return i === nothing ? NaN : t[i] + (t[i+1] - t[i]) * (v - y[i]) / (y[i+1] - y[i])
    end
end

"""
    simulate(study, overrides, actions = []; reltol) -> (sim, failed)

The study's scene with some parameters changed, replaying its log, run to the end. `failed` says
why if it could not be built (then `sim` is nothing), failed, or ran out of wall time.
"""
function simulate(st::Study, overrides, actions = Action[]; reltol = st.reltol)
    st.controller === nothing || return _simulate(st, overrides, actions; reltol)[1:2]
    try
        scene = deepcopy(st.scene)
        apply_overrides!(scene, ["$k=$v" for (k, v) in overrides])
        scene.schedule = sort(vcat(scene.schedule, st.schedule, actions); by = a -> a.t)
        sim = Sim(scene; sample_dt = st.sample_dt, reltol)
        run_until!(sim, st.until; wall_limit = st.wall_limit, work_limit = st.work_limit)
        return sim, !isempty(sim.failed) ? sim.failed : simtime(sim) >= st.until - 1e-6 ? "" :
                    work(sim) >= st.work_limit ? "stopped at the work limit ($(st.work_limit) part-flow evaluations)" : "stopped at the wall-time limit"
    catch err
        return nothing, first(replace(sprint(showerror, err), '\n' => ' ', ',' => ';'), 200)
    end
end

# The same with a controller: (sim, failed, controller). `case` is part of what its noise is seeded from.
function _simulate(st::Study, overrides, actions = Action[]; reltol = st.reltol, case = "")
    c = nothing
    try
        scene_ov, ctl_ov = _split_overrides(overrides)
        scene = deepcopy(st.scene)
        apply_overrides!(scene, ["$k=$v" for (k, v) in scene_ov])
        scene.schedule = sort(vcat(scene.schedule, st.schedule, actions); by = a -> a.t)
        sim = Sim(scene; sample_dt = st.sample_dt, reltol)
        c = _run_controller(st, ctl_ov, scene, sim; seed = _noise_seed(scene_ov, case))
        run_controlled!(sim, c, st.until; wall_limit = st.wall_limit, work_limit = st.work_limit)
        return sim, !isempty(sim.failed) ? sim.failed : simtime(sim) >= st.until - 1e-6 ? "" :
                    work(sim) >= st.work_limit ? "stopped at the work limit ($(st.work_limit) part-flow evaluations)" : "stopped at the wall-time limit", c
    catch err
        return nothing, first(replace(sprint(showerror, err), '\n' => ' ', ',' => ';'), 200), c
    end
end

"Run one member of a study and reduce it to a row of results (column name => value)."
function run_one(st::Study, run::StudyRun; keep = nothing, seen = nothing)
    row = Pair{String,Any}["id" => run.id, "case" => run.case]
    append!(row, [k => get(Dict(run.overrides), k, NaN) for k in st.varied])
    started = Base.time()
    sim, failed = st.controller === nothing ? simulate(st, run.overrides, run.actions) : _simulate(st, run.overrides, run.actions; case = run.case)
    if sim === nothing
        append!(row, [k => NaN for k in vcat("score", first.(_rms_names(st)), ["res." * g.name for g in st.scored], ["res.target." * g.name for g in st.targets],
                                         first.(_offset_names(st)), first.(_target_names(st)), first.(st.metrics))])
        append!(row, ["t_end" => NaN, "failed" => failed, "unconverged" => 0, "mass_error" => NaN, "energy_error" => NaN, "vented_kg" => NaN])
    else
        s = score_run(st, sim)
        push!(row, "score" => study_score(s))
        append!(row, s.rms)
        append!(row, s.z)
        append!(row, s.offsets)
        append!(row, s.targets)
        p_atm = st.scene.settings["P_atm_kPa"] * 1e3
        append!(row, [name => metric(sim, m, p_atm) for (name, m) in st.metrics])
        r = run_summary(sim)
        append!(row, ["t_end" => r.t_s, "failed" => replace(failed, ',' => ';'), "unconverged" => r.unconverged, "mass_error" => r.mass_error,
                      "energy_error" => r.energy_error, "vented_kg" => r.vented_kg])
        seen === nothing || seen(sim)                          # the caller takes what it wants from the run before it is dropped
        keep === nothing || write_csv(keep, sim; columns = unique(vcat([_scored_signal(sim, c.model) for c in st.channels], [_scored_signal(sim, g.model) for g in st.targets], [m["signal"] for (_, m) in st.metrics])))
    end
    st.controller === nothing || append!(row, ["stand_ins" => join(st.controller.stand_ins, ";"), "escape_hatch" => st.controller.escape_hatch ? 1 : 0])
    push!(row, "wall_s" => round(Base.time() - started, digits = 3))
    row
end
_target_names(st::Study) = ["target.$(g.name).$side" => NaN for g in st.targets for side in ("model", "log")]
# The offset columns: one for a channel whose zero is free, one per phase where it is free in each.
_offset_names(st::Study, ci) = (c = st.channels[ci]; c.offset == "free" ? [NaN => "offset.$(c.log)"] : c.offset == "per_phase" ? [NaN => "offset.$(c.log).$p" for (p, _) in st.phases] : Pair{Float64,String}[])
_offset_names(st::Study) = [name => NaN for ci in eachindex(st.channels) for (_, name) in _offset_names(st, ci)]
_rms_names(st::Study) = ["rms.$(c.log).$p" => NaN for c in st.channels for (p, _) in st.phases]

"""
    study(path; resume = false, keep_traces = false, allow_julia_policy = false) -> results file

Run every member of the study in `path`, on as many threads as Julia was started with, and append
one row per run to `<study>_results/results.csv` as it finishes. A run that fails is a row with
the reason in `failed`. With `resume`, runs whose id is already in the file are skipped.

With `[bands]` in the file, also write `bands.csv`: for each time on a grid, the percentiles of
each named signal over the runs (see `percentile_bands`), and `metric_bands.csv`, the same for each
metric. With several cases there is a `bands_<case>.csv` for each: cases are different scenarios,
and are not pooled. Each run is reduced to its signals on the grid as it finishes, so what is held
is runs × grid times × signals numbers and never a whole run.
"""
function study(path; resume = false, keep_traces = false, quiet = false, allow_julia_policy = false)
    st = read_study(path; allow_julia_policy)
    resume && !isempty(st.bands.signals) && error("$path: a study with [bands] signals cannot be resumed: the runs already done left no history to take a band from. Run it without --resume")
    mkpath(st.dir)
    out = joinpath(st.dir, "results.csv")
    done = Set{String}()
    if resume && isfile(out)
        for line in Iterators.drop(eachline(out), 1)
            push!(done, String(first(split(line, ','))))
        end
    else
        rm(out; force = true)
    end
    todo = [r for r in st.runs if !(r.id in done)]
    quiet || println(length(st.runs), " runs", isempty(done) ? "" : " ($(length(done)) already done)", " on ", Threads.nthreads(), " thread(s); results in ", out)
    keep_traces && mkpath(joinpath(st.dir, "traces"))
    lock_ = ReentrantLock()
    header = !isfile(out)
    grids = [zeros(0, length(st.bands.signals)) for _ in st.runs]          # a run that cannot be built stays empty
    seen = isempty(st.bands.signals) ? nothing : (index, sim) -> (grids[index] = _gridded(st, sim))
    Threads.@threads :dynamic for k in eachindex(todo)
        run = todo[k]
        index = findfirst(r -> r.id == run.id, st.runs)
        row = run_one(st, run; keep = keep_traces ? joinpath(st.dir, "traces", "$index.csv") : nothing, seen = seen === nothing ? nothing : sim -> seen(index, sim))
        lock(lock_) do
            open(out, "a") do io
                header && (println(io, join(first.(row), ",")); header = false)
                println(io, join((v isa AbstractFloat ? (isnan(v) ? "" : string(round(v, sigdigits = 7))) : string(v) for v in last.(row)), ","))
            end
        end
    end
    st.mode == "morris" && _screening(st, out, quiet)
    st.bands.on && _write_bands(st, out, grids, quiet)
    out
end

_band_grid(st::Study) = 0:st.bands.grid_s:st.until
# A run's band signals on the study's time grid: one row per grid time the run reached. (A history
# read beyond its end would give its last value for ever, so a run that stopped is cut off here.)
function _gridded(st::Study, sim::Sim)
    grid = _band_grid(st)
    n = count(<=(simtime(sim) + 1e-6), grid)
    [signal_at(sim, key, grid[i]) for i in 1:n, key in st.bands.signals]
end

# The p-th percentile of sorted values, by straight-line interpolation between the two nearest.
function _percentile(v, p)
    h = 1 + (length(v) - 1) * p / 100
    k = floor(Int, h)
    v[k] + (h - k) * (v[min(k + 1, end)] - v[k])
end

"""
    percentile_bands(runs, percentiles) -> (n, bands)

The spread of some signals over several runs. `runs` has one matrix per run, with a row per grid
time and a column per signal; a run that stopped early has fewer rows. `n[i]` is how many runs
reached grid time `i`, and `bands[i, signal, k]` the `percentiles[k]`-th percentile of the signal
over those runs. There is a row for every grid time that at least one run reached.
"""
function percentile_bands(runs, percentiles)
    rows, signals = maximum(r -> size(r, 1), runs; init = 0), isempty(runs) ? 0 : size(runs[1], 2)
    n = [count(r -> size(r, 1) >= i, runs) for i in 1:rows]
    bands = zeros(rows, signals, length(percentiles))
    for i in 1:rows, s in 1:signals
        v = sort!([r[i, s] for r in runs if size(r, 1) >= i])
        bands[i, s, :] = [_percentile(v, p) for p in percentiles]
    end
    n, bands
end

_percent_label(p) = "p" * (isinteger(p) ? string(Int(p)) : string(p))

# After a study with [bands]: bands.csv from the runs' gridded signals, and metric_bands.csv from
# the results file (read back, as for a screening). Each case has its own.
function _write_bands(st::Study, out, grids, quiet)
    cases = unique(r.case for r in st.runs)
    ps = st.bands.percentiles
    cell(x) = isnan(x) ? "" : string(round(x, sigdigits = 7))
    for case in (isempty(st.bands.signals) ? () : cases)
        n, bands = percentile_bands([g for (g, r) in zip(grids, st.runs) if r.case == case], ps)
        file = joinpath(st.dir, length(cases) == 1 ? "bands.csv" : "bands_$case.csv")
        open(file, "w") do io
            println(io, "t_s,n,", join(["$key.$(_percent_label(p))" for key in st.bands.signals for p in ps], ","))
            for i in eachindex(n)
                println(io, _band_grid(st)[i], ",", n[i], ",", join([cell(bands[i, s, k]) for s in eachindex(st.bands.signals) for k in eachindex(ps)], ","))
            end
        end
        quiet || println("bands in ", file)
    end
    isempty(st.metrics) && return
    head, rows = _table(out)
    failed, case_col = findfirst(==("failed"), head), findfirst(==("case"), head)
    file = joinpath(st.dir, "metric_bands.csv")
    open(file, "w") do io
        println(io, "case,metric,n,", join(_percent_label.(ps), ","))
        for case in cases, (name, _) in st.metrics
            c = findfirst(==(name), head)
            # a run that stopped partway has numbers, but not ones to put beside the others
            v = sort!(filter(isfinite, [something(tryparse(Float64, row[c]), NaN) for row in rows if row[case_col] == case && isempty(row[failed])]))
            println(io, case, ",", name, ",", length(v), ",", join([isempty(v) ? "" : cell(_percentile(v, p)) for p in ps], ","))
        end
    end
    quiet || println("metric bands in ", file)
    nothing
end

"""
    elementary_effects(X, y) -> [(mu_star, sigma, n) for each parameter]

Reduce a Morris design to its two numbers per parameter. `X` has one row per run, in trajectory
order (k + 1 rows per trajectory, each differing from the one before in a single column), with each
parameter as its position in its range, 0 to 1. `y` is one output of those runs. An elementary
effect is the change in `y` over one step divided by the step, so it is per whole range of the
parameter and parameters can be compared.

`mu_star` is the mean of the effects' absolute values: how much the parameter matters. `sigma` is
their spread: large beside `mu_star` means the effect depends on where the other parameters are
(an interaction) or is far from a straight line. `n` is how many effects there were; a step with a
failed run at either end is left out.
"""
function elementary_effects(X, y)
    k = size(X, 2)
    effects = [Float64[] for _ in 1:k]
    for start in 1:k+1:size(X, 1), a in start:start+k-1
        j = findfirst(j -> X[a, j] != X[a+1, j], 1:k)
        e = (y[a+1] - y[a]) / (X[a+1, j] - X[a, j])
        isfinite(e) && push!(effects[j], e)
    end
    map(effects) do e
        n = length(e)
        m = n == 0 ? NaN : sum(e) / n
        (mu_star = n == 0 ? NaN : sum(abs, e) / n, sigma = n < 2 ? NaN : sqrt(sum(abs2, e .- m) / (n - 1)), n = n)
    end
end

# After a morris study: read the results back (so that resumed runs count too), and write
# screening.csv with a row for every varied parameter and every scored or measured number.
function _screening(st::Study, out, quiet)
    lines = readlines(out)
    head = split(lines[1], ',')
    rows = Dict(String(first(f)) => f for f in split.(lines[2:end], ','))
    failed = findfirst(==("failed"), head)                      # a run that stopped partway has numbers, but not ones to compare
    X = [(v = Dict(r.overrides)[name]; (lo, hi, logscale) = st.spans[k]; logscale ? log(v / lo) / log(hi / lo) : (v - lo) / (hi - lo)) for r in st.runs, (k, name) in enumerate(st.varied)]
    outputs = [h for h in head if (h == "score" && n_scored(st) > 0) || startswith(h, "res.") || (startswith(h, "target.") && endswith(h, ".model")) || h in first.(st.metrics)]
    file = joinpath(st.dir, "screening.csv")
    open(file, "w") do io
        println(io, "parameter,output,mu_star,sigma,n")
        for o in outputs
            c = findfirst(==(o), head)
            y = [haskey(rows, r.id) && isempty(rows[r.id][failed]) ? something(tryparse(Float64, rows[r.id][c]), NaN) : NaN for r in st.runs]
            ee = elementary_effects(X, y)
            for (k, e) in enumerate(ee)
                println(io, st.varied[k], ",", o, ",", round(e.mu_star, sigdigits = 6), ",", round(e.sigma, sigdigits = 6), ",", e.n)
            end
            if o == "score" && !quiet
                println("Screening: mean |effect| on the score of moving each parameter across its whole range (and the spread of the effects)")
                for k in sortperm([e.mu_star for e in ee]; rev = true)
                    println("  ", rpad(st.varied[k], 34), rpad(round(ee[k].mu_star, sigdigits = 4), 12), "(", round(ee[k].sigma, sigdigits = 4), ", from ", ee[k].n, " effects)")
                end
            end
        end
    end
    quiet || println("screening in ", file)
    file
end

# A results file, or a file of parameter sets, as its header and its rows of cells. (Not `read_csv`:
# a run's id is not a number, and the order of the columns matters here.)
function _table(path)
    lines = filter(!isempty, readlines(path))
    isempty(lines) && error("$path is empty")
    String.(strip.(split(lines[1], ','))), [String.(strip.(split(l, ','))) for l in lines[2:end]]
end

"""
    keep_runs(results; cutoff = 3, nth = 2) -> (; names, sets, cells, runs, kept, scored, ranges)

The runs of a study that its log could not rule out. Every `res.*` column of a `results.csv` is a
difference between model and log in units of how far that reading is trusted (its sigma and any
stated model error together). A run is ruled out when its `nth` largest |res| is over `cutoff`: with
the defaults, when two or more of its scored numbers are more than three of those units out. One is
allowed, because with many scored numbers a single one will stray by chance; allow more (a larger
`nth`) where a study has many. A run that failed or stopped early is ruled out whatever it scored.

`sets` has the kept runs' parameter values (name => value; a parameter the run left alone is not
in its set) and `cells` the same as they stand in the file. `ranges` has, for each varied parameter,
its lowest and highest value among the kept runs and among all of them (`kept_lo`, `kept_hi`, `lo`,
`hi`): a parameter whose kept range is as wide as the sampled one was not determined by this log.
"""
function keep_runs(path; cutoff = 3.0, nth = 2)
    head, rows = _table(path)
    a, b, failed = findfirst(==("case"), head), findfirst(==("score"), head), findfirst(==("failed"), head)
    (a === nothing || b === nothing || failed === nothing) && error("$path is not a study's results.csv (it has no `case`, `score` and `failed` columns)")
    res = findall(startswith("res."), head)
    isempty(res) && error("$path has no res.* columns: its study was not scored against a log, so no run can be ruled out")
    1 <= nth <= length(res) || error("nth = $nth, but $path has $(length(res)) scored numbers: nth must be between 1 and that")
    num(cell) = something(tryparse(Float64, cell), NaN)
    filter!(row -> length(row) == length(head), rows)          # a row cut short by a study that was stopped
    kept = filter(rows) do row
        z = sort!([abs(num(row[c])) for c in res]; rev = true)
        isempty(row[failed]) && !any(isnan, z) && z[nth] <= cutoff
    end
    names = head[a+1:b-1]
    span(rows, c) = (v = filter(isfinite, [num(row[c]) for row in rows]); isempty(v) ? (NaN, NaN) : extrema(v))
    ranges = [(; name, kept_lo = span(kept, a + k)[1], kept_hi = span(kept, a + k)[2], lo = span(rows, a + k)[1], hi = span(rows, a + k)[2]) for (k, name) in enumerate(names)]
    sets = [[name => num(row[a+k]) for (k, name) in enumerate(names) if !isempty(row[a+k])] for row in kept]
    (; names, sets, cells = [row[a+1:b-1] for row in kept], runs = length(rows), kept = length(kept), scored = length(res), ranges)
end

"""
    study_keep(results; cutoff = 3, nth = 2) -> kept file

Write the parameter sets of the runs `keep_runs` keeps to `kept.csv` beside the results (one row per
kept run, one column per varied parameter, as a study with `mode = "sets"` reads them) and print how
many were kept and what range of each parameter is left.
"""
function study_keep(path; cutoff = 3.0, nth = 2, quiet = false)
    k = keep_runs(path; cutoff, nth)
    out = joinpath(dirname(abspath(path)), "kept.csv")
    write(out, join(vcat(join(k.names, ","), join.(k.cells, ",")), "\n") * "\n")
    if !quiet
        println(k.runs, " runs, ", k.kept, " kept (ruled out: failed, or |res| over ", cutoff, " in ", nth, " or more of the ", k.scored, " scored numbers)")
        println("  ", rpad("parameter", 34), rpad("kept", 26), "sampled")
        r4(x) = round(x, sigdigits = 4)
        for r in k.ranges
            println("  ", rpad(r.name, 34), rpad(k.kept == 0 ? "-" : "$(r4(r.kept_lo)) to $(r4(r.kept_hi))", 26), r4(r.lo), " to ", r4(r.hi))
        end
        println("kept sets in ", out)
    end
    out
end
