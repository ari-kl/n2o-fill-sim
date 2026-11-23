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
    channels::Vector{NamedTuple}              # (log column, model signal, sigma in the log's units, weight, y: the log's readings)
    phases::Vector{Pair{String,Tuple{Float64,Float64}}}
    scored::Vector{Scored}                    # every scored number, in the order of the residuals
    metrics::Vector{Pair{String,Dict{String,Any}}}
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
                     "vary", "case", "score", "metrics", "set")
const _SCORE_KEYS = ("channels", "phases", "bin_s")
const _CHANNEL_SCORE_KEYS = ("log", "model", "sigma", "weight")
const _PHASE_KEYS = ("from", "to", "bin_s", "discrepancy")
const _METRIC_FNS = ("final", "max", "min", "mean_over", "time_in_band", "first_crossing")

"""
    read_study(path) -> Study

Read a study file: the scene, the log to replay and score against (optional), what to vary and
how, the cases, and the scores and metrics to record. Unknown keys, parameters and signals named
wrongly are errors here, before anything runs.
"""
read_study(path) = _study(TOML.parsefile(path), dirname(abspath(path)), path, splitext(abspath(path))[1] * "_results", _STUDY_KEYS)

# A study from its parsed file. `path` is only for messages; a fit builds one of these per dataset.
function _study(d, base, path, dir, allowed)
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
    apply_overrides!(deepcopy(scene), ["$k=$(_values(s)[1])" for (k, s) in vary])      # checks every name now
    mode = get(d, "mode", "grid")
    sets = Vector{Pair{String,Float64}}[]
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
    else
        error("$path: mode must be \"grid\", \"oat\" or \"lhs\", not \"$mode\"")
    end
    cases = [(String(c["name"]), [Action(Float64(a["t"]), String(a["part"]), String(a["action"]), Float64(get(a, "value", 0.0))) for a in get(c, "actions", [])])
             for c in get(d, "case", [])]
    isempty(cases) && push!(cases, ("", Action[]))
    for (name, actions) in cases, a in actions
        find_part(scene, a.part) === nothing && error("$path: case $name acts on $(a.part), which is not in the scene")
    end
    label(v) = string(round(v, sigdigits = 6))
    runs = [StudyRun(join(["$k=$(label(v))" for (k, v) in set], ";") * (isempty(name) ? "" : "|" * name), name, set, actions) for set in sets for (name, actions) in cases]
    allunique(r.id for r in runs) || error("$path: two runs would be identical (a value listed twice, or two cases with one name)")
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
        push!(channels, (; log, model = String(c["model"]), sigma = trust(log, get(c, "sigma", NaN)), weight = Float64(get(c, "weight", 1.0)), y = cols[log]))
    end
    isempty(channels) && !haskey(score, "channels") &&
        append!(channels, [(; log = c.name, model = c.target, sigma = trust(c.name, NaN), weight = 1.0, y = cols[c.name]) for c in columns if !isempty(c.target)])
    all(c -> c.sigma > 0, channels) || error("$path: a channel's sigma must be positive")
    bin_s = Dict{String,Float64}()             # phase => bin length, where a phase gives its own
    stated = Dict{String,Any}()                # phase => stated model error: one number, or a table of them by log column
    if haskey(score, "phases")
        sched = sort(vcat(schedule, scene.schedule); by = a -> a.t)
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
    for (k, v) in get(d, "metrics", Dict())
        v isa AbstractDict || error("$path: metric $k must be a table like { fn = \"final\", signal = \"TANK.m\" }")
    end
    metrics = Pair{String,Dict{String,Any}}[String(k) => Dict{String,Any}(v) for (k, v) in get(d, "metrics", Dict())]
    sort!(metrics; by = first)
    for (name, m) in metrics
        get(m, "fn", "") in _METRIC_FNS || error("$path: metric $name needs fn, one of $(join(_METRIC_FNS, ", "))")
        haskey(m, "signal") || error("$path: metric $name needs a signal")
    end
    Study(dir, scene, schedule, t_log, cols, columns, until, Float64(get(d, "sample_dt", 0.5)),
          Float64(get(d, "reltol", 1e-5)), Float64(get(d, "wall_limit", 600.0)), Int(get(d, "work_limit", 5_000_000)), runs, names_v, channels, phases, scored, metrics)
end

"A recorded signal at time `t`, by straight-line interpolation of the history."
function signal_at(sim::Sim, key, t)
    th, y = sim.t_hist, sim.hist[key]
    k = clamp(searchsortedlast(th, t), 1, length(th) - 1)
    length(th) == 1 && return y[1]
    y[k] + (y[k+1] - y[k]) * (clamp(t, th[k], th[k+1]) - th[k]) / (th[k+1] - th[k])
end

# A model signal in the units its log column is in (psi gauge, °C, kg).
_as_logged(role, v, p_atm) = role === :pressure ? (v - p_atm) / PSI : role === :temperature ? v - 273.15 : v

"""
    score_run(study, sim) -> (; rms, residuals, z)

Compare a run with the study's log. The model is interpolated onto the log's times, and times with
no log reading are skipped.

`rms` has one entry per scored channel and phase, in the log's own units (psi, °C, kg): the number
to read.

`z` has one entry per scored number (see `Scored`), with no units: a difference divided by how far
that reading is trusted (`sigma`, widened by any model error stated for the phase), so that a
pressure and a temperature can be weighed against each other. For a channel over a whole phase it is weight · rms / sigma. For a time bin it is
weight · (mean of model − log over the bin) / sigma, with its sign. The one-number `score` is the
root mean square of `z`.

`residuals` is what a fitter minimises. Its sum of squares is the sum of squares of `z`: a bin is
one entry, and a whole phase is one entry per sample, each divided by the square root of the number
of samples, so that a long quiet hold does not outvote a ten-second vent.
"""
function score_run(st::Study, sim::Sim)
    p_atm = st.scene.settings["P_atm_kPa"] * 1e3
    t_end = simtime(sim)
    e2, n = zeros(length(st.channels), length(st.phases)), zeros(Int, length(st.channels), length(st.phases))
    z, residuals = Pair{String,Float64}[], Float64[]
    for (ci, c) in enumerate(st.channels)
        haskey(sim.hist, c.model) || error("no recorded signal named $(c.model). Signals: $(join(sort!(collect(keys(sim.hist))), ", "))")
        role = st.columns[findfirst(col -> col.name == c.log, st.columns)].role
        d = fill(NaN, length(st.t_log))                       # model − log at each scored sample the run reached
        for g in st.scored, i in g.idx
            g.channel == ci && st.t_log[i] <= t_end && (d[i] = _as_logged(role, signal_at(sim, c.model, st.t_log[i]), p_atm) - c.y[i])
        end
        for g in st.scored
            g.channel == ci || continue
            ds = [d[i] for i in g.idx if isfinite(d[i])]
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
    (; rms, residuals, z)
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
n_scored(st::Study) = length(st.scored)
n_residuals(st::Study) = sum(g -> g.binned ? 1 : length(g.idx), st.scored; init = 0)

# A metric's argument in SI: `lo_psig = 600` is a gauge pressure, `value_C = 10` a temperature.
function _si(m, key, p_atm)
    haskey(m, key) && return Float64(m[key])
    haskey(m, key * "_psig") && return Float64(m[key*"_psig"]) * PSI + p_atm
    haskey(m, key * "_C") && return Float64(m[key*"_C"]) + 273.15
    error("metric needs \"$key\" (or $(key)_psig, $(key)_C)")
end

"""
One number from a run: the last value of a signal, its extremes, a mean, the time spent in a band,
or when it first crosses a value. Signals are in SI; `unit = "psig"` or `"C"` converts a result
that is a pressure or a temperature.
"""
function metric(sim::Sim, m, p_atm)
    haskey(sim.hist, m["signal"]) || error("no recorded signal named $(m["signal"])")
    v = _metric(sim, m, p_atm)
    unit = get(m, "unit", "")
    unit == "psig" ? (v - p_atm) / PSI : unit == "C" ? v - 273.15 : v
end
function _metric(sim::Sim, m, p_atm)
    y, t = sim.hist[m["signal"]], sim.t_hist
    fn = m["fn"]
    fn == "final" && return y[end]
    fn == "max" && return maximum(y)
    fn == "min" && return minimum(y)
    # The history is recorded where the solver's steps fall, which can be seconds apart in a quiet
    # stretch. So these treat it as the straight lines between samples, not as the samples alone.
    if fn == "mean_over"
        a, b = Float64(m["from"]), min(Float64(m["to"]), t[end])
        b > a || return NaN
        knots = vcat(a, [x for x in t if a < x < b], b)
        vals = [signal_at(sim, m["signal"], x) for x in knots]
        return sum((knots[i] - knots[i-1]) * (vals[i] + vals[i-1]) / 2 for i in 2:length(knots)) / (b - a)
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

"Run one member of a study and reduce it to a row of results (column name => value)."
function run_one(st::Study, run::StudyRun; keep = nothing)
    row = Pair{String,Any}["id" => run.id, "case" => run.case]
    append!(row, [k => get(Dict(run.overrides), k, NaN) for k in st.varied])
    started = Base.time()
    sim, failed = simulate(st, run.overrides, run.actions)
    if sim === nothing
        append!(row, [k => NaN for k in vcat("score", first.(_rms_names(st)), ["res." * g.name for g in st.scored], first.(st.metrics))])
        append!(row, ["t_end" => NaN, "failed" => failed, "unconverged" => 0, "mass_error" => NaN, "energy_error" => NaN, "vented_kg" => NaN])
    else
        s = score_run(st, sim)
        push!(row, "score" => study_score(s))
        append!(row, s.rms)
        append!(row, s.z)
        p_atm = st.scene.settings["P_atm_kPa"] * 1e3
        append!(row, [name => metric(sim, m, p_atm) for (name, m) in st.metrics])
        r = run_summary(sim)
        append!(row, ["t_end" => r.t_s, "failed" => replace(failed, ',' => ';'), "unconverged" => r.unconverged, "mass_error" => r.mass_error,
                      "energy_error" => r.energy_error, "vented_kg" => r.vented_kg])
        keep === nothing || write_csv(keep, sim; columns = unique(vcat([c.model for c in st.channels], [m["signal"] for (_, m) in st.metrics])))
    end
    push!(row, "wall_s" => round(Base.time() - started, digits = 3))
    row
end
_rms_names(st::Study) = ["rms.$(c.log).$p" => NaN for c in st.channels for (p, _) in st.phases]

"""
    study(path; resume = false, keep_traces = false) -> results file

Run every member of the study in `path`, on as many threads as Julia was started with, and append
one row per run to `<study>_results/results.csv` as it finishes. A run that fails is a row with
the reason in `failed`. With `resume`, runs whose id is already in the file are skipped.
"""
function study(path; resume = false, keep_traces = false, quiet = false)
    st = read_study(path)
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
    Threads.@threads :dynamic for k in eachindex(todo)
        run = todo[k]
        index = findfirst(r -> r.id == run.id, st.runs)
        row = run_one(st, run; keep = keep_traces ? joinpath(st.dir, "traces", "$index.csv") : nothing)
        lock(lock_) do
            open(out, "a") do io
                header && (println(io, join(first.(row), ",")); header = false)
                println(io, join((v isa AbstractFloat ? (isnan(v) ? "" : string(round(v, sigdigits = 7))) : string(v) for v in last.(row)), ","))
            end
        end
    end
    out
end
