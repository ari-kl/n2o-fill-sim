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
    runs::Vector{StudyRun}
    varied::Vector{String}
    channels::Vector{Tuple{String,String,Float64,Float64}}   # (log column, model signal, sigma in the log's units, weight)
    phases::Vector{Pair{String,Tuple{Float64,Float64}}}
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

const _STUDY_KEYS = ("scene", "replay", "manifest", "t0", "until", "sample_dt", "reltol", "wall_limit", "mode", "n", "seed",
                     "vary", "case", "score", "metrics", "set")
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
    channels = Tuple{String,String,Float64,Float64}[(String(c["log"]), String(c["model"]), trust(c["log"], get(c, "sigma", NaN)), Float64(get(c, "weight", 1.0)))
                                                    for c in get(score, "channels", [])]
    isempty(channels) && !haskey(score, "channels") && append!(channels, [(c.name, c.target, trust(c.name, NaN), 1.0) for c in columns if !isempty(c.target)])
    all(c -> c[3] > 0, channels) || error("$path: a channel's sigma must be positive")
    if haskey(score, "phases")
        sched = sort(vcat(schedule, scene.schedule); by = a -> a.t)
        phases = [String(k) => (p isa AbstractVector ? (Float64(p[1]), Float64(p[2])) :
                                (f = _phase_time(p["from"], sched); (f, _phase_time(p["to"], sched; after = f)))) for (k, p) in score["phases"]]
    end
    isempty(phases) && push!(phases, "all" => (0.0, until))
    sort!(phases; by = first)
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
          Float64(get(d, "reltol", 1e-5)), Float64(get(d, "wall_limit", 120.0)), runs, names_v, channels, phases, metrics)
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
    score_run(study, sim) -> (rms, residuals)

Compare a run with the study's log. The model is interpolated onto the log's times, and times with
no log reading are skipped.

`rms` has one entry per scored channel and phase, in the log's own units (psi, °C, kg): the number
to read. `residuals` is what gets added up into the one-number `score`, and what a fitter
minimises. Each difference is divided by how far that reading is trusted (`sigma`), so that a
pressure and a temperature can be weighed against each other, and by the square root of the number
of samples in its phase, so that a long quiet hold does not outvote a ten-second vent. The sum of
squares of `residuals` is then the sum over channels and phases of (weight · rms / sigma)².
"""
function score_run(st::Study, sim::Sim)
    rms, residuals = Pair{String,Float64}[], Float64[]
    p_atm = st.scene.settings["P_atm_kPa"] * 1e3
    t_end = simtime(sim)
    for (logname, model, sigma, weight) in st.channels
        haskey(sim.hist, model) || error("no recorded signal named $model. Signals: $(join(sort!(collect(keys(sim.hist))), ", "))")
        role = st.columns[findfirst(c -> c.name == logname, st.columns)].role
        y = st.cols[logname]
        for (phase, (a, b)) in st.phases
            first_ = length(residuals) + 1
            e2 = 0.0
            for (i, t) in enumerate(st.t_log)
                (a <= t <= b && t <= t_end && isfinite(y[i])) || continue
                d = _as_logged(role, signal_at(sim, model, t), p_atm) - y[i]
                e2 += d^2
                push!(residuals, weight * d / sigma)
            end
            n = length(residuals) - first_ + 1
            n > 0 && (@views residuals[first_:end] ./= sqrt(n))
            push!(rms, "rms.$logname.$phase" => (n == 0 ? NaN : sqrt(e2 / n)))
        end
    end
    rms, residuals
end

"""
The one-number score of a run: the root mean square, over every scored channel and phase, of its
RMS difference in units of that reading's `sigma`. About 1 means the model is, on average, as far
from the log as the readings themselves can be trusted; 3 means three times that.
"""
function study_score(st::Study, rms)
    total, pairs, k = 0.0, 0, 0
    for (_, _, sigma, weight) in st.channels, _ in st.phases
        k += 1
        isnan(last(rms[k])) && continue            # a phase with no log reading on this channel
        total += (weight * last(rms[k]) / sigma)^2
        pairs += 1
    end
    pairs == 0 ? NaN : sqrt(total / pairs)
end

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
        run_until!(sim, st.until; wall_limit = st.wall_limit)
        return sim, !isempty(sim.failed) ? sim.failed : simtime(sim) < st.until - 1e-6 ? "stopped at the wall-time limit" : ""
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
        append!(row, [k => NaN for k in vcat("score", [first(r) for r in _rms_names(st)], first.(st.metrics))])
        append!(row, ["t_end" => NaN, "failed" => failed, "unconverged" => 0, "mass_error" => NaN, "energy_error" => NaN, "vented_kg" => NaN])
    else
        rms, residuals = isempty(st.channels) ? (Pair{String,Float64}[], Float64[]) : score_run(st, sim)
        push!(row, "score" => (isempty(st.channels) ? NaN : study_score(st, rms)))
        append!(row, rms)
        p_atm = st.scene.settings["P_atm_kPa"] * 1e3
        append!(row, [name => metric(sim, m, p_atm) for (name, m) in st.metrics])
        r = run_summary(sim)
        append!(row, ["t_end" => r.t_s, "failed" => replace(failed, ',' => ';'), "unconverged" => r.unconverged, "mass_error" => r.mass_error,
                      "energy_error" => r.energy_error, "vented_kg" => r.vented_kg])
        keep === nothing || write_csv(keep, sim; columns = unique(vcat([c[2] for c in st.channels], [m["signal"] for (_, m) in st.metrics])))
    end
    push!(row, "wall_s" => round(Base.time() - started, digits = 3))
    row
end
_rms_names(st::Study) = ["rms.$(c[1]).$p" => NaN for c in st.channels for (p, _) in st.phases]

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
