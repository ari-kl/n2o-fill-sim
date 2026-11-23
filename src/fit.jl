# Fitting: find the parameter values that bring the model closest to one or several test logs at
# once, and say how well those logs actually determine them. A fit is described by a TOML file;
# see docs/FITTING.md.
#
# Two things matter more than the optimiser. Fitting several logs together is what exposes a model
# that cannot satisfy them all with one set of values. And reporting which parameters the data
# cannot tell apart is what stops a fit from handing back a confident number for something it never
# measured.

using LinearAlgebra: svd, pinv, diag, Diagonal

"One fitted quantity: a scene parameter, its bounds, and the datasets whose scenes it applies to."
struct FitParam
    label::String                 # "VENT_ORIFICE.Cd", or "BOTTLE.fill_kg@hf4" when it belongs to some datasets only
    name::String                  # PART.param or settings.key
    lo::Float64
    hi::Float64
    start::Float64
    log::Bool                     # search in even ratios (for a quantity known only to an order of magnitude)
    datasets::Vector{Int}
end

"A fit file, read and checked."
struct Fit
    name::String
    dir::String
    datasets::Vector{Study}
    names::Vector{String}
    weights::Vector{Float64}
    params::Vector{FitParam}
    explore::Int                  # runs spread over the bounds before the search starts
    seed::Int
    max_evals::Int
    method::String                # "lm" (Levenberg-Marquardt) or "nelder-mead"
    jac_reltol::Float64           # solver tolerance for the sensitivity runs
    jac_step::Float64             # their step, as a fraction of each parameter's range
    held_name::String             # a dataset left out of the fit, to be predicted afterwards ("" if none)
    held::Union{Nothing,Study}
    held_shared::Vector{Int}      # the fitted parameters that apply to it
    held_fixed::Vector{Pair{String,Float64}}   # its own parameters, which stay at the start values given for them
end

const _FIT_KEYS = ("name", "dataset", "parameter", "explore", "seed", "max_evals", "method", "jacobian_reltol", "jacobian_step")
const _DATASET_KEYS = ("name", "weight", "scene", "replay", "manifest", "t0", "until", "sample_dt", "reltol", "wall_limit", "work_limit", "score", "targets", "set")
const _PARAM_KEYS = ("name", "from", "to", "start", "log", "datasets")

# Position of a value within its bounds, 0..1 (on a ratio scale for `log` parameters), and back.
_unit(p::FitParam, v) = p.log ? log(v / p.lo) / log(p.hi / p.lo) : (v - p.lo) / (p.hi - p.lo)
_value(p::FitParam, x) = p.log ? p.lo * (p.hi / p.lo)^x : p.lo + (p.hi - p.lo) * x

"""
    read_fit(path; holdout = nothing) -> Fit

Read a fit file: the datasets (each a scene, a log and what to score, exactly as in a study), and
the parameters to fit with their bounds. A parameter applies to every dataset unless it names the
ones it belongs to; the same name listed twice with different datasets is fitted separately for each.

`holdout` names a dataset to leave out of the fit. Parameters that belong to it alone are not
fitted; shared ones are, on the other datasets.
"""
function read_fit(path; holdout = nothing)
    d = TOML.parsefile(path)
    for k in keys(d)
        k in _FIT_KEYS || error("$path: unknown key \"$k\". Keys: $(join(_FIT_KEYS, ", "))")
    end
    base = dirname(abspath(path))
    dir = splitext(abspath(path))[1] * "_results"
    isempty(get(d, "dataset", [])) && error("$path: no [[dataset]]")
    names = String[String(get(ds, "name", "dataset$k")) for (k, ds) in enumerate(d["dataset"])]
    allunique(names) || error("$path: two datasets share a name")
    datasets = Study[_study(Dict(k => v for (k, v) in ds if !(k in ("name", "weight"))), base, "$path, dataset $(names[k])", dir, _DATASET_KEYS)
                     for (k, ds) in enumerate(d["dataset"])]
    weights = Float64[Float64(get(ds, "weight", 1.0)) for ds in d["dataset"]]
    held, held_shared, held_fixed = nothing, Int[], Pair{String,Float64}[]
    if holdout !== nothing                                     # take it out here, so that everything below sees only the fitted datasets
        h = findfirst(==(holdout), names)
        h === nothing && error("$path: no dataset called \"$holdout\" to hold out. Datasets: $(join(names, ", "))")
        held = datasets[h]
        foreach(v -> deleteat!(v, h), (datasets, names, weights))
        isempty(datasets) && error("$path: holding out $holdout leaves nothing to fit")
    end
    for (k, st) in enumerate(datasets)
        n_scored(st) == 0 && error("$path: dataset $(names[k]) scores nothing (it needs a log and [dataset.score] channels, or a manifest with maps_to)")
    end
    params = FitParam[]
    for p in get(d, "parameter", [])
        for k in keys(p)
            k in _PARAM_KEYS || error("$path: unknown parameter key \"$k\". Keys: $(join(_PARAM_KEYS, ", "))")
        end
        name = String(p["name"])
        lo, hi = Float64(p["from"]), Float64(p["to"])
        hi > lo || error("$path: parameter $name needs from < to")
        logscale = Bool(get(p, "log", false))
        logscale && lo <= 0 && error("$path: parameter $name is on a log scale, so its bounds must be positive")
        listed = haskey(p, "datasets") ? String.(p["datasets"]) : nothing
        on_held = held !== nothing && (listed === nothing || holdout in listed)
        on_held && apply_overrides!(deepcopy(held.scene), ["$name=$lo"])
        listed === nothing || filter!(!=(holdout), listed)
        if listed !== nothing && isempty(listed)               # it belongs to the held-out dataset alone: not fitted
            haskey(p, "start") && push!(held_fixed, name => Float64(p["start"]))
            continue
        end
        which = listed !== nothing ? [something(findfirst(==(n), names), 0) for n in listed] : collect(eachindex(datasets))
        all(>(0), which) || error("$path: parameter $name names a dataset that does not exist. Datasets: $(join(names, ", "))")
        for k in which                                         # the name must exist in each scene it applies to
            apply_overrides!(deepcopy(datasets[k].scene), ["$name=$lo"])
        end
        start = haskey(p, "start") ? Float64(p["start"]) : NaN
        if isnan(start)                                        # default: what the first scene it applies to has now
            s = datasets[which[1]].scene
            m = match(r"^(.+)\.([^.]+)$", name)
            start = m[1] == "settings" ? s.settings[m[2]] : s.parts[find_part(s, m[1])].p[m[2]]
        end
        label = haskey(p, "datasets") ? name * "@" * join(names[which], "+") : name
        push!(params, FitParam(label, name, lo, hi, clamp(start, lo, hi), logscale, which))
        on_held && push!(held_shared, length(params))
    end
    isempty(params) && error("$path: no [[parameter]]")
    allunique(p.label for p in params) || error("$path: a parameter is listed twice for the same datasets")
    for k in eachindex(datasets), name in unique(p.name for p in params)
        count(p -> p.name == name && k in p.datasets, params) <= 1 || error("$path: $name is fitted twice for dataset $(names[k])")
    end
    Fit(String(get(d, "name", splitext(basename(path))[1])), dir, datasets, names, weights,
        params, Int(get(d, "explore", 0)), Int(get(d, "seed", 1)), Int(get(d, "max_evals", 60 * length(params))),
        (m = String(get(d, "method", "lm")); m in ("lm", "nelder-mead") || error("$path: method must be \"lm\" or \"nelder-mead\""); m),
        Float64(get(d, "jacobian_reltol", 1e-7)), Float64(get(d, "jacobian_step", 0.01)), something(holdout, ""), held, held_shared, held_fixed)
end

"""
    evaluate(fit, x; reltol) -> (residuals, readings, failed)

Run every dataset with the parameters at position `x` (each 0..1 within its bounds). `residuals`
is all datasets' residual vectors end to end, each scaled by the square root of its dataset's
weight; its sum of squares is what the fit minimises. `readings` is, for each dataset, the
per-channel, per-phase RMS in physical units, then the offset found for each channel whose zero is
free, then each target's value in the model and in the log. A dataset that fails counts as far from its log, so the search
moves away from it.
"""
function evaluate(fit::Fit, x; reltol = nothing)
    n = length(fit.datasets)
    parts = Vector{Vector{Float64}}(undef, n)
    rms = Vector{Vector{Pair{String,Float64}}}(undef, n)
    failed = fill("", n)
    Threads.@threads :dynamic for k in 1:n
        st = fit.datasets[k]
        overrides = [p.name => _value(p, x[j]) for (j, p) in enumerate(fit.params) if k in p.datasets]
        sim, why = simulate(st, overrides; reltol = something(reltol, st.reltol))
        if sim === nothing || !isempty(why)
            failed[k] = why
            # as bad as 100 sigma in every scored number, in a vector as long as a finished run's
            parts[k] = fill(100.0 * sqrt(fit.weights[k] * n_scored(st) / n_residuals(st)), n_residuals(st))
            rms[k] = vcat(_rms_names(st), _offset_names(st), _target_names(st))
        else
            s = score_run(st, sim)
            rms[k] = vcat(s.rms, s.offsets, s.targets)
            parts[k] = s.residuals .* sqrt(fit.weights[k])
        end
    end
    reduce(vcat, parts), rms, failed
end

# Several positions at once, each of them several datasets run side by side.
_evaluate_all(fit::Fit, xs; reltol = nothing) = fetch.([Threads.@spawn evaluate(fit, x; reltol) for x in xs])

"The fit's one number: as a study's score, over every dataset's scored numbers together."
fit_score(fit::Fit, residuals) = sqrt(sum(abs2, residuals) / sum(n_scored, fit.datasets))

"""
    nelder_mead(f, x0; step, max_evals, xtol, ftol, batch) -> (x, fx, evals)

Minimise `f` over the unit cube from `x0` by the Nelder-Mead simplex method (with the
dimension-dependent coefficients of Gao and Han, 2012, which keep it working beyond a few
parameters). No derivatives are needed, which suits a model whose output has solver noise. Points
are kept inside the cube by clamping.

`batch`, if given, evaluates several points at once (`batch(points) -> values`), and is how spare
processor cores are used. The method's next step depends on the last, so it cannot simply be run
in parallel; but the points it might ask for next (reflected, expanded, contracted outside and
inside) are all known before any of them is evaluated, so with `batch` all four are evaluated
together and the unneeded ones thrown away. The path taken and the answer are exactly those of
the one-at-a-time method, and `evals` counts only the evaluations that method would have made.
"""
function nelder_mead(f, x0; step = 0.1, max_evals = 500, xtol = 2e-3, ftol = 1e-4, batch = nothing)
    n = length(x0)
    α, β, γ, δ = 1.0, 1 + 2 / n, 0.75 - 1 / (2n), 1 - 1 / n
    clampx(x) = clamp.(x, 0.0, 1.0)
    all_at_once(points) = batch === nothing ? [f(x) for x in points] : batch(points)
    X = [copy(x0) for _ in 0:n]
    for k in 1:n
        X[k+1][k] += X[k+1][k] + step <= 1 ? step : -step      # step inwards if the bound is in the way
    end
    F = all_at_once(X)
    evals = n + 1
    while evals < max_evals
        order = sortperm(F)
        X, F = X[order], F[order]
        size = maximum(maximum(abs.(X[k] .- X[1])) for k in 2:n+1)
        (size < xtol && F[end] - F[1] < ftol * max(F[1], 1e-12)) && break
        centre = sum(X[1:n]) ./ n
        xr = clampx(centre .+ α .* (centre .- X[end]))
        xe = clampx(centre .+ β .* (xr .- centre))
        xo = clampx(centre .+ γ .* (xr .- centre))             # contraction on the reflected side
        xi = clampx(centre .+ γ .* (X[end] .- centre))         # contraction on the worst point's side
        ahead = batch === nothing ? nothing : batch([xr, xe, xo, xi])
        fr = ahead === nothing ? f(xr) : ahead[1]; evals += 1
        if fr < F[1]
            fe = ahead === nothing ? f(xe) : ahead[2]; evals += 1
            X[end], F[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < F[n]
            X[end], F[end] = xr, fr
        else
            outside = fr < F[end]
            xc = outside ? xo : xi
            fc = ahead === nothing ? f(xc) : ahead[outside ? 3 : 4]; evals += 1
            if fc < (outside ? fr : F[end])
                X[end], F[end] = xc, fc
            else                                               # shrink towards the best point
                for k in 2:n+1
                    X[k] = clampx(X[1] .+ δ .* (X[k] .- X[1]))
                end
                F[2:end] .= all_at_once(X[2:end])
                evals += n
            end
        end
    end
    k = argmin(F)
    X[k], F[k], evals
end

"""
    levenberg_marquardt(batch, x0; step, max_evals, xtol, ftol) -> (x, S, evals)

Minimise the sum of squares S of a residual vector over the unit cube from `x0`, by the
Levenberg-Marquardt method. `batch(points)` returns the residual vector at each point (or `nothing`
where the model failed).

This is the method for fitting a model to data, and it suits several processor cores: one
iteration needs the residuals at one nearby point per parameter (for the sensitivities), all
independent, then a few trial steps of different lengths, also independent. A simplex search asks
for one point at a time and needs ten or more times as many of them.

The sensitivities are plain differences over `step` of each parameter's range (taken inwards at a
bound): large enough to sit well above the solver's own noise. A parameter at a bound that the
data push further out is held there for that iteration. Four dampings are tried at once, from
nearly a Gauss-Newton step to a short step downhill; the best is taken if it lowers S, and if none
does the damping is raised. If one of the sensitivity runs is itself better than every trial, that
point is taken instead. It stops when a lightly damped step gains less than `ftol` of S and moves
less than `xtol`, twice running, or when no damping helps.
"""
function levenberg_marquardt(batch, x0; step = 0.02, max_evals = 1000, xtol = 1e-3, ftol = 1e-4, report = (args...) -> nothing)
    n = length(x0)
    x = clamp.(x0, 0.0, 1.0)
    r = batch([x])[1]
    r === nothing && error("the model fails at the starting point of the fit")
    S, evals, λ, quiet_steps = sum(abs2, r), 1, 1e-2, 0
    ok(v) = v !== nothing && length(v) == length(r)
    J, moved_since = zeros(length(r), n), true
    nudge, r_nudge, S_nudge = x, r, Inf                        # the best of the sensitivity runs: a point like any other
    while evals + (moved_since ? n : 0) + 4 <= max_evals
        if moved_since                                         # after a step that was not taken, the sensitivities still stand
            hs = [x[k] + step <= 1 ? step : -step for k in 1:n]
            points = [(y = copy(x); y[k] += hs[k]; y) for k in 1:n]
            cols = batch(points)
            evals += n
            fill!(J, 0.0)
            S_nudge = Inf
            for k in 1:n
                ok(cols[k]) || continue
                J[:, k] .= (cols[k] .- r) ./ hs[k]
                Sk = sum(abs2, cols[k])
                Sk < S_nudge && ((nudge, r_nudge, S_nudge) = (points[k], cols[k], Sk))
            end
            moved_since = false
        end
        g, A = J' * r, J' * J
        free = [!((x[k] <= 0 && g[k] > 0) || (x[k] >= 1 && g[k] < 0)) for k in 1:n]
        idx = findall(free)
        isempty(idx) && break
        scale = max.(diag(A)[idx], 1e-12 * maximum(diag(A)) + 1e-300)
        dampings = λ .* (0.1, 1.0, 10.0, 100.0)
        trials = map(dampings) do l
            δ = zeros(n)
            δ[idx] = -((A[idx, idx] + l * Diagonal(scale)) \ g[idx])
            clamp.(x .+ δ, 0.0, 1.0)
        end |> collect
        rs = batch(trials)
        evals += length(trials)
        Ss = [ok(v) ? sum(abs2, v) : Inf for v in rs]
        b = argmin(Ss)
        if S_nudge < min(Ss[b], S)
            # Moving one parameter by its sensitivity step did better than any step worked out from
            # the sensitivities: the model is not close to linear over that distance. Take it.
            x, r, S, moved_since, quiet_steps = nudge, r_nudge, S_nudge, true, 0
            report(evals, S)
        elseif Ss[b] < S
            gain, moved = (S - Ss[b]) / S, maximum(abs.(trials[b] .- x))
            # A short step means the minimum is near only if the step was a nearly free one. A
            # heavily damped step is short because it was held back.
            free_step = dampings[b] <= 1
            x, r, S, λ, moved_since = trials[b], rs[b], Ss[b], max(dampings[b] / 3, 1e-10), true
            report(evals, S)
            quiet_steps = free_step && moved < xtol && (gain < ftol || S < 1e-18) ? quiet_steps + 1 : 0
            (quiet_steps >= 2 || S == 0) && break
        else
            λ *= 1000                                          # none of them helped: shorter steps
            λ > 1e9 && break                                   # not even the shortest: as low as the noise allows
        end
    end
    x, S, evals
end

"""
    identifiability(J, labels) -> (std, corr, statements)

What a set of data can and cannot say about fitted parameters, from the sensitivity `J` of the
residuals to each parameter's position in its range. `std` is each parameter's approximate
standard error (in the same 0..1 units, to be scaled by the misfit), `corr` their correlations,
and `statements` plain sentences naming every pair the data cannot tell apart (|correlation| over
0.9) and every combination of parameters the data do not determine at all (a singular value of
`J` under a thousandth of the largest).
"""
function identifiability(J, labels)
    n = size(J, 2)
    scale = [max(sqrt(sum(abs2, view(J, :, k))), 1e-300) for k in 1:n]
    U = svd(J ./ scale')                                       # columns scaled, so a weak parameter is not mistaken for a null direction
    statements = String[]
    for k in 1:n
        scale[k] < 1e-9 * maximum(scale) && push!(statements, "$(labels[k]) has no measurable effect on any scored reading: these data say nothing about it.")
    end
    for (i, sv) in enumerate(U.S)
        sv < 1e-3 * U.S[1] || continue
        v = U.V[:, i]
        terms = [(v[k], labels[k]) for k in sortperm(abs.(v); rev = true) if abs(v[k]) > 0.2 && scale[k] >= 1e-9 * maximum(scale)]
        length(terms) < 2 && continue
        push!(statements, "The data do not determine the combination " * join(["$(round(c / terms[1][1], digits = 2))·$l" for (c, l) in terms], " + ") *
                          " (moving along it changes no reading)." )
    end
    cov = pinv(J' * J; rtol = 1e-10)
    std = sqrt.(max.(diag(cov), 0.0))
    corr = [std[i] > 0 && std[j] > 0 ? cov[i, j] / (std[i] * std[j]) : 0.0 for i in 1:n, j in 1:n]
    for i in 1:n, j in i+1:n
        abs(corr[i, j]) > 0.9 && push!(statements, "$(labels[i]) and $(labels[j]) are not separately determined (correlation $(round(corr[i, j], digits = 2))): " *
                                                    "a change in one can be made up by the other.")
    end
    std, corr, statements
end

"""
    fit(path; write_scenes = false, holdout = nothing, quiet = false) -> report

Run the fit in `path`. Optionally explore the bounds first with a Latin hypercube, then search from
the best point found (or from the given starts) with Levenberg-Marquardt (or Nelder-Mead, once and
then again from where it stopped). Then run the sensitivity of every residual to every parameter
at a tighter solver tolerance and report what the data determine. Writes into
`<fit file>_results/`: `evals.csv` (every evaluation of the search), `report.txt`, a calibration
record `<date>-<name>.toml`, and with `write_scenes` a copy of each dataset's scene with the fitted
values in place and a note of where each came from.

`holdout` names a dataset to leave out of the fit and then predict with the fitted values: a check
on a log the fit never saw.
"""
function fit(path; write_scenes = false, holdout = nothing, quiet = false)
    ft = read_fit(path; holdout)
    mkpath(ft.dir)
    np = length(ft.params)
    say(args...) = quiet || println(args...)
    log_io = open(joinpath(ft.dir, "evals.csv"), "w")
    println(log_io, "eval,score,", join((p.label for p in ft.params), ","), ",", join(("score." * n for n in ft.names), ","), ",failed")
    evals = Ref(0)
    best = Ref((Inf, Float64[]))
    # Evaluate several positions at once, then score and log them in the order given, so that the
    # log does not depend on which finished first.
    function scored(xs)                                        # (score, residual vector or nothing if a dataset failed) for each
        out = Tuple{Float64,Union{Nothing,Vector{Float64}}}[]
        for (x, (residuals, _, failed)) in zip(xs, _evaluate_all(ft, xs))
            s = fit_score(ft, residuals)
            evals[] += 1
            s < best[][1] && (best[] = (s, copy(x)))
            println(log_io, evals[], ",", round(s, sigdigits = 7), ",", join((round(_value(p, x[j]), sigdigits = 7) for (j, p) in enumerate(ft.params)), ","), ",",
                    join((round(v, sigdigits = 5) for v in _dataset_scores(ft, residuals)), ","), ",", join(filter(!isempty, failed), " | "))
            push!(out, (s, all(isempty, failed) ? residuals : nothing))
        end
        flush(log_io)
        out
    end
    objectives(xs) = first.(scored(xs))
    objective(x) = objectives([x])[1]
    # With more threads than datasets there are cores to spare: use them to look ahead in the search.
    spare = Threads.nthreads() > length(ft.datasets)
    x0 = [_unit(p, p.start) for p in ft.params]
    f0 = objective(x0)
    say("fit \"", ft.name, "\": ", np, " parameters, ", length(ft.datasets), " datasets, ", Threads.nthreads(), " thread(s). Score at the start: ", round(f0, digits = 3))
    if ft.explore > 0                                          # a Latin hypercube over the bounds, to start from the best basin found
        g = _Lcg(UInt64(ft.seed))
        perms = [(p = collect(1:ft.explore); for i in ft.explore:-1:2; j = 1 + floor(Int, _next!(g) * i); p[i], p[j] = p[j], p[i]; end; p) for _ in 1:np]
        points = [[(perms[k][i] - 1 + _next!(g)) / ft.explore for k in 1:np] for i in 1:ft.explore]
        for chunk in Iterators.partition(points, max(1, Threads.nthreads() ÷ 2))
            objectives(collect(chunk))
        end
        say("after ", ft.explore, " exploratory runs: best ", round(best[][1], digits = 3))
    end
    x, fx = best[][2], best[][1]
    if ft.method == "lm"
        x, _, used = levenberg_marquardt(xs -> last.(scored(xs)), x; max_evals = ft.max_evals,
                                         report = (n, S) -> say("  after ", n, " evaluations: score ", round(best[][1], digits = 4)))
        say("search (Levenberg-Marquardt): score ", round(best[][1], digits = 4), " after ", used, " evaluations")
    else
        for (round_, step) in enumerate((0.1, 0.03))           # a second, tighter pass from where the first stopped
            budget = ft.max_evals ÷ 2
            x, fx, used = nelder_mead(objective, x; step, max_evals = budget, batch = spare ? objectives : nothing)
            say("search ", round_, " (Nelder-Mead): score ", round(fx, digits = 4), " after ", used, " evaluations")
        end
    end
    close(log_io)
    x, fx = best[][2], best[][1]
    # Sensitivities at the optimum, with the solver tightened so that differences are signal, not solver noise.
    r0, rms, failed = evaluate(ft, x; reltol = ft.jac_reltol)
    J = zeros(length(r0), np)
    steps = [x[k] + ft.jac_step <= 1 ? ft.jac_step : -ft.jac_step for k in 1:np]
    moved = [(xk = copy(x); xk[k] += steps[k]; xk) for k in 1:np]
    shifted = reduce(vcat, [_evaluate_all(ft, collect(chunk); reltol = ft.jac_reltol) for chunk in Iterators.partition(moved, max(1, Threads.nthreads() ÷ 2))])
    for k in 1:np
        h = steps[k]
        rk = shifted[k][1]
        J[:, k] .= (rk .- r0) ./ h
    end
    labels = [p.label for p in ft.params]
    std, corr, statements = identifiability(J, labels)
    misfit = sqrt(sum(abs2, r0) / max(sum(n_scored, ft.datasets) - np, 1))     # treating each scored number as one observation
    values = [_value(p, x[k]) for (k, p) in enumerate(ft.params)]
    # an approximate standard error in the parameter's own units (a factor, for a log parameter)
    spread = [p.log ? exp(std[k] * misfit * log(p.hi / p.lo)) : std[k] * misfit * (p.hi - p.lo) for (k, p) in enumerate(ft.params)]
    each = _dataset_scores(ft, r0)                             # each fitted dataset's own score
    # The dataset held out, run with the fitted values of the parameters it shares (its own stay
    # where they started), and for comparison with every fitted parameter at its start value.
    held = nothing
    if ft.held !== nothing
        predict(v) = begin
            sim, why = simulate(ft.held, vcat(ft.held_fixed, [ft.params[j].name => v[j] for j in ft.held_shared]))
            (sim === nothing || !isempty(why)) ? (NaN, vcat(_rms_names(ft.held), _offset_names(ft.held), _target_names(ft.held)), why) :
                (s = score_run(ft.held, sim); (study_score(s), vcat(s.rms, s.offsets, s.targets), ""))
        end
        score, readings, why = predict(values)
        held = (; name = ft.held_name, score, start_score = predict([p.start for p in ft.params])[1], readings, failed = why)
    end
    io = IOBuffer()
    println(io, "Fit \"", ft.name, "\"    score ", round(fx, digits = 4), " (", round(f0, digits = 3), " at the start), ", evals[], " evaluations")
    println(io, "\nParameters (± is approximate: one scored number counted as one observation)")
    for (k, p) in enumerate(ft.params)
        at_bound = x[k] < 0.01 ? "  AT ITS LOWER BOUND" : x[k] > 0.99 ? "  AT ITS UPPER BOUND" : ""
        println(io, "  ", rpad(p.label, 34), rpad(round(values[k], sigdigits = 5), 12), p.log ? "×/÷ $(round(spread[k], sigdigits = 3))" : "± $(round(spread[k], sigdigits = 3))",
                "   (bounds ", p.lo, " to ", p.hi, "; started at ", round(p.start, sigdigits = 5), ")", at_bound)
    end
    println(io, "\nHow close, per dataset, channel and phase (RMS in the log's own units; offset: what a channel with a free zero was shifted by; target: a feature in the model and in the log)")
    for (k, name) in enumerate(ft.names)
        println(io, "  ", name, "    score ", round(each[k], sigdigits = 4), isempty(failed[k]) ? "" : "   FAILED: " * failed[k])
        for (key, v) in rms[k]
            println(io, "    ", rpad(replace(key, "rms." => ""), 34), round(v, sigdigits = 4))
        end
    end
    if held !== nothing
        println(io, "\nHeld out: ", held.name, " took no part in the fit. Run with the fitted values of the parameters it shares:")
        println(io, "  score ", round(held.score, sigdigits = 4), "  (", round(held.start_score, sigdigits = 4), " with every parameter at its start value)", isempty(held.failed) ? "" : "   FAILED: " * held.failed)
        for (key, v) in held.readings
            println(io, "    ", rpad(replace(key, "rms." => ""), 34), round(v, sigdigits = 4))
        end
        isempty(ft.held_fixed) || println(io, "  its own parameters stayed at: ", join(["$k = $v" for (k, v) in ft.held_fixed], ", "))
        println(io, "  A fit that scores well on its own logs and badly here has fitted those logs, not the hardware.")
    end
    println(io, "\nWhat these data do and do not determine")
    isempty(statements) ? println(io, "  Every parameter is separately determined (no correlation over 0.9, no undetermined combination).") :
                          foreach(s -> println(io, "  - ", s), statements)
    println(io, "\nSensitivities were taken with steps of ", ft.jac_step, " of each parameter's range at solver tolerance ", ft.jac_reltol, ".")
    report = String(take!(io))
    write(joinpath(ft.dir, "report.txt"), report)
    record = joinpath(ft.dir, "$(_today())-$(ft.name).toml")
    open(record, "w") do out
        println(out, "# Calibration record written by `fit`. Values, what they were fitted to, and how well.")
        println(out, "name = ", repr(ft.name), "\ndate = ", repr(_today()), "\ncommit = ", repr(_commit(dirname(abspath(path)))), "\nmodel_commit = ", repr(_commit(@__DIR__)), "\nscore = ", round(fx, sigdigits = 6))
        println(out, "datasets = ", _list(ft.names))
        for (k, p) in enumerate(ft.params)
            println(out, "\n[[parameter]]\nname = ", repr(p.name), "\ndatasets = ", repr(ft.names[p.datasets]), "\nvalue = ", round(values[k], sigdigits = 7),
                    "\nfrom = ", p.lo, "\nto = ", p.hi, "\nlog = ", p.log, "\nspread = ", round(spread[k], sigdigits = 4))
        end
        for (k, name) in enumerate(ft.names)
            st = ft.datasets[k]
            println(out, "\n[[dataset]]\nname = ", repr(name), "\nweight = ", ft.weights[k], "\nuntil = ", st.until, "\nscore = ", round(each[k], sigdigits = 6))
            println(out, "channels = [", join(["{ log = $(repr(c.log)), model = $(repr(c.model)), sigma = $(c.sigma), weight = $(c.weight), offset = $(repr(c.offset)) }" for c in st.channels], ", "), "]")
            println(out, "phases = { ", join(["$(p) = [$(a), $(b)]" for (p, (a, b)) in st.phases], ", "), " }")
            println(out, "scored = [", join(["{ name = $(repr(g.name)), sigma = $(round(g.sigma, sigdigits = 6)) }" for g in st.scored], ", "), "]")
            isempty(st.targets) || println(out, "targets = [", join(["{ name = $(repr(g.name)), " * join(["$k = $(v isa AbstractString ? repr(v) : v)" for (k, v) in sort!(collect(g.m); by = first)], ", ") * " }" for g in st.targets], ", "), "]")
            for kind in ("rms", "offset", "target")
                println(out, kind, " = { ", join(["$(repr(replace(key, kind * "." => ""))) = $(isnan(v) ? "nan" : round(v, sigdigits = 5))" for (key, v) in rms[k] if startswith(key, kind * ".")], ", "), " }")
            end
        end
        if held !== nothing
            println(out, "\n[holdout]                       # not fitted: predicted with the fitted values")
            println(out, "name = ", repr(held.name), "\nscore = ", isnan(held.score) ? "nan" : round(held.score, sigdigits = 6), "\nstart_score = ", isnan(held.start_score) ? "nan" : round(held.start_score, sigdigits = 6))
            println(out, "readings = { ", join(["$(repr(key)) = $(isnan(v) ? "nan" : round(v, sigdigits = 5))" for (key, v) in held.readings], ", "), " }")
        end
        println(out, "\nidentifiability = ", _list(statements))
    end
    if write_scenes
        for (k, name) in enumerate(ft.names)
            st = ft.datasets[k]
            scene = deepcopy(st.scene)
            apply_overrides!(scene, ["$(p.name)=$(round(values[j], sigdigits = 6))" for (j, p) in enumerate(ft.params) if k in p.datasets])
            file = joinpath(ft.dir, name * ".toml")
            save_scene(file, scene)
            notes = ["#   $(p.name) = $(round(values[j], sigdigits = 6))" for (j, p) in enumerate(ft.params) if k in p.datasets]
            write(file, "# Fitted values written by `fit` \"$(ft.name)\" on $(_today()) [cal $(_today()): $(join(ft.names, ", "))]:\n" * join(notes, "\n") * "\n" * read(file, String))
        end
    end
    say("\n", report, "\nwritten: ", joinpath(ft.dir, "report.txt"), ", ", record)
    (; score = fx, start_score = f0, values = Dict(labels .=> values), spread = Dict(labels .=> spread), statements, rms = Dict(ft.names .=> rms), evals = evals[], dir = ft.dir,
       scores = Dict(ft.names .=> each), holdout = held)
end

# Each dataset's own score from the combined residual vector.
function _dataset_scores(ft::Fit, residuals)
    out, offset = Float64[], 0
    for (k, st) in enumerate(ft.datasets)
        push!(out, sqrt(sum(abs2, view(residuals, offset+1:offset+n_residuals(st))) / (n_scored(st) * ft.weights[k])))
        offset += n_residuals(st)
    end
    out
end

_list(v) = "[" * join(repr.(v), ", ") * "]"        # a list of strings as TOML (an empty one too)
_today() = (t = Libc.strftime("%Y-%m-%d", time()); t)
_commit(dir) = try strip(read(setenv(`git rev-parse --short HEAD`; dir), String)) catch; "unknown" end

"""
    fit_skeleton(scene_path) -> String

A fit file to start from: one parameter for every line of the scene that carries a `[fit]` tag,
with bounds at half and one and a half times its present value.
"""
function fit_skeleton(scene_path)
    io = IOBuffer()
    println(io, "# Draft written by `fit --init`. Add the log to fit to, check the bounds, remove what should stay fixed.")
    println(io, "name = ", repr(splitext(basename(scene_path))[1]), "\nexplore = 40\n\n[[dataset]]\nname = \"log1\"\nscene = ", repr(basename(scene_path)))
    println(io, "# replay = \"log.csv\"\n# t0 = 0.0\n# [dataset.score]\n# channels = [{ log = \"PT202_psig\", model = \"TANK.P\" }]\n# phases = { all = [0, 600] }")
    part = ""
    for line in eachline(scene_path)
        m = match(r"^\s*name\s*=\s*\"([^\"]+)\"", line)
        m === nothing || (part = m[1])
        startswith(strip(line), "[settings]") && (part = "settings")
        occursin("[fit", line) || continue
        kv = match(r"^\s*([A-Za-z0-9_]+)\s*=\s*([-+0-9.eE]+)", line)
        (kv === nothing || isempty(part)) && continue
        v = parse(Float64, kv[2])
        lo, hi = minmax(0.5v, 1.5v)
        println(io, "\n[[parameter]]\nname = ", repr("$part.$(kv[1])"), "\nfrom = ", round(lo, sigdigits = 4), "\nto = ", round(hi, sigdigits = 4), "        # now ", v)
    end
    String(take!(io))
end
