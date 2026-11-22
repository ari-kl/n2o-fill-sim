# Test logs: reading a telemetry CSV, working out what its columns are, and turning it into an
# overlay for the plots and a valve schedule for the simulator.
#
# Column meaning is the trap. The same valves have been logged as `valve_AV201` with 1 = open
# (HF3, HF4) and, in later logs, as `solenoid_AV201` with 1 = "CLOSED". So a 0/1 column is never replayed on a guess: `LogColumn.one_means` has to say.

"""
    read_csv(path; min_dt = 0.0) -> (names, columns)

Read a CSV with a header row. The first column is taken as time; with `min_dt` > 0 rows closer
together than that are skipped, which is how a 200 Hz, 40 MB log is brought down to something to
plot. Anything that is not a number (a blank cell, a device name) reads as NaN.
"""
function read_csv(path; min_dt = 0.0)
    names = String[]
    cols = Vector{Float64}[]
    last = -Inf
    for (i, line) in enumerate(eachline(path))
        if i == 1
            names = String.(strip.(split(line, ',')))
            cols = [Float64[] for _ in names]
            continue
        end
        c = findfirst(',', line)
        c === nothing && continue
        t = tryparse(Float64, SubString(line, 1, c - 1))      # look at the time before splitting the row
        t === nothing && continue
        (t - last >= min_dt || t < last) || continue
        last = t
        f = split(line, ',')
        for (k, col) in enumerate(cols)
            push!(col, k <= length(f) ? something(tryparse(Float64, f[k]), NaN) : NaN)
        end
    end
    names, Dict(zip(names, cols))
end

"""
    schedule_from_log(t, cols, map; t0)

Turn logged 0/1 columns into a schedule. `map` pairs a column name with (part, action when 1,
action when 0), e.g. "AV203" => ("AV203", "open", "close"). The first sample sets every mapped
part's starting state. Missing samples (NaN) are skipped. Times are shifted so `t0` is zero.
"""
function schedule_from_log(t, cols, map; t0 = t[1])
    out = Action[]
    for (col, (part, on, off)) in map
        v = cols[col]
        prev = NaN
        for i in eachindex(t)
            (t[i] < t0 || isnan(v[i]) || v[i] == prev) && continue
            push!(out, Action(max(t[i] - t0, 0.0), part, v[i] > 0.5 ? on : off, 0.0))
            prev = v[i]
        end
    end
    sort!(out; by = a -> a.t)
end

"What one column of a log is taken to be. Filled in by `guess_log_columns`, corrected by the user."
mutable struct LogColumn
    name::String
    role::Symbol                  # :time, :pressure (psi gauge), :mass (kg), :temperature (°C), :valve, :heater, :other
    part::String                  # the scene part a :valve or :heater column drives ("" if none)
    one_means::Symbol             # valves: :open, :closed, or :unknown (must be chosen before replay)
    use::Bool                     # overlay it (readings) / replay it (valve, heater)
    target::String                # the recorded signal a reading is to be compared with, e.g. "TANK.P" ("" if none)
    sigma::Float64                # how far the reading can be trusted, in its own units (psi, °C, kg); NaN: the default for its role
end
LogColumn(name, role, part, one_means, use) = LogColumn(name, role, part, one_means, use, "", NaN)
LogColumn(name, role, part, one_means, use, target) = LogColumn(name, role, part, one_means, use, target, NaN)

"""
How far a reading of each kind is trusted when nothing better is said: 5 psi for a pressure, 1 K
for a temperature, 0.25 kg for a mass. A difference between model and log is counted in these
units, which is what lets a pressure and a temperature be weighed against each other.
"""
default_sigma(role::Symbol) = role === :pressure ? 5.0 : role === :temperature ? 1.0 : role === :mass ? 0.25 : 1.0

_bare(name) = uppercase(replace(replace(name, r"\[.*?\]" => "", r"^(valve|solenoid|relay|kasa)_" => ""), r"[^A-Za-z0-9]" => ""))

"""
    guess_log_columns(names, scene) -> Vector{LogColumn}

A first guess at each column of a log, by its name: units in brackets or a `_psig` / `_kg` suffix
for readings, and 0/1 columns matched to valves and heaters in the scene (`valve_AV-201` and
`solenoid_AV201` both match a part called `AV201`). `valve_` columns and bare valve names are
taken as 1 = open; `solenoid_` and `relay_` columns are left as :unknown.
"""
function guess_log_columns(names, scene::Scene)
    valves = Dict(_bare(p.name) => p.name for p in scene.parts if p.kind == "Valve")
    heated(kind) = (i = findfirst(p -> p.kind == kind && p.p["heater_W"] > 0, scene.parts); i === nothing ? "" : scene.parts[i].name)
    out = LogColumn[]
    for (k, name) in enumerate(names)
        low = lowercase(name)
        col = LogColumn(name, :other, "", :unknown, false)
        if k == 1
            col.role = :time
        elseif occursin("[psi]", low) || endswith(low, "_psig")
            col.role, col.use = :pressure, true
        elseif occursin("[kg]", low) || occursin("[kilograms]", low) || endswith(low, "_kg")
            col.role, col.use = :mass, occursin("fill", low)          # a thrust cell in kg is not tank mass
        elseif occursin("[c]", low) || occursin("[degc]", low) || endswith(low, "_c")
            col.role, col.use = :temperature, occursin("tank", low) || occursin("th", low)   # not every temperature is the tank's
        elseif occursin("heater", low)
            col.role = :heater
            col.part = heated(occursin("supply", low) ? "Bottle" : "Tank")
            col.one_means, col.use = :open, !isempty(col.part)        # for a heater, "open" reads as on
        elseif haskey(valves, _bare(name))
            col.role, col.part = :valve, valves[_bare(name)]
            col.one_means = startswith(low, "solenoid_") || startswith(low, "relay_") ? :unknown : :open
            col.use = col.one_means !== :unknown
        end
        push!(out, col)
    end
    out
end

"The schedule a log's valve and heater columns give, for the columns marked `use` with a known meaning."
function log_schedule(cols, columns::Vector{LogColumn}; t0)
    t = cols[columns[1].name]
    map = Dict{String,Tuple{String,String,String}}()
    for c in columns
        (c.use && !isempty(c.part) && c.one_means !== :unknown) || continue
        on, off = c.role === :heater ? ("on", "off") : c.role === :valve ? ("open", "close") : continue
        map[c.name] = c.one_means === :open ? (c.part, on, off) : (c.part, off, on)
    end
    schedule_from_log(t, cols, map; t0)
end

"Saturation temperature (K) at an absolute pressure, or NaN outside the tabulated range or for a missing reading."
tsat_or_nan(P) = isfinite(P) && SAT_PMIN <= P <= SAT_PMAX ? Tsat(P) : NaN
"Saturation pressure (Pa) at a temperature, or NaN outside the tabulated range or for a missing reading."
psat_or_nan(T) = isfinite(T) && N2O_TTRIP <= T <= SAT_TMAX ? Psat(T) : NaN

"""
The readings a log offers for plotting next to the simulation: (label, kind, values in SI) for
each pressure, mass and temperature column marked `use`. Logged pressures are psi gauge and
temperatures °C.

With `saturation = true` each pressure also gives the saturation temperature it implies and each
temperature the saturation pressure it implies. In a tank holding liquid the two should agree: a
thermometer on the wall that reads the saturation temperature of the tank pressure is following
the liquid, and a pressure above the saturation pressure of the bulk temperature points at a warm
surface layer or a sensor fault.
"""
function log_overlay(cols, columns::Vector{LogColumn}, p_atm; saturation = false)
    out = Tuple{String,Symbol,Vector{Float64}}[]
    for c in columns
        c.use || continue
        v = cols[c.name]
        if c.role === :pressure
            push!(out, (c.name, :pressure, v .* PSI .+ p_atm))
            saturation && push!(out, ("Tsat of " * c.name, :temperature, tsat_or_nan.(v .* PSI .+ p_atm)))
        elseif c.role === :mass
            push!(out, (c.name, :mass, copy(v)))
        elseif c.role === :temperature
            push!(out, (c.name, :temperature, v .+ 273.15))
            saturation && push!(out, ("Psat of " * c.name, :pressure, psat_or_nan.(v .+ 273.15)))
        end
    end
    out
end

# ---------------------------------------------------------------- manifests
#
# A manifest says what a log is, so that nothing has to be guessed and nothing about one test's
# conventions lives in code: which column is time, which boards wrote the rows and how their clocks
# relate, what each channel is and what a 1 means on a valve column, and where the phases of the
# test begin and end. docs/LOGS.md describes the file.

"A log read through a manifest: evenly spaced readings in the units the rest of the code expects."
struct Log
    names::Vector{String}                         # the time column first
    cols::Dict{String,Vector{Float64}}            # psi gauge, kg, °C, 0/1; NaN where there was no sample
    columns::Vector{LogColumn}
    t0::Float64                                   # log time that is simulation time zero
    phases::Dict{String,Tuple{Float64,Float64}}   # seconds after t0
    offsets::Dict{String,Float64}                 # seconds added to each source's clock
    warnings::Vector{String}
end

const _MANIFEST_KEYS = ("file", "time", "t0", "window", "sample_dt", "source_column", "reference_source", "clock_offsets",
                        "debounce_s", "channel", "phases")
const _CHANNEL_KEYS = ("column", "role", "maps_to", "name", "gain", "offset", "one_means", "group", "source", "sigma")
const _ROLES = ("pressure", "mass", "temperature", "valve", "heater")

"Read and check a manifest file. Unknown keys are errors: a misspelt key must not be silently ignored."
function read_manifest(path)
    m = TOML.parsefile(path)
    for k in keys(m)
        k in _MANIFEST_KEYS || error("$path: unknown key \"$k\". Keys: $(join(_MANIFEST_KEYS, ", "))")
    end
    haskey(m, "file") || error("$path: no `file`")
    for c in get(m, "channel", [])
        for k in keys(c)
            k in _CHANNEL_KEYS || error("$path: unknown channel key \"$k\". Keys: $(join(_CHANNEL_KEYS, ", "))")
        end
        haskey(c, "column") || error("$path: a channel has no `column`")
        get(c, "role", "") in _ROLES || error("$path: channel $(c["column"]) needs a role, one of $(join(_ROLES, ", "))")
        c["role"] == "valve" && !(get(c, "one_means", "") in ("open", "closed")) &&
            error("$path: valve channel $(c["column"]) must say one_means = \"open\" or \"closed\"")
        c["role"] in ("valve", "heater") && !haskey(c, "maps_to") && error("$path: channel $(c["column"]) must say which part it maps_to")
    end
    haskey(m, "source_column") && !haskey(m, "reference_source") && error("$path: `source_column` needs a `reference_source`")
    m["file"] = isabspath(m["file"]) ? m["file"] : joinpath(dirname(abspath(path)), m["file"])
    m
end

# Walk the comma-separated fields of `line` once, calling f(k, field) for the wanted column numbers
# (ascending). Splitting every row of a 300 MB log into 60 strings is what makes reading it slow.
function _each_field(f, line, wanted)
    w, col, start = 1, 1, 1
    n = ncodeunits(line)
    while w <= length(wanted)
        stop = something(findnext(',', line, start), n + 1)
        if col == wanted[w]
            f(w, SubString(line, start, prevind(line, stop)))
            w += 1
        end
        stop > n && break
        start, col = stop + 1, col + 1
    end
    nothing
end

"""
    estimate_offsets(path; time, source_column, reference) -> (offsets, drift)

Clock offset of every source that wrote rows to a log, relative to `reference`. Rows are in the
order they arrived, so for each row of another source the last time seen on a reference row is
(to within one reference period) when that row arrived; the median of that minus the row's own
time is the source's offset. `drift` is the difference between that median over the first third
and the last third of the file: near zero if the clocks run at the same rate.
"""
function estimate_offsets(path; time, source_column, reference)
    deltas = Dict{String,Vector{Float64}}()
    wanted, order = Int[], Int[]
    last_ref = NaN
    t, src = Ref(NaN), Ref("")
    for (i, line) in enumerate(eachline(path))
        if i == 1
            h = String.(strip.(split(line, ',')))
            it, is = findfirst(==(time), h), findfirst(==(source_column), h)
            (it === nothing || is === nothing) && error("$path has no column \"$(it === nothing ? time : source_column)\"")
            wanted = sort([it, is]); order = it < is ? [1, 2] : [2, 1]
            continue
        end
        t[] = NaN
        _each_field(line, wanted) do w, field
            order[w] == 1 ? (t[] = something(tryparse(Float64, field), NaN)) : (src[] == field || (src[] = String(field)))
        end
        isnan(t[]) && continue
        if src[] == reference
            last_ref = t[]
        elseif !isnan(last_ref)
            push!(get!(deltas, src[], Float64[]), last_ref - t[])
        end
    end
    med(v) = (s = sort(v); s[(length(s)+1)÷2])
    offsets, drift = Dict(reference => 0.0), Dict(reference => 0.0)
    for (s, d) in deltas
        n = length(d)
        offsets[s] = med(d)
        drift[s] = n >= 30 ? med(d[end-n÷3+1:end]) - med(d[1:n÷3]) : 0.0
    end
    offsets, drift
end

"Remove bounce from a 0/1 series: a change counts only if the new value then lasts at least `hold` seconds."
function debounce(t, v, hold)
    out = similar(v)
    isempty(v) && return out
    state, i, n = v[1], 1, length(v)
    while i <= n
        if v[i] != state
            j = i                                       # the run i..j of the new value
            while j < n && v[j+1] == v[i]
                j += 1
            end
            if j == n || t[j+1] - t[i] >= hold          # it lasts to the end of the data, or long enough
                state = v[i]
            else
                out[i:j] .= state
                i = j + 1
                continue
            end
        end
        out[i] = state
        i += 1
    end
    out
end

# When a phase boundary happens, in seconds after t0: a number, or "<part> open|close[#n]" for the
# n-th time the replayed schedule does that. A part's first action only sets the state it starts
# in, so it is not an event. Without "#n" it is the first such event after `after`.
function _phase_time(spec, schedule; after = -Inf)
    spec isa Real && return Float64(spec)
    m = match(r"^(\S+)\s+(\S+?)(?:#(\d+))?$", strip(String(spec)))
    m === nothing && error("cannot read phase boundary \"$spec\": expected a time or \"PART open\", \"PART close#2\"")
    mine = [a for a in schedule if a.part == m[1]]
    hits = [a.t for a in mine[2:end] if a.action == m[2]]
    if m[3] === nothing
        i = findfirst(>(after), hits)
        i === nothing && error("phase boundary \"$spec\": the log has no such event" * (isfinite(after) ? " after $(after) s" : ""))
        return hits[i]
    end
    n = parse(Int, m[3])
    n <= length(hits) || error("phase boundary \"$spec\": the log has only $(length(hits)) such event(s)")
    hits[n]
end

"""
    load_log(manifest_path) -> Log

Read a log as its manifest describes it: rows from several sources put on one clock, valve and
heater columns debounced, every channel scaled (`gain`, then `offset`) into psi gauge, kg or °C,
channels sharing a `group` averaged, and everything put on an even grid of `sample_dt` seconds
(median of the readings in each interval; the last state for valves and heaters). An interval with
no reading stays NaN: a gap is never filled in.
"""
function load_log(manifest_path)
    m = read_manifest(manifest_path)
    path, warnings = m["file"], String[]
    channels = get(m, "channel", [])
    dt = Float64(get(m, "sample_dt", 0.1))
    multi = haskey(m, "source_column")
    offsets = Dict{String,Float64}()
    if multi
        given = get(m, "clock_offsets", "auto")
        if given == "auto"
            time_name = get(m, "time", first(String.(strip.(split(readline(path), ',')))))
            offsets, drift = estimate_offsets(path; time = time_name, source_column = m["source_column"], reference = m["reference_source"])
            for (s, d) in drift
                abs(d) > 0.5 && push!(warnings, "the clock of $s drifts against $(m["reference_source"]): its offset changes by $(round(d, digits = 2)) s over the log")
            end
        else
            offsets = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in given)
            offsets[m["reference_source"]] = 0.0
        end
    end
    window = haskey(m, "window") ? (Float64(m["window"][1]), Float64(m["window"][2])) : (-Inf, Inf)
    raw_t = [Float64[] for _ in channels]
    raw_v = [Float64[] for _ in channels]
    wanted, what = Int[], Int[]                    # column numbers to read, and what each is: 0 time, -1 source, k channel
    t, off, src = Ref(NaN), Ref(0.0), Ref("")
    for (i, line) in enumerate(eachline(path))
        if i == 1
            h = String.(strip.(split(line, ',')))
            col(name) = something(findfirst(==(name), h), 0)
            it = haskey(m, "time") ? col(m["time"]) : 1
            it == 0 && error("$path has no column \"$(m["time"])\"")
            pairs = Tuple{Int,Int}[(it, 0)]
            if multi
                is = col(m["source_column"])
                is == 0 && error("$path has no column \"$(m["source_column"])\"")
                push!(pairs, (is, -1))
            end
            for (k, c) in enumerate(channels)
                ic = col(c["column"])
                ic == 0 && error("$path has no column \"$(c["column"])\"")
                push!(pairs, (ic, k))
            end
            sort!(pairs)
            allunique(first.(pairs)) || error("$manifest_path: a column is used twice (as time, source or by two channels)")
            wanted, what = first.(pairs), last.(pairs)
            continue
        end
        # time and source can sit after a channel's column, so read those two first
        t[] = NaN
        _each_field(line, wanted) do w, field
            what[w] == 0 ? (t[] = something(tryparse(Float64, field), NaN)) : what[w] == -1 && src[] != field && (src[] = String(field))
        end
        isnan(t[]) && continue
        if multi
            haskey(offsets, src[]) || continue      # a source never seen next to the reference: no clock for it
            off[] = offsets[src[]]
        end
        tr = t[] + off[]
        window[1] <= tr <= window[2] || continue
        _each_field(line, wanted) do w, field
            k = what[w]
            k > 0 || return
            (haskey(channels[k], "source") && src[] != channels[k]["source"]) && return
            v = tryparse(Float64, field)
            v === nothing || (push!(raw_t[k], tr); push!(raw_v[k], v))
        end
    end
    all(isempty, raw_t) && error("$manifest_path: no readings found" * (isfinite(window[1]) ? " in the window $(window)" : ""))
    lo = isfinite(window[1]) ? window[1] : minimum(first(r) for r in raw_t if !isempty(r))
    hi = isfinite(window[2]) ? window[2] : maximum(maximum(r) for r in raw_t if !isempty(r))
    grid = collect(lo:dt:hi)
    hold = Float64(get(m, "debounce_s", 0.0))
    series = Vector{Float64}[]
    for (k, c) in enumerate(channels)
        order = sortperm(raw_t[k])
        ts, vs = raw_t[k][order], raw_v[k][order] .* Float64(get(c, "gain", 1.0)) .+ Float64(get(c, "offset", 0.0))
        state = c["role"] in ("valve", "heater")
        state && hold > 0 && (vs = debounce(ts, vs, hold))
        out = fill(NaN, length(grid))
        a = 1
        for (g, tg) in enumerate(grid)                 # readings in [tg - dt/2, tg + dt/2)
            while a <= length(ts) && ts[a] < tg - dt / 2
                a += 1
            end
            b = a
            while b <= length(ts) && ts[b] < tg + dt / 2
                b += 1
            end
            b > a && (out[g] = state ? vs[b-1] : (s = sort!(vs[a:b-1]); (s[(end+1)÷2] + s[end÷2+1]) / 2))
            a = b
        end
        push!(series, out)
    end
    # Output columns: one per channel, or one per group (the mean of whichever members have a reading).
    names, cols, columns = ["t_s"], Dict("t_s" => grid), [LogColumn("t_s", :time, "", :unknown, false)]
    done = Set{String}()
    for (k, c) in enumerate(channels)
        group = get(c, "group", "")
        name = !isempty(group) ? group : get(c, "name", c["column"])
        name in done && (isempty(group) ? error("$manifest_path: two channels are both called $name") : continue)
        push!(done, name)
        if isempty(group)
            cols[name] = series[k]
        else
            members = [j for (j, d) in enumerate(channels) if get(d, "group", "") == group]
            all(channels[j]["role"] == c["role"] for j in members) || error("$manifest_path: group $group mixes roles")
            cols[name] = [(v = [series[j][g] for j in members if !isnan(series[j][g])]; isempty(v) ? NaN : sum(v) / length(v)) for g in eachindex(grid)]
        end
        role = Symbol(c["role"])
        state = role in (:valve, :heater)
        push!(names, name)
        push!(columns, LogColumn(name, role, state ? c["maps_to"] : "", role === :valve ? Symbol(c["one_means"]) : role === :heater ? :open : :unknown,
                                 true, state ? "" : get(c, "maps_to", ""), Float64(get(c, "sigma", NaN))))
    end
    t0 = Float64(get(m, "t0", grid[1]))
    schedule = log_schedule(cols, columns; t0)
    phases = Dict{String,Tuple{Float64,Float64}}()
    for (name, p) in get(m, "phases", Dict())
        a, b = p isa AbstractVector ? (p[1], p[2]) : (p["from"], p["to"])
        from = _phase_time(a, schedule)
        phases[name] = (from, _phase_time(b, schedule; after = from))
    end
    Log(names, cols, columns, t0, phases, offsets, warnings)
end

"""
    write_log(path, log)

Write a log read through a manifest as a plain CSV on its even grid: the aligned, scaled readings
under their manifest names, and every valve column as 1 = open whatever the original meant by 1.
Missing readings are left blank. The result needs no manifest: `--replay` reads it by name.
"""
function write_log(path, log::Log)
    flip = [c.role === :valve && c.one_means === :closed for c in log.columns]
    open(path, "w") do io
        println(io, join(log.names, ","))
        for i in eachindex(log.cols[log.names[1]])
            for (k, name) in enumerate(log.names)
                v = log.cols[name][i]
                k > 1 && print(io, ",")
                isnan(v) && continue
                flip[k] && (v = 1 - v)
                print(io, log.columns[k].role in (:valve, :heater) ? round(Int, v) : round(v, digits = 3))
            end
            println(io)
        end
    end
end

"""
    suggest_manifest(log_csv, scene) -> String

A draft manifest for a log, from the same name-based guesses the window makes. It is a starting
point to read and correct, not an answer: in particular every valve column's `one_means`.
"""
function suggest_manifest(log_csv, scene::Scene)
    header = String.(strip.(split(readline(log_csv), ',')))
    guesses = guess_log_columns(header, scene)
    io = IOBuffer()
    println(io, "# Draft written by `manifest suggest`. Check every line, above all what a 1 means on each valve column.")
    println(io, "file = ", repr(basename(log_csv)), "\ntime = ", repr(header[1]), "\nsample_dt = 0.1\ndebounce_s = 0.5")
    println(io, "# t0 = 0.0                     # log time that is simulation time zero\n# window = [from, to]          # only this span of log time")
    i = findfirst(h -> lowercase(h) == "source", header)
    if i !== nothing
        println(io, "source_column = ", repr(header[i]), "\nreference_source = \"?\"       # the source whose clock the others are put on\nclock_offsets = \"auto\"")
    end
    for g in guesses
        g.role in (:time, :other) && continue
        println(io, "\n[[channel]]\ncolumn = ", repr(g.name), "\nrole = ", repr(String(g.role)))
        if g.role in (:valve, :heater)
            println(io, "maps_to = ", repr(isempty(g.part) ? "?" : g.part))
            g.role === :valve && println(io, "one_means = ", repr(g.one_means === :unknown ? "?" : String(g.one_means)),
                                         g.one_means === :unknown ? "            # \"open\" or \"closed\": not guessed for this kind of column" : "")
        else
            println(io, "# maps_to = \"TANK.P\"          # the recorded signal to compare it with")
        end
    end
    println(io, "\n# [phases]\n# fill = { from = \"AV201 open\", to = \"AV201 close\" }\n# hold = [170, 1210]            # or seconds after t0")
    String(take!(io))
end
