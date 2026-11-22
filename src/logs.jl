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
    role::Symbol                  # :time, :pressure (psi gauge), :mass (kg), :valve, :heater, :other
    part::String                  # the scene part a :valve or :heater column drives ("" if none)
    one_means::Symbol             # valves: :open, :closed, or :unknown (must be chosen before replay)
    use::Bool                     # overlay it (pressure, mass) / replay it (valve, heater)
end

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

"""
The readings a log offers for plotting next to the simulation: (label, kind, values in SI) for
each pressure and mass column marked `use`. Logged pressures are psi gauge.
"""
function log_overlay(cols, columns::Vector{LogColumn}, p_atm)
    out = Tuple{String,Symbol,Vector{Float64}}[]
    for c in columns
        c.use || continue
        c.role === :pressure && push!(out, (c.name, :pressure, cols[c.name] .* PSI .+ p_atm))
        c.role === :mass && push!(out, (c.name, :mass, copy(cols[c.name])))
    end
    out
end
