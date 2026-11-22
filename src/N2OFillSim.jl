"""
N2OFillSim: a lumped thermo-fluid simulator for filling a hybrid rocket's nitrous oxide run tank.
See README.md for use and docs/MODEL.md for the equations and their sources.
"""
module N2OFillSim
include("n2o.jl")
include("flow.jl")
include("scene.jl")
include("vessel.jl")
include("network.jl")
include("logs.jl")
include("gui.jl")

const USAGE = """
usage: julia --project bin/n2o-fill-sim.jl [scene.toml] [options]
  (no options)            open the editor on the scene (or an empty one)
  --run                   start in run mode
  --replay log.csv        drive the run from a test log (valve and heater columns) and overlay its readings
  --manifest log.toml     the same, for a log described by a manifest (see docs/LOGS.md)
  --t0 SECONDS            log time that corresponds to simulation time zero (default: first row)
  --set NAME=VALUE        change a parameter for this run: PART.param=value or settings.key=value (repeatable)
  --schedule file.csv     add scripted actions from a file with columns t_s,part,action,value
  --headless              no window: run the schedule / replay and print a summary
  --until SECONDS         how long to run headless (default: the end of the log when replaying,
                          otherwise the last action + 60 s)
  --sample-dt SECONDS     headless: shortest interval between recorded samples (default 0.1); samples
                          fall on the solver's steps, so they are further apart when little is changing
  --csv out.csv           headless: write the recorded signals, and the action log as out_events.csv
  --columns A,B,...       with --csv: write only these signals (for example TANK.P,TANK.m)

       julia --project bin/n2o-fill-sim.jl manifest suggest log.csv scene.toml
  prints a draft manifest for a log, to be checked and corrected by hand
       julia --project bin/n2o-fill-sim.jl manifest extract log.toml out.csv
  writes the log as its manifest describes it: one clock, an even grid, valve columns as 1 = open
"""

const _VALUED = ("--replay", "--manifest", "--t0", "--set", "--schedule", "--until", "--sample-dt", "--csv", "--columns", "--shot", "--frames")
const _FLAGS = ("--help", "--run", "--headless")

"Split a command line into positional arguments and options (an option may be given more than once)."
function options(args)
    files, opts, i = String[], Dict{String,Vector{String}}(), 1
    while i <= length(args)
        a = args[i]
        if a in _VALUED
            i < length(args) || error("$a needs a value")
            push!(get!(opts, a, String[]), args[i+1])
            i += 1
        elseif a in _FLAGS
            get!(opts, a, String[])
        elseif startswith(a, "--")
            error("unknown option $a (see --help)")
        else
            push!(files, a)
        end
        i += 1
    end
    files, opts
end

"""
What a finished (or failed) run comes to: how far it got, why it stopped if it failed, how many
junction solves did not settle, the mass and energy balance errors, what was vented, wall time.
"""
function run_summary(sim::Sim, wall_s = NaN)
    m, E = closure(sim)
    (t_s = simtime(sim), failed = sim.failed, unconverged = sim.net.unconverged, mass_error = m, energy_error = E,
     vented_kg = sim.u[end-2], wall_s = wall_s)
end

function main(args = ARGS)
    files, opts = options(args)
    opt(name, default = nothing) = haskey(opts, name) ? opts[name][end] : default
    flag(name) = haskey(opts, name)
    flag("--help") && return print(USAGE)
    if !isempty(files) && files[1] == "manifest"
        if length(files) == 4 && files[2] == "extract"
            log = load_log(files[3])
            foreach(w -> println("warning: ", w), log.warnings)
            write_log(files[4], log)
            return println("wrote ", files[4], ": ", length(log.cols[log.names[1]]), " rows, columns ", join(log.names, ", "))
        end
        (length(files) == 4 && files[2] == "suggest") || error("usage: manifest suggest log.csv scene.toml   or   manifest extract log.toml out.csv")
        return print(suggest_manifest(files[3], load_scene(files[4])))
    end
    path = isempty(files) ? "" : files[1]
    isempty(path) || isfile(path) || flag("--headless") && error("no such scene file: $path")
    scene = isempty(path) || !isfile(path) ? Scene() : load_scene(path)
    apply_overrides!(scene, get(opts, "--set", String[]))
    overlay, log_end = nothing, nothing
    if opt("--replay") !== nothing || opt("--manifest") !== nothing
        if opt("--manifest") !== nothing
            log = load_log(opt("--manifest"))
            foreach(w -> println("warning: ", w), log.warnings)
            names, cols, columns, first_t = log.names, log.cols, log.columns, log.t0
        else
            names, cols = read_csv(opt("--replay"); min_dt = 0.05)
            columns, first_t = guess_log_columns(names, scene), cols[names[1]][1]
        end
        t = cols[names[1]]
        t0 = parse(Float64, opt("--t0", string(first_t)))
        scene.schedule = log_schedule(cols, columns; t0)
        unsure = [c.name for c in columns if c.role === :valve && c.one_means === :unknown]
        isempty(unsure) || println("not replayed, because whether 1 means open or closed is not known: ", join(unsure, ", "), "  (load the log in the window to choose)")
        overlay = (t .- t0, log_overlay(cols, columns, scene.settings["P_atm_kPa"] * 1e3))
        log_end = t[end] - t0
    end
    for file in get(opts, "--schedule", String[])
        append!(scene.schedule, read_events(file))
    end
    sort!(scene.schedule; by = a -> a.t)
    if flag("--headless")
        errors, warnings = check(scene; unattended = true)
        foreach(w -> println("warning: ", w), warnings)
        isempty(errors) || error("cannot run: " * join(errors, "; "))
        sim = Sim(scene; sample_dt = parse(Float64, opt("--sample-dt", "0.1")))
        until = parse(Float64, opt("--until", string(something(log_end, (isempty(scene.schedule) ? 0.0 : scene.schedule[end].t) + 60))))
        wall = @elapsed run_until!(sim, until)
        r = run_summary(sim, wall)
        println("simulated $(round(r.t_s, digits = 1)) s in $(round(wall, digits = 1)) s of wall time")
        isempty(r.failed) || println("FAILED: ", r.failed)
        for n in sim.net.nodes
            n.kind === :vessel || continue
            v = node_report(sim.net, n, sim.u)
            println(rpad(n.name, 10), round(Int, (v.P - sim.net.env.P_atm) / PSI), " psig  ", round(v.T - 273.15, digits = 1), " °C  ",
                    round(v.m, digits = 2), " kg  liquid ", round(Int, 100 * v.liquid_frac), " % of volume, ", round(Int, 100 * (1 - v.x)), " % of mass")
        end
        println("vented ", round(r.vented_kg, digits = 3), " kg; mass/energy balance error ", round.((r.mass_error, r.energy_error), sigdigits = 2),
                r.unconverged > 0 ? "; junction solves that did not settle: $(r.unconverged)" : "")
        if opt("--csv") !== nothing
            columns = opt("--columns") === nothing ? nothing : String.(strip.(split(opt("--columns"), ',')))
            write_csv(opt("--csv"), sim; columns)
            write_events(events_path(opt("--csv")), sim)
        end
        return r
    end
    app = App(scene, path)
    overlay === nothing || set_overlay!(app, overlay...)
    flag("--run") && start_run!(app)
    gui(app; shot = opt("--shot", ""), frames = parse(Int, opt("--frames", "0")))
end

using PrecompileTools: @compile_workload
@compile_workload begin
    s = Scene()
    add_part!(s, "Bottle", 0, 0; name = "B")
    add_part!(s, "Pipe", 3, 0; name = "L")
    add_part!(s, "Valve", 6, 0; name = "V")
    add_part!(s, "Tank", 9, 0; name = "T").p["heater_W"] = 100.0
    add_part!(s, "Valve", 12, 0; name = "W")
    add_part!(s, "Orifice", 14, 0; name = "O")
    add_part!(s, "Atmosphere", 16, 0; name = "A")
    for (a, b) in ((("B", "out"), ("L", "a")), (("L", "b"), ("V", "a")), (("V", "b"), ("T", "bottom")),
                   (("T", "vent"), ("W", "a")), (("W", "b"), ("O", "a")), (("O", "b"), ("A", "a")))
        connect!(s, a, b)
    end
    check(s)
    sim = Sim(s)
    act!(sim, "V", "open"); act!(sim, "W", "open"); act!(sim, "T", "on")
    for _ in 1:20
        advance!(sim, 0.1; max_steps = 2)
    end
    closure(sim)
    for n in sim.net.nodes
        node_report(sim.net, n, sim.u)
    end
    # the file and log handling a headless replay goes through before it reaches the solver
    mktempdir() do dir
        file, log = joinpath(dir, "s.toml"), joinpath(dir, "log.csv")
        save_scene(file, s)
        s2 = load_scene(file)
        apply_overrides!(s2, ["T.fill_kg=0", "settings.T_amb_C=20"])
        write(log, "t_s,PT202_psig,LCFill_kg,V,W\n0.0,1.0,0.0,0,0\n1.0,2.0,0.1,1,0\n2.0,3.0,0.2,1,1\n")
        names, cols = read_csv(log; min_dt = 0.05)
        columns = guess_log_columns(names, s2)
        log_schedule(cols, columns; t0 = 0.0)
        log_overlay(cols, columns, 101325.0)
        run_summary(sim, 0.0)
        write_csv(joinpath(dir, "out.csv"), sim)
        write_events(events_path(joinpath(dir, "out.csv")), sim)
    end
end
# The window code cannot be run here (no display while precompiling), so ask for it to be compiled
# from its signatures instead. Without this the window takes about 15 s to appear.
for f in (frame, draw_canvas, side_panel, plot_panel, shortcuts, build_tab, tanks_tab, signals_tab, view_tab, schedule_tab,
          run_header, browser_popup, log_popup, unsaved_popup, start_run!, survey!, zoom_fit!)
    precompile(f, (App,))
end
precompile(gui, (App,))
precompile(Core.kwcall, (@NamedTuple{shot::String, frames::Int64}, typeof(gui), App))
precompile(main, (Vector{String},))
precompile(draw_part, (Ptr{ig.lib.ImDrawList}, App, Part))
precompile(canvas_edit_clicks, (App, Float32, Float32))
precompile(canvas_run_clicks, (App, Float32, Float32))
precompile(describe, (App, Part))

end
