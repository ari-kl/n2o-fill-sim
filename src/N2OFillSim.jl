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
  --t0 SECONDS            log time that corresponds to simulation time zero (default: first row)
  --headless              no window: run the schedule / replay and print a summary
  --until SECONDS         how long to run headless (default: end of the schedule + 60 s)
  --csv out.csv           headless: write the recorded signals, and the action log as out_events.csv
"""

function main(args = ARGS)
    opt(name, default = nothing) = (i = findfirst(==(name), args); i === nothing ? default : args[i+1])
    flag(name) = name in args
    flag("--help") && return print(USAGE)
    files = [a for (i, a) in enumerate(args) if !startswith(a, "--") && (i == 1 || !(args[i-1] in ("--replay", "--t0", "--until", "--csv", "--shot", "--frames")))]
    path = isempty(files) ? "" : files[1]
    scene = isempty(path) || !isfile(path) ? Scene() : load_scene(path)
    app = App(scene, path)
    if opt("--replay") !== nothing
        names, cols = read_csv(opt("--replay"); min_dt = 0.05)
        columns = guess_log_columns(names, scene)
        t = cols[names[1]]
        t0 = parse(Float64, opt("--t0", string(t[1])))
        scene.schedule = log_schedule(cols, columns; t0)
        unsure = [c.name for c in columns if c.role === :valve && c.one_means === :unknown]
        isempty(unsure) || println("not replayed, because whether 1 means open or closed is not known: ", join(unsure, ", "), "  (load the log in the window to choose)")
        set_overlay!(app, t .- t0, log_overlay(cols, columns, scene.settings["P_atm_kPa"] * 1e3))
    end
    if flag("--headless")
        sim = Sim(scene)
        until = parse(Float64, opt("--until", string((isempty(scene.schedule) ? 0.0 : scene.schedule[end].t) + 60)))
        wall = @elapsed run_until!(sim, until)
        println("simulated $(round(simtime(sim), digits = 1)) s in $(round(wall, digits = 1)) s of wall time")
        isempty(sim.failed) || println("FAILED: ", sim.failed)
        for n in sim.net.nodes
            n.kind === :vessel || continue
            r = node_report(sim.net, n, sim.u)
            println(rpad(n.name, 10), round(Int, (r.P - sim.net.env.P_atm) / PSI), " psig  ", round(r.T - 273.15, digits = 1), " °C  ",
                    round(r.m, digits = 2), " kg  liquid ", round(Int, 100 * r.liquid_frac), " % of volume, ", round(Int, 100 * (1 - r.x)), " % of mass")
        end
        println("vented ", round(sim.u[end-2], digits = 3), " kg; mass/energy balance error ", round.(closure(sim), sigdigits = 2))
        if opt("--csv") !== nothing
            write_csv(opt("--csv"), sim)
            write_events(events_path(opt("--csv")), sim)
        end
        return
    end
    flag("--run") && start_run!(app)
    gui(app; shot = opt("--shot", ""), frames = parse(Int, opt("--frames", "0")))
end

# Compile the simulation path into the package image by running a small fill once at precompile
# time. Without this the first press of Run freezes for several seconds while Julia compiles the
# solver.
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
