# Speed numbers for n2o-fill-sim.   julia --project test/bench.jl [--gui]
# Without --gui it needs no display. With --gui it opens the window for about 15 s.
const T_PROCESS = time()
load_s = @elapsed @eval using N2OFillSim
const F = N2OFillSim
using Statistics: median, quantile

const ROOT = dirname(@__DIR__)
function hf4()
    _, cols = F.read_csv(joinpath(ROOT, "test", "data", "hf4_fill_10hz.csv"))
    s = F.load_scene(joinpath(ROOT, "scenes", "hf4.toml"))
    s.schedule = F.schedule_from_log(cols["t_s"], cols, Dict("AV201" => ("AV201", "open", "close"), "AV202" => ("AV202", "open", "close"),
        "AV203" => ("AV203", "open", "close"), "AV204" => ("AV204", "open", "close"), "heater" => ("TANK", "on", "off")); t0 = 9020.0)
    s, cols
end
ms(x) = round(1e3x, digits = 2)

println("package load                      ", round(load_s, digits = 2), " s")
scene, cols = hf4()
build = @elapsed (sim = F.Sim(scene))
first_step = @elapsed F.advance!(sim, 1 / 60; max_steps = 1)
println("build Sim / first advance         ", ms(build), " ms / ", ms(first_step), " ms")

# Whole replay in one go. Wall time moves by 20 % from run to run, so take the best of five and
# print the work counts, which do not move at all.
replay = Inf
for _ in 1:5
    global sim = F.Sim(scene)
    F.reset_tally!()
    global replay = min(replay, @elapsed F.run_until!(sim, 990.0))
end
c = F.TALLY
println("HF4 replay, 990 s simulated       ", round(replay, digits = 3), " s best of 5   (", sim.integ.stats.naccept, " steps, ", sim.net.unconverged, " unconverged junction solves)")
println("  work: rhs ", c.rhs, "  ODE Jacobians ", c.jac, "  junction Jacobians ", c.newton, "  part flows ", c.edge_flow, "  choke searches ", c.choke, "  eos ", c.eos)
println("  end state: tank ", round((sim.hist["TANK.P"][end] - 101325) / F.PSI, digits = 3), " psig, ", round(sim.hist["TANK.m"][end], digits = 5), " kg")

# One solver step at a time: how bad is the worst one? This is what bounds a frame in the window,
# which never asks for more than one step once its 12 ms are spent.
sim = F.Sim(scene)
steps = Float64[]
while F.simtime(sim) < 990.0 && isempty(sim.failed)
    push!(steps, @elapsed F.advance!(sim, 990.0 - F.simtime(sim); max_steps = 1))
end
order = sortperm(steps; rev = true)
println("single solver steps               median ", ms(median(steps)), " ms, 99th pct ", ms(quantile(steps, 0.99)), " ms, worst ", ms(steps[order[1]]), " ms, next ", ms(steps[order[2]]), " ms")

# Driven exactly the way the window drives it at x1: each frame wants 1/60 s and stops asking
# once 12 ms are spent.
function frame_times(sim, seconds; speed = 1.0)
    out = Float64[]
    for _ in 1:round(Int, 60 * seconds)
        t0 = time_ns(); done = 0.0; want = speed / 60
        while done < want - 1e-9 && (time_ns() - t0) < 12e6 && isempty(sim.failed)
            done += F.advance!(sim, want - done; max_steps = 1)
        end
        push!(out, (time_ns() - t0) / 1e9)
    end
    out
end
frames = frame_times(F.Sim(scene), 330)
println("solver per frame at x1 (0-330 s)  median ", ms(median(frames)), " ms, 99th pct ", ms(quantile(frames, 0.99)), " ms, worst ", ms(maximum(frames)), " ms, frames over 15 ms: ", count(>(0.015), frames))
sim = F.Sim(scene); F.run_until!(sim, 990.0)
hold = @elapsed for _ in 1:5000; F.advance!(sim, 0.1); end
println("recording a quiet hold            ", round(1e6 * hold / 5000, digits = 1), " µs per 0.1 s sample (", round(Int, 0.1 * 5000 / hold), "x real time at most)")

u = copy(sim.u); du = similar(u); J = zeros(length(u), length(u))
F.rhs!(du, u, sim.net, 0.0)
t_rhs = @elapsed for _ in 1:500; F.rhs!(du, u, sim.net, 0.0); end
t_jac = @elapsed for _ in 1:50; F.jac!(J, u, sim.net, 0.0); end
t_refresh = @elapsed for _ in 1:500; F.refresh!(sim); end
a_rhs = @allocated F.rhs!(du, u, sim.net, 0.0)
a_jac = @allocated F.jac!(J, u, sim.net, 0.0)
println("rhs! / jac! / refresh!            ", round(1e6 * t_rhs / 500, digits = 1), " µs / ", round(1e6 * t_jac / 50, digits = 1), " µs / ",
        round(1e6 * t_refresh / 500, digits = 1), " µs;  allocations ", a_rhs, " / ", a_jac, " bytes")

if "--gui" in ARGS
    import CImGui as ig
    scene, cols = hf4()
    app = F.App(scene, joinpath(ROOT, "scenes", "hf4.toml"))
    names, _ = F.read_csv(joinpath(ROOT, "test", "data", "hf4_fill_10hz.csv"))
    F.set_overlay!(app, cols["t_s"] .- 9020.0, F.log_overlay(cols, F.guess_log_columns(names, scene), 101325.0))
    rec = Dict{String,Vector{NTuple{3,Float64}}}()
    note(name) = push!(get!(rec, name, NTuple{3,Float64}[]), (1e3 * unsafe_load(ig.GetIO().DeltaTime), app.frame_ms, app.solver_ms))
    first_frame = Ref(0.0)
    script(app, n) = begin
        n == 2 && (first_frame[] = time() - T_PROCESS)
        n == 10 && (note("warm-up"); delete!(rec, "warm-up"))     # compile the measuring itself before measuring
        20 < n <= 140 && note("editing")
        n == 141 && F.start_run!(app)
        150 < n <= 270 && note("run mode, paused")
        n == 271 && (app.paused = false; app.speed = 100.0)
        if n > 271 && app.speed == 100.0 && F.simtime(app.sim) >= 61.0
            app.speed = 1.0
        end
        n > 271 && app.speed == 1.0 && note(F.simtime(app.sim) < 63.0 ? "x1, idle" : F.simtime(app.sim) < 65.5 ? "x1, fill starting" : "x1, filling")
    end
    F.gui(app; shot = joinpath(tempdir(), "n2o-fill-sim_bench.ppm"), frames = 271 + 120 + 60 * 9, script)
    println("time from process start to the first frame   ", round(first_frame[], digits = 2), " s")
    println("                          frame ms (median / 95th / worst) | our code ms (median / worst) | solver ms (median / worst)")
    for k in ("editing", "run mode, paused", "x1, idle", "x1, fill starting", "x1, filling")
        haskey(rec, k) || continue
        f, o, s = first.(rec[k]), getindex.(rec[k], 2), last.(rec[k])
        println(rpad(k, 22), lpad(length(f), 5), " frames  ", round(median(f), digits = 1), " / ", round(quantile(f, 0.95), digits = 1), " / ", round(maximum(f), digits = 1),
                "  |  ", round(median(o), digits = 2), " / ", round(maximum(o), digits = 1), "  |  ", round(median(s), digits = 2), " / ", round(maximum(s), digits = 1))
    end
end
