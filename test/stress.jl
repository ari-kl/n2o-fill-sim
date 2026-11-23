# Robustness and worst-case speed on scenes nobody tuned.   julia --project test/stress.jl [first last] [-v] [-d] [-p]
#
# 1. A bottle left in the sun until it is supercritical, opened into an empty tank.
# 2. Random plumbing built from default parts, with valves opened and toggled. A quarter of the
#    valves leak when shut (1e-6 to 1e-2 of their area) and three in ten are set part-open (0.001
#    to 1) instead of opened: trickles are what the junction solve handles worst. -p leaves both
#    out, which gives the plumbings as they were before valves could do either.
# Both are driven the way the window drives the solver (one step at a time inside a 12 ms budget),
# so "worst frame" is how long the window would have frozen.
# Because a frame ends on the clock, two runs of the same code do not take the same steps, and the
# count of failed solves moves by a fifth from run to run. To compare two versions of the solver,
# use -d: a frame then ends after 40 solver steps instead, and the failed solves, the time reached and
# the work counts (the costliest step among them) repeat exactly (the times, and "flagged", which
# counts slow frames, still do not).
using N2OFillSim, Random
const F = N2OFillSim
const DETERMINISTIC = "-d" in ARGS
const PLAIN = "-p" in ARGS

"Run `n` window frames at `speed`; returns (wall ms, worst frame ms, frames over 50 ms, part flows in the costliest solver step)."
function frames!(sim, n; speed = 20.0, each = i -> nothing, until = Inf)
    wall = worst = 0.0
    slow = costliest = 0
    for i in 1:n
        (isempty(sim.failed) && sim.t < until) || break
        each(i)
        t0 = time_ns(); done = 0.0; want = speed / 60; steps = 0
        while done < want - 1e-9 && (DETERMINISTIC ? steps < 40 : (time_ns() - t0) < 12e6) && isempty(sim.failed)
            before = F.work(sim)
            done += F.advance!(sim, want - done; max_steps = 1)
            steps += 1; costliest = max(costliest, F.work(sim) - before)
        end
        dt = (time_ns() - t0) / 1e6
        wall += dt; worst = max(worst, dt); dt > 50 && (slow += 1)
    end
    wall, worst, slow, costliest
end

function hot_bottle()
    s = F.Scene()
    F.add_part!(s, "Bottle", 0, 0; name = "B").p["T_C"] = 62.0
    for (i, (kind, name)) in enumerate((("Orifice", "OR"), ("Pipe", "L1"), ("Valve", "V"), ("Pipe", "L2"), ("Tank", "T")))
        F.add_part!(s, kind, 3i, 0; name)
    end
    for (a, b) in ((("B", "out"), ("OR", "a")), (("OR", "b"), ("L1", "a")), (("L1", "b"), ("V", "a")), (("V", "b"), ("L2", "a")), (("L2", "b"), ("T", "bottom")))
        F.connect!(s, a, b)
    end
    s
end

function random_scene(rng)
    s = F.Scene(); n = Ref(0)
    add(kind) = (n[] += 1; F.add_part!(s, kind, 3 * n[], 0; name = kind * string(n[])).name)
    b, t = add("Bottle"), add("Tank")
    chain(from, kinds, to) = begin
        prev = from
        for k in kinds
            p = add(k); F.connect!(s, prev, (p, "a")); prev = (p, "b")
        end
        F.connect!(s, prev, to)
    end
    inline = ["Valve", "Valve", "Pipe", "Orifice", "CheckValve", "Regulator"]
    fill = shuffle(rng, [rand(rng, inline) for _ in 1:rand(rng, 1:4)])
    any(==("Valve"), fill) || rand(rng) < 0.7 && push!(fill, "Valve")
    chain((b, "out"), fill, (t, rand(rng, ["bottom", "bottom", "top"])))
    if rand(rng) < 0.8
        a = add("Atmosphere")
        chain((t, rand(rng, ["vent", "top"])), [rand(rng, ["Valve", "Valve", "Orifice", "Pipe", "ReliefValve", "CheckValve"]) for _ in 1:rand(rng, 1:3)], (a, "a"))
    end
    if rand(rng) < 0.4                                          # a tee off the fill line to a dump
        a = add("Atmosphere"); w = rand(rng, s.wires)
        chain(w.a, [rand(rng, ["Valve", "ReliefValve", "Orifice"])], (a, "a"))
    end
    rand(rng) < 0.3 && F.connect!(s, (add("Gauge"), "a"), rand(rng, s.wires).a)
    s
end

frames!(F.Sim(hot_bottle()), 5)                                 # compile before timing
sim = F.Sim(hot_bottle()); F.act!(sim, "V", "open"); F.reset_tally!()
wall, worst, slow, costliest = frames!(sim, 300; speed = 1e6)
println("hot bottle, 300 frames at max:  wall ", round(Int, wall), " ms  worst frame ", round(Int, worst), " ms  frames over 50 ms ", slow, "  failed solves ", sim.net.unconverged,
        "  costliest step ", costliest, " part flows  reached t=", round(Int, sim.t), " s  (eos ", F.TALLY.eos, ", choke searches ", F.TALLY.choke, ")")

function search(first_seed, last_seed; verbose = false)
    flagged = over_1s = failed_solves = costliest_all = 0
    worst_all = total = 0.0
    for seed in first_seed:last_seed
        rng = MersenneTwister(seed)
        s = random_scene(rng)
        isempty(F.check(s)[1]) || continue
        # Leaks and part openings draw from a stream of their own, so that a seed's plumbing and
        # its toggles are the same with and without them.
        rng2 = MersenneTwister(seed + 1_000_000)
        trickle() = !PLAIN && rand(rng2) < 0.3
        for p in s.parts
            p.kind == "Valve" && !PLAIN && rand(rng2) < 0.25 && (p.p["leak_frac"] = 10.0^rand(rng2, -6:-2))
        end
        sim = F.Sim(s)
        valves = [e.name for e in sim.net.edges if e.kind === :valve]
        wall, worst, slow, costliest = frames!(sim, 2000; each = i -> begin
            i == 30 && foreach(v -> trickle() ? F.act!(sim, v, "set", 10.0^(-3 * rand(rng2))) : F.act!(sim, v, "open"), valves)
            i % 600 == 60 && !isempty(valves) && F.act!(sim, rand(rng, valves), "toggle")
        end, until = 240.0)
        bad = slow > 2 || worst > 400 || !isempty(sim.failed) || sim.t < 240
        flagged += bad; over_1s += worst > 1000; failed_solves += sim.net.unconverged
        worst_all = max(worst_all, worst); total += wall; costliest_all = max(costliest_all, costliest)
        bad && verbose && println("seed ", seed, "  t=", round(sim.t, digits = 1), "  worst ", round(Int, worst), " ms  slow frames ", slow, "  failed solves ", sim.net.unconverged, "  ", join([p.kind for p in s.parts], " "))
    end
    println("random scenes ", first_seed, "-", last_seed, ":  flagged ", flagged, "  with a frame over 1 s ", over_1s, "  worst frame ", round(Int, worst_all), " ms  failed solves ", failed_solves, "  costliest step ", costliest_all, " part flows  total wall ", round(total / 1e3, digits = 1), " s")
end
seeds = filter(a -> !startswith(a, "-"), ARGS)
search((length(seeds) == 2 ? parse.(Int, seeds) : (1, 150))...; verbose = "-v" in ARGS)
