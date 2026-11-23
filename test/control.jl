# Control policies in studies (src/control.jl). Included from runtests.jl, whose `rig` it uses.

# The rig with a pressure gauge on the tank, whose reading a policy may take.
function gauged_rig(; kw...)
    s = rig(; kw...)
    F.add_part!(s, "Gauge", 8, -4; name = "PT")
    F.connect!(s, ("PT", "a"), ("TANK", "top"))
    s
end

# A controller for a scene from a [controller] table written as TOML.
controller(s, body; kw...) = (sim = F.Sim(s); (sim, F._controller(F.TOML.parse(body), s, sim; kw...)))

@testset "control: schedule!, readings and stand-ins, settings by step" begin
    sim = F.Sim(rig())
    F.schedule!(sim, F.Action(5.0, "VENT", "open", 0.0))
    F.schedule!(sim, F.Action(2.0, "FILL", "open", 0.0))
    F.schedule!(sim, F.Action(5.0, "VENT", "close", 0.0))      # after the one already there for 5 s
    @test [(a.t, a.action) for a in sim.schedule] == [(2.0, "open"), (5.0, "open"), (5.0, "close")]
    F.run_until!(sim, 3.0)
    F.schedule!(sim, F.Action(1.0, "FILL", "close", 0.0))      # in the past: applied at once, and what is done stays done
    @test sim.schedule[1].action == "open" && sim.schedule[2].action == "close" && sim.next_action == 2
    F.run_until!(sim, 3.5)
    @test sim.log[end].part == "FILL" && sim.log[end].action == "close" && abs(sim.log[end].t - 3.0) < 0.01

    s = gauged_rig()
    vent(reading) = """
        [[policy]]
        name = "vent"
        kind = "hysteresis"
        reading = "$reading"
        actuator = "VENT"
        open_above_psig = 500
        close_below = 3.0e6
        """
    _, c = controller(s, vent("PT.P_read"))
    @test c.source["PT.P_read"] == "PT.P"                        # a gauge with no sensor parameter reads the truth
    @test c.policies[1].open_above[""] ≈ 500 * F.PSI + 101325 && c.policies[1].close_below[""] == 3.0e6
    @test_throws ErrorException controller(s, vent("TANK.P"))                     # a true signal is not a reading
    _, c = controller(s, "stand_ins = [\"TANK.P\"]\n" * vent("TANK.P"))           # unless it is listed as a stand-in
    @test c.stand_ins == ["TANK.P"] && c.source["TANK.P"] == "TANK.P"
    @test_throws ErrorException controller(s, "stand_ins = [\"TANK.m\"]\n" * vent("PT.P_read"))   # a stand-in nothing reads
    @test_throws ErrorException controller(s, replace(vent("PT.P_read"), "VENT\"" => "PRV\""))     # not a valve
    @test_throws ErrorException controller(s, vent("PT.P_read") * "open_above = 1\n")             # given twice
    @test_throws ErrorException controller(s, vent("PT.P_read") * "opn_above = 1\n")              # no such key
    @test_throws ErrorException controller(s, "noise = { \"PT.P_read\" = { sigma_psig = 1 } }\n" * vent("PT.P_read"))  # a difference has no gauge zero
    s.parts[F.find_part(s, "PT")].p["offset_psi"] = 3.0
    _, c = controller(s, vent("PT.P_read"))
    @test c.source["PT.P_read"] == "PT.P_read"                   # a gauge that has a sensor model is read through it
    @test F._at(F.PerStep("" => 1.0, "hold" => 2.0), "hold") == 2.0 && F._at(F.PerStep("" => 1.0, "hold" => 2.0), "fill") == 1.0
    @test F._fnv1a("") == 0xcbf29ce484222325 && F._fnv1a("a") == 0xaf63dc4c8601ec8c          # the published test values
    @test F._noise_seed(["B" => 2.0, "A" => 1.0], "x") == F._noise_seed(["A" => 1.0, "B" => 2.0], "x") != F._noise_seed(["A" => 1.0, "B" => 2.0], "y")
end

@testset "control: hysteresis with minimum dwells, pulse and wait, heater" begin
    s = gauged_rig()
    sim, c = controller(s, """
        [[policy]]
        name = "vent"
        kind = "hysteresis"
        reading = "PT.P_read"
        actuator = "VENT"
        open_above = 10
        close_below = 5
        min_open_s = 2
        min_shut_s = 3
        """)
    p = c.policies[1]
    cmds() = [(x.t, x.action) for x in c.log]
    at(t, x) = (c.obs["PT.P_read"] = x; F._tick!(p, c, sim, t))
    at(0.0, 11)                  # above: open
    at(1.0, 4)                   # below, but open for only 1 s of the 2 it must stay
    at(2.0, 4)                   # now it may shut
    at(3.0, 11)                  # above again, but shut for only 1 s of the 3
    at(4.0, 11)
    at(5.0, 11)                  # shut for 3 s: open
    at(6.0, 7)                   # inside the band: nothing
    @test [a for (_, a) in cmds()] == ["open", "close", "open"] && [t for (t, _) in cmds()] == [0.0, 2.0, 5.0]
    @test all(a -> a.part == "VENT" && a.policy == "vent", c.log)
    @test [(a.t, a.action) for a in sim.schedule[sim.next_action:end]] == cmds()     # sent to the run with no delay

    sim, c = controller(s, """
        actuation_delay_s = 0.5
        [[policy]]
        name = "pulse"
        kind = "pulse_wait"
        reading = "PT.P_read"
        actuator = "VENT"
        above = 10
        pulse_s = 1
        wait_s = 4
        [[policy]]
        name = "heat"
        kind = "heater_hold"
        heater = "TANK"
        """)
    p, h = c.policies
    at2(t, x) = (c.obs["PT.P_read"] = x; F._tick!(p, c, sim, t); F._tick!(h, c, sim, t))
    for t in 0.0:1.0:8.0
        at2(t, 20)               # held above: a pulse of 1 s, a wait of 4 s, and again
    end
    @test [(x.t, x.action) for x in c.log if x.policy == "pulse"] == [(0.0, "open"), (1.0, "close"), (5.0, "open"), (6.0, "close")]
    @test [x.action for x in c.log if x.policy == "heat"] == ["on"]              # commanded once, not every tick
    @test sim.schedule[sim.next_action].t == 0.5                                   # each command lands half a second after it is sent
end

@testset "control: a tick with no command does not restart the integrator" begin
    # `iter` counts solver steps since the integrator was last (re)started; `naccept` counts them all.
    s = gauged_rig(; tank_fill = 3.0)
    body = """
        dt = 0.2
        [[policy]]
        name = "vent"
        kind = "hysteresis"
        reading = "PT.P_read"
        actuator = "VENT"
        open_above_psig = OPEN
        close_below_psig = 0
        """
    plain = F.Sim(s)
    F.run_until!(plain, 60.0)
    sim, c = controller(s, replace(body, "OPEN" => "5000"))       # never opens: 300 ticks, no command
    F.run_controlled!(sim, c, 60.0)
    @test isempty(c.log) && sim.integ.iter == sim.integ.stats.naccept                  # never restarted
    @test sim.integ.stats.naccept <= plain.integ.stats.naccept + 3                     # and took the steps it takes unticked
    @test F.signal_at(sim, "TANK.P", 60.0) ≈ F.signal_at(plain, "TANK.P", 60.0) rtol = 1e-6
    @test length(sim.t_hist) >= 301 && all(t -> any(x -> abs(x - t) < 1e-9, sim.t_hist), 0.0:0.2:60.0)    # a sample at every tick
    @test all(k -> length(sim.hist[k]) == length(sim.t_hist), c.signals)
    sim, c = controller(s, replace(body, "OPEN" => "100"))        # the same counter does see a command's restart
    F.run_controlled!(sim, c, 20.0)
    @test length(c.log) == 1 && sim.integ.iter < sim.integ.stats.naccept
    @test sim.hist["ctl.vent.cmd"][1] == 1.0 && sim.hist["ctl.vent.enabled"][end] == 1.0
    @test sim.hist["ctl.vent.obs"][end] ≈ sim.hist["PT.P"][end]
end

@testset "control: in a study (columns, overrides, noise, the Julia escape hatch)" begin
    dir = mktempdir()
    F.save_scene(joinpath(dir, "rig.toml"), gauged_rig())
    file(name, body) = (p = joinpath(dir, name); write(p, body); p)
    table(path) = (lines = readlines(path); (split(lines[1], ','), [split(l, ',') for l in lines[2:end]]))
    body = """
        scene = "rig.toml"
        until = 20
        [[case]]
        name = "fill"
        actions = [{ t = 1.0, part = "FILL", action = "open" }]
        [controller]
        dt = 0.5
        noise = { "PT.P_read" = { sigma_psi = 5 } }
        [[controller.policy]]
        name = "vent"
        kind = "hysteresis"
        reading = "PT.P_read"
        actuator = "VENT"
        open_above_psig = 300
        close_below_psig = 250
        [metrics]
        loaded_kg = { fn = "final", signal = "TANK.m" }
        vent_cmd = { fn = "max", signal = "ctl.vent.cmd" }
        [vary]
        "controller.vent.open_above_psig" = [300, 300.001, 1000]
        "BOTTLE.T_C" = [25]
        """
    st = F.read_study(file("s.toml", body))
    @test length(st.runs) == 3 && st.controller.signals == ["ctl.vent.enabled", "ctl.vent.cmd", "ctl.vent.obs"]
    head, rows = table(F.study(file("s.toml", body); quiet = true))
    col(name) = [r[findfirst(==(name), head)] for r in rows]
    @test all(isempty, col("failed")) && col("stand_ins") == ["", "", ""] && col("escape_hatch") == ["0", "0", "0"]
    byid = Dict(r[1] => r for r in rows)
    cmd(id) = byid[id][findfirst(==("vent_cmd"), head)]
    @test cmd("BOTTLE.T_C=25.0;controller.vent.open_above_psig=1000.0|fill") == "0.0" && cmd("BOTTLE.T_C=25.0;controller.vent.open_above_psig=300.0|fill") == "1.0"
    # the same noise whatever the controller's settings: two candidates a hair apart run the same
    sim1, _, c1 = F._simulate(st, st.runs[1].overrides, st.runs[1].actions; case = "fill")
    sim2, _, c2 = F._simulate(st, st.runs[2].overrides, st.runs[2].actions; case = "fill")
    sim3, _, c3 = F._simulate(st, st.runs[1].overrides, st.runs[1].actions; case = "other")
    @test sim1.hist["ctl.vent.obs"] == sim2.hist["ctl.vent.obs"] && [x.t for x in c1.log] == [x.t for x in c2.log] && !isempty(c1.log)
    @test sim1.hist["ctl.vent.obs"] != sim3.hist["ctl.vent.obs"]                  # another case draws other noise
    noisy = sim1.hist["ctl.vent.obs"] .- [F.signal_at(sim1, "PT.P", t) for t in sim1.t_hist]
    @test 2 < sqrt(sum(abs2, noisy) / length(noisy)) / F.PSI < 8
    # mistakes are found when the file is read
    @test_throws ErrorException F.read_study(file("bad1.toml", replace(body, "\"controller.vent.open_above_psig\"" => "\"controller.vnt.open_above_psig\"")))
    @test_throws ErrorException F.read_study(file("bad2.toml", replace(body, "reading = \"PT.P_read\"" => "reading = \"TANK.P\"")))
    @test_throws ErrorException F.read_study(file("bad3.toml", replace(body, "dt = 0.5" => "tick = 0.5")))
    # a policy written in Julia runs only when allowed, and its runs say so
    file("p.jl", """
        (obs, state) -> begin
            x = obs.readings["PT.P_read"]
            open = get(state, "open", false)
            state["open"] = x > obs.settings["above"]
            state["open"] == open ? [] : [("VENT", state["open"] ? "open" : "close")]
        end
        """)
    jl = """
        scene = "rig.toml"
        until = 20
        [[case]]
        name = "fill"
        actions = [{ t = 1.0, part = "FILL", action = "open" }]
        [controller]
        [[controller.policy]]
        name = "mine"
        kind = "julia"
        file = "p.jl"
        reads = ["PT.P_read"]
        actuators = ["VENT"]
        settings = { above = 2.0e6 }
        """
    @test_throws ErrorException F.read_study(file("j.toml", jl))
    head, rows = table(F.study(file("j.toml", jl); quiet = true, allow_julia_policy = true))
    @test isempty(rows[1][findfirst(==("failed"), head)]) && rows[1][findfirst(==("escape_hatch"), head)] == "1"
    st = F.read_study(file("j.toml", jl); allow_julia_policy = true)
    sim, _, c = F._simulate(st, st.runs[1].overrides, st.runs[1].actions; case = "fill")
    @test !isempty(c.log) && c.log[1].action == "open" && F.signal_at(sim, "PT.P", c.log[1].t) > 2.0e6
end

@testset "control: conditions on known readings" begin
    s = gauged_rig()
    sim, c = controller(s, """
        stand_ins = ["TANK.m"]
        [[policy]]
        name = "ready"
        kind = "band_watch"
        reading = "PT.P_read"
        lo_psig = 100
        hi_psig = 200
        hold_s = 3
        [[step]]
        name = "only"
        until = { any = [{ reading = "PT.P_read", ge_psig = 500 }, { rate = "TANK.m", window_s = 2, lt = 0.01 }, { mark = "ready" }] }
        """)
    cond(d) = F._cond(d, s.settings["P_atm_kPa"] * 1e3, "test", String[], c.mark_names)
    c.step, c.t_start["only"] = "only", 10.0
    obs!(p, m = 0.0) = (c.obs["PT.P_read"] = p * F.PSI + 101325; c.obs["TANK.m"] = m)
    obs!(500); @test F._holds(cond(Dict("reading" => "PT.P_read", "ge_psig" => 500)), c, 11.0)
    @test !F._holds(cond(Dict("reading" => "PT.P_read", "gt_psig" => 500)), c, 11.0)
    @test F._holds(cond(Dict("reading" => "PT.P_read", "lt" => 5e6)), c, 11.0)
    @test !F._holds(cond(Dict("elapsed_ge" => 5)), c, 14.0) && F._holds(cond(Dict("elapsed_ge" => 5)), c, 15.0)    # time in the step
    r = cond(Dict("rate" => "TANK.m", "window_s" => 2, "lt" => 0.01))
    rates = [(obs!(0, 0.2t); F._holds(r, c, t)) for t in 0.0:0.5:3.0]        # 0.2 kg/s for 3 s, then it stalls
    @test rates[1] && !any(rates[2:end])                                         # (the first tick has no rate: 0)
    stalled = [(obs!(0, 0.6); F._holds(r, c, t)) for t in 3.5:0.5:6.0]
    @test stalled == [false, false, false, true, true, true]                    # true once the whole window is flat
    h = cond(Dict("held" => Dict("reading" => "PT.P_read", "ge_psig" => 100), "for_s" => 2))
    @test [(obs!(p); F._holds(h, c, t)) for (t, p) in zip(0:1:5, [150, 150, 150, 50, 150, 150])] == [false, false, true, false, false, false]
    m = cond(Dict("mark" => "ready", "plus_s" => 1))
    @test !F._holds(m, c, 5.0)
    c.marks["ready"] = 4.0
    @test !F._holds(m, c, 4.5) && F._holds(m, c, 5.0)
    obs!(0)
    @test F._holds(cond(Dict("any" => [Dict("elapsed_ge" => 100), Dict("mark" => "ready")])), c, 6.0)
    @test !F._holds(cond(Dict("all" => [Dict("elapsed_ge" => 100), Dict("mark" => "ready")])), c, 6.0)
    @test_throws ErrorException cond(Dict("reading" => "PT.P_read", "ge" => 1, "lt" => 2))      # one comparison each
    @test_throws ErrorException cond(Dict("rate" => "TANK.m", "window_s" => 2, "lt_psig" => 1))  # a rate has no gauge zero
    @test_throws ErrorException cond(Dict("mark" => "nope"))
    # a reading in a condition is held to the same rule as a policy's
    @test_throws ErrorException controller(s, "[[step]]\nname = \"a\"\nuntil = { reading = \"TANK.P\", ge = 1 }\n")
    # settings by step: only steps that exist, and a value for each
    vent = "[[policy]]\nname = \"v\"\nkind = \"hysteresis\"\nreading = \"PT.P_read\"\nactuator = \"VENT\"\nclose_below = 1\n"
    steps = "[[step]]\nname = \"a\"\n[[step]]\nname = \"b\"\n"
    _, c2 = controller(s, vent * "open_above = { a = 2, default = 3 }\n" * steps)
    @test F._at(c2.policies[1].open_above, "a") == 2 && F._at(c2.policies[1].open_above, "b") == 3
    @test_throws ErrorException controller(s, vent * "open_above = { a = 2 }\n" * steps)            # no value for b
    @test_throws ErrorException controller(s, vent * "open_above = { a = 2, c = 3, b = 1 }\n" * steps)   # no step c
    @test_throws ErrorException controller(s, vent * "open_above = { a = 2, b = 1 }\n")              # no steps at all
end

@testset "control: a sequence, marks, step and mark columns, windows from marks, a run ending on a step" begin
    dir = mktempdir()
    F.save_scene(joinpath(dir, "rig.toml"), gauged_rig())
    file(name, body) = (p = joinpath(dir, name); write(p, body); p)
    table(path) = (lines = readlines(path); (split(lines[1], ','), [split(l, ',') for l in lines[2:end]]))
    body = """
        scene = "rig.toml"
        until = { step = "hold", edge = "end", max_s = 200 }
        [controller]
        dt = 0.5
        stand_ins = ["TANK.m"]
        [[controller.step]]
        name = "prep"
        duration_s = 2
        [[controller.step]]
        name = "fill"
        on_enter = [{ part = "FILL", action = "open" }]
        on_exit = [{ part = "FILL", action = "close" }]
        until = { reading = "TANK.m", ge = 1.5 }
        operator_delay_s = 1
        timeout_s = 100
        enable = ["ready"]
        [[controller.step]]
        name = "hold"
        duration_s = 20
        [[controller.policy]]
        name = "ready"
        kind = "band_watch"
        enabled = false
        reading = "PT.P_read"
        lo_psig = LO
        hi_psig = 2000
        hold_s = 3
        [metrics]
        m_fill_end = { fn = "final", signal = "TANK.m", to = { step = "fill", edge = "end" } }
        m_ready_5 = { fn = "final", signal = "TANK.m", to = { mark = "ready", plus_s = 5 } }
        p_after_ready = { fn = "min", signal = "TANK.P", from = { mark = "ready" } }
        m_ready_late = { fn = "final", signal = "TANK.m", to = { mark = "ready", plus_s = 1000 } }
        [bands]
        signals = ["TANK.P"]
        align = "ready"
        """
    out = F.study(file("seq.toml", replace(body, "LO" => "100")); quiet = true)
    head, rows = table(out)
    val(name) = (x = rows[1][findfirst(==(name), head)]; isempty(x) ? NaN : parse(Float64, x))
    @test isempty(rows[1][findfirst(==("failed"), head)]) && rows[1][findfirst(==("stand_ins"), head)] == "TANK.m"
    @test val("t_prep.start") == 0 && val("t_prep.end") == 2 && val("t_fill.start") == 2
    @test val("t_hold.start") == val("t_fill.end") && val("t_hold.end") == val("t_fill.end") + 20 == val("t_end")   # the run ends with hold
    @test val("fill.timed_out") == 0 && val("m_fill_end") > 1.5                   # over the line by what came in during the operator's second
    @test val("t_fill.start") < val("t_ready") <= val("t_fill.start") + 30         # watched only from the fill on
    @test val("m_ready_5") > 0 && val("p_after_ready") > 100 * F.PSI && isnan(val("m_ready_late"))   # a window past the run's end has no value
    b = readlines(joinpath(dir, "seq_results", "bands.csv"))
    @test startswith(b[1], "t_s,n,TANK.P.p5") && length(b) - 1 == floor(Int, val("t_end") - val("t_ready")) + 1   # time from the mark
    # a mark that never happens: its columns and its windows are empty, and the run has not failed
    head, rows = table(F.study(file("never.toml", replace(body, "LO" => "1900")); quiet = true))
    @test isempty(rows[1][findfirst(==("failed"), head)]) && isnan(val("t_ready")) && isnan(val("m_ready_5")) && isnan(val("p_after_ready"))
    @test !isnan(val("m_fill_end")) && length(readlines(joinpath(dir, "never_results", "bands.csv"))) == 1
    # a step that times out is marked, and the sequence goes on
    head, rows = table(F.study(file("slow.toml", replace(replace(body, "LO" => "100"), "timeout_s = 100" => "timeout_s = 1")); quiet = true))
    @test val("fill.timed_out") == 1 && val("t_fill.end") == 3 && val("m_fill_end") < 1.5
    # mistakes in events are found when the file is read
    body = replace(body, "LO" => "100")
    @test_throws ErrorException F.read_study(file("bad1.toml", replace(body, "to = { mark = \"ready\", plus_s = 5 }" => "to = { mark = \"redy\" }")))
    @test_throws ErrorException F.read_study(file("bad2.toml", replace(body, "step = \"hold\", edge = \"end\", max_s" => "step = \"hld\", edge = \"end\", max_s")))
    @test_throws ErrorException F.read_study(file("bad3.toml", replace(body, "align = \"ready\"" => "align = \"fill\"")))
end
