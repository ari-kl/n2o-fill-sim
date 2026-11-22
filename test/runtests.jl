using Test, N2OFillSim, CoolProp
const F = N2OFillSim
const N2O = "NitrousOxide"
psia(psig) = (psig + 14.4) * F.PSI          # docs/MODEL.md section 9 uses 14.4 psi for the gauge offset

# A scene built in code: bottle → fill valve → tank, tank vent → valve → orifice → atmosphere,
# relief valve on the tank top.
function rig(; bottle_T = 25.0, tank_fill = 0.0, tank_T = 25.0, vent_mm = 2.0, air = 1.0)
    s = F.Scene()
    F.add_part!(s, "Bottle", 0, 0; name = "BOTTLE").p["T_C"] = bottle_T
    F.add_part!(s, "Valve", 4, 0; name = "FILL")
    t = F.add_part!(s, "Tank", 8, 0; name = "TANK")
    t.p["fill_kg"], t.p["T_C"], t.p["air"] = tank_fill, tank_T, air
    F.add_part!(s, "Valve", 12, 0; name = "VENT")
    F.add_part!(s, "Orifice", 14, 0; name = "VENT_OR").p["d_mm"] = vent_mm
    F.add_part!(s, "ReliefValve", 12, -4; name = "PRV")
    F.add_part!(s, "Atmosphere", 16, 0; name = "ATM")
    F.add_part!(s, "Atmosphere", 16, -4; name = "ATM2")
    for (a, b) in ((("BOTTLE", "out"), ("FILL", "a")), (("FILL", "b"), ("TANK", "bottom")), (("TANK", "vent"), ("VENT", "a")),
                   (("VENT", "b"), ("VENT_OR", "a")), (("VENT_OR", "b"), ("ATM", "a")), (("TANK", "top"), ("PRV", "a")), (("PRV", "b"), ("ATM2", "a")))
        F.connect!(s, a, b)
    end
    s
end
tank(sim) = F.node_report(sim.net, sim.net.nodes[findfirst(n -> n.name == "TANK", sim.net.nodes)], sim.u)
bottle(sim) = F.node_report(sim.net, sim.net.nodes[findfirst(n -> n.name == "BOTTLE", sim.net.nodes)], sim.u)
closes(sim) = all(x -> abs(x) < 2e-4, F.closure(sim))

@testset "N2OFillSim" begin

@testset "properties: saturation line vs CoolProp" begin
    for T in 185.0:4.3:309.0
        st = F.sat_T(T)
        @test st.P ≈ PropsSI("P", "T", T, "Q", 0, N2O) rtol = 1e-5
        @test st.ρl ≈ PropsSI("D", "T", T, "Q", 0, N2O) rtol = 1e-5
        @test st.ρv ≈ PropsSI("D", "T", T, "Q", 1, N2O) rtol = 1e-5
        @test st.hv - st.hl ≈ PropsSI("H", "T", T, "Q", 1, N2O) - PropsSI("H", "T", T, "Q", 0, N2O) rtol = 1e-5
        @test st.cpl ≈ PropsSI("C", "T", T, "Q", 0, N2O) rtol = 1e-4
        @test F.Tsat(st.P) ≈ T atol = 1e-4
    end
end

@testset "properties: single phase vs CoolProp, and state solves" begin
    for (P, T) in ((60e5, 270.0), (100e5, 300.0), (20e5, 300.0), (1e5, 290.0), (40e5, 320.0), (90e5, 330.0), (75e5, 312.0), (300e5, 300.0))
        ρ = PropsSI("D", "P", P, "T", T, N2O)
        e = F.eos(ρ, T)
        @test e.P ≈ P rtol = 1e-8
        @test e.cv ≈ PropsSI("O", "D", ρ, "T", T, N2O) rtol = 1e-8
        @test F.cp_of(e, ρ, T) ≈ PropsSI("C", "D", ρ, "T", T, N2O) rtol = 1e-8
        @test F.state_ρu(ρ, e.u).T ≈ T atol = 1e-6
        @test F.state_Ph(P, e.h).ρ ≈ ρ rtol = 1e-7
        @test F.state_Ps(P, e.s).ρ ≈ ρ rtol = 1e-7
    end
    for T in (190.0, 250.0, 290.0, 305.0, 309.0), x in (0.0, 0.03, 0.5, 0.97, 1.0)
        f = F._mix(F.sat_T(T), x)
        @test F.state_ρu(f.ρ, f.u).x ≈ x atol = 1e-5
        @test F.state_Ph(f.P, f.h).x ≈ x atol = 1e-8
    end
end

@testset "properties: the tables in docs/MODEL.md section 9" begin
    # saturation table: psig → Tsat °C, ρl, ρv, h_fg kJ/kg, cp_l kJ/kg/K
    for (pg, Tc, ρl, ρv, hfg, cp) in ((550, 8.8, 859, 111, 208, 2.52), (700, 18.9, 794, 152, 174, 3.09), (850, 27.4, 719, 208, 135, 4.54),
                                       (950, 32.4, 651, 264, 99, 8.36))
        st = F.sat_P(psia(pg))
        @test st.T - 273.15 ≈ Tc atol = 0.06
        @test st.ρl ≈ ρl atol = 0.6
        @test st.ρv ≈ ρv atol = 0.6
        @test (st.hv - st.hl) / 1e3 ≈ hfg atol = 0.6
        @test st.cpl / 1e3 ≈ cp atol = 0.006
    end
    # fraction of 30 °C bottle liquid that flashes entering a 550 / 700 / 850 psig tank
    hl = F.sat_T(303.15).hl
    @test F.state_Ph(psia(550), hl).x ≈ 0.292 atol = 0.002
    @test F.state_Ph(psia(700), hl).x ≈ 0.203 atol = 0.002
    @test F.state_Ph(psia(850), hl).x ≈ 0.071 atol = 0.002
    # a tank filled at 550 psig goes liquid-full warming to 30 °C above 77 % liquid by volume
    st = F.sat_P(psia(550))
    @test (F.sat_T(303.15).ρl - st.ρv) / (st.ρl - st.ρv) ≈ 0.77 atol = 0.01
    @test F.warm_soak(0.80 * st.ρl * 0.01916 + 0.20 * st.ρv * 0.01916, 0.01916, 303.15)[2] == 1.0
end

@testset "properties: ESDU 91022 transport tables" begin
    @test F.visc_l(273.15) * 1e3 ≈ 0.0774 atol = 1e-4
    @test F.visc_l(293.15) * 1e3 ≈ 0.0570 atol = 1e-4
    @test F.visc_v(273.15) * 1e6 ≈ 15.2 atol = 0.1
    @test F.cond_l(273.15) * 1e3 ≈ 103.0 atol = 0.1
    @test F.cond_v(273.15) * 1e3 ≈ 19.1 atol = 0.1
end

@testset "flow law" begin
    # Subcooled liquid discharging above its vapour pressure is plain incompressible flow.
    up = F.eos(PropsSI("D", "P", 80e5, "T", 280.0, N2O), 280.0)
    @test F.massflux(80e5, up.h, 70e5) ≈ sqrt(2 * PropsSI("D", "P", 80e5, "T", 280.0, N2O) * 10e5) rtol = 1e-3
    # Low-pressure vapour, choked: the ideal-gas critical flux.
    T0, P0 = 300.0, 4e5
    ρ0 = PropsSI("D", "P", P0, "T", T0, N2O)
    γ = PropsSI("C", "P", P0, "T", T0, N2O) / PropsSI("O", "P", P0, "T", T0, N2O)
    ideal = P0 * sqrt(γ / (P0 / ρ0)) * (2 / (γ + 1))^((γ + 1) / (2 * (γ - 1)))
    @test F.massflux(P0, F.eos(ρ0, T0).h, 1e5) ≈ ideal rtol = 0.03
    @test F.massflux(P0, F.eos(ρ0, T0).h, 1.5e5) ≈ F.massflux(P0, F.eos(ρ0, T0).h, 1e5) rtol = 1e-6   # choked
    # Saturated liquid: between HEM and SPI, rising with ΔP, zero at zero, no back flow.
    hl = F.sat_T(293.15).hl
    P1 = F.Psat(293.15)
    g = [F.massflux(P1, hl, P1 - d) for d in (0.0, 10.0, 1e3, 1e5, 1e6, 3e6)]
    @test g[1] == 0 && issorted(g)
    @test F.massflux(P1, hl, P1 + 1e5) == 0
    @test g[5] < sqrt(2 * F.sat_T(293.15).ρl * 1e6)
    @test F.friction_factor(1e5, 1e-4) ≈ 0.0185 atol = 0.001          # Moody chart
end

@testset "scene: edit, check, save and load" begin
    s = rig()
    @test F.check(s) == (String[], String[])
    path = tempname() * ".toml"
    F.save_scene(path, s)
    s2 = F.load_scene(path)
    @test [(p.kind, p.name, p.x, p.y, p.rot, p.p) for p in s2.parts] == [(p.kind, p.name, p.x, p.y, p.rot, p.p) for p in s.parts]
    @test Set((w.a, w.b) for w in s2.wires) == Set((w.a, w.b) for w in s.wires)
    @test s2.settings == s.settings
    F.rename_part!(s2, "FILL", "AV203")
    @test any(w -> w.b == ("AV203", "a"), s2.wires) && F.find_part(s2, "FILL") === nothing
    F.delete_part!(s2, "AV203")
    @test !any(w -> "AV203" in (w.a[1], w.b[1]), s2.wires)
    # wire corners survive a save, a rename and a group move; copies get fresh names and their own wires
    s3 = rig()
    w = s3.wires[findfirst(w -> w.a == ("FILL", "b"), s3.wires)]
    push!(w.via, (6, 2), (7, 5))
    F.save_scene(path, s3)
    @test occursin("via = [[6, 2], [7, 5]]", read(path, String))
    @test F.load_scene(path).wires[findfirst(w -> w.a == ("FILL", "b"), s3.wires)].via == [(6, 2), (7, 5)]
    F.rename_part!(s3, "FILL", "AV1")
    @test s3.wires[findfirst(w -> w.a == ("AV1", "b"), s3.wires)].via == [(6, 2), (7, 5)]
    F.move_parts!(s3, ["AV1", "TANK"], 2, -1)
    @test (s3.parts[F.find_part(s3, "TANK")].x, s3.parts[F.find_part(s3, "TANK")].y) == (10, -1)
    @test s3.wires[findfirst(w -> w.a == ("AV1", "b"), s3.wires)].via == [(8, 1), (9, 4)]        # both ends moved: corners follow
    @test isempty(s3.wires[findfirst(w -> w.b == ("AV1", "a"), s3.wires)].via)
    np, nw = length(s3.parts), length(s3.wires)
    new = F.paste_parts!(s3, F.copy_parts(s3, ["AV1", "TANK", "VENT"]), 0, 10)
    @test length(new) == 3 && length(s3.parts) == np + 3 && length(s3.wires) == nw + 2 && allunique(p.name for p in s3.parts)
    copy_wire = s3.wires[findfirst(w -> w.a == (new[1], "b"), s3.wires)]
    @test copy_wire.b == (new[2], "bottom") && copy_wire.via == [(8, 11), (9, 14)]
    @test s3.parts[F.find_part(s3, new[2])].p == s3.parts[F.find_part(s3, "TANK")].p && s3.parts[F.find_part(s3, new[2])].y == 9
    @test isempty(F.check(s3)[1])
    bad = F.Scene()
    F.add_part!(bad, "Bottle", 0, 0); F.add_part!(bad, "Tank", 4, 0)
    F.connect!(bad, ("BOTTLE1", "out"), ("TANK1", "bottom"))
    @test occursin("wired straight together", F.check(bad)[1][1])
    @test_throws ErrorException F.compile(bad)
    # every example scene loads and compiles
    for f in filter(endswith(".toml"), readdir(joinpath(@__DIR__, "..", "scenes"); join = true))
        sc = F.load_scene(f)
        @test isempty(F.check(sc)[1])
        F.compile(sc)
    end
end

@testset "test logs: reading, column meaning, overlay and schedule" begin
    scene = F.load_scene(joinpath(@__DIR__, "..", "scenes", "hf4.toml"))
    path = joinpath(@__DIR__, "data", "hf4_fill_10hz.csv")
    names, cols = F.read_csv(path)
    _, thin = F.read_csv(path; min_dt = 1.0)
    @test length(cols["t_s"]) > 9000 && 950 < length(thin["t_s"]) < 1100 && minimum(diff(thin["t_s"])) >= 1.0
    columns = F.guess_log_columns(names, scene)
    by(name) = columns[findfirst(c -> c.name == name, columns)]
    @test by("t_s").role === :time && by("PT202_psig").role === :pressure && by("LCFill_kg").role === :mass
    @test (by("AV203").role, by("AV203").part, by("AV203").one_means) == (:valve, "AV203", :open)
    @test (by("heater").role, by("heater").part) == (:heater, "TANK")
    @test by("AV205").role === :other                           # no such part in the scene
    sched = F.log_schedule(cols, columns; t0 = 9020.0)
    @test sched == F.schedule_from_log(cols["t_s"], cols, Dict("AV201" => ("AV201", "open", "close"), "AV202" => ("AV202", "open", "close"),
        "AV203" => ("AV203", "open", "close"), "AV204" => ("AV204", "open", "close"), "heater" => ("TANK", "on", "off")); t0 = 9020.0)
    by("AV204").one_means = :closed                             # the same column read the other way round
    flipped = F.log_schedule(cols, columns; t0 = 9020.0)
    @test [a.action for a in flipped if a.part == "AV204"] == [a.action == "open" ? "close" : "open" for a in sched if a.part == "AV204"]
    overlay = F.log_overlay(cols, columns, 101325.0)
    @test Set(first.(overlay)) == Set(["PT201_psig", "PT202_psig", "LCFill_kg"])
    p202 = overlay[findfirst(o -> o[1] == "PT202_psig", overlay)]
    @test p202[2] === :pressure && p202[3][1] ≈ cols["PT202_psig"][1] * F.PSI + 101325
    # other naming generations, and what must not be guessed
    raw = F.guess_log_columns(["device_timestamp [s]", "LCFILL [kg]", "LCTHRUST [kg]", "PT202 [PSI]", "valve_AV-203", "kasa_Run_Heater", "kasa_Supply_Heater"], scene)
    @test [(c.role, c.part, c.use) for c in raw[2:end]] == [(:mass, "", true), (:mass, "", false), (:pressure, "", true), (:valve, "AV203", true), (:heater, "TANK", true), (:heater, "", false)]
    vec = F.guess_log_columns(["device_timestamp", "source", "PT202 [PSI]", "relay_Safe24", "solenoid_AV203"], scene)
    @test vec[2].role === :other && vec[4].role === :other
    @test (vec[5].part, vec[5].one_means, vec[5].use) == ("AV203", :unknown, false)   # 1 = CLOSED there: never assumed
    # blank cells and text columns read as NaN and are not state changes
    tmp = tempname() * ".csv"
    write(tmp, "t,source,solenoid_AV203,PT202 [PSI]\n0.0,pad,0,1.5\n1.0,pad,,1.6\n2.0,pad,0,\n3.0,pad,1,1.8\n")
    n2, c2 = F.read_csv(tmp)
    @test isnan(c2["source"][1]) && isnan(c2["solenoid_AV203"][2]) && isnan(c2["PT202 [PSI]"][3])
    cols2 = F.guess_log_columns(n2, scene)
    cols2[3].one_means, cols2[3].use = :closed, true
    @test [(a.t, a.action) for a in F.log_schedule(c2, cols2; t0 = 0.0)] == [(0.0, "open"), (3.0, "close")]
end

@testset "log manifests: several boards on one clock, polarity, bounce, groups, gaps, phases" begin
    # A made-up log from three boards. True time runs 0-120 s. Board A (the reference) stamps its
    # rows 1000 s ahead of true time; B's clock is 100.18 s behind A's and C's 1869.08 s behind.
    # Rows are written in the order they arrive.
    dir = mktempdir()
    rows = Tuple{Float64,String}[]
    for k in 0:1200                                              # A, 10 Hz: a pressure
        tt = 0.1k
        push!(rows, (tt, "$(round(1000 + tt, digits = 3)),A,$(round(100 + 2tt, digits = 3)),,,"))
    end
    for k in 0:599                                               # B, 5 Hz: a valve column where 1 means closed
        tt = 0.2k + 0.03
        v = 20 <= tt < 60 ? 0 : 1                                # open from 20 s to 60 s ...
        abs(tt - 60.23) < 0.01 && (v = 0)                        # ... with one bounce as it closes
        push!(rows, (tt, "$(round(1000 + tt - 100.18, digits = 3)),B,,$v,,"))
    end
    for k in 0:239                                               # C, 2 Hz: two thermistors, silent from 80 s to 90 s
        tt = 0.5k + 0.05
        80 <= tt < 90 && continue
        push!(rows, (tt, "$(round(1000 + tt - 1869.08, digits = 3)),C,,,$(round(20 - 0.1tt, digits = 3)),$(round(20.4 - 0.1tt, digits = 3))"))
    end
    sort!(rows; by = first)
    write(joinpath(dir, "log.csv"), "device_timestamp,source,PT [PSI],solenoid_V1,TH1 [C],TH2 [C]\n" * join(last.(rows), "\n") * "\n")
    manifest = joinpath(dir, "log.toml")
    body = """
        file = "log.csv"
        time = "device_timestamp"
        source_column = "source"
        reference_source = "A"
        clock_offsets = "auto"
        sample_dt = 1.0
        debounce_s = 0.5
        t0 = 1010.0
        [[channel]]
        column = "PT [PSI]"
        role = "pressure"
        maps_to = "TANK.P"
        name = "PT_psig"
        [[channel]]
        column = "solenoid_V1"
        role = "valve"
        maps_to = "FILL"
        one_means = "closed"
        [[channel]]
        column = "TH1 [C]"
        role = "temperature"
        group = "TH"
        maps_to = "TANK.T"
        [[channel]]
        column = "TH2 [C]"
        role = "temperature"
        group = "TH"
        [phases]
        fill = { from = "FILL open", to = "FILL close" }
        hold = [60, 100]
        """
    write(manifest, body)
    log = F.load_log(manifest)
    @test abs(log.offsets["B"] - 100.18) < 0.1 && abs(log.offsets["C"] - 1869.08) < 0.1 && log.offsets["A"] == 0
    @test isempty(log.warnings)
    @test log.names == ["t_s", "PT_psig", "solenoid_V1", "TH"] && log.cols["t_s"][1:3] == [1000.0, 1001.0, 1002.0]
    at(name, t) = log.cols[name][findfirst(==(t), log.cols["t_s"])]
    @test abs(at("PT_psig", 1050.0) - 200) < 0.3                  # median of the readings around true time 50 s
    @test abs(at("TH", 1050.0) - (20.2 - 5)) < 0.15               # the two thermistors averaged, on the reference clock
    @test isnan(at("TH", 1085.0)) && !isnan(at("TH", 1079.0)) && !isnan(at("TH", 1091.0))   # a gap stays a gap
    sched = F.log_schedule(log.cols, log.columns; t0 = log.t0)
    @test [(a.part, a.action) for a in sched] == [("FILL", "close"), ("FILL", "open"), ("FILL", "close")]   # no extra pair from the bounce
    @test abs(sched[2].t - 10) <= 1 && abs(sched[3].t - 50) <= 1.5
    @test log.phases["hold"] == (60.0, 100.0) && log.phases["fill"] == (sched[2].t, sched[3].t)
    @test log.columns[2].target == "TANK.P" && log.columns[3].one_means === :closed && log.columns[4].role === :temperature
    over = F.log_overlay(log.cols, log.columns, 101325.0; saturation = true)
    @test Set(first.(over)) == Set(["PT_psig", "Tsat of PT_psig", "TH", "Psat of TH"])
    k = findfirst(==(1050.0), log.cols["t_s"])
    @test over[findfirst(o -> o[1] == "Tsat of PT_psig", over)][3][k] ≈ F.Tsat(200 * F.PSI + 101325) rtol = 1e-3
    @test isnan(over[findfirst(o -> o[1] == "Psat of TH", over)][3][findfirst(==(1085.0), log.cols["t_s"])])
    # written out, it is a plain CSV that replays by name with no manifest: valves as 1 = open
    plain = joinpath(dir, "plain.csv")
    F.write_log(plain, log)
    pn, pc = F.read_csv(plain)
    sc = F.Scene(); F.add_part!(sc, "Valve", 0, 0; name = "solenoid_V1")
    @test pn == log.names && isnan(pc["TH"][findfirst(==(1085.0), pc["t_s"])]) && pc["solenoid_V1"][findfirst(==(1050.0), pc["t_s"])] == 1
    # given offsets are used as given; mistakes in a manifest are errors, not silently ignored
    write(manifest, replace(body, "clock_offsets = \"auto\"" => "clock_offsets = { B = 100.0, C = 1869.0 }"))
    @test F.load_log(manifest).offsets == Dict("A" => 0.0, "B" => 100.0, "C" => 1869.0)
    write(manifest, replace(body, "debounce_s" => "debounce"))
    @test_throws ErrorException F.load_log(manifest)
    write(manifest, replace(body, "one_means = \"closed\"" => ""))
    @test_throws ErrorException F.load_log(manifest)
    write(manifest, replace(body, "column = \"TH2 [C]\"" => "column = \"TH9 [C]\""))
    @test_throws ErrorException F.load_log(manifest)
    # a drifting clock is reported
    drift = [(tt, startswith(split(r, ',')[2], "B") ? replace(r, r"^[^,]+" => string(round(1000 + tt - 100.18 - 0.02tt, digits = 3))) : r) for (tt, r) in rows]
    write(joinpath(dir, "log.csv"), "device_timestamp,source,PT [PSI],solenoid_V1,TH1 [C],TH2 [C]\n" * join(last.(drift), "\n") * "\n")
    write(manifest, body)
    @test any(occursin("drifts", w) for w in F.load_log(manifest).warnings)
    # the draft generator does not guess what a 1 means on a solenoid_ column
    s = F.Scene(); F.add_part!(s, "Valve", 0, 0; name = "V1")
    draft = F.suggest_manifest(joinpath(dir, "log.csv"), s)
    @test occursin("one_means = \"?\"", draft) && occursin("source_column = \"source\"", draft) && occursin("role = \"temperature\"", draft)
    @test F.debounce([0.0, 1, 1.2, 1.4, 3, 4], [0.0, 1, 0, 1, 1, 0], 0.5) == [0.0, 0, 0, 1, 1, 0]
end

@testset "studies: expansion, scoring, metrics, resume, failures, threads" begin
    dir = mktempdir()
    F.save_scene(joinpath(dir, "rig.toml"), rig())
    study(body) = (p = joinpath(dir, "s$(hash(body)).toml"); write(p, body); p)
    table(path) = (lines = readlines(path); (split(lines[1], ','), [split(l, ',') for l in lines[2:end]]))
    # `top` is the part of a study file that must come before any table; `rest` is the tables
    top = "scene = \"rig.toml\"\nuntil = 20\n"
    rest = """
        [[case]]
        name = "fill"
        actions = [{ t = 1.0, part = "FILL", action = "open" }]
        [metrics]
        loaded_kg = { fn = "final", signal = "TANK.m" }
        peak_psig = { fn = "max", signal = "TANK.P", unit = "psig" }
        t_100_psig = { fn = "first_crossing", signal = "TANK.P", value_psig = 100 }
        """
    vary_grid = "[vary]\n\"BOTTLE.T_C\" = [15, 25]\n\"FILL.d_mm\" = { from = 2, to = 6, n = 3 }\n"
    grid = study(top * rest * vary_grid)
    st = F.read_study(grid)
    @test length(st.runs) == 6 && st.varied == ["BOTTLE.T_C", "FILL.d_mm"] && st.runs[1].id == "BOTTLE.T_C=15.0;FILL.d_mm=2.0|fill"
    out = F.study(grid; quiet = true)
    head, rows = table(out)
    col(name) = [r[findfirst(==(name), head)] for r in rows]
    @test length(rows) == 6 && all(isempty, col("failed")) && Set(col("id")) == Set(r.id for r in st.runs)
    num(name, id) = parse(Float64, rows[findfirst(r -> r[1] == id, rows)][findfirst(==(name), head)])
    @test num("loaded_kg", "BOTTLE.T_C=25.0;FILL.d_mm=6.0|fill") > num("loaded_kg", "BOTTLE.T_C=25.0;FILL.d_mm=2.0|fill")      # a bigger valve fills faster
    @test num("peak_psig", "BOTTLE.T_C=25.0;FILL.d_mm=6.0|fill") > num("peak_psig", "BOTTLE.T_C=15.0;FILL.d_mm=6.0|fill") > 100  # a warmer bottle pushes harder
    @test 1 < num("t_100_psig", "BOTTLE.T_C=25.0;FILL.d_mm=6.0|fill") < 20
    # metrics read the history as lines between samples, so sparse samples do not leave holes
    s = rig(); sim = F.Sim(s)
    empty!(sim.t_hist); append!(sim.t_hist, [0.0, 10.0, 20.0]); sim.hist["X"] = [0.0, 100.0, 100.0]
    @test F.metric(sim, Dict("fn" => "mean_over", "signal" => "X", "from" => 4, "to" => 6), 0.0) ≈ 50        # no sample lies in 4-6 s
    @test F.metric(sim, Dict("fn" => "time_in_band", "signal" => "X", "lo" => 50, "hi" => 200), 0.0) ≈ 15    # the second half of the ramp, then the plateau
    @test F.metric(sim, Dict("fn" => "first_crossing", "signal" => "X", "value" => 25), 0.0) ≈ 2.5
    # resume: only what is missing is run again
    lines = readlines(out)
    write(out, join(lines[1:3], "\n") * "\n")
    F.study(grid; resume = true, quiet = true)
    again = readlines(out)
    @test length(again) == 7 && again[2:3] == lines[2:3] && Set(first.(split.(again[2:end], ','))) == Set(r.id for r in st.runs)
    # the other ways of choosing runs
    oat = F.read_study(study(top * "mode = \"oat\"\n" * rest * "[vary]\n\"BOTTLE.T_C\" = [15, 25]\n\"FILL.d_mm\" = [2, 4, 6]\n"))
    @test length(oat.runs) == 6 && isempty(oat.runs[1].overrides) && all(length(r.overrides) == 1 for r in oat.runs[2:end])
    lhs_body = top * "mode = \"lhs\"\nn = 8\nseed = 3\n" * rest * "[vary]\n\"BOTTLE.T_C\" = { from = 10, to = 30 }\n\"FILL.d_mm\" = { from = 1, to = 8, log = true }\n"
    lhs = F.read_study(study(lhs_body))
    temps = sort([Dict(r.overrides)["BOTTLE.T_C"] for r in lhs.runs])
    @test length(lhs.runs) == 8 && all(10 + 2.5(k - 1) <= temps[k] <= 10 + 2.5k for k in 1:8)       # one draw in each of the eight strata
    @test [r.id for r in F.read_study(study(lhs_body)).runs] == [r.id for r in lhs.runs]             # the same seed draws the same runs
    # a run that cannot be made, or fails, is a row that says why; a mistake in the file stops before anything runs
    out = F.study(study(top * rest * "[vary]\n\"TANK.volume_L\" = [19.16, -1]\n"); quiet = true)
    head, rows = table(out)
    failed = [r[findfirst(==("failed"), head)] for r in rows]
    @test length(rows) == 2 && count(isempty, failed) == 1 && any(occursin("must be positive", f) for f in failed)
    @test_throws ErrorException F.read_study(study(top * rest * "[vary]\n\"TANK.volume\" = [1, 2]\n"))      # no such parameter
    @test_throws ErrorException F.read_study(study(top * "sample = 1\n" * rest))                                # no such key
    @test_throws ErrorException F.read_study(study(top * replace(rest, "\"final\"" => "\"last\"")))             # no such metric
    @test_throws ErrorException F.read_study(study(top * rest * "oops = 1\n"))                                  # a key that landed in [metrics]
    # scored against a log: a one-run study gives what a single replay gives
    root = joinpath(@__DIR__, "..")
    one = study("""
        scene = "$(joinpath(root, "scenes", "hf3cf.toml"))"
        replay = "$(joinpath(root, "test", "data", "hf3cf_fill_10hz.csv"))"
        t0 = 15290
        until = 265
        [score]
        channels = [{ log = "PT202_psig", model = "TANK.P" }, { log = "PT201_psig", model = "PT201.P", weight = 0.5 }]
        phases = { sealed = [106.4, 200.8], fill = { from = "AV201 open", to = "AV201 close#3" } }
        """)
    st1 = F.read_study(one)
    @test all(isapprox.(Dict(st1.phases)["fill"], (24.4, 103.4); atol = 1e-6))     # first opening to third closing of AV201
    head, rows = table(F.study(one; quiet = true))
    names, D = F.read_csv(joinpath(root, "test", "data", "hf3cf_fill_10hz.csv"))
    s = F.load_scene(joinpath(root, "scenes", "hf3cf.toml"))
    s.schedule = F.log_schedule(D, F.guess_log_columns(names, s); t0 = 15290.0)
    sim = F.Sim(s; sample_dt = 0.5)
    F.run_until!(sim, 265.0)
    direct = sqrt(sum(((F.signal_at(sim, "TANK.P", t - 15290) - 101325) / F.PSI - D["PT202_psig"][i])^2 for (i, t) in enumerate(D["t_s"]) if 106.4 <= t - 15290 <= 200.8) /
                  count(t -> 106.4 <= t - 15290 <= 200.8, D["t_s"]))
    @test parse(Float64, rows[1][findfirst(==("rms.PT202_psig.sealed"), head)]) ≈ direct rtol = 1e-6
    @test 2 < direct < 3
    # on two threads the results are the same to the last digit
    script = "using N2OFillSim; N2OFillSim.study(ARGS[1]; quiet = true)"
    two = study(top * rest * vary_grid * "# threaded\n")
    run(`$(Base.julia_cmd()) -t 2 --project=$(Base.active_project()) -e $script $two`)
    strip_wall(path) = sort([join(split(l, ',')[1:end-1], ',') for l in readlines(path)[2:end]])
    @test strip_wall(splitext(two)[1] * "_results/results.csv") == strip_wall(splitext(grid)[1] * "_results/results.csv")
end

@testset "command line: overrides, schedule file, columns, defaults, warnings" begin
    dir = mktempdir()
    scene, sched, out = joinpath(dir, "s.toml"), joinpath(dir, "sched.csv"), joinpath(dir, "out.csv")
    F.save_scene(scene, rig())
    write(sched, "t_s,part,action,value\n1.0,FILL,open,0.0\n")
    quiet(args) = redirect_stdout(() -> F.main(args), devnull)
    r = quiet([scene, "--headless", "--set", "BOTTLE.T_C=20", "--set", "settings.T_amb_C=20", "--schedule", sched, "--until", "10",
               "--sample-dt", "0.5", "--csv", out, "--columns", "TANK.P, TANK.m"])
    @test isempty(r.failed) && abs(r.t_s - 10) < 1e-6 && abs(r.mass_error) < 1e-6
    names, cols = F.read_csv(out)
    @test names == ["t_s", "TANK.P", "TANK.m"] && 3 <= length(cols["t_s"]) <= 21 && minimum(diff(cols["t_s"])) >= 0.5 - 1e-6
    @test cols["TANK.m"][end] > 0.5                              # the schedule file opened the fill valve
    @test length(F.read_events(F.events_path(out))) == 1
    @test_throws ErrorException quiet([scene, "--headless", "--set", "TANK.volume=1"])          # no such parameter
    @test_throws ErrorException quiet([scene, "--headless", "--set", "settings.T_amb=1"])       # no such setting
    @test_throws ErrorException quiet([scene, "--headless", "--until", "1", "--csv", out, "--columns", "TANK.nope"])
    @test_throws ErrorException quiet([scene, "--headless", "--bogus"])
    @test_throws ErrorException quiet([joinpath(dir, "missing.toml"), "--headless"])
    # a replay runs to the end of its log unless told otherwise
    root = joinpath(@__DIR__, "..")
    r = quiet([joinpath(root, "scenes", "hf3cf.toml"), "--headless", "--replay", joinpath(root, "test", "data", "hf3cf_fill_10hz.csv"), "--t0", "15290"])
    @test abs(r.t_s - 269.9) < 0.2 && isempty(r.failed)
    # a heater nobody can switch on is worth a warning when nobody is watching, and only then
    s = rig()
    s.parts[F.find_part(s, "TANK")].p["heater_W"] = 100.0
    @test !any(occursin("heater", w) for w in F.check(s)[2])
    @test any(occursin("never switches on", w) for w in F.check(s; unattended = true)[2])
    push!(s.schedule, F.Action(5.0, "TANK", "on", 0.0))
    @test !any(occursin("heater", w) for w in F.check(s; unattended = true)[2])
    push!(s.schedule, F.Action(5.0, "BOTTLE", "on", 0.0))
    @test any(occursin("heater_W is 0", w) for w in F.check(s)[2])
end

@testset "unvented fill stalls below bottle pressure, conserving mass and energy" begin
    sim = F.Sim(rig(bottle_T = 28.0))
    F.act!(sim, "FILL", "open")
    F.run_until!(sim, 200.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    t, b = tank(sim), bottle(sim)
    @test abs(b.P - t.P) < 0.5e5                # equalised
    @test 0.15 < t.liquid_frac < 0.85           # ... with the tank only part full
    @test b.T < 28.0 + 273.15 - 1               # the bottle cooled as it delivered
    @test abs(sim.hist["FILL.mdot"][end]) < 0.01
    @test closes(sim)
end

@testset "wall coupling is set per vessel" begin
    # Cutting the bottle's wall off from its contents leaves the tank's wall coupling alone.
    function run(bottle_mult)
        s = rig()
        s.parts[F.find_part(s, "BOTTLE")].p["h_in_mult"] = bottle_mult
        sim = F.Sim(s)
        F.act!(sim, "FILL", "open")
        F.run_until!(sim, 40.0)
        sim
    end
    a, b = run(1.0), run(0.0)
    drop(sim, key) = sim.hist[key][1] - sim.hist[key][end]
    @test drop(a, "BOTTLE.T") > 1                                 # the bottle cools as it delivers
    @test drop(a, "BOTTLE.Tw") > 10 * abs(drop(b, "BOTTLE.Tw"))   # and takes its wall with it only when coupled
    @test drop(b, "BOTTLE.T") > drop(a, "BOTTLE.T")               # uncoupled, the contents get no heat back from the wall
    @test_throws ErrorException F.apply_overrides!(rig(), ["BOTTLE.h_in_mul=0"])
end

@testset "venting lets more in and costs vented mass" begin
    closed = F.Sim(rig(bottle_T = 28.0)); F.act!(closed, "FILL", "open"); F.run_until!(closed, 200.0)
    vented = F.Sim(rig(bottle_T = 28.0)); F.act!(vented, "FILL", "open"); F.act!(vented, "VENT", "open"); F.run_until!(vented, 200.0)
    @test isempty(vented.failed)
    @test tank(vented).m > tank(closed).m + 0.5
    @test tank(vented).P < tank(closed).P - 1e5
    @test vented.u[end-2] > 0.3           # kg gone out of the vent
    @test closes(vented)
end

@testset "the vent draws vapour below the dip tube and liquid above it" begin
    flow(fill) = (sim = F.Sim(rig(tank_fill = fill, tank_T = 15.0, air = 0.0)); F.act!(sim, "VENT", "open"); F.run_until!(sim, 2.0);
                  (tank(sim).liquid_frac, sim.hist["VENT_OR.mdot"][end]))
    (lf_low, m_low), (lf_high, m_high) = flow(8.0), flow(15.3)      # dip tube at 90 % of the height
    @test lf_low < 0.6 && lf_high > 0.92
    @test m_high > 3 * m_low > 0
end

@testset "vent stuck closed is the same as never venting; fail-safe drops valves to their fail state" begin
    sim = F.Sim(rig(bottle_T = 28.0))
    F.act!(sim, "VENT", "stuck_closed"); F.act!(sim, "VENT", "open"); F.act!(sim, "FILL", "open")
    F.run_until!(sim, 60.0)
    @test sim.u[end-2] < 1e-6
    F.fail_safe!(sim)
    F.run_until!(sim, 62.0)
    @test all(e -> e.opening == 0, filter(e -> e.kind === :valve, sim.net.edges))
end

@testset "vent stuck open: the tank cannot be closed up and keeps losing nitrous" begin
    sim = F.Sim(rig(bottle_T = 28.0))
    F.act!(sim, "VENT", "open"); F.act!(sim, "VENT", "stuck_open"); F.act!(sim, "FILL", "open")
    F.run_until!(sim, 60.0)
    F.act!(sim, "VENT", "close"); F.act!(sim, "FILL", "close")
    F.run_until!(sim, 70.0)
    lost = sim.u[end-2]
    F.run_until!(sim, 300.0)
    @test isempty(sim.failed)
    @test sim.net.edges[findfirst(e -> e.name == "VENT", sim.net.edges)].opening == 1
    @test sim.u[end-2] > lost + 0.5
    @test sim.hist["TANK.P"][end] < sim.hist["TANK.P"][findfirst(>=(70.0), sim.t_hist)] - 3e5
    @test closes(sim)
end

@testset "bottle colder than the tank: flow runs backwards into the bottle" begin
    s = rig(bottle_T = 5.0, tank_fill = 9.0, tank_T = 25.0, air = 0.0)
    s.parts[F.find_part(s, "BOTTLE")].p["fill_kg"] = 10.0
    sim = F.Sim(s)
    m0 = bottle(sim).m
    F.act!(sim, "FILL", "open")
    F.run_until!(sim, 120.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    @test minimum(sim.hist["FILL.mdot"]) < -0.01            # a → b is bottle → tank, so negative
    @test bottle(sim).m > m0 + 1
    @test abs(bottle(sim).P - tank(sim).P) < 1e5
    @test closes(sim)
end

@testset "liquid trapped in a small volume: warming lifts its relief valve" begin
    # The lines themselves hold no fluid, so a trapped section is modelled as a small Tank.
    s = F.Scene()
    s.settings["T_amb_C"], s.settings["solar_W_m2"] = 40.0, 1000.0
    t = F.add_part!(s, "Tank", 0, 0; name = "LINE")
    t.p["volume_L"], t.p["diameter_mm"], t.p["wall_kg"], t.p["absorptivity"], t.p["air"] = 0.3, 12.0, 0.8, 0.9, 0.0
    t.p["T_C"], t.p["fill_kg"] = 15.0, 0.95 * 0.3e-3 * F.sat_T(288.15).ρl
    F.add_part!(s, "ReliefValve", 4, 0; name = "PRV"); F.add_part!(s, "Atmosphere", 7, 0; name = "ATM")
    F.connect!(s, ("LINE", "top"), ("PRV", "a")); F.connect!(s, ("PRV", "b"), ("ATM", "a"))
    sim = F.Sim(s; sample_dt = 1.0)
    F.run_until!(sim, 900.0; wall_limit = 120)
    @test isempty(sim.failed)
    @test maximum(sim.hist["LINE.liquid_frac"]) == 1.0
    @test sim.u[end-2] > 0.005                       # the relief valve passed nitrous
    @test maximum(sim.hist["LINE.P"]) < 1150 * F.PSI + 101325
    @test closes(sim)
end

@testset "valve, heater and ambient state are recorded; the action log repeats the run" begin
    s = rig(bottle_T = 28.0)
    s.parts[F.find_part(s, "TANK")].p["heater_W"] = 200.0
    a = F.Sim(s)
    F.act!(a, "FILL", "open");  F.run_until!(a, 20.0)
    F.act!(a, "VENT", "open");  F.run_until!(a, 26.0)
    F.act!(a, "VENT", "stuck_open"); F.act!(a, "VENT", "close"); F.run_until!(a, 30.0)
    F.act!(a, "VENT", "unstick"); F.act!(a, "TANK", "on"); F.run_until!(a, 60.0)
    h, t = a.hist, a.t_hist
    at(key, time) = h[key][searchsortedfirst(t, time)]
    @test at("FILL.position", 0.0) == 0 && at("FILL.position", 10.0) == 1 && at("FILL.open_cmd", 10.0) == 1
    @test 0 < at("FILL.position", 0.15) < 1                       # part-way through its stroke
    @test at("VENT.position", 28.0) == 1 && at("VENT.open_cmd", 28.0) == 0 && at("VENT.stuck", 28.0) == 1
    @test at("VENT.position", 40.0) == 0 && at("VENT.stuck", 40.0) == 0
    @test at("TANK.heater_on", 25.0) == 0 && at("TANK.heater_on", 40.0) == 1
    @test h["ambient.T"][end] == 298.15 && haskey(h, "PRV.set")
    @test all(v -> length(v) == length(t), values(h))            # every signal has every sample
    csv = tempname() * ".csv"
    F.write_csv(csv, a); F.write_events(F.events_path(csv), a)
    @test occursin("VENT.position", readline(csv))
    events = F.read_events(F.events_path(csv))
    @test [(e.part, e.action) for e in events] == [(e.part, e.action) for e in a.log] && length(events) == 6
    s.schedule = events
    b = F.Sim(s); F.run_until!(b, 60.0)
    @test tank(b).P ≈ tank(a).P rtol = 2e-4
    @test tank(b).m ≈ tank(a).m rtol = 2e-4
    @test b.u[end-2] ≈ a.u[end-2] rtol = 2e-3                    # same mass vented
end

@testset "overfilled tank warming in the sun goes liquid-full and lifts the relief valve" begin
    s = rig(tank_fill = 15.2, tank_T = 12.0, air = 0.0)       # ≈ 93 % liquid at 12 °C
    s.settings["T_amb_C"], s.settings["solar_W_m2"] = 38.0, 1000.0
    s.parts[F.find_part(s, "TANK")].p["absorptivity"] = 0.9
    sim = F.Sim(s; sample_dt = 5.0)
    @test tank(sim).liquid_frac > 0.9
    F.run_until!(sim, 4 * 3600.0; wall_limit = 240)
    @test isempty(sim.failed)
    @test maximum(sim.hist["TANK.liquid_frac"]) == 1.0          # went liquid-full
    @test maximum(sim.hist["PRV.mdot"]) > 0                     # relief valve lifted
    @test maximum(sim.hist["TANK.P"]) < (1100 * F.PSI + 101325) # and held it near its set pressure
    @test sim.u[end-2] > 0.1
    @test closes(sim)
end

@testset "flux has no step where a supercritical supply enters the two-phase region" begin
    e = F.eos(625.0, 61.7 + 273.15)                             # a bottle left in the sun: 146 bar, 62 °C
    G = [F.massflux(P, e.h, P - 30e5) for P in range(80e5, 62e5, length = 181)]   # 0.1 bar apart along its isenthalp
    @test all(isfinite, G) && maximum(abs.(diff(G)) ./ G[2:end]) < 0.02          # it was a 17 % step at 72 bar
end

@testset "a supercritical bottle opened into an empty tank does not stall the junction solve" begin
    # Reported from the window: a bottle left in the sun (62 °C, about 2100 psig) opened through an
    # orifice, hose, valve and hose. Every evaluation "failed" at an imbalance of 1e-6 and paid for
    # the whole fallback: 94 failures and 14 s of wall time for these 5 s. Now 2 failures.
    s = F.Scene()
    F.add_part!(s, "Bottle", 0, 0; name = "B").p["T_C"] = 62.0
    for (i, (kind, name)) in enumerate((("Orifice", "OR"), ("Pipe", "L1"), ("Valve", "V"), ("Pipe", "L2"), ("Tank", "T")))
        F.add_part!(s, kind, 3i, 0; name)
    end
    for (a, b) in ((("B", "out"), ("OR", "a")), (("OR", "b"), ("L1", "a")), (("L1", "b"), ("V", "a")), (("V", "b"), ("L2", "a")), (("L2", "b"), ("T", "bottom")))
        F.connect!(s, a, b)
    end
    sim = F.Sim(s)
    F.act!(sim, "V", "open")
    F.run_until!(sim, 20.0; wall_limit = 60)
    P(name) = F.node_report(sim.net, sim.net.nodes[findfirst(n -> n.name == name, sim.net.nodes)], sim.u).P
    @test isempty(sim.failed) && sim.net.unconverged <= 10
    @test abs(P("B") - P("T")) < 2e5 && P("T") > 60e5              # equalised, both far above the critical pressure
    @test closes(sim)
end

@testset "regulator holds the downstream pressure and flashes the liquid" begin
    s = rig(bottle_T = 28.0)
    F.delete_part!(s, "FILL")
    F.add_part!(s, "Regulator", 3, 0; name = "REG").p["set_psig"] = 700.0
    F.add_part!(s, "Valve", 5, 0; name = "FILL")
    F.add_part!(s, "Gauge", 4, -2; name = "PT")
    F.connect!(s, ("BOTTLE", "out"), ("REG", "a")); F.connect!(s, ("REG", "b"), ("FILL", "a"))
    F.connect!(s, ("FILL", "b"), ("TANK", "bottom")); F.connect!(s, ("PT", "a"), ("REG", "b"))
    sim = F.Sim(s)
    F.act!(sim, "FILL", "open"); F.act!(sim, "VENT", "open")
    F.run_until!(sim, 15.0)
    @test sim.hist["PT.T"][end] < bottle(sim).T - 3              # cold two-phase mixture downstream
    F.run_until!(sim, 120.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    @test maximum(sim.hist["PT.P"]) < 700 * F.PSI + 101325 + 1e3 # never above the set-point
    @test tank(sim).P < 700 * F.PSI + 101325
    @test closes(sim)
    # Behind a shut valve the regulator locks up at its set-point: not at bottle pressure, and not
    # at whatever the line held before. Before the fill valve first opens, and after it shuts mid-fill.
    locked(P) = abs(P - (700 * F.PSI + 101325)) < 1e3
    sim = F.Sim(s)
    F.run_until!(sim, 3.0)
    @test locked(sim.hist["PT.P"][end])
    F.act!(sim, "FILL", "open"); F.act!(sim, "VENT", "open")
    F.run_until!(sim, 10.0)
    @test tank(sim).P < 700 * F.PSI + 101325 - 1e5               # still filling, below the set-point
    F.act!(sim, "FILL", "close")
    F.run_until!(sim, 30.0)
    @test locked(sim.hist["PT.P"][end])
    @test isempty(sim.failed) && sim.net.unconverged == 0
end

# RMS difference, in psi, between a recorded pressure signal and a logged gauge column over [a, b] s of run time.
function log_rms(sim, key, D, col, t0, a, b, p_atm)
    d = [(sim.hist[key][clamp(searchsortedfirst(sim.t_hist, t - t0), 1, length(sim.t_hist))] - p_atm) / F.PSI - D[col][i]
         for (i, t) in enumerate(D["t_s"]) if a <= t - t0 <= b && isfinite(D[col][i])]
    sqrt(sum(abs2, d) / length(d))
end

@testset "HF3 cold flow replay (unvented fill, sealed period, blowdown through an injector)" begin
    # HF4's fitted values with only what was measured changed (see the scene). Numbers on
    # 2026-11-22: sealed-period tank pressure 2.5 psi RMS, fill line 3.2 psi, 4.59 kg loaded
    # against 4.3 kg on the load cell.
    names, D = F.read_csv(joinpath(@__DIR__, "data", "hf3cf_fill_10hz.csv"))
    s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "hf3cf.toml"))
    t0, p_atm = 15290.0, s.settings["P_atm_kPa"] * 1e3
    s.schedule = F.log_schedule(D, F.guess_log_columns(names, s); t0)
    sim = F.Sim(s)
    F.run_until!(sim, 265.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    @test log_rms(sim, "TANK.P", D, "PT202_psig", t0, 106.4, 200.8, p_atm) < 4       # sealed after the last fill pulse
    @test log_rms(sim, "PT201.P", D, "PT201_psig", t0, 106.4, 200.8, p_atm) < 5      # the bottle recovering
    @test 4.4 < sim.hist["TANK.m"][searchsortedfirst(sim.t_hist, 200.0)] < 4.8
    @test closes(sim)
end

@testset "LC26 day 4 replay (vented fill and a 1000 s sealed hold, with thermistors on the tank)" begin
    # Numbers on 2026-11-22: tank pressure 33 psi RMS through the vented fill and 4.9 psi through
    # the hold; the model's liquid 0.7 K below the thermistors on average.
    names, D = F.read_csv(joinpath(@__DIR__, "data", "lc26_day4_fill_1hz.csv"))
    s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "lc26_day4.toml"))
    t0, p_atm = 21950.0, s.settings["P_atm_kPa"] * 1e3
    s.schedule = F.log_schedule(D, F.guess_log_columns(names, s); t0)
    sim = F.Sim(s)
    F.run_until!(sim, 1215.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    @test log_rms(sim, "TANK.P", D, "PT202_psig", t0, 30, 150, p_atm) < 42            # vented fill
    @test log_rms(sim, "TANK.P", D, "PT202_psig", t0, 170, 1210, p_atm) < 7           # sealed hold
    hold = [i for (i, t) in enumerate(D["t_s"]) if 170 <= t - t0 <= 1210 && isfinite(D["TH_tank_C"][i])]
    # the measurement itself: thermistors against the saturation temperature of the logged pressure
    @test maximum(abs(D["TH_tank_C"][i] - (F.Tsat(D["PT202_psig"][i] * F.PSI + p_atm) - 273.15)) for i in hold) < 2
    liquid(i) = sim.hist["TANK.T"][clamp(searchsortedfirst(sim.t_hist, D["t_s"][i] - t0), 1, length(sim.t_hist))] - 273.15
    @test abs(sum(liquid(i) - D["TH_tank_C"][i] for i in hold) / length(hold)) < 2
    @test closes(sim)
end

@testset "HF4 replay (calibration regression)" begin
    _, D = F.read_csv(joinpath(@__DIR__, "data", "hf4_fill_10hz.csv"))
    s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "hf4.toml"))
    t0 = 9020.0
    s.schedule = F.schedule_from_log(D["t_s"], D, Dict("AV201" => ("AV201", "open", "close"), "AV202" => ("AV202", "open", "close"),
        "AV203" => ("AV203", "open", "close"), "AV204" => ("AV204", "open", "close"), "heater" => ("TANK", "on", "off")); t0)
    sim = F.Sim(s)
    wall = @elapsed F.run_until!(sim, 990.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    @test wall < 99                                              # at least 10× faster than real time
    function rms(key, col, a, b)
        e2, n = 0.0, 0
        for (i, td) in enumerate(D["t_s"])
            a <= td <= b || continue
            k = clamp(searchsortedfirst(sim.t_hist, td - t0), 1, length(sim.t_hist))
            e2 += ((sim.hist[key][k] - 101325) / F.PSI - D[col][i])^2
            n += 1
        end
        sqrt(e2 / n)
    end
    @test rms("TANK.P", "PT202_psig", 9160, 9256) < 15           # stall
    @test rms("TANK.P", "PT202_psig", 9256, 9372) < 25           # two vent pulses and rebounds
    @test rms("TANK.P", "PT202_psig", 9387, 9856) < 15           # heater segment
    @test rms("PT201.P", "PT201_psig", 9100, 10000) < 20         # bottle side: cooling and re-warming
    @test 8.8 < sim.hist["TANK.m"][end] < 10.0                   # ≈ 9.3 kg loaded
    @test closes(sim)
end

@testset "HF3 replay (independent check: HF4's fitted values, nothing re-fitted)" begin
    # A different day and a different procedure (one unvented fill, AV203 open before AV201).
    # Only what was measured or is plainly different is changed: the bottle's starting pressure
    # (PT201 = 849 psig) and a cooler day (22 °C is a guess; the log has no temperature).
    _, D = F.read_csv(joinpath(@__DIR__, "data", "hf3_fill_10hz.csv"))
    s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "hf4.toml"))
    s.settings["T_amb_C"] = 22.0
    s.parts[F.find_part(s, "BOTTLE")].p["T_C"] = F.Tsat((849.0 + 14.696) * F.PSI) - 273.15
    s.parts[F.find_part(s, "TANK")].p["T_C"] = 22.0
    t0 = 4520.0
    s.schedule = F.schedule_from_log(D["t_s"], D, Dict("AV201" => ("AV201", "open", "close"), "AV203" => ("AV203", "open", "close"),
        "AV204" => ("AV204", "open", "close")); t0)
    sim = F.Sim(s)
    F.run_until!(sim, 300.0)
    @test isempty(sim.failed) && sim.net.unconverged == 0
    function rms(a, b)
        e2, n = 0.0, 0
        for (i, td) in enumerate(D["t_s"])
            a <= td <= b || continue
            k = clamp(searchsortedfirst(sim.t_hist, td - t0), 1, length(sim.t_hist))
            e2 += ((sim.hist["TANK.P"][k] - 101325) / F.PSI - D["PT202_psig"][i])^2
            n += 1
        end
        sqrt(e2 / n)
    end
    @test rms(4528, 4548) < 45                                   # first 20 s (no pre-charged hose in HF3)
    @test rms(4548, 4620) < 70                                   # fast fill: the model runs ≈ 20 s ahead
    @test rms(4620, 4814) < 12                                   # stall pressure, set by bottle cooling
    @test 8.3 < sim.hist["TANK.m"][end] < 9.5                    # load cell: 9.3 kg including ≈ 0.6 kg of hose load
    @test closes(sim)
end

end
