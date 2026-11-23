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

@testset "surface tension and boiling on a wall" begin
    for T in (200.0, 250.0, 280.0, 300.0)
        @test F.surface_tension(T) ≈ PropsSI("surface_tension", "T", T, "Q", 0, "NitrousOxide") rtol = 2e-3
    end
    @test F.surface_tension(F.N2O_TC + 1) == 0
    # the correlation is pinned at 20 kW/m²: with ΔT = q20 / h20 it must return exactly that
    P = 0.5 * F.N2O_PC
    st = F.sat_P(P)
    L = log(0.5); h20 = exp(0.3092L^3 + 1.649L^2 + 3.641L + 5.272)
    @test 25 < h20 < 40                                           # kW/m²K: a strong coefficient
    @test F.boiling_flux(P, 20 / h20, st) ≈ 20e3 rtol = 1e-9
    @test F.boiling_flux(P, 0.0, st) == 0 && F.boiling_flux(P, -3.0, st) == 0
    @test F.boiling_flux(P, 2 * 20 / h20, st) / 20e3 ≈ 2^(1 / (0.3 * 0.5^0.3)) rtol = 1e-9     # q ∝ ΔT^(1/(1-n))
    @test 2e5 < F.boiling_flux(P, 50.0, st) < 2e6                 # capped at the critical heat flux, a few hundred kW/m²
    @test F.boiling_flux(0.98 * F.N2O_PC, 50.0, F.sat_P(0.98 * F.N2O_PC)) < F.boiling_flux(P, 50.0, st)   # which falls away towards the critical point
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
    # A burst of C's rows stamped with A's clock, in the 10 s where C was otherwise silent: they
    # stay where they are, are counted, and do not spoil C's offset. There are 300 of them against
    # C's 220 ordinary rows, so a median over all of C's rows would give an offset near zero.
    burst = [(tt, "$(round(1000 + tt, digits = 3)),C,,,$(round(20 - 0.1tt, digits = 3)),$(round(20.4 - 0.1tt, digits = 3))") for tt in 80.01 .+ (0:299) ./ 30]
    write(joinpath(dir, "burst.csv"), "device_timestamp,source,PT [PSI],solenoid_V1,TH1 [C],TH2 [C]\n" * join(last.(sort!(vcat(rows, burst); by = first)), "\n") * "\n")
    write(manifest, replace(body, "log.csv" => "burst.csv"))
    log = F.load_log(manifest)
    @test abs(log.offsets["B"] - 100.18) < 0.1 && abs(log.offsets["C"] - 1869.08) < 0.1 && log.offsets["A"] == 0
    @test log.warnings == ["300 rows of C were already on the clock of A and were not shifted"]
    @test abs(at("TH", 1085.0) - (20.2 - 8.5)) < 0.15 && abs(at("TH", 1050.0) - (20.2 - 5)) < 0.15 && log.cols["t_s"][end] == 1120.0
    write(manifest, replace(body, "log.csv" => "burst.csv", "clock_offsets = \"auto\"" => "clock_offsets = { B = 100.18, C = 1869.08 }"))
    log = F.load_log(manifest)                                   # the same with given offsets
    @test length(log.warnings) == 1 && abs(at("TH", 1085.0) - (20.2 - 8.5)) < 0.15 && log.cols["t_s"][end] == 1120.0
    write(manifest, replace(body, "log.csv" => "burst.csv", "sample_dt" => "already_aligned_s = 0\nsample_dt"))
    log = F.load_log(manifest)                                   # the rule switched off: the burst drags C's median to zero
    @test !any(occursin("already", w) for w in log.warnings) && abs(log.offsets["C"]) < 1
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
    # a run that would take too much work is given up on, by count and not by the clock
    out = F.study(study(top * "work_limit = 100\n" * rest); quiet = true)
    head, rows = table(out)
    @test occursin("work limit", rows[1][findfirst(==("failed"), head)]) && parse(Float64, rows[1][findfirst(==("t_end"), head)]) < 20
    @test_throws ErrorException F.read_study(study(top * rest * "[vary]\n\"TANK.volume\" = [1, 2]\n"))      # no such parameter
    @test_throws ErrorException F.read_study(study(top * "sample = 1\n" * rest))                                # no such key
    @test_throws ErrorException F.read_study(study(top * replace(rest, "\"final\"" => "\"last\"")))             # no such metric
    @test_throws ErrorException F.read_study(study(top * rest * "oops = 1\n"))                                  # a key that landed in [metrics]
    # scored against a log: a one-run study gives what a single replay gives
    root = joinpath(@__DIR__, "..")
    fwd(p) = replace(p, '\\' => '/')                              # a Windows path's backslashes are escapes inside a TOML string
    one = study("""
        scene = "$(fwd(joinpath(root, "scenes", "hf3cf.toml")))"
        replay = "$(fwd(joinpath(root, "test", "data", "hf3cf_fill_10hz.csv")))"
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
    val(name) = parse(Float64, rows[1][findfirst(==(name), head)])
    @test val("rms.PT202_psig.sealed") ≈ direct rtol = 1e-6
    @test 2 < direct < 3
    # the score has no units: each RMS over that reading's sigma (5 psi by default for a pressure),
    # every channel and phase counted once, the second channel here at half weight
    ratios = [val("rms.PT202_psig.sealed") / 5, val("rms.PT202_psig.fill") / 5, 0.5 * val("rms.PT201_psig.sealed") / 5, 0.5 * val("rms.PT201_psig.fill") / 5]
    @test val("score") ≈ sqrt(sum(abs2, ratios) / 4) rtol = 1e-5
    # and it is what the residual vector adds up to, so a fitter minimises the same thing
    rms1, res1 = F.score_run(st1, sim)
    @test sqrt(sum(abs2, res1) / 4) ≈ F.study_score(F.score_run(st1, sim)) rtol = 1e-9
    # a sigma given for a channel changes its share of the score and nothing else
    st2 = F.read_study(study(replace(read(one, String), "weight = 0.5" => "weight = 0.5, sigma = 20")))
    @test [c[3] for c in st2.channels] == [5.0, 20.0] && [c[3] for c in st1.channels] == [5.0, 5.0]
    @test F.default_sigma(:temperature) == 1.0 && F.default_sigma(:mass) == 0.25
    # on two threads the results are the same to the last digit
    script = "using N2OFillSim; N2OFillSim.study(ARGS[1]; quiet = true)"
    two = study(top * rest * vary_grid * "# threaded\n")
    run(`$(Base.julia_cmd()) -t 2 --project=$(Base.active_project()) -e $script $two`)
    strip_wall(path) = sort([join(split(l, ',')[1:end-1], ',') for l in readlines(path)[2:end]])
    @test strip_wall(splitext(two)[1] * "_results/results.csv") == strip_wall(splitext(grid)[1] * "_results/results.csv")
end

@testset "fitting: recovers known values from a made-up log, and says what cannot be told apart" begin
    # the pieces
    x, fx, n = F.nelder_mead(v -> (v[1] - 0.3)^2 + 10 * (v[2] - 0.8)^2 + (v[3] - 1.2)^2, [0.5, 0.5, 0.5]; max_evals = 600)
    @test abs(x[1] - 0.3) < 0.01 && abs(x[2] - 0.8) < 0.01 && x[3] == 1.0             # the third wants to leave the cube and is held at its edge
    # looking ahead with spare cores evaluates more points but takes the same path to the same answer
    q(v) = (v[1] - 0.3)^2 + 10 * (v[2] - 0.8)^2 + (v[3] - 1.2)^2
    asked = Ref(0)
    xb, fb, nb = F.nelder_mead(q, [0.5, 0.5, 0.5]; max_evals = 600, batch = pts -> (asked[] += length(pts); q.(pts)))
    @test xb == x && fb == fx && nb == n && asked[] > n
    # Levenberg-Marquardt on a curve fit with known answer (a, b) = (0.3, 0.7) for y = exp(-4 a t) + b t, and a bound in the way
    ts = 0:0.1:2
    resid(v) = [exp(-4 * v[1] * t) + v[2] * t - (exp(-4 * 0.3 * t) + 0.7 * t) for t in ts]
    calls = Ref(0)
    xl, Sl, nl = F.levenberg_marquardt(pts -> (calls[] += length(pts); resid.(pts)), [0.9, 0.1])
    @test abs(xl[1] - 0.3) < 1e-3 && abs(xl[2] - 0.7) < 1e-3 && Sl < 1e-8 && nl == calls[] && nl < 120
    xl, = F.levenberg_marquardt(pts -> [vcat(resid(p), 30 * (p[2] - 1.4)) for p in pts], [0.9, 0.1])    # a strong extra pull on b beyond its range
    @test xl[2] == 1.0 && abs(xl[1] - 0.636) < 0.01                                    # held at the bound; the other at its best value given that (0.636 by grid search)
    xl, = F.levenberg_marquardt(pts -> [p[1] > 0.6 ? nothing : resid(p) for p in pts], [0.5, 0.1])       # the model "fails" in part of the range
    @test abs(xl[1] - 0.3) < 1e-2
    J = [1.0 2.0 0.0; 2.0 4.0 0.0; 0.5 1.0 0.0; 0.0 0.0 0.0]                         # column 2 is twice column 1; column 3 does nothing
    _, _, said = F.identifiability(J .+ [0 0 0; 0 0 0; 0 0 0; 1e-9 0 0], ["a", "b", "c"])
    @test any(occursin("c has no measurable effect", t) for t in said) && any(occursin("a", t) && occursin("b", t) for t in said)
    _, corr, said = F.identifiability([1.0 0.0; 0.0 1.0; 1.0 -1.0], ["a", "b"])
    @test isempty(said) && abs(corr[1, 2]) < 0.9
    # a log made by the model itself with known values, plus half a psi of noise
    dir = mktempdir()
    truth = rig(bottle_T = 27.0)
    F.apply_overrides!(truth, ["FILL.d_mm=4"])
    truth.schedule = [F.Action(1.0, "FILL", "open", 0.0), F.Action(1.0, "VENT", "open", 0.0), F.Action(15.0, "VENT", "close", 0.0)]
    sim = F.Sim(truth)
    F.run_until!(sim, 30.0)
    g = F._Lcg(UInt64(7))
    open(joinpath(dir, "log.csv"), "w") do io
        println(io, "t_s,PT_psig,FILL,VENT")
        for t in 0:0.5:30
            println(io, t, ",", round((F.signal_at(sim, "TANK.P", t) - 101325) / F.PSI + (F._next!(g) - 0.5), digits = 3), ",", t >= 1 ? 1 : 0, ",", 1 <= t < 15 ? 1 : 0)
        end
    end
    F.save_scene(joinpath(dir, "rig.toml"), rig())                # the scene to be fitted has neither value right
    head = """
        name = "made-up"
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        [dataset.score]
        channels = [{ log = "PT_psig", model = "TANK.P", sigma = 0.5 }]
        phases = { venting = [1, 15], closed = [15, 30] }
        [[parameter]]
        name = "FILL.d_mm"
        from = 2
        to = 8
        start = 6
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 22
        """
    write(joinpath(dir, "fit.toml"), head)
    r = F.fit(joinpath(dir, "fit.toml"); write_scenes = true, quiet = true)
    @test abs(r.values["FILL.d_mm"] - 4) < 0.2 && abs(r.values["BOTTLE.T_C"] - 27) < 0.5
    write(joinpath(dir, "fit_nm.toml"), "method = \"nelder-mead\"\n" * head)
    rn = F.fit(joinpath(dir, "fit_nm.toml"); quiet = true)
    @test abs(rn.values["FILL.d_mm"] - r.values["FILL.d_mm"]) < 0.1 && abs(rn.score - r.score) < 0.05      # the two methods agree
    @test r.evals < rn.evals                                      # and the default needs fewer runs
    @test r.score < 1.5 < r.start_score                           # down to about the noise put in
    @test !any(occursin("do not determine", t) for t in r.statements)   # correlated (both set the fill rate) but both recovered
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "made-up.toml"), readdir(r.dir; join = true))))
    @test record["datasets"] == ["run1"] && length(record["parameter"]) == 2 && record["dataset"][1]["channels"][1]["sigma"] == 0.5 && haskey(record, "model_commit")
    fitted = F.load_scene(joinpath(r.dir, "run1.toml"))
    @test abs(fitted.parts[F.find_part(fitted, "FILL")].p["d_mm"] - r.values["FILL.d_mm"]) < 1e-3
    @test occursin("[cal ", readline(joinpath(r.dir, "run1.toml"))) && isfile(joinpath(r.dir, "report.txt")) && countlines(joinpath(r.dir, "evals.csv")) == r.evals + 1
    # a valve's bore and its discharge coefficient only ever act together: the fit must say so, not pick a pair
    write(joinpath(dir, "fit2.toml"), head * "[[parameter]]\nname = \"FILL.Cd\"\nfrom = 0.2\nto = 0.95\n" * "max_evals = 30\n")
    @test_throws ErrorException F.read_fit(joinpath(dir, "fit2.toml"))            # max_evals landed inside the last table: an unknown key
    write(joinpath(dir, "fit2.toml"), "max_evals = 30\n" * head * "[[parameter]]\nname = \"FILL.Cd\"\nfrom = 0.2\nto = 0.95\n")
    r2 = F.fit(joinpath(dir, "fit2.toml"); quiet = true)
    @test any(occursin("FILL.d_mm", t) && occursin("FILL.Cd", t) for t in r2.statements)
    # the same parameter fitted separately per dataset, and mistakes caught before anything runs
    two = replace(head, "[[parameter]]\nname = \"BOTTLE.T_C\"" => "[[dataset]]\nname = \"run2\"\nscene = \"rig.toml\"\nreplay = \"log.csv\"\n[dataset.score]\nchannels = [{ log = \"PT_psig\", model = \"TANK.P\" }]\n[[parameter]]\ndatasets = [\"run1\"]\nname = \"BOTTLE.T_C\"")
    write(joinpath(dir, "fit3.toml"), two * "[[parameter]]\nname = \"BOTTLE.T_C\"\nfrom = 15\nto = 32\ndatasets = [\"run2\"]\n")
    ft = F.read_fit(joinpath(dir, "fit3.toml"))
    @test [p.label for p in ft.params] == ["FILL.d_mm", "BOTTLE.T_C@run1", "BOTTLE.T_C@run2"] && ft.params[2].datasets == [1] && ft.params[1].datasets == [1, 2]
    write(joinpath(dir, "bad.toml"), replace(head, "FILL.d_mm" => "FILL.bore"))
    @test_throws ErrorException F.read_fit(joinpath(dir, "bad.toml"))
    write(joinpath(dir, "bad.toml"), head * "[[parameter]]\nname = \"BOTTLE.T_C\"\nfrom = 1\nto = 2\n")
    @test_throws ErrorException F.read_fit(joinpath(dir, "bad.toml"))
    # a draft fit file from the [fit] tags of a scene
    draft = F.fit_skeleton(joinpath(@__DIR__, "..", "scenes", "hf4.toml"))
    @test occursin("name = \"BOTTLE.fill_kg\"", draft) && occursin("name = \"settings.h_in_mult\"", draft)
end

# A made-up log for the scoring tests: the rig with a 4 mm fill valve and a 27 °C bottle, filled
# with the vent open until 15 s and closed after. The tank starts with 2 kg in it, which makes a
# run several times quicker than starting empty. `extra(t)` is added to the logged pressure, and
# `more` is further columns as name => function of (run, time).
scoring_rig(; kw...) = rig(; tank_fill = 2.0, kw...)
function made_up_log(path; extra = t -> 0.0, more = Pair{String,Function}[], bottle_T = 27.0)
    truth = scoring_rig(; bottle_T)
    F.apply_overrides!(truth, ["FILL.d_mm=4"])
    truth.schedule = [F.Action(1.0, "FILL", "open", 0.0), F.Action(1.0, "VENT", "open", 0.0), F.Action(15.0, "VENT", "close", 0.0)]
    sim = F.Sim(truth)
    F.run_until!(sim, 30.0)
    open(path, "w") do io
        println(io, join(vcat("t_s", "PT_psig", "FILL", "VENT", first.(more)), ","))
        for t in 0:0.5:30
            p = (F.signal_at(sim, "TANK.P", t) - 101325) / F.PSI + extra(t)
            println(io, join(vcat(t, round(p, digits = 3), t >= 1 ? 1 : 0, 1 <= t < 15 ? 1 : 0, [round(f(sim, t), digits = 4) for (_, f) in more]), ","))
        end
    end
    sim
end

@testset "scoring: time bins" begin
    dir = mktempdir()
    made_up_log(joinpath(dir, "log.csv"))
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    study(score) = (p = joinpath(dir, "s$(hash(score)).toml"); write(p, "scene = \"rig.toml\"\nreplay = \"log.csv\"\n[score]\nchannels = [{ log = \"PT_psig\", model = \"TANK.P\" }]\n" * score); p)
    plain = F.read_study(study("phases = { venting = [1, 15], closed = [15, 30] }\n"))
    st = F.read_study(study("bin_s = 5\nphases = { venting = [1, 15], closed = [15, 30] }\n"))
    # a phase is cut into equal bins as near to bin_s long as fit it: three each here
    @test [g.name for g in st.scored] == ["PT_psig.closed.$k" for k in 1:3] ∪ ["PT_psig.venting.$k" for k in 1:3]
    @test F.n_scored(st) == F.n_residuals(st) == 6 && F.n_scored(plain) == 2 && F.n_residuals(plain) == 29 + 31
    sim, failed = F.simulate(st, [])
    s, s0 = F.score_run(st, sim), F.score_run(plain, sim)
    # one residual per bin: the mean of model − log over the bin, with its sign, over sigma (5 psi)
    diff(t) = (F.signal_at(sim, "TANK.P", t) - 101325) / F.PSI - st.cols["PT_psig"][findfirst(==(t), st.t_log)]
    @test s.residuals[1] ≈ sum(diff, 15:0.5:19.5) / 10 / 5 rtol = 1e-9
    @test s.residuals == last.(s.z) && first(s.z[1]) == "res.PT_psig.closed.1"
    @test F.study_score(s) ≈ sqrt(sum(abs2, s.residuals) / 6)
    # the per-phase RMS columns are the same with or without bins
    @test first.(s.rms) == ["rms.PT_psig.closed", "rms.PT_psig.venting"] && last.(s.rms) ≈ last.(s0.rms)
    # a phase can set its own bin length; the other is then scored whole
    own = F.read_study(study("phases = { venting = [1, 15], closed = { from = 15, to = 30, bin_s = 7.5 } }\n"))
    @test [g.name for g in own.scored] == ["PT_psig.closed.1", "PT_psig.closed.2", "PT_psig.venting"] && F.n_residuals(own) == 2 + 29
    # the results file has a column for every scored number
    head = split(readline(F.study(study("bin_s = 5\nphases = { venting = [1, 15], closed = [15, 30] }\n"); quiet = true)), ',')
    @test "res.PT_psig.venting.3" in head && "rms.PT_psig.venting" in head
    @test_throws ErrorException F.read_study(study("bins = 5\n"))                                           # no such key
    @test_throws ErrorException F.read_study(study("phases = { closed = { from = 15, to = 30, bin = 5 } }\n"))
    # a fit on binned residuals finds the values the log was made with
    write(joinpath(dir, "fit.toml"), """
        max_evals = 40
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        [dataset.score]
        channels = [{ log = "PT_psig", model = "TANK.P", sigma = 0.5 }]
        phases = { venting = [1, 15], closed = [15, 30] }
        bin_s = 3.5
        [[parameter]]
        name = "FILL.d_mm"
        from = 2
        to = 8
        start = 6
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 22
        """)
    r = F.fit(joinpath(dir, "fit.toml"); quiet = true)
    @test abs(r.values["FILL.d_mm"] - 4) < 0.1 && abs(r.values["BOTTLE.T_C"] - 27) < 0.3 && r.score < 0.5
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "-fit.toml"), readdir(r.dir; join = true))))
    @test length(record["dataset"][1]["scored"]) == 8
end

@testset "scoring: a model error stated for a phase; a fit is biased without it and recovers with it" begin
    dir = mktempdir()
    # The log is the model's own output at known values, except that during the venting phase it
    # reads 25 psi higher than the model ever can: a stand-in for a mechanism the model lacks.
    made_up_log(joinpath(dir, "log.csv"); extra = t -> 1 <= t < 15 ? 25.0 : 0.0)
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    study(score) = (p = joinpath(dir, "s$(hash(score)).toml"); write(p, "scene = \"rig.toml\"\nreplay = \"log.csv\"\n[score]\n" * score); p)
    one = "channels = [{ log = \"PT_psig\", model = \"TANK.P\" }]\n"
    # the stated error widens sigma for that phase only, in quadrature: sqrt(5² + 12²) = 13
    st = F.read_study(study(one * "phases = { venting = { from = 1, to = 14.9, discrepancy = 12 }, closed = [15, 30] }\n"))
    @test [g.name => g.sigma for g in st.scored] == ["PT_psig.closed" => 5.0, "PT_psig.venting" => 13.0]
    sim, = F.simulate(st, [])
    s = F.score_run(st, sim)
    @test Dict(s.z)["res.PT_psig.venting"] ≈ Dict(s.rms)["rms.PT_psig.venting"] / 13 && Dict(s.z)["res.PT_psig.closed"] ≈ Dict(s.rms)["rms.PT_psig.closed"] / 5
    # it can be given per log column (here two channels read the same column's twin)
    made_up_log(joinpath(dir, "log2.csv"); more = ["PT2_psig" => (sim, t) -> (F.signal_at(sim, "TANK.P", t) - 101325) / F.PSI])
    two = replace(read(study("channels = [{ log = \"PT_psig\", model = \"TANK.P\" }, { log = \"PT2_psig\", model = \"TANK.P\" }]\n" *
                             "phases = { venting = { from = 1, to = 14.9, discrepancy = { PT2_psig = 12 } } }\n"), String), "log.csv" => "log2.csv")
    write(joinpath(dir, "two.toml"), two)
    @test [g.sigma for g in F.read_study(joinpath(dir, "two.toml")).scored] == [5.0, 13.0]
    write(joinpath(dir, "two.toml"), replace(two, "PT2_psig = 12" => "PT3_psig = 12"))
    @test_throws ErrorException F.read_study(joinpath(dir, "two.toml"))                       # not a scored column
    @test_throws ErrorException F.read_study(study(one * "phases = { venting = { from = 1, to = 14.9, discrepancy = -1 } }\n"))
    # End to end. Two values are fitted: the bottle temperature (27 °C in the log), which the
    # closed phase alone is enough to find, and the vent orifice (2 mm), which only acts while
    # venting. A plain fit bends both to soak up the 25 psi it cannot explain.
    fit(phases) = (p = joinpath(dir, "f$(hash(phases)).toml"); write(p, """
        max_evals = 40
        jacobian_reltol = 1e-5
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        set = { "FILL.d_mm" = 4 }
        [dataset.score]
        channels = [{ log = "PT_psig", model = "TANK.P", sigma = 1 }]
        phases = { $phases }
        bin_s = 5
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 22
        [[parameter]]
        name = "VENT_OR.d_mm"
        from = 1
        to = 4
        start = 3
        """); F.fit(p; quiet = true))
    plain = fit("venting = [1, 14.9], closed = [15, 30]")
    @test abs(plain.values["BOTTLE.T_C"] - 27) > 0.5 && abs(plain.values["VENT_OR.d_mm"] - 2) > 0.3 && plain.score > 5
    # With the error stated (25 psi on the venting phase), the fit stops chasing it, and the
    # bottle temperature comes back to within 0.1 °C.
    stated = fit("venting = { from = 1, to = 14.9, discrepancy = 25 }, closed = [15, 30]")
    @test abs(stated.values["BOTTLE.T_C"] - 27) < 0.1 && stated.score < 1
    record = F.TOML.parsefile(only(filter(f -> occursin(r"\d-f\d+\.toml$", f), readdir(stated.dir; join = true))))
    @test record["dataset"][1]["scored"][end]["sigma"] ≈ sqrt(1 + 25^2) rtol = 1e-5
end

@testset "scoring: a channel whose zero is unknown, for the whole run or phase by phase" begin
    dir = mktempdir()
    # A load cell that was not tared (it reads 3.2 kg high), and that steps down by 0.8 kg as the vent valve shuts at 15 s.
    made_up_log(joinpath(dir, "log.csv"); more = ["LC_kg" => (sim, t) -> F.signal_at(sim, "TANK.m", t) + 3.2,
                                                  "LC_step_kg" => (sim, t) -> F.signal_at(sim, "TANK.m", t) + 3.2 - (t >= 15 ? 0.8 : 0.0)])
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig(bottle_T = 27.0))
    study(channel) = (p = joinpath(dir, "s$(hash(channel)).toml"); write(p, """
        scene = "rig.toml"
        replay = "log.csv"
        [set]
        "FILL.d_mm" = 4
        [score]
        channels = [{ $channel, model = "TANK.m" }]
        phases = { venting = [1, 14.9], closed = [15, 30] }
        """); p)
    run(channel) = (st = F.read_study(study(channel)); (st, F.score_run(st, F.simulate(st, [])[1])))
    # scored as it stands, the model at the log's own values looks more than 3 kg out
    st, s = run("log = \"LC_kg\"")
    @test all(>(3.1), last.(s.rms)) && isempty(s.offsets)
    # with the zero free the model lies on the log, and the offset is what the load cell read too much
    st, s = run("log = \"LC_kg\", offset = \"free\"")
    @test all(<(0.01), last.(s.rms)) && first.(s.offsets) == ["offset.LC_kg"] && isapprox(last(s.offsets[1]), -3.2; atol = 0.01)
    # one offset for the run cannot follow a zero that jumps: half the step is left in each phase
    st, s = run("log = \"LC_step_kg\", offset = \"free\"")
    @test all(x -> 0.3 < x < 0.5, last.(s.rms))
    # an offset for each phase can, and both are reported
    st, s = run("log = \"LC_step_kg\", offset = \"per_phase\"")
    @test all(<(0.01), last.(s.rms)) && first.(s.offsets) == ["offset.LC_step_kg.closed", "offset.LC_step_kg.venting"]
    @test isapprox(last(s.offsets[1]), -2.4; atol = 0.01) && isapprox(last(s.offsets[2]), -3.2; atol = 0.01)
    # but then only changes inside a phase are scored: a run that loads the wrong mass by the end of
    # the venting phase and then tracks the log is not told apart by the closed phase
    wrong = F.simulate(st, [])[1]
    wrong.hist["TANK.m"] .+= [t >= 15 ? 1.0 : t / 15 for t in wrong.t_hist]
    sw = Dict(F.score_run(st, wrong).rms)
    @test sw["rms.LC_step_kg.closed"] < 0.01 && sw["rms.LC_step_kg.venting"] > 0.2
    # the offsets are columns of the results
    lines = readlines(F.study(study("log = \"LC_step_kg\", offset = \"per_phase\""); quiet = true))
    @test isapprox(parse(Float64, split(lines[2], ',')[findfirst(==("offset.LC_step_kg.venting"), split(lines[1], ','))]), -3.2; atol = 0.01)
    @test_throws ErrorException F.read_study(study("log = \"LC_kg\", offset = \"yes\""))
    @test_throws ErrorException F.read_study(study("log = \"LC_kg\", ofset = \"free\""))
    @test_throws ErrorException F.read_study(study("log = \"LC_kg\", correct = { by = \"PT_psig\", per_unit = 0.004 }"))
    # a fit on the mass alone, its zero unknown in each phase: the bottle temperature is found from how the mass changes
    write(joinpath(dir, "fit.toml"), """
        max_evals = 25
        jacobian_reltol = 1e-5
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        set = { "FILL.d_mm" = 4 }
        [dataset.score]
        channels = [{ log = "LC_step_kg", model = "TANK.m", sigma = 0.01, offset = "per_phase" }]
        phases = { venting = [1, 14.9], closed = [15, 30] }
        bin_s = 5
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 22
        """)
    r = F.fit(joinpath(dir, "fit.toml"); quiet = true)
    @test abs(r.values["BOTTLE.T_C"] - 27) < 0.2
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "-fit.toml"), readdir(r.dir; join = true))))["dataset"][1]
    @test isapprox(record["offset"]["LC_step_kg.closed"], -2.4; atol = 0.02) && record["channels"][1]["offset"] == "per_phase"
    @test occursin("offset.LC_step_kg.venting", read(joinpath(r.dir, "report.txt"), String))
end

@testset "scoring: targets (a feature of the log against the same feature of the model), slope, windows" begin
    # the pieces: a slope by least squares on unevenly spaced samples, and a window on any metric
    @test F._metric([0.0, 1, 1.5, 4, 10], 3 .+ 2 .* [0.0, 1, 1.5, 4, 10], Dict("fn" => "slope"), 0.0) ≈ 2
    sim = F.Sim(rig())
    empty!(sim.t_hist); append!(sim.t_hist, [0.0, 10.0, 20.0]); sim.hist["X"] = [0.0, 100.0, 100.0]
    m(fn; kw...) = F.metric(sim, Dict{String,Any}("fn" => fn, "signal" => "X", (String(k) => v for (k, v) in kw)...), 0.0)
    @test m("max") == 100 && m("max"; from = 2, to = 5) ≈ 50 && m("min"; from = 5) ≈ 50 && m("final"; to = 5) ≈ 50       # the window's ends are interpolated
    @test m("slope"; from = 0, to = 10) ≈ 10 && m("slope"; from = 10, to = 20) ≈ 0 && 0 < m("slope") < 10
    @test isnan(m("max"; from = 30, to = 40))                                                    # a window the run never reached
    bad = joinpath(mktempdir(), "s.toml")                                                        # a unit shifts a pressure by one atmosphere: no good for a slope
    F.save_scene(joinpath(dirname(bad), "rig.toml"), rig())
    write(bad, "scene = \"rig.toml\"\n[metrics]\nrise = { fn = \"slope\", signal = \"TANK.P\", unit = \"psig\" }\n")
    @test_throws ErrorException F.read_study(bad)
    dir = mktempdir()
    made_up_log(joinpath(dir, "log.csv"); more = ["RAMP_psig" => (sim, t) -> 100 + 2t])
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    study(targets; set = "") = (p = joinpath(dir, "s$(hash((targets, set))).toml"); write(p, """
        scene = "rig.toml"
        replay = "log.csv"
        [set]
        $set
        [score]
        channels = [{ log = "PT_psig", model = "TANK.P" }]
        phases = { venting = [1, 15], closed = [15, 30] }
        [targets]
        $targets
        """); p)
    body = """
        vent_min = { fn = "min", log = "PT_psig", model = "TANK.P", phase = "venting", tol = 5 }
        rise = { fn = "slope", log = "PT_psig", model = "TANK.P", from = "VENT close", to = 30, tol = 0.5 }
        ramp = { fn = "slope", log = "RAMP_psig", model = "TANK.P", from = 15, to = 30, tol = 0.5 }
        """
    st = F.read_study(study(body))
    @test [g.name for g in st.targets] == ["ramp", "rise", "vent_min"] && st.targets[1].value ≈ 2 && st.targets[2].m["from"] == 15
    @test F.n_scored(st) == 2 + 3 && F.n_residuals(st) == 29 + 31 + 3
    run, = F.simulate(st, [])
    s = F.score_run(st, run)
    T, Z = Dict(s.targets), Dict(s.z)
    psig(t) = (F.signal_at(run, "TANK.P", t) - 101325) / F.PSI
    @test T["target.vent_min.model"] ≈ minimum(psig, 1:0.5:15) && T["target.vent_min.log"] ≈ minimum(st.cols["PT_psig"][3:31])
    @test Z["res.target.vent_min"] ≈ (T["target.vent_min.model"] - T["target.vent_min.log"]) / 5
    @test s.residuals[end-2:end] == [Z["res.target.ramp"], Z["res.target.rise"], Z["res.target.vent_min"]]
    @test F.study_score(s) ≈ sqrt(sum(abs2, s.residuals) / 5)
    # a slope is in the log's units per second on both sides (psi/s here, not Pa/s)
    @test T["target.rise.model"] ≈ F._metric(collect(15:0.5:30), psig.(15:0.5:30), Dict("fn" => "slope"), 0.0) && 0.1 < T["target.rise.model"] < 50
    # with the scene at the log's own values, every feature agrees (to a small part of its tolerance:
    # the log was sampled from a run recorded at other instants)
    exact = F.read_study(study(body; set = "\"FILL.d_mm\" = 4\n\"BOTTLE.T_C\" = 27"))
    se = F.score_run(exact, F.simulate(exact, [])[1])
    @test abs(Dict(se.z)["res.target.rise"]) < 0.05 && abs(Dict(se.z)["res.target.vent_min"]) < 0.05
    # a feature the model does not have at all counts as 100 tolerances out: the log passes 700 psig, a run with a 15 °C bottle never does
    cross = "t700 = { fn = \"first_crossing\", log = \"PT_psig\", model = \"TANK.P\", value = 700, tol = 1 }\n"
    st7 = F.read_study(study(cross; set = "\"BOTTLE.T_C\" = 15"))
    @test 1 < st7.targets[1].value < 30 && Dict(F.score_run(st7, F.simulate(st7, [])[1]).z)["res.target.t700"] == 100
    # the results file has both values and the scored number
    head = split(readline(F.study(study(body); quiet = true)), ',')
    @test all(in(head), ["res.target.rise", "target.rise.model", "target.rise.log"])
    for bad in (replace(cross, "value = 700" => "value = 5000"),                                # the log never gets there
                replace(cross, "tol = 1" => "tol = 1, unit = \"psig\""), replace(cross, ", tol = 1" => ""), replace(cross, "first_crossing" => "first"),
                replace(cross, "PT_psig" => "PT9_psig"), replace(cross, "tol = 1" => "tol = 1, phase = \"hold\""), replace(cross, "tol = 1" => "tolerance = 1"))
        @test_throws ErrorException F.read_study(study(bad))
    end
    # in a fit, targets are scored with the residuals and written to the record
    write(joinpath(dir, "fit.toml"), """
        max_evals = 20
        jacobian_reltol = 1e-5
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        set = { "FILL.d_mm" = 4 }
        [dataset.targets]
        rise = { fn = "slope", log = "PT_psig", model = "TANK.P", from = 15, to = 30, tol = 0.05 }
        end_psig = { fn = "final", log = "PT_psig", model = "TANK.P", tol = 1 }
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 22
        """)
    r = F.fit(joinpath(dir, "fit.toml"); quiet = true)
    @test abs(r.values["BOTTLE.T_C"] - 27) < 0.2
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "-fit.toml"), readdir(r.dir; join = true))))["dataset"][1]
    @test Set(t["name"] for t in record["targets"]) == Set(["rise", "end_psig"]) && abs(record["target"]["rise.model"] - record["target"]["rise.log"]) < 0.05
end

@testset "screening: which parameters matter (elementary effects)" begin
    dir = mktempdir()
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    study(top, vary) = (p = joinpath(dir, "s$(hash((top, vary))).toml"); write(p, """
        scene = "rig.toml"
        until = 10
        mode = "morris"
        $top
        [vary]
        $vary
        [[case]]
        name = "fill"
        actions = [{ t = 1.0, part = "FILL", action = "open" }]
        [metrics]
        loaded_kg = { fn = "final", signal = "TANK.m" }
        peak_psig = { fn = "max", signal = "TANK.P", unit = "psig" }
        """); p)
    # the design: n trajectories of k + 1 runs on four levels, each run moving one parameter by two thirds of its range
    unit = "a = { from = 0, to = 1 }\nb = { from = 0, to = 1 }\nc = { from = 0, to = 1 }\n"
    vary = replace(unit, "a =" => "\"BOTTLE.T_C\" =", "b =" => "\"FILL.Cd\" =", "c =" => "\"VENT_OR.Cd\" =")
    st = F.read_study(study("n = 5\nseed = 2", vary))
    X = [Dict(r.overrides)[name] for r in st.runs, name in st.varied]
    @test size(X) == (5 * 4, 3) && all(x -> any(l -> abs(x - l) < 1e-12, (0, 1 / 3, 2 / 3, 1)), X) && allunique(r.id for r in st.runs)
    steps = [X[a+1, :] .- X[a, :] for start in 1:4:20 for a in start:start+2]
    @test all(d -> count(!=(0), d) == 1 && abs(maximum(abs, d) - 2 / 3) < 1e-12, steps)
    @test all(start -> sort([findfirst(!=(0), X[a+1, :] .- X[a, :]) for a in start:start+2]) == [1, 2, 3], 1:4:20)      # every parameter once per trajectory
    @test [r.id for r in F.read_study(study("n = 5\nseed = 2", vary)).runs] == [r.id for r in st.runs]                  # a seed repeats
    @test [r.id for r in F.read_study(study("n = 5\nseed = 3", vary)).runs] != [r.id for r in st.runs]
    # the arithmetic, on an output that is a straight line in each parameter: the effect is its coefficient, with no spread
    ee = F.elementary_effects(X, [3 * X[i, 1] - 2 * X[i, 2] for i in 1:20])
    @test [e.mu_star for e in ee] ≈ [3, 2, 0] && all(e -> e.sigma < 1e-12 && e.n == 5, ee)
    # an interaction shows as spread; a failed run costs the effects on either side of it and nothing else
    ee = F.elementary_effects(X, [X[i, 1] * X[i, 2] for i in 1:20])
    @test ee[1].sigma > 0.1 && ee[3].mu_star == 0
    y = [3 * X[i, 1] for i in 1:20]; y[2] = NaN
    @test sum(e.n for e in F.elementary_effects(X, y)) == 13 && F.elementary_effects(X, y)[1].mu_star ≈ 3
    # on the rig: the bottle's temperature and the fill valve matter to what is loaded; an orifice behind a shut vent valve does not
    file = study("n = 4", "\"BOTTLE.T_C\" = { from = 15, to = 30 }\n\"FILL.d_mm\" = { from = 2, to = 6, log = true }\n\"VENT_OR.d_mm\" = { from = 1, to = 3 }\n")
    out = F.study(file; quiet = true)
    @test countlines(out) == 1 + 4 * 4
    rows = split.(readlines(joinpath(dirname(out), "screening.csv")), ',')
    @test rows[1] == ["parameter", "output", "mu_star", "sigma", "n"] && length(rows) == 1 + 3 * 2
    mu(p, o) = parse(Float64, only(r for r in rows if r[1] == p && r[2] == o)[3])
    @test mu("FILL.d_mm", "loaded_kg") > 0.5 && mu("BOTTLE.T_C", "peak_psig") > 20 && mu("VENT_OR.d_mm", "loaded_kg") < 1e-6
    @test all(r -> r[5] == "4", rows[2:end])
    # resumed, the screening still counts every run
    lines = readlines(out)
    write(out, join(lines[1:6], "\n") * "\n")
    F.study(file; resume = true, quiet = true)
    @test split.(readlines(joinpath(dirname(out), "screening.csv")), ',') == rows
    # a run that stopped partway has numbers in its row, but the steps that touch it are left out
    lines = readlines(out)
    cells = split(lines[2], ','); cells[findfirst(==("failed"), split(lines[1], ','))] = "stopped at the work limit"
    write(out, join(vcat(lines[1], join(cells, ','), lines[3:end]), "\n") * "\n")
    F.study(file; resume = true, quiet = true)
    @test sum(parse(Int, r[5]) for r in split.(readlines(joinpath(dirname(out), "screening.csv")), ',')[2:4]) in (10, 11)      # one or two of the twelve steps
    @test_throws ErrorException F.read_study(study("", "\"BOTTLE.T_C\" = [20]\n"))                                       # a single value is not a range
    @test_throws ErrorException F.read_study(study("", vary * "[[case]]\nname = \"other\"\n"))                           # one case at most
    # against a log, the score and every scored number are screened too, and the ranking is printed
    made_up_log(joinpath(dir, "log.csv"))
    scored = joinpath(dir, "scored.toml")
    write(scored, """
        scene = "rig.toml"
        replay = "log.csv"
        mode = "morris"
        n = 2
        [vary]
        "BOTTLE.T_C" = { from = 15, to = 30 }
        "VENT_OR.d_mm" = { from = 1, to = 3 }
        [score]
        channels = [{ log = "PT_psig", model = "TANK.P" }]
        phases = { venting = [1, 15], closed = [15, 30] }
        """)
    said = mktemp() do path, io
        redirect_stdout(() -> F.study(scored), io); close(io); read(path, String)
    end
    @test occursin("Screening:", said) && occursin("BOTTLE.T_C", said)
    rows = split.(readlines(joinpath(dir, "scored_results", "screening.csv")), ',')
    @test Set(r[2] for r in rows[2:end]) == Set(["score", "res.PT_psig.closed", "res.PT_psig.venting"]) && all(r -> parse(Float64, r[3]) > 0, rows[2:end])
end

@testset "fitting: a log held out and predicted" begin
    dir = mktempdir()
    # two tests on the same hardware (a 4 mm fill valve), with the bottle at 27 °C in one and 22 °C in the other
    made_up_log(joinpath(dir, "a.csv"))
    made_up_log(joinpath(dir, "b.csv"); bottle_T = 22.0)
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    dataset(name) = """
        [[dataset]]
        name = "$name"
        scene = "rig.toml"
        replay = "$name.csv"
        [dataset.score]
        channels = [{ log = "PT_psig", model = "TANK.P", sigma = 0.5 }]
        phases = { venting = [1, 15], closed = [15, 30] }
        bin_s = 5
        """
    bottle(name, start) = "[[parameter]]\nname = \"BOTTLE.T_C\"\nfrom = 15\nto = 32\nstart = $start\ndatasets = [\"$name\"]\n"
    file = joinpath(dir, "fit.toml")
    write(file, "max_evals = 30\njacobian_reltol = 1e-5\n" * dataset("a") * dataset("b") *
                "[[parameter]]\nname = \"FILL.d_mm\"\nfrom = 2\nto = 8\nstart = 6\n" * bottle("a", 24) * bottle("b", 22))
    ft = F.read_fit(file; holdout = "b")
    @test ft.names == ["a"] && [p.label for p in ft.params] == ["FILL.d_mm", "BOTTLE.T_C@a"] && ft.held_shared == [1] && ft.held_fixed == ["BOTTLE.T_C" => 22.0]
    @test_throws ErrorException F.read_fit(file; holdout = "c")
    r = redirect_stdout(() -> F.main(["fit", file, "--holdout", "b", "--write-scenes"]), devnull)
    # b took no part, and its own bottle temperature was never fitted; the valve found from a predicts it
    @test abs(r.values["FILL.d_mm"] - 4) < 0.1 && !haskey(r.values, "BOTTLE.T_C@b")
    @test r.holdout.name == "b" && r.holdout.score < 1 && r.holdout.start_score > 5 * r.holdout.score && isapprox(r.scores["a"], r.score; rtol = 0.05)
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "-fit.toml"), readdir(r.dir; join = true))))
    @test record["datasets"] == ["a"] && record["identifiability"] isa Vector && record["holdout"]["name"] == "b" && isapprox(record["holdout"]["score"], r.holdout.score; rtol = 1e-4)
    @test haskey(record["holdout"]["readings"], "rms.PT_psig.closed") && isapprox(record["dataset"][1]["score"], r.scores["a"]; rtol = 1e-4)
    @test occursin("Held out: b", read(joinpath(r.dir, "report.txt"), String)) && isfile(joinpath(r.dir, "a.toml")) && !isfile(joinpath(r.dir, "b.toml"))
    # held out the other way round, a's bottle stays at its start value (24, not the 27 of its log), and the prediction shows it
    r2 = F.fit(file; holdout = "a", quiet = true)
    @test r2.holdout.score > 3 * r.holdout.score
end

@testset "fitting: the profile of a parameter, and the collinearity of the set" begin
    # the collinearity index: 1 for parameters that act separately, large for a set that cancels
    @test F.collinearity([1.0 0.0; 0.0 2.0; 0.0 0.0]) ≈ 1
    @test F.collinearity([1.0 0.0; 0.0 1.0; 1.0 -1.0]) ≈ sqrt(2)
    @test F.collinearity([1.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 0.0]) ≈ 1                  # a column with no effect is left out
    # three parameters, no two of them alike, of which the third does what the first two do together
    J = [1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 -1.0 0.01]
    _, corr, said = F.identifiability(J, ["a", "b", "c"])
    @test F.collinearity(J) > 20 && any(occursin("nearly collinear", t) for t in said) && F.collinearity(J[:, 1:2]) < 2
    dir = mktempdir()
    made_up_log(joinpath(dir, "log.csv"))
    F.save_scene(joinpath(dir, "rig.toml"), scoring_rig())
    file = joinpath(dir, "fit.toml")
    # the bottle temperature is what the log is about; the tank heater is never switched on, so its power changes nothing
    write(file, """
        max_evals = 20
        jacobian_reltol = 1e-5
        [[dataset]]
        name = "run1"
        scene = "rig.toml"
        replay = "log.csv"
        set = { "FILL.d_mm" = 4 }
        [dataset.score]
        channels = [{ log = "PT_psig", model = "TANK.P", sigma = 0.5 }]
        phases = { venting = [1, 15], closed = [15, 30] }
        bin_s = 5
        [[parameter]]
        name = "BOTTLE.T_C"
        from = 15
        to = 32
        start = 24
        [[parameter]]
        name = "TANK.heater_W"
        from = 0
        to = 200
        start = 50
        """)
    @test_throws ErrorException F.fit(file; profile = "BOTTLE.T", quiet = true)                # no such parameter
    r = redirect_stdout(() -> F.main(["fit", file, "--profile", "BOTTLE.T_C", "--profile-n", "3"]), devnull)
    @test abs(r.values["BOTTLE.T_C"] - 27) < 0.2 && r.collinearity ≈ 1
    # held at 15, 23.5 and 32 °C the score cannot be brought down by the other parameter: the data constrain it.
    # The fitted value is a point of the profile too, and its lowest.
    @test r.profile.values[[1, 2, 4]] == [15.0, 23.5, 32.0] && r.profile.values[3] ≈ r.values["BOTTLE.T_C"] && argmin(r.profile.scores) == 3
    @test !r.profile.flat && all(>(10 * r.score), r.profile.scores[[1, 2, 4]]) && occursin("constrain", r.profile.verdict)
    rows = split.(readlines(joinpath(r.dir, "profile_BOTTLE.T_C.csv")), ',')
    @test rows[1] == ["value", "score", "TANK.heater_W"] && length(rows) == 5 && isapprox(parse(Float64, rows[3][2]), r.profile.scores[2]; rtol = 1e-5)
    record = F.TOML.parsefile(only(filter(f -> endswith(f, "-fit.toml"), readdir(r.dir; join = true))))
    @test haskey(record, "identifiability") && record["profile"]["parameter"] == "BOTTLE.T_C" && record["profile"]["flat"] == false && record["collinearity"] ≈ 1
    @test occursin("Profile of BOTTLE.T_C", read(joinpath(r.dir, "report.txt"), String)) && occursin("Collinearity index", read(joinpath(r.dir, "report.txt"), String))
    @test countlines(joinpath(r.dir, "evals.csv")) == r.evals + 1                              # the profile's runs are not part of the search's log
    # the heater's profile is flat from one bound to the other: the fit says so instead of giving it a value
    r = F.fit(file; profile = "TANK.heater_W", profile_n = 3, quiet = true)
    @test r.profile.flat && occursin("not determined by these data", r.profile.verdict) && maximum(r.profile.scores) < 1.05 * r.score
    @test any(occursin("TANK.heater_W has no measurable effect", t) for t in r.statements)
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

@testset "two-node wall, fill-jet and boiling heat transfer (all off unless asked for)" begin
    function run(sets; until = 60.0)
        s = rig(bottle_T = 28.0, tank_T = 28.0)
        F.apply_overrides!(s, sets)
        sim = F.Sim(s)
        F.act!(sim, "FILL", "open"); F.act!(sim, "VENT", "open")
        F.run_until!(sim, until)
        sim
    end
    last_(sim, k) = sim.hist[k][end]
    plain = run(String[])
    split = run(["TANK.wall_split=1"])
    @test isempty(split.failed) && closes(split) && split.net.unconverged == 0
    @test !haskey(plain.hist, "TANK.Tw_tube") && haskey(split.hist, "TANK.Tw_tube")
    # the wall temperature reported is the mass-weighted mean of the two
    @test last_(split, "TANK.Tw") ≈ 0.5 * last_(split, "TANK.Tw_tube") + 0.5 * last_(split, "TANK.Tw_ends")
    # a thin tube with a tenth of the mass follows the cold liquid; the heavy ends lag behind it
    thin = run(["TANK.wall_split=1", "TANK.wall_tube_frac=0.1", "TANK.wall_link_W_K=1"])
    @test last_(thin, "TANK.Tw_tube") - last_(thin, "TANK.T") < 0.5 * (last_(thin, "TANK.Tw_ends") - last_(thin, "TANK.T"))
    @test closes(thin)
    # tied together by a strong link, the two halves behave as the one wall did
    tied = run(["TANK.wall_split=1", "TANK.wall_link_W_K=1e5"])
    @test abs(last_(tied, "TANK.Tw_tube") - last_(tied, "TANK.Tw_ends")) < 0.01
    @test abs(last_(tied, "TANK.P") - last_(plain, "TANK.P")) < 0.5 * F.PSI && abs(last_(tied, "TANK.m") - last_(plain, "TANK.m")) < 0.01
    # with boiling on, the tube is pinned to the liquid; ends that can only conduct slowly to their wetted face stay warm
    slow = run(["TANK.wall_split=1", "TANK.wall_tube_frac=0.4", "TANK.wall_link_W_K=1", "TANK.h_boil_mult=1", "TANK.wall_ends_W_K=5"]; until = 30.0)
    fast = run(["TANK.wall_split=1", "TANK.wall_tube_frac=0.4", "TANK.wall_link_W_K=1", "TANK.h_boil_mult=1"]; until = 30.0)
    @test last_(slow, "TANK.Tw_tube") - last_(slow, "TANK.T") < 3 && last_(slow, "TANK.Tw_ends") - last_(slow, "TANK.T") > 8
    @test last_(fast, "TANK.Tw_ends") < last_(slow, "TANK.Tw_ends") - 5 && closes(slow) && closes(fast)
    # the fill jet and boiling on a warm wall both pull the wall towards the liquid faster
    gap(sim) = last_(sim, "TANK.Tw") - last_(sim, "TANK.T")
    jet, boil = run(["TANK.h_fill_mult=1"]; until = 30.0), run(["TANK.h_boil_mult=1"]; until = 30.0)
    base = run(String[]; until = 30.0)
    @test gap(jet) < 0.8 * gap(base) && gap(boil) < 0.5 * gap(base) && gap(base) > 2
    @test closes(jet) && closes(boil) && jet.net.unconverged == 0 && boil.net.unconverged == 0
    # none of it allocates in the equations
    allocated(sim) = (du = similar(sim.u); F.rhs!(du, sim.u, sim.net, 0.0); @allocated F.rhs!(du, sim.u, sim.net, 0.0))
    @test allocated(thin) == allocated(plain) && allocated(jet) == allocated(plain)
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

@testset "a junction that no flow responds to does not send the whole junction solve to its fallback" begin
    # For a few milliseconds after the fill valve opens into an air-filled tank, the vent line
    # holds a stream the flow law returns nothing for, so that junction's row of the Jacobian is
    # empty. Left in, the empty row makes the system singular: every solve falls back to the slow
    # method, 400,000 part flows go by in the first half millisecond, and the run takes minutes.
    # The parameter values below are a random draw at which this happens, not a fit.
    names, D = F.read_csv(joinpath(@__DIR__, "data", "lc26_day4_fill_1hz.csv"))
    s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "lc26_day4.toml"))
    F.apply_overrides!(s, ["TANK.wall_split=1", "settings.h_in_mult=1", "settings.interphase_E=418.296516991197", "TANK.h_in_mult=4.522441229963602",
        "BOTTLE.h_in_W_m2K=179.9936819758818", "TANK.wall_tube_frac=0.4050508940519084", "TANK.wall_link_W_K=8.151970808176081",
        "TANK.wall_ends_W_K=27.853527908924153", "TANK.h_boil_mult=0.011543730231007783", "TANK.h_fill_mult=0.08255971377648245",
        "TANK.wall_kg=4.402227895484633", "TANK.heater_W=183.77631878833617", "VENT_ORIFICE.Cd=1.0523450146910518", "settings.h_ext_mult=3.495609596770747",
        "BOTTLE_VALVE.d_mm=4.7083741013324065", "BOTTLE.fill_kg=27.504344745061996", "settings.T_amb_C=15.35230200130054"])
    s.schedule = F.log_schedule(D, F.guess_log_columns(names, s); t0 = 21950.0)
    sim = F.Sim(s; sample_dt = 0.5)
    F.run_until!(sim, 19.1; work_limit = 400_000)                 # the fill valve opens at 19 s
    @test isempty(sim.failed) && F.simtime(sim) > 19.09 && F.work(sim) < 200_000 && closes(sim) && sim.net.unconverged == 0
end

@testset "air in a tank is a gas with heat capacity; valve travel follows simulated time" begin
    names, D = F.read_csv(joinpath(@__DIR__, "data", "lc26_day4_fill_1hz.csv"))
    scene() = (s = F.load_scene(joinpath(@__DIR__, "..", "scenes", "lc26_day4.toml")); s.schedule = F.log_schedule(D, F.guess_log_columns(names, s); t0 = 21950.0); s)
    sim = F.Sim(scene(); sample_dt = 0.01)
    tank = sim.net.nodes[findfirst(n -> n.name == "TANK", sim.net.nodes)]
    F.advance!(sim, 1.0)
    # What an air-filled tank sends down the vent line is close to N2O vapour at its own
    # temperature (not tank pressure over the trace N2O density, which would be about 12 MJ/kg).
    @test abs(tank.port_h[1] - F.vapour_PT(sim.net.env.P_atm, sim.u[tank.i+3]).h) < 2e4
    F.run_until!(sim, 20.0)                                       # the fill valve opens at 19 s
    after = [P for (t, P) in zip(sim.t_hist, sim.hist["TANK.P"]) if t > 19]
    # The first cold vapour has to cool the air too, so the tank does not drop below atmospheric.
    @test minimum(after) > sim.net.env.P_atm - 500 && closes(sim) && isempty(sim.failed)
    # The fill valve takes 0.5 s to open. A tenth of a second in, it is a fifth open, however the
    # run is stepped (and not 20 ms of travel further on for every call).
    opening(piece) = begin
        sim = F.Sim(scene()); F.run_until!(sim, 19.0)
        while F.simtime(sim) < 19.1 - 1e-9; F.advance!(sim, min(piece, 19.1 - F.simtime(sim))); end
        sim.net.edges[findfirst(e -> e.name == "AV201", sim.net.edges)].opening
    end
    @test 0.16 <= opening(0.1) <= 0.24 && 0.16 <= opening(0.001) <= 0.24
end

@testset "a supercritical bottle opened into an empty tank does not stall the junction solve" begin
    # A bottle left in the sun (62 °C, about 2100 psig) opened through an orifice, hose, valve and
    # hose. The imbalance here settles near 1e-6. If every such evaluation counts as failed and
    # pays for the whole fallback, these 5 s take 94 failures and 14 s of wall time. The solve
    # accepts it instead, and 2 failures remain.
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
    # HF4's fitted values with only what was measured changed (see the scene). The model gives
    # sealed-period tank pressure 2.5 psi RMS, fill line 3.2 psi, 4.59 kg loaded against 4.3 kg on
    # the load cell.
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
    # The model gives tank pressure 33 psi RMS through the vented fill and 4.9 psi through the
    # hold; its liquid is 0.7 K below the thermistors on average.
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


@testset "window icon files load at their sizes" begin
    icons = N2OFillSim.window_icons()
    @test size.(icons) == [(16, 16), (32, 32), (48, 48)]
    @test hasmethod(N2OFillSim.GLFW.SetWindowIcon, Tuple{N2OFillSim.GLFW.Window, typeof(icons)})
end

end
