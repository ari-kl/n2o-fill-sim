# The side panel: building the schematic in edit mode, driving and reading the run in run mode.

function param_inputs(values::Dict{String,Float64}, params, app::App)
    for p in params
        v = Ref{Cdouble}(values[p.key])
        ig.SetNextItemWidth(px(110))
        if ig.InputDouble(p.key, v, 0.0, 0.0, "%.4g") && isfinite(v[])
            snapshot!(app)
            values[p.key] = v[]
        end
        tip(p.help)
    end
end

function side_panel(app::App)
    if running(app)
        run_header(app)
        running(app) || return                                  # Stop was pressed
        if ig.BeginTabBar("runtabs")
            for (name, body) in (("Tanks", tanks_tab), ("Signals", signals_tab), ("Schedule", schedule_tab), ("View", view_tab))
                open = ig.BeginTabItem(name)
                mark("tab " * name)
                open && (body(app); ig.EndTabItem())
            end
            ig.EndTabBar()
        end
    else
        button("Simulate") && start_run!(app)
        tip("Check the schematic and start a run (paused). The plumbing cannot be changed while running.")
        if ig.BeginTabBar("edittabs")
            for (name, body) in (("Build", build_tab), ("Settings", a -> param_inputs(a.scene.settings, SETTINGS, a)), ("Schedule", schedule_tab), ("View", view_tab))
                open = ig.BeginTabItem(name)
                mark("tab " * name)
                open && (body(app); ig.EndTabItem())
            end
            ig.EndTabBar()
        end
    end
end

# ---------------------------------------------------------------- edit mode

function build_tab(app::App)
    s = app.scene
    ig.SeparatorText("File")
    wrapped((isempty(app.path) ? "not saved yet" : app.path) * (app.dirty ? "  (changed)" : ""))
    button("Open...") && request!(app, :open)
    tip("Ctrl+O")
    ig.SameLine(); button("Save") && save_scene!(app)
    tip("Ctrl+S")
    ig.SameLine(); button("Save as...") && browse!(app, :save)
    ig.SameLine(); button("New") && request!(app, :new)
    button("Undo") && undo!(app)
    tip("Ctrl+Z")
    ig.SameLine(); button("Redo") && redo!(app)
    tip("Ctrl+Y")
    ig.SameLine(); button("Load test log...") && browse!(app, :log)
    tip("Overlay a recorded test on the plots and, if you want, replay its valve sequence.")

    ig.SeparatorText("Place a part")
    for (i, kind) in enumerate(PALETTE)
        button(kind) && (app.placing = kind; app.wire_from = nothing)
        tip(KINDS[kind].help)
        i % 3 == 0 || i == length(PALETTE) || ig.SameLine()
    end
    wrapped("Click a pin, then another, to wire. Drag parts to move them; shift-click or drag a box to select several. R rotates, Delete removes, Ctrl+C / Ctrl+V / Ctrl+D copy, paste and duplicate. Click a wire to add a corner and drag it; right-click a corner to remove it, or a wire to cut it. Right-drag pans, wheel zooms, F fits.")
    length(app.sel) > 1 && text("$(length(app.sel)) parts selected; showing $(current(app))")

    i = isempty(app.sel) ? nothing : find_part(s, current(app))
    if i !== nothing
        part = s.parts[i]
        ig.SeparatorText(part.kind)
        wrapped(KINDS[part.kind].help)
        ig.SetNextItemWidth(px(160))
        if ig.InputText("name", app.name_buf, 64, ig.ImGuiInputTextFlags_EnterReturnsTrue)
            try
                snapshot!(app)
                rename_part!(s, part.name, str(app.name_buf))
                select!(app, str(app.name_buf))
            catch err
                app.msg = sprint(showerror, err)
            end
        end
        param_inputs(part.p, KINDS[part.kind].params, app)
    end
    errors, warnings = check(s)
    isempty(errors) && isempty(warnings) || ig.SeparatorText("Checks")
    for e in errors
        coloured_wrapped(WARN4, "error: " * e)
    end
    for w in warnings
        coloured_wrapped(AMBER4, "warning: " * w)
    end
end

# ---------------------------------------------------------------- run mode

"Action buttons for a part. `ctx` distinguishes the right-click menu's copies from the side panel's."
function part_actions(app::App, part::Part, ctx = "")
    sim = app.sim
    b(label) = button(label, ctx * label)
    if part.kind == "Valve"
        e = edge_of(app, part.name)
        b(valve_target(e) > 0.5 ? "Close" : "Open") && act!(sim, part.name, "toggle")
        ig.SameLine(); b("Stick open") && act!(sim, part.name, "stuck_open")
        ig.SameLine(); b("Stick closed") && act!(sim, part.name, "stuck_closed")
        e.stuck != 0 && (ig.SameLine(); b("Free") && act!(sim, part.name, "unstick"))
    elseif part.kind in VESSEL_KINDS && part.p["heater_W"] > 0
        on = node_of(app, part.name).heater_on
        b(on ? "Heater off" : "Heater on") && act!(sim, part.name, on ? "off" : "on")
    elseif part.kind in ("Regulator", "ReliefValve")
        ig.SetNextItemWidth(px(110))
        ig.InputDouble("set_psig", app.setpoint, 0.0, 0.0, "%.1f")
        ig.SameLine(); b("Apply") && act!(sim, part.name, "set", app.setpoint[])
    end
end

function run_header(app::App)
    sim, net = app.sim, app.sim.net
    button(app.paused ? "Run" : "Pause") && (app.paused = !app.paused)
    tip("Space")
    ig.SameLine(); button("Step 1 s") && advance!(sim, 1.0)
    ig.SameLine()
    if button("Stop / edit")
        app.sim, app.want_fit = nothing, true
        return
    end
    for (k, sp) in enumerate(SPEEDS)
        ig.RadioButton(speed_label(sp), app.speed == sp) && (app.speed = sp)
        tip("key $k")
        ig.SameLine()
    end
    ig.NewLine()
    text("t = $(round(simtime(sim), digits = 1)) s    running at ×$(round(app.rate, digits = 1))")
    tip("Frame $(round(1e3 * unsafe_load(ig.GetIO().DeltaTime), digits = 1)) ms: solver $(round(app.solver_ms, digits = 1)) ms, building the display $(round(app.frame_ms - app.solver_ms, digits = 1)) ms, the rest is the graphics driver and waiting for the screen.")
    isempty(sim.failed) || coloured_wrapped(WARN4, sim.failed)
    button("De-energise all valves") && fail_safe!(sim)
    tip("Loss of power: every valve goes to its fail state (normally-open ones open).")
    # Valves that are open only because nothing is holding them shut (spring-open valves before
    # they are first commanded, or after a de-energise).
    resting_open = [e.name for e in net.edges if e.kind === :valve && e.normally_open && !e.cmd && e.stuck == 0]
    isempty(resting_open) || coloured_wrapped(AMBER4, "Open because de-energised (spring-open): " * join(resting_open, ", ") * ". Click to command closed.")
end

function tanks_tab(app::App)
    sim, net, set = app.sim, app.sim.net, app.scene.settings
    lo, hi = set["target_lo_psig"], set["target_hi_psig"]
    gauge(P) = (P - p_atm(app)) / PSI                           # the window and MEOP are set in psig
    for (i, n) in enumerate(net.nodes)
        (n.kind === :vessel && app.scene.parts[find_part(app.scene, n.name)].kind == "Tank") || continue
        r = REPORTS[i]
        p = gauge(r.P)
        windowed = hi > lo                                      # a target window is optional
        inwin = windowed && lo <= p <= hi
        coloured(inwin ? OK4 : (0.85f0, 0.85f0, 0.85f0, 1f0), "$(n.name): $(fmt(app, :pressure, r.P))" *
            (windowed ? "  $(inwin ? "in" : p < lo ? "below" : "above") the $(round(Int, lo))-$(round(Int, hi)) psig window" : ""))
        ph = get(sim.hist, n.name * ".P", Float64[])            # time in window, counted once per sample
        seen, in_window = get(app.window_time, n.name, (1, 0.0))
        for k in seen+1:length(ph)
            lo <= gauge(ph[k]) <= hi && (in_window += sim.t_hist[k] - sim.t_hist[k-1])
        end
        app.window_time[n.name] = (max(length(ph), 1), in_window)
        text("  $(fmt(app, :mass, r.m))" * (windowed ? ", time in window $(round(Int, in_window)) s" : ""))
        text("  liquid: $(round(Int, 100 * r.liquid_frac)) % of volume, $(round(Int, 100 * (1 - r.x))) % of mass")
        text("  liquid $(fmt(app, :temperature, r.T))  ullage $(fmt(app, :temperature, r.Tv))  wall $(fmt(app, :temperature, r.Tw))")
        isnan(r.Psat) || text("  liquid's vapour pressure $(fmt(app, :pressure, r.Psat))")
        soaks = [Ts for Ts in (30.0, 35.0) if r.m > 0.05 && Ts + 273.15 > r.T]   # nothing to say about an empty tank, or one already hotter
        isempty(soaks) || text("  if closed and left to warm:")
        for Ts in soaks
            Ps, lf = warm_soak(r.m, n.V, Ts + 273.15)
            bad = gauge(Ps) > set["meop_psig"] || lf >= 1
            coloured(bad ? WARN4 : GREY4, "    at $(Int(Ts)) °C: $(fmt(app, :pressure, Ps))" * (lf >= 1 ? ", LIQUID-FULL" : ", liquid $(round(Int, 100lf)) % of volume"))
        end
        p > set["meop_psig"] && coloured(WARN4, "  ABOVE MEOP ($(round(Int, set["meop_psig"])) psig)")
        r.T > 307.0 && coloured(WARN4, "  near the critical temperature (36.4 °C)")
    end
    cm, ce = closure(sim)
    text("vented $(fmt(app, :mass, sim.u[end-2]))")
    tip("Balance check, should stay near zero: mass $(round(cm, sigdigits = 1)), energy $(round(ce, sigdigits = 1)) (relative).")
    net.unconverged > 0 && coloured(AMBER4, "junction solve failed to converge $(net.unconverged) times")

    ig.SeparatorText("Conditions")
    env = net.env
    amb, sun, wind = Ref(Cfloat(env.T_amb - 273.15)), Ref(Cfloat(env.solar)), Ref(Cfloat(env.wind))
    changed = ig.SliderFloat("ambient °C", amb, -20, 45, "%.0f")
    changed |= ig.SliderFloat("sun W/m²", sun, 0, 1100, "%.0f")
    changed |= ig.SliderFloat("wind m/s", wind, 0, 15, "%.1f")
    if changed
        env.T_amb, env.solar, env.wind = amb[] + 273.15, sun[], wind[]
        _restart!(sim)
    end

    i = isempty(app.sel) ? nothing : find_part(app.scene, current(app))
    if i !== nothing
        ig.SeparatorText(current(app))
        part_actions(app, app.scene.parts[i])
        wrapped(describe(app, app.scene.parts[i]))
    else
        wrapped("Click a valve to open or close it. Click any part to select it; right-click for more actions. Hover for its numbers.")
    end
    ig.SeparatorText("Record")
    if button("Export CSV")
        out = (isempty(app.path) ? "n2o-fill-sim" : splitext(app.path)[1]) * "_run.csv"
        write_csv(out, sim)
        write_events(events_path(out), sim)
        app.msg = "Wrote " * out * " (signals, SI units) and " * basename(events_path(out)) * " (valve and heater actions)"
    end
    ig.SameLine(); button("Load test log...") && browse!(app, :log)
    if button("Keep my clicks as the schedule")
        app.scene.schedule = copy(sim.log)
        app.dirty = true
        app.msg = "$(length(sim.log)) actions stored in the scene; Save to keep them. They replay on the next run."
    end
end

"Tick boxes for every recorded signal that can be plotted, grouped by part, with a filter."
function signals_tab(app::App)
    sim = app.sim
    ig.SetNextItemWidth(-1)
    ig.InputTextWithHint("##filter", "filter, e.g. TANK or mdot", app.filter_buf, length(app.filter_buf))
    needle = lowercase(str(app.filter_buf))
    ig.SmallButton("clear all") && empty!(app.probes)
    groups = Dict{String,Vector{String}}()
    for key in keys(sim.hist)
        signal_kind(key) === :other && continue
        isempty(needle) || occursin(needle, lowercase(key)) || continue
        push!(get!(groups, String(first(rsplit(key, '.'; limit = 2))), String[]), key)
    end
    ig.BeginChild("signals")
    for part in sort!(collect(keys(groups)))
        chosen = count(in(app.probes), groups[part])
        isempty(needle) || ig.SetNextItemOpen(true)
        if ig.TreeNode(part * (chosen > 0 ? "  ($chosen)" : "") * "###" * part)
            for key in sort!(groups[part])
                on = Ref(key in app.probes)
                ig.Checkbox(key, on) && (on[] ? push!(app.probes, key) : delete!(app.probes, key))
            end
            ig.TreePop()
        end
    end
    ig.EndChild()
end

function view_tab(app::App)
    ig.SeparatorText("Colour the plumbing by")
    for (mode, label) in ((:none, "nothing"), (:pressure, "pressure"), (:temperature, "temperature"), (:quality, "vapour share"))
        ig.RadioButton(label, app.colour_by === mode) && (app.colour_by = mode)
        mark("colour " * label)
        ig.SameLine()
    end
    ig.NewLine()
    wrapped("Shown while running. A wire is one colour (no pressure drop along a wire in this model); the step in colour across a part is the drop it causes. Faded means shut off on all sides.")
    if app.colour_by in (:pressure, :temperature)
        kind = colour_kind(app)
        auto = Ref(app.range_auto)
        ig.Checkbox("automatic range (excludes atmosphere)", auto) && (app.range_auto = auto[])
        if !app.range_auto
            lo, hi = Ref{Cdouble}(shown(app, kind, app.range[1])), Ref{Cdouble}(shown(app, kind, app.range[2]))
            ig.SetNextItemWidth(px(90)); a = ig.InputDouble("from", lo, 0.0, 0.0, "%.1f")
            ig.SameLine(); ig.SetNextItemWidth(px(90)); b = ig.InputDouble("to  [" * unit(app, kind) * "]", hi, 0.0, 0.0, "%.1f")
            (a || b) && hi[] > lo[] && (app.range = (to_si(app, kind, lo[]), to_si(app, kind, hi[])))
        end
        if running(app)
            button("Fix the range where it is now") && (app.range = app.shown_range; app.range_auto = false)
            i = isempty(app.sel) ? nothing : find_part(app.scene, current(app))
            if i !== nothing && app.scene.parts[i].kind in FLOW_KINDS
                e = edge_of(app, current(app))
                if button("Range: across " * current(app))
                    a, b = minmax(VALUES[e.a], VALUES[e.b])
                    pad = max(0.15 * (b - a), app.colour_by === :pressure ? 2e3 : 0.2)
                    app.range, app.range_auto = (a - pad, b + pad), false
                end
                tip("Spread the whole colour scale over the drop across the selected part.")
            end
        end
    end
    ig.SeparatorText("Units")
    for (field, label, choices) in UNIT_CHOICES
        text(rpad(label, 12))
        for c in choices
            ig.SameLine()
            ig.RadioButton(UNIT_LABEL[c] * "##" * label, getfield(app.units, field) === c) && setfield!(app.units, field, c)
            mark("unit " * UNIT_LABEL[c])
        end
    end
    wrapped("Units change what is displayed. Scene files keep the unit in each parameter's name, and exported CSV is always SI.")
    ig.SeparatorText("Plots")
    f = Ref(app.follow)
    ig.Checkbox("follow the run", f) && (app.follow = f[])
    mark("follow the run")
    tip("On: the time axis tracks the run. Off: drag and scroll in a plot to pan and zoom; all plots move together.")
    h = Ref(app.plot_frac)
    ig.SetNextItemWidth(px(160))
    ig.SliderFloat("plot height", h, 0.2f0, 0.7f0, "%.2f") && (app.plot_frac = h[]; app.want_fit = true)
end
