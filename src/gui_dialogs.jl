# Dialogs: choosing a file, mapping a test log onto the schematic, editing the schedule.

# ---------------------------------------------------------------- file browser

const BROWSE_TITLE = Dict(:open => "Open a schematic", :save => "Save the schematic as", :log => "Load a test log")

"Open the file browser for `what` (:open, :save or :log), starting beside the current file."
function browse!(app::App, what::Symbol)
    app.browse = what
    here = isempty(app.path) ? "" : dirname(abspath(app.path))
    isdir(app.browse_dir) || (app.browse_dir = pwd())
    isdir(here) && what !== :log && (app.browse_dir = here)
    setbuf!(app.browse_buf, what === :save ? basename(app.path) : "")
    nothing
end

# What happens when a file has been chosen.
function browse_done!(app::App, path)
    what = app.browse
    app.browse = :none
    if what === :open
        open_scene!(app, path)
    elseif what === :save
        save_scene!(app, endswith(lowercase(path), ".toml") ? path : path * ".toml")
    elseif what === :log
        load_log!(app, path)
    end
    nothing
end

function browser_popup(app::App)
    app.browse === :none && return
    title = BROWSE_TITLE[app.browse]
    ig.OpenPopup(title)
    ig.SetNextWindowSize((px(620), px(460)), ig.ImGuiCond_Appearing)
    ig.BeginPopupModal(title) || return
    ext = app.browse === :log ? ".csv" : ".toml"
    wrapped(app.browse_dir)
    button("Up") && (app.browse_dir = dirname(rstrip(app.browse_dir, ['/', '\\'])) |> d -> isempty(d) ? app.browse_dir : d)
    dirs, files = String[], String[]
    try
        for name in readdir(app.browse_dir)
            startswith(name, ".") && continue
            full = joinpath(app.browse_dir, name)
            isdir(full) ? push!(dirs, name) : endswith(lowercase(name), ext) && push!(files, name)
        end
    catch err
        coloured(WARN4, "cannot read this folder")
    end
    chosen = ""
    ig.BeginChild("entries", (0f0, -px(74)), ig.ImGuiChildFlags_Borders)
    for d in dirs
        ig.Selectable("[ " * d * " ]") && (app.browse_dir = joinpath(app.browse_dir, d))
    end
    for f in files
        ig.Selectable(f, f == str(app.browse_buf)) && setbuf!(app.browse_buf, f)
        mark("file " * f)
        ig.IsItemHovered() && ig.IsMouseDoubleClicked(0) && (chosen = f)
    end
    isempty(dirs) && isempty(files) && text("no folders or " * ext * " files here")
    ig.EndChild()
    ig.SetNextItemWidth(-px(60))
    ig.InputText("file", app.browse_buf, length(app.browse_buf), ig.ImGuiInputTextFlags_EnterReturnsTrue) && (chosen = str(app.browse_buf))
    tip("A name in this folder, or a full path. A folder path jumps there.")
    button(app.browse === :save ? "Save" : app.browse === :log ? "Load" : "Open", "browser ok") && (chosen = str(app.browse_buf))
    ig.SameLine()
    button("Cancel", "browser cancel") && (app.browse = :none)
    if !isempty(chosen)
        path = isabspath(chosen) ? chosen : joinpath(app.browse_dir, chosen)
        if isdir(path)
            app.browse_dir = path
            setbuf!(app.browse_buf, "")
        elseif app.browse === :save || isfile(path)
            browse_done!(app, path)
        else
            app.msg = "No such file: " * path
        end
    end
    app.browse === :none && ig.CloseCurrentPopup()
    ig.EndPopup()
end

# ---------------------------------------------------------------- loading a test log

"Read a log and open the mapping dialog on it."
function load_log!(app::App, path)
    try
        names, cols = read_csv(path; min_dt = 0.05)
        isempty(names) && error("the file has no header row")
        app.log_columns = guess_log_columns(names, app.scene)
        app.log_cols, app.log_path, app.log_open = cols, path, true
        app.log_t0[] = suggest_t0(cols, app.log_columns)
    catch err
        app.msg = "Could not read the log: " * first(sprint(showerror, err), 200)
    end
    nothing
end

"Ten seconds before the first mapped valve moves, or the start of the log if none does."
function suggest_t0(cols, columns::Vector{LogColumn})
    t = cols[columns[1].name]
    isempty(t) && return 0.0
    first_move = Inf
    for c in columns
        c.role === :valve || continue
        v = cols[c.name]
        i = findfirst(i -> !isnan(v[i]) && !isnan(v[i-1]) && v[i] != v[i-1], 2:length(v))
        i === nothing || (first_move = min(first_move, t[i+1]))
    end
    isfinite(first_move) ? max(t[1], first_move - 10.0) : t[1]
end

function apply_log!(app::App; replay)
    cols, columns = app.log_cols, app.log_columns
    t = cols[columns[1].name]
    t0 = app.log_t0[]
    set_overlay!(app, t .- t0, log_overlay(cols, columns, p_atm(app)))
    said = "$(length(app.overlay)) readings from $(basename(app.log_path)) are overlaid on the plots."
    if replay
        snapshot!(app)
        app.scene.schedule = log_schedule(cols, columns; t0)
        said *= " Its valve sequence ($(length(app.scene.schedule)) actions) is now the schedule" * (running(app) ? "; stop and Simulate again to replay it." : "; press Simulate to replay it.")
    end
    app.msg, app.log_open = said, false
    nothing
end

function log_popup(app::App)
    app.log_open || return
    ig.OpenPopup("Test log")
    ig.SetNextWindowSize((px(760), px(560)), ig.ImGuiCond_Appearing)
    ig.BeginPopupModal("Test log") || return
    cols, columns = app.log_cols, app.log_columns
    t = cols[columns[1].name]
    wrapped("$(app.log_path):  $(length(t)) rows kept (20 per second at most), $(isempty(t) ? 0 : round(t[end] - t[1], digits = 1)) s, time column `$(columns[1].name)`")
    ig.SetNextItemWidth(px(120))
    ig.InputDouble("log time that becomes t = 0", app.log_t0, 0.0, 0.0, "%.2f")
    isempty(t) || (ig.SameLine(); button("start of log") && (app.log_t0[] = t[1]))
    ig.SameLine(); button("10 s before the first valve moves") && (app.log_t0[] = suggest_t0(cols, columns))
    unsure = false
    flags = ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_BordersInnerV | ig.ImGuiTableFlags_ScrollY
    if ig.BeginTable("columns", 4, flags, (0f0, -px(96)))
        ig.TableSetupColumn("use"); ig.TableSetupColumn("column"); ig.TableSetupColumn("taken as"); ig.TableSetupColumn("a 1 in this column means")
        ig.TableHeadersRow()
        for (k, c) in enumerate(columns)
            c.role in (:time, :other) && continue
            ig.PushID(k)
            ig.TableNextRow()
            ig.TableNextColumn()
            on = Ref(c.use)
            ig.Checkbox("##use", on) && (c.use = on[])
            ig.TableNextColumn(); text(c.name)
            ig.TableNextColumn()
            text(c.role === :pressure ? "pressure, psi gauge: overlay" : c.role === :mass ? "mass, kg: overlay" :
                 isempty(c.part) ? "no matching part" : (c.role === :valve ? "valve " : "heater on ") * c.part)
            ig.TableNextColumn()
            if c.role === :valve && !isempty(c.part)
                ig.RadioButton("open", c.one_means === :open) && (c.one_means = :open)
                mark("open " * c.name)
                ig.SameLine()
                ig.RadioButton("closed", c.one_means === :closed) && (c.one_means = :closed)
                if c.one_means === :unknown
                    ig.SameLine(); coloured(AMBER4, "choose")
                    unsure = true
                end
            elseif c.role === :heater
                text("on")
            end
            ig.PopID()
        end
        ig.EndTable()
    end
    unsure && coloured_wrapped(AMBER4, "Columns marked \"choose\" will not be replayed until you say what a 1 means. Logs with valve_ columns use 1 = open; later recordings with solenoid_ columns write 1 for a reported state of CLOSED.")
    button("Overlay the readings") && apply_log!(app; replay = false)
    ig.SameLine()
    button("Overlay, and replay the valves") && apply_log!(app; replay = true)
    tip("Replaces the scene's schedule with this log's valve and heater sequence (undo brings the old one back).")
    ig.SameLine()
    button("Cancel", "log cancel") && (app.log_open = false)
    app.log_open || ig.CloseCurrentPopup()
    ig.EndPopup()
end

# ---------------------------------------------------------------- schedule editor

const VALVE_ACTIONS = ["open", "close", "stuck_open", "stuck_closed", "unstick"]

# Parts a schedule can act on, and what each can be told to do.
function schedulable(scene::Scene)
    out = Tuple{String,Vector{String}}[]
    for p in scene.parts
        p.kind == "Valve" && push!(out, (p.name, VALVE_ACTIONS))
        p.kind in VESSEL_KINDS && p.p["heater_W"] > 0 && push!(out, (p.name, ["on", "off"]))
        p.kind in ("Regulator", "ReliefValve") && push!(out, (p.name, ["set"]))
    end
    out
end

"The schedule as a table: edit, add, delete, sort and shift scripted actions. Read-only while running."
function schedule_tab(app::App)
    sched = app.scene.schedule
    if running(app)
        sim = app.sim
        text("$(length(sim.schedule)) scripted actions; next is marked. Stop to edit.")
        ig.BeginChild("sched")
        for (i, a) in enumerate(sim.schedule)
            line = "$(lpad(round(a.t, digits = 1), 8)) s   $(rpad(a.part, 14)) $(a.action)" * (a.action == "set" ? " $(a.value) psig" : "")
            i == sim.next_action ? coloured(OK4, "> " * line) : i < sim.next_action ? coloured(GREY4, "  " * line) : text("  " * line)
        end
        ig.EndChild()
        return
    end
    parts = schedulable(app.scene)
    names = first.(parts)
    wrapped("Actions applied at set times during a run. Times are seconds from the start.")
    if isempty(parts)
        text("Nothing to schedule: add a valve, a heater or a regulator first.")
        return
    end
    if button("Add action")
        snapshot!(app)
        push!(sched, Action(isempty(sched) ? 0.0 : sched[end].t + 10, names[1], parts[1][2][1], 0.0))
    end
    ig.SameLine(); button("Sort by time") && (snapshot!(app); sort!(sched; by = a -> a.t))
    ig.SameLine(); button("Clear") && !isempty(sched) && (snapshot!(app); empty!(sched))
    ig.SetNextItemWidth(px(90))
    ig.InputDouble("s", app.shift, 0.0, 0.0, "%.1f")
    ig.SameLine()
    if button("Shift all") && app.shift[] != 0
        snapshot!(app)
        app.scene.schedule = [Action(max(a.t + app.shift[], 0.0), a.part, a.action, a.value) for a in sched]
    end
    sched = app.scene.schedule
    remove = 0
    if ig.BeginTable("sched", 5, ig.ImGuiTableFlags_RowBg | ig.ImGuiTableFlags_ScrollY | ig.ImGuiTableFlags_SizingFixedFit)
        ig.TableSetupColumn("t [s]"); ig.TableSetupColumn("part"); ig.TableSetupColumn("action"); ig.TableSetupColumn("value"); ig.TableSetupColumn("")
        ig.TableHeadersRow()
        for (i, a) in enumerate(sched)
            ig.PushID(i)
            ig.TableNextRow()
            t, part, action, value = a.t, a.part, a.action, a.value
            ig.TableNextColumn()
            tr = Ref{Cdouble}(t)
            ig.SetNextItemWidth(px(74))
            ig.InputDouble("##t", tr, 0.0, 0.0, "%.1f") && (t = max(tr[], 0.0))
            ig.TableNextColumn()
            k = something(findfirst(==(part), names), 0)
            pidx = Ref{Int32}(k - 1)
            ig.SetNextItemWidth(px(118))
            if ig.Combo("##part", pidx, k == 0 ? vcat(names, [part * " (missing)"]) : names) && pidx[] + 1 <= length(names)
                part = names[pidx[]+1]
                k = pidx[] + 1
            end
            ig.TableNextColumn()
            if k != 0
                actions = parts[k][2]
                action in actions || (action = actions[1])
                ai = Ref{Int32}(findfirst(==(action), actions) - 1)
                ig.SetNextItemWidth(px(104))
                ig.Combo("##action", ai, actions) && (action = actions[ai[]+1])
            else
                text(action)
            end
            ig.TableNextColumn()
            if action == "set"
                vr = Ref{Cdouble}(value)
                ig.SetNextItemWidth(px(64))
                ig.InputDouble("##v", vr, 0.0, 0.0, "%.0f") && (value = vr[])
                tip("set-point, psig")
            end
            ig.TableNextColumn()
            ig.SmallButton("x") && (remove = i)
            mark("remove action $i")
            if (t, part, action, value) != (a.t, a.part, a.action, a.value)
                snapshot!(app)
                sched[i] = Action(t, part, action, value)
            end
            ig.PopID()
        end
        ig.EndTable()
    end
    remove != 0 && (snapshot!(app); deleteat!(app.scene.schedule, remove))
    issorted(app.scene.schedule; by = a -> a.t) || coloured_wrapped(AMBER4, "Not in time order: actions are applied in time order whatever the order here. Sort to see them that way.")
end
