# The window: a schematic you edit (place, wire, set parameters) and then run (click valves, watch
# pressures, probe signals onto plots). Dear ImGui redraws everything every frame, so there is no
# widget state to keep in step with the simulation: each frame steps the sim, then draws what is.
#
#   gui.jl         application state, units, the frame, keyboard shortcuts, the window loop
#   gui_canvas.jl  the schematic: drawing, colouring by pressure, mouse handling
#   gui_panels.jl  the side panel: file, palette, properties, run controls, view settings
#   gui_plots.jl   the plots
#   gui_dialogs.jl the file browser, loading a test log, the schedule editor
#
# All scene changes go through scene.jl and all physics through network.jl; these files only draw
# and route input.

import CImGui as ig
import ImPlot
import GLFW
import ModernGL as GL

"Monitor scale factor (1.5 on a 150 % display). Every pixel size goes through `px`."
const UI = Ref(1f0)
px(x) = Float32(x) * UI[]
seconds() = time_ns() / 1e9                    # Base.time() only ticks every ~12 ms on Windows

# ---------------------------------------------------------------- colours

const RGB = NTuple{3,Float32}
u32(c::RGB, a = 255) = ig.IM_COL32(round(Int, c[1]), round(Int, c[2]), round(Int, c[3]), a)
mix(a::RGB, b::RGB, f) = (a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, a[3] + (b[3] - a[3]) * f)
col(r, g, b, a = 255) = ig.IM_COL32(r, g, b, a)
const C_LINE, C_DIM, C_SEL, C_WIRE = col(220, 220, 225), col(120, 120, 130), col(255, 200, 60), col(150, 190, 255)
const C_OPEN, C_SHUT, C_LIQ, C_VAP, C_WARN = col(70, 200, 110), col(205, 80, 80), col(60, 130, 230), col(70, 80, 100), col(255, 90, 60)
const WARN4, OK4, GREY4, AMBER4 = (1f0, 0.4f0, 0.3f0, 1f0), (0.4f0, 0.9f0, 0.5f0, 1f0), (0.7f0, 0.7f0, 0.7f0, 1f0), (1f0, 0.8f0, 0.3f0, 1f0)

# Colour scale for "colour by": blue (low) through green and yellow to red (high).
const SCALE = RGB[(58, 76, 192), (65, 140, 240), (50, 200, 220), (70, 225, 140), (190, 235, 70), (250, 195, 50), (245, 120, 40), (200, 45, 30)]
function scale_colour(f)
    f = clamp(Float32(f), 0f0, 1f0) * (length(SCALE) - 1)
    i = min(floor(Int, f), length(SCALE) - 2)
    mix(SCALE[i+1], SCALE[i+2], f - i)
end

# ---------------------------------------------------------------- units

"Display units. Scene files keep their unit-in-the-name keys and the CSV export stays SI."
mutable struct Units
    P::Symbol                                  # :psig, :psia, :bar, :kPa
    T::Symbol                                  # :C, :K
    m::Symbol                                  # :kg, :lb
    flow::Symbol                               # :gs, :kgs
end
const UNIT_CHOICES = ((:P, "pressure", (:psig, :psia, :bar, :kPa)), (:T, "temperature", (:C, :K)), (:m, "mass", (:kg, :lb)), (:flow, "flow", (:gs, :kgs)))
const UNIT_LABEL = Dict(:psig => "psig", :psia => "psia", :bar => "bar", :kPa => "kPa", :C => "°C", :K => "K", :kg => "kg", :lb => "lb", :gs => "g/s", :kgs => "kg/s")

# A valve or heater trace as the list of its changes: what gets drawn does not grow with the run.
mutable struct Steps
    t::Vector{Float64}
    v::Vector{Float64}
    seen::Int                 # samples of the history already looked at
end

Base.@kwdef mutable struct App
    scene::Scene
    path::String = ""
    sim::Union{Nothing,Sim} = nothing
    paused::Bool = true
    speed::Float64 = 1.0
    rate::Float64 = 0.0                       # simulated seconds per wall second actually achieved
    solver_ms::Float64 = 0.0                  # wall time the last frame spent stepping the simulation
    frame_ms::Float64 = 0.0                   # wall time the last frame spent in our code (solver + building the UI)
    sel::Vector{String} = String[]            # selected part names; the last one is the "current" part
    placing::String = ""
    wire_from::Union{Nothing,Tuple{String,String}} = nothing
    dragging::Bool = false
    drag_last::Tuple{Int,Int} = (0, 0)        # grid cell the pointer was in at the last drag step
    drag_via::Tuple{Int,Int} = (0, 0)         # (wire, corner) being dragged, or (0, 0)
    box_from::Union{Nothing,Tuple{Float64,Float64}} = nothing   # where a selection box was started
    clip::Tuple{Vector{Part},Vector{Wire}} = (Part[], Wire[])   # copied parts and the wires among them
    origin::Tuple{Float32,Float32} = (60f0, 190f0)
    zoom::Float32 = 1f0
    canvas::NTuple{4,Float32} = (0f0, 0f0, 800f0, 500f0)   # where the schematic was drawn last frame: x, y, width, height
    want_fit::Bool = true                     # zoom to fit on the next frame
    undo::Vector{Scene} = Scene[]
    redo::Vector{Scene} = Scene[]
    dirty::Bool = false                       # unsaved changes
    confirm::Symbol = :none                   # what to do once unsaved changes are dealt with (:new, :open, :quit)
    quit::Bool = false
    probes::Set{String} = Set{String}()
    cache::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()   # history in display units, for plotting
    cache_units::NTuple{4,Symbol} = (:psig, :C, :kg, :gs)                  # the units `cache` was built with
    steps::Dict{String,Steps} = Dict{String,Steps}()                       # valve and heater traces as lists of changes
    window_time::Dict{String,Tuple{Int,Float64}} = Dict{String,Tuple{Int,Float64}}()   # per tank: samples counted, seconds in the window
    msg::String = ""
    overlay_t::Vector{Float64} = Float64[]    # a test log shown next to the run: its time base,
    overlay::Vector{Tuple{String,Symbol,Vector{Float64}}} = Tuple{String,Symbol,Vector{Float64}}[]   # and (label, kind, SI values) per reading
    name_buf::Vector{UInt8} = zeros(UInt8, 128)
    filter_buf::Vector{UInt8} = zeros(UInt8, 128)
    setpoint::Base.RefValue{Cdouble} = Ref(700.0)
    units::Units = Units(:psig, :C, :kg, :gs)
    colour_by::Symbol = :none                 # :none, :pressure, :temperature, :quality
    range_auto::Bool = true
    range::Tuple{Float64,Float64} = (1e5, 7e6)         # colour scale limits in SI, used when not automatic
    shown_range::Tuple{Float64,Float64} = (1e5, 7e6)   # the limits actually in use this frame
    follow::Bool = true                       # plots track the run; off lets you zoom and pan
    xlo::Base.RefValue{Cdouble} = Ref(0.0)    # shared time axis of all plots
    xhi::Base.RefValue{Cdouble} = Ref(60.0)
    cursor_t::Float64 = NaN                   # time under the mouse in a plot (NaN if none)
    plot_frac::Float32 = 0.42f0               # share of the window height given to the plots in run mode
    # file browser
    browse::Symbol = :none                    # :none, :open, :save or :log
    browse_dir::String = pwd()
    browse_buf::Vector{UInt8} = zeros(UInt8, 512)
    # a test log being mapped onto the scene (see logs.jl)
    log_open::Bool = false
    log_path::String = ""
    log_cols::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()
    log_columns::Vector{LogColumn} = LogColumn[]
    log_t0::Base.RefValue{Cdouble} = Ref(0.0)
    shift::Base.RefValue{Cdouble} = Ref(0.0)  # schedule editor: seconds to shift every action by
end

App(scene::Scene, path = "") = App(; scene, path)

running(app::App) = app.sim !== nothing
set_overlay!(app::App, t, readings) = (app.overlay_t = collect(t); app.overlay = readings; nothing)
current(app::App) = isempty(app.sel) ? "" : app.sel[end]
select!(app::App, name) = (empty!(app.sel); isempty(name) || (push!(app.sel, name); setbuf!(app.name_buf, name)); nothing)
str(buf) = String(buf[1:something(findfirst(==(0x00), buf), length(buf) + 1)-1])
setbuf!(buf, s) = (fill!(buf, 0x00); copyto!(buf, codeunits(first(s, length(buf) - 1))); buf)

"Call before any change to the scene: makes it undoable and marks the file as modified."
function snapshot!(app::App)
    push!(app.undo, deepcopy(app.scene))
    length(app.undo) > 100 && popfirst!(app.undo)
    empty!(app.redo)
    app.dirty = true
    nothing
end
undo!(app::App) = isempty(app.undo) || (push!(app.redo, app.scene); app.scene = pop!(app.undo); empty!(app.sel); app.dirty = true)
redo!(app::App) = isempty(app.redo) || (push!(app.undo, app.scene); app.scene = pop!(app.redo); empty!(app.sel); app.dirty = true)

# ImGui's Text takes a printf format; these take plain strings.
text(s) = ig.TextUnformatted(s)
coloured(c, s) = ig.TextColored(c, replace(s, "%" => "%%"))
coloured_wrapped(c, s) = (ig.PushStyleColor(ig.ImGuiCol_Text, ig.ImVec4(c...)); wrapped(s); ig.PopStyleColor())
wrapped(s) = ig.TextWrapped(replace(s, "%" => "%%"))
tip(s) = ig.IsItemHovered() && ig.SetTooltip(replace(s, "%" => "%%"))

"Where each named widget was drawn last frame (x0, y0, x1, y1), so a test can click it wherever it landed."
const RECTS = Dict{String,NTuple{4,Float32}}()
mark(label) = (a = ig.GetItemRectMin(); b = ig.GetItemRectMax(); RECTS[label] = (a.x, a.y, b.x, b.y); nothing)
button(label, id = label) = (pressed = ig.Button(label); mark(id); pressed)

# --- values in display units. `kind` is :pressure, :temperature, :mass, :flow, :share, :state or :quality.
p_atm(app::App) = app.scene.settings["P_atm_kPa"] * 1e3
function shown(app::App, kind::Symbol, v)
    u = app.units
    kind === :pressure ? (u.P === :psig ? (v - p_atm(app)) / PSI : u.P === :psia ? v / PSI : u.P === :bar ? v / 1e5 : v / 1e3) :
    kind === :temperature ? (u.T === :C ? v - 273.15 : v) :
    kind === :mass ? (u.m === :kg ? v : v * 2.2046226) :
    kind === :flow ? (u.flow === :gs ? 1e3v : v) : v
end
"Inverse of `shown`: a number typed in display units back to SI."
function to_si(app::App, kind::Symbol, v)
    u = app.units
    kind === :pressure ? (u.P === :psig ? v * PSI + p_atm(app) : u.P === :psia ? v * PSI : u.P === :bar ? v * 1e5 : v * 1e3) :
    kind === :temperature ? (u.T === :C ? v + 273.15 : v) :
    kind === :mass ? (u.m === :kg ? v : v / 2.2046226) :
    kind === :flow ? (u.flow === :gs ? v / 1e3 : v) : v
end
unit(app::App, kind::Symbol) = kind === :pressure ? UNIT_LABEL[app.units.P] : kind === :temperature ? UNIT_LABEL[app.units.T] :
                               kind === :mass ? UNIT_LABEL[app.units.m] : kind === :flow ? UNIT_LABEL[app.units.flow] : kind === :share ? "%" : ""
function fmt(app::App, kind::Symbol, v; with_unit = true)
    x = shown(app, kind, v)
    s = !isfinite(x) ? "-" :
        kind === :pressure ? (app.units.P === :bar ? string(round(x, digits = 2)) : string(round(Int, x))) :
        kind === :temperature ? string(round(x, digits = 1)) :
        kind === :mass ? string(round(x, digits = 2)) :
        kind === :flow ? string(round(x, sigdigits = 3)) : string(round(x, digits = 2))
    with_unit && !isempty(unit(app, kind)) ? s * " " * unit(app, kind) : s
end

"What kind of quantity a recorded signal is, from its name."
function signal_kind(key)
    endswith(key, ".P") || endswith(key, ".Psat") || endswith(key, ".set") ? :pressure :
    endswith(key, ".T") || endswith(key, ".Tv") || endswith(key, ".Tw") ? :temperature :
    endswith(key, ".mdot") ? :flow :
    endswith(key, ".m") ? :mass :
    endswith(key, "_pct") ? :share :
    endswith(key, ".position") || endswith(key, ".open_cmd") || endswith(key, ".heater_on") ? :state : :other
end

include("gui_canvas.jl")
include("gui_panels.jl")
include("gui_plots.jl")
include("gui_dialogs.jl")

# ---------------------------------------------------------------- file actions

function start_run!(app::App)
    errors, _ = check(app.scene)
    if !isempty(errors)
        app.msg = "Cannot run: " * join(errors, "; ")
        return
    end
    try
        app.sim = Sim(app.scene)
        app.paused, app.msg, app.follow, app.want_fit = true, "", true, true     # the plots take space: re-fit
        empty!(app.cache); empty!(app.steps); empty!(app.window_time)
        app.wire_from, app.placing = nothing, ""
        if isempty(app.probes)                              # sensible default traces
            for n in app.sim.net.nodes
                n.kind === :vessel || continue
                push!(app.probes, n.name * ".P", n.name * ".T", n.name * ".m")
                app.scene.parts[find_part(app.scene, n.name)].kind == "Tank" && push!(app.probes, n.name * ".liquid_vol_pct", n.name * ".liquid_mass_pct")
            end
            for g in keys(app.sim.net.gauges)
                push!(app.probes, g * ".P")
            end
            for e in app.sim.net.edges
                e.kind === :valve && push!(app.probes, e.name * ".mdot", e.name * ".position")
            end
        end
    catch err
        app.msg = "Cannot run: " * first(sprint(showerror, err), 300)
    end
    nothing
end

function open_scene!(app::App, path)
    try
        app.scene = load_scene(path)
        app.path, app.msg, app.dirty, app.want_fit = path, "Opened " * path, false, true
        empty!(app.sel); empty!(app.undo); empty!(app.redo); empty!(app.probes)
    catch err
        app.msg = "Open failed: " * first(sprint(showerror, err), 200)
    end
    nothing
end

"Save to `path` (default: where the scene came from). With nowhere to save to, ask."
function save_scene!(app::App, path = app.path)
    isempty(path) && (browse!(app, :save); return false)
    try
        save_scene(path, app.scene)
        app.path, app.msg, app.dirty = path, "Saved " * path, false
        return true
    catch err
        app.msg = "Save failed: " * first(sprint(showerror, err), 200)
        return false
    end
end

function new_scene!(app::App)
    app.scene, app.path, app.dirty = Scene(), "", false
    empty!(app.sel); empty!(app.undo); empty!(app.redo); empty!(app.probes)
    nothing
end

"Do `what` (:new, :open, :quit) now, or ask first if there are unsaved changes."
request!(app::App, what::Symbol) = (app.dirty ? (app.confirm = what) : perform!(app, what); nothing)
function perform!(app::App, what::Symbol)
    app.confirm = :none
    what === :new && new_scene!(app)
    what === :open && browse!(app, :open)             # choosing the file comes after the changes are dealt with
    what === :quit && (app.quit = true)
    nothing
end

function unsaved_popup(app::App)
    app.confirm === :none && return
    ig.OpenPopup("Unsaved changes")
    if ig.BeginPopupModal("Unsaved changes", C_NULL, ig.ImGuiWindowFlags_AlwaysAutoResize)
        text("The schematic has changes that are not saved.")
        what = app.confirm
        if button("Save")
            save_scene!(app) ? perform!(app, what) : (app.confirm = :none)     # nowhere to save yet: the browser takes over
        end
        ig.SameLine()
        button("Discard") && perform!(app, what)
        ig.SameLine()
        button("Cancel") && (app.confirm = :none)
        app.confirm === :none && ig.CloseCurrentPopup()
        ig.EndPopup()
    end
end

# ---------------------------------------------------------------- the frame

const SPEEDS = (1.0, 5.0, 20.0, 100.0, 1e6)        # the last is "as fast as the solver goes", for long holds
speed_label(sp) = sp >= 1e6 ? "max" : "×$(Int(sp))"

"Keyboard shortcuts. They work wherever the mouse is, unless a text field has the keyboard."
function shortcuts(app::App)
    io = ig.GetIO()
    unsafe_load(io.WantTextInput) && return
    ctrl = unsafe_load(io.KeyCtrl)
    pressed(k) = ig.IsKeyPressed(k, false)
    !ctrl && pressed(ig.ImGuiKey_F) && (app.want_fit = true)
    if ctrl && pressed(ig.ImGuiKey_S)
        save_scene!(app)
    elseif ctrl && pressed(ig.ImGuiKey_O)
        request!(app, :open)
    end
    if running(app)
        pressed(ig.ImGuiKey_Space) && (app.paused = !app.paused)
        for (k, key) in enumerate((ig.ImGuiKey_1, ig.ImGuiKey_2, ig.ImGuiKey_3, ig.ImGuiKey_4, ig.ImGuiKey_5))
            pressed(key) && (app.speed = SPEEDS[k])
        end
        pressed(ig.ImGuiKey_Escape) && empty!(app.sel)
        return
    end
    s = app.scene
    if pressed(ig.ImGuiKey_Escape)
        app.wire_from, app.placing = nothing, ""
    elseif !ctrl && pressed(ig.ImGuiKey_R) && !isempty(app.sel)
        snapshot!(app)
        for name in app.sel
            s.parts[find_part(s, name)].rot = mod(s.parts[find_part(s, name)].rot + 1, 4)
        end
    elseif pressed(ig.ImGuiKey_Delete) && !isempty(app.sel)
        snapshot!(app)
        foreach(name -> delete_part!(s, name), app.sel)
        empty!(app.sel)
    elseif ctrl && pressed(ig.ImGuiKey_Z)
        undo!(app)
    elseif ctrl && pressed(ig.ImGuiKey_Y)
        redo!(app)
    elseif ctrl && pressed(ig.ImGuiKey_A)
        app.sel = [p.name for p in s.parts]
    elseif ctrl && pressed(ig.ImGuiKey_C) && !isempty(app.sel)
        app.clip = copy_parts(s, app.sel)
        app.msg = "Copied $(length(app.clip[1])) parts. Ctrl+V pastes them at the pointer."
    elseif ctrl && pressed(ig.ImGuiKey_V) && !isempty(app.clip[1])
        # Paste so that the copied group's top-left part lands under the pointer (or a little
        # down and right of the original, if the pointer is not over the schematic).
        m = ig.GetMousePos()
        x, y, w, h = app.canvas
        x0, y0 = minimum(p.x for p in app.clip[1]), minimum(p.y for p in app.clip[1])
        dx, dy = (x <= m.x <= x + w && y <= m.y <= y + h) ? round.(Int, s2w(app, m.x, m.y)) .- (x0, y0) : (2, 2)
        snapshot!(app)
        app.sel = paste_parts!(s, app.clip, dx, dy)
    elseif ctrl && pressed(ig.ImGuiKey_D) && !isempty(app.sel)
        snapshot!(app)
        app.sel = paste_parts!(s, copy_parts(s, app.sel), 2, 2)
    end
    nothing
end

"One frame: advance the simulation, then draw the panes."
function frame(app::App)
    io = ig.GetIO()
    size = unsafe_load(io.DisplaySize)
    dt = clamp(unsafe_load(io.DeltaTime), 0.001f0, 0.1f0)
    app.solver_ms = 0.0
    if running(app) && !app.paused && isempty(app.sim.failed)
        # Give the solver at most 12 ms of each frame. If that is not enough to keep up (the first
        # seconds of a fill, or a high speed setting) simulated time falls behind rather than the
        # window stuttering; the panel shows the speed actually achieved.
        start, want, done = seconds(), app.speed * dt, 0.0
        while done < want - 1e-9 && seconds() - start < 0.012 && isempty(app.sim.failed)
            done += advance!(app.sim, want - done; max_steps = 1)
        end
        app.solver_ms = 1e3 * (seconds() - start)
        app.rate = 0.9 * app.rate + 0.1 * done / dt
    end
    shortcuts(app)
    side = min(px(400), 0.45f0 * size.x)
    bottom = running(app) ? app.plot_frac * size.y : 0f0
    fixed = ig.ImGuiWindowFlags_NoMove | ig.ImGuiWindowFlags_NoResize | ig.ImGuiWindowFlags_NoCollapse
    ig.SetNextWindowPos((0f0, 0f0)); ig.SetNextWindowSize((size.x - side, size.y - bottom))
    ig.Begin("Schematic" * (running(app) ? "  (running: click valves)" : "  (editing)") * "###schematic", C_NULL, fixed | ig.ImGuiWindowFlags_NoScrollbar | ig.ImGuiWindowFlags_NoScrollWithMouse)
    draw_canvas(app)
    if ig.BeginPopup("partmenu")
        i = find_part(app.scene, current(app))
        i === nothing || (text(current(app)); part_actions(app, app.scene.parts[i], "menu "))
        ig.EndPopup()
    end
    ig.End()
    ig.SetNextWindowPos((size.x - side, 0f0)); ig.SetNextWindowSize((side, size.y))
    ig.Begin(running(app) ? "Run###side" : "Edit###side", C_NULL, fixed)
    side_panel(app)
    isempty(app.msg) || (ig.Separator(); wrapped(app.msg))
    ig.End()
    if running(app)
        ig.SetNextWindowPos((0f0, size.y - bottom)); ig.SetNextWindowSize((size.x - side, bottom))
        ig.Begin("Plots", C_NULL, fixed | ig.ImGuiWindowFlags_NoScrollbar)
        plot_panel(app)
        ig.End()
    end
    unsaved_popup(app)
    browser_popup(app)
    log_popup(app)
    nothing
end

# ---------------------------------------------------------------- the window

# Read the frame just drawn (before it is swapped to the screen) and save it as a binary PPM.
function screenshot(window, path)
    w, h = GLFW.GetFramebufferSize(window)
    buf = Vector{UInt8}(undef, 3w * h)
    GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
    GL.glReadBuffer(GL.GL_BACK)
    GL.glReadPixels(0, 0, w, h, GL.GL_RGB, GL.GL_UNSIGNED_BYTE, buf)
    open(path, "w") do f
        print(f, "P6\n$w $h\n255\n")
        for r in h-1:-1:0
            write(f, @view buf[3w*r+1:3w*(r+1)])
        end
    end
end

window_title(app::App) = "n2o-fill-sim - " * (isempty(app.path) ? "untitled" : basename(app.path)) * (app.dirty ? " *" : "")

"""
    gui(app; shot = "", frames = 0, script = nothing)

Open the window and run until it is closed. With `shot`, save a picture after `frames` frames and
exit; `script(app, frame_number)` can drive the app meanwhile (used for checks).

The loop is ours rather than the library's so that it can sleep when nothing is happening (a
paused schematic should not keep a laptop core busy), ask before closing with unsaved changes,
size itself to the display, and take a screenshot of exactly what was drawn.
"""
function gui(app::App; shot = "", frames = 0, script = nothing)
    GLFW.WindowHint(GLFW.VISIBLE, true)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
    GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 3)
    if Sys.isapple()
        GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE)
        GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, true)
    end
    monitor = GLFW.GetPrimaryMonitor()
    UI[] = clamp(GLFW.GetMonitorContentScale(monitor).xscale, 1f0, 3f0)
    mode = GLFW.GetVideoMode(monitor)
    w = min(round(Int, 1500 * UI[]), round(Int, 0.94 * mode.width))
    h = min(round(Int, 950 * UI[]), round(Int, 0.88 * mode.height))
    title = window_title(app)
    window = GLFW.CreateWindow(w, h, title)
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(1)

    ctx = ig.CreateContext()
    pctx = ImPlot.CreateContext()
    io = ig.GetIO()
    io.IniFilename = C_NULL                       # fixed layout; do not litter imgui.ini files
    ig.StyleColorsDark()
    style = ig.GetStyle()
    ig.ScaleAllSizes(style, UI[])
    style.FontScaleDpi = UI[]
    font = joinpath(pkgdir(ig), "fonts", "Roboto-Medium.ttf")
    isfile(font) && ig.AddFontFromFileTTF(unsafe_load(io.Fonts), font, 15f0)
    ig.lib.ImGui_ImplGlfw_InitForOpenGL(Ptr{ig.lib.GLFWwindow}(window.handle), true)
    ig.lib.ImGui_ImplOpenGL3_Init("#version 330")
    app.origin = (px(60), px(190))
    count, awake, was_running = 0, 10, false
    try
        while !app.quit
            # Sleep until something happens unless the simulation is running or there was input a
            # moment ago (ImGui needs a few frames to settle after an event).
            busy = (running(app) && !app.paused) || awake > 0 || !isempty(shot) || script !== nothing
            busy ? GLFW.PollEvents() : GLFW.WaitEvents(0.25)
            if GLFW.WindowShouldClose(window)
                GLFW.SetWindowShouldClose(window, false)
                isempty(shot) ? request!(app, :quit) : (app.quit = true)
            end
            ig.lib.ImGui_ImplOpenGL3_NewFrame()
            ig.lib.ImGui_ImplGlfw_NewFrame()
            count += 1
            # A test script injects its input here: after the platform layer has queued the real
            # mouse position for this frame, so that the scripted position is the one that counts.
            script === nothing || script(app, count)
            ig.NewFrame()
            started = seconds()
            frame(app)
            app.frame_ms = 1e3 * (seconds() - started)
            # Garbage: drawing makes about 20 kB of it a frame (text, mostly). Left alone it is
            # collected every few hundred frames in a single pause of 25-80 ms, which shows as a
            # dropped frame or five. A full collection where a pause cannot be seen (start-up,
            # entering or leaving run mode) and a small one every two seconds (about half a
            # millisecond) keep it from ever building up.
            mode = running(app)
            (count == 3 || mode != was_running) ? GC.gc() : count % 120 == 0 && GC.gc(false)
            was_running = mode
            delta = unsafe_load(io.MouseDelta)
            active = delta.x != 0 || delta.y != 0 || unsafe_load(io.MouseWheel) != 0 || ig.IsMouseDown(0) || ig.IsMouseDown(1) ||
                     ig.IsMouseDown(2) || unsafe_load(io.WantTextInput) || ig.IsAnyItemActive()
            awake = active ? 10 : max(awake - 1, 0)
            ig.Render()
            fw, fh = GLFW.GetFramebufferSize(window)
            GL.glViewport(0, 0, fw, fh)
            GL.glClearColor(0.08, 0.08, 0.1, 1)
            GL.glClear(GL.GL_COLOR_BUFFER_BIT)
            ig.lib.ImGui_ImplOpenGL3_RenderDrawData(Ptr{Cint}(ig.GetDrawData()))
            if !isempty(shot) && count >= frames
                screenshot(window, shot)
                app.quit = true
            end
            GLFW.SwapBuffers(window)
            window_title(app) == title || (title = window_title(app); GLFW.SetWindowTitle(window, title))
        end
    finally
        ig.lib.ImGui_ImplOpenGL3_Shutdown()
        ig.lib.ImGui_ImplGlfw_Shutdown()
        ImPlot.DestroyContext(pctx)
        ig.DestroyContext(ctx)
        GLFW.DestroyWindow(window)
        app.quit = false
    end
    nothing
end
