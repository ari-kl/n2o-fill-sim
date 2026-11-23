# Drives the real window with injected mouse and key events and checks what happened.
#   julia --project test/gui_input.jl [screenshot.ppm]
# Needs a display. Keep the mouse off the window while it runs (about 5 s): real mouse movement
# is input too and can make a step miss. A failed step is worth re-running before believing it.
using N2OFillSim
const F = N2OFillSim
import CImGui as ig

const ROOT = dirname(@__DIR__)
app = F.App(F.load_scene(joinpath(ROOT, "scenes", "hf4.toml")), "")
steps = Tuple{Int,Function}[]
at(n, f) = push!(steps, (n, f))
io() = ig.GetIO()
# The platform layer reports the real mouse position every frame while the window has focus, so the
# scripted position is re-sent every frame (see `script` below) or the pointer would snap back.
pointer = Ref{Any}(nothing)                                   # (world x, world y, dx, dy) the pointer is held at
move(wx, wy; dx = 0, dy = 0) = (pointer[] = (wx, wy, dx, dy))
function hold()
    p = pointer[]
    p === nothing && return
    if p isa NTuple{2,Float32}                                # a place on the screen, whatever the view does
        ig.AddMousePosEvent(io(), p...)
    elseif p isa String                                       # a widget, by the label it was drawn with
        haskey(F.RECTS, p) || return
        x0, y0, x1, y1 = F.RECTS[p]
        ig.AddMousePosEvent(io(), (x0 + x1) / 2, (y0 + y1) / 2)
    else
        s = F.w2s(app, p[1], p[2])
        ig.AddMousePosEvent(io(), s[1] + p[3], s[2] + p[4])
    end
end
"Click a button, tab or radio button by its label, wherever it was drawn."
function press(n, label)
    at(n, () -> (pointer[] = label)); at(n + 2, () -> ig.AddMouseButtonEvent(io(), 0, true)); at(n + 3, () -> ig.AddMouseButtonEvent(io(), 0, false))
end
function click(n, wx, wy; button = 0, dx = 0, dy = 0)
    at(n, () -> move(wx, wy; dx, dy)); at(n + 1, () -> ig.AddMouseButtonEvent(io(), button, true)); at(n + 2, () -> ig.AddMouseButtonEvent(io(), button, false))
end
key(n, k) = (at(n, () -> ig.AddKeyEvent(io(), k, true)); at(n + 1, () -> ig.AddKeyEvent(io(), k, false)))
ctrl(n, k) = (at(n, () -> ig.AddKeyEvent(io(), ig.ImGuiMod_Ctrl, true)); key(n + 1, k); at(n + 4, () -> ig.AddKeyEvent(io(), ig.ImGuiMod_Ctrl, false)))
failures = String[]
function check(n, what, f)
    at(n, () -> begin
        ok = try f() catch err; false end
        println(ok ? "  ok    " : "  FAIL  ", what)
        if !ok
            push!(failures, what)
            m = ig.GetMousePos()                              # where ImGui thinks the pointer is, for diagnosis
            println("          frame ", n, ": pointer held at ", pointer[], ", ImGui has it at (", round(m.x), ", ", round(m.y), "), canvas ", round.(app.canvas), ", selected ", app.sel)
        end
    end)
end
has(name) = F.find_part(app.scene, name) !== nothing
nparts() = length(app.scene.parts)
nwires() = length(app.scene.wires)
edge(name) = app.sim.net.edges[findfirst(e -> e.name == name, app.sim.net.edges)]
n0, w0 = nparts(), nwires()

# --- edit mode
click(10, 23, -5)                                             # select the PT202 gauge
check(14, "click selects a part", () -> F.current(app) == "PT202")
key(15, ig.ImGuiKey_Delete)
check(18, "Delete removes it and its wire", () -> !has("PT202") && nwires() == w0 - 1 && app.dirty)
ctrl(20, ig.ImGuiKey_Z)
check(26, "Ctrl+Z brings it back", () -> has("PT202") && nwires() == w0)
ctrl(28, ig.ImGuiKey_Y)
check(34, "Ctrl+Y removes it again", () -> !has("PT202"))
ctrl(36, ig.ImGuiKey_Z)
click(44, 13.5, 4; button = 1)                                # right-click the AV201 - DUMP_TEE wire
check(48, "right-click cuts a wire", () -> nwires() == w0 - 1)
ctrl(50, ig.ImGuiKey_Z)
at(58, () -> (app.placing = "Orifice"))
click(59, 30, 6)
check(63, "placing a part from the palette", () -> nparts() == n0 + 1 && app.scene.parts[end].kind == "Orifice" && (app.scene.parts[end].x, app.scene.parts[end].y) == (30, 6))
click(65, 29, 6); click(69, 25, 3)                            # its pin a → TANK.bottom
check(73, "clicking two pins wires them", () -> nwires() == w0 + 1)
at(75, () -> move(30, 6)); at(76, () -> ig.AddMouseButtonEvent(io(), 0, true)); at(78, () -> move(30, 7.4)); at(80, () -> move(30, 8)); at(82, () -> ig.AddMouseButtonEvent(io(), 0, false))
check(84, "dragging moves a part", () -> (app.scene.parts[end].x, app.scene.parts[end].y) == (30, 8))
key(86, ig.ImGuiKey_R)
check(89, "R rotates it", () -> app.scene.parts[end].rot == 1)
key(90, ig.ImGuiKey_Delete)
check(93, "and Delete removes it", () -> nparts() == n0 && nwires() == w0)
key(94, ig.ImGuiKey_F)

# --- run mode
at(100, () -> F.start_run!(app))
check(104, "the run starts paused", () -> app.sim !== nothing && app.paused && isempty(app.msg))
key(105, ig.ImGuiKey_Space)
check(108, "Space runs", () -> !app.paused)
key(109, ig.ImGuiKey_3)
check(112, "3 sets x20", () -> app.speed == 20.0)
click(114, 12, 4); click(120, 21, 4)                          # AV201, AV203
check(126, "left-click commands a valve open", () -> F.valve_target(edge("AV201")) == 1 && F.valve_target(edge("AV203")) == 1 && F.current(app) == "AV203")
at(128, () -> (app.colour_by = :pressure))
check(200, "nitrous is flowing and the tank is filling", () -> edge("AV203").mdot > 0.02 && app.sim.hist["TANK.P"][end] > 3e5)
check(201, "valve position was recorded", () -> app.sim.hist["AV203.position"][end] == 1 && app.sim.hist["AV203.position"][1] == 0)
click(202, 29, -2; button = 1)                                # right-click the vent: its menu opens
press(208, "menu Stick open")
check(214, "the right-click menu sticks a valve", () -> edge("AV204").stuck == 1)
click(216, 5, 9)                                              # click away to close the menu
click(222, 12, 4)
check(230, "clicking again closes the valve", () -> F.valve_target(edge("AV201")) == 0 && isempty(app.sim.failed) && length(app.sim.log) == 4)
press(232, "tab View"); press(238, "colour temperature")
check(244, "the View tab switches what the plumbing is coloured by", () -> app.colour_by === :temperature)
press(246, "tab Tanks"); press(252, "Pause")
check(258, "the Pause button pauses", () -> app.paused)

# --- back to editing: unsaved changes are protected
press(260, "Stop / edit")
check(266, "Stop returns to editing", () -> app.sim === nothing)
at(268, () -> F.request!(app, :new))
check(272, "New with unsaved changes asks first", () -> app.confirm === :new && nparts() == n0)
press(273, "Cancel")
check(279, "Cancel keeps the schematic", () -> app.confirm === :none && nparts() == n0)
at(280, () -> F.request!(app, :new))
press(284, "Discard")
check(290, "Discard starts a new one", () -> nparts() == 0 && !app.dirty)

# --- file browser, test log, schedule editor
at(292, () -> (app.browse_dir = joinpath(ROOT, "scenes"); F.browse!(app, :open)))
press(296, "file regulated_fill.toml"); press(302, "browser ok")
check(308, "the browser opens a schematic", () -> endswith(app.path, "regulated_fill.toml") && has("REG") && app.browse === :none)
at(310, () -> F.load_log!(app, joinpath(ROOT, "test", "data", "hf4_fill_10hz.csv")))
check(314, "a test log opens the mapping dialog", () -> app.log_open && any(c -> c.part == "AV203" && c.one_means === :open, app.log_columns) && 9018 < app.log_t0[] < 9019)
press(316, "Overlay, and replay the valves")
check(322, "its readings are overlaid and its valve sequence becomes the schedule", () -> !app.log_open && length(app.overlay) == 3 && length(app.scene.schedule) > 10 && app.dirty)
press(324, "tab Schedule")
at(330, () -> (global n_actions = length(app.scene.schedule)))
press(331, "Add action")
check(337, "the schedule editor adds an action", () -> length(app.scene.schedule) == n_actions + 1)
press(338, "remove action 1")
check(344, "and removes one", () -> length(app.scene.schedule) == n_actions)

# --- several parts at once, copies, wire corners (still on regulated_fill: REG at (7,4), AV201 at (12,4), AV203 at (21,4))
part(name) = app.scene.parts[F.find_part(app.scene, name)]
shift(n, down) = at(n, () -> ig.AddKeyEvent(io(), ig.ImGuiMod_Shift, down))
press(346, "tab Build")
click(352, 12, 4)
shift(356, true); click(357, 21, 4); shift(361, false)
check(363, "shift-click selects a second part", () -> Set(app.sel) == Set(["AV201", "AV203"]))
at(365, () -> move(21, 4)); at(366, () -> ig.AddMouseButtonEvent(io(), 0, true)); at(368, () -> move(21, 4.6)); at(370, () -> move(21, 6)); at(372, () -> ig.AddMouseButtonEvent(io(), 0, false))
check(374, "dragging one moves the whole selection", () -> (part("AV201").y, part("AV203").y) == (6, 6) && part("AV201").x == 12)
ctrl(376, ig.ImGuiKey_Z)
check(382, "one undo puts both back", () -> (part("AV201").y, part("AV203").y) == (4, 4))
at(384, () -> (app.sel = ["AV201", "DUMP_TEE", "AV202"]; global n_before = (nparts(), nwires())))
ctrl(385, ig.ImGuiKey_D)
check(391, "Ctrl+D duplicates the selection with its wires", () -> (nparts(), nwires()) == n_before .+ (3, 2) && length(app.sel) == 3 && !("AV201" in app.sel))
ctrl(392, ig.ImGuiKey_Z)
at(398, () -> (app.sel = ["PT202"]))
ctrl(399, ig.ImGuiKey_C)
at(404, () -> move(30, 9)); ctrl(405, ig.ImGuiKey_V)
check(411, "Ctrl+C then Ctrl+V pastes at the pointer", () -> nparts() == n_before[1] + 1 && (app.scene.parts[end].kind, app.scene.parts[end].x, app.scene.parts[end].y) == ("Gauge", 30, 9))
ctrl(412, ig.ImGuiKey_Z)
at(418, () -> (empty!(app.sel); move(18, 1))); at(419, () -> ig.AddMouseButtonEvent(io(), 0, true)); at(421, () -> move(20, 3)); at(423, () -> move(23.5, 5.5)); at(425, () -> ig.AddMouseButtonEvent(io(), 0, false))
check(427, "dragging a box over empty space selects what is inside", () -> Set(app.sel) == Set(["QD", "AV203"]))
wire_ab(a, b) = app.scene.wires[findfirst(w -> w.a == a && w.b == b, app.scene.wires)]
at(429, () -> move(24, 4)); at(430, () -> ig.AddMouseButtonEvent(io(), 0, true)); at(432, () -> move(24, 5)); at(434, () -> move(24, 7)); at(436, () -> ig.AddMouseButtonEvent(io(), 0, false))
check(438, "clicking a wire adds a corner that can be dragged", () -> wire_ab(("AV203", "b"), ("TANK", "bottom")).via == [(24, 7)])
click(440, 24, 7; button = 1)
check(444, "right-clicking the corner removes it", () -> isempty(wire_ab(("AV203", "b"), ("TANK", "bottom")).via))

# The pointer is sent before this frame's steps as well as after. With the window focused and the
# real mouse parked outside it, the platform layer queues the real position every frame; a button
# event queued straight after that would be taken at the real position, and the click would miss.
script(app, n) = (hold(); foreach(((k, f),) -> k == n && f(), steps); hold())
# --- units, plot follow, quitting with unsaved changes
press(446, "tab View"); press(452, "unit bar")
check(458, "the View tab switches units", () -> app.units.P === :bar && endswith(F.fmt(app, :pressure, 3e5), " bar"))
press(460, "follow the run")
check(466, "'follow the run' can be turned off", () -> !app.follow)
at(468, () -> F.request!(app, :quit))
check(472, "quitting with unsaved changes asks first", () -> app.confirm === :quit && !app.quit)
press(473, "Cancel")
check(479, "and Cancel keeps the window open", () -> app.confirm === :none && !app.quit)

at(480, () -> (app.sel = ["PT201", "PT202"]; global n_before_delete = nparts()))
key(481, ig.ImGuiKey_Delete)
check(486, "Delete removes every selected part", () -> nparts() == n_before_delete - 2 && !has("PT201") && !has("PT202"))

# --- a second run (regulated_fill, less its two gauges): what a click shows while paused
at(490, () -> F.start_run!(app))
click(494, 12, 4)                                             # AV201, with time standing still
check(498, "a valve clicked while paused is marked as told to open before it has moved", () -> app.paused && F.simtime(app.sim) == 0 &&
      edge("AV201").opening == 0 && F.pending(edge("AV201")) == 1 && F.pending(edge("AV203")) == 0)
at(500, () -> (app.speed = 5.0; global t_before_step = F.simtime(app.sim)))
press(501, "Step 1 s")                                        # the button acts on frame 504
# A frame is at most 0.1 s long to the simulation, so at x5 the first one does half the step at most.
check(506, "Step 1 s is worked off over several frames", () -> !app.paused && 0 < app.step_left < 1 && t_before_step < F.simtime(app.sim) < t_before_step + 1)
check(700, "and ends paused, one second on", () -> app.paused && app.step_left == 0 && abs(F.simtime(app.sim) - t_before_step - 1) < 1e-6)

# --- back in the editor: the right button pans when dragged and acts only when it is not
press(702, "Stop / edit")
at(708, () -> (move(13.5, 4); global view_before = (app.origin, nwires(), length(app.undo), F.wire_at(app, 13.5, 4))))    # on the AV201 - DUMP_TEE wire
at(710, () -> ig.AddMouseButtonEvent(io(), 1, true))
# The view moves under a dragging pointer, so the drag is given in screen positions.
at(711, () -> (global drag_from = F.w2s(app, 13.5, 4); pointer[] = drag_from .+ (25f0, 0f0)))
at(713, () -> (pointer[] = drag_from .+ (60f0, 10f0)))
at(715, () -> ig.AddMouseButtonEvent(io(), 1, false))
check(718, "a right-drag that starts on a wire pans and leaves the wire", () -> view_before[4] !== nothing && nwires() == view_before[2] && length(app.undo) == view_before[3] && app.origin[1] > view_before[1][1] + 30)
at(719, () -> (app.origin = view_before[1]))

# --- a number typed into a parameter box is one undo step, and undo keeps the part selected
const PARAM = F.KINDS["Valve"].params[1].key
type(n, s) = at(n, () -> ig.AddInputCharactersUTF8(io(), s))
press(720, "tab Build")
click(726, 12, 4)                                             # AV201
at(730, () -> (global edit_before = (length(app.undo), part("AV201").p[PARAM])))
press(731, "param " * PARAM)
type(736, "3"); type(738, "1"); type(740, "7")                # three changes of value: 3, 31, 317
key(742, ig.ImGuiKey_Enter)
check(746, "a typed parameter is one undo step", () -> part("AV201").p[PARAM] == 317 && length(app.undo) == edit_before[1] + 1)
ctrl(748, ig.ImGuiKey_Z)
check(754, "undoing it keeps the part selected", () -> part("AV201").p[PARAM] == edit_before[2] && app.sel == ["AV201"])
shift(756, true); ctrl(757, ig.ImGuiKey_Z); shift(762, false)
check(764, "Ctrl+Shift+Z redoes", () -> part("AV201").p[PARAM] == 317)
ctrl(766, ig.ImGuiKey_Z)

shot = isempty(ARGS) ? joinpath(tempdir(), "n2o-fill-sim_gui_input.ppm") : ARGS[1]
F.gui(app; shot, frames = parse(Int, get(ENV, "N2O_FILL_SIM_TEST_FRAMES", "774")), script)
println(isempty(failures) ? "all window checks passed" : "FAILED: " * join(failures, "; "))
println("screenshot: ", shot)
exit(isempty(failures) ? 0 : 1)
