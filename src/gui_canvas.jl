# The schematic canvas: where things are, how they are drawn, and what the mouse does.

grid(app::App) = px(26) * app.zoom              # pixels per grid unit
w2s(app::App, x, y) = (app.origin[1] + Float32(x) * grid(app), app.origin[2] + Float32(y) * grid(app))
s2w(app::App, sx, sy) = ((sx - app.origin[1]) / grid(app), (sy - app.origin[2]) / grid(app))

# Part-local offset → screen, applying the part's rotation.
function local_pt(app::App, part::Part, dx, dy)
    for _ in 1:mod(part.rot, 4)
        dx, dy = -dy, dx
    end
    w2s(app, part.x + dx, part.y + dy)
end

half_extent(kind) = kind == "Tank" ? (0.9, 2.6) : kind == "Bottle" ? (0.8, 1.9) : kind == "Junction" ? (0.3, 0.3) : (0.9, 0.6)

function part_at(app::App, wx, wy)
    for part in Iterators.reverse(app.scene.parts)
        hx, hy = half_extent(part.kind)
        isodd(part.rot) && ((hx, hy) = (hy, hx))
        cy = part.kind == "Bottle" ? part.y + 0.2 : part.y
        abs(wx - part.x) <= hx && abs(wy - cy) <= hy && return part
    end
    nothing
end

"""
How near the pointer must be, in grid units, to take a pin or a wire corner: 0.35 of a grid square,
but never less than 8 pixels on the screen. Zoomed right out, 0.35 of a square is 3 pixels, too
small a thing to hit.
"""
grab_radius(app::App) = max(0.35, px(8) / grid(app))

# The nearest within reach, not the first: zoomed out, the reach of neighbouring pins overlaps.
# Zoomed out, the reach of a part's own pins also covers nearly all of its body, which would leave
# nowhere to click to select or drag it. So over a body, beyond the 0.35 a pin always has, a pin
# is taken only where the pointer is nearer to it than to the middle of the part.
function pin_at(app::App, wx, wy)
    best, found = grab_radius(app), nothing
    over = part_at(app, wx, wy)
    over === nothing || (best = min(best, max(0.35, hypot(wx - over.x, wy - over.y))))
    for part in app.scene.parts, pin in KINDS[part.kind].pins
        px_, py_ = pin_pos(part, pin)
        d = hypot(wx - px_, wy - py_)
        d < best && ((best, found) = (d, (part.name, pin.name)))
    end
    found
end

"""
Call `f(x1, y1, x2, y2, leg)` for each straight run of a wire, in grid units: from its first pin
through any corners placed on it to its second pin, horizontal first and then vertical between
consecutive stops. `leg` counts the stretches between corners. Stops and returns true as soon as
`f` returns true. Nothing is allocated: this runs for every wire on every frame.
"""
function each_run(f, app::App, w::Wire)
    pa, pb = find_part(app.scene, w.a[1]), find_part(app.scene, w.b[1])
    (pa === nothing || pb === nothing) && return false
    A, B = app.scene.parts[pa], app.scene.parts[pb]
    ia = findfirst(p -> p.name == w.a[2], KINDS[A.kind].pins)
    ib = findfirst(p -> p.name == w.b[2], KINDS[B.kind].pins)
    (ia === nothing || ib === nothing) && return false
    prev, last = pin_pos(A, KINDS[A.kind].pins[ia]), pin_pos(B, KINDS[B.kind].pins[ib])
    for leg in 1:length(w.via)+1
        (x1, y1), (x2, y2) = prev, leg <= length(w.via) ? w.via[leg] : last
        if x1 != x2 && y1 != y2
            (f(x1, y1, x2, y1, leg) === true || f(x2, y1, x2, y2, leg) === true) && return true
        else
            f(x1, y1, x2, y2, leg) === true && return true
        end
        prev = (x2, y2)
    end
    false
end

"The wire under a point, as (index of the wire, which leg between its pins and corners), or `nothing`."
function wire_at(app::App, wx, wy)
    near(ax, ay, bx, by) = (ax == bx ? abs(wx - ax) < 0.25 && min(ay, by) - 0.1 <= wy <= max(ay, by) + 0.1 :
                                       abs(wy - ay) < 0.25 && min(ax, bx) - 0.1 <= wx <= max(ax, bx) + 0.1)
    hit = Ref(0)
    for (i, w) in enumerate(app.scene.wires)
        each_run((x1, y1, x2, y2, leg) -> near(x1, y1, x2, y2) && (hit[] = leg; true), app, w) && return (i, hit[])
    end
    nothing
end

"The wire corner under a point (the nearest within reach), as (wire index, corner index), or `nothing`."
function via_at(app::App, wx, wy)
    best, found = grab_radius(app), nothing
    for (i, w) in enumerate(app.scene.wires), (k, v) in enumerate(w.via)
        d = hypot(wx - v[1], wy - v[2])
        d < best && ((best, found) = (d, (i, k)))
    end
    found
end

"Choose zoom and offset so the whole schematic is visible, with room for the labels."
function zoom_fit!(app::App)
    app.want_fit = false
    isempty(app.scene.parts) && return
    x0 = y0 = Inf
    x1 = y1 = -Inf
    for p in app.scene.parts
        hx, hy = half_extent(p.kind)
        x0, x1 = min(x0, p.x - hx - 1.5), max(x1, p.x + hx + (p.kind in VESSEL_KINDS ? 6.5 : 3.0))
        y0, y1 = min(y0, p.y - hy - 3.2), max(y1, p.y + hy + (p.kind in VESSEL_KINDS ? 3.0 : 1.8))   # room for the colour bar above, read-outs below
    end
    _, _, w, h = app.canvas
    app.zoom = clamp(min(w / ((x1 - x0) * px(26)), h / ((y1 - y0) * px(26))), 0.35f0, 1.6f0)
    cx, cy = app.canvas[1] + w / 2, app.canvas[2] + h / 2
    app.origin = (cx - Float32((x0 + x1) / 2) * grid(app), cy - Float32((y0 + y1) / 2) * grid(app))
    nothing
end

# ---------------------------------------------------------------- live values and "colour by"

edge_of(app::App, name) = (k = findfirst(e -> e.name == name, app.sim.net.edges); k === nothing ? nothing : app.sim.net.edges[k])
node_index(app::App, name) = something(findfirst(n -> n.name == name, app.sim.net.nodes), 0)
node_of(app::App, name) = (k = node_index(app, name); k == 0 ? nothing : app.sim.net.nodes[k])

# Per-frame scratch, filled by `survey!`: one report per node, its colour-by value, whether it is
# cut off from everything by closed valves.
const REPORTS = Any[]
const VALUES = Float64[]
const CUT_OFF = Bool[]
const OFF_SCALE = Bool[]      # at atmospheric pressure: drawn grey, not as the bottom of the scale

colour_kind(app::App) = app.colour_by === :pressure ? :pressure : app.colour_by === :temperature ? :temperature : :quality

"Evaluate every node once for this frame and work out the colour scale limits."
function survey!(app::App)
    net, u = app.sim.net, app.sim.u
    empty!(REPORTS); empty!(VALUES); empty!(CUT_OFF); empty!(OFF_SCALE)
    lo, hi = Inf, -Inf
    for (i, n) in enumerate(net.nodes)
        r = node_report(net, n, u)
        v = app.colour_by === :temperature ? r.T : app.colour_by === :quality ? r.x : r.P
        cut = n.kind === :junction && all(k -> flow_area(net.edges[k], net.nodes[net.edges[k].a].P, net.nodes[net.edges[k].b].P) <= 0, net.adj[i])
        push!(REPORTS, r); push!(VALUES, v); push!(CUT_OFF, cut)
        # The atmosphere, and anything vented down to it, is left out of the automatic range:
        # with them in, a fill's 5 psi drops would all be one colour against 0-950 psig.
        vented = n.kind === :atm || (app.colour_by === :pressure && r.P < net.env.P_atm + 2e4)
        push!(OFF_SCALE, vented && app.range_auto)
        cut || vented || (lo = min(lo, v); hi = max(hi, v))
    end
    if app.colour_by === :quality
        lo, hi = 0.0, 1.0
    elseif !app.range_auto
        lo, hi = app.range
    end
    isfinite(lo) || ((lo, hi) = (0.0, 1.0))
    pad = app.colour_by === :pressure ? 2e3 : app.colour_by === :temperature ? 0.2 : 0.0
    hi - lo < 2pad && ((lo, hi) = ((lo + hi) / 2 - pad, (lo + hi) / 2 + pad))
    app.shown_range = (lo, hi)
    nothing
end

"Colour of node `i` on the current scale; `nothing` when not colouring."
function tint(app::App, i)
    (app.colour_by === :none || i == 0) && return nothing
    OFF_SCALE[i] && return C_DIM
    lo, hi = app.shown_range
    u32(scale_colour((VALUES[i] - lo) / (hi - lo)), CUT_OFF[i] ? 90 : 255)
end
pin_tint(app::App, part, pin) = tint(app, get(app.sim.net.pin_node, (part, pin), 0))

function colour_bar(dl, app::App)
    x, y, w, h = app.canvas
    bw, bh, n = px(220), px(10), 44
    x0, y0 = x + px(14), y + px(26)
    for k in 0:n-1
        ig.AddRectFilled(dl, (x0 + bw * k / n, y0), (x0 + bw * (k + 1) / n + 1f0, y0 + bh), u32(scale_colour((k + 0.5) / n)))
    end
    lo, hi = app.shown_range
    kind = colour_kind(app)
    label(v) = kind === :quality ? string(round(Int, 100v), " %") : fmt(app, kind, v; with_unit = false)
    title = kind === :quality ? "vapour share of mass" : string(app.colour_by, " [", unit(app, kind), "]", app.range_auto ? "  (auto range)" : "")
    ig.AddText(dl, (x0, y0 - px(17)), C_LINE, title)
    ig.AddText(dl, (x0, y0 + bh + px(2)), C_DIM, label(lo))
    mid = label((lo + hi) / 2)
    ig.AddText(dl, (x0 + bw / 2 - ig.CalcTextSize(mid).x / 2, y0 + bh + px(2)), C_DIM, mid)
    ig.AddText(dl, (x0 + bw - ig.CalcTextSize(label(hi)).x, y0 + bh + px(2)), C_DIM, label(hi))
end

# ---------------------------------------------------------------- drawing

"Which way a valve has been told to go and has not yet got to: +1 opening, -1 closing, 0 if it is where it was told to be."
pending(e) = Int(sign(valve_target(e) - e.opening))

function draw_part(dl, app::App, part::Part)
    L(dx, dy) = local_pt(app, part, dx, dy)
    z = UI[] * app.zoom
    selected = part.name in app.sel
    c = selected ? C_SEL : C_LINE
    line(a, b, color = c, t = 2f0) = ig.AddLine(dl, L(a...), L(b...), color, t * UI[])
    k = part.kind
    live = running(app)
    if k in VESSEL_KINDS
        hx, hy = k == "Tank" ? (0.8, 2.5) : (0.7, 1.7)
        cy = k == "Tank" ? 0.0 : 0.3
        frac = 0.0
        outline = c
        if live
            i = node_index(app, part.name)
            frac = clamp(REPORTS[i].liquid_frac, 0.0, 1.0)
            selected || (outline = something(tint(app, i), c))
        elseif part.p["fill_kg"] > 0
            T = part.p["T_C"] + 273.15
            ρ = part.p["fill_kg"] / (part.p["volume_L"] * 1e-3)
            frac = T < SAT_TMAX ? clamp((ρ - sat_T(T).ρv) / (sat_T(T).ρl - sat_T(T).ρv), 0.0, 1.0) : 1.0
        end
        # Vessels are drawn upright whatever the rotation, so the liquid level reads as a level.
        x0, y0 = w2s(app, part.x - hx, part.y + cy - hy)
        x1, y1 = w2s(app, part.x + hx, part.y + cy + hy)
        ig.AddRectFilled(dl, (x0, y0), (x1, y1), C_VAP, 8f0 * z)
        frac > 0 && ig.AddRectFilled(dl, (x0, y1 - (y1 - y0) * Float32(frac)), (x1, y1), C_LIQ, 8f0 * z)
        ig.AddRect(dl, (x0, y0), (x1, y1), outline, 8f0 * z, (app.colour_by === :none ? 2f0 : 3f0) * UI[])
        for pin in KINDS[k].pins                               # stubs from the body to each pin
            sx, sy = pin_pos(part, pin)
            ig.AddLine(dl, w2s(app, clamp(sx, part.x - hx, part.x + hx), clamp(sy, part.y + cy - hy, part.y + cy + hy)), w2s(app, sx, sy), outline, 2f0 * UI[])
        end
        if k == "Tank"                                         # dip tube mark
            yv = part.y + cy + hy - 2hy * part.p["vent_height_frac"]
            ig.AddLine(dl, w2s(app, part.x + hx - 0.35, yv), w2s(app, part.x + hx, yv), C_DIM, 2f0 * UI[])
            ig.AddLine(dl, w2s(app, part.x + hx - 0.35, yv), w2s(app, part.x + hx - 0.35, part.y - 2), C_DIM, 2f0 * UI[])
        end
    elseif k in FLOW_KINDS
        # When colouring, each side of a flow part takes the colour of what it is connected to, so
        # the step in colour across the part is the pressure drop it causes.
        ca = live ? pin_tint(app, part.name, "a") : nothing
        cb = live ? pin_tint(app, part.name, "b") : nothing
        tinted = ca !== nothing && cb !== nothing
        la, lb = something(ca, c), something(cb, c)
        thick = tinted ? 3f0 : 2f0
        # Valves show their position; check valves, relief valves and regulators whether anything
        # is passing. `shut` is how far closed the part is drawn, 0 to 1 (nothing has a position
        # while editing: drawn filled, in grey).
        state, shut = C_DIM, 1f0
        if live && k == "Valve"
            e = edge_of(app, part.name)
            shut = e.opening > 0.99 ? 0f0 : e.opening < 0.01 ? 1f0 : 1f0 - Float32(e.opening)
            state = shut == 0 ? C_OPEN : shut == 1 ? C_SHUT : col(230, 180, 60)
            e.stuck != 0 && ig.AddText(dl, L(-0.5, 0.75), C_WARN, "STUCK")
        elseif live && k != "Pipe" && k != "Orifice"
            shut = abs(edge_of(app, part.name).mdot) > 1e-6 ? 0f0 : 1f0
            state = shut == 0 ? C_OPEN : C_SHUT
        end
        fa, fb = tinted ? (la, lb) : (state, state)
        if k == "Pipe"
            line((-1, 0), (-0.7, 0), la, thick); line((0.7, 0), (1, 0), lb, thick)
            pts = ((-0.7, 0.0), (-0.45, -0.22), (-0.15, 0.22), (0.15, -0.22), (0.45, 0.22), (0.7, 0.0))
            for i in 1:5                                       # the drop is spread along a pipe
                seg = tinted ? u32(mix(scale_rgb(app, part.name, "a"), scale_rgb(app, part.name, "b"), (i - 0.5f0) / 5)) : c
                line(pts[i], pts[i+1], seg, 3f0)
            end
        elseif k == "Orifice"
            line((-1, 0), (0, 0), la, thick); line((0, 0), (1, 0), lb, thick)
            line((0, -0.5), (0, -0.12), c, 3f0); line((0, 0.12), (0, 0.5), c, 3f0)
        else
            line((-1, 0), (-0.6, 0), la, thick); line((0.6, 0), (1, 0), lb, thick)
            # Open and shut differ by shape as well as colour, so the state reads without telling
            # red from green: shut is filled with a bar across the flow path, open is an outline
            # only, part-way is an outline with the fill closing in on the throat. The outline thins
            # when zoomed out so that the smallest open symbol still has a hole in it.
            thin = clamp(2f0 * z, 1f0, thick * UI[])
            if k == "CheckValve"
                ig.AddTriangle(dl, L(-0.6, -0.45), L(-0.6, 0.45), L(0.5, 0), fa, thin)
                shut > 0 && ig.AddTriangleFilled(dl, L(-0.6, -0.45), L(-0.6, 0.45), L(0.5, 0), fa)
                line((0.55, -0.45), (0.55, 0.45), tinted ? lb : c)
            else
                for (side, f) in ((-1, fa), (1, fb))
                    ig.AddTriangle(dl, L(0.6side, -0.4), L(0.6side, 0.4), L(0, 0), f, thin)
                    shut > 0 && ig.AddTriangleFilled(dl, L(0.6side * shut, -0.4shut), L(0.6side * shut, 0.4shut), L(0, 0), f)
                end
                live && shut == 1 && line((0, -0.58), (0, 0.58), c, 2f0)
            end
            stem = tinted && live ? state : c                  # the actuator shows open / shut when the body shows pressure
            if k == "Valve"
                # A click changes the command at once, but the body shows where the valve is, and
                # that only moves with simulated time: paused, a click would seem to do nothing.
                # So the actuator takes the colour of where the valve has been told to go and carries
                # an arrow until it gets there: outline pointing away = opening, filled pointing at
                # the body = closing (the same shapes as open and shut).
                want = live ? pending(edge_of(app, part.name)) : 0
                if want != 0
                    stem = want > 0 ? C_OPEN : C_SHUT
                    want > 0 ? ig.AddTriangle(dl, L(-0.28, -0.85), L(0.28, -0.85), L(0, -1.3), stem, thin) :
                               ig.AddTriangleFilled(dl, L(-0.28, -1.3), L(0.28, -1.3), L(0, -0.85), stem)
                end
                heavy = tinted || want != 0 ? 3f0 : 2f0
                line((0, 0), (0, -0.7), stem, heavy); line((-0.3, -0.7), (0.3, -0.7), stem, heavy)
                part.p["normally_open"] > 0.5 && ig.AddText(dl, L(0.35, -1.0), C_DIM, "NO")
            elseif k == "Regulator"
                line((0, 0), (0, -0.5), stem); ig.AddCircle(dl, L(0, -0.75), 0.27f0 * grid(app), stem, 0, 2f0 * UI[])
                line((0.2, 0.5), (0.8, 0.5)); line((0.8, 0.5), (0.6, 0.35)); line((0.8, 0.5), (0.6, 0.65))
            elseif k == "ReliefValve"
                line((0, 0), (0, -0.3), stem); line((0, -0.3), (-0.2, -0.45), stem); line((-0.2, -0.45), (0.2, -0.6), stem)
                line((0.2, -0.6), (-0.2, -0.75), stem); line((-0.2, -0.75), (0, -0.9), stem)
            end
        end
    elseif k == "Atmosphere"
        line((-1, 0), (-0.2, 0)); line((-0.2, -0.5), (-0.2, 0.5)); line((0.05, -0.33), (0.05, 0.33)); line((0.3, -0.16), (0.3, 0.16))
    elseif k == "Junction"
        ig.AddCircleFilled(dl, L(0, 0), 0.16f0 * grid(app), live ? something(pin_tint(app, part.name, "a"), c) : c)
    elseif k == "Gauge"
        g = live ? something(pin_tint(app, part.name, "a"), c) : c
        line((0, 1), (0, 0.45), g); ig.AddCircle(dl, L(0, 0), 0.45f0 * grid(app), selected ? C_SEL : g, 0, 2f0 * UI[])
        ig.AddText(dl, (L(0, 0)[1] - px(4), L(0, 0)[2] - px(8)), c, "P")
    end
    # Label, with live readings when running.
    lx, ly = k in VESSEL_KINDS ? w2s(app, part.x + 1.05, part.y - (k == "Tank" ? 0.4 : 1.2)) :
             k == "Gauge" ? w2s(app, part.x + 0.6, part.y - 0.65) :
             k == "Valve" && mod(part.rot, 4) == 1 ? w2s(app, part.x + 0.75, part.y + 0.3) :   # its actuator and arrow point this way: keep clear
             k == "Valve" && mod(part.rot, 4) == 2 ? w2s(app, part.x - 0.9, part.y + 1.4) :    # likewise, pointing down
             isodd(part.rot) ? w2s(app, part.x + 0.75, part.y - 0.3) : w2s(app, part.x - 0.9, part.y + 0.55)
    k == "Junction" || ig.AddText(dl, (lx, ly), selected ? C_SEL : C_DIM, part.name)
    live || return
    sim, net = app.sim, app.sim.net
    put(row, s, color = C_LINE) = ig.AddText(dl, (lx, ly + px(16) * row), color, s)
    if k in VESSEL_KINDS
        r = REPORTS[node_index(app, part.name)]
        over = (r.P - p_atm(app)) / PSI > app.scene.settings["meop_psig"]
        put(1, fmt(app, :pressure, r.P), over ? C_WARN : C_LINE)
        put(2, fmt(app, :temperature, r.T) * " liquid")
        put(3, fmt(app, :mass, r.m))
        # Two different things: how much of the vessel the liquid occupies, and how much of the
        # nitrous in it is liquid. Venting a closed tank lowers both.
        put(4, r.liquid_frac >= 0.999 ? "LIQUID-FULL" : "liquid $(round(Int, 100 * r.liquid_frac)) % of volume", r.liquid_frac >= 0.999 ? C_WARN : C_LINE)
        put(5, "liquid $(round(Int, 100 * (1 - r.x))) % of mass")
        node_of(app, part.name).heater_on && put(6, "HEATER ON", col(255, 160, 60))
    elseif k == "Gauge"
        g = get(net.gauges, part.name, 0)
        g == 0 || put(1, fmt(app, :pressure, REPORTS[g].P))
    elseif k in FLOW_KINDS
        e = edge_of(app, part.name)
        abs(e.mdot) > 3e-4 && k != "Pipe" && k != "Orifice" && put(1, "$(e.mdot > 0 ? ">" : "<") " * fmt(app, :flow, abs(e.mdot)), col(140, 220, 255))
    end
    nothing
end

# Scale colour of a pin's node as RGB (for blending along a pipe).
function scale_rgb(app::App, part, pin)
    i = get(app.sim.net.pin_node, (part, pin), 0)
    lo, hi = app.shown_range
    i == 0 || OFF_SCALE[i] ? RGB((120, 120, 130)) : scale_colour((VALUES[i] - lo) / (hi - lo))
end

"What the hover tooltip says about a part."
function describe(app::App, part::Part)
    k = part.kind
    lines = [part.name * "  (" * k * ")"]
    if !running(app)
        for p in KINDS[k].params[1:min(end, 5)]
            push!(lines, "$(p.key) = $(part.p[p.key])")
        end
        return join(lines, "\n")
    end
    net = app.sim.net
    if k in VESSEL_KINDS
        r = REPORTS[node_index(app, part.name)]
        push!(lines, "pressure " * fmt(app, :pressure, r.P) * (isnan(r.Psat) ? "" : ",  liquid's vapour pressure " * fmt(app, :pressure, r.Psat)))
        push!(lines, "liquid " * fmt(app, :temperature, r.T) * ",  ullage " * fmt(app, :temperature, r.Tv) * ",  wall " * fmt(app, :temperature, r.Tw))
        push!(lines, fmt(app, :mass, r.m) * ",  liquid $(round(Int, 100 * r.liquid_frac)) % of volume, $(round(Int, 100 * (1 - r.x))) % of mass")
    elseif k in FLOW_KINDS
        e = edge_of(app, part.name)
        ra, rb = REPORTS[e.a], REPORTS[e.b]
        up = e.mdot >= 0 ? ra : rb
        push!(lines, "flow " * fmt(app, :flow, abs(e.mdot)) * (e.mdot > 0 ? "  a to b" : e.mdot < 0 ? "  b to a" : ""))
        push!(lines, "a: " * fmt(app, :pressure, ra.P) * "   b: " * fmt(app, :pressure, rb.P) * "   drop $(round((ra.P - rb.P) / PSI, digits = 1)) psi")
        abs(e.mdot) > 0 && push!(lines, "upstream " * fmt(app, :temperature, up.T) * ", $(round(Int, 100 * up.x)) % vapour by mass")
        k == "Valve" && push!(lines, "position $(round(Int, 100 * e.opening)) % open" * (e.stuck != 0 ? " (stuck)" : ""))
    else
        i = k == "Atmosphere" ? node_index(app, "atmosphere") : get(net.pin_node, (part.name, "a"), 0)
        if i != 0
            r = REPORTS[i]
            push!(lines, fmt(app, :pressure, r.P) * ",  " * fmt(app, :temperature, r.T) * ",  $(round(Int, 100 * r.x)) % vapour by mass")
            CUT_OFF[i] && push!(lines, "shut off on all sides: holding its last pressure")
        end
    end
    join(lines, "\n")
end

function draw_canvas(app::App)
    dl = ig.GetWindowDrawList()
    p0 = ig.GetCursorScreenPos()
    avail = ig.GetContentRegionAvail()
    app.canvas = (p0.x, p0.y, max(avail.x, px(50)), max(avail.y, px(50)))
    app.want_fit && zoom_fit!(app)
    ig.InvisibleButton("canvas", (app.canvas[3], app.canvas[4]), ig.ImGuiButtonFlags_MouseButtonLeft | ig.ImGuiButtonFlags_MouseButtonRight)
    hovered = ig.IsItemHovered()
    io = ig.GetIO()
    mouse = ig.GetMousePos()
    wx, wy = s2w(app, mouse.x, mouse.y)
    live = running(app)
    live && survey!(app)
    ig.PushClipRect(dl, (p0.x, p0.y), (p0.x + avail.x, p0.y + avail.y), true)
    step = grid(app)                                           # grid dots
    if step > px(10)
        x = p0.x + mod(app.origin[1] - p0.x, step)
        while x < p0.x + avail.x
            y = p0.y + mod(app.origin[2] - p0.y, step)
            while y < p0.y + avail.y
                ig.AddRectFilled(dl, (x - 1f0, y - 1f0), (x + 1f0, y + 1f0), col(60, 60, 70))
                y += step
            end
            x += step
        end
    end
    # The pin or wire corner a click would take now is drawn lit, so that a target grown to a
    # usable size is not a surprise. A pin wins over a corner, as it does on the click.
    aim = hovered && !live && isempty(app.placing) && !ig.IsMouseDown(0)
    hot_pin = aim ? pin_at(app, wx, wy) : nothing
    hot_via = aim && hot_pin === nothing ? via_at(app, wx, wy) : nothing
    for (i, w) in enumerate(app.scene.wires)
        c = live ? something(pin_tint(app, w.a[1], w.a[2]), C_WIRE) : C_WIRE
        t = (live && app.colour_by !== :none ? 3f0 : 2f0) * UI[]
        each_run((x1, y1, x2, y2, _) -> (ig.AddLine(dl, w2s(app, x1, y1), w2s(app, x2, y2), c, t); false), app, w)
        live && continue
        for (k, v) in enumerate(w.via)                         # corners can be dragged while editing
            x, y = w2s(app, v...)
            hot = hot_via == (i, k)
            r = px(hot ? 5 : 3)
            ig.AddRectFilled(dl, (x - r, y - r), (x + r, y + r), hot ? C_SEL : C_WIRE)
        end
    end
    for part in app.scene.parts
        draw_part(dl, app, part)
    end
    if !live
        for part in app.scene.parts, pin in KINDS[part.kind].pins
            ig.AddCircle(dl, w2s(app, pin_pos(part, pin)...), px(4), app.wire_from == (part.name, pin.name) ? C_SEL : col(120, 200, 140), 0, 1.5f0 * UI[])
            hot_pin == (part.name, pin.name) && ig.AddCircleFilled(dl, w2s(app, pin_pos(part, pin)...), px(5), C_SEL)
        end
        if app.wire_from !== nothing
            i = find_part(app.scene, app.wire_from[1])
            if i !== nothing
                part = app.scene.parts[i]
                pin = KINDS[part.kind].pins[findfirst(p -> p.name == app.wire_from[2], KINDS[part.kind].pins)]
                ig.AddLine(dl, w2s(app, pin_pos(part, pin)...), (mouse.x, mouse.y), C_SEL, 1.5f0 * UI[])
            end
        end
        isempty(app.placing) || ig.AddText(dl, (mouse.x + px(12), mouse.y + px(8)), C_SEL, "place " * app.placing)
        if app.box_from !== nothing                            # selection box being dragged out
            x, y = w2s(app, app.box_from...)
            ig.AddRect(dl, (min(x, mouse.x), min(y, mouse.y)), (max(x, mouse.x), max(y, mouse.y)), C_SEL, 0f0, 1f0 * UI[])
        end
    elseif app.colour_by !== :none
        colour_bar(dl, app)
    end
    ig.PopClipRect(dl)

    hovered || return
    wheel = unsafe_load(io.MouseWheel)
    if wheel != 0                                              # zoom about the cursor
        z = clamp(app.zoom * (1 + 0.12f0 * wheel), 0.35f0, 3f0)
        app.origin = (mouse.x - Float32(wx) * px(26) * z, mouse.y - Float32(wy) * px(26) * z)
        app.zoom = z
    end
    if ig.IsMouseDragging(1) || ig.IsMouseDragging(2)          # pan
        d = unsafe_load(io.MouseDelta)
        app.origin = (app.origin[1] + d.x, app.origin[2] + d.y)
    end
    over = part_at(app, wx, wy)
    over !== nothing && !app.dragging && !ig.IsMouseDown(0) && isempty(app.placing) && app.wire_from === nothing && ig.SetTooltip(replace(describe(app, over), "%" => "%%"))
    live ? canvas_run_clicks(app, wx, wy) : canvas_edit_clicks(app, wx, wy)
end

# The right button both pans (drag) and acts (cut a wire, remove a corner, open a part's menu).
# Acting on the press made every pan that happened to start on a wire cut it, so a right-click
# counts when the button comes up without having been dragged anywhere.
right_clicked() = ig.IsMouseReleased(1) && !ig.IsMouseDragPastThreshold(1)

function canvas_edit_clicks(app::App, wx, wy)
    s = app.scene
    gx, gy = round(Int, wx), round(Int, wy)
    shift = unsafe_load(ig.GetIO().KeyShift)
    if ig.IsMouseClicked(0)
        app.dragging, app.drag_via, app.box_from, app.drag_last = false, (0, 0), nothing, (gx, gy)
        pin = pin_at(app, wx, wy)
        if !isempty(app.placing)
            snapshot!(app)
            select!(app, add_part!(s, app.placing, gx, gy).name)
            app.placing = ""
        elseif pin !== nothing
            if app.wire_from === nothing
                app.wire_from = pin
            else
                snapshot!(app)
                connect!(s, app.wire_from, pin)
                app.wire_from = nothing
            end
        else
            app.wire_from = nothing
            via, part, wire = via_at(app, wx, wy), part_at(app, wx, wy), wire_at(app, wx, wy)
            if via !== nothing
                app.drag_via = via                             # an existing corner: drag it
            elseif part !== nothing
                if shift                                       # shift-click adds to or removes from the selection
                    part.name in app.sel ? filter!(!=(part.name), app.sel) : push!(app.sel, part.name)
                elseif part.name in app.sel                    # already selected: keep the group, make this the current part
                    filter!(!=(part.name), app.sel); push!(app.sel, part.name)
                else
                    select!(app, part.name)
                end
                isempty(app.sel) || setbuf!(app.name_buf, current(app))
            elseif wire !== nothing                            # a click on a wire puts a corner there, ready to drag
                snapshot!(app)
                insert!(s.wires[wire[1]].via, wire[2], (gx, gy))
                app.drag_via, app.dragging = wire, true
            else                                               # empty space: start a selection box
                shift || empty!(app.sel)
                app.box_from = (wx, wy)
            end
        end
    elseif ig.IsMouseDragging(0)
        if app.drag_via != (0, 0)
            w = s.wires[app.drag_via[1]]
            if w.via[app.drag_via[2]] != (gx, gy)
                app.dragging || snapshot!(app)
                app.dragging = true
                w.via[app.drag_via[2]] = (gx, gy)
            end
        elseif app.box_from === nothing && !isempty(app.sel) && app.wire_from === nothing
            dx, dy = gx - app.drag_last[1], gy - app.drag_last[2]
            if (dx, dy) != (0, 0)                              # the whole selection moves together
                app.dragging || snapshot!(app)
                app.dragging = true
                move_parts!(s, app.sel, dx, dy)
                app.drag_last = (gx, gy)
            end
        end
    elseif ig.IsMouseReleased(0)
        if app.box_from !== nothing
            x0, x1 = minmax(app.box_from[1], wx)
            y0, y1 = minmax(app.box_from[2], wy)
            for p in s.parts
                x0 <= p.x <= x1 && y0 <= p.y <= y1 && !(p.name in app.sel) && push!(app.sel, p.name)
            end
            isempty(app.sel) || setbuf!(app.name_buf, current(app))
        end
        app.dragging, app.drag_via, app.box_from = false, (0, 0), nothing
    elseif right_clicked()
        via, wire = via_at(app, wx, wy), wire_at(app, wx, wy)
        if app.wire_from !== nothing || !isempty(app.placing)
            app.wire_from, app.placing = nothing, ""
        elseif via !== nothing                                 # right-click a corner: remove it
            snapshot!(app)
            deleteat!(s.wires[via[1]].via, via[2])
        elseif wire !== nothing && part_at(app, wx, wy) === nothing
            snapshot!(app)
            disconnect!(s, s.wires[wire[1]])
        end
    end
    nothing
end

function canvas_run_clicks(app::App, wx, wy)
    (ig.IsMouseClicked(0) || right_clicked()) || return
    part = part_at(app, wx, wy)
    select!(app, part === nothing ? "" : part.name)
    part === nothing && return
    if ig.IsMouseClicked(0) && part.kind == "Valve"
        act!(app.sim, part.name, "toggle")
    elseif !ig.IsMouseClicked(0)
        ig.OpenPopup("partmenu")
    end
    nothing
end
