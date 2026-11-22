# The scene: what is on the schematic and how it is wired. This is the file format, the editor's
# data model and the input to the simulator (network.jl), so there is exactly one description of a
# plumbing system. Parameter names carry their unit (`d_mm`, `set_psig`); conversion to SI happens
# in one place, `compile` in network.jl.

import TOML

"One parameter of a part kind: key (with unit suffix), default, and a one-line explanation."
struct Param
    key::String
    default::Float64
    help::String
end

"A connection point on a part, at (dx, dy) grid units from the part centre before rotation."
struct Pin
    name::String
    dx::Int
    dy::Int
end

struct Kind
    pins::Vector{Pin}
    params::Vector{Param}
    help::String
end

const _VESSEL_PARAMS = [
    Param("volume_L", 19.16, "internal volume"),
    Param("diameter_mm", 143.0, "inside diameter (vertical cylinder; sets level and wall areas)"),
    Param("fill_kg", 0.0, "N2O mass at the start"),
    Param("T_C", 25.0, "starting temperature of contents and wall"),
    Param("air", 1.0, "1 = the rest of the volume starts as air at atmospheric pressure, 0 = N2O only"),
    Param("wall_kg", 6.0, "mass of wall metal that follows the contents' temperature"),
    Param("wall_cp_J_kgK", 900.0, "wall specific heat (aluminium 900, steel 470)"),
    Param("absorptivity", 0.4, "fraction of solar flux absorbed"),
    Param("emissivity", 0.8, "infrared emissivity of the outside surface"),
    Param("shade", 0.0, "fraction of sunlight blocked, 0..1"),
    Param("heater_W", 0.0, "heater power delivered to the wall when switched on"),
]

_two_pin(params, help) = Kind([Pin("a", -1, 0), Pin("b", 1, 0)], params, help)
const _RESTRICTION = [
    Param("d_mm", 6.0, "bore diameter"),
    Param("Cd", 0.7, "discharge coefficient (calibration knob)"),
]

"Every kind of part that can be placed, in palette order."
const KINDS = Dict{String,Kind}(
    "Bottle" => Kind([Pin("out", 0, -2)],
        [Param("volume_L", 40.0, "internal (water) volume"),
         Param("diameter_mm", 230.0, "inside diameter"),
         Param("fill_kg", 25.0, "N2O mass at the start"),
         Param("T_C", 25.0, "starting temperature"),
         Param("draw_height_frac", 0.0, "height of the outlet inside the bottle: 0 = liquid draw (siphon tube or inverted), 1 = vapour draw"),
         Param("wall_kg", 30.0, "cylinder metal mass"),
         Param("wall_cp_J_kgK", 470.0, "wall specific heat (aluminium 900, steel 470)"),
         Param("absorptivity", 0.5, "fraction of solar flux absorbed"),
         Param("emissivity", 0.8, "infrared emissivity"),
         Param("shade", 0.0, "fraction of sunlight blocked, 0..1"),
         Param("heater_W", 0.0, "heater power when switched on")],
        "Supply cylinder: a two-phase vessel that cools as it delivers."),
    "Tank" => Kind([Pin("top", 0, -3), Pin("vent", 1, -2), Pin("bottom", 0, 3)],
        vcat(_VESSEL_PARAMS, [Param("vent_height_frac", 0.9, "height of the `vent` pin's dip-tube inlet, as a fraction of tank height; liquid reaches the vent at this level")]),
        "Run tank. `bottom` draws liquid, `top` draws ullage, `vent` switches from vapour to liquid when the level reaches the dip tube."),
    "Pipe" => _two_pin(
        [Param("L_m", 1.0, "length"),
         Param("id_mm", 10.0, "inside diameter"),
         Param("rough_mm", 0.015, "wall roughness"),
         Param("K", 0.0, "sum of minor-loss coefficients (bends, fittings)")],
        "Tube or hose: wall friction plus minor losses. It holds no fluid of its own (see docs/MODEL.md)."),
    "Valve" => _two_pin(
        vcat(_RESTRICTION, [
            Param("Cv", 0.0, "US flow coefficient; if > 0 it replaces d_mm and Cd"),
            Param("normally_open", 0.0, "1 = open when de-energised (fail state)"),
            Param("travel_s", 0.3, "time to stroke fully")]),
        "On/off valve. Click to command it in run mode."),
    "Orifice" => _two_pin(_RESTRICTION, "Fixed restriction (orifice, quick-disconnect, fitting)."),
    "Regulator" => _two_pin(
        vcat(_RESTRICTION, [
            Param("set_psig", 700.0, "outlet pressure it closes at"),
            Param("droop_psi", 30.0, "outlet pressure drop below set-point at which it is fully open")]),
        "Pressure-reducing regulator, flow a → b only. Saturated liquid flashes across it."),
    "CheckValve" => _two_pin(
        vcat(_RESTRICTION, [Param("crack_psi", 1.0, "cracking pressure")]),
        "Check valve, flow a → b only."),
    "ReliefValve" => _two_pin(
        vcat(_RESTRICTION, [Param("set_psig", 1000.0, "inlet pressure at which it starts to open (fully open 10 % above)")]),
        "Relief valve, flow a → b only. Wire b to Atmosphere."),
    "Atmosphere" => Kind([Pin("a", -1, 0)], Param[], "Ambient boundary. Everything vented ends up here."),
    "Junction" => Kind([Pin("a", 0, 0)], Param[], "A dot to join several wires (tee, cross)."),
    "Gauge" => Kind([Pin("a", 0, 1)], Param[], "Shows pressure and temperature of whatever it is wired to (PT / TC)."),
)
const PALETTE = ["Bottle", "Tank", "Pipe", "Valve", "Orifice", "Regulator", "CheckValve", "ReliefValve", "Atmosphere", "Junction", "Gauge"]

const SETTINGS = [
    Param("T_amb_C", 25.0, "ambient air temperature"),
    Param("solar_W_m2", 0.0, "solar flux on exposed surfaces (clear summer noon ≈ 900)"),
    Param("wind_m_s", 1.0, "wind speed"),
    Param("P_atm_kPa", 101.325, "local atmospheric pressure; gauge readings are relative to this"),
    Param("target_lo_psig", 0.0, "bottom of a tank pressure window to judge the run against (leave both at 0 for none)"),
    Param("target_hi_psig", 0.0, "top of that window"),
    Param("meop_psig", 1000.0, "maximum expected operating pressure (warning threshold)"),
    Param("interphase_E", 700.0, "boiling enhancement on liquid-to-surface heat transfer (Zilliac-Karabeyoglu E; fitted, literature 240-1300)"),
    Param("cond_tau_s", 0.5, "time constant for supersaturated ullage vapour to condense"),
    Param("h_in_mult", 1.0, "multiplier on wall-to-fluid natural convection inside vessels"),
    Param("h_ext_mult", 1.0, "multiplier on outside convection (0 = perfectly insulated)"),
]

mutable struct Part
    kind::String
    name::String
    x::Int
    y::Int
    rot::Int                               # quarter turns clockwise
    p::Dict{String,Float64}
end

"""
A connection between two pins. `via` lists grid points the drawn wire passes through, in order
from `a` to `b`; it only affects how the wire is drawn. Between consecutive points the wire runs
horizontally, then vertically.
"""
struct Wire
    a::Tuple{String,String}                # (part name, pin name)
    b::Tuple{String,String}
    via::Vector{Tuple{Int,Int}}
end
Wire(a, b) = Wire(a, b, Tuple{Int,Int}[])

"A scripted action: at time `t`, `open`/`close` a valve, switch a heater `on`/`off`, or `set` a regulator."
struct Action
    t::Float64
    part::String
    action::String
    value::Float64
end

mutable struct Scene
    settings::Dict{String,Float64}
    parts::Vector{Part}
    wires::Vector{Wire}
    schedule::Vector{Action}
end

Scene() = Scene(Dict(p.key => p.default for p in SETTINGS), Part[], Wire[], Action[])

defaults(kind) = Dict(p.key => p.default for p in KINDS[kind].params)
find_part(s::Scene, name) = findfirst(p -> p.name == name, s.parts)

"Pin position on the canvas, in grid units, after the part's rotation."
function pin_pos(part::Part, pin::Pin)
    dx, dy = pin.dx, pin.dy
    for _ in 1:mod(part.rot, 4)
        dx, dy = -dy, dx
    end
    part.x + dx, part.y + dy
end

# ---------------------------------------------------------------- edit operations
# Plain functions on the scene, so the editor is testable without a window.

function add_part!(s::Scene, kind, x, y; name = nothing)
    haskey(KINDS, kind) || error("unknown part kind $kind")
    if name === nothing
        prefix = Dict("Bottle" => "BOTTLE", "Tank" => "TANK", "Pipe" => "L", "Valve" => "V", "Orifice" => "OR",
                      "Regulator" => "REG", "CheckValve" => "CV", "ReliefValve" => "PRV", "Atmosphere" => "ATM",
                      "Junction" => "J", "Gauge" => "PT")[kind]
        n = 1
        while find_part(s, "$prefix$n") !== nothing
            n += 1
        end
        name = "$prefix$n"
    end
    find_part(s, name) === nothing || error("a part named $name already exists")
    push!(s.parts, Part(kind, name, x, y, 0, defaults(kind)))
    s.parts[end]
end

function delete_part!(s::Scene, name)
    filter!(p -> p.name != name, s.parts)
    filter!(w -> w.a[1] != name && w.b[1] != name, s.wires)
    filter!(a -> a.part != name, s.schedule)
    s
end

function rename_part!(s::Scene, old, new)
    old == new && return s
    (isempty(new) || occursin('.', new)) && error("a part name cannot be empty or contain '.'")
    find_part(s, new) === nothing || error("a part named $new already exists")
    s.parts[find_part(s, old)].name = new
    fix(pin) = pin[1] == old ? (new, pin[2]) : pin
    s.wires = [Wire(fix(w.a), fix(w.b), w.via) for w in s.wires]
    s.schedule = [a.part == old ? Action(a.t, new, a.action, a.value) : a for a in s.schedule]
    s
end

function connect!(s::Scene, a::Tuple{String,String}, b::Tuple{String,String})
    a == b && return s
    any(w -> (w.a == a && w.b == b) || (w.a == b && w.b == a), s.wires) || push!(s.wires, Wire(a, b))
    s
end

disconnect!(s::Scene, w::Wire) = (filter!(x -> x !== w && !(x.a == w.a && x.b == w.b), s.wires); s)

"""
    copy_parts(scene, names) -> (parts, wires)

Independent copies of the named parts and of the wires that run between them, for pasting.
"""
function copy_parts(s::Scene, names)
    keep = Set(names)
    [deepcopy(p) for p in s.parts if p.name in keep], [deepcopy(w) for w in s.wires if w.a[1] in keep && w.b[1] in keep]
end

"""
    paste_parts!(scene, (parts, wires), dx, dy) -> names of the new parts

Add copies of `parts`, moved by (dx, dy) grid units and given fresh names, with the wires between
them reconnected to the copies.
"""
function paste_parts!(s::Scene, clip, dx, dy)
    parts, wires = clip
    renamed = Dict{String,String}()
    for p in parts
        q = add_part!(s, p.kind, p.x + dx, p.y + dy)
        q.rot, q.p = p.rot, copy(p.p)
        renamed[p.name] = q.name
    end
    for w in wires
        connect!(s, (renamed[w.a[1]], w.a[2]), (renamed[w.b[1]], w.b[2]))
        append!(s.wires[end].via, [(x + dx, y + dy) for (x, y) in w.via])
    end
    [renamed[p.name] for p in parts]
end

"Move the named parts by (dx, dy) grid units, taking the corners of wires between them along."
function move_parts!(s::Scene, names, dx, dy)
    keep = Set(names)
    for p in s.parts
        p.name in keep && ((p.x, p.y) = (p.x + dx, p.y + dy))
    end
    for w in s.wires
        w.a[1] in keep && w.b[1] in keep && (w.via .= [(x + dx, y + dy) for (x, y) in w.via])
    end
    s
end

# ---------------------------------------------------------------- nets

"""
Group pins into nets (sets of pins joined by wires). Returns a Dict from (part, pin) to net number
and the number of nets. Every pin of every part gets a net, wired or not.
"""
function nets(s::Scene)
    pins = [(p.name, pin.name) for p in s.parts for pin in KINDS[p.kind].pins]
    id = Dict(pin => i for (i, pin) in enumerate(pins))
    parent = collect(1:length(pins))
    root(i) = (while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
    for w in s.wires
        haskey(id, w.a) && haskey(id, w.b) || continue
        parent[root(id[w.a])] = root(id[w.b])
    end
    number = Dict{Int,Int}()
    out = Dict{Tuple{String,String},Int}()
    for pin in pins
        out[pin] = get!(number, root(id[pin]), length(number) + 1)
    end
    out, length(number)
end

const FLOW_KINDS = ("Pipe", "Valve", "Orifice", "Regulator", "CheckValve", "ReliefValve")
const VESSEL_KINDS = ("Bottle", "Tank")

"Problems that stop a scene from running (errors) or that are probably mistakes (warnings)."
function check(s::Scene)
    errors, warnings = String[], String[]
    net, n = nets(s)
    kind = Dict(p.name => p.kind for p in s.parts)
    members = [Tuple{String,String}[] for _ in 1:n]
    for (pin, k) in net
        push!(members[k], pin)
    end
    any(p -> p.kind in VESSEL_KINDS, s.parts) || push!(errors, "there is no Bottle or Tank")
    for m in members
        vessels = unique(first.(filter(pin -> kind[pin[1]] in VESSEL_KINDS, m)))
        atm = any(pin -> kind[pin[1]] == "Atmosphere", m)
        length(vessels) > 1 && push!(errors, "$(join(vessels, " and ")) are wired straight together; put a valve, pipe or orifice between them")
        !isempty(vessels) && atm && push!(errors, "$(vessels[1]) is wired straight to Atmosphere; put a valve or orifice between them")
        !isempty(vessels) && count(pin -> kind[pin[1]] == vessels[1], m) > 1 &&
            push!(warnings, "two pins of $(vessels[1]) are wired together")
        if length(m) == 1 && kind[m[1][1]] in FLOW_KINDS
            push!(warnings, "$(m[1][1]).$(m[1][2]) is not connected (treated as a capped dead end)")
        end
    end
    for p in s.parts, q in KINDS[p.kind].params
        v = get(p.p, q.key, q.default)
        isfinite(v) || push!(errors, "$(p.name).$(q.key) is not a number")
        q.key in ("volume_L", "diameter_mm", "d_mm", "id_mm", "L_m", "wall_cp_J_kgK") && v <= 0 &&
            push!(errors, "$(p.name).$(q.key) must be positive")
    end
    for a in s.schedule
        find_part(s, a.part) === nothing && push!(warnings, "schedule refers to missing part $(a.part)")
    end
    errors, warnings
end

# ---------------------------------------------------------------- file I/O

function load_scene(path)
    d = TOML.parsefile(path)
    s = Scene()
    for (k, v) in get(d, "settings", Dict())
        s.settings[k] = Float64(v)
    end
    for t in get(d, "part", [])
        kind, name = t["kind"], t["name"]
        haskey(KINDS, kind) || error("$path: unknown part kind \"$kind\" for $name")
        part = add_part!(s, kind, Int(get(t, "x", 0)), Int(get(t, "y", 0)); name)
        part.rot = Int(get(t, "rot", 0))
        for (k, v) in t
            k in ("kind", "name", "x", "y", "rot") && continue
            haskey(part.p, k) || error("$path: $name has no parameter \"$k\"")
            part.p[k] = Float64(v)
        end
    end
    pin(str) = (i = findlast('.', str); (String(str[1:i-1]), String(str[i+1:end])))
    for w in get(d, "wire", [])
        n = length(s.wires)
        connect!(s, pin(w["from"]), pin(w["to"]))
        length(s.wires) > n && append!(s.wires[end].via, [(Int(v[1]), Int(v[2])) for v in get(w, "via", [])])
    end
    for a in get(d, "schedule", [])
        push!(s.schedule, Action(Float64(a["t"]), a["part"], a["action"], Float64(get(a, "value", 0.0))))
    end
    sort!(s.schedule; by = a -> a.t)
    s
end

_num(v) = v == round(v) && abs(v) < 1e9 ? string(Int(v), ".0") : string(v)

function save_scene(path, s::Scene)
    open(path, "w") do io
        println(io, "# n2o-fill-sim scene. Units are in the parameter names. See README.md for the format.\n")
        println(io, "[settings]")
        for p in SETTINGS
            println(io, rpad(p.key, 16), " = ", rpad(_num(s.settings[p.key]), 10), " # ", p.help)
        end
        for part in s.parts
            println(io, "\n[[part]]")
            println(io, "kind = \"", part.kind, "\"\nname = \"", part.name, "\"")
            println(io, "x = ", part.x, "\ny = ", part.y, "\nrot = ", part.rot)
            for p in KINDS[part.kind].params
                println(io, p.key, " = ", _num(part.p[p.key]))
            end
        end
        for w in s.wires
            println(io, "\n[[wire]]\nfrom = \"", w.a[1], ".", w.a[2], "\"\nto = \"", w.b[1], ".", w.b[2], "\"")
            isempty(w.via) || println(io, "via = [", join(["[$x, $y]" for (x, y) in w.via], ", "), "]")
        end
        for a in s.schedule
            println(io, "\n[[schedule]]\nt = ", _num(a.t), "\npart = \"", a.part, "\"\naction = \"", a.action, "\"")
            a.action == "set" && println(io, "value = ", _num(a.value))
        end
    end
end
