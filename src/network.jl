# Scene → equations → time stepping.
#
# The network is a graph. Vessels (bottles, tanks) hold fluid and carry the state. Flow parts are
# edges: each turns the pressures either side of it and the upstream enthalpy into a mass flow.
# Where wires join without a vessel there is a junction. A junction holds no fluid: its pressure is
# whatever makes the flows into it balance, found afresh every time the equations are evaluated,
# exactly as a circuit simulator finds node voltages. Fluid crossing a junction keeps its enthalpy
# (streams mix), so flashing through a regulator or a long hose needs no special treatment.
#
# Junctions were first given their own volume and fluid state, which would also show a hose
# charging and liquid trapped between closed valves. That made the equations unsolvable in
# practice: millilitres of near-incompressible liquid sitting on the two-phase boundary. To study
# trapped liquid, put a small Tank where the trapped volume is.

using LinearAlgebra: lu!, ldiv!
using OrdinaryDiffEqRosenbrock: Rodas5P
using SciMLBase: ODEFunction, ODEProblem, init, step!, reinit!, add_tstop!, successful_retcode, ReturnCode

const PSI = 6894.757293168
const EDGE_KIND = Dict("Pipe" => :pipe, "Valve" => :valve, "Orifice" => :orifice, "Regulator" => :regulator,
                       "CheckValve" => :check, "ReliefValve" => :relief)

mutable struct Edge
    name::String
    kind::Symbol                  # :pipe, :valve, :orifice, :regulator, :check, :relief
    a::Int                        # node index at pin a, and which of that node's ports
    pa::Int
    b::Int
    pb::Int
    area::Float64                 # Cd·A fully open, m² (for a pipe: bore area)
    set::Float64                  # Pa: regulator / relief set-point (absolute), check cracking ΔP
    band::Float64                 # Pa over which it goes from shut to fully open
    L::Float64                    # pipe only
    D::Float64
    relrough::Float64
    K::Float64
    opening::Float64              # 0..1, valves
    travel::Float64               # s to stroke
    normally_open::Bool
    cmd::Bool                     # energised (commanded away from its fail state)
    stuck::Int                    # 0 = healthy, 1 = stuck open, 2 = stuck closed
    mdot::Float64                 # kg/s, positive a → b
    # Memo of the last upstream state resolved for this part and its choking point (see hem_flux).
    m_P1::Float64
    m_h1::Float64
    m_up::Fluid
    m_Pc::Float64
    m_G::Float64
    m_Pu::Float64
end

mutable struct Net
    nodes::Vector{Node}
    edges::Vector{Edge}
    env::Env
    nfluid::Int                   # states belonging to vessels; three bookkeeping states follow
    junctions::Vector{Int}        # node indices
    adj::Vector{Vector{Int}}      # edges touching each node
    gauges::Dict{String,Int}      # Gauge part name → node index
    pin_node::Dict{Tuple{String,String},Int}
    # work space for the junction solve
    order::Vector{Int}
    res::Vector{Float64}
    p0::Vector{Float64}
    step::Vector{Float64}
    jac::Matrix{Float64}
    row::Vector{Int}              # node index → row in the junction system (0 if not a junction)
    lu::Any                       # factorised junction Jacobian from the last full Newton step
    unconverged::Int              # evaluations in which the junction solve did not converge
    # which parts can carry flow and which junctions need solving right now (see plan!)
    live::Vector{Bool}
    active::Vector{Int}
    pruned::Vector{Int}
    feeder::Vector{Int}
    deg::Vector{Int}
    stack::Vector{Int}
    stamp::UInt
end

"Is this valve meant to be open right now?"
valve_target(e::Edge) = e.stuck == 1 ? 1.0 : e.stuck == 2 ? 0.0 : (e.cmd != e.normally_open ? 1.0 : 0.0)

function Env(settings::Dict)
    g(k) = settings[k]
    Env(g("T_amb_C") + 273.15, g("solar_W_m2"), g("wind_m_s"), g("P_atm_kPa") * 1e3, g("interphase_E"),
        g("cond_tau_s"), g("h_in_mult"), g("h_ext_mult"))
end

"""
    compile(scene) -> (net, u0)

Turn a scene into a network and its initial state. Each net of wires becomes one node: the vessel
whose pin is on it, the atmosphere, or a junction. Throws if `check(scene)` finds errors.
"""
function compile(s::Scene)
    errors, _ = check(s)
    isempty(errors) || error(join(errors, "; "))
    env = Env(s.settings)
    net_of, nnet = nets(s)
    nodes = Node[]
    node_of_net = zeros(Int, nnet)
    port_of_net = ones(Int, nnet)             # which vessel port a net lands on
    atm = 0
    for p in s.parts
        p.kind in VESSEL_KINDS || continue
        n = Node(p.name, :vessel)
        g(k) = p.p[k]
        n.V, n.D = g("volume_L") * 1e-3, g("diameter_mm") * 1e-3
        n.H = n.V / (π * n.D^2 / 4)
        n.Cw = max(g("wall_kg") * g("wall_cp_J_kgK"), 1.0)
        n.A_in = π * n.D * n.H + π * n.D^2 / 2
        n.absorb, n.emiss, n.shade, n.heater_W = g("absorptivity"), g("emissivity"), g("shade"), g("heater_W")
        push!(nodes, n)
        for pin in KINDS[p.kind].pins
            height = p.kind == "Bottle" ? g("draw_height_frac") : pin.name == "bottom" ? 0.0 : pin.name == "vent" ? g("vent_height_frac") : 1.0
            push!(n.ports, height); push!(n.port_h, 0.0); push!(n.port_fl, 0.0); push!(n.port_cover, 0.0)
            k = net_of[(p.name, pin.name)]
            node_of_net[k], port_of_net[k] = length(nodes), length(n.ports)
        end
    end
    for p in s.parts
        p.kind == "Atmosphere" || continue
        if atm == 0
            push!(nodes, Node("atmosphere", :atm))
            atm = length(nodes)
            push!(nodes[atm].port_h, 0.0); push!(nodes[atm].port_fl, 0.0)
        end
        node_of_net[net_of[(p.name, "a")]] = atm
    end
    for p in s.parts, pin in KINDS[p.kind].pins          # remaining nets with a flow part: junctions
        k = net_of[(p.name, pin.name)]
        (p.kind in FLOW_KINDS && node_of_net[k] == 0) || continue
        n = Node("net$k", :junction)
        n.P = n.P_hold = env.P_atm
        n.hl = vapour_PT(env.P_atm, env.T_amb).h
        push!(nodes, n)
        node_of_net[k] = length(nodes)
    end
    pin_node = Dict(pin => node_of_net[k] for (pin, k) in net_of if node_of_net[k] != 0)

    edges = Edge[]
    for p in s.parts
        p.kind in FLOW_KINDS || continue
        g(k) = p.p[k]
        ka, kb = net_of[(p.name, "a")], net_of[(p.name, "b")]
        e = Edge(p.name, EDGE_KIND[p.kind], node_of_net[ka], port_of_net[ka], node_of_net[kb], port_of_net[kb],
                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.3, false, false, 0, 0.0,
                 NaN, NaN, Fluid(NaN, NaN, NaN, NaN, NaN, NaN, NaN, false), NaN, NaN, Inf)
        if p.kind == "Pipe"
            e.L, e.D, e.K = g("L_m"), g("id_mm") * 1e-3, g("K")
            e.relrough = g("rough_mm") * 1e-3 / e.D
            e.area = π * e.D^2 / 4
        else
            e.area = g("Cd") * π / 4 * (g("d_mm") * 1e-3)^2
            if p.kind == "Valve"
                g("Cv") > 0 && (e.area = g("Cv") * CV_AREA)
                e.normally_open, e.travel = g("normally_open") > 0.5, max(g("travel_s"), 0.0)
                e.opening = valve_target(e)
            elseif p.kind == "Regulator"
                e.set, e.band = g("set_psig") * PSI + env.P_atm, max(g("droop_psi"), 0.1) * PSI
            elseif p.kind == "CheckValve"
                e.set, e.band = g("crack_psi") * PSI, 0.2 * g("crack_psi") * PSI + 500.0
            elseif p.kind == "ReliefValve"
                e.set, e.band = g("set_psig") * PSI + env.P_atm, 0.1 * g("set_psig") * PSI + 500.0
            end
        end
        push!(edges, e)
    end

    i = 1
    for n in nodes
        n.kind === :vessel || continue
        n.i = i
        i += nstates(n)
    end
    junctions = [k for k in eachindex(nodes) if nodes[k].kind === :junction]
    row = zeros(Int, length(nodes))
    row[junctions] .= 1:length(junctions)
    adj = [[k for k in eachindex(edges) if edges[k].a == j || edges[k].b == j] for j in eachindex(nodes)]
    gauges = Dict(p.name => get(pin_node, (p.name, "a"), 0) for p in s.parts if p.kind == "Gauge")
    net = Net(nodes, edges, env, i - 1, junctions, adj, gauges, pin_node, collect(eachindex(nodes)),
              zeros(length(junctions)), zeros(length(junctions)), zeros(length(junctions)),
              zeros(length(junctions), length(junctions)), row, nothing, 0,
              falses(length(edges)) |> collect, Int[], Int[], zeros(Int, length(nodes)), zeros(Int, length(nodes)), Int[], UInt(1))

    u0 = zeros(i + 2)
    for n in nodes
        n.kind === :vessel || continue
        # Start every vessel as one mixed volume at T; switch_mode! then splits it into zones if
        # it has both liquid and ullage. A bottle filled to its limit on a hot day legitimately
        # starts liquid-full.
        p = s.parts[find_part(s, n.name)]
        T = p.p["T_C"] + 273.15
        air = get(p.p, "air", 0.0) > 0.5
        m = max(p.p["fill_kg"], (air ? 0.001 : 1.0) * vapour_PT(env.P_atm, T).ρ * n.V)   # never a vacuum
        ρ = m / n.V
        if T < SAT_TMAX && sat_T(T).ρv < ρ < sat_T(T).ρl
            st = sat_T(T)
            x = (1 / ρ - 1 / st.ρl) / (1 / st.ρv - 1 / st.ρl)
            u0[n.i+1] = m * (st.ul + x * (st.uv - st.ul))
        else
            u0[n.i+1] = m * eos(ρ, T).u
            air && ρ < N2O_RHOC && (u0[n.i+5] = max(env.P_atm - eos(ρ, T).P, 0.0) * n.V / (R_AIR * T))
        end
        u0[n.i], u0[n.i+4] = m, T
        switch_mode!(n, u0)
    end
    net, u0
end

# ---------------------------------------------------------------- flows

"Effective flow area of an edge, given the pressures at its a and b ends."
function flow_area(e::Edge, Pa, Pb)
    e.kind === :valve && return e.area * e.opening
    e.kind === :regulator && return e.area * clamp((e.set - Pb) / e.band, 0.0, 1.0)
    e.kind === :check && return e.area * clamp((Pa - Pb - e.set) / e.band, 0.0, 1.0)
    e.kind === :relief && return e.area * clamp((Pa - e.set) / e.band, 0.0, 1.0)
    e.area
end

"""
Mass flow a → b through an edge (negative for b → a), given the pressure at each end and the
enthalpy each end would send if it were upstream.
"""
function edge_flow(e::Edge, Pa, Pb, ha, hb)
    TALLY.edge_flow += 1
    forward = Pa >= Pb
    forward || e.kind in (:valve, :orifice, :pipe) || return 0.0     # the rest are one-way
    A = flow_area(e, Pa, Pb)
    A <= 0 && return 0.0
    P1, P2, h1 = forward ? (Pa, Pb, ha) : (Pb, Pa, hb)
    if P1 != e.m_P1 || h1 != e.m_h1                 # a new upstream state: resolve it once and remember
        e.m_P1, e.m_h1, e.m_up, e.m_Pc, e.m_Pu = P1, h1, state_Ph(P1, h1, e.m_up), NaN, Inf
    end
    up = e.m_up
    G = massflux(up, P1, P2, e)
    if e.kind === :pipe && G > 0
        # Friction as a discharge coefficient: ΔP = (1 + K + f·L/D)·G²/2ρ, so cd = (1 + K + f·L/D)^-½
        # with f depending on the Reynolds number cd·G·D/μ. Substituting back converges by only a
        # factor of ten a time for a long hose, so after two substitutions the secant rule takes
        # over: four evaluations of f reach 1e-12, where eight substitutions reached 1e-8.
        R0 = G * e.D / visc_mix(up)
        φ(c) = 1 / sqrt(1 + e.K + friction_factor(c * R0, e.relrough) * e.L / e.D)
        c0 = 1 / sqrt(1 + e.K + 0.02 * e.L / e.D)
        cd = φ(c0)
        F0 = cd - c0
        for _ in 1:3
            F1 = φ(cd) - cd
            F1 == F0 && break                           # already exact
            c0, cd = cd, cd - F1 * (cd - c0) / (F1 - F0)
            F0 = F1
        end
        G *= cd
    end
    forward ? A * G : -A * G
end

sent_h(n::Node, port) = n.kind === :junction ? n.hl : n.port_h[port]

"""
Decide which parts can carry flow right now and which junctions actually need solving.

A junction with only one open path is a dead end: nothing can flow into it, so its pressure is
simply that of whatever it hangs off, and the part leading to it carries exactly zero. Removing
dead ends (repeatedly, since removing one can expose another) leaves the junctions that sit on a
real flow path; only those go to Newton. With every valve shut there may be none at all. Besides
the saving, this makes dead-end flows exactly zero instead of a trickle of solver noise.

Sets `net.live` (per edge), `net.active` (junction node indices to solve), `net.row`, and
`net.pruned` / `net.feeder` (dead-end junctions in the order removed, and the edge each hangs by;
0 for a junction cut off on every side). Returns true if the plan differs from last time.
"""
function plan!(net::Net)
    nodes, edges = net.nodes, net.edges
    stamp = UInt(0)
    for (k, e) in enumerate(edges)
        net.live[k] = flow_area(e, nodes[e.a].P, nodes[e.b].P) > 0
        net.live[k] || (e.mdot = 0.0)                   # shut: no flow, and it will not be visited again
        stamp = hash(net.live[k], stamp)
    end
    changed = stamp != net.stamp
    net.stamp = stamp
    empty!(net.pruned); empty!(net.active)
    fill!(net.row, 0)
    deg, stack = net.deg, net.stack
    empty!(stack)
    for j in net.junctions
        deg[j] = count(k -> net.live[k], net.adj[j])
        deg[j] <= 1 && push!(stack, j)
    end
    while !isempty(stack)
        j = pop!(stack)
        deg[j] < 0 && continue                          # already removed
        k = findfirst(k -> net.live[k], net.adj[j])
        net.feeder[j] = k === nothing ? 0 : net.adj[j][k]
        deg[j] = -1
        push!(net.pruned, j)
        k === nothing && continue
        e = edges[net.adj[j][k]]
        net.live[net.adj[j][k]] = false
        e.mdot = 0.0
        other = e.a == j ? e.b : e.a
        if nodes[other].kind === :junction && deg[other] > 0
            deg[other] -= 1
            deg[other] <= 1 && push!(stack, other)
        end
    end
    for j in net.junctions
        deg[j] >= 0 || continue
        push!(net.active, j)
        net.row[j] = length(net.active)
    end
    changed
end

# Dead-end junctions take the pressure of what they hang off (in reverse order of removal, so that
# what they hang off has already been set). One cut off on every side keeps its last pressure.
function settle_dead_ends!(net::Net)
    for j in Iterators.reverse(net.pruned)
        k = net.feeder[j]
        if k == 0
            net.nodes[j].P = net.nodes[j].P_hold
        else
            e = net.edges[k]
            net.nodes[j].P = _rest_pressure(e, j, net.nodes[e.a == j ? e.b : e.a].P, net.nodes[j].P_hold)
        end
    end
    nothing
end

# Where a dead end stops taking or giving flow through `e`, whose far end is at `P`. Through a
# pipe, valve or orifice that is when the pressures are equal. A part that shuts on the dead end's
# own pressure stops earlier: a regulator locks up at its set-point (and, being one-way, never
# pulls its outlet down), a check valve stops one cracking pressure short, a relief valve reseats
# at its set pressure.
function _rest_pressure(e::Edge, j, P, P_hold)
    e.kind === :regulator && e.b == j && return max(P_hold, min(P, e.set))
    e.kind === :check && return e.b == j ? P - e.set : P + e.set
    e.kind === :relief && e.a == j && return max(P, e.set)
    P
end

"""
One sweep at the current pressures: the flow in every live edge, every active junction's mixture
enthalpy and its mass imbalance (left in `dm_l`). Nodes are visited from high pressure to low so
that a junction's enthalpy is known before anything flows out of it.
"""
function sweep!(net::Net)
    nodes = net.nodes
    sortperm!(net.order, nodes; by = n -> -n.P)
    for j in net.active
        n = nodes[j]
        n.dm_l = n.min_l = n.Hin_l = 0.0
    end
    for i in net.order
        n = nodes[i]
        n.kind === :junction && net.row[i] != 0 && n.min_l > 0 && (n.hl = n.Hin_l / n.min_l)
        for k in net.adj[i]
            net.live[k] || continue
            e = net.edges[k]
            other = e.a == i ? e.b : e.a
            o = nodes[other]
            # each edge is handled once, from its higher-pressure end
            (n.P > o.P || (n.P == o.P && i < other)) || continue
            m = edge_flow(e, nodes[e.a].P, nodes[e.b].P, sent_h(nodes[e.a], e.pa), sent_h(nodes[e.b], e.pb))
            e.mdot = m
            flow = abs(m)
            n.kind === :junction && (n.dm_l -= flow)
            if o.kind === :junction
                o.dm_l += flow
                o.min_l += flow
                o.Hin_l += flow * sent_h(n, e.a == i ? e.pa : e.pb)
            end
        end
    end
    nothing
end

# Sweep, then collect the active junctions' imbalances. Returns the worst one and the largest flow.
function _residual!(net::Net)
    sweep!(net)
    scale = worst = 0.0
    for (k, e) in enumerate(net.edges)
        net.live[k] && (scale = max(scale, abs(e.mdot)))
    end
    for (q, j) in enumerate(net.active)
        net.res[q] = net.nodes[j].dm_l
        worst = max(worst, abs(net.res[q]))
    end
    worst, scale
end

"""
Find the junction pressures at which every junction's flows balance.

Dead ends are removed first (`plan!`). The rest is Newton's method from the last solution, with
the Jacobian from differencing each edge with respect to its two end pressures. Flow goes like
√ΔP, on which a full Newton step overshoots by exactly its own length and cycles for ever, so a
step is halved if the imbalance does not fall. Parts whose opening depends on pressure (relief
and check valves, regulators) can open or shut as the pressures move, so the plan is checked
again afterwards and the solve repeated if it changed.
"""
function solve_junctions!(net::Net)
    for round in 1:4
        plan!(net) && (net.lu = nothing)                # a different set of unknowns: old factorisation is useless
        _solve_active!(net)
        settle_dead_ends!(net)
        # Did the solve open anything that was planned as shut? Then its flow was left out: go
        # again. (The other way round needs nothing: a part that the solved pressures shut, such
        # as a regulator reaching lock-up, was given its flow at those pressures, which is zero.)
        stable = true
        for (k, e) in enumerate(net.edges)
            now_open = flow_area(e, net.nodes[e.a].P, net.nodes[e.b].P) > 0
            was_open = net.live[k] || e.mdot == 0 && _hangs_dead(net, k)
            now_open && !was_open && (stable = false)
        end
        stable && return
    end
    net.unconverged += 1                                # still opening and shutting parts after four rounds
    nothing
end

# Is edge k one that plan! switched off as the feeder of a dead end (rather than because it is shut)?
_hangs_dead(net::Net, k) = any(j -> net.feeder[j] == k, net.pruned)

function _solve_active!(net::Net)
    na = length(net.active)
    na == 0 && return sweep!(net)
    nodes, r, P0 = net.nodes, view(net.res, 1:na), net.p0
    J = view(net.jac, 1:na, 1:na)
    # Flow runs downhill, so no junction can end up outside the range spanned by the vessels, the
    # atmosphere and the pressures the junctions are holding.
    Plo, Phi = Inf, -Inf
    for n in nodes
        p = n.kind === :junction ? n.P_hold : n.P
        Plo, Phi = min(Plo, p), max(Phi, p)
    end
    for j in net.active
        nodes[j].P = clamp(nodes[j].P, Plo, Phi)
    end
    worst, scale = _residual!(net)
    # First reuse the Jacobian factorised last time. Successive evaluations (the ODE solver's
    # stages, its finite-difference columns) differ only slightly, so this usually converges in
    # two or three sweeps without rebuilding the Jacobian, which is the expensive part.
    if net.lu !== nothing
        for _ in 1:8
            converged = worst < 1e-9 * scale + 1e-10
            r .= .-r
            ldiv!(net.lu, r)
            all(isfinite, r) || break
            for (q, j) in enumerate(net.active)
                nodes[j].P = clamp(nodes[j].P + r[q], Plo, Phi)
            end
            before = worst
            worst, scale = _residual!(net)
            converged && return                        # that was the extra step past the tolerance
            worst > 0.3 * before && break              # not contracting: the Jacobian is stale
        end
    end
    for iter in 1:40
        converged = worst < 1e-9 * scale + 1e-10       # kg/s; 1e-10 is a third of a milligram an hour
        TALLY.newton += 1
        fill!(J, 0.0)
        for (k, e) in enumerate(net.edges)
            net.live[k] || continue
            ra, rb = net.row[e.a], net.row[e.b]
            (ra == 0 && rb == 0) && continue
            Pa, Pb = nodes[e.a].P, nodes[e.b].P
            ha, hb = sent_h(nodes[e.a], e.pa), sent_h(nodes[e.b], e.pb)
            δ = 1e-6 * max(Pa, Pb) + 1.0
            # Only junction ends are unknowns: an end at a vessel or the atmosphere needs no
            # derivative. The downstream end goes first, because moving it leaves the upstream
            # state (and with it the remembered flash and choking point) untouched.
            m0 = edge_flow(e, Pa, Pb, ha, hb)
            first_b = Pa >= Pb
            da = db = 0.0
            for end_b in (first_b, !first_b)
                if end_b
                    rb != 0 && (db = (edge_flow(e, Pa, Pb + δ, ha, hb) - m0) / δ)
                else
                    ra != 0 && (da = (edge_flow(e, Pa + δ, Pb, ha, hb) - m0) / δ)
                end
            end
            ra != 0 && (J[ra, ra] -= da)
            rb != 0 && (J[rb, rb] += db)
            ra != 0 && rb != 0 && (J[ra, rb] -= db; J[rb, ra] += da)
        end
        r .= .-r
        net.lu = lu!(J; check = false)
        ldiv!(net.lu, r)                               # r now holds the Newton step
        big = maximum(abs, r)
        isfinite(big) || break
        big > 2e6 && (r .*= 2e6 / big)                 # no pressure moves more than 20 bar at once
        # Pressures that have stopped moving are converged, whatever the imbalance says: the flow
        # laws carry round-off of their own, and near it the tolerance above cannot be met.
        converged |= big < 1e-2
        for (q, j) in enumerate(net.active)
            P0[q] = nodes[j].P
            net.step[q] = r[q]
        end
        for λ in (1.0, 0.5, 0.25)
            for (q, j) in enumerate(net.active)
                nodes[j].P = clamp(P0[q] + λ * net.step[q], Plo, Phi)
            end
            w2, scale = _residual!(net)
            # Once converged, take the one extra full step regardless: it brings the result down
            # to round-off, which matters because the ODE solver differentiates through this solve.
            (converged || w2 < worst || λ == 0.25) && (worst = w2; break)
        end
        converged && return
    end
    loose = 1e-6 * scale + 1e-10                       # one part in a million of the largest flow
    worst < loose && return                            # as good as the flow laws' round-off allows
    # Newton did not settle (a kink in a flow law: a check valve cracking, a restriction starting
    # to choke next to a closing valve). Fall back to something that cannot fail: balance each
    # junction in turn with its neighbours held, and repeat. Each junction's imbalance falls
    # monotonically as its own pressure rises, so its root can always be bracketed and found.
    for pass in 1:60
        moved = 0.0
        for j in net.active
            before = nodes[j].P
            _balance_one!(net, j, Plo, Phi)
            moved = max(moved, abs(nodes[j].P - before))
        end
        worst, scale = _residual!(net)                 # also brings the junction enthalpies up to date
        (worst < 1e-9 * scale + 1e-10 || moved < 1e-3 || pass >= 3 && worst < loose) && return
    end
    net.unconverged += 1
    nothing
end

# Mass flowing into junction j minus mass flowing out, at the present pressures. As in `sweep!`,
# what flows out carries the mixture enthalpy of what flows in, so the inflows are taken first.
function _imbalance(net::Net, j)
    nodes = net.nodes
    n = nodes[j]
    flow_in = H = flow_out = 0.0
    for inward in (true, false)
        for k in net.adj[j]
            net.live[k] || continue
            e = net.edges[k]
            other = e.a == j ? e.b : e.a
            o = nodes[other]
            (o.P > n.P || (o.P == n.P && other < j)) == inward || continue
            m = abs(edge_flow(e, nodes[e.a].P, nodes[e.b].P, sent_h(nodes[e.a], e.pa), sent_h(nodes[e.b], e.pb)))
            inward ? (flow_in += m; H += m * sent_h(o, e.a == other ? e.pa : e.pb)) : (flow_out += m)
        end
        inward && flow_in > 0 && (n.hl = H / flow_in)
    end
    flow_in - flow_out
end

# Set junction j to the pressure at which it balances, its neighbours held. Only its own parts are
# evaluated. Walk from where it is in the direction that reduces the imbalance, with a growing
# stride, until the imbalance changes sign; then close in by false position (Illinois variant),
# with a plain bisection every third step so that a kinked flow law cannot slow it down.
function _balance_one!(net::Net, j, Plo, Phi)
    n = net.nodes[j]
    a = n.P = clamp(n.P, Plo, Phi)
    fa = _imbalance(net, j)
    (fa == 0 || !isfinite(fa)) && return
    b, fb, stride = a, fa, 1e3
    while (fb > 0) == (fa > 0)
        (b >= Phi && fb > 0 || b <= Plo && fb < 0) && return        # balanced nowhere in range: rest at the limit
        a, fa = b, fb
        b = n.P = clamp(b + (fb > 0 ? stride : -stride), Plo, Phi)
        stride *= 4
        fb = _imbalance(net, j)
        (fb == 0 || !isfinite(fb)) && return
    end
    side = 0
    for it in 1:60
        abs(b - a) < 0.01 && break
        c = it % 3 == 0 ? (a + b) / 2 : (a * fb - b * fa) / (fb - fa)
        min(a, b) < c < max(a, b) || (c = (a + b) / 2)
        n.P = c
        fc = _imbalance(net, j)
        fc == 0 && return
        if (fc > 0) == (fa > 0)
            a, fa = c, fc
            side == 1 && (fb /= 2)
            side = 1
        else
            b, fb = c, fc
            side == 2 && (fa /= 2)
            side = 2
        end
    end
    n.P = abs(fa) < abs(fb) ? a : b
    nothing
end

"The right-hand side of the ODE system. The last three states total vented mass, vented enthalpy and heat taken in."
function rhs!(du, u, net::Net, t)
    TALLY.rhs += 1
    nodes, env = net.nodes, net.env
    for n in nodes
        evaluate!(n, u, env)
    end
    solve_junctions!(net)
    for n in nodes
        n.kind === :vessel && (n.dm_l = n.min_l = n.Hin_l = 0.0)
    end
    vent_m = vent_h = 0.0
    for e in net.edges
        m = e.mdot
        m == 0 && continue
        up, kup, dn, kdn = m > 0 ? (nodes[e.a], e.pa, nodes[e.b], e.pb) : (nodes[e.b], e.pb, nodes[e.a], e.pa)
        h = sent_h(up, kup)
        outflow!(up, kup, abs(m), u)
        inflow!(dn, kdn, abs(m), h)
        dn.kind === :atm && (vent_m += abs(m); vent_h += abs(m) * h)
        up.kind === :atm && (vent_m -= abs(m); vent_h -= abs(m) * h)
    end
    q = 0.0
    for n in nodes
        q += derivative!(du, n, u, env)
    end
    du[end-2], du[end-1], du[end] = vent_m, vent_h, q
    nothing
end

"""
Jacobian of `rhs!` by forward differences with a relative step (the library default uses an
absolute floor that is far too coarse for small states). The bookkeeping states do not feed back.
"""
function jac!(J, u, net::Net, t)
    TALLY.jac += 1
    n = length(u)
    f0, f1 = zeros(n), zeros(n)
    rhs!(f0, u, net, t)
    fill!(J, 0.0)
    for j in 1:net.nfluid
        uj = u[j]
        δ = 1e-7 * abs(uj) + 1e-11
        u[j] = uj + δ
        rhs!(f1, u, net, t)
        u[j] = uj
        @views J[:, j] .= (f1 .- f0) ./ δ
    end
    nothing
end

# ---------------------------------------------------------------- running

"""
A running simulation: the network, the integrator, the remaining schedule and the recorded history.

`t` and `u` are "now": the time and state the caller has advanced to. The integrator is allowed to
run ahead of now (it takes whatever step size the physics allows, up to seconds when nothing is
happening) and `u` is read off its interpolant. Forcing a solver step to end at every screen frame
instead cost 60 steps per simulated second and made the window lag. Anything that changes the
equations (a valve command, a heater, a scheduled action) first winds the integrator back to now.
"""
mutable struct Sim
    net::Net
    integ::Any
    t::Float64
    u::Vector{Float64}
    schedule::Vector{Action}
    next_action::Int
    log::Vector{Action}                       # every valve/heater action applied, scripted or clicked
    m0::Float64                               # N2O mass and energy at t = 0, for the closure check
    E0::Float64
    t_hist::Vector{Float64}
    hist::Dict{String,Vector{Float64}}
    sample_dt::Float64
    last_sample::Float64
    failed::String                            # non-empty once the solver has given up
end

function Sim(scene::Scene; sample_dt = 0.1, reltol = 1e-5)
    net, u0 = compile(scene)
    f = ODEFunction(rhs!; jac = jac!, tgrad = (dT, u, p, t) -> fill!(dT, 0.0))
    prob = ODEProblem(f, u0, (0.0, Inf), net)
    integ = init(prob, Rodas5P(); reltol, abstol = 1e-8, save_everystep = false, save_start = false,
                 save_end = false, maxiters = typemax(Int), dtmax = 5.0)
    m0, E0 = totals(net, u0)
    sim = Sim(net, integ, 0.0, copy(u0), sort(scene.schedule; by = a -> a.t), 1, Action[], m0, E0, Float64[], Dict{String,Vector{Float64}}(), sample_dt, -Inf, "")
    refresh!(sim)
    sample!(sim)
    sim
end

simtime(sim::Sim) = sim.t

"Evaluate flows and junction pressures at the present state, for display and recording."
function refresh!(sim::Sim)
    net, u = sim.net, sim.u
    for n in net.nodes
        evaluate!(n, u, net.env)
    end
    solve_junctions!(net)
    for j in net.junctions
        net.nodes[j].P_hold = net.nodes[j].P
    end
    nothing
end

function totals(net::Net, u)
    m = E = 0.0
    for n in net.nodes
        mi, Ei = inventory(n, u)
        m += mi; E += Ei
    end
    m + u[end-2], E + u[end-1] - u[end]
end

"""
Relative mass and energy imbalance since t = 0: what is in the vessels plus what was vented, minus
heat taken in, compared with the start. Both should stay near zero; the energy figure is scaled by
the starting mass times 200 kJ/kg.
"""
function closure(sim::Sim)
    m, E = totals(sim.net, sim.u)
    (m - sim.m0) / max(sim.m0, 1e-3), (E - sim.E0) / (max(sim.m0, 1e-3) * 2e5)
end

"Restart the integrator from now. Needed whenever the equations change under it."
_restart!(sim::Sim) = reinit!(sim.integ, copy(sim.u); t0 = sim.t, tf = Inf, erase_sol = true, reset_dt = true)

"Wind the integrator back to now if it has run ahead."
_rewind!(sim::Sim) = sim.integ.t > sim.t && _restart!(sim)

"Apply an action (`open`, `close`, `toggle`, `deenergise`, `stuck_open`, `stuck_closed`, `unstick`, `on`, `off`, `set`) to a part by name."
function act!(sim::Sim, part, action, value = 0.0; record = true)
    net = sim.net
    ie = findfirst(e -> e.name == part, net.edges)
    in_ = findfirst(n -> n.name == part, net.nodes)
    _rewind!(sim)
    if ie !== nothing
        e = net.edges[ie]
        if e.kind === :valve
            action == "open" && (e.cmd = !e.normally_open)
            action == "close" && (e.cmd = e.normally_open)
            action == "toggle" && (e.cmd = !e.cmd)
            action == "deenergise" && (e.cmd = false)
            action == "stuck_open" && (e.stuck = 1)
            action == "stuck_closed" && (e.stuck = 2)
            action == "unstick" && (e.stuck = 0)
        elseif e.kind in (:regulator, :relief) && action == "set"
            e.set = value * PSI + net.env.P_atm
            _restart!(sim)
        else
            return false
        end
    elseif in_ !== nothing && action in ("on", "off")
        net.nodes[in_].heater_on = action == "on"
        _restart!(sim)
    else
        return false
    end
    record && push!(sim.log, Action(simtime(sim), part, action, value))
    true
end

"Drop every valve to its de-energised (fail) state, as on loss of power."
fail_safe!(sim::Sim) = foreach(e -> e.kind === :valve && act!(sim, e.name, "deenergise"), sim.net.edges)

_valves_moving(sim::Sim) = any(e -> e.kind === :valve && e.opening != valve_target(e), sim.net.edges)

# Move valves towards their targets by `dt` seconds of travel.
function _stroke!(sim::Sim, dt)
    for e in sim.net.edges
        e.kind === :valve || continue
        target = valve_target(e)
        step = e.travel > 0 ? dt / e.travel : 1.0
        e.opening = target > e.opening ? min(e.opening + step, target) : max(e.opening - step, target)
        # Ten steps of 0.1 do not sum to exactly 1, and a valve left 1e-16 open is "open" to the
        # junction solve with a flow too small to pin any pressure down. Arrive exactly.
        abs(e.opening - target) < 1e-9 && (e.opening = target)
    end
end

# Make `t` now: read the state off the integrator (which is at or past t), refresh, record if due.
function _arrive!(sim::Sim, t)
    integ = sim.integ
    t >= integ.t ? copyto!(sim.u, integ.u) : integ(sim.u, t)
    sim.t = min(t, integ.t)
    refresh!(sim)
    sim.t - sim.last_sample >= sim.sample_dt - 1e-9 && sample!(sim)
    nothing
end

"""
    advance!(sim, dt; max_steps = 2000) -> seconds actually advanced

Move now forward by up to `dt` seconds: run scheduled actions when they fall due, stroke valves,
step the equations as far as needed, switch vessel modes, and record samples. It returns early if
it uses `max_steps` solver steps (so a caller drawing frames never hangs), and if the solver fails
`sim.failed` says why and the simulation stops advancing.
"""
function advance!(sim::Sim, dt; max_steps = 2000)
    isempty(sim.failed) || return 0.0
    integ = sim.integ
    t_start = sim.t
    target = t_start + dt
    steps = 0
    while sim.t < target - 1e-9
        # Things that happen now: scheduled actions falling due, valves part-way through a stroke.
        due = sim.next_action <= length(sim.schedule) && sim.schedule[sim.next_action].t <= sim.t + 1e-9
        while sim.next_action <= length(sim.schedule) && sim.schedule[sim.next_action].t <= sim.t + 1e-9
            a = sim.schedule[sim.next_action]
            act!(sim, a.part, a.action, a.value)
            sim.next_action += 1
        end
        moving = _valves_moving(sim)
        if moving                                # valve areas change in 20 ms slices
            _stroke!(sim, 0.02)
            _restart!(sim)
        elseif due
            _rewind!(sim)
        end
        # The integrator must stop exactly where the equations next change; otherwise it is free.
        event = moving ? sim.t + 0.02 : Inf
        sim.next_action <= length(sim.schedule) && (event = min(event, sim.schedule[sim.next_action].t))
        t_stop = min(target, max(event, sim.t + 1e-6))
        isfinite(event) && event > integ.t && add_tstop!(integ, max(event, integ.t + 1e-9))
        while integ.t < t_stop - 1e-12
            try
                step!(integ)
            catch err
                sim.failed = "solver error at t = $(round(integ.t, digits = 3)) s: " * first(sprint(showerror, err), 300)
                return sim.t - t_start
            end
            steps += 1
            if !successful_retcode(integ.sol.retcode) && integ.sol.retcode != ReturnCode.Default
                sim.failed = "solver stopped at t = $(round(integ.t, digits = 3)) s: $(integ.sol.retcode)"
                return sim.t - t_start
            end
            switched = false
            for n in sim.net.nodes
                switched |= switch_mode!(n, integ.u)
            end
            if switched                          # the state was rewritten: now jumps to this step's end
                _arrive!(sim, integ.t)
                _restart!(sim)
                t_stop = min(t_stop, sim.t)
                break
            end
            integ.t <= t_stop + 1e-12 && _arrive!(sim, integ.t)
            if steps >= max_steps
                _arrive!(sim, min(integ.t, t_stop))
                return sim.t - t_start
            end
        end
        sim.t < t_stop && _arrive!(sim, t_stop)
    end
    sim.t - t_start
end

"Run until simulation time `t`, or until the solver fails or `wall_limit` seconds have passed."
function run_until!(sim::Sim, t; wall_limit = Inf)
    start = Base.time()
    while simtime(sim) < t - 1e-9 && isempty(sim.failed) && Base.time() - start < wall_limit
        advance!(sim, t - simtime(sim))
    end
    sim
end

# ---------------------------------------------------------------- reading results

"""
Current values for one node, in SI. Vessels report both zones; a junction reports the fluid
passing through it. Call after `refresh!` (which `advance!` does).
"""
function node_report(net::Net, n::Node, u)
    env = net.env
    if n.kind === :atm
        return (; P = env.P_atm, T = env.T_amb, m = 0.0, level = 0.0, liquid_frac = 0.0, Tv = env.T_amb, Tw = env.T_amb, Psat = NaN, x = 1.0)
    elseif n.kind === :junction
        f = state_Ph(n.P, n.hl)
        return (; P = n.P, T = f.T, m = 0.0, level = 0.0, liquid_frac = 0.0, Tv = f.T, Tw = f.T, Psat = NaN, x = f.x)
    elseif n.strat
        s = strat(n, u)
        return (; P = s.P, T = s.Tl, m = s.ml + s.mv, level = s.level, liquid_frac = s.Vl / n.V, Tv = s.Tv, Tw = u[n.i+4],
                Psat = Psat(s.Tl), x = s.mv / (s.ml + s.mv))
    end
    f = homog(n, u)
    lf = f.twophase ? (1 - f.x) * f.ρ / sat_T(f.T).ρl : f.x >= 1 ? 0.0 : 1.0
    (; P = f.P + air_pressure(n, u, f), T = f.T, m = u[n.i], level = lf * n.H, liquid_frac = lf, Tv = f.T, Tw = u[n.i+4],
     Psat = f.T < SAT_TMAX ? Psat(f.T) : NaN, x = f.x)
end

function sample!(sim::Sim)
    net, u = sim.net, sim.u
    t = simtime(sim)
    push!(sim.t_hist, t)
    put(k, v) = push!(get!(sim.hist, k, Float64[]), v)
    for n in net.nodes
        n.kind === :atm && continue
        r = node_report(net, n, u)
        put(n.name * ".P", r.P); put(n.name * ".T", r.T)
        if n.kind === :vessel
            put(n.name * ".m", r.m); put(n.name * ".level", r.level); put(n.name * ".liquid_frac", r.liquid_frac)
            put(n.name * ".liquid_vol_pct", 100 * r.liquid_frac); put(n.name * ".liquid_mass_pct", 100 * (1 - r.x))
            put(n.name * ".Tv", r.Tv); put(n.name * ".Tw", r.Tw); put(n.name * ".Psat", r.Psat)
        end
    end
    for e in net.edges
        put(e.name * ".mdot", e.mdot)
        # Recorded by meaning, not as a bare 0/1: position is how far open the valve actually is,
        # open_cmd is what it is being told to do, stuck says why the two may differ.
        if e.kind === :valve
            put(e.name * ".position", e.opening)
            put(e.name * ".open_cmd", e.cmd != e.normally_open ? 1.0 : 0.0)
            put(e.name * ".stuck", Float64(e.stuck))
        elseif e.kind in (:regulator, :relief)
            put(e.name * ".set", e.set)
        end
    end
    for n in net.nodes
        n.kind === :vessel && n.heater_W > 0 && put(n.name * ".heater_on", n.heater_on ? 1.0 : 0.0)
    end
    put("ambient.T", net.env.T_amb); put("ambient.solar", net.env.solar); put("ambient.wind", net.env.wind)
    for (g, k) in net.gauges
        k == 0 && continue
        r = node_report(net, net.nodes[k], u)
        put(g * ".P", r.P); put(g * ".T", r.T)
    end
    put("vented.m", u[end-2])
    sim.last_sample = t
    nothing
end

"Warm-soak check: pressure and liquid volume fraction if N2O mass m sat closed in volume V at T."
function warm_soak(m, V, T)
    ρ = m / V
    if T < SAT_TMAX
        st = sat_T(T)
        st.ρv < ρ < st.ρl && return st.P, (ρ - st.ρv) / (st.ρl - st.ρv)
    end
    eos(ρ, T).P, ρ > N2O_RHOC ? 1.0 : 0.0
end

"""
Write the timestamped action log (every valve, heater and set-point action applied, scripted or
clicked) as CSV: `t_s,part,action,value`. With the scene file this is enough to repeat the run.
"""
function write_events(path, sim::Sim)
    open(path, "w") do io
        println(io, "t_s,part,action,value")
        for a in sim.log
            println(io, a.t, ",", a.part, ",", a.action, ",", a.value)
        end
    end
end

"Read an action log written by `write_events` back as a schedule."
function read_events(path)
    out = Action[]
    for line in Iterators.drop(eachline(path), 1)
        f = split(line, ',')
        length(f) == 4 && push!(out, Action(parse(Float64, f[1]), String(f[2]), String(f[3]), parse(Float64, f[4])))
    end
    sort!(out; by = a -> a.t)
end

"Path of the events file that goes with a signal CSV: `run.csv` → `run_events.csv`."
events_path(csv) = splitext(csv)[1] * "_events.csv"

"Write the recorded history as CSV (SI units, one column per signal)."
function write_csv(path, sim::Sim)
    keys_ = sort!(collect(keys(sim.hist)))
    open(path, "w") do io
        println(io, "t_s,", join(keys_, ","))
        for i in eachindex(sim.t_hist)
            print(io, sim.t_hist[i])
            for k in keys_
                print(io, ",", sim.hist[k][i])
            end
            println(io)
        end
    end
end
