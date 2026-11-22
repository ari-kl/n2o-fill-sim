# The plots. One plot per kind of quantity that has a signal ticked, sharing one time axis.

const PLOT_KINDS = ((:pressure, "Pressure"), (:temperature, "Temperature"), (:flow, "Mass flow"), (:mass, "Mass"),
                    (:share, "Liquid share"), (:state, "Valves and heaters"))
const MAX_POINTS = 4000       # per line; longer histories are drawn every k-th sample

# Draw every k-th sample straight from the arrays (stride is in bytes), without copying them.
strided(k) = (base = ImPlot.ImPlotSpec(); ImPlot.ImPlotSpec(ntuple(i -> fieldname(ImPlot.ImPlotSpec, i) === :Stride ? Int32(8k) : getfield(base, i), fieldcount(ImPlot.ImPlotSpec))...))

function plot_panel(app::App)
    sim = app.sim
    n = length(sim.t_hist)
    u = app.units
    tag = (u.P, u.T, u.m, u.flow)                 # the cache holds display units
    tag == app.cache_units || (empty!(app.cache); app.cache_units = tag)
    active = Tuple{Symbol,String,Vector{String}}[]
    for (kind, title) in PLOT_KINDS
        keys_ = sort!([k for k in app.probes if haskey(sim.hist, k) && signal_kind(k) === kind])
        has_overlay = any(o -> o[2] === kind, app.overlay)
        (isempty(keys_) && !has_overlay) || push!(active, (kind, title, keys_))
    end
    if isempty(active)
        text("Nothing to plot. Tick signals in the Signals tab.")
        return
    end
    t_now = simtime(sim)
    if app.follow
        app.xlo[], app.xhi[] = 0.0, max(t_now, 10.0) * 1.02
    end
    np = length(active)
    rows = np <= 3 ? 1 : 2
    cols = cld(np, rows)
    avail = ig.GetContentRegionAvail()
    w = (avail.x - px(8) * cols) / cols
    h = (avail.y - px(6) * rows) / rows
    stride = max(1, cld(n, MAX_POINTS))
    spec = strided(stride)
    cursor = NaN
    k_cursor = isfinite(app.cursor_t) ? clamp(searchsortedfirst(sim.t_hist, app.cursor_t), 1, n) : 0
    for (i, (kind, title, keys_)) in enumerate(active)
        label = isempty(unit(app, kind)) ? title : title * " [" * unit(app, kind) * "]"
        fit = kind === :state ? ImPlot.ImPlotAxisFlags_None : ImPlot.ImPlotAxisFlags_AutoFit | ImPlot.ImPlotAxisFlags_RangeFit
        if ImPlot.BeginPlot(label * "###" * String(kind), "t [s]", "", ig.ImVec2(w, h); y_flags = fit)
            ImPlot.SetupAxisLinks(ImPlot.ImAxis_X1, app.xlo, app.xhi)
            # 0 = shut or off, 1 = open or on: a little room either side, or the traces sit on the frame
            kind === :state && ImPlot.SetupAxisLimits(ImPlot.ImAxis_Y1, -0.1, 1.1, ImPlot.ImPlotCond_Always)
            ImPlot.SetupLegend(ImPlot.ImPlotLocation_SouthEast)
            for key in keys_
                c = get!(app.cache, key, Float64[])
                length(c) > n && empty!(c)
                for j in length(c)+1:n
                    push!(c, shown(app, kind, sim.hist[key][j]))
                end
                # With the mouse over a plot, every legend entry shows its value at that time.
                name = k_cursor == 0 ? key : key * ":  " * string(round(c[k_cursor], sigdigits = 4))
                if kind === :state
                    st = get!(() -> Steps(Float64[], Float64[], 0), app.steps, key)
                    st.seen > n && (empty!(st.t); empty!(st.v); st.seen = 0)
                    for j in st.seen+1:n
                        (isempty(st.v) || c[j] != st.v[end]) && (push!(st.t, sim.t_hist[j]); push!(st.v, c[j]))
                    end
                    st.seen = n
                    push!(st.t, sim.t_hist[n]); push!(st.v, c[n])      # hold the last value up to now
                    ImPlot.PlotStairs(name * "###" * key, st.t, st.v, length(st.t))
                    pop!(st.t); pop!(st.v)
                else
                    ImPlot.PlotLine(name * "###" * key, sim.t_hist, c, cld(n, stride), spec)
                end
            end
            last = searchsortedlast(app.overlay_t, t_now)       # logged data, up to now
            first_ = searchsortedfirst(app.overlay_t, 0.0)
            if last > first_
                for (name, okind, ys) in app.overlay
                    okind === kind || continue
                    c = get!(() -> [shown(app, kind, y) for y in ys], app.cache, "log " * name)    # converted once per unit choice
                    ImPlot.PlotLine("log " * name, view(app.overlay_t, first_:last), view(c, first_:last), last - first_ + 1)
                end
            end
            isfinite(app.cursor_t) && ImPlot.PlotInfLines("##cursor", [app.cursor_t])
            ImPlot.IsPlotHovered() && (cursor = ImPlot.GetPlotMousePos().x)
            ImPlot.EndPlot()
        end
        i % cols == 0 || i == np || ig.SameLine()
    end
    app.cursor_t = cursor
    nothing
end
