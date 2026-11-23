# The plots. One plot per kind of quantity that has a signal ticked, sharing one time axis.

const PLOT_KINDS = ((:pressure, "Pressure"), (:temperature, "Temperature"), (:flow, "Mass flow"), (:mass, "Mass"),
                    (:share, "Liquid share"), (:state, "Valves and heaters"))
const MAX_POINTS = 4000       # per line; longer histories are drawn every k-th sample

# The default plot style with some fields changed. `Stride` is in bytes: 8k draws every k-th
# sample straight from the arrays, without copying them.
function plot_spec(; changes...)
    base = ImPlot.ImPlotSpec()
    ImPlot.ImPlotSpec(ntuple(fieldcount(ImPlot.ImPlotSpec)) do i
        f = fieldname(ImPlot.ImPlotSpec, i)
        haskey(changes, f) ? convert(fieldtype(ImPlot.ImPlotSpec, i), changes[f]) : getfield(base, i)
    end...)
end
strided(k) = plot_spec(Stride = 8k)

# A legend entry carries the series' value: at the time under the mouse, or now. The number is
# padded to a fixed width with figure spaces (as wide as a digit), because the legend stands beside
# the plot and takes its width from the longest entry: unpadded, the plot would change width with
# every digit, under the mouse that is choosing the time.
legend_label(key, x) = key * ":  " * lpad(string(round(x, sigdigits = 4)), 8, '\u2007')

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
    value_w, log_w = ig.CalcTextSize(legend_label("", 0.0)).x, ig.CalcTextSize("log ").x
    k_cursor = isfinite(app.cursor_t) ? clamp(searchsortedfirst(sim.t_hist, app.cursor_t), 1, n) : 0
    for (i, (kind, title, keys_)) in enumerate(active)
        label = isempty(unit(app, kind)) ? title : title * " [" * unit(app, kind) * "]"
        fit = kind === :state ? ImPlot.ImPlotAxisFlags_None : ImPlot.ImPlotAxisFlags_AutoFit | ImPlot.ImPlotAxisFlags_RangeFit
        if ImPlot.BeginPlot(label * "###" * String(kind), "t [s]", "", ig.ImVec2(w, h); y_flags = fit)
            ImPlot.SetupAxisLinks(ImPlot.ImAxis_X1, app.xlo, app.xhi)
            # 0 = shut or off, 1 = open or on: a little room either side, or the traces sit on the frame
            kind === :state && ImPlot.SetupAxisLimits(ImPlot.ImAxis_Y1, -0.1, 1.1, ImPlot.ImPlotCond_Always)
            # The legend stands beside the plot, not on it: inside, it sat on the newest data. Only
            # when that would leave almost no plot (a narrow window, long names) does it go back
            # inside, in the corner with the oldest data.
            legend_w = 0f0
            for key in keys_
                legend_w = max(legend_w, ig.CalcTextSize(key).x + value_w)
            end
            for (name, okind, _) in app.overlay
                okind === kind && (legend_w = max(legend_w, ig.CalcTextSize(name).x + log_w))
            end
            beside = w - legend_w - px(34) >= px(150)             # 34: the colour patch and the legend's padding
            ImPlot.SetupLegend(beside ? ImPlot.ImPlotLocation_NorthEast : ImPlot.ImPlotLocation_NorthWest,
                               beside ? ImPlot.ImPlotLegendFlags_Outside : ImPlot.ImPlotLegendFlags_None)
            for key in keys_
                c = get!(app.cache, key, Float64[])
                length(c) > n && empty!(c)
                for j in length(c)+1:n
                    push!(c, shown(app, kind, sim.hist[key][j]))
                end
                # Inside the plot the value is shown only while the mouse picks a time, to cover less.
                name = beside || k_cursor != 0 ? legend_label(key, c[k_cursor == 0 ? n : k_cursor]) : key
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
                # Measured lines carry dots, so that they differ from simulated ones by more than
                # hue: 16 to 32 of them however many samples are drawn. The spacing is a power of
                # two, so that as a replay adds samples the dots stay on the samples they were on
                # and every other one drops out at a doubling, rather than all of them creeping.
                count = last - first_ + 1
                every = prevpow(2, max(1, count ÷ 16))
                dots = plot_spec(Stride = 8every, Marker = ImPlot.ImPlotMarker_Circle, MarkerSize = px(2.5), MarkerFillColor = ig.ImVec4(0.08f0, 0.08f0, 0.1f0, 1f0))
                for (name, okind, ys) in app.overlay
                    okind === kind || continue
                    c = get!(() -> [shown(app, kind, y) for y in ys], app.cache, "log " * name)    # converted once per unit choice
                    ImPlot.PlotLine("log " * name, view(app.overlay_t, first_:last), view(c, first_:last), count)
                    ImPlot.PlotScatter("log " * name, view(app.overlay_t, first_:last), view(c, first_:last), cld(count, every), dots)
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
