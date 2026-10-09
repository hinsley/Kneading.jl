using CairoMakie

function read_table(path)
    lines = readlines(path)
    header = split(first(lines), '\t')
    return [Dict(zip(header, split(line, '\t'; keepempty = true))) for line in lines[2:end] if !isempty(line)]
end

function grid(rows, value; seed = nothing)
    shifts = sort!(unique(parse(Float64, row["shift"]) for row in rows))
    currents = sort!(unique(parse(Float64, row["current"]) for row in rows))
    values = fill(NaN, length(shifts), length(currents))
    xs = Dict(x => j for (j, x) in enumerate(shifts))
    ys = Dict(y => i for (i, y) in enumerate(currents))
    for row in rows
        isnothing(seed) || row["seed"] == seed || continue
        values[xs[parse(Float64, row["shift"])], ys[parse(Float64, row["current"])]] = value(row)
    end
    return 1000 .* shifts, 1000 .* currents, values
end

resting(row) = row["status"] == "numerical_failure"

function word_color(row)
    row["status"] == "complete" || resting(row) || return NaN
    return Float64(mod(hash((row["transition_word"], resting(row))), 4096))
end

function tail_period(row; tail = 40, longest = 16)
    word = row["transition_word"]
    row["status"] == "complete" && length(word) >= tail || return NaN
    w = word[end-tail+1:end]
    for p in 1:longest
        all(w[k] == w[k-p] for k in p+1:tail) && return p
    end
    return NaN
end

function events_before_rest(row; slow = 2.0)
    resting(row) || return NaN
    events = parse(Float64, row["events"])
    gap = events == 1 ? parse(Float64, row["last_event_time"]) : parse(Float64, row["longest_return"])
    return events > 0 && gap > slow ? events - 1 : events
end

longest_return(row) = log10(parse(Float64, row["longest_return"]))

const word_palette = [RGBf(0.15 + 0.7 * mod(0.618k, 1), 0.2 + 0.6 * mod(0.381k + 0.3, 1),
    0.25 + 0.65 * mod(0.7k + 0.6, 1)) for k in 0:4095]

function equilibrium_curves!(axis, equilibria; linewidth = 2.5)
    x, y, counts = grid(equilibria, row -> length(split(row["types"], ',')))
    _, _, upper = grid(equilibria, row -> parse(Int, first(split(row["types"], ',')[end])))
    contour!(axis, x, y, counts; levels = [2.0], color = :black, linewidth)
    contour!(axis, x, y, upper; levels = [1.0], color = :black, linewidth, linestyle = :dash)
    return axis
end

function homoclinic_curves!(axis, orbits; color = :red, linewidth = 1.6)
    x, y, steps = grid(orbits, events_before_rest; seed = "real_saddle")
    finite = filter(isfinite, steps)
    isempty(finite) && return axis
    interior = copy(steps)
    for (column, original) in zip(eachrow(interior), eachrow(steps)), k in 1:length(column)
        any(isnan, original[k:min(k + 2, end)]) && (column[k] = NaN)
    end
    levels = collect(0.5:1:min(maximum(finite), 14))
    contour!(axis, x, y, interior; levels, color, linewidth)
    return axis
end

function heatmap_axis(position, title, x, y, values, equilibria, orbits; colormap, colorrange = nothing,
    label = nothing, homoclinics = true, kw...)
    axis = Axis(position[1, 1]; title, xlabel = "V_K2 shift (mV)", ylabel = "I_app (pA)",
        xgridvisible = false, ygridvisible = false, kw...)
    range_kw = isnothing(colorrange) ? (;) : (; colorrange)
    plot = heatmap!(axis, x, y, values; colormap, nan_color = :white, interpolate = false, range_kw...)
    equilibrium_curves!(axis, equilibria)
    homoclinics && homoclinic_curves!(axis, orbits)
    limits!(axis, first(x), last(x), first(y), last(y))
    isnothing(label) || Colorbar(position[1, 2], plot; label)
    return axis
end

function plot_leech(directory; output = joinpath(directory, "figures"))
    orbits = read_table(joinpath(directory, "orbits.tsv"))
    equilibria = read_table(joinpath(directory, "equilibria.tsv"))
    mkpath(output)
    x, y, words = grid(orbits, word_color; seed = "saddle_focus")
    _, _, periods = grid(orbits, tail_period; seed = "saddle_focus")
    _, _, focus_return = grid(orbits, longest_return; seed = "saddle_focus")
    _, _, saddle_events = grid(orbits, events_before_rest; seed = "real_saddle")
    _, _, saddle_return = grid(orbits, longest_return; seed = "real_saddle")
    word_map = cgrad(word_palette; categorical = true)
    paths = String[]

    figure = Figure(size = (1500, 1300), fontsize = 18, backgroundcolor = :white)
    Label(figure[0, 1:2], "Leech heart interneuron: flow kneading at m_K2 minima, $(length(x)) × $(length(y)) points";
        fontsize = 24, font = :bold, tellwidth = false)
    heatmap_axis(figure[1, 1], "Saddle-focus unstable manifold: transition words", x, y, words, equilibria, orbits;
        colormap = word_map, colorrange = (-0.5, 4095.5))
    heatmap_axis(figure[1, 2], "Saddle-focus unstable manifold: period of the word tail", x, y, periods, equilibria, orbits;
        colormap = cgrad(:turbo, 16; categorical = true), colorrange = (0.5, 16.5), label = "events per period")
    heatmap_axis(figure[2, 1], "Saddle-focus unstable manifold: log10 longest return (s)", x, y, focus_return,
        equilibria, orbits; colormap = :magma, label = "log10 s")
    heatmap_axis(figure[2, 2], "Real-saddle separatrix: events before rest", x, y, saddle_events, equilibria, orbits;
        colormap = cgrad(:viridis, 21; categorical = true), colorrange = (-0.5, 20.5), label = "m_K2 minima")
    Label(figure[3, 1:2], "solid: fold of equilibria (three equilibria below)  ·  dashed: Andronov–Hopf  ·  " *
        "red: homoclinic orbits of the middle saddle (steps of the real-saddle staircase)  ·  white: no word";
        fontsize = 16, tellwidth = false, color = :gray30)
    path = joinpath(output, "leech-heart-interneuron-kneading.png")
    save(path, figure; px_per_unit = 1)
    push!(paths, abspath(path))

    figure = Figure(size = (1500, 650), fontsize = 18, backgroundcolor = :white)
    heatmap_axis(figure[1, 1], "Real-saddle separatrix: events before rest", x, y, saddle_events, equilibria, orbits;
        colormap = cgrad(:viridis, 21; categorical = true), colorrange = (-0.5, 20.5), label = "m_K2 minima")
    heatmap_axis(figure[1, 2], "Real-saddle separatrix: log10 longest return (s)", x, y, saddle_return, equilibria,
        orbits; colormap = :magma, label = "log10 s")
    path = joinpath(output, "leech-heart-interneuron-saddle.png")
    save(path, figure; px_per_unit = 1)
    push!(paths, abspath(path))

    figure = Figure(size = (1600, 1650), fontsize = 30, backgroundcolor = :white)
    axis = heatmap_axis(figure[1, 1:5], "", x, y, words, equilibria, orbits; colormap = word_map,
        colorrange = (-0.5, 4095.5), xlabelsize = 34, ylabelsize = 34, homoclinics = false)
    homoclinic_curves!(axis, orbits; linewidth = 2.5)
    for (text, position) in slide_labels
        text!(axis, position...; text, color = :white, strokecolor = :black, strokewidth = 2,
            fontsize = 30, font = :bold, align = (:center, :center))
    end
    traces = isfile(joinpath(directory, "traces.tsv")) ? read_table(joinpath(directory, "traces.tsv")) : []
    labels = unique(row["label"] for row in traces)
    for (k, label) in enumerate(labels)
        rows = filter(row -> row["label"] == label, traces)
        shift, current = 1000 .* parse.(Float64, (rows[1]["shift"], rows[1]["current"]))
        scatter!(axis, [shift], [current]; color = :white, strokecolor = :black, strokewidth = 2, markersize = 26)
        text!(axis, shift, current; text = label, fontsize = 20, font = :bold, align = (:center, :center))
        inset = Axis(figure[2, k]; title = label, titlesize = 26, xlabel = "t (s)", ylabel = k == 1 ? "V (mV)" : "",
            xgridvisible = false, ygridvisible = false, yticklabelsvisible = k == 1, xlabelsize = 24, ylabelsize = 24,
            xticklabelsize = 20, yticklabelsize = 20)
        lines!(inset, parse.(Float64, getindex.(rows, "time")), 1000 .* parse.(Float64, getindex.(rows, "voltage"));
            color = :black, linewidth = 1.2)
        ylims!(inset, -55, 45)
    end
    rowsize!(figure.layout, 1, Relative(0.72))
    Label(figure[3, 1:5], "colors: transition words of the critical orbit on the saddle-focus unstable manifold  ·  " *
        "white: no word\nsolid: fold  ·  dashed: Andronov–Hopf  ·  red: homoclinic orbits to the middle saddle";
        fontsize = 22, tellwidth = false, color = :gray30)
    path = joinpath(output, "leech-heart-interneuron-slide.png")
    save(path, figure; px_per_unit = 1)
    push!(paths, abspath(path))
    foreach(p -> println("Saved ", p), paths)
    return paths
end

const slide_labels = [
    ("small\noscillations", (-25.8, 4.0)),
    ("bursting:\nspike adding", (-19.0, 8.0)),
    ("tonic\nspiking", (-10.0, 2.0)),
    ("Andronov–Hopf", (-21.0, 24.0)),
    ("stable equilibrium", (-24.5, 31.0)),
    ("fold", (-26.5, -8.0)),
    ("rest after N spikes:\nhomoclinics to the saddle", (-19.0, -24.0)),
]

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    plot_leech(isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "leech-heart-interneuron") : ARGS[1])
