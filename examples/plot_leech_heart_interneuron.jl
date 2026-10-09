using CairoMakie

function read_table(path)
    lines = readlines(path)
    header = split(first(lines), '\t')
    rows = [Dict(zip(header, split(line, '\t'; keepempty = true))) for line in lines[2:end] if !isempty(line)]
    return rows
end

function grid(rows, key; seed = nothing)
    selected = isnothing(seed) ? rows : filter(row -> row["seed"] == seed, rows)
    shifts = sort!(unique(parse(Float64, row["shift"]) for row in rows))
    currents = sort!(unique(parse(Float64, row["current"]) for row in rows))
    values = fill(NaN, length(shifts), length(currents))
    xs = Dict(x => j for (j, x) in enumerate(shifts))
    ys = Dict(y => i for (i, y) in enumerate(currents))
    for row in selected
        values[xs[parse(Float64, row["shift"])], ys[parse(Float64, row["current"])]] = key(row)
    end
    return 1000 .* shifts, 1000 .* currents, values
end

function tail_period(word; tail = 48, longest = 24)
    length(word) < tail && return NaN
    w = word[end-tail+1:end]
    for p in 1:longest
        all(w[k] == w[k-p] for k in p+1:tail) && return p
    end
    return NaN
end

word_hash(word) = isempty(word) ? NaN : Float64(mod(hash(word), 4096))

function equilibrium_curves!(axis, equilibria)
    count_x, count_y, counts = grid(equilibria, row -> length(split(row["types"], ',')))
    _, _, upper = grid(equilibria, row -> parse(Int, first(split(row["types"], ',')[end])))
    contour!(axis, count_x, count_y, counts; levels = [2.0], color = :black, linewidth = 2.5)
    contour!(axis, count_x, count_y, upper; levels = [1.0], color = :black, linewidth = 2.5, linestyle = :dash)
    return axis
end

function panel(figure, position, title, x, y, values, equilibria; colormap, colorrange = nothing, label = "", categorical = false)
    axis = Axis(figure[position...]; title, xlabel = "V_K2 shift (mV)", ylabel = "I_app (pA)",
        xgridvisible = false, ygridvisible = false)
    range_kw = isnothing(colorrange) ? (;) : (; colorrange)
    plot = heatmap!(axis, x, y, values; colormap, nan_color = :white, interpolate = false, range_kw...)
    equilibrium_curves!(axis, equilibria)
    limits!(axis, first(x), last(x), first(y), last(y))
    categorical || Colorbar(figure[position[1], position[2] + 1], plot; label)
    return axis
end

function plot_leech(directory; output = joinpath(directory, "figures"))
    orbits = read_table(joinpath(directory, "orbits.tsv"))
    equilibria = read_table(joinpath(directory, "equilibria.tsv"))
    mkpath(output)
    focus_word = row -> row["status"] == "complete" ? word_hash(row["transition_word"]) : NaN
    x, y, words = grid(orbits, focus_word; seed = "saddle_focus")
    _, _, periods = grid(orbits, row -> row["status"] == "complete" ? tail_period(row["transition_word"]) : NaN; seed = "saddle_focus")
    _, _, focus_return = grid(orbits, row -> log10(parse(Float64, row["longest_return"])); seed = "saddle_focus")
    _, _, saddle_events = grid(orbits, row -> parse(Float64, row["events"]); seed = "real_saddle")
    _, _, saddle_return = grid(orbits, row -> log10(parse(Float64, row["longest_return"])); seed = "real_saddle")
    complete = count(row -> row["seed"] == "saddle_focus" && row["status"] == "complete", orbits)
    focus_points = count(row -> row["seed"] == "saddle_focus", orbits)
    palette = [RGBf(0.15 + 0.7 * mod(0.618k, 1), 0.2 + 0.6 * mod(0.381k + 0.3, 1), 0.25 + 0.65 * mod(0.7k + 0.6, 1)) for k in 0:4095]
    paths = String[]

    figure = Figure(size = (1500, 1250), fontsize = 18, backgroundcolor = :white)
    Label(figure[0, 1:4], "Leech heart interneuron: flow kneading of m_K2 minima, $(length(x)) × $(length(y))";
        fontsize = 24, font = :bold, tellwidth = false)
    panel(figure, (1, 1), "Saddle-focus W^u: transition words ($complete/$focus_points complete)", x, y, words, equilibria;
        colormap = cgrad(palette; categorical = true), colorrange = (-0.5, 4095.5), categorical = true)
    panel(figure, (1, 3), "Saddle-focus W^u: period of word tail", x, y, periods, equilibria;
        colormap = cgrad(:turbo, 12; categorical = true), colorrange = (0.5, 12.5), label = "events per period")
    panel(figure, (2, 1), "Saddle-focus W^u: log10 longest return (s)", x, y, focus_return, equilibria;
        colormap = :magma, label = "log10 s")
    panel(figure, (2, 3), "Real saddle W^u: events before rest", x, y, saddle_events, equilibria;
        colormap = cgrad(:viridis, 16; categorical = true), colorrange = (-0.5, 15.5), label = "m_K2 minima")
    Label(figure[3, 1:4], "solid: fold of equilibria (three equilibria below)   ·   dashed: Andronov–Hopf of the depolarized equilibrium   ·   white: no word";
        fontsize = 16, tellwidth = false, color = :gray35)
    path = joinpath(output, "leech-heart-interneuron-kneading.png")
    save(path, figure; px_per_unit = 1)
    push!(paths, abspath(path))

    figure = Figure(size = (1500, 650), fontsize = 18, backgroundcolor = :white)
    panel(figure, (1, 1), "Real saddle W^u: log10 longest return (s)", x, y, saddle_return, equilibria;
        colormap = :magma, label = "log10 s")
    panel(figure, (1, 3), "Saddle-focus W^u: transition words", x, y, words, equilibria;
        colormap = cgrad(palette; categorical = true), colorrange = (-0.5, 4095.5), categorical = true)
    path = joinpath(output, "leech-heart-interneuron-saddle.png")
    save(path, figure; px_per_unit = 1)
    push!(paths, abspath(path))
    foreach(p -> println("Saved ", p), paths)
    return paths
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    plot_leech(isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "leech-heart-interneuron") : ARGS[1])
