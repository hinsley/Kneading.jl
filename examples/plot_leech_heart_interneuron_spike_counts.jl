include(joinpath(@__DIR__, "plot_leech_heart_interneuron.jl"))

counts_of(row) = isempty(row["counts"]) ? Int[] : parse.(Int, split(row["counts"], ','))

count_word_color(row) = row["status"] == "complete" ? Float64(mod(hash(row["counts"]), 4096)) : NaN

function count_outcome(row)
    row["status"] == "complete" && return NaN
    row["status"] == "rest" && return 1.0
    row["status"] == "no_spikes" && return 2.0
    return 3.0
end

function final_count(row; tail = 6)
    counts = counts_of(row)
    row["status"] == "complete" && length(counts) >= tail || return NaN
    window = counts[end-tail+1:end]
    return all(==(first(window)), window) ? Float64(first(window)) : 0.0
end

first_count(row) = isempty(row["counts"]) ? NaN : Float64(first(counts_of(row)))

const outcome_colors = cgrad([RGBf(0.75, 0.75, 0.75), RGBf(0.35, 0.35, 0.35), RGBf(1.0, 0.8, 0.8)]; categorical = true)
const count_colors = cgrad(vcat([RGBf(0, 0, 0)], [cgrad(:turbo)[k / 15] for k in 1:15]); categorical = true)

function plot_spike_counts(directory; output = joinpath(directory, "figures"))
    rows = [merge(row, Dict("seed" => "saddle_focus")) for row in read_table(joinpath(directory, "spike-counts.tsv"))]
    orbits = read_table(joinpath(directory, "orbits.tsv"))
    equilibria = read_table(joinpath(directory, "equilibria.tsv"))
    mkpath(output)
    x, y, words = grid(rows, count_word_color)
    _, _, outcomes = grid(rows, count_outcome)
    _, _, final = grid(rows, final_count)
    _, _, initial = grid(rows, first_count)
    figure = Figure(size = (1500, 1300), fontsize = 18, backgroundcolor = :white)
    Label(figure[0, 1:2], "Leech heart interneuron: spikes per burst along the saddle-focus critical orbit, " *
        "$(length(x)) × $(length(y)) points"; fontsize = 22, font = :bold, tellwidth = false)
    axis = heatmap_axis(figure[1, 1], "Spike-count words (first 14 bursts)", x, y, words, equilibria, orbits;
        colormap = cgrad(word_palette; categorical = true), colorrange = (-0.5, 4095.5))
    heatmap!(axis, x, y, outcomes; colormap = outcome_colors, colorrange = (0.5, 3.5), nan_color = :transparent)
    equilibrium_curves!(axis, equilibria)
    homoclinic_curves!(axis, orbits)
    heatmap_axis(figure[1, 2], "Spikes per burst at the end (0: not periodic)", x, y, final, equilibria, orbits;
        colormap = count_colors, colorrange = (-0.5, 15.5), label = "spikes per burst")
    heatmap_axis(figure[2, 1], "Spikes in the first burst", x, y, initial, equilibria, orbits;
        colormap = cgrad(:viridis, 40; categorical = true), colorrange = (0.5, 40.5), label = "spikes")
    _, _, orientation = grid(orbits, word_color; seed = "saddle_focus")
    heatmap_axis(figure[2, 2], "Orientation words, same orbits", x, y, orientation, equilibria, orbits;
        colormap = cgrad(word_palette; categorical = true), colorrange = (-0.5, 4095.5))
    Label(figure[3, 1:2], "light gray: ends at rest  ·  dark gray: no spikes above 0 mV  ·  pink: unfinished after 60 s  ·  " *
        "white: no seed\nsolid: fold  ·  dashed: Andronov–Hopf  ·  red: homoclinic orbits to the middle saddle";
        fontsize = 16, tellwidth = false, color = :gray30)
    path = joinpath(output, "leech-heart-interneuron-spike-counts.png")
    save(path, figure; px_per_unit = 1)
    println("Saved ", abspath(path))
    return abspath(path)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    plot_spike_counts(isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "leech-heart-interneuron") : ARGS[1])
