using CairoMakie

function read_rossler_scan(path)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("the scan file is empty"))
    header = split(first(lines), '\t')
    columns = Dict(name => i for (i, name) in enumerate(header))
    required = ("c", "a", "status", "raw_code", "transition_code", "raw_length", "transition_length")
    all(name -> haskey(columns, name), required) ||
        throw(ArgumentError("the scan must contain c, a, status, codes, and word lengths"))
    rows = [split(line, '\t'; keepempty = true) for line in Iterators.drop(lines, 1) if !isempty(line)]
    isempty(rows) && throw(ArgumentError("the scan has no parameter points"))
    c_values = sort!(unique(parse(Float64, row[columns["c"]]) for row in rows))
    a_values = sort!(unique(parse(Float64, row[columns["a"]]) for row in rows))
    length(rows) == length(c_values) * length(a_values) ||
        throw(ArgumentError("the scan must contain a complete rectangular parameter grid"))
    c_indices = Dict(c => j for (j, c) in enumerate(c_values))
    a_indices = Dict(a => i for (i, a) in enumerate(a_values))
    raw = fill(NaN, length(a_values), length(c_values))
    transition = copy(raw)
    seen = falses(size(raw))
    statuses = Dict{String,Int}()
    for row in rows
        i = a_indices[parse(Float64, row[columns["a"]])]
        j = c_indices[parse(Float64, row[columns["c"]])]
        seen[i, j] && throw(ArgumentError("the scan contains a duplicated parameter point"))
        seen[i, j] = true
        status = row[columns["status"]]
        statuses[status] = get(statuses, status, 0) + 1
        if status == "complete"
            parse(Int, row[columns["raw_length"]]) == 8 &&
                parse(Int, row[columns["transition_length"]]) == 7 ||
                throw(ArgumentError("this Rössler plot requires 8 raw and 7 transition signs"))
            raw_code = parse(Int, row[columns["raw_code"]])
            transition_code = parse(Int, row[columns["transition_code"]])
            0 <= raw_code <= 255 && 0 <= transition_code <= 127 ||
                throw(ArgumentError("a complete word code is outside its encoding range"))
            raw[i, j] = raw_code
            transition[i, j] = transition_code
        end
    end
    return (; c_values, a_values, raw, transition, statuses)
end

function rossler_word_color(code)
    hue = mod(137code, 256) / 256
    saturation = 0.55 + 0.28 * ((code & 3) / 3)
    value = 0.76 + 0.18 * (((code >> 2) & 3) / 3)
    chroma = value * saturation
    x = chroma * (1 - abs(mod(6hue, 2) - 1))
    sector = floor(Int, 6hue)
    channels = ((chroma, x, 0.0), (x, chroma, 0.0), (0.0, chroma, x),
        (0.0, x, chroma), (x, 0.0, chroma), (chroma, 0.0, x))[sector + 1]
    rgb = round.(Int, 255 .* (channels .+ value .- chroma)) ./ 255
    return RGBf(rgb...)
end

function plot_rossler_scan(path; output_directory = joinpath(dirname(path), "figures"))
    data = read_rossler_scan(path)
    mkpath(output_directory)
    paths = String[]
    total = length(data.raw)
    completed = get(data.statuses, "complete", 0)
    grid = "$(length(data.c_values)) × $(length(data.a_values))"
    for (kind, bits, values, description) in (
        ("raw", 8, data.raw, "1 = positive tangent component; 0 = negative"),
        ("transition", 7, data.transition, "1 = preserved sign; 0 = reversed sign"),
    )
        last_code = 2^bits - 1
        palette = [rossler_word_color(kind == "raw" ? code : 128 | code) for code in 0:last_code]
        figure = Figure(size = (1400, 1050), fontsize = 20, backgroundcolor = :white)
        Label(figure[1, 1:2], "Rössler y-minima · $bits $kind signs";
            fontsize = 28, font = :bold, tellwidth = false)
        Label(figure[2, 1:2], "b = 0.3   ·   $grid parameter points   ·   initial critical event included";
            fontsize = 18, tellwidth = false, color = :gray35)
        axis = Axis(figure[3, 1]; xlabel = "c", ylabel = "a", xgridvisible = false,
            ygridvisible = false, xticks = 2:1:7, yticks = 0.30:0.05:0.55)
        plot = heatmap!(axis, data.c_values, data.a_values, permutedims(values);
            colormap = cgrad(palette; categorical = true),
            colorrange = (-0.5, last_code + 0.5), nan_color = :white, interpolate = false)
        limits!(axis, first(data.c_values), last(data.c_values), first(data.a_values), last(data.a_values))
        ticks = unique(vcat(collect(0:32:last_code), [last_code]))
        Colorbar(figure[3, 2], plot; label = "Binary word encoded as an integer", ticks)
        Label(figure[4, 1:2], description; fontsize = 18, tellwidth = false)
        Label(figure[5, 1:2], "$completed / $total complete words · white = incomplete or failed";
            fontsize = 17, tellwidth = false, color = :gray35)
        filename = joinpath(output_directory, "rossler-flow-kneading-$kind-$bits.png")
        save(filename, figure; px_per_unit = 1)
        push!(paths, abspath(filename))
        println("Saved ", abspath(filename))
    end
    return paths
end

function main()
    input = isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "rossler-flow-kneading.tsv") : ARGS[1]
    destination = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(input), "figures")
    return plot_rossler_scan(input; output_directory = destination)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
