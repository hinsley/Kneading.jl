package_root = dirname(@__DIR__)
package_root in LOAD_PATH || pushfirst!(LOAD_PATH, package_root)

using CairoMakie
using Kneading.Diagrams

include("chebyshev_cubic_scan.jl")
using .ChebyshevCubicScan

function environment_integer(name::String, default::Int)
    return parse(Int, get(ENV, name, string(default)))
end

function main()
    grid_size = environment_integer("CHEBYSHEV_GRID_SIZE", 1000)
    iterates = environment_integer("CHEBYSHEV_ITERATES", 20)
    output_path = get(
        ENV,
        "CHEBYSHEV_OUTPUT",
        joinpath(@__DIR__, "chebyshev_cubic_kneading_diagram.png"),
    )

    compute_seconds = @elapsed diagram = compute_chebyshev_cubic_diagram(
        grid_size = grid_size,
        iterates = iterates,
    )
    render_seconds = @elapsed save_kneading_contours(
        output_path,
        diagram;
        colors = Dict(
            :left_critical => RGBf(0.85, 0.15, 0.12),
            :right_critical => RGBf(0.10, 0.30, 0.90),
        ),
        opacity = iterate -> Float64(iterate - 1)^(-1.2),
        xticks = -2.0:0.5:2.0,
        yticks = -2.0:0.5:2.0,
    )

    println("grid_size=$grid_size")
    println("iterates=$iterates")
    println("layers=$(length(diagram.layers))")
    println("compute_seconds=$(round(compute_seconds; digits = 6))")
    println("render_seconds=$(round(render_seconds; digits = 6))")
    println("output=$output_path")
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
