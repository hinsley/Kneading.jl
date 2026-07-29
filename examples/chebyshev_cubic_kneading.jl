package_root = dirname(@__DIR__)
package_root in LOAD_PATH || pushfirst!(LOAD_PATH, package_root)

using CairoMakie
using Kneading.Diagrams

@inline function chebyshev_cubic(u, v, x)
    return ((u - v) / 2) * (4x^3 - 3x) + (u + v) / 2
end

function main()
    grid_size = parse(Int, get(ENV, "CHEBYSHEV_GRID_SIZE", "1000"))
    iterates = parse(Int, get(ENV, "CHEBYSHEV_ITERATES", "20"))
    output_path = get(
        ENV,
        "CHEBYSHEV_OUTPUT",
        joinpath(@__DIR__, "chebyshev_cubic_kneading_diagram.png"),
    )

    parameter_values = collect(range(-2.0, 2.0; length = grid_size))
    plane = ParameterPlane(
        parameter_values,
        parameter_values;
        xname = "𝑢",
        yname = "𝑣",
    )
    critical_points = (-0.5, 0.5)
    diagram = KneadingDiagram(
        plane;
        metadata = (
            family = :chebyshev_cubic,
            critical_points,
            iterates,
        ),
    )
    orbit_values = Array{Float64}(undef, grid_size, grid_size, 2)
    orbit_values[:, :, 1] .= critical_points[1]
    orbit_values[:, :, 2] .= critical_points[2]

    for iterate in 2:iterates
        scan_plane!(orbit_values, plane) do values, u, v
            values[1] = chebyshev_cubic(u, v, values[1])
            values[2] = chebyshev_cubic(u, v, values[2])
        end
        for (source_index, source) in enumerate((:left_critical, :right_critical))
            for level in critical_points
                add_level_contours!(
                    diagram,
                    view(orbit_values, :, :, source_index);
                    source,
                    iterate,
                    level,
                )
            end
        end
    end

    save_kneading_contours(
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

    println("Saved $(length(diagram.layers)) contour layers to $output_path")
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
