module ChebyshevCubicScan

using Kneading.Diagrams

export chebyshev_cubic,
    chebyshev_critical_points,
    compute_chebyshev_cubic_diagram

chebyshev_critical_points() = (-0.5, 0.5)

@inline function chebyshev_cubic(u::Real, v::Real, x::Real)
    return ((u - v) / 2) * (4 * x^3 - 3 * x) + (u + v) / 2
end

function _advance_critical_orbits!(
    orbit_values::Array{Float64,3},
    plane::ParameterPlane,
)
    for x_index in eachindex(plane.x)
        u = plane.x[x_index]
        for y_index in eachindex(plane.y)
            v = plane.y[y_index]
            orbit_values[y_index, x_index, 1] =
                chebyshev_cubic(u, v, orbit_values[y_index, x_index, 1])
            orbit_values[y_index, x_index, 2] =
                chebyshev_cubic(u, v, orbit_values[y_index, x_index, 2])
        end
    end
    return orbit_values
end

function compute_chebyshev_cubic_diagram(;
    grid_size::Integer = 1000,
    iterates::Integer = 20,
    u_range::Tuple{<:Real,<:Real} = (-2.0, 2.0),
    v_range::Tuple{<:Real,<:Real} = (-2.0, 2.0),
)
    grid_size >= 2 || throw(ArgumentError("grid_size must be at least two"))
    iterates >= 2 || throw(ArgumentError("iterates must be at least two"))
    u_range[1] < u_range[2] || throw(ArgumentError("the u range must be ordered"))
    v_range[1] < v_range[2] || throw(ArgumentError("the v range must be ordered"))

    u_values = collect(range(Float64(u_range[1]), Float64(u_range[2]); length = grid_size))
    v_values = collect(range(Float64(v_range[1]), Float64(v_range[2]); length = grid_size))
    plane = ParameterPlane(u_values, v_values; xname = "𝑢", yname = "𝑣")
    critical_points = chebyshev_critical_points()
    diagram = KneadingDiagram(
        plane;
        metadata = (
            family = :chebyshev_cubic,
            critical_points = critical_points,
            iterates = Int(iterates),
        ),
    )

    orbit_values = Array{Float64}(undef, length(v_values), length(u_values), 2)
    orbit_values[:, :, 1] .= critical_points[1]
    orbit_values[:, :, 2] .= critical_points[2]

    for iterate in 2:Int(iterates)
        _advance_critical_orbits!(orbit_values, plane)
        for (source_index, source) in enumerate((:left_critical, :right_critical))
            field = view(orbit_values, :, :, source_index)
            for level in critical_points
                add_level_contours!(
                    diagram,
                    field;
                    source,
                    iterate,
                    level,
                )
            end
        end
    end
    return diagram
end

end
