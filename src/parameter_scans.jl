"""
    ParameterPlane(x, y; xname="x", yname="y")

A rectangular parameter grid. Scalar fields on the grid use Julia matrix
indexing: the first dimension selects rows corresponding to `y` values, and
the second dimension selects columns corresponding to `x` values, as in
`field[y_index, x_index]`.
"""
struct ParameterPlane{TX<:Real,TY<:Real}
    x::Vector{TX}
    y::Vector{TY}
    xname::String
    yname::String

    function ParameterPlane(
        x::AbstractVector{TX},
        y::AbstractVector{TY};
        xname::AbstractString = "x",
        yname::AbstractString = "y",
    ) where {TX<:Real,TY<:Real}
        length(x) >= 2 || throw(ArgumentError("the x axis needs at least two points"))
        length(y) >= 2 || throw(ArgumentError("the y axis needs at least two points"))
        all(isfinite, x) || throw(ArgumentError("the x axis must be finite"))
        all(isfinite, y) || throw(ArgumentError("the y axis must be finite"))
        all(x[index] < x[index + 1] for index in 1:(length(x) - 1)) ||
            throw(ArgumentError("the x axis must be strictly increasing"))
        all(y[index] < y[index + 1] for index in 1:(length(y) - 1)) ||
            throw(ArgumentError("the y axis must be strictly increasing"))
        return new{TX,TY}(collect(x), collect(y), String(xname), String(yname))
    end
end

"""
A line segment on a parameter plane.
"""
struct ContourSegment{T<:AbstractFloat}
    x1::T
    y1::T
    x2::T
    y2::T
end

"""
A set of contour segments from one source, iterate, and scalar level.
"""
struct ContourLayer{T<:AbstractFloat}
    source::Symbol
    iterate::Int
    level::T
    segments::Vector{ContourSegment{T}}
end

"""
A parameter plane, its contour layers, and scan-specific metadata.
"""
struct KneadingDiagram{P,M}
    plane::P
    layers::Vector{ContourLayer{Float64}}
    metadata::M
end

function KneadingDiagram(plane::ParameterPlane; metadata = NamedTuple())
    return KneadingDiagram(
        plane,
        ContourLayer{Float64}[],
        metadata,
    )
end

"""
    scan_plane!(sampler!, destination, plane)

Evaluate `sampler!` at every point of `plane`. The sampler receives a view of
`destination[y_index, x_index, :]` followed by the `x` and `y` parameter
values. `destination` is a three-dimensional output array whose first
dimension is `y` rows, second is `x` columns, and third stores the fields at
each parameter pair. Sampler calls can run concurrently and in any order.
`sampler!` must be thread-safe and must not mutate shared state outside its
provided destination slice.
"""
function scan_plane!(
    sampler!,
    destination::AbstractArray{T,3},
    plane::ParameterPlane,
) where {T}
    size(destination, 1) == length(plane.y) ||
        throw(DimensionMismatch("the first destination dimension must match the y axis"))
    size(destination, 2) == length(plane.x) ||
        throw(DimensionMismatch("the second destination dimension must match the x axis"))

    Threads.@threads for x_index in eachindex(plane.x) # Column
        x = plane.x[x_index]
        for y_index in eachindex(plane.y) # Row
            y = plane.y[y_index]
            sampler!(view(destination, y_index, x_index, :), x, y)
        end
    end
    return destination
end

"""
    scan_continuation!(
        step!,
        state_factory,
        destination,
        plane;
        along=:columns,
        direction=:forward,
    )

Scan each row or column as an independent continuation line. `along` must be
`:rows` or `:columns`, and defaults to `:columns`. Column continuation keeps
`x` fixed and traverses `y`. Row continuation keeps `y` fixed and traverses
`x`. `direction` must be `:forward` for increasing parameter indices or
`:backward` for decreasing parameter indices.

For each line, `state_factory(line_index, initial_x, initial_y)` is called once
at the first point in the selected direction. `line_index` is the `x` index
for a column or the `y` index for a row. The factory must return a separate
mutable state for that line. At each point, `step!` is called as
`step!(destination_slice, state, x, y)`. It may mutate the destination slice
and the line state.

Lines may run concurrently and in any order. Each line stays on one threaded
loop iteration, and its points run sequentially in the selected direction.
Both callbacks must be thread-safe and must not mutate shared state outside
the provided destination slice and line state. Destination slices use
`destination[y_index, x_index, :]`.
"""
function scan_continuation!(
    step!,
    state_factory,
    destination::AbstractArray{T,3},
    plane::ParameterPlane;
    along = :columns,
    direction = :forward,
) where {T}
    along in (:rows, :columns) ||
        throw(ArgumentError("along must be :rows or :columns"))
    direction in (:forward, :backward) ||
        throw(ArgumentError("direction must be :forward or :backward"))
    size(destination, 1) == length(plane.y) ||
        throw(DimensionMismatch("the first destination dimension must match the y axis"))
    size(destination, 2) == length(plane.x) ||
        throw(DimensionMismatch("the second destination dimension must match the x axis"))

    x_indices = direction === :forward ?
        (firstindex(plane.x):lastindex(plane.x)) :
        (lastindex(plane.x):-1:firstindex(plane.x))
    y_indices = direction === :forward ?
        (firstindex(plane.y):lastindex(plane.y)) :
        (lastindex(plane.y):-1:firstindex(plane.y))

    if along === :rows
        initial_x_index = first(x_indices)
        Threads.@threads for y_index in eachindex(plane.y)
            y = plane.y[y_index]
            state = state_factory(
                y_index,
                plane.x[initial_x_index],
                y,
            )
            for x_index in x_indices
                x = plane.x[x_index]
                step!(
                    view(destination, y_index, x_index, :),
                    state,
                    x,
                    y,
                )
            end
        end
    else
        initial_y_index = first(y_indices)
        Threads.@threads for x_index in eachindex(plane.x)
            x = plane.x[x_index]
            state = state_factory(
                x_index,
                x,
                plane.y[initial_y_index],
            )
            for y_index in y_indices
                y = plane.y[y_index]
                step!(
                    view(destination, y_index, x_index, :),
                    state,
                    x,
                    y,
                )
            end
        end
    end
    return destination
end

"""
Examine one edge of a marching-squares cell. Return the linearly interpolated
contour zero-crossing point when there is one; otherwise, return `nothing`.
"""
function _edge_point(
    first_value::Float64,
    second_value::Float64,
    first_point::NTuple{2,Float64},
    second_point::NTuple{2,Float64},
)
    if first_value == 0.0 && second_value == 0.0
        return nothing
    elseif first_value == 0.0
        return first_point
    elseif second_value == 0.0
        return second_point
    elseif signbit(first_value) == signbit(second_value)
        return nothing
    end

    fraction = first_value / (first_value - second_value)
    x = first_point[1] + fraction * (second_point[1] - first_point[1])
    y = first_point[2] + fraction * (second_point[2] - first_point[2])
    return (x, y)
end

"""
The marching-squares lookup table that returns the edge crossings to connect
for a cell sign pattern.
"""
function _edge_pairs(case_index::Int, center_value::Float64)
    case_index == 1 && return ((4, 1),) # left-bottom
    case_index == 2 && return ((1, 2),) # bottom-right
    case_index == 3 && return ((4, 2),) # left-right
    case_index == 4 && return ((2, 3),) # right-top
    if case_index == 5
        return center_value >= 0.0 ?
            ((1, 2), (3, 4)) : # bottom-right; top-left
            ((4, 1), (2, 3)) # left-bottom; right-top
    end
    case_index == 6 && return ((1, 3),) # bottom-top
    case_index == 7 && return ((3, 4),) # top-left
    case_index == 8 && return ((3, 4),) # top-left
    case_index == 9 && return ((1, 3),) # bottom-top
    if case_index == 10
        return center_value >= 0.0 ?
            ((4, 1), (2, 3)) : # left-bottom; right-top
            ((1, 2), (3, 4)) # bottom-right; top-left
    end
    case_index == 11 && return ((2, 3),) # right-top
    case_index == 12 && return ((4, 2),) # left-right
    case_index == 13 && return ((1, 2),) # bottom-right
    case_index == 14 && return ((4, 1),) # left-bottom
    return ()
end

"""
    level_contours(plane, values; level=0)

Return linearly interpolated marching-squares segments for one scalar level.
Cells with nonfinite values are skipped.
"""
function level_contours(
    plane::ParameterPlane,
    values::AbstractMatrix{<:Real};
    level::Real = 0.0,
)
    size(values) == (length(plane.y), length(plane.x)) ||
        throw(DimensionMismatch("the dimensions of values must match the parameter plane"))

    contour_level = Float64(level)
    segments = ContourSegment{Float64}[]
    sizehint!(segments, 4 * (length(plane.x) + length(plane.y)))

    for x_index in 1:(length(plane.x) - 1)
        x0 = Float64(plane.x[x_index])
        x1 = Float64(plane.x[x_index + 1])
        for y_index in 1:(length(plane.y) - 1)
            y0 = Float64(plane.y[y_index])
            y1 = Float64(plane.y[y_index + 1])

            bottom_left = Float64(values[y_index, x_index]) - contour_level
            bottom_right = Float64(values[y_index, x_index + 1]) - contour_level
            top_right = Float64(values[y_index + 1, x_index + 1]) - contour_level
            top_left = Float64(values[y_index + 1, x_index]) - contour_level
            corner_values = (bottom_left, bottom_right, top_right, top_left)
            all(isfinite, corner_values) || continue

            case_index =
                (bottom_left >= 0.0 ? 1 : 0) +
                (bottom_right >= 0.0 ? 2 : 0) +
                (top_right >= 0.0 ? 4 : 0) +
                (top_left >= 0.0 ? 8 : 0)
            (case_index == 0 || case_index == 15) && continue

            corners = (
                (x0, y0),
                (x1, y0),
                (x1, y1),
                (x0, y1),
            )
            edge_corners = ((1, 2), (2, 3), (3, 4), (4, 1))
            edge_points = Vector{Union{Nothing,NTuple{2,Float64}}}(undef, 4)
            for edge in 1:4
                first_corner, second_corner = edge_corners[edge]
                edge_points[edge] = _edge_point(
                    corner_values[first_corner],
                    corner_values[second_corner],
                    corners[first_corner],
                    corners[second_corner],
                )
            end

            center_value = sum(corner_values) / 4
            for (first_edge, second_edge) in _edge_pairs(case_index, center_value)
                first_point = edge_points[first_edge]
                second_point = edge_points[second_edge]
                (isnothing(first_point) || isnothing(second_point)) && continue
                push!(
                    segments,
                    ContourSegment(
                        first_point[1],
                        first_point[2],
                        second_point[1],
                        second_point[2],
                    ),
                )
            end
        end
    end
    return segments
end

"""
    boolean_contours(plane, values)

Return marching-squares segments along boundaries between `false` and `true`
grid values.
"""
function boolean_contours(
    plane::ParameterPlane,
    values::AbstractMatrix{Bool},
)
    return level_contours(plane, values; level = 0.5)
end

"""
    add_level_contours!(diagram, values; source, iterate, level)

Extract one contour level from `values` and append it to `diagram`.
"""
function add_level_contours!(
    diagram::KneadingDiagram,
    values::AbstractMatrix{<:Real};
    source::Symbol,
    iterate::Integer,
    level::Real,
)
    iterate >= 1 || throw(ArgumentError("the iterate must be positive"))
    contour_level = Float64(level)
    segments = level_contours(diagram.plane, values; level = contour_level)
    push!(
        diagram.layers,
        ContourLayer(source, Int(iterate), contour_level, segments),
    )
    return diagram
end

"""
    plot_kneading_contours(diagram; kwargs...)

Plot a kneading diagram. Load CairoMakie to activate this method.
"""
function plot_kneading_contours(args...; kwargs...)
    throw(ArgumentError("load CairoMakie before plotting a kneading diagram"))
end

"""
    save_kneading_contours(path, diagram; kwargs...)

Plot and save a kneading diagram. Load CairoMakie to activate this method.
"""
function save_kneading_contours(args...; kwargs...)
    throw(ArgumentError("load CairoMakie before saving a kneading diagram"))
end
