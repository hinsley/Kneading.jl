### Directions and Orientations

@enum LapOrientation::Int8 Increasing = 1 Decreasing = -1

# The two unit tangent directions at a point of the interval.
@enum TangentDirection::Int8 LeftSide = -1 RightSide = 1

### Directed Points

"""A point whose direction breaks ties in direct comparisons."""
struct DirectedPoint{T<:Real} <: Real
    value::T
    direction::TangentDirection
end

# Arithmetic removes the direction.
Base.convert(::Type{T}, point::DirectedPoint) where {T<:Real} =
    convert(T, point.value)
Base.promote_rule(
    ::Type{DirectedPoint{T}},
    ::Type{S},
) where {T<:Real,S<:Real} = promote_type(T, S)
Base.float(point::DirectedPoint) = float(point.value)
Base.:+(point::DirectedPoint) = +point.value
Base.:-(point::DirectedPoint) = -point.value
Base.abs(point::DirectedPoint) = abs(point.value)
Base.inv(point::DirectedPoint) = inv(point.value)
Base.:+(left::DirectedPoint, right::DirectedPoint) = left.value + right.value
Base.:-(left::DirectedPoint, right::DirectedPoint) = left.value - right.value
Base.:*(left::DirectedPoint, right::DirectedPoint) = left.value * right.value
Base.:/(left::DirectedPoint, right::DirectedPoint) = left.value / right.value
Base.:^(base::DirectedPoint, exponent::Integer) = base.value ^ exponent
Base.:^(base::DirectedPoint, exponent::DirectedPoint) = base.value ^ exponent.value

# Ordered comparisons use the direction only when the values are equal.
Base.:(==)(left::DirectedPoint, right::DirectedPoint) =
    left.value == right.value && left.direction == right.direction
Base.:(==)(::DirectedPoint, ::Real) = false
Base.:(==)(::Real, ::DirectedPoint) = false
Base.:(==)(::DirectedPoint, ::Base.AbstractIrrational) = false
Base.:(==)(::Base.AbstractIrrational, ::DirectedPoint) = false

Base.:<(left::DirectedPoint, right::DirectedPoint) =
    left.value < right.value ||
    (left.value == right.value &&
     left.direction == LeftSide &&
     right.direction == RightSide)
Base.:<(point::DirectedPoint, value::Real) =
    point.value < value || (point.value == value && point.direction == LeftSide)
Base.:<(value::Real, point::DirectedPoint) =
    value < point.value || (value == point.value && point.direction == RightSide)

Base.:<=(left::DirectedPoint, right::DirectedPoint) = left < right || left == right
Base.:<=(point::DirectedPoint, value::Real) = point < value
Base.:<=(value::Real, point::DirectedPoint) = value < point

### Lap Partitions

"""An interval domain and its ordered interior partition points."""
struct LapPartition{T<:Real}
    domain::Tuple{T,T}
    points::Vector{T}

    """Create a lap partition from ordered interior points."""
    function LapPartition(
        domain::Tuple{T,T},
        points::AbstractVector{T},
    ) where {T<:Real}
        left, right = domain
        left < right || throw(ArgumentError("the domain endpoints must be correctly ordered"))

        values = collect(points)
        for index in eachindex(values)
            left < values[index] < right ||
                throw(ArgumentError("partition points must be inside the domain"))
            index == firstindex(values) && continue
            values[index - 1] < values[index] ||
                throw(ArgumentError("partition points must be strictly increasing"))
        end
        return new{T}(domain, values)
    end
end

### Partitioned Interval Maps

"""A piecewise-monotone interval map with one orientation for each lap."""
struct PartitionedIntervalMap{F,T<:Real}
    f::F
    partition::LapPartition{T}
    orientations::Vector{LapOrientation}
end

"""Return the itinerary column for a tangent direction."""
_direction_index(direction::TangentDirection) = direction == LeftSide ? 1 : 2

"""Evaluate the selected one-sided value of a map."""
function _one_sided_value(f, x::Real, direction::TangentDirection)
    output = f(DirectedPoint(x, direction))
    output isa DirectedPoint && return output.value
    output isa Real || throw(ArgumentError("an interval map must return a real number"))
    return output
end

"""Infer each lap orientation from its endpoint images."""
function _infer_orientations(f, partition::LapPartition)
    boundaries = [partition.domain[1]; partition.points; partition.domain[2]]
    orientations = LapOrientation[]

    for lap in 1:(length(boundaries) - 1)
        left_image = _one_sided_value(f, boundaries[lap], RightSide)
        right_image = _one_sided_value(f, boundaries[lap + 1], LeftSide)

        left_image == right_image &&
            throw(ArgumentError("lap $lap has equal endpoint images"))
        push!(orientations, left_image < right_image ? Increasing : Decreasing)
    end
    return orientations
end

"""Create a partitioned map and infer any omitted orientations."""
function PartitionedIntervalMap(
    f,
    partition::LapPartition;
    orientations = nothing,
)
    if isnothing(orientations)
        resolved_orientations = _infer_orientations(f, partition)
    else
        resolved_orientations = LapOrientation[orientation for orientation in orientations]
    end

    length(resolved_orientations) == length(partition.points) + 1 ||
        throw(ArgumentError("there must be one orientation for each lap"))
    return PartitionedIntervalMap(f, partition, resolved_orientations)
end

### Kneading Data

"""The underlying 1D map and the two finite lap itineraries associated with
each partition point, represented as a two-column matrix."""
struct KneadingData{M}
    interval_map::M
    itineraries::Matrix{Vector{Int}}
end

"""Return the itinerary for one directed partition point."""
Base.getindex(data::KneadingData, point::Integer, direction::TangentDirection) =
    data.itineraries[point, _direction_index(direction)]

"""Return the lap that contains a directed point."""
function _lap_index(partition::LapPartition, x::Real, direction::TangentDirection)
    left, right = partition.domain
    left <= x <= right || throw(DomainError(x, "an image is outside the domain"))
    x == left && direction == LeftSide &&
        throw(DomainError(x, "the left side is outside the domain"))
    x == right && direction == RightSide &&
        throw(DomainError(x, "the right side is outside the domain"))

    index = searchsortedfirst(partition.points, x)
    if index <= length(partition.points) && partition.points[index] == x
        return direction == LeftSide ? index : index + 1
    end
    return index
end

"""Propagate a tangent direction through one lap."""
function _next_direction(
    direction::TangentDirection,
    orientation::LapOrientation,
)
    orientation == Increasing && return direction
    return direction == LeftSide ? RightSide : LeftSide
end

"""Compute one finite lap itinerary."""
function _itinerary(
    interval_map::PartitionedIntervalMap,
    point::Real,
    direction::TangentDirection,
    iterations::Int,
)
    lap = _lap_index(interval_map.partition, point, direction)
    itinerary = Vector{Int}(undef, iterations + 1)

    for step in eachindex(itinerary)
        itinerary[step] = lap
        step == lastindex(itinerary) && break

        point = _one_sided_value(interval_map.f, point, direction)
        direction = _next_direction(direction, interval_map.orientations[lap])
        lap = _lap_index(interval_map.partition, point, direction)
    end
    return itinerary
end

"""Compute both one-sided lap itineraries through `iterations` images."""
function kneading_data(interval_map::PartitionedIntervalMap, iterations::Integer)
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative"))
    isempty(interval_map.partition.points) &&
        throw(ArgumentError("the partition must have at least one point"))

    count = Int(iterations)
    point_count = length(interval_map.partition.points)
    itineraries = Matrix{Vector{Int}}(undef, point_count, 2)

    for index in 1:point_count
        point = interval_map.partition.points[index]
        itineraries[index, _direction_index(LeftSide)] =
            _itinerary(interval_map, point, LeftSide, count)
        itineraries[index, _direction_index(RightSide)] =
            _itinerary(interval_map, point, RightSide, count)
    end
    return KneadingData(interval_map, itineraries)
end
