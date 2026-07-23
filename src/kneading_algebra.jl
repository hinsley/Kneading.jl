### Algebraic Types

"""
A formal power series represented by coefficients listed in ascending order, starting
with the constant term.

The two multiplicands must have the same number of coefficients. Multiplication does
not calculate higher-degree terms.
"""
struct PowerSeriesJet{T<:Integer}
    coefficients::Vector{T}

    function PowerSeriesJet(coefficients::AbstractVector{T}) where {T<:Integer}
        isempty(coefficients) && throw(ArgumentError("coefficients must not be empty"))
        return new{T}(collect(coefficients))
    end
end

"""An exact finite polynomial in ascending coefficient order."""
struct Polynomial{T}
    coefficients::Vector{T}

    function Polynomial(coefficients::AbstractVector{T}) where {T}
        values = collect(coefficients)
        if isempty(values)
            isconcretetype(T) || throw(ArgumentError(
                "an empty coefficient vector must have a concrete element type",
            ))
            push!(values, zero(T))
        end
        while length(values) > 1 && iszero(values[end])
            pop!(values)
        end
        return new{T}(values)
    end
end

### Coefficient Arithmetic

Base.eltype(::Type{PowerSeriesJet{T}}) where {T} = T
Base.eltype(::PowerSeriesJet{T}) where {T} = T
Base.length(series::PowerSeriesJet) = length(series.coefficients)
Base.iszero(series::PowerSeriesJet) = all(iszero, series.coefficients)
Base.:-(series::PowerSeriesJet) = PowerSeriesJet(-series.coefficients)

Base.eltype(::Type{Polynomial{T}}) where {T} = T
Base.eltype(::Polynomial{T}) where {T} = T
Base.length(polynomial::Polynomial) = length(polynomial.coefficients)
Base.iszero(polynomial::Polynomial) = all(iszero, polynomial.coefficients)
Base.:-(polynomial::Polynomial) = Polynomial(-polynomial.coefficients)

function Base.zero(series::PowerSeriesJet{T}) where {T}
    return PowerSeriesJet(zeros(T, length(series)))
end

function Base.one(series::PowerSeriesJet{T}) where {T}
    coefficients = zeros(T, length(series))
    coefficients[1] = one(T)
    return PowerSeriesJet(coefficients)
end

Base.zero(polynomial::Polynomial{T}) where {T} = Polynomial(T[zero(T)])
Base.one(polynomial::Polynomial{T}) where {T} = Polynomial(T[one(T)])

function _same_precision(left::PowerSeriesJet, right::PowerSeriesJet)
    length(left) == length(right) ||
        throw(DimensionMismatch("jets must have the same number of coefficients"))
    return length(left)
end

function Base.:+(left::PowerSeriesJet{T}, right::PowerSeriesJet{T}) where {T}
    _same_precision(left, right)
    return PowerSeriesJet(left.coefficients + right.coefficients)
end

function Base.:-(left::PowerSeriesJet{T}, right::PowerSeriesJet{T}) where {T}
    _same_precision(left, right)
    return PowerSeriesJet(left.coefficients - right.coefficients)
end

function Base.:*(left::PowerSeriesJet{T}, right::PowerSeriesJet{T}) where {T}
    coefficient_count = _same_precision(left, right)
    coefficients = zeros(T, coefficient_count)

    for left_index in eachindex(left.coefficients)
        maximum_right_index = coefficient_count - left_index + 1
        for right_index in 1:maximum_right_index
            coefficients[left_index + right_index - 1] +=
                left.coefficients[left_index] * right.coefficients[right_index]
        end
    end
    return PowerSeriesJet(coefficients)
end

function Base.:+(left::Polynomial{T}, right::Polynomial{T}) where {T}
    coefficients = zeros(T, max(length(left), length(right)))
    for index in eachindex(coefficients)
        if index <= length(left)
            coefficients[index] += left.coefficients[index]
        end
        if index <= length(right)
            coefficients[index] += right.coefficients[index]
        end
    end
    return Polynomial(coefficients)
end

function Base.:-(left::Polynomial{T}, right::Polynomial{T}) where {T}
    coefficients = zeros(T, max(length(left), length(right)))
    for index in eachindex(coefficients)
        if index <= length(left)
            coefficients[index] += left.coefficients[index]
        end
        if index <= length(right)
            coefficients[index] -= right.coefficients[index]
        end
    end
    return Polynomial(coefficients)
end

function Base.:*(left::Polynomial{T}, right::Polynomial{T}) where {T}
    coefficients = zeros(T, length(left) + length(right) - 1)
    for left_index in eachindex(left.coefficients)
        for right_index in eachindex(right.coefficients)
            coefficients[left_index + right_index - 1] +=
                left.coefficients[left_index] * right.coefficients[right_index]
        end
    end
    return Polynomial(coefficients)
end

function Base.:^(polynomial::Polynomial, exponent::Integer)
    exponent >= 0 || throw(ArgumentError("a polynomial exponent must be nonnegative"))
    result = one(polynomial)
    factor = polynomial
    remaining = exponent

    while remaining > 0
        if isodd(remaining)
            result *= factor
        end
        remaining >>= 1
        if remaining > 0
            factor *= factor
        end
    end
    return result
end

function _constant_coefficient(series)
    return series.coefficients[1]
end

function _assert_unit_constant(series)
    constant = _constant_coefficient(series)
    abs(constant) == one(constant) ||
        throw(ArgumentError("a Bareiss divisor must have constant coefficient ±1"))
    return nothing
end

function _divide_exact(
    numerator::PowerSeriesJet{T},
    denominator::PowerSeriesJet{T},
) where {T<:Integer}
    coefficient_count = _same_precision(numerator, denominator)
    _assert_unit_constant(denominator)
    quotient = zeros(T, coefficient_count)
    constant = _constant_coefficient(denominator)

    for degree in 0:(coefficient_count - 1)
        value = numerator.coefficients[degree + 1]
        for denominator_degree in 1:degree
            value -= denominator.coefficients[denominator_degree + 1] *
                quotient[degree - denominator_degree + 1]
        end
        quotient[degree + 1] = div(value, constant)
    end
    return PowerSeriesJet(quotient)
end

function _divide_exact(
    numerator::Polynomial{T},
    denominator::Polynomial{T},
) where {T<:Integer}
    iszero(denominator) && throw(DivideError())
    iszero(numerator) && return zero(denominator)
    length(numerator) < length(denominator) &&
        throw(ArgumentError("polynomial division is not exact"))

    remainder = copy(numerator.coefficients)
    quotient = zeros(T, length(numerator) - length(denominator) + 1)

    while !(length(remainder) == 1 && iszero(remainder[1])) &&
          length(remainder) >= length(denominator)
        shift = length(remainder) - length(denominator)
        coefficient, scalar_remainder =
            divrem(remainder[end], denominator.coefficients[end])
        iszero(scalar_remainder) ||
            throw(ArgumentError("polynomial division is not exact"))
        quotient[shift + 1] = coefficient

        for index in eachindex(denominator.coefficients)
            remainder[shift + index] -= coefficient * denominator.coefficients[index]
        end
        while length(remainder) > 1 && iszero(remainder[end])
            pop!(remainder)
        end
    end

    all(iszero, remainder) || throw(ArgumentError("polynomial division is not exact"))
    return Polynomial(quotient)
end

"""Convert a finite jet to a polynomial with zero unstored coefficients."""
polynomial_approximation(jet::PowerSeriesJet) = Polynomial(jet.coefficients)

### Kneading Matrices

"""A kneading matrix and the orientation of each of its columns."""
struct KneadingMatrix{S} <: AbstractMatrix{S}
    entries::Matrix{S}
    orientations::Vector{LapOrientation}

    function KneadingMatrix(
        entries::AbstractMatrix{S},
        orientations,
    ) where {S}
        size(entries, 1) > 0 || throw(ArgumentError("a kneading matrix needs a row"))
        size(entries, 2) == size(entries, 1) + 1 ||
            throw(DimensionMismatch("a kneading matrix must have one more column than row"))
        length(orientations) == size(entries, 2) ||
            throw(DimensionMismatch("each kneading-matrix column needs an orientation"))
        values = LapOrientation[orientation for orientation in orientations]
        return new{S}(Matrix(entries), values)
    end
end

Base.size(matrix::KneadingMatrix) = size(matrix.entries)
Base.getindex(matrix::KneadingMatrix, indices...) = getindex(matrix.entries, indices...)

"""Construct the finite kneading matrix from one-sided lap itineraries."""
function kneading_matrix(data::KneadingData)
    partition_point_count = size(data.itineraries, 1)
    lap_count = partition_point_count + 1
    coefficient_count = length(data[1, LeftSide])

    entries = [
        PowerSeriesJet(zeros(BigInt, coefficient_count))
        for _ in 1:partition_point_count, _ in 1:lap_count
    ]

    for point_index in 1:partition_point_count
        for side in (LeftSide, RightSide)
            itinerary = data[point_index, side]
            length(itinerary) == coefficient_count ||
                throw(DimensionMismatch("all itineraries must have the same length"))

            increment_sign = side == LeftSide ? -1 : 1
            orientation_sign = 1
            for coefficient_index in eachindex(itinerary)
                lap = itinerary[coefficient_index]
                entries[point_index, lap].coefficients[coefficient_index] +=
                    increment_sign * orientation_sign
                orientation_sign *= Int(data.interval_map.orientations[lap])
            end
        end
    end

    return KneadingMatrix(entries, copy(data.interval_map.orientations))
end

"""Check the constant-coefficient pattern required by a kneading matrix."""
function _validate_constant_skeleton(matrix::KneadingMatrix)
    for row in axes(matrix, 1)
        for column in axes(matrix, 2)
            expected = if column == row
                -1
            elseif column == row + 1
                1
            else
                0
            end
            _constant_coefficient(matrix[row, column]) == expected ||
                throw(ArgumentError("the kneading matrix has an invalid constant term"))
        end
    end
    return nothing
end

"""Convert all entries to `BigInt` jets with the shortest stored length."""
function _common_precision_entries(matrix::KneadingMatrix{<:PowerSeriesJet})
    coefficient_count = minimum(length, matrix.entries)
    return map(matrix.entries) do entry
        coefficients = BigInt.(entry.coefficients[1:coefficient_count])
        PowerSeriesJet(coefficients)
    end
end

### Determinant Algorithms

# Bareiss uses exact division, including truncated-series division, in O(n^3) ring operations.
function _bareiss_determinant(matrix::AbstractMatrix)
    row_count, column_count = size(matrix)
    row_count == column_count || throw(DimensionMismatch("the matrix must be square"))
    row_count > 0 || throw(ArgumentError("the matrix must not be empty"))

    work = copy(matrix)
    previous_pivot = one(work[1, 1])

    # By the definition of the kneading matrix, each Bareiss pivot
    # has constant coefficient 1 or -1.
    for pivot_index in 1:(row_count - 1)
        pivot = work[pivot_index, pivot_index]
        _assert_unit_constant(pivot)

        for row in (pivot_index + 1):row_count
            for column in (pivot_index + 1):column_count
                numerator = pivot * work[row, column] -
                    work[row, pivot_index] * work[pivot_index, column]
                work[row, column] = _divide_exact(numerator, previous_pivot)
            end
        end
        previous_pivot = pivot
    end

    determinant = work[end, end]
    _assert_unit_constant(determinant)
    return determinant
end

"""Compute a dot product with ring addition and multiplication."""
function _ring_dot(left, right)
    value = zero(left[1])
    for index in eachindex(left, right)
        value += left[index] * right[index]
    end
    return value
end

"""Multiply a matrix and vector with ring addition and multiplication."""
function _ring_matvec(matrix, vector)
    result = [zero(vector[1]) for _ in axes(matrix, 1)]
    for row in axes(matrix, 1)
        for column in axes(matrix, 2)
            result[row] += matrix[row, column] * vector[column]
        end
    end
    return result
end

# Berkowitz computes the characteristic polynomial with only additions and
# multiplications in O(n^4) ring operations.
function _berkowitz_characteristic_polynomial(matrix::AbstractMatrix)
    order = size(matrix, 1)
    order == size(matrix, 2) || throw(DimensionMismatch("the matrix must be square"))
    order > 0 || throw(ArgumentError("the matrix must not be empty"))
    order == 1 && return [one(matrix[1, 1]), -matrix[1, 1]]

    trailing_matrix = matrix[2:end, 2:end]
    trailing_coefficients = _berkowitz_characteristic_polynomial(trailing_matrix)
    toeplitz_column = Vector{eltype(matrix)}(undef, order + 1)
    toeplitz_column[1] = one(matrix[1, 1])
    toeplitz_column[2] = -matrix[1, 1]

    row = collect(matrix[1, 2:end])
    powered_column = collect(matrix[2:end, 1])
    for index in 3:(order + 1)
        toeplitz_column[index] = -_ring_dot(row, powered_column)
        if index <= order
            powered_column = _ring_matvec(trailing_matrix, powered_column)
        end
    end

    coefficients = [zero(matrix[1, 1]) for _ in 1:(order + 1)]
    for output_index in 0:order
        for trailing_index in 0:min(output_index, order - 1)
            coefficients[output_index + 1] +=
                toeplitz_column[output_index - trailing_index + 1] *
                trailing_coefficients[trailing_index + 1]
        end
    end
    return coefficients
end

function _berkowitz_determinant(matrix::AbstractMatrix)
    coefficients = _berkowitz_characteristic_polynomial(matrix)
    determinant = isodd(size(matrix, 1)) ? -coefficients[end] : coefficients[end]
    return determinant
end

### Kneading Determinants

# Keep P/(1 - εt) unexpanded: roots in (0,1) are numerator roots.
"""A normalized kneading determinant `P(t) / (1 - εt)`."""
struct KneadingDeterminant{P}
    numerator::P
    epsilon::LapOrientation
    deleted_lap::Int

    function KneadingDeterminant(
        numerator::P,
        epsilon::LapOrientation,
        deleted_lap::Integer,
    ) where {P}
        _constant_coefficient(numerator) == one(_constant_coefficient(numerator)) ||
            throw(ArgumentError("a kneading numerator must start with 1"))
        deleted_lap >= 1 || throw(ArgumentError("the deleted lap must be positive"))
        return new{P}(numerator, epsilon, Int(deleted_lap))
    end
end

function polynomial_approximation(
    determinant::KneadingDeterminant{<:PowerSeriesJet},
)
    numerator = polynomial_approximation(determinant.numerator)
    return KneadingDeterminant(
        numerator,
        determinant.epsilon,
        determinant.deleted_lap,
    )
end

function Base.denominator(determinant::KneadingDeterminant)
    T = eltype(determinant.numerator)
    return Polynomial(T[one(T), -convert(T, Int(determinant.epsilon))])
end

function _minor_determinant(matrix, backend)
    # A modular/CRT backend can replace BigInt arithmetic if coefficients grow.
    backend == :bareiss && return _bareiss_determinant(matrix)
    backend == :berkowitz && return _berkowitz_determinant(matrix)
    throw(ArgumentError("the supported backends are :bareiss and :berkowitz"))
end

"""Compute a normalized kneading determinant from finite kneading data."""
function kneading_determinant(
    matrix::KneadingMatrix{<:PowerSeriesJet};
    mode = :jet,
    deleted_lap::Integer = 1,
    backend = :bareiss,
)
    _validate_constant_skeleton(matrix)
    lap_count = size(matrix, 2)
    1 <= deleted_lap <= lap_count || throw(BoundsError(matrix, (:, deleted_lap)))
    columns = [column for column in 1:lap_count if column != deleted_lap]
    jet_entries = _common_precision_entries(matrix)

    # Deleting lap j leaves diagonal -1 before j and +1 after j at t = 0.
    # The bidiagonal blocks disconnect, so every leading minor has constant ±1.

    working_entries = if mode == :jet
        jet_entries
    elseif mode == :polynomial
        # Above the jet order, this surrogate can depend on the deleted lap.
        map(polynomial_approximation, jet_entries)
    else
        throw(ArgumentError("the supported modes are :jet and :polynomial"))
    end

    numerator = _minor_determinant(working_entries[:, columns], backend)
    iseven(deleted_lap) && (numerator = -numerator)
    _constant_coefficient(numerator) == one(_constant_coefficient(numerator)) ||
        throw(AssertionError("a normalized kneading numerator must start with 1"))

    return KneadingDeterminant(
        numerator,
        matrix.orientations[deleted_lap],
        Int(deleted_lap),
    )
end

### Display

function _show_series(io::IO, coefficients; remainder_order = nothing)
    printed = false
    for degree in 0:(length(coefficients) - 1)
        coefficient = coefficients[degree + 1]
        iszero(coefficient) && continue
        negative = coefficient < zero(coefficient)
        magnitude = abs(coefficient)

        if printed
            print(io, negative ? " - " : " + ")
        elseif negative
            print(io, "-")
        end

        if degree == 0 || magnitude != one(magnitude)
            print(io, magnitude)
            degree > 0 && print(io, "*")
        end
        if degree > 0
            print(io, "t")
            degree > 1 && print(io, "^", degree)
        end
        printed = true
    end

    printed || print(io, "0")
    isnothing(remainder_order) || print(io, " + O(t^", remainder_order, ")")
end

Base.show(io::IO, polynomial::Polynomial) =
    _show_series(io, polynomial.coefficients)

Base.show(io::IO, jet::PowerSeriesJet) =
    _show_series(io, jet.coefficients; remainder_order = length(jet))

function Base.show(io::IO, determinant::KneadingDeterminant)
    print(io, "(")
    show(io, determinant.numerator)
    print(io, ") / (")
    show(io, denominator(determinant))
    print(io, ")")
end
