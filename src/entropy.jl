### Topological Entropy

# For positive entropy, h = -log(r), where r is the smallest numerator root in
# (0, 1). A finite series jet is converted to a polynomial by setting its
# unstored coefficients to zero. The calculation then uses exact rational
# arithmetic to remove repeated factors, count and isolate every root with a
# Sturm sequence, refine the smallest root, and apply -log with outward
# rounding. The reported interval applies to the finite polynomial surrogate,
# not its unknown tail.
#
# Main call path:
# 1. entropy_estimate obtains the numerator polynomial, using
#    polynomial_approximation when the numerator is a jet.
# 2. real_root_intervals isolates its roots in (0, 1).
# 3. _refine_root narrows the interval containing the smallest root.
# 4. _entropy_interval converts the root bounds to entropy bounds.
#
# Newton's method needs a suitable initial guess and can converge to the wrong
# root. A bracketed search needs a known sign-changing interval and can miss a
# repeated root. Sturm isolation uses exact root counts to identify the smallest
# root before its interval is narrowed.

### Result Types

"""An exact rational interval containing one real root."""
struct RootInterval
    lower::Rational{BigInt}
    upper::Rational{BigInt}
end

"""A finite-surrogate entropy estimate and its certified interval."""
struct EntropyEstimate
    root_interval::RootInterval
    entropy_interval::Tuple{BigFloat,BigFloat}
    estimate::BigFloat
end

### Polynomial Operations

"""Evaluate a polynomial at a point with Horner's method."""
function _evaluate(polynomial::Polynomial, point)
    value = zero(point)
    for coefficient in Iterators.reverse(polynomial.coefficients)
        value = value * point + coefficient
    end
    return value
end

"""Return the formal derivative of a polynomial."""
function _derivative(polynomial::Polynomial{T}) where {T}
    length(polynomial) == 1 && return zero(polynomial)
    coefficients = [
        (degree * polynomial.coefficients[degree + 1])
        for degree in 1:(length(polynomial) - 1)
    ]
    return Polynomial(coefficients)
end

"""Convert all coefficients to exact `Rational{BigInt}` values."""
function _rational_polynomial(polynomial::Polynomial)
    coefficients = Rational{BigInt}[
        convert(Rational{BigInt}, coefficient)
        for coefficient in polynomial.coefficients
    ]
    return Polynomial(coefficients)
end

"""Divide two rational polynomials and return their quotient and remainder."""
function _polynomial_divrem(
    numerator::Polynomial{T},
    denominator::Polynomial{T},
) where {T<:Rational}
    iszero(denominator) && throw(DivideError())
    if length(numerator) < length(denominator)
        return zero(numerator), numerator
    end

    remainder = copy(numerator.coefficients)
    quotient = zeros(T, length(numerator) - length(denominator) + 1)

    while !(length(remainder) == 1 && iszero(remainder[1])) &&
          length(remainder) >= length(denominator)
        shift = length(remainder) - length(denominator)
        coefficient = remainder[end] / denominator.coefficients[end]
        quotient[shift + 1] = coefficient

        for index in eachindex(denominator.coefficients)
            remainder[shift + index] -= coefficient * denominator.coefficients[index]
        end
        while length(remainder) > 1 && iszero(remainder[end])
            pop!(remainder)
        end
    end

    return Polynomial(quotient), Polynomial(remainder)
end

"""Normalize a nonzero polynomial to leading coefficient one."""
function _monic(polynomial::Polynomial)
    iszero(polynomial) && return polynomial
    leading_coefficient = polynomial.coefficients[end]
    return Polynomial(polynomial.coefficients / leading_coefficient)
end

"""Return the monic greatest common divisor of two polynomials."""
function _polynomial_gcd(left::Polynomial, right::Polynomial)
    while !iszero(right)
        _, remainder = _polynomial_divrem(left, right)
        left, right = right, remainder
    end
    return _monic(left)
end

"""Remove repeated factors from a polynomial."""
function _square_free_part(polynomial::Polynomial)
    derivative = _derivative(polynomial)
    iszero(derivative) && return _monic(polynomial)
    common_factor = _polynomial_gcd(polynomial, derivative)
    quotient, remainder = _polynomial_divrem(polynomial, common_factor)
    iszero(remainder) || throw(AssertionError("square-free division must be exact"))
    return _monic(quotient)
end

"""Remove all copies of a root at an interval endpoint."""
function _remove_endpoint_root(polynomial::Polynomial{T}, endpoint::T) where {T}
    factor = Polynomial(T[-endpoint, one(T)])
    while _evaluate(polynomial, endpoint) == zero(T)
        polynomial, remainder = _polynomial_divrem(polynomial, factor)
        iszero(remainder) || throw(AssertionError("endpoint division must be exact"))
    end
    return polynomial
end

"""Construct the Sturm sequence of a polynomial."""
function _sturm_sequence(polynomial::Polynomial)
    sequence = [polynomial, _derivative(polynomial)]
    while !iszero(sequence[end])
        _, remainder = _polynomial_divrem(sequence[end - 1], sequence[end])
        iszero(remainder) && break
        push!(sequence, -remainder)
    end
    return sequence
end

### Sturm Root Isolation

"""Count nonzero sign changes in a Sturm sequence at a point."""
function _sign_variations(sequence, point)
    changes = 0
    previous_sign = 0

    for polynomial in sequence
        value = _evaluate(polynomial, point)
        current_sign = value < 0 ? -1 : value > 0 ? 1 : 0
        current_sign == 0 && continue
        if previous_sign != 0 && current_sign != previous_sign
            changes += 1
        end
        previous_sign = current_sign
    end
    return changes
end

"""Count distinct roots in an open interval from a Sturm sequence."""
_root_count(sequence, left, right) =
    _sign_variations(sequence, left) - _sign_variations(sequence, right)

"""Choose an interval split that is not a polynomial root."""
function _nonroot_split(polynomial, left, right)
    split = (left + right) / 2
    while iszero(_evaluate(polynomial, split))
        split = (left + split) / 2
    end
    return split
end

"""Append intervals that each contain one distinct root."""
function _isolate_roots!(intervals, polynomial, sequence, left, right, count)
    count == 0 && return intervals
    if count == 1
        push!(intervals, RootInterval(left, right))
        return intervals
    end

    split = _nonroot_split(polynomial, left, right)
    left_count = _root_count(sequence, left, split)
    _isolate_roots!(intervals, polynomial, sequence, left, split, left_count)
    _isolate_roots!(intervals, polynomial, sequence, split, right, count - left_count)
    return intervals
end

"""Isolate all distinct polynomial roots in an open rational interval."""
function real_root_intervals(
    polynomial::Polynomial;
    interval = (BigInt(0) // 1, BigInt(1) // 1),
)
    iszero(polynomial) && throw(ArgumentError("the zero polynomial has no isolated roots"))
    left = convert(Rational{BigInt}, interval[1])
    right = convert(Rational{BigInt}, interval[2])
    left < right || throw(ArgumentError("the interval endpoints must be ordered"))

    square_free = _square_free_part(_rational_polynomial(polynomial))
    square_free = _remove_endpoint_root(square_free, left)
    square_free = _remove_endpoint_root(square_free, right)
    length(square_free) == 1 && return RootInterval[]

    sequence = _sturm_sequence(square_free)
    count = _root_count(sequence, left, right)
    intervals = RootInterval[]
    return _isolate_roots!(intervals, square_free, sequence, left, right, count)
end

### Root Refinement

"""Narrow an isolated root interval to the requested relative precision."""
function _refine_root(polynomial::Polynomial, interval::RootInterval, bits::Integer)
    bits >= 1 || throw(ArgumentError("bits must be positive"))
    rational_polynomial = _square_free_part(_rational_polynomial(polynomial))
    left, right = interval.lower, interval.upper
    rational_polynomial = _remove_endpoint_root(rational_polynomial, left)
    rational_polynomial = _remove_endpoint_root(rational_polynomial, right)
    left_value = _evaluate(rational_polynomial, left)
    right_value = _evaluate(rational_polynomial, right)
    sign(left_value) == -sign(right_value) ||
        throw(AssertionError("an isolated simple root must be bracketed by a sign change"))

    scale = BigInt(1) << bits
    while left != right && (right - left) * scale > left
        midpoint = (left + right) / 2
        midpoint_value = _evaluate(rational_polynomial, midpoint)
        if iszero(midpoint_value)
            left = midpoint
            right = midpoint
        elseif sign(midpoint_value) == sign(left_value)
            left = midpoint
            left_value = midpoint_value
        else
            right = midpoint
            right_value = midpoint_value
        end
    end
    return RootInterval(left, right)
end

### Entropy Bounds

"""Convert root bounds to outward-rounded entropy bounds."""
function _entropy_interval(root::RootInterval, bits::Integer)
    precision = Int(bits) + 32
    return setprecision(BigFloat, precision) do
        upper_root = BigFloat(root.upper, RoundUp)
        lower_entropy = -setrounding(BigFloat, RoundUp) do
            log(upper_root)
        end

        lower_root = BigFloat(root.lower, RoundDown)
        upper_entropy = -setrounding(BigFloat, RoundDown) do
            log(lower_root)
        end
        (lower_entropy, upper_entropy)
    end
end

### Entropy Estimates

"""Estimate entropy from the smallest polynomial root in `(0, 1)`."""
function _entropy_estimate(polynomial::Polynomial; bits::Integer)
    intervals = real_root_intervals(polynomial)
    isempty(intervals) && return nothing
    root_interval = _refine_root(polynomial, first(intervals), bits)
    entropy_interval = _entropy_interval(root_interval, bits)
    estimate = (entropy_interval[1] + entropy_interval[2]) / 2
    return EntropyEstimate(root_interval, entropy_interval, estimate)
end

"""Estimate entropy from a finite-jet kneading determinant."""
function entropy_estimate(
    determinant::KneadingDeterminant{<:PowerSeriesJet};
    bits::Integer = 128,
)
    approximation = polynomial_approximation(determinant)
    return _entropy_estimate(approximation.numerator; bits)
end

"""Estimate entropy from a polynomial kneading determinant."""
function entropy_estimate(
    determinant::KneadingDeterminant{<:Polynomial};
    bits::Integer = 128,
)
    # Sturm isolation certifies roots of this finite surrogate.
    return _entropy_estimate(determinant.numerator; bits)
end
