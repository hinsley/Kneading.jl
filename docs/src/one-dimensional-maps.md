# [One-dimensional maps](@id Kneading.OneDimensionalMaps)

`Kneading.OneDimensionalMaps` represents piecewise-monotone self-maps of an
interval and derives finite kneading data, determinants, and entropy
estimates.

```@meta
CurrentModule = Kneading.OneDimensionalMaps
```

## Example

The following map has three monotone laps:

```math
f(x) =
\begin{cases}
x/2, & 0 \leq x < 1/3, \\
3x-1, & 1/3 \leq x < 2/3, \\
3-3x, & 2/3 \leq x \leq 1.
\end{cases}
```

```@example interval-map
using Kneading.OneDimensionalMaps

function f(x)
    if x < 1 // 3
        return x / 2
    elseif x < 2 // 3
        return 3x - 1
    else
        return 3 - 3x
    end
end

partition = LapPartition((0 // 1, 1 // 1), [1 // 3, 2 // 3])
interval_map = PartitionedIntervalMap(f, partition)
data = kneading_data(interval_map, 10)
matrix = kneading_matrix(data)
determinant = kneading_determinant(matrix)
entropy = entropy_estimate(determinant)

(determinant = determinant, entropy = entropy.estimate)
```

The returned entropy interval certifies the finite determinant surrogate. It
does not certify the unknown infinite tail of the kneading series.

## Directions and orientations

`LapOrientation` has the values `Increasing` and `Decreasing`.
`TangentDirection` has the values `LeftSide` and `RightSide`; these values
select the one-sided orbit at a partition point.

## Interval maps and kneading data

```@docs
LapPartition
PartitionedIntervalMap
KneadingData
kneading_data
```

## Kneading algebra

```@docs
PowerSeriesJet
Polynomial
KneadingMatrix
kneading_matrix
KneadingDeterminant
kneading_determinant
polynomial_approximation
```

## Entropy estimates

```@docs
RootInterval
EntropyEstimate
real_root_intervals
entropy_estimate
```
