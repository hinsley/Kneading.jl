# Kneading.jl
Symbolic dynamics for one-dimensional maps

The current intention is to support the calculation of topological entropy for piecewise-continuous, piecewise-monotone self-maps of an interval, as well as to support the rendering of kneading diagrams for both 1D maps and flows with 1D return map reductions (e.g., the Lorenz family).

Currently, this package does not support weighted kneading theory or
generalized topological pressures.

# Package namespaces

The package is divided into two public namespaces:

- `Kneading.OneDimensionalMaps` provides interval maps, kneading algebra, and
  entropy estimates.
- `Kneading.Diagrams` provides parameter-plane scans, contour geometry, and
  plotting hooks.

# One-dimensional example

The snippet below uses kneading itineraries truncated at 10 iterates to calculate an approximation of the topological entropy for the one-dimensional map defined piecewise by

$$
f(x) =
\begin{cases}
x/2, & 0 \leq x < 1/3, \\
3x-1, & 1/3 \leq x < 2/3, \\
3-3x, & 2/3 \leq x \leq 1.
\end{cases}
$$

```jl
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

domain = (0 // 1, 1 // 1)
partition_points = [1 // 3, 2 // 3]
partition = LapPartition(domain, partition_points)

interval_map = PartitionedIntervalMap(f, partition)

@show interval_map.orientations # LapOrientation[Increasing, Increasing, Decreasing]

data = kneading_data(interval_map, 10)

@show size(data.itineraries)    # (2, 2)
@show length(data[1, LeftSide]) # 11

matrix = kneading_matrix(data)

@show size(matrix)         # (2, 3)
@show eltype(matrix)       # PowerSeriesJet{BigInt}
@show length(matrix[1, 1]) # 11

determinant = kneading_determinant(matrix)

@show determinant # (1 - 2*t + O(t^11)) / (1 - t)

# We can force the determinant calculation to use polynomials instead of
# power-series jets so that higher-degree terms are retained in multiplications.
#=
polynomial_determinant = kneading_determinant(
    matrix;
    mode = :polynomial,
)
# Its higher coefficients can depend on which lap is deleted.
=#

entropy = entropy_estimate(determinant)

@show entropy.estimate # 0.693147180559945...
@show entropy.entropy_interval # This interval certifies the finite surrogate, not its unknown tail.
```

# Parameter-plane scans

Use the diagram API with:

```jl
using Kneading.Diagrams
```

`ParameterPlane` defines a rectangular parameter grid. `scan_plane!`
evaluates one or more scalar fields on that grid. `level_contours` extracts
interpolated contour segments, and `KneadingDiagram` collects those segments
with their source and iterate metadata.

These operations do not depend on how a scalar field was produced. The field
can come from an interval map, a projected return map, or stored numerical
data.

Loading CairoMakie activates `plot_kneading_contours` and
`save_kneading_contours`. CairoMakie is not required for scanning or contour
extraction. See [the Chebyshev cubic example](examples/README.md) for a complete
scan and plot.
