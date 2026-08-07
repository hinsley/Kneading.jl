# [Kneading diagrams](@id Kneading.Diagrams)

`Kneading.Diagrams` separates parameter-plane sampling from contour
extraction. A sampled scalar field can come from an interval map, a return-map
reduction of a flow, or stored numerical data.

```@meta
CurrentModule = Kneading.Diagrams
```

## Parameter-plane scan

`scan_plane!` evaluates one or more scalar fields at every Cartesian grid
point. The first array dimension corresponds to the plane's `y` values and the
second corresponds to its `x` values.

```@example diagrams
using Kneading.Diagrams

plane = ParameterPlane(range(-1, 1; length = 21), range(-1, 1; length = 21))
samples = zeros(length(plane.y), length(plane.x), 1)

scan_plane!(samples, plane) do output, x, y
    output[1] = x^2 + y^2
end

segments = level_contours(plane, @view(samples[:, :, 1]); level = 0.5)
length(segments)
```

Each independent point in `scan_plane!` may run concurrently. Use
`scan_continuation!` when state must instead advance sequentially along each
row or column.

## Planes and scans

```@docs
ParameterPlane
scan_plane!
scan_continuation!
```

## Contours and diagrams

```@docs
ContourSegment
ContourLayer
KneadingDiagram
level_contours
boolean_contours
add_level_contours!
```

## Plotting

Load CairoMakie to activate the plotting extension:

```julia
using CairoMakie
using Kneading.Diagrams

plot_kneading_contours(diagram)
save_kneading_contours("diagram.png", diagram)
```

```@docs
plot_kneading_contours
save_kneading_contours
```

The complete Chebyshev cubic scan is available in the
[examples directory](https://github.com/hinsley/Kneading.jl/tree/main/examples).
