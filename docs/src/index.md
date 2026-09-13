# Kneading.jl

```@meta
CurrentModule = Kneading
```

Kneading.jl provides symbolic-dynamics tools for maps and flow-derived return
data. Its current scope includes kneading data, finite kneading determinants,
certified finite-surrogate entropy estimates, real-saddle and saddle-focus
initialization, flow-normal tangent integration, extremum-event orientation
words, parameter-plane scans, and kneading-diagram contours.

The package does not currently support weighted kneading theory or generalized
topological pressures.

```@docs
Kneading
```

## Installation

Kneading.jl is currently installed directly from GitHub:

```julia
import Pkg
Pkg.add(url = "https://github.com/hinsley/Kneading.jl")
```

## Package structure

The package has five public namespaces:

- [`Kneading.OneDimensionalMaps`](@ref) contains interval maps, kneading
  algebra, and entropy estimates.
- [`Kneading.Diagrams`](@ref) contains parameter-plane scans, contour
  geometry, and optional CairoMakie plotting methods.
- [`Kneading.RealSaddleInitialization`](@ref) validates a real saddle and
  constructs the displaced orbit and tangent seed.
- [`Kneading.FlowNormalTangents`](@ref) integrates one unit tangent transverse
  to the current flow.
- [`Kneading.FlowKneading`](@ref) captures orientation words along critical
  orbits and continues their initializations across parameter planes.

```@docs
OneDimensionalMaps
Diagrams
RealSaddleInitialization
FlowNormalTangents
FlowKneading
```

## Tutorials

- [Calculate a Chebyshev cubic kneading diagram](@ref chebyshev-cubic-tutorial) follows both critical
  orbits across a two-parameter plane, extracts critical-relation contours,
  and renders the finished diagram.
- [Flow kneading](@ref flow-kneading) includes executable Lorenz and Rössler
  examples, initialization diagnostics, and a Rössler diagram workflow.

```@contents
Pages = [
    "tutorials/chebyshev-cubic-kneading-diagram.md",
    "one-dimensional-maps.md",
    "flow-kneading.md",
    "diagrams.md",
]
Depth = 2
```

## Source

Kneading.jl is developed on
[GitHub](https://github.com/hinsley/Kneading.jl).
