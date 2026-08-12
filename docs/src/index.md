# Kneading.jl

```@meta
CurrentModule = Kneading
```

Kneading.jl provides symbolic-dynamics tools for one-dimensional maps. Its
current scope includes kneading data, finite kneading determinants, certified
finite-surrogate entropy estimates, parameter-plane scans, and kneading
diagram contours.

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

The package has two public namespaces:

- [`Kneading.OneDimensionalMaps`](@ref) contains interval maps, kneading
  algebra, and entropy estimates.
- [`Kneading.Diagrams`](@ref) contains parameter-plane scans, contour
  geometry, and optional CairoMakie plotting methods.

```@docs
OneDimensionalMaps
Diagrams
```

## Tutorials

- [Calculate a Chebyshev cubic kneading diagram](@ref chebyshev-cubic-tutorial) follows both critical
  orbits across a two-parameter plane, extracts critical-relation contours,
  and renders the finished diagram.

```@contents
Pages = [
    "tutorials/chebyshev-cubic-kneading-diagram.md",
    "one-dimensional-maps.md",
    "diagrams.md",
]
Depth = 2
```

## Source

Kneading.jl is developed on
[GitHub](https://github.com/hinsley/Kneading.jl).
