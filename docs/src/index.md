# Kneading.jl

Kneading.jl provides symbolic-dynamics tools for one-dimensional maps. Its
current scope includes kneading data, finite kneading determinants, certified
finite-surrogate entropy estimates, parameter-plane scans, and kneading
diagram contours.

The package does not currently support weighted kneading theory or generalized
topological pressures.

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

```@contents
Pages = ["one-dimensional-maps.md", "diagrams.md"]
Depth = 2
```

## Source

Kneading.jl is developed on
[GitHub](https://github.com/hinsley/Kneading.jl).
