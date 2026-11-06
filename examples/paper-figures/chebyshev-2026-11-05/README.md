# Chebyshev cubic: maximal critical-orbit Lyapunov exponent and kneading contours

The figure is generated with CairoMakie, with no title and LaTeXStrings for
all mathematical axis labels and ticks. The parameter scan has a resolution of 3000 × 3000. The final PNG is
3000 pixels wide, with its height fitted tightly to the plot and colorbar.
Ticks and labels are enlarged by approximately 50% relative to the first export. The scan covers the closed
rectangle u ∈ [0,1], v ∈ [-1,0], including its edges.

## Calculation

f(u,v,x) = (u-v)/2 * (4x³-3x) + (u+v)/2,
f′(u,v,x) = (u-v)(6x²-3/2).

Follow the two critical orbits starting at c₋ = -1/2 and c₊ = 1/2.
Their first images are exactly u and v. Discard 8,192 further iterations,
then average log |f′| over 16,384 iterates, separately for each orbit.
The background is max(λ₋, λ₊). Never include the zero derivative at the
initial critical point in the Lyapunov average. These are finite-time
estimates, in natural logarithms per iterate. A derivative equal to zero
is retained as -∞, without an artificial logarithm floor.

The whole rectangle has u,v ∈ [-1,1], so the critical orbits remain in the
forward-invariant interval [-1,1]. No escape mask is needed here.

Computation uses Float64 on the NVIDIA GeForce RTX 4070 Ti. See
validation.txt for analytic tests, independent CPU comparisons, and
sampled checks with twice the transient and averaging length.

## Contours and color scale

Contours use Kneading.jl's actual level_contours / add_level_contours!
implementation and its documented Chebyshev cubic tutorial. Orbit entries
2 through 20 represent fⁿ(cₛ)=cₜ for n=1,...,19 and cₛ,cₜ∈{-1/2,+1/2}.
There are 76 layers, drawn in descending orbit-entry order so earlier
iterates appear on top. Styles match the documented example:

- Left critical orbit: RGB (0.85, 0.15, 0.12).
- Right critical orbit: RGB (0.10, 0.30, 0.90).
- Opacity: n^(-1.2), where n is the number of map applications.
- Constant line width 2 at the 1000-unit logical figure size, as in the
  package's CairoMakie extension; export uses px_per_unit=3.

The source page uses tanh(λ/0.5) saturation, with green #24c77a at the
negative limit, black at zero, and gold #f2c94c at the positive limit.
The colorbar uses the same transformation and is labeled in λ values,
not transformed values. The finite ticks are -0.5, -0.25, 0, 0.25, and 0.5;
the asymptotic endpoints are left unlabeled.

## Reproduce

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --threads=4 --project=. compute.jl
julia --startup-file=no --threads=4 --project=. render.jl
```

Project.toml and Manifest.toml preserve the package environment. The local
Kneading.jl checkout is pinned to f2d9add0cfc3651fc38dfee928f08ef240d9e89f.
The .jld2 files retain both Lyapunov fields, their maximum, parameter
coordinates, calculation settings, and the extracted contour geometry.

Sources:
- https://hinsley.github.io/dynamical-systems/projects/chebyshev-cubic-kneading-diagram.html
- https://github.com/hinsley/Kneading.jl/blob/f2d9add0cfc3651fc38dfee928f08ef240d9e89f/docs/src/tutorials/chebyshev-cubic-kneading-diagram.md
- https://github.com/hinsley/Kneading.jl/blob/f2d9add0cfc3651fc38dfee928f08ef240d9e89f/ext/KneadingCairoMakieExt.jl

The fetched source page is preserved as source-page.html. Its Lyapunov
burn-in and averaging lengths are not exposed in that HTML, so this scan
uses the explicit settings above rather than assuming the original lengths.
