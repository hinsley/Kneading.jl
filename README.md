# Kneading.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://hinsley.github.io/Kneading.jl/dev/)

Kneading theory for one-dimensional maps and symbolic analysis of flows with approximately one-dimensional return maps.

## Kneading for one-dimensional maps

For piecewise-continuous, piecewise-monotone self-maps of an interval, `Kneading.OneDimensionalMaps` supports:

- Lap partitions, increasing and decreasing branches, and one-sided itineraries at partition points.
- Finite kneading data, kneading matrices, and kneading determinants.
- Truncated power-series and polynomial arithmetic for kneading algebra.
- Topological entropy estimates from finite determinant approximations.

Weighted kneading theory and generalized topological pressures are not currently supported.

See the [one-dimensional maps documentation](https://hinsley.github.io/Kneading.jl/dev/one-dimensional-maps/).

## Flow kneading

For autonomous ODE systems whose attractors admit approximately one-dimensional return maps, `Kneading.FlowKneading` computes symbolic orientation words by transporting a tangent vector normal to the flow and sampling it at selected events. Support includes:

- Critical-orbit initialization from real saddles with one unstable direction.
- Saddle-focus initialization and parameter continuation of smooth return-map critical points, with seed-refinement checks.
- Variational tangent integration with projection and normalization.
- Capture of local maxima or minima of a state variable, with optional event-acceptance filters.
- Component or observable-direction signs, orientation-preservation and reversal words, event states, and return times.
- Parameter-plane scans with continuation, incomplete-word reporting, and export of scan results.

The lower-level `Kneading.RealSaddleInitialization` and
`Kneading.FlowNormalTangents` namespaces also expose orbit-seeding and
tangent-integration tools independently of symbolic encoding.

See the [flow-kneading documentation](https://hinsley.github.io/Kneading.jl/dev/flow-kneading/).

## Kneading diagrams and parameter scans

`Kneading.Diagrams` provides rectangular parameter grids, independent and continuation-based scans, and contour extraction from scalar or Boolean fields. These tools accept data from interval maps, flow-derived return maps, or other numerical calculations.

Contour layers retain their source and iterate metadata. Optional CairoMakie integration provides plotting and figure export; scanning and contour extraction do not require a plotting dependency.

See the [kneading-diagram documentation](https://hinsley.github.io/Kneading.jl/dev/diagrams/).
