# Kneading diagram examples

`Kneading.Diagrams` exports `ParameterPlane`, `scan_plane!`,
`level_contours`, `add_level_contours!`, and `KneadingDiagram`. Loading
CairoMakie activates `plot_kneading_contours` and
`save_kneading_contours`.

`chebyshev_cubic_kneading.jl` uses the public package APIs to scan the two
critical orbits of

$$
f_{u,v}(x)
=
\frac{u-v}{2}\left(4x^3-3x\right)
+
\frac{u+v}{2}.
$$

The critical points are $-1/2$ and $1/2$. The example draws a contour whenever
an iterate of either critical point crosses either critical point. Red curves
come from the left critical orbit. Blue curves come from the right critical
orbit.

The same package API accepts scalar fields produced by other computations. The
scan and plotting code does not depend on the Chebyshev family.

Install the plotting dependency:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

Generate the full `1000 x 1000`, 20-iterate diagram:

```bash
julia --project=examples examples/chebyshev_cubic_kneading.jl
```

Run a smaller scan:

```bash
CHEBYSHEV_GRID_SIZE=200 CHEBYSHEV_ITERATES=8 \
    julia --project=examples examples/chebyshev_cubic_kneading.jl
```

Run the dependency-free example tests from the repository root:

```bash
julia --project=. examples/test_chebyshev_cubic_kneading.jl
```

## Rössler flow kneading

`rossler_kneading_diagram.jl` uses the public saddle-focus initializer and
flow-kneading scan API for the origin-fixed Rössler system. From the repository
root, run:

```sh
julia --project=. --threads=auto examples/rossler_kneading_diagram.jl
julia --project=examples examples/plot_rossler_kneading.jl output/rossler-flow-kneading.tsv
```

The default grid is 256 by 256. Set `ROSSLER_RESOLUTION=32` for a quick trial
and `ROSSLER_OUTPUT` to choose the TSV destination. Plotting uses the optional
environment installed above and writes raw-sign and transition-word images
under `output/figures/`, leaving incomplete words white.

The example reproduces the reference protocol: fixed event index 4, RK4 step
0.05, and eight raw signs including the initial critical event. New calculations
can instead use adaptive integration and automatic seed refinement, as shown in
the [flow-kneading guide](../docs/src/flow-kneading.md).

Run the checked-in reference fixtures independently with:

```sh
julia --project=. examples/verify_rossler_reference.jl
```
