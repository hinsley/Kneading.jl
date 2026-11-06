# Chebyshev cubic orbit diagram at v = -1/2

Map: `f(x; u,v) = (u-v)/2 * (4x^3 - 3x) + (u+v)/2`.

The plot uses 6,001 equally spaced values of u from 0 to 1, inclusive.
Each parameter is independently initialized at x = -1/2; 32,768 iterates
are discarded, then 1,024 consecutive iterates are recorded in Float64.
The other critical point +1/2 maps directly to -1/2 on this slice, so its
forward orbit has the same tail. No continuation between parameters is used.
As with any finite orbit diagram, very long transients can remain near
neutral bifurcations. The exact endpoint u=1 sends the critical orbit to
the fixed point x=1; it is intentionally retained.

The computation implements the fixed-initial-state algorithm described at
https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/orbitdiagram/
using direct scalar iteration, without installing ChaosTools.

Run from this directory:
```sh
julia --project=. --threads=8 orbit_diagram_v_minus_half.jl
```

Outputs: CairoMakie PNG (3000 x 1800), PDF, and JLD2 data containing
the parameter vector, state matrix (saved iterate x parameter), and settings.
The plot has LaTeX axis and tick labels, no title, and black scatter points.

Validation checks: all samples finite and in [-1/2,1]; exact superstable
2-cycle at u=1/2; endpoint fixed orbit at u=1; critical image identities
verified to floating-point tolerance.
