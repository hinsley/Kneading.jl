# Topological entropy, Figure 6 parameter region

**Current figure:** adaptively refined; see `HTOP-refinement.md` and `entropy-refinement-final.json`. The final data are in `entropy_boundary_refined_final_3000.npz`. A total of 186,583 cells were recomputed; only two retain the 0.01 fallback. The code below describes the initial degree-64 scan retained for provenance. `./replot_HTOP.sh` uses the current refined display cache and does not repeat these computations.

New computation on 2026-11-05. Region u in [0,0.99], v in [-0.99,0], 3000 x 3000 samples uniform in (-log(1-u), log(1+v)). Original-coordinate axis labels and ticks match the lower panel of Figure 6.

Method: https://hinsley.github.io/dynamical-systems/projects/chebyshev-entropy-method.html
Both one-sided germs at each critical point; signed itinerary differences; first-column kneading minor divided formally by (1-t), truncated at degree 64 (and 32 for sensitivity). Entropy estimate is minus the natural logarithm of the smallest positive real root in (0,1); no root means a zero finite-order estimate. This region has invariant interval [-1,1], so escape continuation is unnecessary. Uses Float64 itineraries without fast-math or fused operations, exact integer polynomial coefficients, exact-integer initial Bernstein numerators, long-double subdivision/bracketing, and Arb root isolation for ambiguous cases. This is a finite-order estimate, not a certified bound on the true entropy.

`entropy_boundary.cpp` and `entropy_boundary.so` compute the coefficients and roots. `finish_entropy.py` applies Arb fallbacks, validates sampled production results independently, and saves all raw estimates and masks. `entropy-reference-polynomials.txt` contains 109 cases whose entire coefficient arrays matched Kneading.jl exactly, including critical hits. `entropy-validation.json` records numerical validation, sensitivity and fill counts. The Chebyshev endpoint recovers log(3); the superstable period-two example gives zero.

The website masks cells with |h64-h32|>0.01. In this new figure those cells, plus any unresolved cells, receive DISPLAY value 0.01 as requested. 85,155 of 9,000,000 cells are filled; none remain unresolved after Arb. Actual computed zeros are preserved. Raw h64/h32, failure/sensitivity/fill masks and display values are retained in `entropy_boundary_log_3000.npz`. The fill is not interpreted as computed positive entropy.

Colors: linear black to #f2c94c across [0,log(3)]. Red/blue contours are exactly the cached 19 iterates per critical orbit from Figure 6: opacity n^(-1.6), linewidth 2*n^(-1.1), later iterates drawn first.

Fast replotting: `render_HTOP_boundary.jl` renders the heatmap and axes using CairoMakie. `render_HTOP_contour_overlay.jl` renders the transparent contour layer ONCE. `compose_HTOP.py` overlays the cached 2007-square contour image in the identical transformed-coordinate rectangle, avoiding repeat contour calculations or drawing. Output: `HTOP_boundary_log_kneading.png`. Re-run only the base renderer and compositor for label/colorbar changes. Preserve `HTOP_contours_opacity16_width11.png` for reuse.

Replot command: `./replot_HTOP.sh` (uses cached entropy and contour layers).

Kneading.jl reference commit: f2d9add0cfc3651fc38dfee928f08ef240d9e89f.

Validation: 600 independent Arb comparisons, maximum absolute entropy difference 4.7573e-13; 871 distinct ambiguous polynomials resolved by Arb. Output label in the paper: `fig:cubic_htop_zoom`.
