# Entropy edge refinement (2026-11-05)

The original image filled 85,155 degree-32/64-sensitive cells at 0.01. This refinement targets those cells and tightens the sensitivity threshold to 0.001, giving 173,916 cells. Eight additional zero outliers on the symmetry diagonal are checked independently.

The map, grid, natural-log entropy convention and kneading-determinant construction are unchanged. The critical itineraries and determinant truncations are extended adaptively: 256, 512, 1024, 2048, 4096, 8192 and longer for the difficult cases. The comparison tolerance is 0.001, a numerical diagnostic, not a certified bound on the true entropy.

Two numerical issues require particular treatment:

* For numerically eventually periodic kneading rows, both Float64 and independent long-double critical itineraries are run for 32,768 steps. Each signed-symbol tail must repeat at least eight times in the final 8192 entries, and both precisions must agree on the entire prefix, start and period. The periodic series are then summed rationally, the exact integer determinant numerator is simplified against its denominator, and Arb isolates its roots. This is the same kneading determinant with a numerically verified periodic continuation; the inferred continuation is not a rigorous theorem about the exact real orbit. Thirty-two sampled rational determinants reproduced every coefficient of the independently computed degree-512 series exactly.
* On the symmetry line $v=-u$, the kneading rows have the exact form $(A,B,C)$ and $(-C,-B,-A)$, so $D=B(C-A)/(1-t)$. The two factors are truncated separately at degrees 256 and 512 and their roots isolated with Arb. This algebraically identical determinant representation avoids losing a nearly repeated root when the product is expanded and truncated. The degree comparison uses the largest entropy from either factor. Seven sampled factorizations reproduced all 513 coefficients of the expanded determinant exactly. Adjacent-degree expanded-polynomial checks were also explored and cached, but are not used for the final symmetric estimates because both adjacent truncations can miss the limiting root.

The raw intermediate arrays, root results, periodic-tail information, masks and actual per-cell orders are retained. Unresolved cells are kept explicitly flagged rather than described as converged. The original degree-64 image and raw data remain available for comparison.

The existing CairoMakie base renderer and cached transparent contour layer are reused. No contour geometry or styling is changed.

Final result: 173,924 cells recomputed (including eight additional zero outliers); 173,922 pass the tightened 0.001 agreement check. Two cells remain flagged after 32,768 terms, at (u-index,v-index)=(271,1827) and (1172,2728), zero-based, and retain the display fallback 0.01. Their last finite estimates, about 0.00187973, are preserved in the raw refined array. There is no smoothing or spatial interpolation of entropy values.

Final cache: `entropy_boundary_refined_final_3000.npz`. `method` records direct extended truncation (1), rational periodic tails (2), or factored symmetric determinant (3); the eight extra points are stored separately. `entropy-refinement-final.json` gives counts. Replotting uses the same `entropy_display_3000.bin`, CairoMakie renderer, and transparent contour cache as before.

## Boundary-zero audit

Visual QA exposed short-order zero/zero agreements that the original difference mask could not flag. A second pass recomputed 15,613 zero boundary cells (zeros adjacent to values above 0.01, plus nonperiodic zeros from the first refinement). Numerically periodic tails were checked independently as above; other cells required longer truncations, with zero/zero agreement not accepted before degree 4096. All cells in this second pass ultimately passed the 0.001 diagnostic. Four additional symmetric cases used the exact determinant factorization; one nonsymmetric parameter and its reflected conjugate required 32,768 terms. The final image retains only the two unresolved cells identified in the first pass.

`merge_boundary_zeros.py` assembles the final raw/display fields. It preserves the first pass in `entropy_boundary_refined_primary_3000.npz` and stores boundary-zero IDs and orders separately in the final NPZ. `entropy-refinement-final.json` is the current aggregate summary. The old root-order comparison experiments are retained for audit, not used as the final symmetry estimator.
