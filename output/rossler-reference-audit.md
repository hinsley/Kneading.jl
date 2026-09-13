# Rössler reference audit

The independent regression check uses the previously generated origin-fixed Rössler scan, not newly generated expected values.

Source: `/Users/carterhinsley/Documents/Dev/MultimodalMaps.jl/flow_folding/results/rossler_y_minima_critical_orbit_scan_256/coarse_scan.tsv`.

Source SHA-256: `fa30cfeb7b97223a8447d2115ca1115675ed568ccd128b162dd4f729d7ac84d0`.

## Reference protocol

- System: `dx = -y - z`, `dy = x + a*y`, `dz = 0.3*x + z*(x-c)`.
- Grid: 256 values of `c` from 2 to 7, and 256 values of `a` from 0.30 to 0.55.
- Critical point: fourth accepted local minimum of `y` from a saddle-focus launch near the origin. No seed-radius refinement was used in this historical run.
- Critical-point solve: first-order forward sensitivities, event-time correction, capped Newton steps with a central finite-difference slope; scan and bisection fallback. The launch guard is half the linearized spiral period. The criticality tolerance is `1e-6`, and the ODE absolute and relative tolerances are `1e-9`.
- Kneading orbit: start at the critical point, include its sign at time zero, then capture seven more local minima of `y`, without discarding a transient. Fixed-step RK4 uses `dt = 0.05` and linearly interpolated events; the tangent is projected orthogonal to the flow and normalized at every step and event.
- Stop after eight raw signs, time 2000, or a state component exceeding `1e6`. Positive tangent components encode as one. Adjacent equal signs encode positive monotonicity, producing seven transition signs from eight raw signs.

The reference contains 65,536 successful critical-point initializations, 48,420 complete words, and 17,116 incomplete orbits. The historical status calls every incomplete orbit `orbit_max_time`, even when a state bound ended integration earlier. The new API reports the actual stop condition, so comparison uses the completeness mask rather than requiring identical status strings.

## Independent 25-point check

`test/fixtures/rossler_reference.tsv` contains exact historical rows at the Cartesian product of grid indices 1, 64, 128, 192, and 256. The verifier performs a fresh initialization at every point without using the stored radius, state, or tangent as an initial guess.

Run:

```sh
julia --project=. --startup-file=no examples/verify_rossler_reference.jl
```

Verified results:

- 355 assertions passed.
- 25/25 raw words or incomplete prefixes matched exactly.
- 25/25 completeness masks matched.
- Maximum critical-state error: `5.86875e-7` in Euclidean norm.
- Maximum oriented unit-tangent error: `1.14075e-9`.
- Five second-order-sensitivity Newton corrections, initialized with a perturbation of the newly computed seed coordinate, produced the same words and masks. Maximum state difference from finite-difference Newton: `6.52811e-7`.
- A separate integration-only check starting from the stored critical states and tangents reproduced all 25 historical words or prefixes. This isolates the orbit integration from the independent initialization check.

These checks reproduce the historical fixed-index numerical protocol. They do not substitute for the default seed-radius refinement test, nor establish that the historical diagram is invariant under time-step refinement.

## Critical-return orientation sensitivity

At `a = 0.3`, `c = 5.5`, the default refined initializer accepts `M = 4` after checking the `M = 5` candidate. The raw and transition words were checked at three tolerance levels and with adaptive integration and two RK4 step sizes. Sign thresholds were set to zero only for this diagnostic, to retain the first eight signs for comparison.

| Tolerance level | Criticality tolerance | Initialization absolute/relative tolerances | Orbit absolute/relative tolerances | Adaptive raw word | Adaptive transition word |
| --- | --- | --- | --- | --- | --- |
| Default | `1e-8` | `1e-11`, `1e-11` | `1e-10`, `1e-9` | `01000110` | `0011010` |
| Tight | `1e-10` | `1e-12`, `1e-12` | `1e-12`, `1e-11` | `00111001` | `1011010` |
| Tighter | `1e-12` | `1e-13`, `1e-13` | `1e-13`, `1e-12` | `00111001` | `1011010` |

For all three initialization tolerance levels, RK4 with either `dt = 0.05` or `dt = 0.025` produced raw word `00111001` and transition word `1011010`. Across all nine cases, `transition_word[2:end]` was exactly `011010`.

The default adaptive case complements every raw sign after the initial critical event relative to the other eight cases; only the first transition differs. The first three normalized tangent components have absolute values above `0.9999` in every case. A small-component sign guard therefore does not reveal this orientation sensitivity: normalization can hide the strong contraction of the unnormalized tangent near the critical return.

This is evidence of a stable six-transition suffix at this parameter point, not proof that all parameter points have the same property. A historical fixed-step reproduction and numerical convergence of every symbol are distinct checks.

## Full-grid comparison

The completed new scan is `output/rossler-flow-kneading.tsv`.

New scan SHA-256: `3a0e2df11605cd0db8a7af5cf914c5919f19b92709303ad547fd907a9e59bda4`.

It was generated through `examples/rossler_kneading_diagram.jl` using the public `SaddleFocusInitializer`, `scan_flow_kneading`, and `write_flow_scan` interfaces, with fresh numerical initialization and no stored reference states or words used as scan input.

The completed comparison reports:

- 65,536/65,536 parameter coordinates matched.
- 48,420 complete words in both scans, with no new initialization failures.
- 65,536/65,536 completeness masks matched.
- All 17,116 incomplete prefixes matched exactly.
- 48,406/48,420 complete raw words matched exactly; 14 differed.
- Combining complete words and incomplete prefixes, 65,522/65,536 grid entries matched exactly, or approximately 99.979%.
- The seven-transition words likewise differ at those same 14 complete grid points.
- Maximum full critical-state difference: `5.93608277e-7`; RMS difference: `1.77839576e-7`.

All 14 differing grid cells have a differently labeled immediate neighbor in the historical scan, placing them on resolved symbolic-region boundaries. At 12 of these cells, the new word occurs in an immediate historical neighbor; at the other two it does not. This observation does not identify either finite-precision word as mathematically exact.

At the differing points, critical-state differences range from `5.12098e-9` to `1.94998e-8`. The historical criticality residual magnitudes range from `3.27572e-10` to `1.21451e-9`, while the new residual magnitudes range from `1.57647e-8` to `6.27678e-8`; all satisfy the configured `1e-6` residual tolerance.

For an independent isolation check, the new library integrated each of the 14 historical critical states and section tangents with the same fixed-step orbit options. It reproduced all 14 historical words exactly. At these points, this separates differences in the freshly computed initializations from the orbit integrator and event encoding. No words, masks, or initial states were patched into the new scan to force agreement.

Exact differing coordinates, words, residuals, and state errors are recorded in [rossler-reference-mismatches.tsv](rossler-reference-mismatches.tsv).

Repeat the comparison:

```sh
julia --project=. --startup-file=no examples/verify_rossler_reference.jl --full output/rossler-flow-kneading.tsv /Users/carterhinsley/Documents/Dev/MultimodalMaps.jl/flow_folding/results/rossler_y_minima_critical_orbit_scan_256/coarse_scan.tsv
```

The comparator reports exact word/prefix agreement, completeness-mask agreement, bitwise disagreements, and critical-state errors. It does not silently treat numerical disagreements as matches.
