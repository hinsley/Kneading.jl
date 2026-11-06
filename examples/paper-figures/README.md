# Paper figure source archive

Source backup for the figures developed in the paper-editing chat through
2026-11-05. The computation, plotting, composition, validation, and retained
style variants are included, not just the final renderers. Original working
scripts are preserved byte-for-byte; `SOURCE-MANIFEST.json` records their
SHA-256 hashes and the small supporting input files.

## Figure index

All paths below are in `chebyshev-2026-11-05/` unless noted otherwise.
Figure numbers describe the paper on the snapshot date; use the LaTeX labels
as the durable identifiers.

| Figure or variant | Calculation | Rendering and assembly |
| --- | --- | --- |
| ZILESmax, original linear axes | `compute.jl` | `render_ZILES_max.jl`; original `render.jl` |
| ZILESmin, linear axes | Same cached critical-orbit exponents | `render_ZILES_min.jl` |
| ZILESmin, logarithmic axes near zero | `compute_ZILESmin_log.jl` | `render_ZILESmin_log.jl` |
| ZILESmin, logarithmic axes near u=1, v=-1 | `compute_ZILESmin_boundary_log.jl` | `render_ZILESmin_boundary_log*.jl`; final style is `render_ZILESmin_boundary_log_tapered_opacity09.jl` |
| Orbit diagram at v=-1/2, linear axes | `orbit_diagram_v_minus_half.jl` | Same script |
| Orbit diagram, logarithmic u and 10x retained iterates | `orbit_diagram_v_minus_half_log_10x.jl` | Same script |
| Figure 6 (`fig:ZILESmin`), orbit panel above Lyapunov/kneading panel | Boundary-log and orbit calculations above | `render_Figure6_orbit_upper_panel.jl`, `compose_Figure6.py`; cached-layout path: `render_Figure6_orbit_axis_frame.jl`, `refresh_Figure6_layout.py` |
| Figure 7 (`fig:cubic_htop_zoom`), refined topological entropy | `entropy_boundary.cpp`, `finish_entropy.py`, refinement/periodic-tail/symmetry scripts | `render_HTOP_boundary.jl`, `render_HTOP_contour_overlay.jl`, `compose_HTOP.py`, `replot_HTOP.sh` |
| Figure 8 (`fig:cubic_critical_cobweb`), 5,000-iterate critical cobweb | `render_cubic_critical_cobweb.jl`, including 256-bit orbit computation | Same script; use its PNG in the paper |
| Earlier contour-rule explanatory illustration | `../contour-rule-explanation/draw_contour_example.py` | Same script; this was a discussion illustration, not a final paper figure |

The earlier contour diagram examples already live in this repository's
`examples/` directory. The interactive cubic website and its symbolic orbit
finder are separately versioned in
[hinsley.github.io](https://github.com/hinsley/hinsley.github.io/tree/main/dynamical-systems/projects).

## Environment and first run

The archived Julia environment uses Julia 1.13.3 and pins package versions in
`Manifest.toml`. The `Kneading.jl` path dependency refers to commit
`f2d9add0cfc3651fc38dfee928f08ef240d9e89f`. Restore that dependency without
changing the surrounding repository checkout:

```sh
cd examples/paper-figures/chebyshev-2026-11-05
bash setup_environment.sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

The original Lyapunov scans require a CUDA-capable GPU and were calculated on
an RTX 4070 Ti. Rendering and the cobweb do not require GPU computation. C++
entropy helpers require GCC with OpenMP and Boost development headers:

```sh
bash build_entropy.sh
```

The cobweb is self-contained after environment setup:

```sh
julia --project=. render_cubic_critical_cobweb.jl
```

Its 2700 x 2700 PNG is the paper asset. A vector PDF is also exported, but
some PDF viewers thicken the densely overlaid sub-pixel lines; the PNG
preserves the intended thin, translucent appearance.

## Calculation and cache notes

Read `README.md` in the snapshot directory for the original linear scan,
`ZILES.md` for the subsequent variants, `Figure6-stacked.md` for exact panel
alignment and fast layout refreshes, and `HTOP.md` / `HTOP-refinement.md` for
the entropy computation and its limitations. Those documents are historical
working notes: later sections can supersede earlier descriptions.

For example, from the snapshot directory:

```sh
julia --project=. --threads=4 compute_ZILESmin_boundary_log.jl
julia --project=. --threads=4 render_ZILESmin_boundary_log_tapered_opacity09.jl
julia --project=. --threads=8 orbit_diagram_v_minus_half_log_10x.jl
julia --project=. --threads=8 render_Figure6_orbit_upper_panel.jl
python compose_Figure6.py
```

The initial entropy sweep uses the archived exact grid input:

```sh
./entropy_boundary 3000
python finish_entropy.py
```

The final entropy image additionally uses the adaptive refinement scripts
and assembly steps described in `HTOP-refinement.md`. These are archived
research scripts with staged cache dependencies, not a newly tested
one-command reconstruction pipeline. The cached-layout Figure 6 variant
also uses an earlier rendered orbit interior; its recipe is documented in
`Figure6-stacked.md`. Do not confuse an initial entropy scan with the final
refined data, or a full orbit redraw with the pixel-preserving layout edits.

Large computed arrays and rendered caches remain in the original local
working directory (`/home/carter/Documents/chebyshev-figure`). They are not
included in this **source-code** backup. `LOCAL-CACHE-INVENTORY.json` lists
their names and sizes. Small exact grid and refinement-selection inputs are
included here, along with the saved validation reports. Generated arrays,
binary libraries, downloaded Boost headers, and screenshots are excluded
from Git. No manuscript source snapshots are included.

The explanatory contour-rule script is an unchanged historical source: it
uses DejaVu Sans from the Linux system font directory and an absolute output
path. Adjust its final `im.save(...)` argument to run it elsewhere.

## Verify the backup

From this directory:

```sh
python3 verify_archive.py
```

This checks every archived file against its original working-file hash. It
does not rerun the expensive numerical scans.
