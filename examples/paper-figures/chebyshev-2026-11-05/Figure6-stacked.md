# Figure 6: orbit diagram above ZILESmin

Overleaf project 6abfd08d4a9a262f96c0c49e now includes
`figures/Figure6_ZILESmin_with_orbit.png` at full text width, under
the existing `fig:ZILESmin` label. Compiled and verified on page 10.
The caption was minimally extended to identify the upper orbit panel
and lower Lyapunov/kneading panel; surrounding prose was preserved.

Reproduce:
```sh
julia --project=. --threads=8 render_Figure6_orbit_upper_panel.jl
python3 compose_Figure6.py
```

The upper panel reads the existing 6,001 x 10,240 Float64 orbit data,
uses black markers of size 0.8 at opacity 0.025 (reduced from 0.10),
and a 669 x 185 logical-pixel plot area. Its horizontal bounds align
to the lower panel within pixel rounding, using the same transformed
coordinate -log(1-u) and parameter ticks. It has a linear state axis.

The lower panel reuses `ZILESmin_boundary_log_tapered_opacity09.png`
without resampling. Python/Pillow removes the outer margins and joins
the data-area boundaries directly, preserving the full v=0 tick label.
The orbit panel shows state-axis ticks -0.5, 0, 0.5, and 1 and label x on the
left. It hides its horizontal ticks and label; the lower panel provides
the shared parameter-axis ticks. The colorbar label is translated 105 pixels
left by compose_Figure6.py, closer to the bar, without resampling.
Final image: 2910 x 2901 pixels. The single shared u axis is at the bottom.
The former gapped, opacity-0.10 image and renderer are saved separately.

Original ZILESmin image and all prior orbit plots, sources, and data remain
available. To revert the paper's figure, restore its includegraphics path
to `figures/ZILESmin.png` and remove the new upper/lower caption wording.

For this layout-only update, render_Figure6_orbit_right_axis_frame.jl
rendered the new CairoMakie axes; the unchanged data interior
(443,47,2447,599) was copied from Figure6_orbit_upper_panel_no_ticks.png
into that frame, preserving every data pixel. The full upper renderer also
contains the updated right-axis settings for future regeneration.

## Fast layout updates

Run `python3 refresh_Figure6_layout.py` to redraw only the CairoMakie axes
and reuse cached data pixels. Edit `render_Figure6_orbit_axis_frame.jl` for
axis formatting, keeping the full upper renderer consistent. This avoids
loading or redrawing 61,450,240 orbit points for label/tick changes.
Raw orbit, Lyapunov, and contour data remain on disk. The cached image is
valid for the existing data limits, viewport, resolution, and marker style;
changes to those require refreshing the relevant rendered data panel.

The upper state range is now [-0.8, 1.3], adding balanced internal padding
above and below the orbit data, so the -0.5 tick remains centered on its value without colliding with the
lower v=0 label. The full renderer uses these limits. Fast layout refresh
affinely rescales the cached data interior from its original [-0.52,1.02]
range, with Lanczos interpolation; no orbit calculations are repeated.
The original bottom spine is excluded from the data crop. The full renderer
is retained for exact marker-size regeneration if required.

The lower v=0 tick label is shifted downward by 30 output pixels (10
logical pixels) in compose_Figure6.py; its tick and the data are unchanged.

Colorbar refinement: moved lambda_min another 45 pixels left (105 total),
then trimmed 90 pixels from the right edge, retaining about 60 pixels of
white margin beyond the label. Data pixels are unchanged by this step.
