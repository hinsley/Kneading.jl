# Orbit diagram at v = -1/2, log-scaled near u = 1

This variant preserves the original linear-axis plot and its source/data.
It uses the same coordinate transformation and cutoff as ZILESmin:
horizontal coordinate -log(1-u), with 0 <= u <= 0.99.
Tick labels remain actual parameter values: 0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99.

6,001 parameter samples are uniformly spaced in the transformed coordinate.
At each sample, start from the critical point x=-1/2, discard 32,768 iterates,
then retain 10,240 consecutive iterates (10 times the original 1,024).
The other critical point maps to -1/2 in one step on this parameter slice,
so the same orbit tail covers both critical points.
All 61,450,240 points are plotted with the original marker size and opacity.

The computation uses Float64 direct scalar iteration; rendering uses CairoMakie.
The PNG is 3000 x 1800. JLD2 stores the original u values, transformed values,
full state matrix (iterate x parameter), and computation settings.

Run from this directory:
```sh
julia --project=. --threads=8 orbit_diagram_v_minus_half_log_10x.jl
```
