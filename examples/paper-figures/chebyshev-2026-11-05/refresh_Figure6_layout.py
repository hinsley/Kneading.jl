"""Refresh axes from CairoMakie while reusing unchanged orbit/heatmap pixels.
Uses the cached orbit raster with an affine vertical mapping to the new limits.
Data, horizontal limits, and output resolution stay unchanged.
"""
from pathlib import Path
import subprocess
import sys
from PIL import Image
root = Path(__file__).resolve().parent
if "--reuse-frame" not in sys.argv:
    subprocess.run(["julia", "--project=.", "render_Figure6_orbit_axis_frame.jl"], cwd=root, check=True)
frame = Image.open(root / "Figure6_orbit_axis_frame.png").convert("RGB")
cached = Image.open(root / "Figure6_orbit_upper_panel_no_ticks.png").convert("RGB")
box = (443, 47, 2447, 598)
# Original raster covered x-state [-.52, 1.02]; the new axis covers [-.8, 1.3], with balanced internal padding.
# Affine mapping preserves the original coordinates and adds blank space below.
# Include the original subpixel boundary positions in the mapping.
old_top, old_bottom = 45.0, 600.0
old_span, new_span = 1.54, 2.1
scale = old_span / new_span
top_padding = (1.3 - 1.02) / new_span * (old_bottom - old_top)
new_top = old_top + top_padding + (box[1] - old_top) * scale
new_bottom = old_top + top_padding + (box[3] - old_top) * scale
height = round(new_bottom) - round(new_top)
interior = cached.crop(box).resize((box[2]-box[0], height), Image.Resampling.LANCZOS)
frame.paste(interior, (box[0], round(new_top)))
frame.save(root / "Figure6_orbit_upper_panel.png")
subprocess.run(["python3", "compose_Figure6.py"], cwd=root, check=True)
