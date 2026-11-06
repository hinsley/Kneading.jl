"""Stack CairoMakie panels without resampling the ZILESmin raster."""
from pathlib import Path
from PIL import Image

root = Path(__file__).resolve().parent
upper = Image.open(root / "Figure6_orbit_upper_panel.png").convert("RGB")
# Upper axis bottom is y=600 px; lower axis top is y=72 px.
# Remove outside margins so these data-area boundaries meet exactly.
# Horizontal decorations are disabled in the CairoMakie renderer.
upper = upper.crop((0, 0, upper.width, 600))
lower = Image.open(root / "ZILESmin_boundary_log_tapered_opacity09.png").convert("RGB")
# Move only the colorbar label 105 pixels closer; keep the original raster.
label_box = (2820, 950, 3000, 1200)
label = lower.crop(label_box)
lower.paste("white", label_box)
lower.paste(label, (2715, 950))
# Nudge only the v=0 label down by 30 output pixels (10 logical pixels).
# Its tick and all data stay at their original coordinates.
zero_box = (360, 35, 412, 108)
zero_label = lower.crop(zero_box)
lower.paste("white", zero_box)
lower.paste(zero_label, (360, 65))
assert upper.width == lower.width == 3000
combined = Image.new("RGB", (lower.width, upper.height + lower.height - 72), "white")
combined.paste(upper, (0, 0))
combined.paste(lower.crop((0,72,lower.width,lower.height)), (0, upper.height))
# Preserve the full v=0 label straddling the touching boundary, while
# keeping the upper panel's data intact underneath the lower top margin.
combined.paste(lower.crop((0,35,438,72)), (0,upper.height-72+35))
combined.paste(lower.crop((2451,0,3000,72)), (2451,upper.height-72))
# Trim the right margin, retaining 60 pixels beyond the moved label.
combined = combined.crop((0, 0, 2910, combined.height))
combined.save(root / "Figure6_ZILESmin_with_orbit.png")
print("Saved", combined.size)
