"""Stack CairoMakie panels without resampling the ZILESmin raster."""
from pathlib import Path
from PIL import Image

root = Path(__file__).resolve().parent
upper = Image.open(root / "Figure6_orbit_upper_panel.png").convert("RGB")
# Keep a single shared u label beneath the lower panel. The upper panel's
# tick labels end above this crop; only its unused bottom label area is removed.
upper = upper.crop((0, 0, upper.width, 750))
lower = Image.open(root / "ZILESmin_boundary_log_tapered_opacity09.png").convert("RGB")
assert upper.width == lower.width == 3000
combined = Image.new("RGB", (lower.width, upper.height + lower.height), "white")
combined.paste(upper, (0, 0))
combined.paste(lower, (0, upper.height))
combined.save(root / "Figure6_ZILESmin_with_orbit.png")
print("Saved", combined.size)
