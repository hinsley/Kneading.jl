from pathlib import Path
from PIL import Image
root=Path(__file__).resolve().parent
base=Image.open(root/'HTOP_boundary_base.png').convert('RGBA')
contours=Image.open(root/'HTOP_contours_opacity16_width11.png').convert('RGBA')
assert base.size==(2910,2373)
assert contours.size==(2007,2007)
# Shared transformed-coordinate rectangle: x=[147,816], y=[98,767] points at 3 px/point.
base.alpha_composite(contours,(441,72))
base.convert('RGB').save(root/'HTOP_boundary_log_kneading.png')
base.thumbnail((1000,1000))
base.convert('RGB').save(root/'HTOP_preview.png')
