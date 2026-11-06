#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
test -s HTOP_contours_opacity16_width11.png
julia --project=. render_HTOP_boundary.jl
uv run --with pillow python compose_HTOP.py
