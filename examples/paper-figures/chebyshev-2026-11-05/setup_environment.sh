#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
kneading_revision=f2d9add0cfc3651fc38dfee928f08ef240d9e89f
if [[ ! -d Kneading.jl ]]; then
    git clone https://github.com/hinsley/Kneading.jl.git Kneading.jl
    git -C Kneading.jl checkout --detach "$kneading_revision"
fi
if [[ $(git -C Kneading.jl rev-parse HEAD) != "$kneading_revision" ]] ||
   [[ -n $(git -C Kneading.jl status --porcelain) ]]; then
    echo 'Existing Kneading.jl dependency differs from the archived revision; leaving it untouched.' >&2
    exit 1
fi
julia --project=. -e 'using Pkg; Pkg.instantiate()'
