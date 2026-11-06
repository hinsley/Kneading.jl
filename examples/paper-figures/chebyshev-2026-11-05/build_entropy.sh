#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
# Requires GCC/OpenMP and Boost development headers. Never enable fast-math:
# changes in critical itineraries can alter the kneading coefficients.
compiler=${CXX:-g++}
flags=(-std=c++17 -O3 -ffp-contract=off -fopenmp)
"$compiler" "${flags[@]}" entropy_boundary.cpp -o entropy_boundary
"$compiler" "${flags[@]}" -shared -fPIC entropy_bridge.cpp -o entropy_boundary.so
"$compiler" "${flags[@]}" -shared -fPIC entropy_refine.cpp -o entropy_refine.so
"$compiler" "${flags[@]}" -shared -fPIC periodic_tail.cpp -o periodic_tail.so
"$compiler" "${flags[@]}" -shared -fPIC symmetric_entropy.cpp -o symmetric_entropy.so
