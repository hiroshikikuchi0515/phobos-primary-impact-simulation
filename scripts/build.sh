#!/usr/bin/env bash
set -euo pipefail
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$project_root/build"
g++ -O3 -std=c++17 -fopenmp -Wall -Wextra \
  "$project_root/src/phobos_impact_source_sphere_production_v5.cpp" \
  -o "$project_root/build/phobos_source_sphere_v5"
echo "Compiled: $project_root/build/phobos_source_sphere_v5"
"$project_root/build/phobos_source_sphere_v5" --help
