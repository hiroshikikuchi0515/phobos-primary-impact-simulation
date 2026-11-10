#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
MODE="${1:-help}"
BASE=figures_for_paper_v004_pg1_v1
case "$MODE" in
  inventory)
    python3 tools/index_existing_figures_v1.py --root .. --out "$BASE/00_inventory"
    echo "Open: $BASE/00_inventory/START_HERE_v1.txt"
    echo "Preview: $BASE/00_inventory/existing_figures_contact_v1.jpg"
    ;;
  check)
    python3 tools/make_v004_pg1_impact_maps_v1.py --dry-run
    ;;
  maps)
    python3 tools/make_v004_pg1_impact_maps_v1.py
    echo "New maps: $BASE/02_pg1_maps"
    ;;
  collect)
    if [[ ! -f "$BASE/02_pg1_maps/figure_provenance_v1.json" ]]; then
      echo "New map output missing; run 'maps' first." >&2; exit 2
    fi
    out="$BASE/paper_figures_review_v1.tar.gz"
    if [[ -e "$out" ]];then echo "Refusing to overwrite $out" >&2;exit 3;fi
    members=(02_pg1_maps)
    [[ ! -d "$BASE/00_inventory" ]] || members+=(00_inventory)
    tar -czf "$out" -C "$BASE" "${members[@]}"
    ls -lh "$out"
    ;;
  *)
    echo 'Usage: bash tools/run_paper_figures_v1.sh {inventory|check|maps|collect}'
    echo '  inventory: find historical figures, plot scripts, and make contact sheet'
    echo '  check: verify that PG ON/v004 26-column event catalogs exist (no reading)'
    echo '  maps: generate Main & Supplementary maps from two 7.953-km/s cases'
    echo '  collect: bundle small figures, PNG/PDF, gridded CSV, and provenance for ChatGPT'
    exit 2;;
esac
