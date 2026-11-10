# As-used production and analysis scripts (archival copies)

These six files were copied unchanged from `05_production_v004_pg1/tools/` and
`05_production_v004_pg1/scripts/` in the research workspace. They document the
specific 36-case campaign and map/masking workflow. They are **not** a standalone
runnable checkout as published here. In particular:

- `run_lwss_campaign_v1.py` resolves `BASE` from its original parent directory,
  expects `src/` and the v004 OBJ model, and writes simulation outputs there.
- `make_v004_pg1_impact_maps_v1.py` requires full raw 26-column ON/OFF/THROUGH
  event files and the matching OBJ; neither is redistributable by assumption.
- `rebuild_exact_four_sector_from_raw_v2.py` also requires a harmonized,
  licensed observed crater-coordinate catalog (not in this repository).
- Several as-run scripts refer to exact workspace subdirectories that are
  **not** guaranteed to exist in this GitHub checkout.

To run these scripts, use an appropriately reconstructed local
`05_production_v004_pg1/` tree with the original `src/`, input shape and
catalogs. Do not edit the archive in place; use a versioned reproduction wrapper
when a fully portable workflow is developed.

The C++ source used for production is in the repository's existing top-level
`src/` directory. It is not duplicated here to avoid divergent copies.
