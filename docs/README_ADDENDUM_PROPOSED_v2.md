# README追記案（公開前にREADME.mdへ統合）

## Archived 36-case campaign and derived model tables

The original production campaign covered four orbital radii (9,376, 12,000,
15,000, 20,000 km) and nine Mars-centric asymptotic projectile speeds
(0.5, 1, 2, 3, 4, 7.95298, 20, 40, and 60 km/s). The archived manifests record
10,000,000 source proposals per case, a 20-km injection sphere, an Ernst v004
shape model and enabled Martian and Phobos gravitational terms.

For inspection without distributing raw collision events, the repository includes:

- `data/derived/manifest_36cases_portable_v2.csv` — all 36 run configurations
  and provenance (absolute local file paths removed).
- `data/derived/sector_summary_36cases_v2.csv` — 36-case crater-production and
  flux summary data.
- `data/derived/model_global_four_sector_36cases_v2.csv` — 144 global
  (36 cases × 4 quadrants) model-sector rows.
- `data/derived/model_288_reference.csv` — 288 rows including both global and
  outside-Stickney25 masks.
- `workflow/as_used/` — unchanged copies of archived campaign/analysis scripts.
- `docs/FINAL36_CODE_AUDIT_JP_v2.md` — detailed data and provenance audit.

Validate the distributed summaries with:

```bash
python3 scripts/audit_final36_aggregates_v2.py
```

The command verifies CSV consistency, **not** reruns of the numerical
trajectories. It also does not regenerate the full set of manuscript figures.
Raw event catalogs, the external shape model, and licensed observed crater
coordinates remain outside this staging repository. Consult the release
checklist and data-access requirements before claiming full reproducibility.

**Status: PRE-RELEASE STAGING COPY.** This software and supplementary tables
must be verified against the final submitted manuscript before public release.
