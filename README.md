# Phobos primary-impact simulations

Source code and selected derived research data supporting the manuscript:

> Kikuchi, H., *Primary-impact asymmetries on Phobos: orbital-distance and speed dependence, Mars shielding, and the limitations of crater-record inference*, prepared for *Icarus* (2026).

**Status:** Research software and derived-data release candidate. A versioned GitHub release and Zenodo DOI have not yet been issued. This repository documents and supports verification of important manuscript results, but it is **not a self-contained, end-to-end reproduction environment for every figure**.

## Scientific scope

The production campaign contains 36 Monte Carlo cases: four orbital radii (9,376, 12,000, 15,000 and 20,000 km) and nine isotropic Mars-centric asymptotic speeds (0.5, 1, 2, 3, 4, 7.95298, 20, 40 and 60 km/s), with 10 million proposals per case. It includes Mars gravity and shielding, Phobos point-mass gravity during trajectory integration, and the Ernst SPC v004 irregular Phobos shape model. The relative crater-production proxy uses q = 2.5 and mu = gamma = 0.55.

Spatial comparisons use four body-fixed quadrants (leading/Mars-facing LM, leading/anti-Mars LA, trailing/Mars-facing TM and trailing/anti-Mars TA). The Stickney25 mask excludes a 25-degree angular cap centered at east longitude -50 degrees and latitude +1.3 degrees. The two overlapping observational crater catalogs (Salamunićcar et al., 2014; Karachevtseva et al., 2014) are analyzed separately. Nominal goodness-of-fit p values are descriptive fixed-model tests and do not estimate an actual impactor population or reconstruct orbital history.

## Repository contents

- `src/phobos_impact_source_sphere_production_v5.cpp` and `src/phobos_impact_production_v3.cpp`: C++17/OpenMP Monte Carlo implementation, including paired shielding ON/OFF/THROUGH. Both sources are required to compile.
- `scripts/build.sh`: Linux build script.
- `data/derived/`: the 36-case run manifest and production summaries; 288 masked sector records; 144 global model-sector records; 144 model-catalog comparison tests; observational quadrant and Stickney counts; reported two-speed mixture thresholds; 18 OBJ source-radius trials; and 60 spherical benchmark runs. See [the table index](data/derived/README.md).
- `analysis/`, `figures/archived_v004/` and `workflow/as_used/`: selected analysis, historical figure-generation and as-used campaign scripts. Several require external data or the original project directory layout, and are **not** a portable workflow for every final manuscript figure.
- `docs/`: source-code provenance, audits, release checks and noted reproducibility limitations.
- `data/external/README.md`: external third-party inputs and their access considerations.

## Build

A Linux C++17 toolchain with OpenMP is required:

```bash
bash scripts/build.sh
./build/phobos_source_sphere_v5 --help
```

Building does not itself reproduce the manuscript campaign. Full event-level reruns require separately acquired shape data, the original numerical settings, the indicated random seeds and a compatible environment. The original campaign used eight OpenMP threads and static scheduling. Consult the [36-case manifest](data/derived/manifest_36cases_portable_v2.csv), [code audit](docs/FINAL36_CODE_AUDIT_JP_v2.md) and [workflow notes](workflow/as_used/README.md).

## Scope and limits of reproducibility

The public tables allow inspection and independent checking of the 36-case model summaries, sector shares, observed regional counts, exploratory nominal model-catalog comparisons, published rounded mixture thresholds and numerical validation statistics. The `sphere_6groups_summary_v2.csv` file contains **arithmetic means of per-run ratios**, while Supplementary Table S10 reports **ratios based on pooled weighted production rates**; small differences are expected.

The original 26-column ON/OFF/THROUGH impact-event catalogs are retained in the author's research archive and **are not included**. The Ernst v004 OBJ mesh and original third-party crater catalogs are **not redistributed**. This archive is not claimed to recreate every published figure from a clean checkout; see the [figure provenance audit](docs/FIGURE_PROVENANCE_v3_JP.md).

## External sources

- Ernst et al. (2023), *High-resolution shape models of Phobos and Deimos from stereophotoclinometry*, https://doi.org/10.1186/s40623-023-01814-7. The original mesh is not included; the run manifest records its filename and SHA-256.
- Salamunićcar et al. (2014), PH9224GT crater catalog, https://doi.org/10.1016/j.asr.2013.11.006.
- Karachevtseva et al. (2014), *The Phobos information system*, https://doi.org/10.1016/j.pss.2013.12.015.

Derived observational counts are provided for interpreting manuscript comparisons; they do not confer redistribution rights for original crater coordinates or catalog records.

## Citation and reuse

See [CITATION.cff](CITATION.cff) for citation metadata. Once release `v1.0.0` is archived on Zenodo, cite its **version-specific DOI** rather than only the moving `main` branch. The software release is versioned as `1.0.0` in `CITATION.cff` in preparation for the GitHub Release; the DOI is intentionally omitted until assigned.

**Licensing by material type:**

- **Author-owned software source** (including C++, Python and shell scripts): [MIT License](LICENSE), Copyright (c) 2026 Hiroshi Kikuchi. The MIT grant applies only to contributions for which the author can grant rights.
- **Original author-created derived numeric tables** in `data/derived/`: [Creative Commons Attribution 4.0 International](data/derived/LICENSE) (CC BY 4.0), to the extent of rights held by the author. The published data provide scientific summary values, not the original impact-event records or external inputs.
- **Third-party inputs**, including the original Ernst v004 shape model and the Salamunićcar and Karachevtseva crater catalogs: **not included and not relicensed**. Their citation and any original access/use conditions remain applicable. Derived counts from those catalogs do not grant rights in their underlying records.

Research conclusions and any third-party work remain attributable to their original sources. Availability under an open license does not mean the simulation is a complete reproduction environment for all article figures.

## Funding

JAXA SSF Program Japan Grant Number JPJXSSF24MX17002.
