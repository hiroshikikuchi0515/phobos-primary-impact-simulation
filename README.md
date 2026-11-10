# Phobos primary-impact simulations

Research software associated with a manuscript in preparation for *Icarus* on the spatial distribution of primary impacts on Phobos, including Martian gravitational deflection, shielding by Mars, and the irregular shape of Phobos.

**Status: PRE-RELEASE STAGING COPY.** The files in this starter directory are supplied to assist repository setup; they have not yet been vetted as the complete, self-contained research release. Do not label the repository or these files as a finished reproduction package until the checks in `docs/RELEASE_CHECKLIST.md` are complete.

## What is here

- `src/phobos_impact_source_sphere_production_v5.cpp`: finite source-sphere Monte Carlo implementation (research source; source-sphere weighting, paired Mars shielding ON/OFF).
- `src/phobos_impact_production_v3.cpp`: base source included by the v5 program. It is needed to compile v5.
- `scripts/build.sh`: build command for GCC/Clang with C++17 and OpenMP.
- `analysis/`: selected independent benchmark and validation/rotation-analysis scripts. These are archival scripts, not yet a fully packaged workflow.
- `data/derived/model_288_reference.csv`: model comparison table from project results, for provenance/inspection. Verify the final table and its column definitions before releasing.
- `data/external/`: explanation of separately obtained input assets.
- `docs/`: reproducibility, publication, and release notes.
- `CITATION.cff`: citation metadata for the code, initially without a DOI.

## Build

Prerequisites: C++17 compiler with OpenMP support. For example, on Ubuntu:

```bash
bash scripts/build.sh
```

To see all command-line options:

```bash
./build/phobos_source_sphere_v5 --help
```

## Model scope

The production source uses a uniform-area injection sphere weighted by the inward flux of the Mars-only local Liouville distribution. For the production case, Mars gravity and point-mass Phobos gravity are enabled and the Phobos surface is represented with an external Ernst v004 OBJ shape model. The validation cases include a sphere target with Phobos gravity disabled. Detailed assumptions, test conditions, and numerical uncertainties must be described in conjunction with the associated article.

This repository does not automatically distribute the third-party crater catalogs or the Phobos shape model. Their sources, permission terms, and required versions should be documented separately.

## Reproducibility status

The code archive should eventually include precise input versions, model parameters, run manifests, numerical outputs needed to regenerate article figures, and working end-to-end analysis instructions. A code repository alone is **not** equivalent to release of full observation data or mesh inputs.

## Citing the code

Use GitHub's **Cite this repository** control, populated from `CITATION.cff`. A versioned DOI can be added later by archiving a GitHub Release with Zenodo. Do not invent a DOI in a manuscript before one is issued.

## License

**To be selected by the author after checking ownership and third-party dependencies.** No license is granted by the presence of the current files. Do not assume materials taken from other teams or missions may be redistributed.
