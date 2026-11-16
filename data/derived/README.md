# Derived tables and provenance

These are **derived numeric summaries**, not the original 26-column event files, full crater catalogs or Ernst v004 OBJ mesh. Row counts exclude headers.

| File | Rows | Contents |
| --- | ---: | --- |
| `manifest_36cases_portable_v2.csv` | 36 | Run parameters, seeds and input mesh hash; referenced raw event files are not deposited. |
| `sector_summary_36cases_v2.csv` | 36 | Per-case global L/T and M/A, shielding and sector results. |
| `model_288_reference.csv` | 288 | 36 cases x four sectors x two masks (Global and Outside Stickney25), with model areas, impact counts and weighted crater-production shares. |
| `model_global_four_sector_36cases_v2.csv` | 144 | 36 cases x four sectors **of the model**, not 144 model-catalog tests. |
| `observed_four_quadrant_counts_v1.csv` | 4 | Derived LM/LA/TM/TA observed counts for two catalogs and two masks, transcribed from manuscript v31 Table 2. |
| `observed_Stickney25_disjoint_counts_v2.csv` | 6 | Observational crater-size groups and Stickney25 cap counts. |
| `catalog_gof_144cases_v1.csv` | 144 | 36 model cases x two catalogs x two masks, with unrotated G and nominal chi-square p (df=3), extracted from the rotation analysis table. |
| `rotation_144case_robustness_v30.csv` | 144 | Baseline test and classification sensitivity under 18 small diagnostic coordinate rotations. |
| `table4_speed_mixture_v1.csv` | 12 | Approximate rounded thresholds from manuscript v31 Table 4 for 7.95298 km/s plus a 20, 40 or 60 km/s fast component. `None` means no joint non-rejection was found on the scanned grid, **not zero**. |
| `obj_18runs_validation_v1.csv` | 18 | Individual v004 irregular-body, Phobos-gravity-ON source-radius validation trials. |
| `sphere_60runs_normalized_v2.csv` | 60 | Individual spherical-target, Phobos-gravity-OFF benchmark runs: two speeds x three source radii x ten seeds. |
| `sphere_6groups_summary_v2.csv` | 6 | Arithmetic means of individual-run metrics and their sample SD / standard error, **not** pooled weighted ratios. |

**Interpretation.** Four quadrants LM/LA/TM/TA refer to leading/trailing and Mars-facing/anti-Mars sectors. Stickney25 excludes the 25-degree angular cap centered at east longitude -50 degrees and latitude +1.3 degrees. The observational catalog counts are derived from Salamunićcar et al. (2014) and Karachevtseva et al. (2014); the catalogs share individual craters and are **not independent samples**. Tests report nominal, unadjusted p values for fixed model fields, not fitted population constraints.

**Important limitations.** The sphere summary's ratio columns are unweighted averages of per-run ratios; Supplementary Table S10 reports ratios based on pooled weighted production rates, so minor differences are expected. Numeric speeds near 7.95298 km/s may have small floating-point representation differences across source files. The original event catalogs, shape model and source crater catalogs are not distributed here. The root README and `docs/FIGURE_PROVENANCE_v3_JP.md` document the limits of figure-level reproduction.

**License.** Original author-owned contributions in the derived tables are distributed under [CC BY 4.0](LICENSE), with attribution to Hiroshi Kikuchi. The license does not grant any rights in the underlying third-party crater catalogs, Ernst shape models or other external inputs. See the [root README](../../README.md) for the software/data license split and citations to the originating crater catalogs.
