# code/legacy/ — archived items (not part of the `run_all.py` pipeline)

These two files are the authors' **original master scripts**, taken verbatim from
`Tcell/` (only a "archive inclusion note" header was prepended and server absolute
paths were wrapped in `translate_path()`; the body logic is unchanged).

## Why they are included

`section6`'s `07_01_prognosis_main.R` **reads** the tables below, but `07_01`–`07_05`
only read and never write — the repository previously had no script that could produce
them. These two files are the producers:

| Intermediate table | `yuhou.R` (25-gene version) | `yuhou2.R` (23-gene version) |
|---|---|---|
| `LASSO_prognostic_formula.csv` | L289 | L80 |
| `GSE10846_risk_scores.csv` | L304 | L96 |
| `Cox_univariable_multivariable_results.csv` | **L957** | L411 |
| `Calibration_data_1_3_5yr.csv` | **L1571** | L528 |
| `Immune_infiltration_scores.csv` | **L1731** | — |
| `Immune_RiskScore_correlations.csv` | **L1734** | — |
| `All_cohorts_risk_scores.csv` (5 cohorts) | **L3356** | — |
| `MultiCohort_KM_statistics.csv` | **L3404** (5 cohorts) | from L1995 (6 cohorts) |
| `TCGA_DLBC_risk_scores.csv` | **L3414** | — |
| `GSE87371/11318/181063_risk_scores.csv` | **L3421/3426/3431** | — |

> **How this was established**: in the data root, `DLBCL_prognosis/tables/`
> `Cox_univariable_multivariable_results.csv` has the columns
> `Variable,label,type,HR,CI_low,CI_high,pval,pval_str,hr_str`, which matches
> `yuhou.R` L957's `select(Variable, label, type=section, HR, CI_low, CI_high,
> pval, pval_str, hr_str)` **character for character**; the 5 cohorts in
> `All_cohorts_risk_scores.csv` (412/221/199/882/45) are listed explicitly in
> `yuhou.R` L3350-3351. `LASSO_prognostic_formula.csv` contains 15 selected genes,
> all of which lie inside the 23-gene candidate set of `07_02_lasso_cv_curves.R`.

## Relationship between the two

| | `yuhou.R` | `yuhou2.R` |
|---|---|---|
| Gene set | 25 (HMGB1 + HAVCR2 + 23) | **23** (without HMGB1/HAVCR2) = the candidate set of `07_02` |
| Output directory | `/mnt/results/DLBCL_prognosis/` → `data/DLBCL_prognosis/` | `/mnt/results/DLBCL_prognosis_v2/` → same (merged by LEGACY_MAP) |
| Cohorts | **5** (GSE10846/87371/11318/181063/TCGA) = those used in the paper | 6 (one extra: **NCICCR-DLBCL**) |
| Relation to the paper | **Produces the tables actually present in the data root** | Parallel branch; its LASSO part is already covered by `07_02` |

## Blocks that can be skipped / are not required (line numbers = original line numbers)

`yuhou.R`
- L225-268 — exploration: switching to all 420 samples, trying elastic net (α=0.5), trying plain LASSO
- L374-442 / L455-541 / L542-720 — **three successive redraw attempts** of the CV curve and coefficient trajectory; the later ones supersede the earlier
  (the corresponding figures are now produced by `07_01`'s `A0*`)
- L2064-2205 — one-off variable/object reload debugging block
- L2897-2955 — GSE87371 directionality investigation (**but the `cens_os` direction fix at L2992 is required — do not skip it**)

`yuhou2.R`
- L685-1150 — the **NCICCR-DLBCL** cohort (not used in the paper; includes the GDC download and a comparison of two scaling schemes)
- L1538-1690 — a TCGA implementation via TCGAbiolinks (`yuhou.R` uses the GDC API; either one is enough)
- L2020-2394 — repeated recolouring/relayout redraws of the 6-cohort KM plot

## Running notes

- ⚠️ **These two files are hybrid R + embedded Python working files**: the body contains
  Python fragments such as `import requests` / `pd.read_csv(...)` (`yuhou.R` from about
  L2371, `yuhou2.R` from about L681 handle TCGA/GDC interaction). The file therefore
  **cannot be `source()`d as a whole** — `Rscript -e 'parse("<file>")'` errors at the
  Python blocks (**this is true of the original, not something introduced by archiving**).
  It was evidently executed manually, block by block. The purpose of including them is to
  **preserve the producing logic and parameters**, not to provide a one-click entry point.
- They **depend on upstream products**: `GSE10846_series_matrix.txt.gz` (publicly
  downloadable from GEO, see `data/README.md`), `LASSO_prognostic_formula.csv`, and so on.
  They are not a "download the repository and run" entry point.
- The original used `/workspace` as a staging directory; this maps to
  `<data root>/external/workspace`, so it is worth creating that directory before the
  first run.
- Required R packages: `GEOquery`, `Biobase`, `hgu133plus2.db`, `illuminaHumanv4.db`,
  `glmnet`, `survival`, `survminer`, `rms`, `ggplot2`, `patchwork`, `cowplot`,
  `dplyr`, `stringr`, `httr`/`jsonlite` (TCGA GDC API); `yuhou2.R` additionally needs
  `reticulate`, or manual block-by-block execution of its Python sections.

## Changes made when archiving (auditable)

Relative to the originals under `Tcell/`, the archived copies differ in **only** two
ways, verified line by line with a level diff:

1. an "archive inclusion note" comment block was prepended to the file header
   (29 lines for `yuhou.R`, 25 lines for `yuhou2.R`);
2. server absolute path literals were wrapped in `translate_path()`
   (80 occurrences in `yuhou.R`, 33 in `yuhou2.R`).

Other than that, **not a single line of body logic was changed**.
