# DLBCL HMGB1–HAVCR2–MAFB multi-omics analysis pipeline

**Multi-omics analysis code** for the association of HMGB1, HAVCR2 and MAFB with T-cell
function in the diffuse large B-cell lymphoma (DLBCL) tumour microenvironment.
Organised by the Results subsections of the paper:

| Part | Content | Figures |
|---|---|---|
| `section1_恶性B细胞与髓系通讯` | scRNA-seq QC/integration/annotation, malignancy calling, LIANA cell–cell communication | Figure 1 (a–h), S1 |
| `section2_髓系轨迹与LA_TAM终末态` | Myeloid re-subtyping, PAGA/CellRank trajectory, LA-TAM scoring | Figure 2 (a–f), S2, S8 |
| `section3_HAVCR2-MAFB调控模块` | pySCENIC regulon activity, RSS receptor prioritisation, pseudotime correlation | Figure 3 (a–g), S3, S9 |
| `section4_MAFB计算机扰动` | CellOracle: GRN → perturbation simulation → fate probability → sensitivity | Figure 4 (a–h) |
| `section5_空间转录组验证` | GeoMx: Q3 normalisation → paired ROIs → DEA → GSEA → spatial co-variation | Figure 5 (a–f), S5 |
| `section6_MAFB预后评分模型` | Bulk multi-cohort LASSO-Cox, KM/calibration, CIBERSORT immune infiltration | Figure 6 (a–g), S6 + supplementary tables |

> ⚠ **This version contains only the code for Results 1–6.** The paper removed
> "HAVCR2 structure and druggability" (formerly Figure 7 / S7: GDSC IC50, Boltz-2,
> molecular docking, MD), so that part is not included in this repository.
>
> ⚠ This study is **entirely computational**; no experimental validation is included.
> Terms such as "axis / driver / targeting" are **hypothesis-generating** statements —
> please read them as the original text's "candidate / associated with".

> ℹ **Language note.** Documentation files are provided in English; the original
> Chinese versions are kept alongside as `*.zh-CN.md` (e.g. `README.zh-CN.md`).
> Source code and its comments remain in Chinese and have **not** been modified, so
> that the translation cannot affect any numerical result.

---

## Directory layout

```
.
├── run_all.py              # one-shot entry point (--stage 1..6 / --list / --check / --dry-run)
├── config/                 # path resolution (paths.py/paths.R), plotting conventions, render wrapper, parameters
├── section1..section6/     # one directory per Results subsection: README.md + code/<module>/
│   └── code/config/        # per-section copy of config/ (so a section can be run on its own;
│                           #   identical in content to the root config/)
├── tools/                  # figure reproduction pipeline (Fig1/Fig3 one-click rebuild, generic composer, panel extractors)
├── docs/                   # panel ↔ code mapping, version audit, reproduction guide and status report
├── 00_download/            # public data download script
├── environment.yml         # main conda environment (+ environment_celloracle/pyscenic.yml)
├── install_R_packages.R    # R dependencies
└── data/                   # data directory placeholder (not distributed with the repository, see data/README.md)
```

---

## Quick start

```bash
# 1) Environment
conda env create -f environment.yml && conda activate dlbcl
Rscript install_R_packages.R
#   CellOracle / pySCENIC have demanding dependencies; separate environments are recommended:
#   conda env create -f environment_celloracle.yml
#   conda env create -f environment_pyscenic.yml

# 2) Data (~8 GB, not included in this repository) -- see data/README.md
export DLBCL_DATA_ROOT=/path/to/data      # Windows: set DLBCL_DATA_ROOT=D:\bulk-download
python 00_download/download_all.py --list
python 00_download/download_all.py --all

# 3) One-shot reproduction
python run_all.py                 # everything (section1 -> section6)
python run_all.py --stage 4       # only the CellOracle perturbation
python run_all.py --list          # list each section's scripts and readiness
python run_all.py --check         # input-data pre-check only
python run_all.py --dry-run       # print the commands that would be executed
```

`DLBCL_DATA_ROOT` is the only path you need to configure. Every script resolves paths
through `config/paths.py` (on the R side, `config/paths.R`). The historical absolute
paths kept in the code (`D:\bulk-download\…`, `/mnt/results/…`) are mapped at runtime
to the current data root by `config.paths.translate()` — **no manual edits required**.

### ⚠ Known environment pitfall: legacy h5ad files cannot be read

The h5ad files distributed with the project were written by an older anndata version;
their `uns/log1p` carries `encoding_type='null'`, and anndata 0.11+ removed the
corresponding read method, so it raises
`IORegistryError: No read method registered for IOSpec(encoding_type='null', ...)`.
Run `python tools/fix_h5ad.py` to repair this: it **copies** the h5ad into
`results/data_fixed/` and drops that key (**the original data is not modified**), and
`config/paths.py` will automatically prefer the repaired copy.

---

## From code to figures

Main figures are produced by panel-level scripts and then assembled. One-click rebuild
and composition:

```bash
python tools/rebuild_fig1.py --fast                 # Figure 1
python tools/rebuild_fig3.py                        # Figure 3
python tools/fig_compose_pdf.py results/assembled/figN.json    # N = 2,4,5,6
```

| Figure | Main plotting script | Data source |
|---|---|---|
| Fig 1 | `section1/code/10_figures/10_00_umap_cnv_myeloid_panels.py` + `tools/rebuild_fig1.py` | `section1/code/01_*`, `02_*` |
| Fig 2 | `section2/code/03_trajectory_PAGA/03_01_paga_trajectory.py` | `section2/code/01_*` |
| Fig 3 | `section3/code/04_pySCENIC/04_01_pyscenic_downstream_stats_figures.py` | pySCENIC output CSV |
| Fig 4 | `section4/code/05_CellOracle/05_90_reproduce_figures_local.py` | `section4/code/step*` |
| Fig 5 | `section5/code/06_spatial_GSE232853/06_01_spatial_analysis.py` | `data/GSE232853_v2/` |
| Fig 6 | `section6/code/07_bulk_prognosis/07_01_prognosis_main.R`, `07_02_lasso_cv_curves.R` | `data/DLBCL_prognosis/` |

The mapping between panels and final figures is in `docs/figure_panel_map.md` and
`docs/figure_code_map_final.md`; per-section reproduction status is in
`docs/复现状态报告_20260913.md` (Chinese; an English version is in progress).

---

## Figure conventions

Journal page metrics: `\textwidth = 160 mm`, `\textheight = 216 mm`.
The older scripts opened canvases as large as 13×13 in or even 20×6.5 in while using
the default 10 pt font, so after layout scaling the on-figure text was only 1–2 pt —
**that**, not pixel count, is the real reason "the figures look blurry".

The conventions are now centralised in `config/plot_style.py`
(on the R side, `config/plot_style.R`):

```python
from config.plot_style import apply_paper_style, panel_figsize, save_figure

apply_paper_style()                                    # font from 8 pt, minimum 6 pt
fig, axes = plt.subplots(2, 3, figsize=panel_figsize(3, 2))
save_figure(fig, FIG_DIR / "Fig1")                     # writes PDF (vector) + PNG (600 dpi)
```

**Vector PDF is preferred** — it stays sharp at any magnification.

To retrofit existing scripts in bulk (without editing them, by running a wrapper):

```bash
python tools/run_with_retrofit.py <script.py> <display_width_in> <font_pt>
Rscript tools/run_with_retrofit.R <script.R> <display_width_in> <font_pt>
```

How it works is described in `config/retrofit.py` / `config/retrofit.R`: it intercepts
the canvas size passed to `plt.subplots` / `ggsave` and scales the font size
accordingly, redirecting all output to `results/` — **without overwriting the original
images**.

---

## Where the analysis parameters live

Analysis code and hyper-parameters are **embedded in the pipeline scripts** (not held in
a central config), so opening a script shows them directly:
QC/doublet removal (scrublet threshold), normalisation and HVG (`n_top_genes`),
integration (harmonypy), clustering (leiden resolution), PAGA/DPT parameters, pySCENIC
downstream statistics, CellOracle per-step parameters (step00–step09: GRN filtering,
perturbation simulation, sensitivity grid), spatial analysis (Q3 correction), and
LASSO/Cox (λ grid, cross-validation folds).

Randomness: `seed = 42` (see `config/params.yaml`). The thresholds of interest to
reviewers (LIANA `magnitude_rank < 0.05`, GRN edge counts, HVG counts, …) are
collected in `config/params.yaml`, with accompanying sensitivity-analysis scripts
(`step08*`).

---

## Known gaps (stated honestly)

> This repository contains **code only**. "Download the repository → figures come out"
> does **not** hold: most inputs read by the scripts are **analysis intermediates**
> (~8 GB) that exist in no public download source. The gaps are listed by nature below
> (audited 2026-09-14; each item traces back to a file path).

### A. Data you must supply yourself (not in the repository)

| Category | Content | How to obtain |
|---|---|---|
| Public raw data | GSE182434 raw count matrix; GSE10846/GSE87371/GSE11318/GSE181063 series matrices; GeoMx GSE232853 raw RLT matrix; GDSC2 training matrix | `python 00_download/download_all.py --list` (the cohort list printed by that script is **not fully consistent with what the code actually uses**; treat the `GSE*` identifiers appearing in each section's code as authoritative) |
| Analysis intermediates | `adata_processed.h5ad`, `adata_mac_*`, pySCENIC outputs, `celloracle0331/`, `GSE232853_v2/`'s `expression_raw_filtered.csv`/`Q3_normalization_stats.csv`/`task2_*`, `DLBCL_prognosis/tables/*.csv`, `GSE10846_sur_model.Rdata`, etc. | **Must be provided by the authors** (cannot be downloaded from GEO) |
| Require registration | CIBERSORT's `LM22.txt`; CellOracle's `base_GRN_human_promoter.csv`; `Macro_mono_cellmarker.xlsx` | see `data/README.md` |

### B. Upstream scripts (originally 4 missing, **all now recovered**)

1. ~~**pySCENIC core**~~ **recovered**: `section3` originally held only the downstream
   statistics script `04_01_pyscenic_downstream_stats_figures.py`. The core was
   recovered from the Biomni analysis records and is now included as three files under
   `section3/code/04_pySCENIC/`:
   `04_00a_setup_scenic_db.sh` (reference database download + the three official CLI
   commands), `04_00b_build_loom.py` (builds the 327-cell × 9,118-gene loom), and
   `04_00c_pyscenic_run_and_export.py` (GRNBoost2 → RcisTarget NES ≥ 3.0 → AUCell, plus
   export of the `data_*.csv` consumed downstream). Parameters, database URLs and random
   seeds are **transcribed verbatim** from the original records.
   See `section3/code/04_pySCENIC/README_recovered.md`.
2. ~~**LIANA core**~~ **recovered**: `section1/code/02_cellcomm_LIANA/02_00_liana_run.py`
   is named "run" but is in fact a **plotting script** (it reads three result tables).
   Two core scripts are now included:
   `02_00a_liana_rank_aggregate_14celltypes.py` (all cells, 14 cell types) and
   `02_00b_liana_rank_aggregate_13subtypes.py` (13 subtypes, including the
   Malignant/Normal B split rule). Parameters
   `expr_prop=0.1 / min_cells=5 / n_perms=100 / seed=42` are transcribed verbatim.
   See `section1/code/02_cellcomm_LIANA/README_recovered.md`.
   In addition, the **malignancy calling** step (inferCNVpy →
   `cnv/malignancy_classification.csv`) and the **differential-communication volcano**
   (`cellcomm/volcano_data_malignant_vs_normal_monomac.csv`) were recovered as
   `code/01_scRNA_GSE182434/01_00_cnv_malignancy_classification.py` and
   `code/02_cellcomm_LIANA/02_00c_volcano_malignant_vs_normal.py`.
3. ~~**Spatial transcriptomics upstream**~~ **recovered**: the 11 CSVs read by
   `06_01_spatial_analysis.py` in `section5` now have producers — five scripts were added
   under `section5/code/06_spatial_GSE232853/`:
   `06_00a` (download the GEO raw matrix + sample metadata + strict QC),
   `06_00b` (GeoMx Q3 normalisation + Harmony batch correction + PCA/UMAP),
   `06_00c` (two DEA variants + paired ROI table + LA_TAM grouping table),
   `06_00d` (GSEA prerank + ssGSEA),
   `06_00e` (MAFB target-gene scoring + spatial co-variation master table).
   The raw data is publicly downloadable from GEO
   (`GSE232853_Processed_data_CD20_CD68_final.csv.gz`).
   See `section5/code/06_spatial_GSE232853/README_recovered.md`.
4. ~~**Several prognostic intermediate tables**~~ **recovered**: the
   `Cox_univariable_multivariable_results.csv`, `Calibration_data_1_3_5yr.csv`,
   `Immune_infiltration_scores.csv`, `All_cohorts_risk_scores.csv`,
   `GSE87371/GSE11318/GSE181063/TCGA_*_risk_scores.csv` and others read by
   `07_01_prognosis_main.R` are produced by the authors' original master scripts, now
   included as `section6/code/legacy/yuhou.R` (25-gene version, **which produces the
   tables actually present in the data root**) and
   `section6/code/legacy/yuhou2.R` (23-gene parallel branch, containing the NCICCR
   cohort that the paper does not use).
   The table ↔ line-number cross-reference and the list of "exploratory blocks that can
   be skipped" are in `section6/code/legacy/README.md`.
   `07_02_lasso_cv_curves.R` additionally produces `LASSO_prognostic_formula.csv` and
   `GSE10846_risk_scores.csv`.

### C. Known code-level issues (fixed or documented in the archive copy)

1. ✅ **`07_05_cibersort.R` originally lacked path migration** (no
   `source(config/paths.R)` in the source project, and three bare absolute paths) →
   the archive copy injects it and wraps the paths in `translate_path()`.
2. ✅ **`config/paths.py` / `paths.R` `LEGACY_MAP` was missing 2 prefixes**
   (`D:/BadiduNetdiskDownload/R/Tcell`, `…/R/GSE10846survive`) → now mapped to
   `data/external/Tcell` and `data/external/GSE10846survive`.
3. ⚠️ **Recomputing CIBERSORT requires `exp.txt`**: line 117 of `07_05` calls
   `cibersort(lm22f, "exp.txt", perm = 1000, QN = T)`, but `exp.txt` is **not generated
   by the script and is absent from both the repository and the data**; there is an
   `if (!file.exists("GSE10846_ciber.Rdata"))` short-circuit, i.e. it only runs when
   `GSE10846_ciber.Rdata` has been cached. To genuinely recompute, supply your own
   `exp.txt` (a gene × sample expression matrix) under `data/external/Tcell/`.
   ⚠ Conversely: **CIBERSORT's `LM22.txt` is read from the `extdata/` of the R package
   `CIBERSORT`** (`system.file("extdata","LM22.txt", package="CIBERSORT")`), not from the
   working directory — so that R package must be installed first
   (`install_R_packages.R` cannot install it automatically; it must be obtained
   manually).
4. ⚠️ **`00_download/download_all.py` cohort list once disagreed with reality**: as
   shipped, the script listed GSE32918 / GSE4475 / GSE23501, whereas `section6`'s code
   actually uses GSE10846 / GSE87371 / GSE11318 / GSE181063 / TCGA-DLBC.
   **The archive copy has been corrected** (the old identifiers survive only as an
   explanatory comment), and the wrong claim that GEO provides `adata_processed.h5ad`
   has been fixed.
5. ⚠️ **`/tmp/…` intermediate paths**: `section4`'s `step*` scripts use
   `/tmp/oracle_after_02B.pkl` and similar for inter-step staging, which resolves to a
   literal `\tmp` path on non-Unix systems (Windows). This is fine on Linux/macOS.
6. ⚠️ **In the source project, `config/paths.py` relied on a Windows junction**
   (`section*/code/config` was a directory junction pointing at the root `config/`).
   The archive has materialised these into real copies and switched `paths.py` to
   "self-locate the repository root upwards", which is friendlier for cross-platform use
   and for git.

### D. Other

- The audit and mapping documents under `docs/` record version differences encountered
  while tidying up (for example, some figures' panel sources come from an earlier
  version in the authors' asset library). They are for traceability only and do not
  affect the pipeline.
- Fig 4's panel sources are vector PDFs written directly by matplotlib; if you render
  the same SVG on another platform, note that semi-transparent fills (`fill-opacity`)
  may be lost.

---

## Data availability

All raw data come from public databases (GEO / GDC / GDSC / MSigDB); the list is
available via `python 00_download/download_all.py --list`.
This repository **contains** neither raw data nor large intermediate products.

---

## Citation and licence

- Code licence: **MIT** (see `LICENSE`)
- Citation metadata: see `CITATION.cff`
- If you use this code, please cite the original paper (details to be added).
