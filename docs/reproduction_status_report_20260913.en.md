# Reproduction status report (2026-09-13)

> **Archive note**: this report was written on 2026-09-13 and is the original record of
> that date. Its content about Figure 7 / `run_all.py --stage 7` (HAVCR2 structure and
> druggability) corresponds to the part later removed from the paper; **this archive
> (section1–6) does not contain that part's code** and the material is retained for
> traceability only.

> This report is based on **measurement**: it records the results of running each figure's
> reproduction chain inside the archive repository `DLBCL_HMGB1_HAVCR2_MAFB/` with
> `DLBCL_DATA_ROOT=<repo>/00_共用/数据`.
> Goal: to state clearly in one place what can already be reproduced with one command,
> what is missing, and how to fix it.

---

## 1. Summary

| Main figure | Reproduction entry point | Status | Fidelity |
|---|---|---|---|
| **Fig1** | `python tools/rebuild_fig1.py [--fast]` | ✅ runs | pixel difference vs the final version **0.039%** (only a tiny difference in 1f) |
| **Fig2** | `run_all.py --stage 2` + compose with `fig2.json` | 🟡 scripts complete, awaiting a real run for verification | — |
| **Fig3** | `python tools/rebuild_fig3.py` | ✅ runs | consistent with the latest final version (174 × 351 mm) |
| **Fig4** | `05_90_reproduce_figures_local.py` | ✅ runs | built-in numerical checks all pass; of 23 figures, **0 differ by >0.5%** |
| **Fig5** | `run_all.py --stage 5` + `fig5.json` | ✅ runs | data complete (all gaps are figures only) |
| **Fig6** | `Rscript 07_01_prognosis_main.R` / `07_02_lasso_cv_curves.R` | 🟡 environment/data ready, awaiting a real run | R packages 22/23 |
| **Fig7** | `run_all.py --stage 7` + `fig7.json` | ✅ runs | see §5 |

**Core conclusion: the plotting chain for all 7 main figures is usable; the archived data
is self-contained (no dependency on `D:/bulk-download`).**

---

## 2. Environment (measured, not a design draft)

`00_共用/环境/*.yml` describes **how to build the environment on a new computer**; the
environments actually available on this machine are:

| Interpreter | Version | Packages present | Purpose |
|---|---|---|---|
| **System Python310** (primary)<br>`C:\Users\admin\AppData\Local\Programs\Python\Python310\python.exe` | 3.10.11 | **30/36**: scanpy, anndata, infercnvpy, harmonypy, liana, celltypist, cellrank, gseapy, adjustText, pypdf, pypdfium2, reportlab, openpyxl, matplotlib, seaborn, statsmodels, scrublet, leidenalg, h5py, loompy … | **Reproduce all plotting and post-processing** |
| conda `scRNA`<br>`D:\ProgramFiles\Miniconda\envs\scRNA\python.exe` | 3.10.20 | 22 (no infercnvpy/liana/celltypist/pypdf) | Fallback |
| **R** `D:\R\R-4.4.2\bin\Rscript.exe` | 4.4.2 | **22/23** core packages | Fig6 prognostic modelling |
| WorkBuddy venv | 3.13.14 | pypdf / pypdfium2 / reportlab / matplotlib | Figure composition, rendering |

> ⚠ The three yml environments (dlbcl / celloracle / pyscenic) have **not** been created
> on this machine. They are only needed to **re-run the upstream analysis** (pySCENIC to
> build the GRN, CellOracle to build the oracle), because the intermediate results are
> already archived.

**Missing packages (needed only to "re-run the analysis")**:
`celloracle`, `pyscenic` (+ arboreto/ctxcore), `decoupler`, `pymol` (none of these in
Python310), and R's `estimate`.

---

## 3. Fixes made in this round (3, all side effects of tidying up)

### 3.1 `tools/04_01_gate_runner.py` had a stale path ✅ fixed

After tidying, the main script moved from `04_pySCENIC/` to
`section3_HAVCR2-MAFB调控模块/code/04_pySCENIC/`, so the hard-coded `MAIN` path in the
runner could no longer find the file. It now resolves **new location first → old location
as fallback → global best-effort search**.

### 3.2 The archived `GSE182434/scenic/` was missing 33 PNGs ✅ restored (15.7 MB)

`need_redraw()` in `04_01_pyscenic_downstream_stats_figures.py` decides whether to skip a
section by testing whether `OUT_DIR/{name}.png` exists. Tidying had cleared all the PNGs,
so every gate was judged "needs redrawing", which then executed a **historically broken
block** inside the script (line 214 uses `ax_heat`, which is not defined until line 1454)
→ `NameError`. The script's own comment states explicitly that it relies on "skip if the
figure already exists" to bypass that code. **Restoring the PNGs restores the original
gating logic.**

### 3.3 The archived `celloracle0331/` was missing 94 files ✅ restored (30.2 MB)

What was missing was exactly the `figures/` **SVG/PNG sources for all of Fig4's panels**
(the `OUT_DIR` of `05_90_reproduce_figures_local.py`). The originals are preserved in
`00_共用/历史版本/fig4原产物_20260913/`.

> **Shared lesson**: the rule used when tidying up to "exclude figures" was too
> aggressive. PNG/SVG inside a data directory is not necessarily "pure output" — some of
> it acts as a **cache marker** or as a **panel source**, and deleting it breaks
> reproducibility. ⇒ When cleaning a data directory in future, check each directory for
> scripts that decide based on the presence of an output.

---

## 4. Overview of archive gaps (`00_共用/数据/` vs the original sources)

**1017** files are missing, broken down as:

| Category | Count | Handling |
|---|---|---|
| Figures (png/svg/pdf) | 847 | mostly **non-blocking** (the scripts regenerate them); the 2 essential cases have been restored |
| Data-type files (csv/h5ad/npy/json/txt…) | 153 | checked item by item, see below |
| Archives (zip) | 3 | deliberately excluded ✅ |
| **`project0428/`** (75) | 75 | **unrelated to this paper** (GSE253902/GSE314596 gate analysis), deliberately excluded ✅ |
| `GSE232853/` (original version, 19) | 19 | the paper uses the `_v2` reprocessed version, deliberately excluded ✅ |
| `GSE182434/celloracle_paper_final_20260328_160843/` (54) | 54 | an older run; the paper uses the **0331** version, not restored for now |

**Conclusion: none of the data-type gaps block reproduction.**

---

## 5. Reproduction notes per figure

### Fig1 (`tools/rebuild_fig1.py`)
- Panels: a = truncated `10_00` script; b = hand-written dot plot; c = secondary output of the same run as a; d/e = static assets; f = volcano extracted from `02_02` and re-rendered
- Finalised parameters: 1a uses `pad_pt=10, th=252` (prevents top clipping); `--fast` reuses the existing 1a/1c
- Measured: 493 × 608 pt, **0.039% pixel difference** vs `终版_20260912/Fig1_细胞图谱与CNV.pdf`

### Fig3 (`tools/rebuild_fig3.py`)
- Gating: `PYSCENIC_FORCE_GATE=<output name>` re-renders only the specified section
- 3d uses the single-section extractor `fig3d_rss_extract.py` (the gate version in the main file has a different layout from the final version)
- Measured: 493 × 996 pt (174 × 351 mm); consistent in both layout and content with the 2330 version in the final directory (that final version is a snapshot from 09-12 23:33, and `fig3.json` was adjusted again at 09-13 00:56)

### Fig4 (`section4.../05_CellOracle/05_90_reproduce_figures_local.py`)
- Reads the CSVs needed by the 23 figures under `celloracle0331/`
- **The script has built-in numerical checks**, and all of them pass in practice:
  `SNR Mono=0.6561 / IFN_TAM=2.0028 / LA_TAM=7.1464`, matching the reference values exactly
- Comparison of outputs against the original outputs: **0 figures differ by >0.5%**; most are pixel-identical. Only 3 (`fig_effect_size_summary_dotplot`, `fig_horizontal_heatmap_all_deg`, `fig_merged_ora_bubble_all_cells`) differ in size — the originals were produced through `run_with_retrofit.py 3.2 7` scaling (960 px wide), whereas a bare run gives the native large figure. **Reproducing those three requires the retrofit route.**

### Fig5 (`06_01_spatial_analysis.py`)
- **Successfully run, exit=0.** Output goes to `00_共用/数据/GSE232853_v2/figures_reproduced/`
  (including `fig2_violin_..._Q3corrected`, `S5C_2x2_stratified_violin_Q3corrected`,
  `expression_q3_lognorm.csv` 16560 × 572, `Table_S8_*`, `fig2_three_gene_summary.csv`).
- Data is complete: the 62 gaps in `GSE232853_v2` are **all figures**, and the 26 gaps in
  `spatial_analysis` are **all figures**.

### Fig6 (R)
- Environment: R 4.4.2, **22 of 23** core packages installed, only `estimate` missing.
- Data is complete: the 63 gaps in `DLBCL_prognosis/` are **all figures under `figures/`**;
  the key `*_data.csv` files (LASSO / Cox / calibration / multi-cohort KM / immune
  infiltration) and `GSE10846_sur_model.Rdata` are **all present**.
- ⚠ **Parts that cannot be reproduced locally**: CIBERSORT depends on `LM22.txt` +
  `CIBERSORT.R`, which are not available locally (only the calling script
  `07_05_cibersort.R` exists in the Zenodo snapshot). They must be downloaded after
  registering at <https://ciberx.stanford.edu/>, **for academic use only**, and are
  impractical to redistribute.

### Fig7
- Scripts are complete (`08_01_ic50_and_boltz2.py`, 2786 lines; most data hard-coded plus
  2 CSVs read; no network dependency).
- A measured run **had not finished after >26 minutes** (no errors); the cause is
  unverified (possibly reading the 260 MB
  `GDSC2_RNAseq_log2TPM_training_cells.csv`, or large-scale regression/permutation tests;
  block buffering of the log may also be hiding progress).
- **c/d/g/h are PyMOL-rendered bitmaps and f needs GROMACS output**; these are external-tool
  products and cannot be re-exported as vectors from within Python.

---

## 6. How to reproduce

```bash
# 0) Data root: resolved automatically to 00_共用/数据/ by default, no setup needed
cd <repo>/DLBCL_HMGB1_HAVCR2_MAFB

# 1) Full analysis (by section)
python run_all.py --check          # input pre-check
python run_all.py --list           # list each section's scripts
python run_all.py --stage 4        # section4 only
python run_all.py                  # everything

# 2) One-command rebuild of the main figures
python tools/rebuild_fig1.py --fast
python tools/rebuild_fig3.py
python tools/fig_compose_pdf.py results/assembled/figN.json   # N=2,4,5,6,7

# 3) Re-render a single panel (shrink the canvas, keep the font size)
python tools/run_with_retrofit.py <script.py> <display_width_in> <font_pt> --out <output dir>
```

> Choosing the Python interpreter: `run_all.py` uses `python` from PATH — make sure it is
> the **system Python310**. Or:
> `export PATH="/c/Users/admin/AppData/Local/Programs/Python/Python310:$PATH"`

---

## 7. Remaining to-dos

1. **Real runs of Fig2 / Fig6 for verification** (scripts and data are both ready; no
   panel-by-panel comparison has been done yet)
2. **Install the missing packages** (only needed to re-run the pySCENIC / CellOracle upstream)
3. **Fig6 CIBERSORT resources**: `LM22.txt` + `CIBERSORT.R` (academic registration required)
4. **Fig7 external tools**: PyMOL (bitmaps for c/d/g/h), GROMACS (the RMSD trajectory in f)
5. **Fig2a/2b panel verification**: the output of `10_00_umap_cnv_myeloid_panels.py` needs
   to be aligned with the final version
6. **Clean up redundancy**: 1740 groups of duplicate files (6.5 GB) in the project, awaiting
   the user's decision on the authoritative paths
