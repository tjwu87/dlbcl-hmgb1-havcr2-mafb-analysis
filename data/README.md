# Data directory (placeholder — not distributed with the repository)

This repository **contains code only**. All raw data come from public databases
(GEO / GDC / GDSC / MSigDB) and total roughly 8 GB. Download it yourself and place it
here, or point an environment variable at a location you already have:

```bash
# Option 1: put it here (recommended)
#   Unpack/download each dataset under data/, using directory names that match what
#   `run_all.py --check` reports:
#   data/GSE182434/       single-cell discovery cohort
#   data/GSE232853_v2/    GeoMx spatial transcriptomics (reprocessed version)
#   data/celloracle0331/  CellOracle perturbation outputs
#   data/DLBCL_prognosis/ bulk prognostic modelling
#   data/GDSC2/           drug sensitivity
#   data/spatial_analysis/

# Option 2: point at any existing directory
export DLBCL_DATA_ROOT=/path/to/data      # Linux / macOS
set    DLBCL_DATA_ROOT=D:\bulk-download   # Windows

# See which datasets are needed
python 00_download/download_all.py --list
```

`config/paths.py` resolves in this order: `DLBCL_DATA_ROOT` → `<repo>/data/` →
`data/` inside the repository.

Reference files that require registration and must be obtained manually (impractical to
distribute with the package):

| File | Place under | Notes |
|---|---|---|
| `GSE10846_series_matrix.txt.gz` | `data/external/GSE10846survive/` | Publicly downloadable from GEO; the R scripts in `section6` look for it at this path |
| `GSE10846_ciber.Rdata` | `data/external/Tcell/` | **Cache** of the CIBERSORT result (with it the script runs; delete it and a recompute is required) |
| `immu_check_point.txt` | `data/external/Tcell/` | Immune-checkpoint gene list (the script reads it from the working directory) |
| `DLBCL_MRGs_Cluster_PD1.csv` | `data/external/Tcell/` | TIDE input table |
| `exp.txt` | `data/external/Tcell/` | **Required only to recompute CIBERSORT** (gene × sample expression matrix; the script does not generate it) |
| the `CIBERSORT` R package | R library | Provides `cibersort()` and the bundled `extdata/LM22.txt`; register at <https://ciberx.stanford.edu/> (academic use only) |
| `base_GRN_human_promoter.csv` | `data/GSE182434/celloracle_rerun_20260324_035651/` | CellOracle prior GRN, see <https://github.com/morris-lab/CellOracle> |
| `Macro_mono_cellmarker.xlsx` | `data/` | Myeloid marker table (used by Fig1a/1c and the myeloid panels) |

> The two `data/external/` subdirectory names correspond to the `external/Tcell` and
> `external/GSE10846survive` entries in `LEGACY_MAP` in
> `config/paths.py` / `paths.R` — the historical absolute paths in the scripts
> (`D:/BadiduNetdiskDownload/R/…`) are translated here at runtime, so no code changes
> are needed.
