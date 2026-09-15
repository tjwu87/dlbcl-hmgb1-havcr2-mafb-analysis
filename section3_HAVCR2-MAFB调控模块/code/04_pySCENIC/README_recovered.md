# code/04_pySCENIC/ — upstream core scripts (recovered)

`04_01_pyscenic_downstream_stats_figures.py` only does **statistics and plotting**; it
reads `<data root>/GSE182434/scenic/data_*.csv`. Those CSVs are produced by the pySCENIC
core, for which the repository originally had no code. The three scripts below were
recovered from the Biomni platform's analysis records; parameters and database URLs are
transcribed **verbatim**.

| Script | Purpose |
| --- | --- |
| `04_00a_setup_scenic_db.sh` | Downloads the three pySCENIC reference databases (TF list / hg38 motif rankings / motif–TF annotations) and documents the three official CLI commands in comments |
| `04_00b_build_loom.py` | Builds the loom input (Mono + IFN_TAM + LA_TAM = 327 cells × 9,118 genes, filtered to genes expressed in ≥ 5% of cells) |
| `04_00c_pyscenic_run_and_export.py` | GRNBoost2 → co-expression modules → RcisTarget (NES ≥ 3.0) → AUCell, plus export of the downstream `data_*.csv` files |

## Environment

`pyscenic`, `arboreto`, `ctxcore`, `loompy`, plus `scanpy`/`anndata` to read h5ad.

## Why the Python API rather than the CLI

The Biomni records show that the CLI's dask backend was unstable in the container
(`pyscenic ctx` produced zero bytes of output for a long time, and `aggregate_func` is
incompatible with the newer dask `from_delayed`). The final run used the Python API with
`client_or_address='custom_multiprocessing'` to bypass dask; the threshold parameters are
identical to the CLI. The equivalent CLI commands are in the header comment of `04_00a`.

## Prerequisites

`04_00b` needs an h5ad with annotated myeloid subtypes (a `mac_subtype` column; cells
outside Mono/IFN_TAM/LA_TAM are excluded). That annotation is a single-cell upstream
step and is outside the scope of this repository.
