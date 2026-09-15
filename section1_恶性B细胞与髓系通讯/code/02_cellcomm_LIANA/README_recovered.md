# code/02_cellcomm_LIANA/ — upstream core scripts (recovered)

`02_00_liana_run.py` / `02_01` / `02_02` are all scripts that **read CSV and draw
figures**; the LIANA `rank_aggregate` computation itself was originally absent from the
repository. The scripts below were recovered from the Biomni platform's analysis
records; the parameters (`expr_prop=0.1`, `min_cells=5`, `n_perms=100`, `seed=42`) and
the 13-subtype construction rules are transcribed **verbatim**.

| Script | Output |
| --- | --- |
| `01_00_cnv_malignancy_classification.py` (in `code/01_scRNA_GSE182434/`) | `cnv/malignancy_classification.csv` — malignancy vs normal B-cell calls by inferCNVpy |
| `02_00a_liana_rank_aggregate_14celltypes.py` | `liana_results_full.csv`, `liana_results_significant.csv`, `interaction_weight_matrix.csv` |
| `02_00b_liana_rank_aggregate_13subtypes.py` | `liana_results_13subtypes_full.csv`, `liana_results_13subtypes_significant.csv` |
| `02_00c_volcano_malignant_vs_normal.py` | `volcano_data_malignant_vs_normal_monomac.csv` — differential communication malignant vs normal B → macrophages |

The 13 subtypes = T/NK/pDC mapped from `cell_type` + B cells split into
Malignant/Normal B using `cnv/malignancy_classification.csv` + myeloid cells split into
Mono/DC_1/DC_2/IFN_TAM/LA_TAM by taking the argmax of five marker-set scores.

## Environment

`scanpy`, `anndata`, `liana` (≥ 1.1), plus `infercnvpy` for the CNV step. The input is
`adata_processed.h5ad`; its `.X` must be log1p-normalised expression (LIANA treats it as
natural log-normalised).
