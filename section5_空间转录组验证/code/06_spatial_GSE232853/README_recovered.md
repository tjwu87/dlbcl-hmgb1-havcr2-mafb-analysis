# code/06_spatial_GSE232853/ — upstream core scripts (recovered)

`06_01_spatial_analysis.py` is a **plotting + statistics** script that reads 11 CSVs,
but the code producing those CSVs was originally absent from the repository. The five
scripts below were recovered from the Biomni platform's analysis records
(`01_对话记录/14_空间转录组分析.md`); parameters and thresholds are transcribed
**verbatim**.

| Script | Purpose | Output |
| --- | --- | --- |
| `06_00a_download_and_build_meta.py` | Download the GEO raw matrix + parse column names to build metadata + strict QC (drops 6 `Full ROI` samples) | `meta_df.csv`, `meta_filtered.csv`, `expression_raw_filtered.csv` |
| `06_00b_q3_norm_harmony_umap.py` | GeoMx Q3 normalisation + Harmony batch correction + PCA/UMAP | `expression_q3_lognorm.csv`, `Q3_normalization_stats.csv`, `task2_PCA_coordinates.csv`, `task2_UMAP_coordinates.csv`, `task2_highly_variable_genes.csv`, `GSE232853_adata_v2.h5ad` |
| `06_00c_dea_and_roi_tables.py` | Two DEA variants + paired ROI table + grouping table | `fig2_wilcoxon_*.csv`, `fig3_DEA_CD68_scanpy.csv`, `task4_DEA_CD68_DLBCL_vs_Normal.csv`, `fig4_paired_ROI_16gene_LATAM.csv`, `fig4_correlation_stats.csv`, `fig5_DLBCL_LATAM_group_assignment.csv` |
| `06_00d_gsea_ssgsea.py` | GSEA prerank + ssGSEA (both with MSigDB Hallmark 2020) | `fig5_GSEA_Hallmark_*.csv`, `fig5_GSEA_ranking_CD20_DLBCL.csv`, `fig6_ssGSEA_Hallmark_CD20_DLBCL.csv` |
| `06_00e_target_gene_score_master.py` | MAFB target-gene scoring on CD68+ ROIs + spatial co-variation master table | `fig7_master_spatial_coevolution.csv` |

## Data source (public on GEO, no registration)

```
GSE232853_Processed_data_CD20_CD68_final.csv.gz
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232853/suppl/
```

16,560 genes × 578 ROI samples. Adjacent ROIs on the same slide are masked with CD20 /
CD68 antibodies respectively, so each ROI yields a pair of CD20+ and CD68+ samples
(572 valid samples: CD20 253 + CD68 319).

## Key parameters (transcribed from the Biomni records)

| Step | Parameters |
| --- | --- |
| Q3 normalisation | `X / Q3_i * median(Q3)`, with Q3 = 75th percentile of each sample's non-zero genes; then `log1p` |
| HVG / PCA | `n_top_genes=2000` (seurat); `n_comps=50` |
| Harmony | `run_harmony(..., 'Slide_ID', max_iter_harmony=20, random_state=42)` |
| UMAP | `n_neighbors=15, n_pcs=30, min_dist=0.35, spread=1.2, random_state=42` |
| DEA | gene-wise `mannwhitneyu` + BH correction; the other variant is `sc.tl.rank_genes_groups(method='wilcoxon')` |
| LA_TAM_Score | mean of 16 genes (`FTL` is absent from the panel and was dropped) |
| Rich / Poor | median split of LA_TAM_Score within DLBCL samples (57 ROIs each) |
| GSEA prerank | `MSigDB_Hallmark_2020, min_size=10, max_size=500, permutation_num=1000, seed=42` |
| ssGSEA | `sample_norm_method='rank', scale=True`, taking the `NES` from the result matrix |

## Pitfalls (from the original records — following them avoids the same errors)

1. **Harmony output shape**: inside the container `ho.Z_corr` already comes back as
   `(samples, PCs)`; use it directly, **do not `.T` it**. An early transposition caused a
   shape mismatch and downstream errors. Also convert the PCA matrix to C-contiguous
   (`np.ascontiguousarray`) before handing it to `harmonypy`.
2. **`sc.get.rank_genes_groups_df` column names**: the returned columns must be renamed
   manually to `['gene','scores','logfoldchanges','pvals','pvals_adj','pts','pts_rest']`.
3. **ssGSEA input orientation**: gene × sample (i.e. transpose the expression matrix), and
   the `res2d` output needs
   `pivot(index='Name', columns='Term', values='NES')` to become samples × pathways.
4. **The original records use two output directories** (`GSE232853/` and
   `GSE232853_v2/`). This repository standardises on `GSE232853_v2/` (the one actually
   used in the paper and the one present in the data root).

## Environment

`scanpy`, `anndata`, `harmonypy`, `gseapy`, `seaborn`, `pandas`, `numpy`, `scipy`,
`statsmodels`. On its first run `gseapy` downloads the MSigDB Hallmark gene sets
(network access required).
