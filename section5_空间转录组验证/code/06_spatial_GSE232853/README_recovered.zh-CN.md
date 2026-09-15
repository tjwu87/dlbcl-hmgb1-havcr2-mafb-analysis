# code/06_spatial_GSE232853/ —— 上游本体脚本（找回件）

`06_01_spatial_analysis.py` 是**出图 + 统计**脚本，它读 11 张 CSV，但这些 CSV 的
产出代码原先不在仓库里。下面 5 个脚本从 Biomni 平台的分析记录
（`01_对话记录/14_空间转录组分析.md`）中找回，参数与阈值均为**原文照录**。

| 脚本 | 作用 | 产出 |
| --- | --- | --- |
| `06_00a_download_and_build_meta.py` | 下载 GEO 原始矩阵 + 解析列名建 meta + 严格 QC（剔除 6 个 Full ROI） | `meta_df.csv`、`meta_filtered.csv`、`expression_raw_filtered.csv` |
| `06_00b_q3_norm_harmony_umap.py` | GeoMx Q3 归一化 + Harmony 批次校正 + PCA/UMAP | `expression_q3_lognorm.csv`、`Q3_normalization_stats.csv`、`task2_PCA_coordinates.csv`、`task2_UMAP_coordinates.csv`、`task2_highly_variable_genes.csv`、`GSE232853_adata_v2.h5ad` |
| `06_00c_dea_and_roi_tables.py` | 差异表达（两套口径）+ 配对 ROI 表 + 分组表 | `fig2_wilcoxon_*.csv`、`fig3_DEA_CD68_scanpy.csv`、`task4_DEA_CD68_DLBCL_vs_Normal.csv`、`fig4_paired_ROI_16gene_LATAM.csv`、`fig4_correlation_stats.csv`、`fig5_DLBCL_LATAM_group_assignment.csv` |
| `06_00d_gsea_ssgsea.py` | GSEA prerank + ssGSEA（均用 MSigDB Hallmark 2020） | `fig5_GSEA_Hallmark_*.csv`、`fig5_GSEA_ranking_CD20_DLBCL.csv`、`fig6_ssGSEA_Hallmark_CD20_DLBCL.csv` |
| `06_00e_target_gene_score_master.py` | CD68+ ROI 的 MAFB 靶基因打分 + 空间共变主表 | `fig7_master_spatial_coevolution.csv` |

## 数据来源（GEO 公开，无需注册）

```
GSE232853_Processed_data_CD20_CD68_final.csv.gz
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232853/suppl/
```

16 560 基因 × 578 ROI 样本。同一张切片相邻 ROI 分别用 CD20 / CD68 抗体掩膜，
因此每个 ROI 有 CD20+ 与 CD68+ 两个配对样本（572 个有效样本：CD20 253 + CD68 319）。

## 关键参数（照 Biomni 原文）

| 步骤 | 参数 |
| --- | --- |
| Q3 归一化 | `X / Q3_i * median(Q3)`，Q3 = 每样本非零基因 75 分位；再 `log1p` |
| HVG / PCA | `n_top_genes=2000`（seurat）；`n_comps=50` |
| Harmony | `run_harmony(..., 'Slide_ID', max_iter_harmony=20, random_state=42)` |
| UMAP | `n_neighbors=15, n_pcs=30, min_dist=0.35, spread=1.2, random_state=42` |
| DEA | `mannwhitneyu` 逐基因 + BH 校正；另一套为 `sc.tl.rank_genes_groups(method='wilcoxon')` |
| LA_TAM_Score | 16 基因（`FTL` 不在 panel 故剔除）的均值 |
| Rich / Poor | 在 DLBCL 样本内按 LA_TAM_Score 中位数二分（各 57 个 ROI） |
| GSEA prerank | `MSigDB_Hallmark_2020, min_size=10, max_size=500, permutation_num=1000, seed=42` |
| ssGSEA | `sample_norm_method='rank', scale=True`，取结果矩阵的 `NES` |

## 踩坑记录（照原文，读者照做可避）

1. **Harmony 输出的形状**：`ho.Z_corr` 在容器里返回的已经是 `(样本, PC)`，
   直接用即可，**不要 `.T`**。早期误转置导致形状错乱、后续报错。
   喂给 `harmonypy` 前要把 PCA 矩阵转成 C 连续（`np.ascontiguousarray`）。
2. **`sc.get.rank_genes_groups_df` 的列名**：返回列需手动重命名为
   `['gene','scores','logfoldchanges','pvals','pvals_adj','pts','pts_rest']`。
3. **ssGSEA 输入方向**：gene × sample（即表达矩阵要转置），输出 `res2d`
   需 `pivot(index='Name', columns='Term', values='NES')` 才是 样本 × 通路。
4. **原文中的输出目录有两个**（`GSE232853/` 与 `GSE232853_v2/`）。本仓库统一到
   `GSE232853_v2/`（即论文实际使用、也是数据根里现存的那一份）。

## 环境

`scanpy`、`anndata`、`harmonypy`、`gseapy`、`seaborn`、`pandas`、`numpy`、`scipy`、`statsmodels`。
`gseapy` 首次运行会自动下载 MSigDB Hallmark 基因集（需要网络）。
