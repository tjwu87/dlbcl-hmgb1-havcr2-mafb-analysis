# code/02_cellcomm_LIANA/ —— 上游本体脚本（找回件）

`02_00_liana_run.py` / `02_01` / `02_02` 都是**读 CSV 出图**的脚本，
真正计算细胞通讯的 LIANA `rank_aggregate` 原先不在仓库里。下面这些脚本
从 Biomni 平台的分析记录中找回，参数（`expr_prop=0.1`、`min_cells=5`、
`n_perms=100`、`seed=42`）与 13 亚型构造规则均为原文照录。

| 脚本 | 产出 |
| --- | --- |
| `01_00_cnv_malignancy_classification.py`（在 `code/01_scRNA_GSE182434/`） | `cnv/malignancy_classification.csv` —— inferCNVpy 恶性/正常 B 细胞判定 |
| `02_00a_liana_rank_aggregate_14celltypes.py` | `liana_results_full.csv`、`liana_results_significant.csv`、`interaction_weight_matrix.csv` |
| `02_00b_liana_rank_aggregate_13subtypes.py` | `liana_results_13subtypes_full.csv`、`liana_results_13subtypes_significant.csv` |
| `02_00c_volcano_malignant_vs_normal.py` | `volcano_data_malignant_vs_normal_monomac.csv` —— 恶性 vs 正常 B → 巨噬细胞的差异通讯 |

13 亚型 = T/NK/pDC 由 `cell_type` 映射 ＋ B 细胞按 `cnv/malignancy_classification.csv`
拆 Malignant/Normal B ＋ 髓系按 5 组 marker 打分取最大分成 Mono/DC_1/DC_2/IFN_TAM/LA_TAM。

## 环境

`scanpy`、`anndata`、`liana`（≥1.1），CNV 那步另需 `infercnvpy`。输入为
`adata_processed.h5ad`，其 `.X` 须为 log1p 归一化表达（LIANA 按 natural
log-normalized 处理）。
