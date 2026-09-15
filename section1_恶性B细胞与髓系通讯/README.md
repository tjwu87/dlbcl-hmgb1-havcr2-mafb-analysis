# 恶性B细胞与髓系通讯

> **论文对应**：Results 1 — 恶性 B 细胞与髓系细胞通讯  
> **图件**：Figure 1（a–h）、Figure S1

## 分析内容

GSE182434 scRNA-seq：QC → Harmony 整合 → 聚类 → CellTypist/CellVote 注释 → inferCNVpy 恶性判定 → LIANA 细胞通讯 → 火山图

## 目录

```
section1_恶性B细胞与髓系通讯/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/01_scRNA_GSE182434/01_01_qc_integration.py` | |
| `code/01_scRNA_GSE182434/01_03_annotation_concordance.py` | |
| `code/01_scRNA_GSE182434/01_04_supplementary_panels.py` | |
| `code/02_cellcomm_LIANA/02_00_liana_run.py` | |
| `code/02_cellcomm_LIANA/02_01_fig1de_aggregate_circle_heatmap.py` | |
| `code/02_cellcomm_LIANA/02_02_liana_figures_13subtypes_volcano.py` | |
| `code/02_cellcomm_LIANA/02_03_liana_supplementary.py` | |
| `code/02_cellcomm_LIANA/02_00a_liana_rank_aggregate_14celltypes.py` | 上游本体（找回件）：LIANA 全细胞 14 类型 `rank_aggregate`，产出 `liana_results_full/significant.csv` 与 `interaction_weight_matrix.csv` |
| `code/02_cellcomm_LIANA/02_00b_liana_rank_aggregate_13subtypes.py` | 上游本体（找回件）：13 亚型 `rank_aggregate`，产出 `liana_results_13subtypes_full/significant.csv`（论文主图所用）—— 见 `code/02_cellcomm_LIANA/README_recovered.md` |
| `code/10_figures/10_00_umap_cnv_myeloid_panels.py` | |

## 输入数据（公用，位于 `data/`）

- `GSE182434/adata_processed.h5ad`
- `GSE182434/cnv/adata_cnv.h5ad`
- `GSE182434/cnv/malignancy_classification.csv`
- `GSE182434/cellcomm/liana_results_full.csv`
- `GSE182434/cellcomm/liana_results_13subtypes_full.csv`
- `GSE182434/cellcomm/volcano_data_malignant_vs_normal_monomac.csv`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 1
```

## 环境状态

- 需要：scanpy, anndata, harmonypy, celltypist, GEOparse, infercnvpy, liana, decoupler
- 本机 `scRNA` 环境已具备：scanpy, anndata, harmonypy
- 尚缺：celltypist, GEOparse, infercnvpy, liana, decoupler

## 产出图件

- `figures/论文成图/Fig1.png`
- `figures/论文成图/FigS1.png`
- `figures/论文成图/Figure S1.png`
- `figures/论文成图/Figure1.pdf`
- `figures/论文成图/Figure1.png`
- `figures/重排版/Fig1.png`
- `figures/面板原件/1A.png`
- `figures/面板原件/1B.png`
- `figures/面板原件/1C.png`
- `figures/面板原件/1D.png`
- `figures/面板原件/1E.png`
- `figures/面板原件/1F.png`
- `figures/面板原件/1G.png`
- `figures/面板原件/1H.png`
- `figures/面板原件/Figure S1-1.png`
- `figures/面板原件/Figure S1-2.png`
- `figures/面板原件/Figure1.1.png`
- `figures/面板原件/S1.png`
- `figures/面板原件/celltype_composition_per_sample2.png`
- `figures/面板原件/chromosome_heatmap.png`
- `figures/面板原件/circle_CD8T_receiver.png`
- `figures/面板原件/circle_CD8T_sender.png`
- `figures/面板原件/circle_IFNTAM_receiver.png`
- `figures/面板原件/circle_IFNTAM_sender.png`
- `figures/面板原件/circle_LATAM_receiver.png`
- `figures/面板原件/circle_LATAM_sender.png`
- `figures/面板原件/circle_MalignantB_receiver.png`
- `figures/面板原件/circle_MalignantB_sender.png`
- `figures/面板原件/circle_plot_13subtypes.png`
- `figures/面板原件/circle_plot_all.png`
- `figures/面板原件/circle_plot_malignant_vs_normal.png`
- `figures/面板原件/cnv_score_distribution.png`
- `figures/面板原件/dotplot_malignant_normal_toplr.png`
- `figures/面板原件/heatmap_13subtypes.png`
- `figures/面板原件/heatmap_7ct_focused_w.png`
- `figures/面板原件/heatmap_7subtypes.png`
- `figures/面板原件/heatmap_malignant_vs_normal.png`
- `figures/面板原件/malignancy_bcell_umap_composition.png`
- `figures/面板原件/metabolic_heatmap2.png`
- `figures/面板原件/metabolic_violin2.png`
- `figures/面板原件/network_MAFB_tree_organic.png`
- `figures/面板原件/supp_A_pseudotime_distribution.png`
- `figures/面板原件/supp_B_paga_connectivity_heatmap.png`
- `figures/面板原件/supp_C_sample_composition.png`
- `figures/面板原件/supp_D_cellcomm_strength.png`
- `figures/面板原件/umap_celltypes2.png`
- `figures/面板原件/umap_overview2.png`
- `figures/面板原件/volcano_malignant_vs_normal_monomac.png`
- `figures/面板原件/volcano_receptor_pseudotime_corr.png`
