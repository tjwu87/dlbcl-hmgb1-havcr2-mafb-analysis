# 空间转录组验证

> **论文对应**：Results 5 — 空间转录组验证（GeoMx）  
> **图件**：Figure 5（a–f）、Figure S5

## 分析内容

GSE232853 GeoMx：Q3 归一化 → ROI/AOI 表达 → 配对 ROI 比较 → DEA → GSEA/ssGSEA → 空间共变

## 目录

```
section5_空间转录组验证/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/06_spatial_GSE232853/06_01_spatial_analysis.py` | |
| `code/06_spatial_GSE232853/06_00a_download_and_build_meta.py` | 上游本体（找回件）：下载 GSE232853 GeoMx 原始矩阵、构建样本元数据、严格 QC |
| `code/06_spatial_GSE232853/06_00b_q3_norm_harmony_umap.py` | 上游本体（找回件）：GeoMx Q3 归一化 + Harmony 批次校正 + PCA/UMAP |
| `code/06_spatial_GSE232853/06_00c_dea_and_roi_tables.py` | 上游本体（找回件）：两套 DEA 口径 + 配对 ROI 表 + LA_TAM 分组表 |
| `code/06_spatial_GSE232853/06_00d_gsea_ssgsea.py` | 上游本体（找回件）：GSEA prerank + ssGSEA（MSigDB Hallmark 2020） |
| `code/06_spatial_GSE232853/06_00e_target_gene_score_master.py` | 上游本体（找回件）：CD68+ ROI 的 MAFB 靶基因打分 + 空间共变主表 —— 见 `code/06_spatial_GSE232853/README_recovered.md` |

## 输入数据（公用，位于 `data/`）

- `GSE232853_v2/expression_q3_lognorm.csv`
- `GSE232853_v2/meta_filtered.csv`
- `GSE232853_v2/fig4_paired_ROI_16gene_LATAM.csv`
- `GSE232853_v2/fig5_GSEA_Hallmark_LATAM_Rich_vs_Poor_CD20.csv`
- `GSE232853_v2/fig7_master_spatial_coevolution.csv`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 5
```

## 环境状态

- 需要：scanpy, anndata, gseapy, decoupler
- 本机 `scRNA` 环境已具备：scanpy, anndata, gseapy
- 尚缺：decoupler

## 产出图件

- `figures/论文成图/Fig5.png`
- `figures/论文成图/FigS5.png`
- `figures/论文成图/Figure S5.png`
- `figures/论文成图/Figure5.pdf`
- `figures/论文成图/Figure5.png`
- `figures/重排版/Fig5.png`
- `figures/面板原件/Figure S5 - 副本.png`
- `figures/面板原件/Figure S5.1.png`
- `figures/面板原件/Figure5.1.png`
- `figures/面板原件/S5.png`
- `figures/面板原件/S5A_ROI_AOI_design_schematic.png`
- `figures/面板原件/S5B_PCA_UMAP_by_sample.png`
- `figures/面板原件/S5C_2x2_stratified_violin.png`
- `figures/面板原件/S5C_2x2_stratified_violin_Q3corrected.png`
- `figures/面板原件/fig1_UMAP_dual.png`
- `figures/面板原件/fig2_violin_HMGB1_HAVCR2_MAFB.png`
- `figures/面板原件/fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected.png`
- `figures/面板原件/fig2_violin_HMGB1_HAVCR2_MAFB_v2.png`
- `figures/面板原件/fig3_volcano_CD68_DLBCL_vs_Normal.png`
- `figures/面板原件/fig4_scatter_HMGB1_HAVCR2_LATAM.png`
- `figures/面板原件/fig4_scatter_v3_16gene_LATAM.png`
- `figures/面板原件/fig5_GSEA_joyplot.png`
- `figures/面板原件/fig5_GSEA_joyplot_v2.png`
- `figures/面板原件/fig7_spatial_coevolution_split_heatmap.png`
- `figures/面板原件/fig9_violin_HMGB1_TGS_Normal_vs_DLBCL.png`
- `figures/面板原件/task4_CD68_volcano_DLBCL_vs_Normal.png`
- `figures/面板原件/task5_dual_chamber_heatmap.png`
- `figures/面板原件/task6_GSEA_ridge_plot.png`
