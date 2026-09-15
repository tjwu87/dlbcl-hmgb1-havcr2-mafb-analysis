# MAFB计算机扰动

> **论文对应**：Results 4 — MAFB 计算机扰动（CellOracle）  
> **图件**：Figure 4（a–h）

## 分析内容

CellOracle：base GRN → oracle 构建 → links 过滤 → GRN 拟合 → MAFB KO/OE/Null 扰动 → 命运概率 → 敏感性分析

## 目录

```
section4_MAFB计算机扰动/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/05_CellOracle/05_90_reproduce_figures_local.py` | |
| `code/05_CellOracle/05_91_grn_edge_sensitivity.py` | |
| `code/05_CellOracle/step00A_env.py` | |
| `code/05_CellOracle/step00B_freeze_params.py` | |
| `code/05_CellOracle/step01_preprocess.py` | |
| `code/05_CellOracle/step02AB_oracle_build.py` | |
| `code/05_CellOracle/step02C_get_links.py` | |
| `code/05_CellOracle/step02DE_filter_audit.py` | |
| `code/05_CellOracle/step02F_fit_grn.py` | |
| `code/05_CellOracle/step03_null.py` | |
| `code/05_CellOracle/step03_perturbation.py` | |
| `code/05_CellOracle/step04A_fate_calc.py` | |
| `code/05_CellOracle/step04B_fate_stats.py` | |
| `code/05_CellOracle/step04_5_consistency_check.py` | |
| `code/05_CellOracle/step04_scoring_figures.py` | |
| `code/05_CellOracle/step05_scenic_targets.py` | |
| `code/05_CellOracle/step06_program_score.py` | |
| `code/05_CellOracle/step07C_ko_response_enrichment.py` | |
| `code/05_CellOracle/step07_enrichment.py` | |
| `code/05_CellOracle/step08AB_sensitivity.py` | |
| `code/05_CellOracle/step08C_sensitivity_hvg4000.py` | |
| `code/05_CellOracle/step08D_sensitivity_n30.py` | |
| `code/05_CellOracle/step08E_k_sensitivity.py` | |
| `code/05_CellOracle/step08_summary.py` | |
| `code/05_CellOracle/step09_export_data.py` | |
| `code/05_CellOracle/step09_fig23_deg_ora.py` | |

## 输入数据（公用，位于 `data/`）

- `celloracle0331/oracle_WT.pkl`
- `celloracle0331/oracle_MAFB_KO.pkl`
- `celloracle0331/oracle_MAFB_OE.pkl`
- `celloracle0331/oracle_Null.pkl`
- `GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv`
- `GSE182434/celloracle_paper_final_20260331_085223/`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 4
```

## 环境状态

- 需要：celloracle, pyscenic, ctxcore
- 本机 `scRNA` 环境已具备：（无）
- 尚缺：celloracle, pyscenic, ctxcore

## 产出图件

- `figures/论文成图/Fig4.png`
- `figures/论文成图/FigS4.png`
- `figures/论文成图/Figure4.pdf`
- `figures/论文成图/Figure4.png`
- `figures/论文成图/FigureS4.png`
- `figures/重排版/Fig4.png`
- `figures/面板原件/Fig4new.png`
- `figures/面板原件/FigA_MAFB_KO_mirror_bubble_plot.png`
- `figures/面板原件/FigB_MAFB_TF_target_network.png`
- `figures/面板原件/FigB_streamline_WT_vs_KO.png`
- `figures/面板原件/FigC_perturbation_density_map.png`
- `figures/面板原件/FigVF_OE_B_streamline_WT_vs_OE.png`
- `figures/面板原件/FigVF_WT_KO_OE_vector_fields.png`
- `figures/面板原件/Figure4.1.png`
- `figures/面板原件/Figure4.2.png`
- `figures/面板原件/Figure4.3.png`
- `figures/面板原件/Figure4.5.png`
- `figures/面板原件/Figure4.6.png`
- `figures/面板原件/GSE182434_celloracle_fig01_mafb_quiver_vector_field (2).png`
- `figures/面板原件/GSE182434_celloracle_fig02_mafb_streamline_vector_field (1).png`
- `figures/面板原件/GSE182434_celloracle_fig03_mafb_differential_vector_field (2).png`
- `figures/面板原件/GSE182434_celloracle_fig04_mafb_velocity_density_per_subtype (1).png`
- `figures/面板原件/GSE182434_celloracle_fig05_mafb_target_violin_plots (1).png`
- `figures/面板原件/GSE182434_celloracle_fig06_mafb_butterfly_ko_oe (1).png`
- `figures/面板原件/GSE182434_celloracle_fig07_mafb_enrichment_bubble (1).png`
- `figures/面板原件/GSE182434_celloracle_fig08_mafb_regulatory_network (1).png`
- `figures/面板原件/GSE182434_celloracle_fig09_pseudotime_mafb_expression (1).png`
- `figures/面板原件/GSE182434_celloracle_fig10_mafb_fate_probability_matrix (1).png`
- `figures/面板原件/GSE182434_celloracle_fig11_mafb_target_dotplot (1).png`
- `figures/面板原件/GSE182434_celloracle_fig12_mafb_ko_waterfall (1).png`
- `figures/面板原件/GSE182434_celloracle_fig13_mafb_net_velocity_ko_oe (1).png`
- `figures/面板原件/GSE182434_celloracle_fig14_mafb_mechanism_summary (1).png`
- `figures/面板原件/S4.png`
- `figures/面板原件/celloracle_MAFB_KO_fate_probability_plot.png`
- `figures/面板原件/celloracle_MAFB_KO_p05_sig_violin.png`
- `figures/面板原件/celloracle_MAFB_KO_projection_violin.png`
- `figures/面板原件/fig01_mafb_quiver_vector_field.png`
- `figures/面板原件/fig02_mafb_streamline_vector_field.png`
- `figures/面板原件/fig03_mafb_differential_vector_field.png`
- `figures/面板原件/fig04_mafb_velocity_density_per_subtype.png`
- `figures/面板原件/fig05_mafb_target_violin_plots.png`
- `figures/面板原件/fig08_mafb_regulatory_network_FINAL_UPDATED.png`
- `figures/面板原件/fig09_mafb_expression_umap.png`
- `figures/面板原件/fig10_fate_score_violin_wt_vs_ko.png`
- `figures/面板原件/fig12_mafb_ko_waterfall.png`
- `figures/面板原件/fig15_program_score_violin.png`
- `figures/面板原件/fig16_program_score_paired_violin.png`
- `figures/面板原件/fig18_grn_coefficient_heatmap.png`
- `figures/面板原件/fig_effect_size_summary_dotplot.png`
- `figures/面板原件/fig_horizontal_heatmap_all_deg.png`
- `figures/面板原件/fig_merged_ora_bubble_all_cells.png`