# MAFB预后评分模型

> **论文对应**：Results 6 — MAFB 预后评分模型  
> **图件**：Figure 6（a–g）、Figure S6、Table S5C / S19 / S20 / S23 / S24

## 分析内容

bulk 多队列（GSE10846/GSE87371/GSE11318/GSE181063/TCGA）：LASSO-Cox 建模 → 风险评分 → KM/TimeROC/校准 → CIBERSORT 免疫浸润

## 目录

```
section6_MAFB预后评分模型/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/07_bulk_prognosis/07_01_prognosis_main.R` | |
| `code/07_bulk_prognosis/07_02_lasso_cv_curves.R` | |
| `code/07_bulk_prognosis/07_03_immune_deconvolution.R` | |
| `code/07_bulk_prognosis/07_04_supplementary_tables.R` | |
| `code/07_bulk_prognosis/07_05_cibersort.R` | |
| `code/10_figures/supp_tables/S24_BtoMac_ranking.py` | |
| `code/10_figures/supp_tables/mono_5gene_vs_16gene_score.py` | |
| `code/10_figures/supp_tables/tableS20.R` | |
| `code/10_figures/supp_tables/tableS23_ARI.py` | |
| `code/10_figures/supp_tables/tableS5C_GRN_5000.py` | |
| `code/legacy/yuhou.R` | 归档件（不进入流程）：25 基因版原始母脚本，**产出 07_01 所读、但本仓库原先无处产出的中间表**（Cox 表 / 校准数据 / 免疫浸润 / 多队列与 TCGA 风险评分）—— 见 `code/legacy/README.md` |
| `code/legacy/yuhou2.R` | 归档件（不进入流程）：23 基因最终模型平行支线，含论文未使用的 NCICCR-DLBCL 队列 |

## 输入数据（公用，位于 `data/`）

- `GSE10846_expression_gene_level.csv`
- `DLBCL_prognosis/GSE10846_sur_model.Rdata`
- `DLBCL_prognosis/tables/`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 6
```

## 环境状态

- 需要：R: survival, glmnet, timeROC, survminer, CIBERSORT
- 本机 `scRNA` 环境已具备：（无）
- 尚缺：R: survival, glmnet, timeROC, survminer, CIBERSORT

## 产出图件

- `figures/论文成图/Fig6.png`
- `figures/论文成图/FigS6.png`
- `figures/论文成图/Figure S6.png`
- `figures/论文成图/Figure6.pdf`
- `figures/论文成图/Figure6.png`
- `figures/重排版/Fig6.png`
- `figures/面板原件/B1_RiskScore_tripanel.png`
- `figures/面板原件/B3_TimeROC.png`
- `figures/面板原件/B5_Calibration_ggplot.png`
- `figures/面板原件/Fig1A_LASSO_CV_curve.png`
- `figures/面板原件/Fig1B_LASSO_coef_trajectory.png`
- `figures/面板原件/Fig1B_LASSO_coef_trajectory2.png`
- `figures/面板原件/Fig2_Cox_forest_plot.png`
- `figures/面板原件/Fig3A_Nomogram.png`
- `figures/面板原件/Fig3B_Calibration_curves.png`
- `figures/面板原件/Fig4_Immune_infiltration.png`
- `figures/面板原件/Fig5_MultiCohort_KM_curves.png`
- `figures/面板原件/Figure S6.1.png`
- `figures/面板原件/RADAR.png`
- `figures/面板原件/S6.png`
- `figures/面板原件/TIDE.png`
- `figures/面板原件/checkpoints (2).png`
- `figures/面板原件/checkpoints.png`
- `figures/面板原件/ciber.png`
