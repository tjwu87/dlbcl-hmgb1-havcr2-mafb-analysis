# MAFB prognostic scoring model

> **Paper mapping**: Results 6 — MAFB prognostic scoring model  
> **Figures**: Figure 6 (a–g), Figure S6, Table S5C / S19 / S20 / S23 / S24

## What this part does

Bulk multi-cohort (GSE10846/GSE87371/GSE11318/GSE181063/TCGA): LASSO-Cox modelling → risk score → KM/TimeROC/calibration → CIBERSORT immune infiltration

## Layout

```
section6_MAFB预后评分模型/
├── code/            analysis scripts (see below)
├── data/            intermediates specific to this part (for comparison / reuse)
├── figures/         paper figures and panels
└── README.md
```

## Scripts

| Script | Notes |
| --- | --- |
| `code/07_bulk_prognosis/07_01_prognosis_main.R` |  |
| `code/07_bulk_prognosis/07_02_lasso_cv_curves.R` |  |
| `code/07_bulk_prognosis/07_03_immune_deconvolution.R` |  |
| `code/07_bulk_prognosis/07_04_supplementary_tables.R` |  |
| `code/07_bulk_prognosis/07_05_cibersort.R` |  |
| `code/10_figures/supp_tables/S24_BtoMac_ranking.py` |  |
| `code/10_figures/supp_tables/mono_5gene_vs_16gene_score.py` |  |
| `code/10_figures/supp_tables/tableS20.R` |  |
| `code/10_figures/supp_tables/tableS23_ARI.py` |  |
| `code/10_figures/supp_tables/tableS5C_GRN_5000.py` |  |
| `code/legacy/yuhou.R` | Archived item (not part of the pipeline): original 25-gene master script, **which produces the intermediate tables read by `07_01` that previously had no producer in this repository** (Cox table / calibration data / immune infiltration / multi-cohort and TCGA risk scores) — see `code/legacy/README.md` |
| `code/legacy/yuhou2.R` | Archived item (not part of the pipeline): parallel branch of the 23-gene final model, containing the NCICCR-DLBCL cohort that the paper does not use |

## Input data (shared, under `data/`)

- `GSE10846_expression_gene_level.csv`
- `DLBCL_prognosis/GSE10846_sur_model.Rdata`
- `DLBCL_prognosis/tables/`

## How to run

```bash
# No data root needed: config/paths.py resolves to data/ by default
cd DLBCL_HMGB1_HAVCR2_MAFB
# whole pipeline
python run_all.py
# this part only
python run_all.py --stage 6
```

## Environment status

- Required: R: survival, glmnet, timeROC, survminer, CIBERSORT
- Already available in the local `scRNA` environment: (none)
- Still missing: R: survival, glmnet, timeROC, survminer, CIBERSORT

## Output figures

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
