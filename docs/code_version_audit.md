# 代码版本审计报告

生成时间：2026-09-11 01:17　（工具：`tools/audit_code_versions.py`）

**判读**：「产物命中」= 脚本声明写出的文件在项目目录中被找到。
命中数高且产物时间新 → 该版本产出了实际结果（大概率是论文所用版本）。
命中数为 0 → 早期探索脚本，或输出位于服务器（`/mnt/results`）未在本次拷贝范围内。

| 脚本 | 脚本时间 | 声明产物 | 命中 | 最新产物时间 | 产物示例 |
|---|---|---:|---:|---|---|
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09_figures_complete.py` | 1980-01-01 | 48 | **46** | 2026-04-03 | fig01_mafb_quiver_vector_field.png、fig01_mafb_quiver_vector_field.svg、fig02_mafb_streamline_vector_field.png |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09C_figures_13_23.py` | 1980-01-01 | 22 | **20** | 2026-04-03 | fig13_mafb_net_velocity_ko_oe.png、fig13_mafb_net_velocity_ko_oe.svg、fig14_mafb_mechanism_summary.png |
| `yuhou.R` | 2026-03-29 | 14 | **13** | 2026-03-30 | All_cohorts_risk_scores.csv、Calibration_data_1_3_5yr.csv、Fig1A_LASSO_CV_curve.png |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09B_figures_07_12.py` | 1980-01-01 | 12 | **12** | 2026-04-03 | fig07_mafb_enrichment_bubble.png、fig07_mafb_enrichment_bubble.svg、fig08_mafb_regulatory_network.png |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09A_figures_01_06.py` | 1980-01-01 | 12 | **12** | 2026-04-03 | fig01_mafb_quiver_vector_field.png、fig01_mafb_quiver_vector_field.svg、fig02_mafb_streamline_vector_field.png |
| `yuhou2.R` | 2026-03-30 | 13 | **9** | 1980-01-01 | Calibration_data_1_3_5yr.csv、Fig1A_LASSO_CV_curve.png、Fig1B_LASSO_coef_trajectory.png |
| `mono_mac_5subtypes.py` | 2026-04-06 | 8 | **8** | 2026-04-06 | monomac_paga_trajectory.png、monomac_paga_trajectory.svg、monomac_paga_trajectory_rawpt.png |
| `FigS2补图.py` | 2026-09-09 | 7 | **7** | 2026-09-09 | S8_MAFB_CD163_check.csv、S8_leave_one_patient_out.csv、S8_perpatient.csv |
| `tableS23_ARI.py` | 2026-09-10 | 6 | **6** | 2026-09-10 | ARI_crosstab_CellType_vs_celltype.csv、ARI_crosstab_CellType_vs_leiden06.csv、TableS23a_crosstab_CellType_vs_leiden06.csv |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09_export_data.py` | 1980-01-01 | 5 | **5** | 1980-01-01 | cell_fate_scores_MAFB_KO.csv、cell_metadata.csv、imputed_count_WT.csv |
| `FigS8.py` | 2026-09-08 | 4 | **4** | 2026-09-09 | S8_MAFB_CD163_check.csv、S8_sensitivity.csv、S8_signature_filtering.csv |
| `PAOTU.py` | 2026-04-21 | 4 | **4** | 2026-04-20 | mac_m1_m2_violin2.png、mac_subtype_dotplot2.png、mac_subtype_umap2.png |
| `S24.py` | 2026-09-10 | 3 | **3** | 2026-09-10 | TableS24_BtoMac_ranking_sensitivity.csv、TableS24_clean.csv、liana_top20_under_0.1.csv |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09_fig23_deg_ora.py` | 1980-01-01 | 2 | **2** | 2026-04-03 | fig23_mafb_ko_ora_enrichment_bubble.png、fig23_mafb_ko_ora_enrichment_bubble.svg |
| `tableS5C_GRN_5000.py` | 2026-09-10 | 1 | **1** | 2026-09-10 | TableS5f_GRN_edge_cap_concordance.csv |
| `tableS19.py` | 2026-09-06 | 6 | **1** | 2026-09-06 | FigS3_receptor_prioritization.png |
| `10846_LASSO.r` | 2026-03-13 | 0 | **0** | — | — |
| `10846lasso-2.r` | 2026-03-30 | 0 | **0** | — | — |
| `4SUPP.py` | 2026-03-20 | 0 | **0** | — | — |
| `All_in_one0205.rmd` | 2026-02-08 | 0 | **0** | — | — |
| `All_in_one1124.Rmd` | 2026-01-01 | 9 | **0** | — | — |
| `B_cell_liana.py` | 2026-03-24 | 0 | **0** | — | — |
| `Bottom_left_axis.R` | 2025-05-14 | 0 | **0** | — | — |
| `GDSC_11_Complexes\02_resistance_E98D\pymol_WT_vs_E98D_overlay.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\02_resistance_E98D\pymol_WT_vs_E98D_overlay2.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\03_V75M\pymol_V75M_steric_occupancy.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\DABRAFENIB\4.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\DABRAFENIB\fir.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\DABRAFENIB\sec.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\DABRAFENIB\thir.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\2.py` | 2026-03-10 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\3.py` | 2026-03-10 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\boltz.py` | 2026-03-10 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\complex.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\fir.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\make_publication_figures.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\make_publication_figures1.py` | 2026-03-09 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\make_publication_figures2.py` | 2026-03-09 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\make_publication_figures3.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\render_docking.py` | 2026-03-10 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\sec.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\vina.py` | 2026-03-27 | 0 | **0** | — | — |
| `GDSC_11_Complexes\crizotinib\visualize_complex.py` | 2026-03-10 | 0 | **0** | — | — |
| `GDSC_11_Complexes\三个分子对接及耐药位点预测\02_resistance_E98D\pymol_WT_vs_E98D_overlay.py` | 1980-01-01 | 0 | **0** | — | — |
| `GDSC_11_Complexes\三个分子对接及耐药位点预测\03_V75M\pymol_V75M_steric_occupancy.py` | 2026-03-27 | 0 | **0** | — | — |
| `GSE182434_singlecell.py` | 2026-03-24 | 0 | **0** | — | — |
| `GSE232853.py` | 2026-04-07 | 2 | **0** | — | — |
| `PAGA.py` | 2026-03-21 | 0 | **0** | — | — |
| `WODEYUHOU.rMD` | 2026-03-29 | 0 | **0** | — | — |
| `boltz-2.py` | 2026-04-03 | 7 | **0** | — | — |
| `bulk-download\03_V75M\pymol_V75M_steric_occupancy.py` | 2026-03-27 | 0 | **0** | — | — |
| `bulk-download\DLBCL_prognosis\tables\LASSO.R` | 2026-03-29 | 0 | **0** | — | — |
| `bulk-download\DLBCL_prognosis\tables\LASSO复现.r` | 2026-03-30 | 0 | **0** | — | — |
| `bulk-download\DLBCL_prognosis\tables\bupao.r` | 2026-09-06 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step00A_env.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step00B_freeze_params.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step01_preprocess.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step02AB_oracle_build.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step02C_get_links.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step02DE_filter_audit.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step02F_fit_grn.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step03_null.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step03_perturbation.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step04A_fate_calc.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step04B_fate_stats.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step04_5_consistency_check.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step04_scoring_figures.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step05_scenic_targets.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step06_program_score.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step07C_ko_response_enrichment.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step07_enrichment.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step08AB_sensitivity.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step08C_sensitivity_hvg4000.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step08D_sensitivity_n30.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step08E_k_sensitivity.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step08_summary.py` | 1980-01-01 | 0 | **0** | — | — |
| `bulk-download\GSE182434\celloracle_paper_final_20260331_085223\scripts\step09_figures.py` | 1980-01-01 | 0 | **0** | — | — |
| `celloracle0324.py` | 2026-03-24 | 0 | **0** | — | — |
| `celloracle0331.py` | 2026-09-10 | 2 | **0** | — | — |
| `cibersort.r` | 2026-03-30 | 0 | **0** | — | — |
| `liana.py` | 2026-04-04 | 0 | **0** | — | — |
| `mycolors.R` | 2025-05-14 | 0 | **0** | — | — |
| `paga2.py` | 2026-04-06 | 0 | **0** | — | — |
| `pyscenic-suoxiao.py` | 2026-03-24 | 0 | **0** | — | — |
| `pyscenic.py` | 2026-03-22 | 0 | **0** | — | — |
| `tableS20.R` | 2026-09-06 | 0 | **0** | — | — |
| `xibaotongxunchayihuoshantu.py` | 2026-02-27 | 0 | **0** | — | — |
| `生图\03_V75M\pymol_V75M_steric_occupancy.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\DABRAFENIB\4.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\DABRAFENIB\fir.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\DABRAFENIB\sec.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\DABRAFENIB\thir.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\2.py` | 2026-03-10 | 0 | **0** | — | — |
| `生图\crizotinib\3.py` | 2026-03-10 | 0 | **0** | — | — |
| `生图\crizotinib\boltz.py` | 2026-03-10 | 0 | **0** | — | — |
| `生图\crizotinib\complex.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\fir.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\make_publication_figures.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\make_publication_figures1.py` | 2026-03-09 | 0 | **0** | — | — |
| `生图\crizotinib\make_publication_figures2.py` | 2026-03-09 | 0 | **0** | — | — |
| `生图\crizotinib\make_publication_figures3.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\render_docking.py` | 2026-03-10 | 0 | **0** | — | — |
| `生图\crizotinib\sec.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\vina.py` | 2026-03-27 | 0 | **0** | — | — |
| `生图\crizotinib\visualize_complex.py` | 2026-03-10 | 0 | **0** | — | — |
| `髓系细胞5基因_16基因打分对比.py` | 2026-09-07 | 0 | **0** | — | — |

---

## 同名/近名脚本分组（需人工确认最终版）

### ``

- `GDSC_11_Complexes\DABRAFENIB\4.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\2.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\3.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `生图\DABRAFENIB\4.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\2.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `生图\crizotinib\3.py`　脚本 2026-03-10　命中 0/0　最新产物 —

### `all_in_one`

- `All_in_one0205.rmd`　脚本 2026-02-08　命中 0/0　最新产物 —
- `All_in_one1124.Rmd`　脚本 2026-01-01　命中 0/9　最新产物 —

### `boltz`

- `GDSC_11_Complexes\crizotinib\boltz.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `boltz-2.py`　脚本 2026-04-03　命中 0/7　最新产物 —
- `生图\crizotinib\boltz.py`　脚本 2026-03-10　命中 0/0　最新产物 —

### `celloracle`

- `celloracle0324.py`　脚本 2026-03-24　命中 0/0　最新产物 —
- `celloracle0331.py`　脚本 2026-09-10　命中 0/2　最新产物 —

### `complex`

- `GDSC_11_Complexes\crizotinib\complex.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\complex.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `fir`

- `GDSC_11_Complexes\DABRAFENIB\fir.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\fir.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\DABRAFENIB\fir.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\fir.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `make_publication_figures`

- `GDSC_11_Complexes\crizotinib\make_publication_figures.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\make_publication_figures1.py`　脚本 2026-03-09　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\make_publication_figures2.py`　脚本 2026-03-09　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\make_publication_figures3.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\make_publication_figures.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\make_publication_figures1.py`　脚本 2026-03-09　命中 0/0　最新产物 —
- `生图\crizotinib\make_publication_figures2.py`　脚本 2026-03-09　命中 0/0　最新产物 —
- `生图\crizotinib\make_publication_figures3.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `paga`

- `PAGA.py`　脚本 2026-03-21　命中 0/0　最新产物 —
- `paga2.py`　脚本 2026-04-06　命中 0/0　最新产物 —

### `pymol_5m_steric_occupancy`

- `GDSC_11_Complexes\03_V75M\pymol_V75M_steric_occupancy.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\三个分子对接及耐药位点预测\03_V75M\pymol_V75M_steric_occupancy.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `bulk-download\03_V75M\pymol_V75M_steric_occupancy.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\03_V75M\pymol_V75M_steric_occupancy.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `pymol_wt_vs_e98d_overlay`

- `GDSC_11_Complexes\02_resistance_E98D\pymol_WT_vs_E98D_overlay.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\02_resistance_E98D\pymol_WT_vs_E98D_overlay2.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\三个分子对接及耐药位点预测\02_resistance_E98D\pymol_WT_vs_E98D_overlay.py`　脚本 1980-01-01　命中 0/0　最新产物 —

### `render_docking`

- `GDSC_11_Complexes\crizotinib\render_docking.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `生图\crizotinib\render_docking.py`　脚本 2026-03-10　命中 0/0　最新产物 —

### `sec`

- `GDSC_11_Complexes\DABRAFENIB\sec.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `GDSC_11_Complexes\crizotinib\sec.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\DABRAFENIB\sec.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\sec.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `tables`

- `tableS19.py`　脚本 2026-09-06　命中 1/6　最新产物 2026-09-06
- `tableS20.R`　脚本 2026-09-06　命中 0/0　最新产物 —

### `thir`

- `GDSC_11_Complexes\DABRAFENIB\thir.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\DABRAFENIB\thir.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `vina`

- `GDSC_11_Complexes\crizotinib\vina.py`　脚本 2026-03-27　命中 0/0　最新产物 —
- `生图\crizotinib\vina.py`　脚本 2026-03-27　命中 0/0　最新产物 —

### `visualize_complex`

- `GDSC_11_Complexes\crizotinib\visualize_complex.py`　脚本 2026-03-10　命中 0/0　最新产物 —
- `生图\crizotinib\visualize_complex.py`　脚本 2026-03-10　命中 0/0　最新产物 —

### `yuhou`

- `yuhou.R`　脚本 2026-03-29　命中 13/14　最新产物 2026-03-30
- `yuhou2.R`　脚本 2026-03-30　命中 9/13　最新产物 1980-01-01
