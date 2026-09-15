# 主图 ↔ 子图 ↔ 代码 权威映射表

> 判定依据（三重交叉验证，非猜测）：
> 1. **直接阅读原始主图**（助手已可读图），逐 panel 记录内容；
> 2. **panel 文字导出**（`docs/panel_text_dump.md`，从 SVG 提取可读文字）；
> 3. **渲染核对**——把候选 SVG 渲染成图，与原图逐 panel 目视比对（`tools/panel_contact_sheet.py`）。
>
> 置信度：★★★ = 渲染后目视确认；★★ = 文件名/文字高度吻合；★ = 按分析流程推断。

生成时间：2026-09-11

---

## Fig1 — 单细胞图谱 + 细胞通讯（6 panel）

| panel | 内容 | 产出文件（仓库内相对路径） | 产出脚本 | 置信 |
|---|---|---|---|---|
| **1a** | UMAP「Refined Cell Types」+「Tissue (DLBCL vs Tonsil)」+ 类型图例 | `10_00_all_panels_GSE182434/umap_celltypes2.svg` | `10_figures/10_00_all_panels_GSE182434.py` | ★★★ |
| **1b** | Marker 点图：14 细胞类型 × 基因（色=Fraction of cells，径=Mean expression） | `01_01_qc_integration/marker_dotplot.svg` | `01_scRNA_GSE182434/01_01_qc_integration.py` | ★★★ |
| **1c** | inferCNV 推断拷贝数热图（行分 Normal B / Normal B (DLBCL) / Malignant B） | `10_00_all_panels_GSE182434/chromosome_heatmap.svg` | `10_figures/10_00_all_panels_GSE182434.py` | ★★★ |
| **1d** | 细胞互作环形网络（Malignant vs Normal B） | `02_02_liana_visualization/circle_plot_malignant_vs_normal_v2.svg` | `02_cellcomm_LIANA/02_02_liana_visualization.py` | ★★★ |
| **1e** | 互作强度热图（source × target，14 细胞类型） | `02_02_liana_visualization/interaction_heatmap.svg` | 同上 | ★★★ |
| **1f** | 差异互作散点（ΔMagnitude vs −log10 p；Malignant 82 / Normal 29） | `02_02_liana_visualization/volcano_malignant_vs_normal_monomac.svg` | 同上 | ★★★ |

---

## Fig2 — 髓系亚群 + 拟时序（6 panel）

| panel | 内容 | 产出文件 | 产出脚本 | 置信 |
|---|---|---|---|---|
| **2a** | UMAP 髓系亚型（Mono / DC_1 / LA_TAM / IFN_TAM / DC_2） | `01_02_myeloid_subtypes/*`（子图名待渲染核对） | `01_scRNA_GSE182434/01_02_myeloid_subtypes.py` | ★★ |
| **2b** | 点图：5 亚型 × 基因（FCN1/S100A8/C1QA/CD74/CCL18/MARCO/TREM2/SPP1…） | 同上 | 同上 | ★★ |
| **2c** | 细胞互作网络：1 张全亚型环形 + 4 张以 Malignant B / LA_TAM / IFN_TAM / Normal B 为中心的辐射图 | `02_02_liana_visualization/circle_plot_13subtypes.svg` + `circle_*_sender.svg`（4 张） | `02_cellcomm_LIANA/02_02_liana_visualization.py` | ★★★（环形）/ ★★（辐射） |
| **2d** | 2 张 UMAP（亚型 / TF differentiation），带 PT 标注 | `03_01_paga_trajectory/monomac_paga_trajectory.svg`、`monomac_trajectory_directed.svg` | `03_trajectory_PAGA/03_01_paga_trajectory.py` | ★★ |
| **2e** | 拟时序热图（Mono/IFN_TAM/LA_TAM/DC_2/DC_1 × 基因） | `03_01_paga_trajectory/monomac_gene_dynamics.svg` 或 `scenic_tf_heatmap.svg` | 同上 | ★ |
| **2f** | 4 张富集条形图：KEGG/GO-BP × IFN_TAM/LA_TAM up | `03_01_paga_trajectory/enrichment_barplot_IFN_TAM_vs_LA_TAM.svg` | 同上 | ★★ |

---

## Fig3 — 调控子（regulon）分析（7 panel）

| panel | 内容 | 产出文件 | 产出脚本 | 置信 |
|---|---|---|---|---|
| **3a** | 热图：受体基因（按峰值表达排序）× 细胞（Mono/IFN_TAM/LA_TAM），标 ▲正相关 ▼负相关 | `04_01_pyscenic_regulons/receptor_pseudotime_heatmap_path.svg` | `04_pySCENIC/04_01_pyscenic_regulons.py`（亦见于 `03_01`） | ★★★ |
| **3b** | 4 折线图：CD74 (r=0.418) / AXL (0.277) / HLA-DPA1 (0.218) / HAVCR2 (0.232) vs 拟时序 | `04_01_pyscenic_regulons/receptor_pseudotime_lineplots_path.svg` | 同上 | ★★★ |
| **3c** | 热图：regulon 活性（MAFB(+)·IRF1(+)·STAT1(+)…）× Mono/IFN_TAM/LA_TAM | `04_01_pyscenic_regulons/Fig1_regulon_activity_heatmap_paper.svg` | 同上 | ★★★ |
| **3d** | 散点「Regulon Specificity Score」：RSS-IFN_TAM vs RSS-LA_TAM（LA_TAM 特异 22 / IFN_TAM 特异 7 / 共享 10） | `04_01_pyscenic_regulons/Fig2_RSS_scatter_paper.svg` | 同上 | ★★★ |
| **3e** | Sankey（7 层）：Malignant B Cell → Secreted Ligand → Mono/Mac Receptor → Signaling Pathway → Core TF → Top5 Target → TAM Subtype | `04_01_pyscenic_regulons/sankey_7layer_TF_targets_v3_paper.svg` | 同上 | ★★★ |
| **3f** | 4×2 折线图：ETV5/MITF/IRF8/IRF7 vs MAFB_KO；ETV5/BHLHE40/MAFB/IRF8 vs AXL | `04_01_pyscenic_regulons/lineplot_LA_TAM_TF_receptor_2x4_paper.svg` | 同上 | ★★★ |
| **3g** | 富集点图：GO BP / KEGG / Reactome（MAFB regulon targets） | `04_01_pyscenic_regulons/MAFB_gseapy_enrichment_dotplot_paper.svg` | 同上 | ★★★ |

**未使用**：`Fig3_TF_target_network_paper.svg`、`pseudotime_TF_receptor_dynamics_paper.svg`（原图未采用，可作补充）。

---

## Fig4 — CellOracle 扰动（8 panel）★ 全部确认

| panel | 内容 | 产出文件 | 置信 |
|---|---|---|---|
| **4a** | 3 张 UMAP：Cell Subtypes / MAFB Expression / MAFB KO – Velocity Magnitude | `05_90_reproduce_figures_local/fig09_mafb_expression_umap.svg` | ★★★ |
| **4b** | 3 小提琴：Mono / IFN_TAM / LA_TAM 的 LA_TAM fate score，WT vs MAFB KO | `05_90_reproduce_figures_local/fig10_fate_score_violin_wt_vs_ko.svg` | ★★★ |
| **4c** | 3 张 UMAP：WT (Baseline) / MAFB KO / Null (Randomized GRN) | `05_90_reproduce_figures_local/fig01_mafb_quiver_vector_field.svg` | ★★★ |
| **4d** | MAFB 调控网络（Top15 靶基因 + 共调控因子） | `05_90_reproduce_figures_local/fig08_mafb_regulatory_network_FINAL_LOGIC_FIXED.svg` | ★★★ |
| **4e** | 热图：All Significant Down-regulated Genes (MAFB KO vs WT) | `05_90_reproduce_figures_local/fig_horizontal_heatmap_all_deg.svg` | ★★★ |
| **4f** | 3 小提琴：MAFB KO Effect on LA_TAM Program Score（Δ=−0.036/−0.167/−0.000） | `05_90_reproduce_figures_local/fig16_program_score_paired_violin.svg` | ★★★ |
| **4g** | 点图：Down-regulated genes（GO BP / KEGG / MSigDB Hallmark） | `05_90_reproduce_figures_local/fig_merged_ora_bubble_all_cells.svg` | ★★★ |
| **4h** | 散点：MAFB KO Effect-Size（Fate/Program score + SNR） | `05_90_reproduce_figures_local/fig_effect_size_summary_dotplot.svg` | ★★★ |

> 产出脚本：`05_CellOracle/05_90_reproduce_figures_local.py`

---

## Fig5 — 空间转录组验证（6 panel）★ 全部确认

| panel | 内容 | 产出文件 | 置信 |
|---|---|---|---|
| **5a** | 2 张 UMAP：Cell Compartment (CD20/CD68)、Tissue Type (DLBCL/Normal) | `06_01_spatial_analysis/fig1_UMAP_dual.svg` | ★★★ |
| **5b** | 火山图：CD68⁺ Macrophages DLBCL vs Normal（标出 MAFB / HAVCR2 / HMGB1） | `06_01_spatial_analysis/fig3_volcano_DEA_CD68.svg` | ★★★ |
| **5c** | 3 小提琴：HMGB1 / HAVCR2 / MAFB（CD20 vs CD68） | `06_01_spatial_analysis/fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected.svg` | ★★★ |
| **5d** | 双腔热图：CD20⁺ vs CD68⁺ 配对 ROI（DLBCL n=114 / Normal n=114） | `06_01_spatial_analysis/task5_dual_chamber_heatmap.svg` | ★★★ |
| **5e** | 2 散点：B-cell HMGB1 vs Macrophage HAVCR2 / vs LA_TAM Score | `06_01_spatial_analysis/fig4_scatter_HMGB1_HAVCR2_LATAM.svg` | ★★★ |
| **5f** | GSEA ridge plot（Depleted ↔ Enriched in LA_TAM_Rich，18 条 Hallmark） | `06_01_spatial_analysis/fig5_GSEA_joyplot.svg` | ★★★ |

> 产出脚本：`06_spatial_GSE232853/06_01_spatial_analysis.py`

---

## Fig6 — 预后模型（8 panel）

| panel | 内容 | 产出文件 | 产出脚本 | 置信 |
|---|---|---|---|---|
| **6a** | LASSO 系数路径 + 非零系数个数 | `07_02_lasso_cv_curves/Fig1B_LASSO_coef_trajectory.pdf`、`Fig1A_LASSO_CV_curve.pdf` | `07_bulk_prognosis/07_02_lasso_cv_curves.R` | ★★★ |
| **6b** | 6 张 KM 曲线（5 队列验证） | `07_01_prognosis_main/A5_MultiCohort_KM.pdf` | `07_bulk_prognosis/07_01_prognosis_main.R` | ★★★ |
| **6c** | 森林图：Univariable / Multivariable Cox | `07_01_prognosis_main/A1_Cox_forest.pdf` | 同上 | ★★★ |
| **6d** | 校准曲线：1/3/5 年 OS | `07_01_prognosis_main/A3_Calibration_curves.pdf` | 同上 | ★★★ |
| **6e** | 雷达图：22 种免疫细胞浸润 | `07_01_prognosis_main/A4_Immune_infiltration.pdf` | 同上 | ★★ |
| **6f** | 4 散点（H2/M1 Macrophage、CD8 T、B Cell vs Risk Score）+ Spearman rho 条形图 | 待确认（候选：`A4_Immune_infiltration.pdf` 内含；或 `07_03_immune_deconvolution.R` 产出） | `07_bulk_prognosis/07_03_immune_deconvolution.R` | ★ |
| **6g** | 免疫检查点小提琴（20 基因） | 待确认 | 同上 | ★ |
| **6h** | ICB Response Rate + TIDE / T-cell Exclusion | 待确认 | 同上 | ★ |

> ⚠ `A2_Nomogram.pdf` 不在 Fig6 中原图未见 → 应为**补充图**。

---

## Fig7 — 药物敏感 + 分子对接（10 panel）

| panel | 内容 | 产出文件 | 产出脚本 | 置信 |
|---|---|---|---|---|
| **7a** | Tim-3 Docking Score vs GDSC2 IC50 散点 + Top-10 候选药条形 | `08_01_ic50_and_boltz2/top200_gdsc2_matched.svg`、`Figure_IC50_Raincloud_HAVCR2_docking.svg` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7b** | 10 种药物 IC50 小提琴（High vs Low Risk） | 同上（`Figure_IC50_Raincloud_HAVCR2_docking.svg`） | 同上 | ★★ |
| **7c** | 分子对接结构图（2 行概览+放大） | `09_01_pymol_crizotinib/docking_overview.png`、`binding_site_closeup.png` | `09_structure_docking_MD/09_01_pymol_crizotinib.py` | ★★ |
| **7d** | 分子对接结构图（含残基标注） | `09_02_pymol_V75M/*` 等 | `09_structure_docking_MD/09_0*_pymol_*.py` | ★ |
| **7e** | Confidence metrics（pTM/ipTM/pLDDT）+ Binding pocket pLDDT 热图 | `08_01_ic50_and_boltz2/pocket_pLDDT_heatmap_2drugs.png` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7f** | RMSD 轨迹 vs Time (ns) | 需 GROMACS 输出 | 外部工具（MD） | — |
| **7g** | 结构对比（V73M vs WT） | `09_02_pymol_V75M/*` | `09_structure_docking_MD/09_02_pymol_V75M.py` | ★ |
| **7h** | 结构 / 配体视图 | `09_0*_pymol_*.py` 产出 | 同上 | ★ |
| **7i** | pLDDT / 结合口袋 RMSD + 分子内/间接触散点 | `08_01_ic50_and_boltz2/pLDDT_profile_2drugs.png`、`delta_pLDDT_profile.png` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7j** | 3 张 PAE 热图（WT / V73M / ΔPAE） | `08_01_ic50_and_boltz2/pocket_delta_heatmap.png` | 同上 | ★★ |

> ⚠ Fig7 的 c/d/g/h 为 PyMOL 渲染的**位图**，依赖外部工具，无法重出矢量；f 需 GROMACS。

---

## 汇总

| 主图 | panel 数 | 已定位 | 完整度 |
|---|---|---|---|
| Fig1 | 6 | 6 | **100%** |
| Fig2 | 6 | 4 | 67% |
| Fig3 | 7 | 7 | **100%** |
| Fig4 | 8 | 8 | **100%** |
| Fig5 | 6 | 6 | **100%** |
| Fig6 | 8 | 5 | 63% |
| Fig7 | 10 | 8 | 80%（含位图，部分依赖外部工具） |

## 已知偏差（需留意）

1. **`_paper.svg` 为原始文件**（与 `bulk-download/GSE182434/scenic/` 下字节完全一致），
   即渲染出来的就是论文所用版本，可放心使用。
2. **Fig3c 行数差异**：`Fig1_regulon_activity_heatmap_paper.svg` 含 26 个 regulon，
   原图仅展示自 `SPI1(+)` 起的 15 个（疑似原图做过裁剪）。数值一致。
3. **Fig3g 列数差异**：当前产出为 GO BP / KEGG / Reactome 三列，原图为四列
   （多一列 Metabolism）。内容同源，属脚本版本差异。
4. **Fig1c 标题与轴标签重叠**（`chromosome_heatmap.svg` 渲染可见），重排时需修正。
