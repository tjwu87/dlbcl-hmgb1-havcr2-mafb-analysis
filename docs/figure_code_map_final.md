# 主图 + 附图 子图构成 ↔ 代码位置（逐张核对版）

> 核对方法（三重）：① 直接阅读 `论文结果/主图/*.png` 逐 panel 记录内容；
> ② 从候选 SVG **提取文字**比对（基因名/数值/类型数）；③ 渲染候选图与原图目视比对。
>
> 状态标记：✅已验证（渲染或数值比对一致）｜🟡 高度吻合（文字/文件名吻合，未渲染比对）
> ｜⚠️ 需重写（当前代码产不出该 panel）｜❌ 缺数据/缺脚本
>
> 更新：2026-09-11。此前版本的错误（1d/1e 配错脚本）已在本版改正。

---

## Fig1 — 单细胞图谱 + 细胞通讯（6 panel）

| panel | 内容（原图实测） | 代码位置 | 状态 |
|---|---|---|---|
| 1a | UMAP「Refined Cell Types」+「Tissue (DLBCL vs Tonsil)」+ 14 类图例（CD8 T 5516…pDC/Other 41） | `10_figures/10_00_all_panels_GSE182434.py` → `umap_celltypes2`（原始 **18×7 in**） | ✅ 已跑出 |
| 1b | Marker 点图：14 类 × ~25 基因（CD3D/CD8A/CD4/FOXP3/MS4A1/CD38/MKI67/LYZ/GNLY/IGHG1/BCL6…），色=Fraction of cells，径=Mean expression | `01_scRNA_GSE182434/01_01_qc_integration.py` → `marker_dotplot` | ⏳ 脚本可跑（GEO 原始文件已就位并补好路径映射），正在运行 |
| 1c | inferCNV 拷贝数热图：行分 Normal B / Normal B (DLBCL) / Malignant B，列 chr1–22 | `10_figures/10_00_all_panels_GSE182434.py` → `chromosome_heatmap`（原始 **20×11 in**） | ✅ 已跑出 |
| 1d | 环形通讯网络：**9 种细胞类型**（Malignant B, Normal B, CD4 T, CD8 T, Tregs, TFH, NK, Mono/Mac, pDC/Other） | `02_cellcomm_LIANA/02_01_liana_malignant_vs_normal.py` → `draw_circle_comm_v3()` → `circle_plot_all` / `Malignant_vs_Normal_B_cell_circle_plot`（原始 13×13 in） | ⚠️ 脚本 bug（`fmt` 未定义）已修；重跑时因与 01_01 同时运行**内存不足**，需单独重跑 |
| 1e | 9×9 交互权重热图（值 0–18.89，行/列同 9 类） | 同上 → `heatmap_malignant_vs_normal`（原始 12×10 in） | ⚠️ 同上 |
| 1f | 差异互作散点：ΔMagnitude vs −log10(p)，Total 239 / Malignant-enriched 82 / Normal-enriched 29 | `02_cellcomm_LIANA/02_02_liana_visualization.py` → `volcano_malignant_vs_normal_monomac`（原始 12×8.5 in） | ✅ 数值完全吻合 |

> ❌ 之前错配记录（避免再犯）：1d 不是 `circle_plot_malignant_vs_normal_v2.svg`（那是 02_02 的版本），
> 1e 不是 `interaction_heatmap.svg`（那是 **14 类**，原图是 **9 类**）。
> 数据佐证：`cellcomm/liana_results_full.csv` = 14 类、`liana_results_13subtypes_full.csv` = 13 类，
> 9 类 = 14 类聚合 + B 细胞按恶性/正常拆分（正是 02_01 做的事）。

---

## Fig2 — 髓系亚群 + 拟时序（6 panel）

| panel | 内容（原图实测） | 代码位置 | 状态 |
|---|---|---|---|
| 2a | UMAP 髓系 5 亚型（Mono/DC_1/DC_2/IFN_TAM/LA_TAM） | ⚠️ **无脚本**。产物 `bulk-download/GSE182434/mono_mac/mac_subtype_umap2.png`（位图）存在 | ❌ 需按 `adata_mac_subtyped.h5ad` 重写 |
| 2b | 点图：5 亚型 × 25 基因（FCN1, S100A9, S100A8, S100A4, APOBEC3A, LTB, CLEC10A, CD1C, JAML, CD1E, APOC1, CCL18, APOE, CTSD, MT1H, MT1G, CCL8, CCL2, MT1X, DNASE1L3, CST3, CLEC9A, SNX3…） | ⚠️ 同上：`mac_subtype_dotplot2.png` 无脚本。基因集与 `01_02_myeloid_subtypes.py` 的 `subtype_markers` 一致 | ❌ 需重写 |
| 2c | 辐射式通讯网络（Malignant B 为枢纽，13 类） | `02_cellcomm_LIANA/02_02_liana_visualization.py` → `circle_plot_13subtypes` + `circle_MalignantB_sender/IFNTAM_sender/LATAM_sender/LATAM_receiver` | 🟡 |
| 2d | 左：PAGA 图（Mono n=127 / DC_1 n=84 / DC_2 n=173 / IFN_TAM n=65 / LA_TAM n=80，边权 0.27/0.84/1.00…）；右：DPT Pseudotime UMAP（Mono PT=0.20…LA_TAM PT=0.40） | `mono_mac_5subtypes.py`（仓库根目录，**尚未迁入仓库**）→ `monomac_paga_trajectory` / `monomac_trajectory_directed`；仓库副本 `03_trajectory_PAGA/03_01_paga_trajectory.py` 亦产出 | 🟡 |
| 2e | 拟时序热图：行基因 FCN1, S100A9, S100A4, MT1H, CCL8, CCL2, CXCL8, CCL19, CTSD, CLEC9A, CD1C, CLEC10A, LTB，按 Mono/IFN_TAM/LA_TAM/DC_2/DC_1 分组 | `03_trajectory_PAGA/03_01_paga_trajectory.py` → `monomac_gene_dynamics` | 🟡 |
| 2f | 4 张富集条形图：KEGG/GO-BP × (IFN_TAM up, LA_TAM up)，n=13/9/8/5 | `03_trajectory_PAGA/03_01_paga_trajectory.py` → `enrichment_barplot_IFN_TAM_vs_LA_TAM`（当前只出 2 张，需确认是否含 GO-BP） | ⚠️ 数量待核 |

---

## Fig3 — 调控子分析（7 panel）✅ 全部验证

| panel | 内容 | 代码位置 | 状态 |
|---|---|---|---|
| 3a | 受体基因拟时序热图（Mono/IFN_TAM/LA_TAM，▲正相关 ▼负相关） | `04_pySCENIC/04_01_pyscenic_regulons.py` → `receptor_pseudotime_heatmap_path` | ✅ 渲染比对一致 |
| 3b | 4 折线：CD74 r=0.418 / AXL 0.277 / HLA-DPA1 0.218 / HAVCR2 0.232 | 同上 → `receptor_pseudotime_lineplots_path` | ✅ |
| 3c | regulon 活性热图（MAFB(+), IRF1(+), STAT1(+)…26 行，原图展示自 SPI1(+) 起 15 行） | 同上 → `Fig1_regulon_activity_heatmap_paper` | ✅（原图疑似裁剪过） |
| 3d | RSS 散点：RSS-IFN_TAM vs RSS-LA_TAM（特异 22/7，共享 10） | 同上 → `Fig2_RSS_scatter_paper` | ✅ |
| 3e | 7 层 Sankey：Malignant B Cell → Secreted Ligand → Mono/Mac Receptor → Signaling Pathway → Core TF → Top5 Target → TAM Subtype | 同上 → `sankey_7layer_TF_targets_v3_paper` | ✅ |
| 3f | 4×2 折线：ETV5/MITF/IRF8/IRF7 vs MAFB_KO；ETV5/BHLHE40/MAFB/IRF8 vs AXL | 同上 → `lineplot_LA_TAM_TF_receptor_2x4_paper` | ✅ |
| 3g | 富集点图 GO-BP/KEGG/Reactome（MAFB regulon targets） | 同上 → `MAFB_gseapy_enrichment_dotplot_paper` | ✅（原图为 4 列含 Metabolism，现 3 列） |

---

## Fig4 — CellOracle 扰动（8 panel）✅ 全部验证

产出脚本统一为 `05_CellOracle/05_90_reproduce_figures_local.py`。

| panel | 内容 | 文件 | 状态 |
|---|---|---|---|
| 4a | 3 UMAP：Cell Subtypes / MAFB Expression / MAFB KO Velocity | `fig09_mafb_expression_umap` | ✅ |
| 4b | 3 小提琴：LA_TAM fate score，WT vs MAFB KO | `fig10_fate_score_violin_wt_vs_ko` | ✅ |
| 4c | 3 UMAP：WT (Baseline) / MAFB KO / Null (Randomized GRN) | `fig01_mafb_quiver_vector_field` | ✅ |
| 4d | MAFB 调控网络（Top15 靶基因 + 共调控因子） | `fig08_mafb_regulatory_network_FINAL_LOGIC_FIXED` | ✅ |
| 4e | 热图：All Significant Down-regulated Genes (KO vs WT) | `fig_horizontal_heatmap_all_deg` | ✅ |
| 4f | 3 小提琴：MAFB KO Effect on LA_TAM Program Score（Δ=−0.036/−0.167/−0.000） | `fig16_program_score_paired_violin` | ✅ Δ 值吻合 |
| 4g | ORA 点图（GO-BP/KEGG/MSigDB Hallmark，n=86 DEGs） | `fig_merged_ora_bubble_all_cells` | ✅ |
| 4h | 效应量散点（Fate/Program score + SNR） | `fig_effect_size_summary_dotplot` | ✅ |

---

## Fig5 — 空间转录组验证（6 panel）✅ 全部验证

产出脚本统一为 `06_spatial_GSE232853/06_01_spatial_analysis.py`。

| panel | 内容 | 文件 | 状态 |
|---|---|---|---|
| 5a | 2 UMAP：Cell Compartment (CD20 n=253/CD68 n=319)、Tissue (DLBCL n=250/Normal n=322) | `fig1_UMAP_dual` | ✅ |
| 5b | 火山图 CD68⁺ DLBCL vs Normal（Up 75 / Down 72） | `fig3_volcano_DEA_CD68` | ✅ |
| 5c | 3 小提琴 HMGB1/HAVCR2/MAFB（padj 2.48e-14 / 6.11e-74 / 1.45e-75） | `fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected` | ✅ |
| 5d | 双腔热图 DLBCL(n=114)/Normal(n=114) | `task5_dual_chamber_heatmap` | ✅ |
| 5e | 2 散点（B-cell HMGB1 vs Macrophage HAVCR2 r=−0.473；vs LA_TAM r=+0.852） | `fig4_scatter_HMGB1_HAVCR2_LATAM` | ✅ |
| 5f | GSEA ridge（18 条 Hallmark，NES −1.88…3.68） | `fig5_GSEA_joyplot` | ✅ |

---

## Fig6 — 预后模型（8 panel）

| panel | 内容（原图实测） | 代码位置 | 状态 |
|---|---|---|---|
| 6a | LASSO 系数路径（Number of Non-zero Genes 20→0；基因 BASP1/AAK1/CHMP4B/CHCHD10/CDK19/IFIT3/NAE1/MFGE8/PYROXD2/PRICKLE1/THUMPD3/TUBB1/UPPS1） | `07_bulk_prognosis/07_02_lasso_cv_curves.R` → `Fig1B_LASSO_coef_trajectory`；`07_01` → `A0b_LASSO_coef` | 🟡 |
| 6b | **5 张** KM：GSE10846 (Training, n=412)、GSE39571 (n=552)、GSE31312 (n=299)、GSE23501 (n=80)、TCGA-DLBC (n=48) | `07_01_prognosis_main.R` → `A5_MultiCohort_KM` | ⚠️ **队列对不上**：现产出为 GSE10846(181)/GSE23501(60)/GSE32918(49)/GSE31312(470)/TCGA(48)，原图是 GSE39571(552)、GSE31312(299) → 需核对验证集清单 |
| 6c | 森林图 Univariable/Multivariable（Risk Score, Age, Sex, Stage, IPI） | `07_01` → `A1_Cox_forest` | ✅ |
| 6d | 校准曲线 1/3/5 年 OS（MAE 0.042/0.061/0.01） | `07_01` → `A3_Calibration_curves` | 🟡 |
| 6e | 雷达图 22 种免疫细胞（Low vs High Risk，含显著性星号） | `07_03_immune_deconvolution.R` | ⚠️ 该脚本**无 savefig 调用**，且依赖 `D:/BadiduNetdiskDownload/R/GSE10846survive/GSE10846_series_matrix.txt.gz`（本地无）→ 需补齐/重写 |
| 6f | 4 散点（M2 Macrophage r=0.255 / CD8⁺ T r=0.008 / M1 r=−0.135 / B Cell r=0.236）+ Spearman 横向条形（GSE10846, n=412） | 同上 | ⚠️ 同上 |
| 6g | 21 个免疫检查点小提琴（CD200R1, CD274, CD276, CD40, CD48, CTLA4, ENTPD1, HAVCR1, ICOS, ICOSLG, LAIR1, LAG3, LGALS9, PDCD1, TIGIT, TMIGD2, TNFRSF14/18/9, TNFSF14/15/18/9…） | 同上 | ⚠️ 同上 |
| 6h | ICB Response Rate（66/34% vs 48/54%）+ TIDE/Dysfunction/Exclusion 小提琴 + TIDE Score 瀑布 | 同上 | ⚠️ 同上 |
| — | `A2_Nomogram`、`B1_RiskScore_tripanel`、`B2_KM_training`、`B3_TimeROC`、`B5_Calibration_*` 均不在 Fig6 中 → 属**补充图** | `07_01` | — |

---

## Fig7 — 药物敏感 + 分子对接（10 panel）

| panel | 内容 | 代码位置 | 状态 |
|---|---|---|---|
| 7a | Tim-3 Docking Score vs GDSC2 IC50 + Top-10 候选药条形 | `08_drug_GDSC/08_01_ic50_and_boltz2.py` → `top200_gdsc2_matched`（源文件在 `bulk-download/GDSC2/results/`，已复制入仓） | 🟡 |
| 7b | 10 种药物 IC50 小提琴（High vs Low Risk） | 同上 → `Figure_IC50_Raincloud_HAVCR2_docking` | 🟡 |
| 7c | 分子对接结构图（概览+放大） | `09_structure_docking_MD/09_01_pymol_crizotinib.py` → `docking_overview` / `binding_site_closeup`（**位图**） | ❌ PyMOL 位图，无法矢量重出 |
| 7d | 分子对接结构图（残基 GLN169/ASN260/LEU253/THR108/ASP103/ARG98/ASP104） | `09_0*_pymol_*.py` | ❌ 同上 |
| 7e | Confidence metrics 三角（pTM/ipTM/pLDDT）+ Binding pocket pLDDT 热图 | `08_01` → `pocket_pLDDT_heatmap_2drugs.png` | 🟡 位图 |
| 7f | RMSD 轨迹 vs Time (ns) | **需 GROMACS** | ❌ 无脚本无数据 |
| 7g | 结构对比 V73M vs WT | `09_02_pymol_V75M.py` | ❌ 位图 |
| 7h | 结构/配体视图（配体 A3P0） | `09_0*_pymol_*.py` | ❌ 位图 |
| 7i | pLDDT 曲线 + 分子内/间接触散点 | `08_01` → `pLDDT_profile_2drugs.png`、`delta_pLDDT_profile.png` | 🟡 位图 |
| 7j | 3 张 PAE 热图（WT/V73M/ΔPAE） | `08_01` → `pocket_delta_heatmap.png` | 🟡 位图 |

---

## 附图（`论文结果/附图/` 共 17 个文件，含新旧两套命名）

| 附图 | 内容 | 代码位置 | 状态 |
|---|---|---|---|
| FigS1 | QC + UMAP overview（Tissue/Patient/Leiden res=0.6/Cell Type） | `10_figures/10_00_all_panels_GSE182434.py` → `umap_overview2`、`celltype_composition_per_sample2`、`celltype_mean_proportion2`、`marker_genes2`；另 `01_scRNA_GSE182434/01_04_supplementary_panels.py` | 🟡 |
| FigS2 | 髓群相关补充 | `10_figures/FigS2/10_10_FigS2_supplementary.py`（产出文件名却是 `FigS8.svg`，**命名错位需修**） | ⚠️ |
| FigS3 | 受体优先级排序 | `10_figures/FigS8/10_12_tableS19.py` → `FigS3_receptor_prioritization` | 🟡 |
| FigS4 | 未定位 | — | ❌ |
| FigS5 | 空间转录组补充（ROI/AOI 设计、PCA/UMAP 批次检查、组织×腔室小提琴） | `06_spatial_GSE232853/06_01_spatial_analysis.py` → `S5A_ROI_AOI_design_schematic`、`S5B_PCA_UMAP_by_sample`、`S5C_2x2_stratified_violin_Q3corrected` | ✅ |
| FigS6 | 未定位 | — | ❌ |
| FigS7 | 未定位 | — | ❌ |
| FigS8 | LAM/Lipid-laden macrophage 特征打分（Jaitin 2019 / Kloosterman 2024，AUROC 0.85/0.79） | `10_figures/FigS8/10_11_FigS8_supplementary.py`；另有根目录 `FigS8.py`、`FigS2补图.py` | 🟡 |
| 其他 | `supp_A_pseudotime_distribution` / `supp_B_paga_connectivity_heatmap` / `supp_C_sample_composition` / `supp_D_cellcomm_strength` | `4SUPP.py`（仓库根目录，**尚未迁入仓库**） | 🟡 |

---

## 汇总

| 主图 | panel 数 | 已验证 | 可复现 | 需重写/缺 |
|---|---|---|---|---|
| Fig1 | 6 | 4（1a/1c/1f 已跑出，1b 运行中） | 6 | 1d/1e 需单独重跑（内存） |
| Fig2 | 6 | 3 | 4 | **2a/2b 需重写**（无脚本） |
| Fig3 | 7 | 7 | 7 | — |
| Fig4 | 8 | 8 | 8 | — |
| Fig5 | 6 | 6 | 6 | — |
| Fig6 | 8 | 4 | 4 | **6e–6h 需补 07_03**（无 savefig、缺 GEO 矩阵）；6b 队列需核对 |
| Fig7 | 10 | 2 | 6 | c/d/g/h 位图、f 需 GROMACS |
| 附图 | 8 组 | 1 | 5 | S4/S6/S7 未定位 |

## 阻塞项（需要你决定/提供）

1. **Fig2a/2b**：无脚本。我可以按 `adata_mac_subtyped.h5ad` + `01_02` 的 marker 基因表重写，工作量小。
2. **Fig6e–h**：`07_03_immune_deconvolution.R` 没有出图代码，且依赖 `GSE10846_series_matrix.txt.gz`（本地没有）。请提供该文件，或授权我按现有 `A4_Immune_infiltration` 的思路重写。
3. **Fig6b 队列对不上**：原图是 GSE39571(552)/GSE31312(299)，现脚本是 GSE10846(181)/GSE32918(49)。请确认哪套是论文用的。
4. **FigS4/S6/S7**：本地没有对应的代码或产物，需要你指认内容。
5. **Fig7 c/d/g/h**：PyMOL 位图。如果要矢量，需要在 PyMOL 里重新导出（脚本可给，但渲染要 PyMOL 环境）。


---

## 【决定性证据】`Tcell/生图/` —— 拼装原图的素材库（265 文件：225 png + 18 py + 10 pse）

作者把拼主图用的 panel 原件都存在这里，**文件名即出处**。已核对 Fig1：

| 生图文件 | 对应子图 | 尺寸(px) | 产出脚本 |
|---|---|---|---|
| `1A.png` | **Fig1a** 双 UMAP「Cell Type Annotation」 | 4157×2073 | `10_00_all_panels_GSE182434.py` → `umap_celltypes2` |
| `1B.png` | **Fig1b** marker 点图（14 类 × 基因，含星号） | 2792×839 | `01_01_qc_integration.py` → `marker_dotplot` |
| `1C.png` | **Fig1c** inferCNV 热图 | 2850×1670 | `10_00_all_panels_GSE182434.py` → `chromosome_heatmap` |
| `1E.png` | **Fig1d** 环形网络（9 类，Malignant vs Normal B） | 1643×1719 | `02_01_liana_malignant_vs_normal.py` → `circle_plot_malignant_vs_normal` |
| `1F.png` | **Fig1e** 9×9 权重热图（Top 300，值 2.99/6.85/18.89…） | 1381×1186 | `02_01_liana_malignant_vs_normal.py` → `heatmap_malignant_vs_normal` |
| `1H.png` | **Fig1f** 差异互作火山（Target: Monocytes/Macrophages） | 1778×1324 | `02_02_liana_visualization.py` → `volcano_malignant_vs_normal_monomac` |
| `1D.png` | ~~代谢通路热图（AUCell）~~ | 2928×1300 | **未用于 Fig1**（备选） |
| `1G.png` | ~~Top LR 圆形点图~~ | 2959×1993 | **未用于 Fig1**（备选） |
| `3A.png` | **Fig3a** 受体基因拟时序热图 | 1920×1724 | `04_01_pyscenic_regulons.py` → `receptor_pseudotime_heatmap_path` |

> **结论**：我上一版修正后的 1d/1e 映射（→ `02_01_liana_malignant_vs_normal.py`）
> 被 `生图/circle_plot_malignant_vs_normal.png` 与 `heatmap_malignant_vs_normal.png`
> 的文件名**直接证实**（这正是 02_01 的输出文件名）。

### 生图里其他可用线索（未逐一核对）

- Fig4：`GSE182434_celloracle_fig01…fig14`、`fig01…fig23` 系列 → `05_90_reproduce_figures_local.py`
- Fig6：`Fig1A_LASSO_CV_curve` / `Fig1B_LASSO_coef_trajectory` / `Fig2_Cox_forest_plot` /
  `Fig3A_Nomogram` / `Fig3B_Calibration_curves` / `Fig4_Immune_infiltration` /
  `Fig5_MultiCohort_KM_curves` / `B1_RiskScore_tripanel` / `B3_TimeROC` / `B5_Calibration_ggplot` /
  `RADAR.png`（6e）/ `checkpoints.png`（6g）/ `TIDE.png`（6h）/ `ciber.png`
- Fig7：`PAE_heatmaps_Criz_Dab` / `RMSD_Combined(_paper)` / `RMSF(_paper)` / `Rg_Complex_Comparison` /
  `SASA_*` / `Hbond_*` / `boltz*.png` / `havcr2_screening_setup_top200_gdsc2_matched` /
  子目录 `crizotinib/`（18 个 py：make_publication_figures*.py、render_docking.py、vina.py 等）、
  `03_V75M/`、`DABRAFENIB/`（PyMOL 脚本 + .pse 会话 + .cif 结构）
- 附图：`S1.png…S10.png`、`Figure S1-1/S1-2`、`Figure S5.1`、`Figure S6.1`、
  `supp_A/B/C/D_*`（→ `4SUPP.py`）、`S5A/S5B/S5C`（→ `06_01`）
- 髓系：`mac_subtype_umap2` / `mac_subtype_dotplot2`（2a/2b）、`monomac_paga_trajectory(_rawpt)` /
  `monomac_trajectory_directed(_rawpt)`（2d）、`monomac_gene_dynamics`（2e）、
  `enrichment_barplot_IFN_TAM_vs_LA_TAM`（2f）、`volcano_IFN_TAM_vs_LA_TAM` / `volcano_IFN_TAM_vs_rest`
- pySCENIC：`Fig1_regulon_activity_heatmap(_paper)` / `Fig2_RSS_scatter` / `sankey_7layer_TF_targets_v3` /
  `lineplot_LA_TAM_TF_receptor_2x4` / `MAFB_gseapy_enrichment_dotplot`
- 空间：`fig1_UMAP_dual` / `fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected` / `fig3_volcano_CD68_DLBCL_vs_Normal` /
  `fig4_scatter_HMGB1_HAVCR2_LATAM` / `fig5_GSEA_joyplot` / `task5_dual_chamber_heatmap` / `fig7_spatial_coevolution_split_heatmap`
