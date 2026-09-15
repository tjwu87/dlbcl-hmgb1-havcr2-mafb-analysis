# Authoritative main-figure ↔ panel ↔ code mapping

> Basis for these judgements (triple cross-validation, not guesswork):
> 1. **Reading the original main figures directly** (the assistant can read images), recording each panel's content;
> 2. **Exporting panel text** (`docs/panel_text_dump.md`, extracting readable text from the SVGs);
> 3. **Render-and-compare** — rendering candidate SVGs and comparing them panel by panel against the original (`tools/panel_contact_sheet.py`).
>
> Confidence: ★★★ = visually confirmed after rendering; ★★ = filename/text matches closely; ★ = inferred from the analysis pipeline.

Generated: 2026-09-11

---

## Fig1 — single-cell atlas + cell communication (6 panels)

| Panel | Content | Producing file (repo-relative) | Producing script | Confidence |
|---|---|---|---|---|
| **1a** | UMAP "Refined Cell Types" + "Tissue (DLBCL vs Tonsil)" + type legend | `10_00_all_panels_GSE182434/umap_celltypes2.svg` | `10_figures/10_00_all_panels_GSE182434.py` | ★★★ |
| **1b** | Marker dot plot: 14 cell types × genes (colour = fraction of cells, size = mean expression) | `01_01_qc_integration/marker_dotplot.svg` | `01_scRNA_GSE182434/01_01_qc_integration.py` | ★★★ |
| **1c** | inferCNV copy-number heatmap (rows grouped Normal B / Normal B (DLBCL) / Malignant B) | `10_00_all_panels_GSE182434/chromosome_heatmap.svg` | `10_figures/10_00_all_panels_GSE182434.py` | ★★★ |
| **1d** | Circular cell–cell interaction network (Malignant vs Normal B) | `02_02_liana_visualization/circle_plot_malignant_vs_normal_v2.svg` | `02_cellcomm_LIANA/02_02_liana_visualization.py` | ★★★ |
| **1e** | Interaction-strength heatmap (source × target, 14 cell types) | `02_02_liana_visualization/interaction_heatmap.svg` | same as above | ★★★ |
| **1f** | Differential interaction scatter (ΔMagnitude vs −log10 p; Malignant 82 / Normal 29) | `02_02_liana_visualization/volcano_malignant_vs_normal_monomac.svg` | same as above | ★★★ |

---

## Fig2 — myeloid subtypes + pseudotime (6 panels)

| Panel | Content | Producing file | Producing script | Confidence |
|---|---|---|---|---|
| **2a** | UMAP of myeloid subtypes (Mono / DC_1 / LA_TAM / IFN_TAM / DC_2) | `01_02_myeloid_subtypes/*` (sub-figure names pending render check) | `01_scRNA_GSE182434/01_02_myeloid_subtypes.py` | ★★ |
| **2b** | Dot plot: 5 subtypes × genes (FCN1/S100A8/C1QA/CD74/CCL18/MARCO/TREM2/SPP1…) | same as above | same as above | ★★ |
| **2c** | Communication network: one all-subtype circle + four radial plots centred on Malignant B / LA_TAM / IFN_TAM / Normal B | `02_02_liana_visualization/circle_plot_13subtypes.svg` + `circle_*_sender.svg` (4) | `02_cellcomm_LIANA/02_02_liana_visualization.py` | ★★★ (circle) / ★★ (radial) |
| **2d** | 2 UMAPs (subtype / TF differentiation) with PT annotation | `03_01_paga_trajectory/monomac_paga_trajectory.svg`, `monomac_trajectory_directed.svg` | `03_trajectory_PAGA/03_01_paga_trajectory.py` | ★★ |
| **2e** | Pseudotime heatmap (Mono/IFN_TAM/LA_TAM/DC_2/DC_1 × genes) | `03_01_paga_trajectory/monomac_gene_dynamics.svg` or `scenic_tf_heatmap.svg` | same as above | ★ |
| **2f** | 4 enrichment bar plots: KEGG/GO-BP × IFN_TAM/LA_TAM up | `03_01_paga_trajectory/enrichment_barplot_IFN_TAM_vs_LA_TAM.svg` | same as above | ★★ |

---

## Fig3 — regulon analysis (7 panels)

| Panel | Content | Producing file | Producing script | Confidence |
|---|---|---|---|---|
| **3a** | Heatmap: receptor genes (ordered by peak expression) × cells (Mono/IFN_TAM/LA_TAM), with ▲ positive ▼ negative correlation | `04_01_pyscenic_regulons/receptor_pseudotime_heatmap_path.svg` | `04_pySCENIC/04_01_pyscenic_regulons.py` (also in `03_01`) | ★★★ |
| **3b** | 4 line plots: CD74 (r=0.418) / AXL (0.277) / HLA-DPA1 (0.218) / HAVCR2 (0.232) vs pseudotime | `04_01_pyscenic_regulons/receptor_pseudotime_lineplots_path.svg` | same as above | ★★★ |
| **3c** | Heatmap: regulon activity (MAFB(+)·IRF1(+)·STAT1(+)…) × Mono/IFN_TAM/LA_TAM | `04_01_pyscenic_regulons/Fig1_regulon_activity_heatmap_paper.svg` | same as above | ★★★ |
| **3d** | RSS scatter: RSS-IFN_TAM vs RSS-LA_TAM (22 LA_TAM-specific / 7 IFN_TAM-specific / 10 shared) | `04_01_pyscenic_regulons/Fig2_RSS_scatter_paper.svg` | same as above | ★★★ |
| **3e** | Sankey (7 layers): Malignant B Cell → Secreted Ligand → Mono/Mac Receptor → Signaling Pathway → Core TF → Top5 Target → TAM Subtype | `04_01_pyscenic_regulons/sankey_7layer_TF_targets_v3_paper.svg` | same as above | ★★★ |
| **3f** | 4×2 line plots: ETV5/MITF/IRF8/IRF7 vs MAFB_KO; ETV5/BHLHE40/MAFB/IRF8 vs AXL | `04_01_pyscenic_regulons/lineplot_LA_TAM_TF_receptor_2x4_paper.svg` | same as above | ★★★ |
| **3g** | Enrichment dot plot: GO BP / KEGG / Reactome (MAFB regulon targets) | `04_01_pyscenic_regulons/MAFB_gseapy_enrichment_dotplot_paper.svg` | same as above | ★★★ |

**Not used**: `Fig3_TF_target_network_paper.svg`, `pseudotime_TF_receptor_dynamics_paper.svg` (not adopted in the original figure; usable as supplementary).

---

## Fig4 — CellOracle perturbation (8 panels) ★ all confirmed

| Panel | Content | Producing file | Confidence |
|---|---|---|---|
| **4a** | 3 UMAPs: Cell Subtypes / MAFB Expression / MAFB KO – Velocity Magnitude | `05_90_reproduce_figures_local/fig09_mafb_expression_umap.svg` | ★★★ |
| **4b** | 3 violins: LA_TAM fate score in Mono / IFN_TAM / LA_TAM, WT vs MAFB KO | `05_90_reproduce_figures_local/fig10_fate_score_violin_wt_vs_ko.svg` | ★★★ |
| **4c** | 3 UMAPs: WT (Baseline) / MAFB KO / Null (Randomized GRN) | `05_90_reproduce_figures_local/fig01_mafb_quiver_vector_field.svg` | ★★★ |
| **4d** | MAFB regulatory network (Top15 target genes + co-regulators) | `05_90_reproduce_figures_local/fig08_mafb_regulatory_network_FINAL_LOGIC_FIXED.svg` | ★★★ |
| **4e** | Heatmap: All Significant Down-regulated Genes (MAFB KO vs WT) | `05_90_reproduce_figures_local/fig_horizontal_heatmap_all_deg.svg` | ★★★ |
| **4f** | 3 violins: MAFB KO Effect on LA_TAM Program Score (Δ = −0.036/−0.167/−0.000) | `05_90_reproduce_figures_local/fig16_program_score_paired_violin.svg` | ★★★ |
| **4g** | Dot plot: Down-regulated genes (GO BP / KEGG / MSigDB Hallmark) | `05_90_reproduce_figures_local/fig_merged_ora_bubble_all_cells.svg` | ★★★ |
| **4h** | Scatter: MAFB KO Effect-Size (Fate/Program score + SNR) | `05_90_reproduce_figures_local/fig_effect_size_summary_dotplot.svg` | ★★★ |

> Producing script: `05_CellOracle/05_90_reproduce_figures_local.py`

---

## Fig5 — spatial transcriptomics validation (6 panels) ★ all confirmed

| Panel | Content | Producing file | Confidence |
|---|---|---|---|
| **5a** | 2 UMAPs: Cell Compartment (CD20/CD68), Tissue Type (DLBCL/Normal) | `06_01_spatial_analysis/fig1_UMAP_dual.svg` | ★★★ |
| **5b** | Volcano: CD68⁺ Macrophages DLBCL vs Normal (MAFB / HAVCR2 / HMGB1 annotated) | `06_01_spatial_analysis/fig3_volcano_DEA_CD68.svg` | ★★★ |
| **5c** | 3 violins: HMGB1 / HAVCR2 / MAFB (CD20 vs CD68) | `06_01_spatial_analysis/fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected.svg` | ★★★ |
| **5d** | Dual-chamber heatmap: CD20⁺ vs CD68⁺ paired ROIs (DLBCL n=114 / Normal n=114) | `06_01_spatial_analysis/task5_dual_chamber_heatmap.svg` | ★★★ |
| **5e** | 2 scatters: B-cell HMGB1 vs Macrophage HAVCR2 / vs LA_TAM Score | `06_01_spatial_analysis/fig4_scatter_HMGB1_HAVCR2_LATAM.svg` | ★★★ |
| **5f** | GSEA ridge plot (Depleted ↔ Enriched in LA_TAM_Rich, 18 Hallmark sets) | `06_01_spatial_analysis/fig5_GSEA_joyplot.svg` | ★★★ |

> Producing script: `06_spatial_GSE232853/06_01_spatial_analysis.py`

---

## Fig6 — prognostic model (8 panels)

| Panel | Content | Producing file | Producing script | Confidence |
|---|---|---|---|---|
| **6a** | LASSO coefficient paths + number of non-zero coefficients | `07_02_lasso_cv_curves/Fig1B_LASSO_coef_trajectory.pdf`, `Fig1A_LASSO_CV_curve.pdf` | `07_bulk_prognosis/07_02_lasso_cv_curves.R` | ★★★ |
| **6b** | 6 KM curves (5-cohort validation) | `07_01_prognosis_main/A5_MultiCohort_KM.pdf` | `07_bulk_prognosis/07_01_prognosis_main.R` | ★★★ |
| **6c** | Forest plot: univariable / multivariable Cox | `07_01_prognosis_main/A1_Cox_forest.pdf` | same as above | ★★★ |
| **6d** | Calibration curves: 1/3/5-year OS | `07_01_prognosis_main/A3_Calibration_curves.pdf` | same as above | ★★★ |
| **6e** | Radar plot: 22 immune cell infiltrates | `07_01_prognosis_main/A4_Immune_infiltration.pdf` | same as above | ★★ |
| **6f** | 4 scatters (H2/M1 Macrophage, CD8 T, B Cell vs Risk Score) + Spearman rho bar chart | To confirm (candidates: contained in `A4_Immune_infiltration.pdf`; or produced by `07_03_immune_deconvolution.R`) | `07_bulk_prognosis/07_03_immune_deconvolution.R` | ★ |
| **6g** | Immune-checkpoint violin plot (20 genes) | To confirm | same as above | ★ |
| **6h** | ICB Response Rate + TIDE / T-cell Exclusion | To confirm | same as above | ★ |

> ⚠ `A2_Nomogram.pdf` does not appear in Fig6 → it should be a **supplementary figure**.

---

## Fig7 — drug sensitivity + molecular docking (10 panels)

> Not included in this release: the paper removed this Results subsection. The mapping
> is retained here for traceability of the original figure set.

| Panel | Content | Producing file | Producing script | Confidence |
|---|---|---|---|---|
| **7a** | Tim-3 Docking Score vs GDSC2 IC50 scatter + Top-10 candidate drug bars | `08_01_ic50_and_boltz2/top200_gdsc2_matched.svg`, `Figure_IC50_Raincloud_HAVCR2_docking.svg` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7b** | IC50 violins for 10 drugs (High vs Low Risk) | same as above | same as above | ★★ |
| **7c** | Docking structure figure (2-row overview + close-up) | `09_01_pymol_crizotinib/docking_overview.png`, `binding_site_closeup.png` | `09_structure_docking_MD/09_01_pymol_crizotinib.py` | ★★ |
| **7d** | Docking structure figure (with residue labels) | `09_02_pymol_V75M/*` etc. | `09_structure_docking_MD/09_0*_pymol_*.py` | ★ |
| **7e** | Confidence metrics (pTM/ipTM/pLDDT) + binding-pocket pLDDT heatmap | `08_01_ic50_and_boltz2/pocket_pLDDT_heatmap_2drugs.png` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7f** | RMSD trajectory vs time (ns) | requires GROMACS output | external tool (MD) | — |
| **7g** | Structural comparison (V73M vs WT) | `09_02_pymol_V75M/*` | `09_structure_docking_MD/09_02_pymol_V75M.py` | ★ |
| **7h** | Structure / ligand view | produced by `09_0*_pymol_*.py` | same as above | ★ |
| **7i** | pLDDT / binding-pocket RMSD + intra-/inter-molecular contact scatter | `08_01_ic50_and_boltz2/pLDDT_profile_2drugs.png`, `delta_pLDDT_profile.png` | `08_drug_GDSC/08_01_ic50_and_boltz2.py` | ★★ |
| **7j** | 3 PAE heatmaps (WT / V73M / ΔPAE) | `08_01_ic50_and_boltz2/pocket_delta_heatmap.png` | same as above | ★★ |

> ⚠ Fig7's c/d/g/h are **bitmaps** rendered by PyMOL and depend on an external tool, so they
> cannot be re-exported as vectors; f requires GROMACS.

---

## Summary

| Main figure | Panels | Located | Completeness |
|---|---|---|---|
| Fig1 | 6 | 6 | **100%** |
| Fig2 | 6 | 4 | 67% |
| Fig3 | 7 | 7 | **100%** |
| Fig4 | 8 | 8 | **100%** |
| Fig5 | 6 | 6 | **100%** |
| Fig6 | 8 | 5 | 63% |
| Fig7 | 10 | 8 | 80% (includes bitmaps; some depend on external tools) |

## Known discrepancies (worth noting)

1. **The `_paper.svg` files are the original files** (byte-identical to those under
   `bulk-download/GSE182434/scenic/`), i.e. what they render is the version used in the
   paper — safe to use.
2. **Fig3c row-count difference**: `Fig1_regulon_activity_heatmap_paper.svg` contains 26
   regulons, while the original figure shows only the 15 starting from `SPI1(+)`
   (the original was apparently cropped). The values agree.
3. **Fig3g column-count difference**: the current output has three columns
   (GO BP / KEGG / Reactome), while the original has four (plus a Metabolism column).
   Same underlying content; a script-version difference.
4. **Fig1c title overlaps the axis labels** (visible when rendering
   `chromosome_heatmap.svg`); to be fixed during re-layout.
