# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""STEP 00B: Freeze main analysis parameters. k_main=20 per user instruction."""
import os, sys, pandas as pd
from datetime import datetime

t0 = datetime.now()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUT_TSV = os.path.join(OUTDIR, "00_audit", "00B_frozen_parameters.tsv")

params = [
    ("cell_input_file", "/tmp/adata_mac_annotated_fixed.h5ad", "input", "yes", "no",
     "Fixed h5ad: uns/log1p/base null removed for anndata 0.10.9 compat",
     "Original: /mnt/results/GSE182434/trajectory/adata_mac_annotated.h5ad"),
    ("subset_subtypes", "Mono,IFN_TAM,LA_TAM", "preprocessing", "yes", "no",
     "TAM-oriented branch; DC_1/DC_2 excluded", "272 cells: Mono=127, IFN_TAM=65, LA_TAM=80"),
    ("HVG_layer_priority", "counts (fallback: raw_count)", "preprocessing", "yes", "no",
     "seurat_v3 requires counts-like layer; log1p_norm forbidden as HVG input (H1)", "counts confirmed PRESENT"),
    ("HVG_flavor", "seurat_v3", "preprocessing", "yes", "no",
     "Standard HVG for scRNA-seq; CellOracle best practice", ""),
    ("HVG_n_top_genes_main", "3000", "preprocessing", "yes", "no",
     "Primary HVG count; fixed per plan H2", ""),
    ("HVG_n_top_genes_sensitivity", "4000,5000", "preprocessing", "no", "yes",
     "Sensitivity only; must not change main analysis (H3)", "STEP 08C (4000), Fig21 (5000)"),
    ("expression_matrix_name", "log1p_norm", "oracle", "yes", "no",
     "CellOracle uses log1p_norm for PCA and imputation", "log1p_norm confirmed PRESENT"),
    ("n_pca_components", "50", "oracle", "yes", "no",
     "Standard PCA dims for CellOracle; sufficient for 272-cell dataset", ""),
    ("knn_imputation_k_candidate_set", "15,20,30,40", "oracle", "no", "yes",
     "All candidates per H15; main k selected on stability not result strength",
     "272 cells; k=15/20 preserve local structure; k=30/40 stronger smoothing"),
    ("knn_imputation_k_main", "20", "oracle", "yes", "no",
     "k=20 selected per user instruction: appropriate local smoothing for 272-cell dataset; "
     "preserves local trajectory structure; k=15/30/40 used as sensitivity (STEP 08E)",
     "Set to 20 per user instruction (hixing set) before STEP 02"),
    ("knn_imputation_k_sensitivity", "15,30,40", "oracle", "no", "yes",
     "Remaining candidates; used in STEP 08E k-sensitivity analysis",
     "Must not overwrite main analysis results"),
    ("knn_main_selection_rule",
     "k=20 selected per user instruction: local smoothing for 272-cell dataset; "
     "k=15/30/40 as sensitivity",
     "oracle", "yes", "no",
     "Rule: select k satisfying (1) Mono fate KO↓ stable, (2) LA_TAM highest SNR, "
     "(3) no direction flip at adjacent k, (4) no excessive fragmentation",
     "k=20 chosen per user instruction before STEP 02"),
    ("cluster_column", "mac_subtype", "oracle", "yes", "no",
     "obs column for GRN cluster fitting", ""),
    ("filter_p_main", "0.001", "grn_filter", "yes", "no",
     "Strict p-value threshold for main GRN edges", ""),
    ("filter_threshold_main", "2000", "grn_filter", "yes", "no",
     "Global top-N edge threshold for main analysis (strict_direct definition)",
     "threshold_number = global truncation, not per-TF (G6, H6)"),
    ("filter_p_supported", "0.001", "grn_filter", "no", "yes",
     "Same p threshold, relaxed threshold_number for supported_direct", ""),
    ("filter_threshold_supported", "5000", "grn_filter", "no", "yes",
     "Relaxed global top-N for supported_direct: present at t=5000, absent at t=2000",
     "supported_direct ≠ confirmed direct binding (H7)"),
    ("transition_n_neighbors_main", "200", "perturbation", "yes", "no",
     "n_neighbors for estimate_transition_prob; covers most of 272-cell dataset", ""),
    ("transition_n_neighbors_sensitivity", "30", "perturbation", "no", "yes",
     "Sensitivity only (STEP 08D)", ""),
    ("sigma_corr", "0.05", "perturbation", "yes", "no",
     "Bandwidth for calculate_embedding_shift; CellOracle default", ""),
    ("null_model", "use_randomized_GRN=True", "perturbation", "yes", "no",
     "Randomized GRN as noise baseline for SNR", ""),
    ("fate_definition", "sum of transition_prob toward LA_TAM neighbors",
     "fate", "yes", "no",
     "Unique formal fate score definition (H9)", ""),
    ("WT_baseline_rule", "per-subtype imputed mean", "fate", "yes", "no",
     "WT baseline must be per-subtype, not global mean (H11)", ""),
    ("OE_cascade_rule",
     "candidate_1=2*max; candidate_2=p95*1.5; candidate_3=max+1*std",
     "perturbation", "yes", "no",
     "Cascade test for OE value; first passing candidate used (H13)",
     "Each candidate tested on fresh Oracle copy; failed copy discarded"),
    ("random_seed", "42", "global", "yes", "no",
     "Fixed seed for all random processes where API supports it", ""),
    ("random_seed_applicability",
     "knn_imputation: NO (no seed param); "
     "estimate_transition_prob: YES (random_seed param); "
     "simulate_shift: PARTIAL (use_randomized_GRN flag only)",
     "global", "yes", "no",
     "Verified in STEP 00A; limitation noted in final report",
     "knn_imputation randomness cannot be fixed; minor batch variation possible"),
]

header = ["parameter_name","value","scope","is_primary_analysis",
          "allowed_to_vary_in_sensitivity","rationale","notes"]
with open(OUT_TSV, "w") as f:
    f.write("\t".join(header) + "\n")
    for row in params:
        f.write("\t".join(str(x) for x in row) + "\n")

print(f"[PASS] Frozen parameters saved: {OUT_TSV}")
print(f"       {len(params)} parameters")

# Verify key values
df = pd.read_csv(OUT_TSV, sep="\t", index_col="parameter_name")
for k in ["knn_imputation_k_main","knn_imputation_k_sensitivity","filter_threshold_main"]:
    print(f"  {k}: {df.loc[k,'value']}")

print(f"\n[PASS] 00B complete. Elapsed: {(datetime.now()-t0).total_seconds():.1f}s")
print("  knn_imputation_k_main = 20")
