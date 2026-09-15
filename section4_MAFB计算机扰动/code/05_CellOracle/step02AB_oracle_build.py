# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 02A: Prepare Oracle input (verify adata_hvg_counts3000.h5ad)
STEP 02B: Init Oracle + PCA + knn_imputation(k=20)
约定2: 双备份 oracle_checkpoint_02B.pkl → /tmp/ + $OUTDIR/02_oracle_grn/
"""
import os, sys, shutil, pickle, time
import numpy as np, pandas as pd
import anndata as ad
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"
def elapsed_min(): return (datetime.now()-t0).total_seconds()/60

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

params_df = pd.read_csv(os.path.join(OUTDIR,"00_audit","00B_frozen_parameters.tsv"),
                        sep="\t", index_col="parameter_name")
def P(k): return params_df.loc[k,"value"]

OUTDIR_02 = os.path.join(OUTDIR, "02_oracle_grn")
GRN_PATH  = translate("/mnt/results/GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv")
H5AD_PATH = os.path.join(OUTDIR, "01_preprocessing", "adata_hvg_counts3000.h5ad")
K_MAIN    = int(P("knn_imputation_k_main"))   # 20
N_PCA     = int(P("n_pca_components"))         # 50
EXPR_MAT  = P("expression_matrix_name")        # log1p_norm
CLUSTER_COL = P("cluster_column")              # mac_subtype

# Checkpoint paths (约定2)
TMP_CKPT  = "/tmp/oracle_after_02B.pkl"
MNTCKPT   = os.path.join(OUTDIR_02, "oracle_checkpoint_02B.pkl")

lines=[]; log=lines.append
log("="*60); log("STEP 02A+02B: Oracle Build")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log(f"Params: k={K_MAIN}, n_pca={N_PCA}, expr_mat={EXPR_MAT}")
log("="*60)

# ── Check if checkpoint already exists (约定2 recovery) ───────
if os.path.exists(TMP_CKPT):
    log(f"\n[INFO] Checkpoint found at {TMP_CKPT} — loading instead of rebuilding")
    with open(TMP_CKPT,"rb") as f:
        oracle = pickle.load(f)
    log(f"  [PASS] Oracle loaded from checkpoint. shape={oracle.adata.shape}")
    log(f"  [PASS] 02A+02B skipped (checkpoint recovery)")
    # Write audit
    audit_path = os.path.join(OUTDIR_02, "02A_input_audit.txt")
    with open(audit_path,"w") as f: f.write("\n".join(lines)+"\n")
    audit_b = os.path.join(OUTDIR_02, "02B_oracle_build_log.txt")
    with open(audit_b,"w") as f: f.write("\n".join(lines)+"\n")
    print("\n".join(lines))
    print(f"\n02A+02B COMPLETE (checkpoint recovery) — {elapsed()}")
    sys.exit(0)

# ── 02A: Verify input ─────────────────────────────────────────
log("\n--- 02A: Prepare Oracle Input ---")
assert os.path.exists(H5AD_PATH), f"FAIL: {H5AD_PATH} not found"
adata = ad.read_h5ad(H5AD_PATH)
log(f"  [PASS] adata_hvg_counts3000.h5ad: shape={adata.shape}")
log(f"         layers: {list(adata.layers.keys())}")

# Ensure raw_count
if "raw_count" not in adata.layers:
    if "counts" in adata.layers:
        adata.layers["raw_count"] = adata.layers["counts"].copy()
        log("  [PASS] counts → raw_count copied")
    else:
        log("  [FAIL] No counts/raw_count"); sys.exit(1)
else:
    log("  [PASS] raw_count already present")

# Verify mac_subtype
log(f"  mac_subtype counts: {dict(adata.obs['mac_subtype'].value_counts())}")

# ── CRITICAL: set X = raw_count (integer counts) before Oracle import ──
# import_anndata_as_raw_count requires X to be non-negative raw counts
import scipy.sparse as sp
raw = adata.layers["raw_count"]
if sp.issparse(raw):
    adata.X = raw.copy()
else:
    import scipy.sparse as sp2
    adata.X = sp2.csr_matrix(raw.copy())
X_check = adata.X.toarray() if sp.issparse(adata.X) else adata.X
assert X_check.min() >= 0, f"FAIL: X still has negatives after raw_count swap: min={X_check.min()}"
log(f"  [PASS] X set to raw_count (min={X_check.min():.0f}, max={X_check.max():.0f})")

# Save 02A audit
audit_a = os.path.join(OUTDIR_02, "02A_input_audit.txt")
with open(audit_a,"w") as f: f.write("\n".join(lines)+"\n")
log(f"  [PASS] 02A done [{elapsed()}]")

# ── 02B: Init Oracle + PCA + imputation ──────────────────────
log("\n--- 02B: Oracle Init + PCA + Imputation ---")
log(f"  Parameters from frozen TSV:")
log(f"    expression_matrix_name = {EXPR_MAT}")
log(f"    n_pca_components = {N_PCA}")
log(f"    knn_imputation_k = {K_MAIN}")
log(f"    cluster_column = {CLUSTER_COL}")

import celloracle as co
log(f"  celloracle version: {co.__version__}")

# Init Oracle
oracle = co.Oracle()
log(f"  [PASS] co.Oracle() instantiated [{elapsed()}]")

# Import adata
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name=CLUSTER_COL,
    embedding_name="X_umap"
)
log(f"  [PASS] import_anndata_as_raw_count done [{elapsed()}]")

# Import base GRN
grn = pd.read_csv(GRN_PATH, index_col=0).reset_index()
log(f"  [PASS] Base GRN loaded: shape={grn.shape}")
oracle.import_TF_data(TF_info_matrix=grn)
log(f"  [PASS] import_TF_data done [{elapsed()}]")

# PCA
oracle.perform_PCA()
log(f"  [PASS] perform_PCA(n_components={N_PCA}) done [{elapsed()}]")

# knn_imputation
log(f"  Running knn_imputation(k={K_MAIN})...")
t_imp = time.time()
oracle.knn_imputation(k=K_MAIN)
log(f"  [PASS] knn_imputation(k={K_MAIN}) done in {time.time()-t_imp:.1f}s [{elapsed()}]")

# Verify imputed_count
assert "imputed_count" in oracle.adata.layers, "FAIL: imputed_count layer missing"
log(f"  [PASS] imputed_count layer present")
log(f"  Oracle adata shape: {oracle.adata.shape}")

# ── 约定2: 双备份 ─────────────────────────────────────────────
log(f"\n--- 约定2: Dual Checkpoint Save ---")
# /tmp backup
with open(TMP_CKPT,"wb") as f: pickle.dump(oracle, f)
log(f"  [PASS] /tmp checkpoint: {TMP_CKPT} ({os.path.getsize(TMP_CKPT)//1024//1024}MB)")
# $OUTDIR persistent backup
shutil.copy(TMP_CKPT, MNTCKPT)
log(f"  [PASS] /mnt checkpoint: {MNTCKPT} ({os.path.getsize(MNTCKPT)//1024//1024}MB)")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 02A+02B COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if elapsed_min()>45 else 'NOT triggered'}")
log(f"55min checkpoint: {'YES' if elapsed_min()>55 else 'NOT triggered'}")
log("stall: NO")
log("="*60)

audit_b = os.path.join(OUTDIR_02, "02B_oracle_build_log.txt")
with open(audit_b,"w") as f: f.write("\n".join(lines)+"\n")
with open(audit_a,"w") as f: f.write("\n".join(lines)+"\n")
print("\n".join(lines))
print(f"\n02A+02B COMPLETE — {total:.1f}s")
