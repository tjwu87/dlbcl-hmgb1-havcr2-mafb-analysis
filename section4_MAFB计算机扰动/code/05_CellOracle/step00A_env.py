# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 00A: Environment & Input Validation
Reads OUTDIR from tmp_manifest/outdir.txt
All checks must PASS before proceeding.
"""
import sys, os, importlib, inspect, glob
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR_FILE = os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")
with open(OUTDIR_FILE) as f:
    OUTDIR = f.read().strip()

AUDIT_FILE = os.path.join(OUTDIR, "00_audit", "00_environment_audit.txt")
lines = []; log = lines.append; FAIL_FLAGS = []

log("=" * 60)
log("STEP 00A: Environment & Input Audit")
log(f"Start: {t0.isoformat()}")
log(f"OUTDIR: {OUTDIR}")
log(f"Python: {sys.version.split()[0]}")
log("=" * 60)

# ── 1. Package versions ───────────────────────────────────────
log("\n--- 1. Package Versions ---")
for label, pkg in [("celloracle","celloracle"),("scanpy","scanpy"),("anndata","anndata"),
                    ("numpy","numpy"),("pandas","pandas"),("scipy","scipy"),
                    ("sklearn","sklearn"),("matplotlib","matplotlib"),
                    ("seaborn","seaborn"),("gseapy","gseapy")]:
    try:
        mod = importlib.import_module(pkg)
        ver = getattr(mod, "__version__", "unknown")
        log(f"  [PASS] {label}: {ver}")
    except ImportError as e:
        log(f"  [FAIL] {label}: {e}")
        FAIL_FLAGS.append(f"{label} import failed")

# ── 2. celloracle functional ──────────────────────────────────
log("\n--- 2. celloracle functional ---")
try:
    import celloracle as co
    log(f"  [PASS] import celloracle as co: OK ({co.__version__})")
    try:
        o = co.Oracle(); del o
        log("  [PASS] co.Oracle() instantiation: OK")
    except Exception as e:
        log(f"  [FAIL] co.Oracle(): {e}"); FAIL_FLAGS.append("Oracle() failed")
except Exception as e:
    log(f"  [FAIL] import celloracle: {e}"); FAIL_FLAGS.append("celloracle import failed"); co=None

# ── 3. Base GRN ───────────────────────────────────────────────
log("\n--- 3. Base GRN ---")
import pandas as pd
GRN_PATH = None
for p in [translate("/mnt/results/GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv"),
          translate("/mnt/results/GSE182434/celloracle_paper_final_20260331_034717/02_oracle_grn/base_GRN_human_promoter.csv")]:
    if os.path.exists(p):
        grn = pd.read_csv(p, index_col=0).reset_index()
        log(f"  [PASS] GRN loaded from local cache: {p}")
        log(f"         shape={grn.shape}, cols[:4]={grn.columns.tolist()[:4]}")
        GRN_PATH = p; break
if not GRN_PATH:
    log("  [FAIL] Base GRN not found in any candidate path")
    FAIL_FLAGS.append("Base GRN not found")

# ── 4. h5ad ───────────────────────────────────────────────────
log("\n--- 4. h5ad Input ---")
import anndata as ad, numpy as np
H5AD_FIXED = "/tmp/adata_mac_annotated_fixed.h5ad"
adata = None
if os.path.exists(H5AD_FIXED):
    try:
        adata = ad.read_h5ad(H5AD_FIXED)
        log(f"  [PASS] Fixed h5ad readable: shape={adata.shape}")
        log(f"         h5ad in use: {H5AD_FIXED}")
        log(f"         Fix applied: uns/log1p/base null removed")
    except Exception as e:
        log(f"  [FAIL] Fixed h5ad read error: {e}"); FAIL_FLAGS.append("h5ad read failed")
else:
    log(f"  [FAIL] {H5AD_FIXED} not found — run setup_env.sh first")
    FAIL_FLAGS.append("h5ad not found")

# ── 5. obs['mac_subtype'] ─────────────────────────────────────
log("\n--- 5. obs['mac_subtype'] ---")
if adata is not None:
    if "mac_subtype" in adata.obs.columns:
        vc = adata.obs["mac_subtype"].value_counts()
        log("  [PASS] mac_subtype present")
        for k,v in vc.items(): log(f"         {k}: {v}")
    else:
        log("  [FAIL] mac_subtype missing"); FAIL_FLAGS.append("mac_subtype missing")

# ── 6. Layers ─────────────────────────────────────────────────
log("\n--- 6. Layers ---")
if adata is not None:
    log(f"  Available layers: {list(adata.layers.keys())}")
    if "counts" in adata.layers:
        log("  [PASS] counts layer: PRESENT")
    elif "raw_count" in adata.layers:
        log("  [WARN] counts absent; raw_count PRESENT (will use raw_count)")
    else:
        log("  [FAIL] Neither counts nor raw_count"); FAIL_FLAGS.append("No counts layer")
    if "log1p_norm" in adata.layers:
        log("  [PASS] log1p_norm layer: PRESENT")
    else:
        log("  [FAIL] log1p_norm missing"); FAIL_FLAGS.append("log1p_norm missing")

# ── 7. UMAP ───────────────────────────────────────────────────
log("\n--- 7. UMAP Embedding ---")
if adata is not None:
    log(f"  obsm keys: {list(adata.obsm.keys())}")
    if "X_umap_paga" in adata.obsm:
        log("  [PASS] X_umap_paga: PRESENT → will use X_umap_paga")
    elif "X_umap" in adata.obsm:
        log("  [PASS] X_umap: PRESENT → will use X_umap")
    else:
        log("  [FAIL] No UMAP found"); FAIL_FLAGS.append("No UMAP")

# ── 8. dpt_pseudotime ────────────────────────────────────────
log("\n--- 8. dpt_pseudotime ---")
if adata is not None:
    if "dpt_pseudotime" in adata.obs.columns:
        dpt = adata.obs["dpt_pseudotime"]
        log(f"  [PASS] dpt_pseudotime present, range=[{dpt.min():.4f},{dpt.max():.4f}], nan={dpt.isna().sum()}")
    else:
        log("  [WARN] dpt_pseudotime absent — will compute in STEP 01D")

# ── 9. 23 pySCENIC genes ─────────────────────────────────────
log("\n--- 9. 23 pySCENIC genes in var_names ---")
SCENIC = ["CEP192","CHCHD10","NAGK","PPRC1","BASP1","MTR","TIMM17A","AAK1",
          "MGAT1","USP12","FAM20A","APOL3","GSDMD","ELMO1","PICALM","RASSF4",
          "IFIT3","NAE1","NAP1L4","PPP2R5A","MRGBP","PPP1R18","CTSL"]
n_present = 0
if adata is not None:
    for g in SCENIC:
        s = "PRESENT" if g in adata.var_names else "ABSENT"
        if g in adata.var_names: n_present += 1
        log(f"    {g}: {s}")
    log(f"  Total: {n_present}/23")

# ── 10. Existing celloracle_* dirs ───────────────────────────
log("\n--- 10. Existing celloracle_* directories ---")
for pat in [translate("/mnt/results/GSE182434/celloracle_*/"),
            translate("/mnt/results/GSE182434/celloracle_paper_final_*/")]:
    for d in sorted(glob.glob(pat)):
        log(f"  EXISTS: {os.path.basename(d.rstrip('/'))}")

# ── 11. Seed support ─────────────────────────────────────────
log("\n--- 11. Seed parameter support ---")
if co is not None:
    for mname in ["knn_imputation","estimate_transition_prob","simulate_shift"]:
        try:
            sig = inspect.signature(getattr(co.Oracle, mname))
            params = list(sig.parameters.keys())
            seed_p = [p for p in params if "seed" in p.lower() or "random" in p.lower()]
            log(f"  {mname}: seed/random params = {seed_p if seed_p else 'NONE — API does not support seed'}")
        except Exception as e:
            log(f"  {mname}: ERROR {e}")

# ── 12. Current OUTDIR ────────────────────────────────────────
log(f"\n--- 12. Current OUTDIR ---")
log(f"  {OUTDIR}")
log(f"  Timestamp: {t0.strftime('%Y%m%d_%H%M%S')}")

# ── Final verdict ─────────────────────────────────────────────
log(f"\n{'='*60}")
log(f"Elapsed: {elapsed()}")
if FAIL_FLAGS:
    log("RESULT: FAIL")
    for f in FAIL_FLAGS: log(f"  FAIL: {f}")
    log("STOPPING")
else:
    log("RESULT: ALL PASS")
log("="*60)

with open(AUDIT_FILE, "w") as f: f.write("\n".join(lines)+"\n")
print(f"Saved: {AUDIT_FILE}")
print("\n".join(lines))
if FAIL_FLAGS: sys.exit(1)
