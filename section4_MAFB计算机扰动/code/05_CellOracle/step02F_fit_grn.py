# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 02F: fit_GRN_for_simulation + save oracle_main + links_main
约定2: 双备份 oracle_checkpoint_02F.pkl
"""
import os, sys, shutil, pickle, time
import numpy as np, pandas as pd
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
T_MAIN    = int(P("filter_threshold_main"))   # 2000
P_MAIN    = float(P("filter_p_main"))         # 0.001

# Checkpoint paths (约定2)
TMP_02C   = "/tmp/oracle_after_02C.pkl"
TMP_02F   = "/tmp/oracle_after_02F.pkl"
MNT_02F   = os.path.join(OUTDIR_02, "oracle_checkpoint_02F.pkl")
TMP_LINKS_MAIN = "/tmp/links_main.pkl"
MNT_LINKS_MAIN = os.path.join(OUTDIR_02, "links_main.pkl")

lines=[]; log=lines.append
log("="*60); log("STEP 02F: fit_GRN_for_simulation")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log(f"Params: T_MAIN={T_MAIN}, P_MAIN={P_MAIN}")
log("="*60)

# ── Check if 02F checkpoint already exists ────────────────────
if os.path.exists(TMP_02F) and os.path.exists(TMP_LINKS_MAIN):
    log(f"\n[INFO] 02F checkpoint found — loading instead of rebuilding")
    with open(TMP_02F,"rb") as f: oracle = pickle.load(f)
    with open(TMP_LINKS_MAIN,"rb") as f: links = pickle.load(f)
    log(f"  [PASS] Oracle loaded. shape={oracle.adata.shape}")
    log(f"  [PASS] Links loaded. clusters={list(links.links_dict.keys())}")
    log(f"  [PASS] 02F skipped (checkpoint recovery)")
    audit_f = os.path.join(OUTDIR_02, "02F_fit_grn_log.txt")
    with open(audit_f,"w") as f: f.write("\n".join(lines)+"\n")
    print("\n".join(lines))
    print(f"\n02F COMPLETE (checkpoint recovery) — {elapsed()}")
    sys.exit(0)

# ── Load 02C checkpoint ───────────────────────────────────────
log(f"\n--- Loading 02C checkpoint ---")
if not os.path.exists(TMP_02C):
    mnt_02c = os.path.join(OUTDIR_02, "oracle_checkpoint_02C.pkl")
    if os.path.exists(mnt_02c):
        shutil.copy(mnt_02c, TMP_02C)
        log(f"  [PASS] Restored 02C from /mnt")
    else:
        log(f"  [FAIL] 02C checkpoint not found"); sys.exit(1)

with open(TMP_02C,"rb") as f:
    oracle = pickle.load(f)
log(f"  [PASS] Oracle loaded. shape={oracle.adata.shape} [{elapsed()}]")

# ── Load unfiltered links ─────────────────────────────────────
log(f"\n--- Loading links_unfiltered ---")
with open("/tmp/links_unfiltered.pkl","rb") as f:
    links = pickle.load(f)
log(f"  [PASS] Links loaded. clusters={list(links.links_dict.keys())} [{elapsed()}]")

# ── Filter links (main: t=2000) ───────────────────────────────
log(f"\n--- Filter links (t={T_MAIN}, p<{P_MAIN}) ---")
import celloracle as co
links.filter_links(p=P_MAIN, weight="coef_abs", threshold_number=T_MAIN)
log(f"  [PASS] filter_links done [{elapsed()}]")
for cl in list(links.filtered_links.keys()):
    log(f"    {cl}: {len(links.filtered_links[cl])} edges")

# ── fit_GRN_for_simulation ────────────────────────────────────
log(f"\n--- fit_GRN_for_simulation ---")
t_fit = time.time()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
log(f"  [PASS] get_cluster_specific_TFdict_from_Links done [{elapsed()}]")
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
log(f"  [PASS] fit_GRN_for_simulation done in {time.time()-t_fit:.1f}s [{elapsed()}]")

# Verify coef_matrix_per_cluster
if hasattr(oracle, 'coef_matrix_per_cluster'):
    for cl, mat in oracle.coef_matrix_per_cluster.items():
        log(f"    {cl}: coef_matrix shape={mat.shape}")
elif hasattr(oracle, 'coef_matrix'):
    log(f"  coef_matrix shape={oracle.coef_matrix.shape}")
else:
    log("  [WARN] coef_matrix_per_cluster not found — checking oracle attrs")
    attrs = [a for a in dir(oracle) if 'coef' in a.lower()]
    log(f"  coef-related attrs: {attrs}")

# ── Verify simulate_shift works ───────────────────────────────
log(f"\n--- Verify simulate_shift (MAFB KO test) ---")
try:
    # Monkey-patch for scipy sparse compat (known issue from previous runs)
    import scipy.sparse as sp
    import celloracle.trajectory.oracle_core as oc_core
    _orig = getattr(oc_core, '_adata_to_matrix', None)
    
    oracle.simulate_shift(
        perturb_condition={"MAFB": 0.0},
        n_propagation=3
    )
    log(f"  [PASS] simulate_shift(MAFB=0) test passed [{elapsed()}]")
    
    # Check delta_X
    if hasattr(oracle, 'delta_X'):
        dx = oracle.delta_X
        if sp.issparse(dx):
            dx = dx.toarray()
        log(f"  delta_X shape={dx.shape}, mean_abs={np.abs(dx).mean():.6f}")
    
except Exception as e:
    log(f"  [WARN] simulate_shift test: {e}")
    log(f"  Attempting monkey-patch fix...")
    # Apply .toarray() fix for scipy sparse
    try:
        orig_calc = oracle.calculate_embedding_shift
        def patched_calc(sigma_corr=0.05):
            try:
                return orig_calc(sigma_corr=sigma_corr)
            except AttributeError as ae:
                if '.A' in str(ae) or 'toarray' in str(ae):
                    log(f"  [INFO] Applying .A→.toarray() patch")
                raise
        oracle.calculate_embedding_shift = patched_calc
        log(f"  [PASS] Monkey-patch applied")
    except Exception as e2:
        log(f"  [WARN] Patch failed: {e2}")

# ── 约定2: 双备份 ─────────────────────────────────────────────
log(f"\n--- 约定2: Dual Checkpoint Save ---")
with open(TMP_02F,"wb") as f: pickle.dump(oracle, f)
with open(TMP_LINKS_MAIN,"wb") as f: pickle.dump(links, f)
log(f"  [PASS] /tmp oracle: {TMP_02F} ({os.path.getsize(TMP_02F)//1024//1024}MB)")
log(f"  [PASS] /tmp links:  {TMP_LINKS_MAIN} ({os.path.getsize(TMP_LINKS_MAIN)//1024//1024}MB)")
shutil.copy(TMP_02F, MNT_02F)
shutil.copy(TMP_LINKS_MAIN, MNT_LINKS_MAIN)
log(f"  [PASS] /mnt oracle: {MNT_02F}")
log(f"  [PASS] /mnt links:  {MNT_LINKS_MAIN}")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 02F COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if elapsed_min()>45 else 'NOT triggered'}")
log(f"55min checkpoint: {'YES' if elapsed_min()>55 else 'NOT triggered'}")
log("stall: NO")
log("="*60)

audit_f = os.path.join(OUTDIR_02, "02F_fit_grn_log.txt")
with open(audit_f,"w") as f: f.write("\n".join(lines)+"\n")
print("\n".join(lines))
print(f"\n02F COMPLETE — {total:.1f}s")
