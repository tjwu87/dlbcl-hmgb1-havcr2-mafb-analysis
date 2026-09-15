# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 02C: get_links (unfiltered) for all 3 clusters
约定2: 双备份 oracle_checkpoint_02C.pkl + links_unfiltered.pkl
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
SEED      = int(P("random_seed"))          # 42
SIGMA     = float(P("sigma_corr"))         # 0.05

# Checkpoint paths (约定2)
TMP_02B   = "/tmp/oracle_after_02B.pkl"
TMP_02C   = "/tmp/oracle_after_02C.pkl"
MNT_02C   = os.path.join(OUTDIR_02, "oracle_checkpoint_02C.pkl")
TMP_LINKS = "/tmp/links_unfiltered.pkl"
MNT_LINKS = os.path.join(OUTDIR_02, "links_unfiltered.pkl")

lines=[]; log=lines.append
log("="*60); log("STEP 02C: get_links (unfiltered)")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log(f"Params: sigma_corr={SIGMA}, seed={SEED}")
log("="*60)

# ── Check if 02C checkpoint already exists ────────────────────
if os.path.exists(TMP_02C) and os.path.exists(TMP_LINKS):
    log(f"\n[INFO] 02C checkpoint found — loading instead of rebuilding")
    with open(TMP_02C,"rb") as f: oracle = pickle.load(f)
    with open(TMP_LINKS,"rb") as f: links = pickle.load(f)
    log(f"  [PASS] Oracle loaded. shape={oracle.adata.shape}")
    log(f"  [PASS] Links loaded. clusters={list(links.links_dict.keys())}")
    log(f"  [PASS] 02C skipped (checkpoint recovery)")
    audit_c = os.path.join(OUTDIR_02, "02C_get_links_log.txt")
    with open(audit_c,"w") as f: f.write("\n".join(lines)+"\n")
    print("\n".join(lines))
    print(f"\n02C COMPLETE (checkpoint recovery) — {elapsed()}")
    sys.exit(0)

# ── Load 02B checkpoint ───────────────────────────────────────
log(f"\n--- Loading 02B checkpoint ---")
if not os.path.exists(TMP_02B):
    # Restore from /mnt
    mnt_02b = os.path.join(OUTDIR_02, "oracle_checkpoint_02B.pkl")
    if os.path.exists(mnt_02b):
        shutil.copy(mnt_02b, TMP_02B)
        log(f"  [PASS] Restored 02B from /mnt")
    else:
        log(f"  [FAIL] 02B checkpoint not found"); sys.exit(1)

with open(TMP_02B,"rb") as f:
    oracle = pickle.load(f)
log(f"  [PASS] Oracle loaded. shape={oracle.adata.shape} [{elapsed()}]")

# ── 02C: get_links ────────────────────────────────────────────
log(f"\n--- 02C: get_links ---")
import celloracle as co
log(f"  celloracle version: {co.__version__}")

t_links = time.time()
links = oracle.get_links(
    cluster_name_for_GRN_unit=oracle.cluster_column_name,
    alpha=10,
    verbose_level=10,
    test_mode=False
)
log(f"  [PASS] get_links done in {time.time()-t_links:.1f}s [{elapsed()}]")

# Verify
clusters = list(links.links_dict.keys())
log(f"  Clusters: {clusters}")
for cl in clusters:
    df = links.links_dict[cl]
    log(f"    {cl}: {len(df)} edges, cols={list(df.columns)}")

# Save coef CSVs (unfiltered)
coef_dir = os.path.join(OUTDIR_02, "coef_unfiltered")
os.makedirs(coef_dir, exist_ok=True)
for cl in clusters:
    df = links.links_dict[cl]
    out = os.path.join(coef_dir, f"coef_unfiltered_{cl}.csv")
    df.to_csv(out, index=False)
    log(f"  [PASS] Saved {out}")

# ── 约定2: 双备份 ─────────────────────────────────────────────
log(f"\n--- 约定2: Dual Checkpoint Save ---")
# Save oracle + links to /tmp
with open(TMP_02C,"wb") as f: pickle.dump(oracle, f)
with open(TMP_LINKS,"wb") as f: pickle.dump(links, f)
log(f"  [PASS] /tmp oracle: {TMP_02C} ({os.path.getsize(TMP_02C)//1024//1024}MB)")
log(f"  [PASS] /tmp links:  {TMP_LINKS} ({os.path.getsize(TMP_LINKS)//1024//1024}MB)")
# Copy to /mnt
shutil.copy(TMP_02C, MNT_02C)
shutil.copy(TMP_LINKS, MNT_LINKS)
log(f"  [PASS] /mnt oracle: {MNT_02C}")
log(f"  [PASS] /mnt links:  {MNT_LINKS}")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 02C COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if elapsed_min()>45 else 'NOT triggered'}")
log(f"55min checkpoint: {'YES' if elapsed_min()>55 else 'NOT triggered'}")
log("stall: NO")
log("="*60)

audit_c = os.path.join(OUTDIR_02, "02C_get_links_log.txt")
with open(audit_c,"w") as f: f.write("\n".join(lines)+"\n")
print("\n".join(lines))
print(f"\n02C COMPLETE — {total:.1f}s")
