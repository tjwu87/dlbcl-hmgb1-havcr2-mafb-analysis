# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 03 Null 条件补跑
Null 定义：oracle 副本上调用 simulate_shift(perturb_condition={}, use_randomized_GRN=True)
其余参数与 KO 完全相同：n_neighbors=200, knn_random=True, sampled_fraction=1, sigma_corr=0.05
来源：从 02_oracle_grn/oracle_checkpoint_02F.pkl 加载独立副本
保存：03_perturbation/Null/oracle_Null.pkl, delta_X_Null.csv, delta_embedding_Null.csv
"""
import os, sys, shutil, pickle, copy, time
import numpy as np, pandas as pd
import scipy.sparse as sp
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

OUTDIR_NULL = os.path.join(OUTDIR, "03_perturbation", "Null")
os.makedirs(OUTDIR_NULL, exist_ok=True)

SEED    = int(P("random_seed"))                  # 42
SIGMA   = float(P("sigma_corr"))                 # 0.05
N_NEIGH = int(P("transition_n_neighbors_main"))  # 200
N_PROP  = 3  # CellOracle default

def log(msg):
    print(msg, flush=True)

log("="*60)
log("STEP 03 Null: simulate_shift(perturb_condition={}, use_randomized_GRN=True)")
log(f"Start: {t0.isoformat()}")
log(f"OUTDIR: {OUTDIR}")
log(f"Params: seed={SEED}, sigma={SIGMA}, n_neighbors={N_NEIGH}, n_prop={N_PROP}")
log("="*60)

# ── Check if already done ─────────────────────────────────────
out_pkl = os.path.join(OUTDIR_NULL, "oracle_Null.pkl")
out_dx  = os.path.join(OUTDIR_NULL, "delta_X_Null.csv")
out_de  = os.path.join(OUTDIR_NULL, "delta_embedding_Null.csv")

# NOTE: always rerun — perturb_condition corrected to {"MAFB": 0}
# (previous run used {} which gave delta_X=0; overwrite old files)

# ── Load oracle_main (02F checkpoint) ─────────────────────────
log(f"\n--- Loading oracle_main (02F checkpoint) ---")
TMP_02F = "/tmp/oracle_after_02F.pkl"
if not os.path.exists(TMP_02F):
    mnt_02f = os.path.join(OUTDIR, "02_oracle_grn", "oracle_checkpoint_02F.pkl")
    if os.path.exists(mnt_02f):
        shutil.copy(mnt_02f, TMP_02F)
        log(f"  [PASS] Restored 02F from /mnt")
    else:
        log(f"  [FAIL] 02F checkpoint not found"); sys.exit(1)

with open(TMP_02F,"rb") as f:
    oracle_main = pickle.load(f)
log(f"  [PASS] Oracle loaded. shape={oracle_main.adata.shape} [{elapsed()}]")

# ── Null: deepcopy + simulate_shift(use_randomized_GRN=True) ──
log(f"\n--- Null: deepcopy + simulate_shift ---")
oracle_null = copy.deepcopy(oracle_main)
log(f"  deepcopy done [{elapsed()}]")

oracle_null.simulate_shift(
    perturb_condition={"MAFB": 0},
    use_randomized_GRN=True,
    n_propagation=N_PROP
)
log(f"  [PASS] simulate_shift(perturb_condition={{MAFB:0}}, use_randomized_GRN=True) done [{elapsed()}]")

# delta_X from layers
delta_X = oracle_null.adata.layers['delta_X']
if sp.issparse(delta_X):
    delta_X = delta_X.toarray()
log(f"  delta_X shape={delta_X.shape}, mean_abs={np.abs(delta_X).mean():.6f}")

# ── estimate_transition_prob ──────────────────────────────────
log(f"\n--- estimate_transition_prob ---")
oracle_null.estimate_transition_prob(
    n_neighbors=N_NEIGH,
    knn_random=True,
    sampled_fraction=1,
    random_seed=SEED
)
log(f"  [PASS] estimate_transition_prob(n_neighbors={N_NEIGH}) done [{elapsed()}]")

# ── calculate_embedding_shift ─────────────────────────────────
log(f"\n--- calculate_embedding_shift ---")
oracle_null.calculate_embedding_shift(sigma_corr=SIGMA)
log(f"  [PASS] calculate_embedding_shift(sigma={SIGMA}) done [{elapsed()}]")

delta_emb = oracle_null.delta_embedding
if sp.issparse(delta_emb):
    delta_emb = delta_emb.toarray()
log(f"  delta_embedding shape={delta_emb.shape}, mean_abs={np.abs(delta_emb).mean():.6f}")

# ── Save results ──────────────────────────────────────────────
log(f"\n--- Saving results ---")
cell_ids = oracle_null.adata.obs_names.tolist()
gene_ids = oracle_null.adata.var_names.tolist()

# oracle pkl
with open(out_pkl,"wb") as f:
    pickle.dump(oracle_null, f)
log(f"  [PASS] oracle_Null.pkl saved ({os.path.getsize(out_pkl)//1024//1024}MB)")

# delta_X
dx_df = pd.DataFrame(delta_X, index=cell_ids, columns=gene_ids)
dx_df.to_csv(out_dx)
log(f"  [PASS] delta_X_Null.csv saved: shape={dx_df.shape}")

# delta_embedding
de_df = pd.DataFrame(delta_emb, index=cell_ids, columns=['UMAP1','UMAP2'])
de_df.to_csv(out_de)
log(f"  [PASS] delta_embedding_Null.csv saved: shape={de_df.shape}")

# ── Final report ──────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 03 Null COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if elapsed_min()>45 else 'NOT triggered'}")
log(f"55min checkpoint: {'YES' if elapsed_min()>55 else 'NOT triggered'}")
log(f"stall: NO")
log("="*60)
log(f"\n>>> KEY RESULT: delta_embedding mean_abs = {np.abs(delta_emb).mean():.6f}")
log(f">>> (compare: MAFB_KO=0.454095, MAFB_OE=0.284353, WT≈0)")
log(f"\n停止等待确认。")
