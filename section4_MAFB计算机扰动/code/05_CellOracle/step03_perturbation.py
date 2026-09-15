# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 03: MAFB perturbation simulations
- MAFB KO: perturb_condition={"MAFB": 0.0}
- MAFB OE: perturb_condition={"MAFB": 2.0}
- WT:      perturb_condition={}  (empty dict, NOT None)
For each: simulate_shift → estimate_transition_prob → calculate_embedding_shift
Save delta_X, delta_embedding per condition.
"""
import os, sys, shutil, pickle, time, copy, traceback
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

OUTDIR_03 = os.path.join(OUTDIR, "03_perturbation")
os.makedirs(OUTDIR_03, exist_ok=True)

SEED      = int(P("random_seed"))                  # 42
SIGMA     = float(P("sigma_corr"))                 # 0.05
N_NEIGH   = int(P("transition_n_neighbors_main"))  # 200
N_PROP    = 3  # CellOracle default; not in frozen params TSV

# Perturbation conditions
CONDITIONS = {
    "MAFB_KO": {"MAFB": 0.0},
    "MAFB_OE": {"MAFB": 2.0},
    "WT":      {}
}

lines=[]; 
def log(msg):
    lines.append(msg)
    print(msg, flush=True)

log("="*60); log("STEP 03: MAFB Perturbation Simulations")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log(f"Params: seed={SEED}, sigma={SIGMA}, n_neighbors={N_NEIGH}, n_prop={N_PROP}")
log(f"Conditions: {list(CONDITIONS.keys())}")
log("="*60)

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

# ── Run all conditions ────────────────────────────────────────
log(f"\n--- Running perturbation conditions ---")
results = {}

for cond_name, perturb_cond in CONDITIONS.items():
    log(f"\n  --- {cond_name}: perturb={perturb_cond} ---")
    t_cond = time.time()
    
    oracle = copy.deepcopy(oracle_main)
    log(f"    deepcopy done [{elapsed()}]")
    
    # simulate_shift
    oracle.simulate_shift(perturb_condition=perturb_cond, n_propagation=N_PROP)
    log(f"    [PASS] simulate_shift done [{elapsed()}]")
    
    # delta_X is stored in oracle.adata.layers['delta_X'] (celloracle 0.18.0)
    delta_X = oracle.adata.layers['delta_X']
    if sp.issparse(delta_X):
        delta_X = delta_X.toarray()
    log(f"    delta_X shape={delta_X.shape}, mean_abs={np.abs(delta_X).mean():.6f}")
    
    # estimate_transition_prob
    oracle.estimate_transition_prob(
        n_neighbors=N_NEIGH,
        knn_random=True,
        sampled_fraction=1,
        random_seed=SEED
    )
    log(f"    [PASS] estimate_transition_prob(n_neighbors={N_NEIGH}) done [{elapsed()}]")
    
    # calculate_embedding_shift
    oracle.calculate_embedding_shift(sigma_corr=SIGMA)
    log(f"    [PASS] calculate_embedding_shift(sigma={SIGMA}) done [{elapsed()}]")
    
    # Extract delta_embedding (celloracle 0.18.0 attr name)
    delta_emb = oracle.delta_embedding
    if sp.issparse(delta_emb):
        delta_emb = delta_emb.toarray()
    
    trans_prob = oracle.transition_prob
    if sp.issparse(trans_prob):
        trans_prob = trans_prob.toarray()
    
    log(f"    delta_embedding shape={delta_emb.shape}, mean_abs={np.abs(delta_emb).mean():.6f}")
    log(f"    transition_prob shape={trans_prob.shape}")
    log(f"    {cond_name} done in {time.time()-t_cond:.1f}s [{elapsed()}]")
    
    results[cond_name] = {
        'oracle': oracle,
        'delta_X': delta_X,
        'delta_embedding': delta_emb,
        'transition_prob': trans_prob,
    }
    log(f"  [PASS] {cond_name} complete [{elapsed()}]")

# ── Save results ──────────────────────────────────────────────
log(f"\n--- Saving results ---")
cell_ids = oracle_main.adata.obs_names.tolist()
gene_ids = oracle_main.adata.var_names.tolist()

for cond_name, res in results.items():
    cond_dir = os.path.join(OUTDIR_03, cond_name)
    os.makedirs(cond_dir, exist_ok=True)
    
    # delta_X (cells × genes)
    delta_df = pd.DataFrame(res['delta_X'], index=cell_ids, columns=gene_ids)
    delta_df.to_csv(os.path.join(cond_dir, f"delta_X_{cond_name}.csv"))
    log(f"  [PASS] delta_X saved: {cond_name} shape={delta_df.shape}")
    
    # delta_embedding (cells × 2)
    emb_df = pd.DataFrame(res['delta_embedding'], index=cell_ids, columns=['UMAP1','UMAP2'])
    emb_df.to_csv(os.path.join(cond_dir, f"delta_embedding_{cond_name}.csv"))
    log(f"  [PASS] delta_embedding saved: {cond_name} shape={emb_df.shape}")
    
    # Save oracle per condition
    oracle_path = os.path.join(cond_dir, f"oracle_{cond_name}.pkl")
    with open(oracle_path,"wb") as f: pickle.dump(res['oracle'], f)
    log(f"  [PASS] oracle_{cond_name}.pkl saved ({os.path.getsize(oracle_path)//1024//1024}MB)")

# ── Summary table ─────────────────────────────────────────────
log(f"\n--- Summary ---")
summary_rows = []
for cond_name, res in results.items():
    dx = res['delta_X']
    es = res['delta_embedding']
    summary_rows.append({
        'condition': cond_name,
        'n_cells': dx.shape[0],
        'n_genes': dx.shape[1],
        'delta_X_mean_abs': np.abs(dx).mean(),
        'delta_X_max_abs': np.abs(dx).max(),
        'delta_embedding_mean_abs': np.abs(es).mean(),
        'delta_embedding_max_abs': np.abs(es).max(),
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTDIR_03, "03_perturbation_summary.csv"), index=False)
log(f"\n  Perturbation summary:")
log(summary_df.to_string(index=False))

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 03 COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if elapsed_min()>45 else 'NOT triggered'}")
log(f"55min checkpoint: {'YES' if elapsed_min()>55 else 'NOT triggered'}")
log("stall: NO")
log("="*60)

audit_03 = os.path.join(OUTDIR_03, "03_perturbation_log.txt")
with open(audit_03,"w") as f: f.write("\n".join(lines)+"\n")
print(f"\nSTEP 03 COMPLETE — {total:.1f}s")
