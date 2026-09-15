# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>


import os, pickle
import numpy as np
import pandas as pd

OUTDIR = translate("/mnt/results/GSE182434/celloracle_paper_final_20260331_085223")
DOUT   = f"{OUTDIR}/09_figures/data"
os.makedirs(DOUT, exist_ok=True)

print("Loading oracle objects...")
def load_pkl(path):
    with open(path, "rb") as f:
        return pickle.load(f)

oracle_wt   = load_pkl(f"{OUTDIR}/03_perturbation/WT/oracle_WT.pkl")
oracle_ko   = load_pkl(f"{OUTDIR}/03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl")
oracle_oe   = load_pkl(f"{OUTDIR}/03_perturbation/MAFB_OE/oracle_MAFB_OE.pkl")
oracle_null = load_pkl(f"{OUTDIR}/03_perturbation/Null/oracle_Null.pkl")
print("All oracle objects loaded")

adata = oracle_wt.adata
gene_names  = np.array(adata.var_names)
subtypes    = adata.obs["mac_subtype"].values
umap_coords = adata.obsm["X_umap"]
cell_ids    = np.array(adata.obs_names)

# ── 1. Cell metadata (UMAP + subtype) ─────────────────────────────────────────
df_meta = pd.DataFrame({
    "cell_id":  cell_ids,
    "subtype":  subtypes,
    "UMAP1":    umap_coords[:, 0],
    "UMAP2":    umap_coords[:, 1],
})
df_meta.to_csv(f"{DOUT}/cell_metadata.csv", index=False)
print(f"Saved cell_metadata.csv ({len(df_meta)} rows)")

# ── 2. log1p_norm expression (for Fig05, Fig09, Fig11) ────────────────────────
X_log = adata.layers["log1p_norm"]
if hasattr(X_log, "toarray"): X_log = X_log.toarray()
df_log = pd.DataFrame(X_log, index=cell_ids, columns=gene_names)
df_log.to_csv(f"{DOUT}/log1p_norm_expression.csv")
print(f"Saved log1p_norm_expression.csv ({df_log.shape})")

# ── 3. imputed_count WT (for Fig23 DEG) ───────────────────────────────────────
X_imp = oracle_wt.adata.layers["imputed_count"]
if hasattr(X_imp, "toarray"): X_imp = X_imp.toarray()
df_imp = pd.DataFrame(X_imp, index=cell_ids, columns=gene_names)
df_imp.to_csv(f"{DOUT}/imputed_count_WT.csv")
print(f"Saved imputed_count_WT.csv ({df_imp.shape})")

# ── 4. simulated_count KO (for Fig23 DEG) ─────────────────────────────────────
X_sim = oracle_ko.adata.layers["simulated_count"]
if hasattr(X_sim, "toarray"): X_sim = X_sim.toarray()
df_sim = pd.DataFrame(X_sim, index=cell_ids, columns=gene_names)
df_sim.to_csv(f"{DOUT}/simulated_count_KO.csv")
print(f"Saved simulated_count_KO.csv ({df_sim.shape})")

# ── 5. cell_fate_scores_MAFB_KO.csv (missing file) ────────────────────────────
# Reconstruct from oracle_ko obs if fate score is stored there
obs_ko = oracle_ko.adata.obs
fate_cols = [c for c in obs_ko.columns if "fate" in c.lower() or "score" in c.lower()]
print(f"KO obs fate columns: {fate_cols}")
if fate_cols:
    df_fate_ko = obs_ko[fate_cols].copy()
    df_fate_ko.index.name = "cell_id"
    df_fate_ko.to_csv(f"{OUTDIR}/04_fate/cell_fate_scores_MAFB_KO.csv")
    print(f"Saved cell_fate_scores_MAFB_KO.csv ({df_fate_ko.shape})")
else:
    # Try to get from WT oracle (fate scores computed on WT adata)
    obs_wt = oracle_wt.adata.obs
    fate_cols_wt = [c for c in obs_wt.columns if "fate" in c.lower() or "score" in c.lower()]
    print(f"WT obs fate columns: {fate_cols_wt}")
    # Load from existing cell_fate_scores_WT.csv as reference
    df_fate_wt = pd.read_csv(f"{OUTDIR}/04_fate/cell_fate_scores_WT.csv", index_col=0)
    print(f"WT fate scores columns: {df_fate_wt.columns.tolist()}")

print("Data export complete")
