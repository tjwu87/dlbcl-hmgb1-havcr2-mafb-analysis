# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 04: Perturbation scoring + figures
- Fate bias score (LA_TAM fate) per condition
- UMAP quiver plots (delta_embedding arrows)
- Gene-level response: top up/down genes per condition vs WT
- Per-cluster delta_X summary
"""
import os, sys, shutil, pickle, copy
import numpy as np, pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

params_df = pd.read_csv(os.path.join(OUTDIR,"00_audit","00B_frozen_parameters.tsv"),
                        sep="\t", index_col="parameter_name")
def P(k): return params_df.loc[k,"value"]

OUTDIR_03 = os.path.join(OUTDIR, "03_perturbation")
OUTDIR_04 = os.path.join(OUTDIR, "04_figures")
os.makedirs(OUTDIR_04, exist_ok=True)

CLUSTER_COL = P("cluster_column")  # mac_subtype
CONDITIONS = ["MAFB_KO", "MAFB_OE", "WT"]
CLUSTER_COLORS = {"Mono": "#4878CF", "IFN_TAM": "#D65F5F", "LA_TAM": "#6ACC65"}

lines=[]; 
def log(msg):
    lines.append(msg)
    print(msg, flush=True)

log("="*60); log("STEP 04: Perturbation Scoring + Figures")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log("="*60)

# ── Load oracle_main for cell metadata ────────────────────────
log(f"\n--- Loading oracle_main ---")
TMP_02F = "/tmp/oracle_after_02F.pkl"
if not os.path.exists(TMP_02F):
    shutil.copy(os.path.join(OUTDIR,"02_oracle_grn","oracle_checkpoint_02F.pkl"), TMP_02F)
with open(TMP_02F,"rb") as f:
    oracle_main = pickle.load(f)
log(f"  [PASS] Oracle loaded. shape={oracle_main.adata.shape} [{elapsed()}]")

# Cell metadata
obs = oracle_main.adata.obs.copy()
umap = oracle_main.adata.obsm['X_umap']
cell_ids = oracle_main.adata.obs_names.tolist()
gene_ids = oracle_main.adata.var_names.tolist()
clusters = obs[CLUSTER_COL].values
log(f"  Clusters: {dict(pd.Series(clusters).value_counts())}")

# ── Load perturbation results ─────────────────────────────────
log(f"\n--- Loading perturbation results ---")
delta_X_dict = {}
delta_emb_dict = {}
trans_prob_dict = {}

for cond in CONDITIONS:
    cond_dir = os.path.join(OUTDIR_03, cond)
    dx = pd.read_csv(os.path.join(cond_dir, f"delta_X_{cond}.csv"), index_col=0)
    de = pd.read_csv(os.path.join(cond_dir, f"delta_embedding_{cond}.csv"), index_col=0)
    delta_X_dict[cond] = dx.values
    delta_emb_dict[cond] = de.values
    log(f"  [PASS] {cond}: delta_X={dx.shape}, delta_emb={de.shape}")

# ── 04A: Fate bias score (LA_TAM fate) ───────────────────────
log(f"\n--- 04A: Fate Bias Score (LA_TAM fate) ---")
# Load transition_prob from oracle pkl
fate_scores = {}
for cond in CONDITIONS:
    cond_dir = os.path.join(OUTDIR_03, cond)
    oracle_pkl = os.path.join(cond_dir, f"oracle_{cond}.pkl")
    with open(oracle_pkl,"rb") as f:
        oracle_c = pickle.load(f)
    
    trans_prob = oracle_c.transition_prob
    if sp.issparse(trans_prob):
        trans_prob = trans_prob.toarray()
    
    # Fate score = sum of transition_prob toward LA_TAM cells
    la_tam_mask = (clusters == "LA_TAM")
    fate_score = trans_prob[:, la_tam_mask].sum(axis=1)
    fate_scores[cond] = fate_score
    log(f"  {cond}: LA_TAM fate score mean={fate_score.mean():.4f}, std={fate_score.std():.4f}")

# Save fate scores
fate_df = pd.DataFrame(fate_scores, index=cell_ids)
fate_df['cluster'] = clusters
fate_df['dpt_pseudotime'] = obs['dpt_pseudotime'].values
fate_df.to_csv(os.path.join(OUTDIR_04, "04A_fate_scores_LA_TAM.csv"))
log(f"  [PASS] Fate scores saved [{elapsed()}]")

# Per-cluster fate score summary
log(f"\n  Per-cluster fate score summary:")
for cl in ["Mono", "IFN_TAM", "LA_TAM"]:
    mask = (clusters == cl)
    row = {"cluster": cl}
    for cond in CONDITIONS:
        row[f"{cond}_mean"] = fate_scores[cond][mask].mean()
        row[f"{cond}_std"] = fate_scores[cond][mask].std()
    log(f"    {cl}: KO={row['MAFB_KO_mean']:.4f}±{row['MAFB_KO_std']:.4f}, "
        f"OE={row['MAFB_OE_mean']:.4f}±{row['MAFB_OE_std']:.4f}, "
        f"WT={row['WT_mean']:.4f}±{row['WT_std']:.4f}")

# ── 04B: UMAP quiver plots ────────────────────────────────────
log(f"\n--- 04B: UMAP Quiver Plots ---")

def plot_quiver(umap_coords, delta_emb, clusters, cluster_colors, title, outpath, scale=1.0):
    fig, ax = plt.subplots(figsize=(7, 6))
    # Background scatter
    for cl, color in cluster_colors.items():
        mask = (clusters == cl)
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   c=color, s=15, alpha=0.6, label=cl, zorder=2)
    # Quiver arrows
    ax.quiver(umap_coords[:, 0], umap_coords[:, 1],
              delta_emb[:, 0], delta_emb[:, 1],
              alpha=0.5, scale=scale, scale_units='xy', angles='xy',
              width=0.003, color='black', zorder=3)
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.legend(loc='upper right', fontsize=9, markerscale=1.5)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    return outpath

for cond in ["MAFB_KO", "MAFB_OE"]:
    de = delta_emb_dict[cond]
    # Normalize arrows for visibility
    arrow_scale = np.percentile(np.abs(de), 95) * 10 + 1e-9
    out_png = os.path.join(OUTDIR_04, f"04B_quiver_{cond}.png")
    out_svg = os.path.join(OUTDIR_04, f"04B_quiver_{cond}.svg")
    plot_quiver(umap, de, clusters, CLUSTER_COLORS,
                f"MAFB Perturbation: {cond}\n(delta_embedding arrows)",
                out_png, scale=arrow_scale)
    plot_quiver(umap, de, clusters, CLUSTER_COLORS,
                f"MAFB Perturbation: {cond}\n(delta_embedding arrows)",
                out_svg, scale=arrow_scale)
    log(f"  [PASS] Quiver plot saved: {cond} [{elapsed()}]")

# WT baseline (should be ~0 arrows)
de_wt = delta_emb_dict["WT"]
out_wt = os.path.join(OUTDIR_04, "04B_quiver_WT.png")
plot_quiver(umap, de_wt, clusters, CLUSTER_COLORS,
            "WT baseline (delta_embedding ≈ 0)", out_wt, scale=1.0)
log(f"  [PASS] WT quiver saved [{elapsed()}]")

# ── 04C: Gene-level response (KO vs WT, OE vs WT) ────────────
log(f"\n--- 04C: Gene-level Response ---")

wt_dx = delta_X_dict["WT"]  # should be 0

for cond in ["MAFB_KO", "MAFB_OE"]:
    dx = delta_X_dict[cond]
    diff = dx - wt_dx  # delta_X relative to WT
    
    # Per-gene mean delta across all cells
    gene_mean = diff.mean(axis=0)
    gene_abs  = np.abs(diff).mean(axis=0)
    
    gene_df = pd.DataFrame({
        'gene': gene_ids,
        'mean_delta_vs_WT': gene_mean,
        'mean_abs_delta_vs_WT': gene_abs,
    }).sort_values('mean_abs_delta_vs_WT', ascending=False)
    
    out_gene = os.path.join(OUTDIR_04, f"04C_gene_response_{cond}_vs_WT.csv")
    gene_df.to_csv(out_gene, index=False)
    log(f"  [PASS] Gene response saved: {cond} vs WT")
    
    # Top 20 up/down
    top_up   = gene_df.nlargest(20, 'mean_delta_vs_WT')[['gene','mean_delta_vs_WT']].to_string(index=False)
    top_down = gene_df.nsmallest(20, 'mean_delta_vs_WT')[['gene','mean_delta_vs_WT']].to_string(index=False)
    log(f"\n  {cond} vs WT — Top 20 UP genes:\n{top_up}")
    log(f"\n  {cond} vs WT — Top 20 DOWN genes:\n{top_down}")

# ── 04D: Per-cluster delta_X heatmap (top 30 genes) ──────────
log(f"\n--- 04D: Per-cluster delta_X summary ---")

for cond in ["MAFB_KO", "MAFB_OE"]:
    dx = delta_X_dict[cond]
    rows = []
    for cl in ["Mono", "IFN_TAM", "LA_TAM"]:
        mask = (clusters == cl)
        cl_mean = dx[mask].mean(axis=0)
        rows.append(cl_mean)
    
    cl_df = pd.DataFrame(rows, index=["Mono","IFN_TAM","LA_TAM"], columns=gene_ids)
    
    # Top 30 genes by max abs across clusters
    top_genes = cl_df.abs().max(axis=0).nlargest(30).index.tolist()
    cl_top = cl_df[top_genes]
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(14, 4))
    sns.heatmap(cl_top, cmap='RdBu_r', center=0, ax=ax,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'mean delta_X'})
    ax.set_title(f"{cond}: Per-cluster mean delta_X (top 30 genes)", fontsize=12)
    ax.set_xlabel("Gene"); ax.set_ylabel("Cluster")
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.tight_layout()
    out_hm = os.path.join(OUTDIR_04, f"04D_cluster_deltaX_heatmap_{cond}.png")
    plt.savefig(out_hm, dpi=150, bbox_inches='tight')
    plt.close()
    log(f"  [PASS] Cluster heatmap saved: {cond} [{elapsed()}]")
    
    # Save full per-cluster table
    cl_df.to_csv(os.path.join(OUTDIR_04, f"04D_cluster_deltaX_{cond}.csv"))

# ── 04E: Fate score violin plot ───────────────────────────────
log(f"\n--- 04E: Fate score violin plot ---")

fate_long = []
for cond in CONDITIONS:
    for cl in ["Mono", "IFN_TAM", "LA_TAM"]:
        mask = (clusters == cl)
        for v in fate_scores[cond][mask]:
            fate_long.append({'condition': cond, 'cluster': cl, 'fate_score': v})
fate_long_df = pd.DataFrame(fate_long)

fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
for i, cl in enumerate(["Mono", "IFN_TAM", "LA_TAM"]):
    ax = axes[i]
    sub = fate_long_df[fate_long_df['cluster']==cl]
    sns.violinplot(data=sub, x='condition', y='fate_score', ax=ax,
                   palette=['#E74C3C','#3498DB','#95A5A6'],
                   order=["MAFB_KO","MAFB_OE","WT"], inner='box')
    ax.set_title(f"{cl}", fontsize=12)
    ax.set_xlabel("Condition"); ax.set_ylabel("LA_TAM fate score" if i==0 else "")
    ax.tick_params(axis='x', rotation=30)

plt.suptitle("LA_TAM Fate Score by Condition and Cluster", fontsize=13, y=1.02)
plt.tight_layout()
out_violin = os.path.join(OUTDIR_04, "04E_fate_score_violin.png")
plt.savefig(out_violin, dpi=150, bbox_inches='tight')
plt.close()
log(f"  [PASS] Fate score violin saved [{elapsed()}]")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 04 COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

audit_04 = os.path.join(OUTDIR_04, "04_scoring_figures_log.txt")
with open(audit_04,"w") as f: f.write("\n".join(lines)+"\n")
print(f"\nSTEP 04 COMPLETE — {total:.1f}s")
