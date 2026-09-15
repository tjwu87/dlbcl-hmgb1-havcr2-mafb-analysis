# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python
# step09_fig23_deg_ora.py
# Fig23: MAFB KO vs WT DEG (imputed_count vs simulated_count) + ORA enrichment bubble
# Data: oracle_WT.adata.layers['imputed_count'] vs oracle_KO.adata.layers['simulated_count']
# Test: Wilcoxon signed-rank (paired per cell), BH correction
# Foreground: padj < 0.05 & |log2FC| > 0.3, per subtype, up/down
# Background: HVG3000 all genes
# Plot style: Fig07 enrichment bubble framework

import os, sys, time, warnings, pickle
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
warnings.filterwarnings("ignore")

START = time.time()
OUTDIR = translate("/mnt/results/GSE182434/celloracle_paper_final_20260331_085223")
OUT    = f"{OUTDIR}/09_figures"
DOUT   = f"{OUTDIR}/09_figures/data"
os.makedirs(OUT,  exist_ok=True)
os.makedirs(DOUT, exist_ok=True)

def elapsed():
    return f"{(time.time()-START)/60:.1f}min"

print(f"[{elapsed()}] step09_fig23 started")

# ── Load oracle objects ───────────────────────────────────────────────────────
def load_pkl(path):
    with open(path, "rb") as f:
        return pickle.load(f)

print(f"[{elapsed()}] Loading oracle WT...")
oracle_wt = load_pkl(f"{OUTDIR}/03_perturbation/WT/oracle_WT.pkl")
print(f"[{elapsed()}] Loading oracle KO...")
oracle_ko = load_pkl(f"{OUTDIR}/03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl")

adata_wt = oracle_wt.adata
adata_ko = oracle_ko.adata

gene_names  = np.array(adata_wt.var_names)
subtypes    = adata_wt.obs["mac_subtype"].values
subtypes_order = ["Mono", "IFN_TAM", "LA_TAM"]

# Extract expression matrices
X_wt = adata_wt.layers["imputed_count"]
X_ko = adata_ko.layers["simulated_count"]
if hasattr(X_wt, "toarray"): X_wt = X_wt.toarray()
if hasattr(X_ko, "toarray"): X_ko = X_ko.toarray()

print(f"[{elapsed()}] WT imputed_count: {X_wt.shape}, KO simulated_count: {X_ko.shape}")
print(f"  WT range: [{X_wt.min():.4f}, {X_wt.max():.4f}]")
print(f"  KO range: [{X_ko.min():.4f}, {X_ko.max():.4f}]")

# ── Per-subtype Wilcoxon signed-rank DEG ─────────────────────────────────────
print(f"[{elapsed()}] Computing per-subtype DEG (Wilcoxon signed-rank + BH)...")

all_deg_rows = []
deg_per_subtype = {}  # {subtype: {"up": [genes], "down": [genes]}}

for st in subtypes_order:
    mask = subtypes == st
    n_cells = mask.sum()
    print(f"  {st}: {n_cells} cells")

    wt_sub = X_wt[mask, :]   # (n_cells, n_genes)
    ko_sub = X_ko[mask, :]

    # log2FC: mean(KO) - mean(WT) in log space (already log1p-like imputed)
    mean_wt = wt_sub.mean(0)
    mean_ko = ko_sub.mean(0)
    log2fc  = mean_ko - mean_wt   # difference in log-space ≈ log2FC proxy

    # Wilcoxon signed-rank test per gene (paired: each cell is a pair)
    pvals = np.ones(len(gene_names))
    for g in range(len(gene_names)):
        wt_g = wt_sub[:, g]
        ko_g = ko_sub[:, g]
        diff = ko_g - wt_g
        if np.any(diff != 0):
            try:
                _, pvals[g] = stats.wilcoxon(wt_g, ko_g, alternative="two-sided",
                                              zero_method="wilcox")
            except Exception:
                pvals[g] = 1.0

    # BH correction
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")

    # Build DEG table
    df_st = pd.DataFrame({
        "gene":    gene_names,
        "subtype": st,
        "mean_wt": mean_wt,
        "mean_ko": mean_ko,
        "log2FC":  log2fc,
        "p_value": pvals,
        "padj":    padj,
    })
    df_st["direction"] = "ns"
    df_st.loc[(df_st["padj"] < 0.05) & (df_st["log2FC"] >  0.3), "direction"] = "up"
    df_st.loc[(df_st["padj"] < 0.05) & (df_st["log2FC"] < -0.3), "direction"] = "down"

    n_up   = (df_st["direction"] == "up").sum()
    n_down = (df_st["direction"] == "down").sum()
    print(f"    up={n_up}, down={n_down}, total_sig={(df_st['padj']<0.05).sum()}")

    all_deg_rows.append(df_st)
    deg_per_subtype[st] = {
        "up":   df_st[df_st["direction"] == "up"]["gene"].tolist(),
        "down": df_st[df_st["direction"] == "down"]["gene"].tolist(),
    }

df_deg_all = pd.concat(all_deg_rows, ignore_index=True)
deg_path = f"{DOUT}/fig23_deg_imputed_vs_simulated.csv"
df_deg_all.to_csv(deg_path, index=False)
print(f"[{elapsed()}] DEG table saved: {deg_path} ({len(df_deg_all)} rows)")

# ── ORA enrichment per subtype × direction ───────────────────────────────────
print(f"[{elapsed()}] Running ORA enrichment (gseapy enrichr)...")

gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]
background = gene_names.tolist()

all_enr_rows = []

for st in subtypes_order:
    for direction in ["up", "down"]:
        fg_genes = deg_per_subtype[st][direction]
        print(f"  {st} {direction}: {len(fg_genes)} foreground genes")

        if len(fg_genes) < 3:
            print(f"    SKIP: foreground < 3 genes")
            continue

        for gs in gene_sets:
            try:
                enr = gp.enrichr(
                    gene_list=fg_genes,
                    gene_sets=gs,
                    background=background,
                    outdir=None,
                    verbose=False,
                )
                df_enr = enr.results.copy()
                if len(df_enr) == 0:
                    continue
                df_enr["subtype"]   = st
                df_enr["direction"] = direction
                df_enr["gene_set"]  = gs
                df_enr["n_fg_genes"] = len(fg_genes)
                df_enr["n_bg_genes"] = len(background)
                all_enr_rows.append(df_enr)
            except Exception as e:
                print(f"    ERROR {gs}: {e}")

if all_enr_rows:
    df_enr_all = pd.concat(all_enr_rows, ignore_index=True)
    enr_path = f"{DOUT}/fig23_ora_enrichment_results.csv"
    df_enr_all.to_csv(enr_path, index=False)
    print(f"[{elapsed()}] Enrichment saved: {enr_path} ({len(df_enr_all)} rows)")
    print(f"  Columns: {df_enr_all.columns.tolist()}")
    sig = df_enr_all[df_enr_all["Adjusted P-value"] < 0.05]
    print(f"  Significant (padj<0.05): {len(sig)} terms")
    print(sig.groupby(["subtype","direction","gene_set"]).size().to_string())
else:
    print("WARNING: No enrichment results — all foreground sets < 3 genes")
    df_enr_all = pd.DataFrame()

print(f"[{elapsed()}] ORA complete")

# ── Fig23: Enrichment bubble plot (Fig07 framework) ──────────────────────────
print(f"[{elapsed()}] Drawing Fig23...")

if len(df_enr_all) == 0:
    print("No enrichment data — creating placeholder Fig23")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5,
            "No significant ORA terms\n(foreground sets < 3 genes per subtype/direction)",
            ha="center", va="center", fontsize=14, transform=ax.transAxes)
    ax.set_title("Fig23: MAFB KO vs WT — ORA Enrichment\n(imputed_count vs simulated_count)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig23_mafb_ko_ora_enrichment_bubble.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/fig23_mafb_ko_ora_enrichment_bubble.svg", bbox_inches="tight")
    plt.close()
    print(f"[{elapsed()}] Fig23 placeholder saved")
else:
    # Select top terms per subtype × direction × gene_set
    top_terms = []
    for st in subtypes_order:
        for direction in ["up", "down"]:
            for gs in gene_sets:
                sub = df_enr_all[
                    (df_enr_all["subtype"] == st) &
                    (df_enr_all["direction"] == direction) &
                    (df_enr_all["gene_set"] == gs)
                ].copy()
                sub = sub[sub["Adjusted P-value"] < 0.05]
                if len(sub) == 0:
                    continue
                sub = sub.nsmallest(5, "Adjusted P-value")
                sub["panel"] = f"{st}\n({direction})"
                top_terms.append(sub)

    if not top_terms:
        print("No significant terms at padj<0.05 — relaxing to padj<0.2")
        for st in subtypes_order:
            for direction in ["up", "down"]:
                for gs in gene_sets:
                    sub = df_enr_all[
                        (df_enr_all["subtype"] == st) &
                        (df_enr_all["direction"] == direction) &
                        (df_enr_all["gene_set"] == gs)
                    ].copy()
                    if len(sub) == 0:
                        continue
                    sub = sub.nsmallest(3, "Adjusted P-value")
                    sub["panel"] = f"{st}\n({direction})"
                    top_terms.append(sub)

    df_plot = pd.concat(top_terms, ignore_index=True)
    df_plot["-log10_padj"] = -np.log10(df_plot["Adjusted P-value"] + 1e-300)
    df_plot["overlap_frac"] = df_plot["n_fg_genes"] / df_plot["n_bg_genes"]
    df_plot["Term_short"] = df_plot["Term"].apply(
        lambda x: x[:55] + "..." if len(x) > 55 else x)

    # One panel per subtype × direction
    panels = [f"{st}\n({d})" for st in subtypes_order for d in ["up", "down"]]
    panels_present = [p for p in panels if p in df_plot["panel"].values]

    n_panels = len(panels_present)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    lib_colors = {
        "GO_Biological_Process_2023": "#4DBBD5",
        "KEGG_2021_Human":            "#E64B35",
        "MSigDB_Hallmark_2020":       "#00A087",
    }
    subtype_colors = {"Mono": "#E64B35", "IFN_TAM": "#4DBBD5", "LA_TAM": "#00A087"}

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 7, nrows * 5 + 1),
                             squeeze=False)
    fig.patch.set_facecolor("white")

    for idx, panel in enumerate(panels_present):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df_plot[df_plot["panel"] == panel].copy()
        sub = sub.sort_values("-log10_padj", ascending=True)

        for gs in gene_sets:
            gs_sub = sub[sub["gene_set"] == gs]
            if len(gs_sub) == 0:
                continue
            color = lib_colors.get(gs, "#888888")
            ax.scatter(
                gs_sub["-log10_padj"],
                gs_sub["Term_short"],
                s=gs_sub["overlap_frac"] * 3000 + 30,
                c=color, alpha=0.85, linewidths=0.5,
                edgecolors="white", label=gs.replace("_2023","").replace("_2021_Human","").replace("_2020",""),
                zorder=3
            )

        ax.axvline(x=-np.log10(0.05), color="gray", linestyle="--",
                   linewidth=1, alpha=0.7)
        ax.set_xlabel("-log10(Adjusted P-value)", fontsize=10)
        ax.set_title(panel.replace("\n", " — "), fontsize=11, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.3)
        if idx == 0:
            ax.legend(title="Gene Set", fontsize=8, loc="lower right", framealpha=0.8)

    # Turn off empty panels
    for idx in range(len(panels_present), nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    plt.suptitle(
        "Fig23: MAFB KO vs WT — ORA Enrichment per Subtype\n"
        "(WT: imputed_count | KO: simulated_count | Wilcoxon signed-rank, BH | "
        "padj<0.05, |log2FC|>0.3 | background=HVG3000)",
        fontsize=12, fontweight="bold", y=1.01
    )
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig23_mafb_ko_ora_enrichment_bubble.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{OUT}/fig23_mafb_ko_ora_enrichment_bubble.svg", bbox_inches="tight")
    plt.close()
    print(f"[{elapsed()}] Fig23 saved ({len(df_plot)} terms plotted)")

print(f"\n[{elapsed()}] step09_fig23 COMPLETE")
