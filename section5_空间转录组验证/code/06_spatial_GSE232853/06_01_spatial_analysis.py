# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

"""
reproduce_all_figures.py
Local reproduction of Fig1–Fig9 + S5A/S5B/S5C
GSE232853 · GeoMx DSP · DLBCL macrophage spatial transcriptomics
Requirements: pandas numpy matplotlib seaborn scipy statsmodels gseapy
"""

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import TwoSlopeNorm, ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.stats import (mannwhitneyu, pearsonr, spearmanr,
                          gaussian_kde, zscore, t as t_dist)
from statsmodels.stats.multitest import multipletests
import seaborn as sns
warnings.filterwarnings("ignore")

OUT = translate(r"D:\bulk-download\GSE232853_v2\figures_reproduced")
os.makedirs(OUT, exist_ok=True)

# ─────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────
MASK_COLORS   = {"CD20": "#2166AC", "CD68": "#D6604D"}
TISSUE_COLORS = {"DLBCL": "#B2182B", "Normal": "#4393C3"}

def sig_stars(p):
    if p < 0.0001: return "****"
    elif p < 0.001: return "***"
    elif p < 0.01:  return "**"
    elif p < 0.05:  return "*"
    return "ns"

def regression_ci(x, y, x_line, ci=0.95):
    n = len(x)
    m, b = np.polyfit(x, y, 1)
    y_hat  = m * x + b
    se     = np.sqrt(np.sum((y - y_hat)**2) / (n - 2))
    ss_x   = np.sum((x - np.mean(x))**2)
    t_val  = t_dist.ppf((1 + ci) / 2, df=n - 2)
    y_line = m * x_line + b
    se_ln  = se * np.sqrt(1/n + (x_line - np.mean(x))**2 / ss_x)
    return y_line, t_val * se_ln

def clean_term(t):
    return (t.replace("_", " ").title()
             .replace("Emt","EMT").replace("Pi3K","PI3K")
             .replace("Mtorc1","mTORC1").replace("Il-","IL-")
             .replace("Il-6","IL-6").replace("Tnf-Alpha","TNF-α")
             .replace("Dna","DNA").replace("G2-M","G2/M")
             .replace("  "," "))

def savefig(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {name}")

# ─────────────────────────────────────────────
# Load shared data
# ─────────────────────────────────────────────
print("Loading data...")

# Meta (572 filtered ROIs, v2)
meta = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/meta_filtered.csv"))
# Normalise column names defensively
meta.columns = [c.strip() for c in meta.columns]

# UMAP (572 ROIs, from v2 pipeline)
umap_df = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/task2_UMAP_coordinates.csv"))
umap_df.columns = [c.strip() for c in umap_df.columns]

# PCA (578 samples, for S5B)
pca_df = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/task2_PCA_coordinates.csv"))
pca_df.columns = [c.strip() for c in pca_df.columns]



# Wilcoxon stats for Fig2
wilcox2 = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig2_wilcoxon_HMGB1_HAVCR2_MAFB.csv"))

# DEA for Fig3
dea = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig3_DEA_CD68_scanpy.csv"))
dea.columns = [c.strip() for c in dea.columns]

# Paired ROI table for Fig4
paired = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig4_paired_ROI_16gene_LATAM.csv"))

# Correlation stats for Fig4
corr4 = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig4_correlation_stats.csv"))

# GSEA results for Fig5
gsea_df = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig5_GSEA_Hallmark_LATAM_Rich_vs_Poor_CD20.csv"))

# DLBCL LA_TAM group assignment for Fig6
group_df = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig5_DLBCL_LATAM_group_assignment.csv"))
group_df.columns = [c.strip() for c in group_df.columns]
roi_to_group = group_df.set_index("Unique_ROI")["LA_TAM_Group"].to_dict()

# ssGSEA scores for Fig6
scores_wide = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig6_ssGSEA_Hallmark_CD20_DLBCL.csv"), index_col=0)
scores_wide["LA_TAM_Group"] = scores_wide.index.map(roi_to_group)
path_cols   = [c for c in scores_wide.columns if c != "LA_TAM_Group"]
rich_scores = scores_wide[scores_wide["LA_TAM_Group"]=="LA_TAM_Rich"][path_cols].astype(float)
poor_scores = scores_wide[scores_wide["LA_TAM_Group"]=="LA_TAM_Poor"][path_cols].astype(float)

# Master spatial co-evolution table for Fig7/8/9
master = pd.read_csv(translate("D:/bulk-download/GSE232853_v2/fig7_master_spatial_coevolution.csv"))

print("  All files loaded.\n")

# ═══════════════════════════════════════════════════════════════
# FIG 1 — Dual UMAP (Mask | Tissue_Type)
# ═══════════════════════════════════════════════════════════════
print("Fig1: UMAP...")

# Detect UMAP column names flexibly
u1_col = [c for c in umap_df.columns if "umap" in c.lower() and "1" in c][0]
u2_col = [c for c in umap_df.columns if "umap" in c.lower() and "2" in c][0]
mask_col   = [c for c in umap_df.columns if "mask" in c.lower()][0]
tissue_col = [c for c in umap_df.columns if "tissue" in c.lower()][0]

fig, axes = plt.subplots(1, 2, figsize=(14, 9))
fig.patch.set_facecolor("#FFFFFF")

for ax, (col, cmap_dict, title) in zip(axes, [
    (mask_col,   MASK_COLORS,   "Cell Compartment"),
    (tissue_col, TISSUE_COLORS, "Tissue Type"),
]):
    ax.set_facecolor("#FFFFFF")
    ax.scatter(umap_df[u1_col], umap_df[u2_col],
               c="#DDDDDD", s=8, alpha=0.3, linewidths=0, zorder=1, rasterized=True)
    for grp in sorted(umap_df[col].unique()):
        idx = umap_df[col] == grp
        color = cmap_dict.get(grp, "#888888")
        ax.scatter(umap_df.loc[idx, u1_col], umap_df.loc[idx, u2_col],
                   c=color, s=18, alpha=0.75, linewidths=0, zorder=2, rasterized=True)
    handles = [mpatches.Patch(color=cmap_dict.get(g,"#888"),
                              label=f"{g}  (n={(umap_df[col]==g).sum()})")
               for g in sorted(umap_df[col].unique())]
    # 版面窄（3.12 in）时图例放坐标区内会盖住 UMAP，改为外置下方两列
    ax.legend(handles=handles, fontsize=10, framealpha=0.92,
              edgecolor="#CCCCCC", loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=2, borderaxespad=0.0)
    ax.set_xlabel("UMAP 1", fontsize=11); ax.set_ylabel("UMAP 2", fontsize=11)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=8)
    ax.tick_params(labelsize=9)
    for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")

n_rois = len(umap_df)
fig.suptitle(f"GeoMx DSP — UMAP  |  GSE232853 · DLBCL · {n_rois} ROIs",
             fontsize=12, fontweight="bold", y=0.99)
plt.tight_layout()
savefig(fig, "fig1_UMAP_dual")



# ═══════════════════════════════════════════════════════════════
# FIG 3 — Volcano: CD68+ DLBCL vs Normal DEA
# ═══════════════════════════════════════════════════════════════
print("Fig3: Volcano DEA...")

# Detect column names flexibly
lfc_col  = [c for c in dea.columns if "logfold" in c.lower() or "log2fc" in c.lower()
            or "lfc" in c.lower() or "logfc" in c.lower()][0]
pval_col = [c for c in dea.columns if "padj" in c.lower() or "fdr" in c.lower()
            or "pval" in c.lower() or "p_val" in c.lower()][0]
gene_col = [c for c in dea.columns if "gene" in c.lower() or "names" in c.lower()][0]

dea = dea.dropna(subset=[lfc_col, pval_col])
dea["neg_log10_p"] = -np.log10(dea[pval_col].clip(lower=1e-300))
dea["color"] = "#AAAAAA"
dea.loc[(dea[lfc_col] >  1) & (dea[pval_col] < 0.05), "color"] = "#B2182B"
dea.loc[(dea[lfc_col] < -1) & (dea[pval_col] < 0.05), "color"] = "#4393C3"

n_up   = ((dea[lfc_col] >  1) & (dea[pval_col] < 0.05)).sum()
n_dn   = ((dea[lfc_col] < -1) & (dea[pval_col] < 0.05)).sum()

fig, ax = plt.subplots(figsize=(9, 7))
fig.patch.set_facecolor("#FFFFFF"); ax.set_facecolor("#FFFFFF")
ax.scatter(dea[lfc_col], dea["neg_log10_p"], c=dea["color"],
           s=6, alpha=0.55, linewidths=0, rasterized=True)
ax.axvline( 1, color="#B2182B", lw=1.0, ls="--", alpha=0.6)
ax.axvline(-1, color="#4393C3", lw=1.0, ls="--", alpha=0.6)
ax.axhline(-np.log10(0.05), color="#888888", lw=1.0, ls="--", alpha=0.6)

# Label top genes
top_up = dea[dea["color"]=="#B2182B"].nlargest(8, "neg_log10_p")
top_dn = dea[dea["color"]=="#4393C3"].nlargest(8, "neg_log10_p")
for _, r in pd.concat([top_up, top_dn]).iterrows():
    ax.text(r[lfc_col], r["neg_log10_p"]+0.3, r[gene_col],
            fontsize=6.5, ha="center", color=r["color"], fontweight="bold")

ax.text(0.97, 0.97, f"Up in DLBCL: {n_up}", transform=ax.transAxes,
        ha="right", va="top", fontsize=9, color="#B2182B", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#B2182B", alpha=0.85))
ax.text(0.03, 0.97, f"Down in DLBCL: {n_dn}", transform=ax.transAxes,
        ha="left", va="top", fontsize=9, color="#4393C3", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#4393C3", alpha=0.85))

ax.set_xlabel("log2 Fold Change (DLBCL vs Normal)", fontsize=11)
ax.set_ylabel("−log10(padj)", fontsize=11)
ax.set_title("Differential Expression: CD68+ Macrophages\nDLBCL vs Normal · GSE232853",
             fontsize=11, fontweight="bold", pad=10)
ax.tick_params(labelsize=9)
for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")
plt.tight_layout()
savefig(fig, "fig3_volcano_DEA_CD68")

# ────────────────────────────────────────────────────────────
# TABLE 4 — S11 Key DEA Genes (HAVCR2, MAFB, HMGB1)
#   CD68+ DLBCL vs Normal
#   Values: gene, log2FC (2dp), padj, Direction
# ────────────────────────────────────────────────────────────
print("\n[4/5] Exporting S11 key DEA genes...")

_lfc_col  = [c for c in dea.columns if any(k in c.lower() for k in
             ["logfold", "log2fc", "logfc", "lfc"])][0]
_pval_col = [c for c in dea.columns if any(k in c.lower() for k in
             ["padj", "pval_adj", "p_val_adj", "fdr", "pvals_adj"])][0]
_gene_col = [c for c in dea.columns if any(k in c.lower() for k in
             ["gene", "names"])][0]

_key_genes = ["HAVCR2", "MAFB", "HMGB1"]
_dea_sub = dea[dea[_gene_col].isin(_key_genes)].copy()

if len(_dea_sub) == 0:
    # Try case-insensitive match
    _dea_sub = dea[dea[_gene_col].str.upper().isin([g.upper() for g in _key_genes])].copy()

_dea_sub["log2FC_2dp"]    = _dea_sub[_lfc_col].apply(lambda x: round(float(x), 2))
_dea_sub["padj_formatted"] = _dea_sub[_pval_col].apply(
    lambda x: "< 0.001" if float(x) < 0.001 else f"{float(x):.3e}")
_dea_sub["Direction"] = _dea_sub[_lfc_col].apply(
    lambda x: "Up_in_DLBCL" if float(x) > 0 else "Down_in_DLBCL")

_dea_export_cols = [_gene_col, _lfc_col, "log2FC_2dp",
                    _pval_col, "padj_formatted", "Direction"]
_dea_export_cols = [c for c in _dea_export_cols if c in _dea_sub.columns]
_dea_sub[_dea_export_cols].to_csv(
    os.path.join(OUT, "dea_key_genes_for_results.csv"), index=False)
print(f"  ✓ dea_key_genes_for_results.csv  ({len(_dea_sub)} genes found)")

if len(_dea_sub) < len(_key_genes):
    _found = _dea_sub[_gene_col].tolist()
    _missing_genes = [g for g in _key_genes if g not in _found]
    print(f"  [WARNING] Not found in DEA: {_missing_genes}")
    print(f"  DEA gene column sample: {dea[_gene_col].head(10).tolist()}")
# ═══════════════════════════════════════════════════════════════
# FIG 4 — Scatter: HMGB1 vs HAVCR2 & LA_TAM_Score
# ═══════════════════════════════════════════════════════════════
print("Fig4: Scatter paired ROI...")

# Detect column names
hmgb1_col  = [c for c in paired.columns if "hmgb1" in c.lower()][0]
havcr2_col = [c for c in paired.columns if "havcr2" in c.lower()][0]
latam_col  = [c for c in paired.columns if "latam" in c.lower() or "la_tam" in c.lower()][0]
tissue_c   = [c for c in paired.columns if "tissue" in c.lower()][0]

b_hmgb1  = paired[hmgb1_col].values
m_havcr2 = paired[havcr2_col].values
la_tam   = paired[latam_col].values
tissue   = paired[tissue_c].values
dlbcl_idx = tissue == "DLBCL"
norm_idx  = tissue == "Normal"

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor("#FFFFFF")

panel_cfg = [
    (axes[0], m_havcr2,
     "B cell HMGB1 (Q3 log-norm)",
     "Macrophage HAVCR2 / TIM-3 (Q3 log-norm)",
     "Fig A: B-cell HMGB1 vs Macrophage HAVCR2"),
    (axes[1], la_tam,
     "B cell HMGB1 (Q3 log-norm)",
     "Macrophage LA_TAM_Score\n(mean of 16 genes: PTGDS, CCL18, APOE, CHI3L1…)",
     "Fig B: B-cell HMGB1 vs LA_TAM_Score"),
]

for ax, y_all, xlabel, ylabel, title in panel_cfg:
    ax.set_facecolor("#FFFFFF")
    x_line = np.linspace(b_hmgb1.min()-0.1, b_hmgb1.max()+0.1, 200)
    stat_y = {"DLBCL": 0.97, "Normal": 0.83}
    for tissue_label, color in TISSUE_COLORS.items():
        idx = dlbcl_idx if tissue_label=="DLBCL" else norm_idx
        x_t = b_hmgb1[idx]; y_t = y_all[idx]
        ax.scatter(x_t, y_t, c=color, s=22, alpha=0.60, linewidths=0.3,
                   edgecolors="white", label=f"{tissue_label} (n={idx.sum()})",
                   zorder=3, rasterized=True)
        y_line, ci_band = regression_ci(x_t, y_t, x_line)
        ax.plot(x_line, y_line, color=color, lw=2.0, zorder=4)
        ax.fill_between(x_line, y_line-ci_band, y_line+ci_band,
                        color=color, alpha=0.12, zorder=2)
        r, pr   = pearsonr(x_t, y_t)
        rho, ps = spearmanr(x_t, y_t)
        ax.text(0.03, stat_y[tissue_label],
                f"{tissue_label}\nr = {r:+.3f},  p = {pr:.2e}\nρ = {rho:+.3f}",
                transform=ax.transAxes, fontsize=8.5, va="top", color=color,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.28", fc="white",
                          ec=color, alpha=0.88, linewidth=0.9))
    ax.set_xlabel(xlabel, fontsize=10); ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=8)
    ax.legend(fontsize=9, framealpha=0.9, edgecolor="#CCCCCC",
              loc="lower right", markerscale=1.3)
    ax.tick_params(labelsize=9)
    for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")

fig.suptitle("Paired ROI Co-expression Analysis · 228 Matched CD20+/CD68+ ROIs\n"
             "GSE232853 · Q3-Normalized · 95% CI regression bands",
             fontsize=9, fontweight="bold", y=1.04)
plt.tight_layout()
savefig(fig, "fig4_scatter_HMGB1_HAVCR2_LATAM")


print("\n[1/5] Exporting Fig4/5E correlation stats...")

_hmgb1_col  = [c for c in paired.columns if "hmgb1"  in c.lower()][0]
_havcr2_col = [c for c in paired.columns if "havcr2" in c.lower()][0]
_latam_col  = [c for c in paired.columns if "latam"  in c.lower() or "la_tam" in c.lower()][0]
_tissue_c   = [c for c in paired.columns if "tissue" in c.lower()][0]

_b_hmgb1  = paired[_hmgb1_col].astype(float).values
_m_havcr2 = paired[_havcr2_col].astype(float).values
_la_tam   = paired[_latam_col].astype(float).values
_tissue   = paired[_tissue_c].astype(str).values

_dlbcl = _tissue == "DLBCL"
_norm  = _tissue == "Normal"

corr_rows = []
for _aname, _y in [("HMGB1_vs_HAVCR2",       _m_havcr2),
                   ("HMGB1_vs_LA_TAM_score",  _la_tam)]:
    for _tname, _idx in [("All",    np.ones(len(paired), dtype=bool)),
                         ("DLBCL",  _dlbcl),
                         ("Normal", _norm)]:
        _x = _b_hmgb1[_idx]
        _y_sub = _y[_idx]
        _n = int(_idx.sum())
        if _n < 3:
            continue
        _r,   _pr  = pearsonr(_x, _y_sub)
        _rho, _ps  = spearmanr(_x, _y_sub)
        corr_rows.append({
            "analysis":     _aname,
            "tissue":       _tname,
            "n":            _n,
            "pearson_r":    round(float(_r),   4),
            "pearson_p":    float(_pr),
            "spearman_rho": round(float(_rho), 4),
            "spearman_p":   float(_ps),
            "x_col":        _hmgb1_col,
            "y_col":        _havcr2_col if _aname == "HMGB1_vs_HAVCR2" else _latam_col,
        })

_corr_df = pd.DataFrame(corr_rows)
_corr_df.to_csv(os.path.join(OUT, "fig4_scatter_correlation_values.csv"), index=False)
print(f"  ✓ fig4_scatter_correlation_values.csv  ({len(_corr_df)} rows)")

# ────────────────────────────────────────────────────────────
# TABLE 2 — S9 Paired ROI Summary Stats
#   HMGB1 / HAVCR2 / MAFB: CD20 vs CD68 within paired ROIs
#   Groups: All / DLBCL / Normal
#   Values: n, median, q25, q75, median_diff, MWU p
# ────────────────────────────────────────────────────────────
print("\n[2/5] Exporting S9 paired ROI summary stats...")

# Auto-detect CD20/CD68 paired columns for each gene
def _find_paired_col(df, mask, gene):
    """Find column like CD20_HMGB1 or HMGB1_CD20 (case-insensitive)."""
    for c in df.columns:
        cl = c.lower()
        if mask.lower() in cl and gene.lower() in cl:
            return c
    return None

_s9_genes = ["HMGB1", "HAVCR2", "MAFB"]
_s9_rows  = []
_s9_missing = []

for _gene in _s9_genes:
    _col_cd20 = _find_paired_col(paired, "B_", _gene)
    _col_cd68 = _find_paired_col(paired, "Mac", _gene)

    if _col_cd20 is None or _col_cd68 is None:
        _s9_missing.append(_gene)
        continue

    for _tname, _tmask in [("All",    np.ones(len(paired), dtype=bool)),
                           ("DLBCL",  _dlbcl),
                           ("Normal", _norm)]:
        _sub = paired[_tmask][[_col_cd20, _col_cd68, _tissue_c]].dropna(
            subset=[_col_cd20, _col_cd68])
        _n = len(_sub)
        if _n == 0:
            continue
        _a = _sub[_col_cd20].astype(float).values
        _b = _sub[_col_cd68].astype(float).values
        _diff = _b - _a
        try:
            _, _wp = mannwhitneyu(_a, _b, alternative="two-sided")
        except Exception:
            _wp = float("nan")
        _s9_rows.append({
            "gene":                    _gene,
            "tissue":                  _tname,
            "n_pairs":                 _n,
            "median_CD20":             round(float(np.median(_a)), 4),
            "q25_CD20":                round(float(np.percentile(_a, 25)), 4),
            "q75_CD20":                round(float(np.percentile(_a, 75)), 4),
            "median_CD68":             round(float(np.median(_b)), 4),
            "q25_CD68":                round(float(np.percentile(_b, 25)), 4),
            "q75_CD68":                round(float(np.percentile(_b, 75)), 4),
            "median_diff_CD68_minus_CD20": round(float(np.median(_diff)), 4),
            "MWU_p":                   float(_wp),
            "col_CD20":                _col_cd20,
            "col_CD68":                _col_cd68,
        })

if _s9_missing:
    print(f"  [WARNING] Could not find paired columns for: {_s9_missing}")
    print(f"  Available paired columns: {[c for c in paired.columns if 'CD20' in c or 'CD68' in c]}")

if _s9_rows:
    _s9_df = pd.DataFrame(_s9_rows)
    _s9_df.to_csv(os.path.join(OUT, "paired_roi_summary_stats.csv"), index=False)
    print(f"  ✓ paired_roi_summary_stats.csv  ({len(_s9_df)} rows)")
else:
    print("  [SKIP] No paired summary rows generated — check column names above")


# ═══════════════════════════════════════════════════════════════
# FIG 5 — GSEA Joyplot (prerank, Hallmark)
# ═══════════════════════════════════════════════════════════════
print("Fig5: GSEA Joyplot...")

sig    = gsea_df[gsea_df["FDR q-val"] < 0.25].copy()
top_up = sig[sig["NES"] > 0].nlargest(10, "NES")
top_dn = sig[sig["NES"] < 0].nsmallest(8,  "NES")
plot_df = pd.concat([top_dn, top_up]).sort_values("NES", ascending=True).reset_index(drop=True)
plot_df["Term_clean"] = plot_df["Term"].apply(clean_term)
n_paths = len(plot_df)

norm_nes = Normalize(vmin=-2.5, vmax=4.0)
cmap_nes = plt.cm.RdBu_r
rng      = np.random.default_rng(42)

def fdr_stars(q):
    if q < 0.001: return "***"
    elif q < 0.01: return "**"
    elif q < 0.05: return "*"
    return "·"

fig_h = max(6.0, n_paths * 0.34)
fig, ax = plt.subplots(figsize=(15, fig_h))
fig.patch.set_facecolor("#FFFFFF"); ax.set_facecolor("#FFFFFF")
overlap = 0.55
x_range = np.linspace(-4.5, 6.5, 600)

for i, row in plot_df.iterrows():
    y_base = i * (1 - overlap)
    nes    = row["NES"]
    fdr    = max(row["FDR q-val"], 1e-10)
    color  = cmap_nes(norm_nes(nes))
    spread = np.clip(0.55 + 0.35*(-np.log10(fdr)/10), 0.30, 0.90)
    samples = rng.normal(loc=nes, scale=spread, size=800)
    from scipy.stats import gaussian_kde as gkde
    kde     = gkde(samples, bw_method=0.35)
    density = kde(x_range); density = density / density.max() * 0.85
    y_vals  = y_base + density

    ax.fill_between(x_range, y_base, y_vals, color=color, alpha=0.75, zorder=i+2)
    ax.plot(x_range, y_vals, color=color, lw=0.8, alpha=0.95, zorder=i+3)
    ax.axhline(y_base, color="#CCCCCC", lw=0.4, zorder=i+1)

    peak_idx = np.argmin(np.abs(x_range - nes))
    ax.vlines(nes, y_base, y_base + density[peak_idx],
              color="white", lw=1.5, zorder=i+4, alpha=0.9)
    ax.text(nes, y_base + density[peak_idx] + 0.10, f"{nes:.2f}",
            fontsize=11, color=color, ha="center", va="bottom", fontweight="bold",
            zorder=i+6,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
    ax.text(6.2, y_base + 0.28, fdr_stars(fdr),
            fontsize=13, color="#333333", va="center", ha="right", zorder=i+5)
    ax.text(-4.35, y_base + 0.28, row["Term_clean"],
            fontsize=12, va="center", ha="left", color="#111111",
            fontweight="bold" if abs(nes) > 2.5 else "normal", zorder=i+5)

ax.axvline(0, color="#555555", lw=1.2, ls="-", zorder=1, alpha=0.6)
y_top = n_paths * (1 - overlap) + 0.7
ax.annotate("", xy=(-2.5, y_top), xytext=(-0.3, y_top),
            arrowprops=dict(arrowstyle="->", color="#2980B9", lw=1.5))
ax.text(-2.6, y_top, "Depleted in LA_TAM_Rich",
        fontsize=14, color="#2980B9", style="italic", ha="right", va="center")
ax.annotate("", xy=(3.5, y_top), xytext=(0.3, y_top),
            arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.5))
ax.text(3.6, y_top, "Enriched in LA_TAM_Rich",
        fontsize=14, color="#C0392B", style="italic", ha="left", va="center")
ax.text(6.2, -0.15, "FDR: ***​ <0.001  ​** <0.01  * <0.05  · <0.25",
        fontsize=13, color="#666666", ha="right", va="bottom", style="italic")

sm = ScalarMappable(cmap=cmap_nes, norm=norm_nes)
sm.set_array([])
cbar_ax = fig.add_axes([0.92, 0.18, 0.016, 0.32])
cbar = plt.colorbar(sm, cax=cbar_ax, label="NES")
cbar.ax.tick_params(labelsize=8)

ax.set_xlim(-4.5, 6.5); ax.set_ylim(-0.4, y_top + 0.5)
ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=13, labelpad=8)
ax.set_yticks([])
ax.set_title("GSEA Hallmark Pathways — LA_TAM_Rich vs LA_TAM_Poor"
             "CD20⁺ Tumor B-cell Compartment · DLBCL ROIs"
             "Ranking: mean expression difference (Rich − Poor) · 1000 permutations",
             fontsize=10, fontweight="bold", pad=12)
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_edgecolor("#CCCCCC")
ax.tick_params(axis="x", labelsize=9)
plt.tight_layout(rect=[0, 0, 0.91, 1])
savefig(fig, "fig5_GSEA_joyplot")

# ────────────────────────────────────────────────────────────
# TABLE 3 — S10 GSEA Top Pathways
#   LA_TAM_Rich vs LA_TAM_Poor, Hallmark
#   All significant pathways (FDR < 0.25), sorted by NES
#   Values: Term, NES (2dp), FDR q-val, Direction
# ────────────────────────────────────────────────────────────
print("\n[3/5] Exporting S10 GSEA top pathways...")

_gsea_sig = gsea_df[gsea_df["FDR q-val"] < 0.25].copy()
_gsea_up  = _gsea_sig[_gsea_sig["NES"] > 0].nlargest(10, "NES")
_gsea_dn  = _gsea_sig[_gsea_sig["NES"] < 0].nsmallest(8,  "NES")
_gsea_out = pd.concat([_gsea_dn, _gsea_up]).sort_values("NES", ascending=False).reset_index(drop=True)

_gsea_out["Term_clean"]      = _gsea_out["Term"].apply(clean_term)
_gsea_out["NES_2dp"]         = _gsea_out["NES"].apply(lambda x: round(x, 2))
_gsea_out["FDR_formatted"]   = _gsea_out["FDR q-val"].apply(
    lambda q: "< 0.001" if q < 0.001 else f"{q:.3f}")
_gsea_out["Direction"]       = _gsea_out["NES"].apply(
    lambda x: "Enriched_LA_TAM_Rich" if x > 0 else "Enriched_LA_TAM_Poor")

_gsea_export_cols = ["Term", "Term_clean", "NES", "NES_2dp",
                     "FDR q-val", "FDR_formatted", "Direction"]
_gsea_export_cols = [c for c in _gsea_export_cols if c in _gsea_out.columns]
_gsea_out[_gsea_export_cols].to_csv(
    os.path.join(OUT, "gsea_top_pathways_for_results.csv"), index=False)
print(f"  ✓ gsea_top_pathways_for_results.csv  ({len(_gsea_out)} pathways)")

# ═══════════════════════════════════════════════════════════════
# FIG 6 — ssGSEA Joyplot (per-ROI scores, Rich vs Poor)
# ═══════════════════════════════════════════════════════════════
print("Fig6: ssGSEA Joyplot...")

# Rank pathways by mean difference
diff_rows = []
for path in path_cols:
    r = rich_scores[path].values.astype(float)
    p = poor_scores[path].values.astype(float)
    stat, pval = mannwhitneyu(r, p, alternative="two-sided")
    diff_rows.append({"Pathway": path,
                      "mean_Rich": float(np.mean(r)),
                      "mean_Poor": float(np.mean(p)),
                      "mean_diff": float(np.mean(r) - np.mean(p)),
                      "MWU_p": pval})
diff_df = pd.DataFrame(diff_rows)
_, padj, _, _ = multipletests(diff_df["MWU_p"], method="fdr_bh")
diff_df["padj"] = padj
diff_df = diff_df.sort_values("mean_diff", ascending=False)

top_up6 = diff_df.head(8)["Pathway"].tolist()
top_dn6 = diff_df.tail(7)["Pathway"].tolist()
selected = top_dn6 + top_up6
labels   = [clean_term(p) for p in selected]
n6       = len(selected)

COLORS6  = {"LA_TAM_Rich": "#C0392B", "LA_TAM_Poor": "#4393C3"}
overlap6 = 0.52

all_vals6 = []
for path in selected:
    all_vals6.extend(rich_scores[path].values.astype(float).tolist())
    all_vals6.extend(poor_scores[path].values.astype(float).tolist())
x_lo6 = np.percentile(all_vals6, 0.5)  - (np.percentile(all_vals6,99.5)-np.percentile(all_vals6,0.5))*0.12
x_hi6 = np.percentile(all_vals6, 99.5) + (np.percentile(all_vals6,99.5)-np.percentile(all_vals6,0.5))*0.12
x_range6 = np.linspace(x_lo6, x_hi6, 600)
diff_idx6 = diff_df.set_index("Pathway")

fig, ax = plt.subplots(figsize=(13, n6 * 0.80 + 2))
fig.patch.set_facecolor("#FFFFFF"); ax.set_facecolor("#FFFFFF")

for i, (path, label) in enumerate(zip(selected, labels)):
    y_base = i * (1 - overlap6)
    for grp, color in COLORS6.items():
        vals = (rich_scores[path].values.astype(float) if grp=="LA_TAM_Rich"
                else poor_scores[path].values.astype(float))
        if np.std(vals) < 1e-6: continue
        kde     = gaussian_kde(vals, bw_method=0.30)
        density = kde(x_range6); density = density / density.max() * 0.80
        y_vals  = y_base + density
        ax.fill_between(x_range6, y_base, y_vals, color=color, alpha=0.45, zorder=i*2+2)
        ax.plot(x_range6, y_vals, color=color, lw=1.0, alpha=0.9, zorder=i*2+3)
        med     = np.clip(np.median(vals), x_lo6, x_hi6)
        med_idx = np.argmin(np.abs(x_range6 - med))
        ax.vlines(med, y_base, y_base + density[med_idx],
                  color=color, lw=1.8, ls="--", zorder=i*2+4, alpha=0.9)
    ax.axhline(y_base, color="#DDDDDD", lw=0.5, zorder=i*2+1)
    ax.text(x_lo6 - 0.005, y_base + 0.25, label, fontsize=8.5, va="center",
            ha="right", color="#111111",
            fontweight="bold" if abs(diff_idx6.loc[path,"mean_diff"]) > 0.02 else "normal")
    md  = diff_idx6.loc[path,"mean_diff"]
    pv  = diff_idx6.loc[path,"padj"]
    star = "**" if pv<0.01 else ("*" if pv<0.05 else ("·" if pv<0.25 else ""))
    col_ann = "#C0392B" if md > 0 else "#4393C3"
    ax.text(x_hi6 - 0.002, y_base + 0.25, f"Δ={md:+.3f}{star}",
            fontsize=7.5, va="center", ha="right", color=col_ann, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))

ax.axvline(0, color="#888888", lw=1.0, ls="-", alpha=0.6, zorder=1)
handles6 = [mpatches.Patch(color=COLORS6["LA_TAM_Rich"], alpha=0.7,
                            label=f"LA_TAM_Rich (n={len(rich_scores)})"),
            mpatches.Patch(color=COLORS6["LA_TAM_Poor"], alpha=0.7,
                            label=f"LA_TAM_Poor (n={len(poor_scores)})"),
            plt.Line2D([0],[0], color="#555555", lw=1.5, ls="--", label="Median")]
ax.legend(handles=handles6, loc="lower right", fontsize=9,
          framealpha=0.92, edgecolor="#CCCCCC")
y_top6 = n6 * (1 - overlap6) + 0.5
mid6   = (x_lo6 + x_hi6) / 2
ax.text(mid6 - (x_hi6-x_lo6)*0.25, y_top6, "← Depleted in Rich",
        fontsize=8.5, color="#4393C3", style="italic", ha="center")
ax.text(mid6 + (x_hi6-x_lo6)*0.20, y_top6, "Enriched in Rich →",
        fontsize=8.5, color="#C0392B", style="italic", ha="center")
ax.set_xlim(x_lo6-0.01, x_hi6+0.01); ax.set_ylim(-0.3, y_top6+0.4)
ax.set_xlabel("ssGSEA Normalized Enrichment Score (NES)", fontsize=10, labelpad=8)
ax.set_yticks([])
ax.set_title("ssGSEA Hallmark Pathway Scores per ROI"
             "CD20⁺ Tumor B-cells · DLBCL · LA_TAM_Rich vs LA_TAM_Poor"
             "Δ = mean(Rich) − mean(Poor)  |  ** padj<0.01  * padj<0.05  · padj<0.25",
             fontsize=10, fontweight="bold", pad=12)
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_edgecolor("#CCCCCC")
ax.tick_params(axis="x", labelsize=9)
plt.tight_layout()
savefig(fig, "fig6_ssGSEA_joyplot")

# ═══════════════════════════════════════════════════════════════
# FIG 7 — Spatial Co-evolution Split Heatmap (DLBCL | Normal)
# ═══════════════════════════════════════════════════════════════
print("Fig7: Spatial co-evolution heatmap...")

feat_cols7 = ["CD20_Ligand","CD68_Receptor","CD68_TF","CD68_Target_Gene_Score"]
master_z   = master.copy()
for col in feat_cols7:
    master_z[col] = zscore(master[col].values)

dlbcl_z  = master_z[master_z["Tissue_Type"]=="DLBCL"].sort_values("CD20_Ligand").reset_index(drop=True)
normal_z = master_z[master_z["Tissue_Type"]=="Normal"].sort_values("CD20_Ligand").reset_index(drop=True)
n_d, n_n = len(dlbcl_z), len(normal_z)

mac_rows7      = ["CD68_Receptor","CD68_TF","CD68_Target_Gene_Score"]
row_labels_mac = ["HAVCR2\n(Receptor)", "\nMAFB(TF)", "Target Gene\nScore (23g)"]
row_labels_tum = ["HMGB1\n(Ligand)"]
norm7  = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
CMAP7  = "RdBu_r"

fig = plt.figure(figsize=(17, 6))
fig.patch.set_facecolor("#FFFFFF")
gs7 = gridspec.GridSpec(4, 4,
    height_ratios=[0.06, 3, 0.15, 1],
    width_ratios=[n_d, 0.5, n_n, 1],
    hspace=0.04, wspace=0.03, figure=fig)

for gs_col, (df_z, tissue, n_rois) in enumerate(
        [(dlbcl_z,"DLBCL",n_d),(normal_z,"Normal",n_n)]):
    c  = 0 if tissue=="DLBCL" else 2
    tc = TISSUE_COLORS[tissue]
    ax_bar   = fig.add_subplot(gs7[0, c])
    ax_mac   = fig.add_subplot(gs7[1, c])
    ax_div   = fig.add_subplot(gs7[2, c])
    ax_tumor = fig.add_subplot(gs7[3, c])

    ax_bar.imshow(np.zeros((1,n_rois)), aspect="auto",
                  cmap=ListedColormap([tc]), vmin=0, vmax=1, interpolation="nearest")
    ax_bar.set_xticks([]); ax_bar.set_yticks([])
    ax_bar.set_title(f"{tissue}  (n={n_rois})", fontsize=16,
                     fontweight="bold", color=tc, pad=5)

    mat_mac = df_z[mac_rows7].T.values.astype(float)
    im = ax_mac.imshow(mat_mac, aspect="auto", cmap=CMAP7, norm=norm7,
                       interpolation="nearest")
    ax_mac.set_yticks(range(3))
    ax_mac.set_yticklabels(row_labels_mac if c==0 else [""]*3, fontsize=12, rotation=90, va='center')
    if c==0: ax_mac.set_ylabel("Macrophage\nCompartment (CD68+)",
                               fontsize=11, fontweight="bold", labelpad=4)
    ax_mac.set_xticks([])
    for ri, feat in enumerate(mac_rows7):
        r_v, p_v = pearsonr(df_z["CD20_Ligand"].values, df_z[feat].values)
        p_str = f"p={p_v:.1e}" if p_v >= 1e-4 else "p<1e-4"
        ax_mac.text(n_rois*0.98, ri, f"r={r_v:+.2f} {p_str}",
                    fontsize=13, va="center", ha="right", color="white",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc="#00000066", ec="none"))

    ax_div.set_facecolor("#1A1A1A"); ax_div.set_xticks([]); ax_div.set_yticks([])
    if c==0:
        ax_div.text(0.5, 0.5, "── Compartment Boundary ──",
                    transform=ax_div.transAxes, ha="center", va="center",
                    fontsize=11, color="white", fontweight="bold")

    mat_tum = df_z[["CD20_Ligand"]].T.values.astype(float)
    ax_tumor.imshow(mat_tum, aspect="auto", cmap=CMAP7, norm=norm7,
                    interpolation="nearest")
    ax_tumor.set_yticks([0])
    ax_tumor.set_yticklabels(row_labels_tum if c==0 else [""], fontsize=11, rotation=90, va='center')
    if c==0: ax_tumor.set_ylabel("Tumor B-cell\nCompartment (CD20+)",
                                 fontsize=11, fontweight="bold", labelpad=4)
    ax_tumor.set_xticks([])
    ax_tumor.set_xlabel("ROIs sorted by HMGB1 ↑  (Low → High)",
                        fontsize=12, labelpad=5, color=tc)
    r_tgs, p_tgs = pearsonr(df_z["CD20_Ligand"].values,
                             df_z["CD68_Target_Gene_Score"].values)
    ax_tumor.text(0.5, -0.42,
                  f"HMGB1 ↔ Target Gene Score:  r={r_tgs:+.3f},  p={p_tgs:.2e}",
                  transform=ax_tumor.transAxes, ha="center", fontsize=13,
                  color=tc, fontweight="bold",
                  bbox=dict(boxstyle="round,pad=0.25", fc="white",
                            ec=tc, alpha=0.9, linewidth=1.0))

ax_cbar7 = fig.add_subplot(gs7[1:4, 3])
sm7 = cm.ScalarMappable(cmap=CMAP7, norm=norm7); sm7.set_array([])
plt.colorbar(sm7, cax=ax_cbar7, label="Z-score")
ax_cbar7.tick_params(labelsize=8)
fig.suptitle("Spatial Co-evolution Heatmap: DLBCL vs Normal (Side-by-Side)"
             "HMGB1 (Ligand) → HAVCR2 (Receptor) → MAFB (TF) → Target Gene Score (23 genes) "
             "GSE232853 · 228 Paired ROIs · Z-scored on pooled data · Sorted by HMGB1 within group",
             fontsize=10, fontweight="bold", y=1.03)
savefig(fig, "fig7_spatial_coevolution_split_heatmap")

# ═══════════════════════════════════════════════════════════════
# FIG 8 — 2D Density Contour: HMGB1 vs Target Gene Score
# ═══════════════════════════════════════════════════════════════
print("Fig8: 2D Density contour...")

master_r = master.rename(columns={"CD20_Ligand":"HMGB1",
                                   "CD68_Target_Gene_Score":"Target Gene Score"})
dlbcl8  = master_r[master_r["Tissue_Type"]=="DLBCL"]
normal8 = master_r[master_r["Tissue_Type"]=="Normal"]

fig, ax = plt.subplots(figsize=(9, 7))
fig.patch.set_facecolor("#FFFFFF"); ax.set_facecolor("#FFFFFF")
x_all8  = master_r["HMGB1"].values
x_line8 = np.linspace(x_all8.min()-0.1, x_all8.max()+0.1, 300)

for tissue, df, color in [("Normal", normal8, TISSUE_COLORS["Normal"]),
                           ("DLBCL",  dlbcl8,  TISSUE_COLORS["DLBCL"])]:
    x = df["HMGB1"].values; y = df["Target Gene Score"].values
    ax.scatter(x, y, c=color, s=18, alpha=0.35, linewidths=0,
               zorder=2, rasterized=True)
    xy  = np.vstack([x, y])
    kde = gaussian_kde(xy, bw_method=0.35)
    xi  = np.linspace(x.min()-0.2, x.max()+0.2, 120)
    yi  = np.linspace(y.min()-0.02, y.max()+0.02, 120)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = kde(np.vstack([Xi.ravel(), Yi.ravel()])).reshape(Xi.shape)
    ax.contourf(Xi, Yi, Zi, levels=6, colors=[color]*6, alpha=0.08, zorder=3)
    ax.contour(Xi, Yi, Zi, levels=5, colors=[color],
               linewidths=[0.6,0.9,1.2,1.5,1.8], alpha=0.85, zorder=4)
    y_line8, ci_band8 = regression_ci(x, y, x_line8)
    ax.plot(x_line8, y_line8, color=color, lw=2.2, zorder=5)
    ax.fill_between(x_line8, y_line8-ci_band8, y_line8+ci_band8,
                    color=color, alpha=0.12, zorder=3)
    r, pr   = pearsonr(x, y)
    rho, ps = spearmanr(x, y)
    y_pos = 0.97 if tissue=="Normal" else 0.82
    ax.text(0.03, y_pos,
            f"{tissue}  (n={len(x)}) r = {r:+.3f},  p = {pr:.2e} ρ = {rho:+.3f},  p = {ps:.2e}",
            transform=ax.transAxes, fontsize=9, va="top", color=color,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec=color, alpha=0.90, linewidth=1.0))

ax.set_xlabel("B-cell HMGB1 expression (Q3 log-norm)", fontsize=11)
ax.set_ylabel("Macrophage Target Gene Score (mean of 23 genes, Q3 log-norm)", fontsize=11)
ax.set_title("2D Density Contour: HMGB1 vs Target Gene Score"
             "GSE232853 · 228 Paired CD20+/CD68+ ROIs · DLBCL vs Normal",
             fontsize=11, fontweight="bold", pad=10)
handles8 = [mpatches.Patch(color=TISSUE_COLORS[t], alpha=0.75, label=t)
            for t in ["DLBCL","Normal"]]
ax.legend(handles=handles8, fontsize=10, framealpha=0.92,
          edgecolor="#CCCCCC", loc="lower right")
ax.tick_params(labelsize=9)
for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")
plt.tight_layout()
savefig(fig, "fig8_2D_density_HMGB1_TGS")

# ═══════════════════════════════════════════════════════════════
# FIG 9 — Violin: HMGB1 & TGS, Normal vs DLBCL
# ═══════════════════════════════════════════════════════════════
print("Fig9: Violin HMGB1 & TGS Normal vs DLBCL...")

master9 = master.rename(columns={"CD20_Ligand":"HMGB1",
                                  "CD68_Target_Gene_Score":"Target Gene Score"})
FEATURES9 = ["HMGB1", "Target Gene Score"]
YLABELS9  = {"HMGB1": "HMGB1 expression (CD20+ ROI, Q3 log-norm)",
             "Target Gene Score": "Target Gene Score (CD68+ ROI, mean 23 genes, Q3 log-norm)"}
ORDER9 = ["Normal", "DLBCL"]

fig, axes = plt.subplots(1, 2, figsize=(10, 6))
fig.patch.set_facecolor("#FFFFFF")

for ax, feat in zip(axes, FEATURES9):
    ax.set_facecolor("#FFFFFF")
    groups9 = [master9.loc[master9["Tissue_Type"]==t, feat].dropna().values
               for t in ORDER9]
    parts = ax.violinplot(groups9, positions=[1,2], widths=0.55,
                          showmedians=False, showextrema=False)
    for pc, tissue in zip(parts["bodies"], ORDER9):
        pc.set_facecolor(TISSUE_COLORS[tissue]); pc.set_edgecolor("white"); pc.set_alpha(0.65)
    rng9 = np.random.default_rng(42)
    for i, (grp, tissue) in enumerate(zip(groups9, ORDER9), start=1):
        q1, med, q3 = np.percentile(grp, [25,50,75])
        ax.plot([i,i], [q1,q3], color="black", lw=2.5, zorder=4)
        ax.scatter(i, med, color="white", s=55, zorder=5,
                   edgecolors="black", linewidths=1.2)
        jitter = rng9.uniform(-0.12, 0.12, size=len(grp))
        ax.scatter(i+jitter, grp, color=TISSUE_COLORS[tissue],
                   s=12, alpha=0.45, linewidths=0, zorder=3)
    stat9, p9 = mannwhitneyu(groups9[0], groups9[1], alternative="two-sided")
    y_max9 = max(np.max(g) for g in groups9)
    y_top9 = y_max9 * 1.05
    ax.plot([1,1,2,2], [y_top9, y_top9*1.02, y_top9*1.02, y_top9], color="black", lw=1.2)
    ax.text(1.5, y_top9*1.025, f"{sig_stars(p9)} p={p9:.2e}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks([1,2]); ax.set_xticklabels(ORDER9, fontsize=11)
    ax.set_ylabel(YLABELS9[feat], fontsize=10)
    ax.set_title(feat, fontsize=12, fontweight="bold", pad=8)
    ax.tick_params(axis="y", labelsize=9)
    for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")
    ax.set_xlim(0.4, 2.6)
    for i, (grp, tissue) in enumerate(zip(groups9, ORDER9), start=1):
        ax.text(i, ax.get_ylim()[0] - (ax.get_ylim()[1]-ax.get_ylim()[0])*0.06,
                f"n={len(grp)}", ha="center", va="top", fontsize=8.5,
                color=TISSUE_COLORS[tissue], fontweight="bold")

handles9 = [mpatches.Patch(color=TISSUE_COLORS[t], alpha=0.75, label=t) for t in ORDER9]
fig.legend(handles=handles9, loc="upper center", ncol=2, fontsize=10,
           framealpha=0.9, edgecolor="#CCCCCC", bbox_to_anchor=(0.5, 1.01))
fig.suptitle("HMGB1 & Target Gene Score: Normal vs DLBCL"
             "GSE232853 · 228 Paired CD20+/CD68+ ROIs",
             fontsize=11, fontweight="bold", y=1.06)
plt.tight_layout()
savefig(fig, "fig9_violin_HMGB1_TGS_Normal_vs_DLBCL")


# ═══════════════════════════════════════════════════════════════
# S5A — GeoMx ROI / AOI Design Schematic
# ═══════════════════════════════════════════════════════════════
print("S5A: GeoMx ROI/AOI composition and analysis schema...")

# ── 从 meta 提取数据组成信息 ──────────────────────────────────────────
n_total    = len(meta)
mask_col   = [c for c in meta.columns if "mask"   in c.lower()][0]
tissue_col = [c for c in meta.columns if "tissue" in c.lower()][0]
n_cd20     = (meta[mask_col]   == "CD20").sum()
n_cd68     = (meta[mask_col]   == "CD68").sum()
n_dlbcl    = (meta[tissue_col] == "DLBCL").sum()
n_normal   = (meta[tissue_col] == "Normal").sum()
n_paired   = len(paired)

# ── canvas ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(18, 11))
fig.patch.set_facecolor("#FFFFFF")
ax.set_facecolor("#FFFFFF")
W, H = 18, 11
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

# ── 颜色 ──────────────────────────────────────────────────────────────
C_DLBCL    = TISSUE_COLORS.get("DLBCL",  "#C0392B")
C_NORMAL   = TISSUE_COLORS.get("Normal", "#2980B9")
C_CD20     = MASK_COLORS.get("CD20",    "#2471A3")
C_CD68     = MASK_COLORS.get("CD68",    "#CB4335")
C_PAIRED   = "#7B2D8B"
C_ANALYSIS = ["#B2182B", "#D6604D", "#4393C3", "#2166AC"]

# ── helpers ───────────────────────────────────────────────────────────
def rbox(x, y, w, h, ec, fc, lw=1.8, ls="-", zorder=2, alpha=1.0):
    ax.add_patch(plt.Rectangle(
        (x, y), w, h, linewidth=lw,
        edgecolor=ec, facecolor=fc,
        linestyle=ls, zorder=zorder, alpha=alpha))

def arr(x0, y0, x1, y1, color="#666666", lw=1.4, rad=0.0):
    """x0,y0 = 起点（尾部）  x1,y1 = 终点（箭头头部）"""
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(
                    arrowstyle="->", color=color, lw=lw,
                    connectionstyle=f"arc3,rad={rad}"))

def label(x, y, txt, fs=9, fw="normal", color="#222222",
          ha="center", va="center", ls="normal", zorder=3):
    ax.text(x, y, txt, ha=ha, va=va, fontsize=fs,
            fontweight=fw, color=color, style=ls,
            linespacing=1.45, zorder=zorder)

# ══════════════════════════════════════════════════════════════════════
# TITLE
# ══════════════════════════════════════════════════════════════════════
label(W/2, 10.70,
      "ROI/AOI Composition and Analysis Framework",
      fs=16, fw="bold", color="#111111")

# ══════════════════════════════════════════════════════════════════════
# ROW 1 — 总 cohort 框（居中）
# ══════════════════════════════════════════════════════════════════════
R1W, R1H = 4.0, 1.6
R1X = W/2 - R1W/2
R1Y = 8.45
rbox(R1X, R1Y, R1W, R1H, "#444444", "#F2F2F2", lw=2.0)
cx1 = W/2

label(cx1, R1Y+R1H-0.30, "GeoMx DSP Cohort  (GSE232853)",
      fs=14, fw="bold", color="#222222")
label(cx1, R1Y+R1H-0.72, f"Total AOIs  n = {n_total}",
      fs=13, fw="bold", color="#333333")
label(cx1, R1Y+0.44,
      "Q3 normalization\nHarmony batch correction\n18,677 genes (WTA)",
      fs=11, color="#777777", ls="italic")

# ══════════════════════════════════════════════════════════════════════
# ROW 2 — 三个分组框
#   左:  Tissue Type  (DLBCL / Normal)
#   中:  Cell Compartment  (CD68 / CD20)
#   右:  Paired ROIs
# ══════════════════════════════════════════════════════════════════════
R2H = 1.90
R2Y = 6.10

# —— 左框：Tissue Type ————————————————————————————————————————————————
LX, LW = 0.5, 5.0
rbox(LX, R2Y, LW, R2H, "#555555", "#FAFAFA", lw=1.5)
label(LX+LW/2, R2Y+R2H-0.28, "Tissue Type",
      fs=14, fw="bold", color="#333333")

rbox(LX+0.20, R2Y+0.50, 2.1, 0.85, C_DLBCL, "#FDECEA", lw=1.4)
label(LX+0.20+1.05, R2Y+0.50+0.55,
      "DLBCL", fs=13, fw="bold", color=C_DLBCL)
label(LX+0.20+1.05, R2Y+0.50+0.20,
      f"n = {n_dlbcl}\n64 patients · 87 biopsies",
      fs=10, color=C_DLBCL)

rbox(LX+2.70, R2Y+0.50, 2.1, 0.85, C_NORMAL, "#EBF5FB", lw=1.4)
label(LX+2.70+1.05, R2Y+0.50+0.55,
      "Normal (RLT)", fs=13, fw="bold", color=C_NORMAL)
label(LX+2.70+1.05, R2Y+0.50+0.20,
      f"n = {n_normal}\n24 tonsil samples",
      fs=10, color=C_NORMAL)

# —— 中框：Cell Compartment ───────────────────────────────────────────
MX, MW = 6.50, 5.0
rbox(MX, R2Y, MW, R2H, "#555555", "#FAFAFA", lw=1.5)
label(MX+MW/2, R2Y+R2H-0.28, "Cell Compartment",
      fs=14, fw="bold", color="#333333")

rbox(MX+0.20, R2Y+0.50, 2.1, 0.85, C_CD68, "#FDECEA", lw=1.4)
label(MX+0.20+1.05, R2Y+0.50+0.55,
      "CD68+ AOI", fs=13, fw="bold", color=C_CD68)
label(MX+0.20+1.05, R2Y+0.50+0.20,
      f"n = {n_cd68}\nMacrophage compartment",
      fs=10, color=C_CD68)

rbox(MX+2.70, R2Y+0.50, 2.1, 0.85, C_CD20, "#EBF5FB", lw=1.4)
label(MX+2.70+1.05, R2Y+0.50+0.55,
      "CD20+ AOI", fs=13, fw="bold", color=C_CD20)
label(MX+2.70+1.05, R2Y+0.50+0.20,
      f"n = {n_cd20}\nB-cell compartment",
      fs=10, color=C_CD20)

# —— 右框：Paired ROIs ────────────────────────────────────────────────
RX, RW = 12.50, 5.0
rbox(RX, R2Y, RW, R2H, C_PAIRED, "#F5EEF8", lw=2.0, ls="--")
label(RX+RW/2, R2Y+R2H-0.28, "Paired ROIs",
      fs=15, fw="bold", color=C_PAIRED)
label(RX+RW/2, R2Y+R2H-0.68, f"n = {n_paired}",
      fs=14, fw="bold", color=C_PAIRED)
label(RX+RW/2, R2Y+0.88,
      "Same tissue section with both",
      fs=11, color="#555555")
label(RX+RW/2, R2Y+0.56,
      "CD20+ and CD68+ AOIs matched by Unique_ROI",
      fs=11, color="#555555")
label(RX+RW/2, R2Y+0.22,
      "Used for cross-compartment correlation analyses",
      fs=11, color="#888888", ls="italic")

# ── Arrows: Row1 → Row2（全部直线，rad=0.0）──────────────────────────
arr(cx1, R1Y,        LX+LW/2,  R2Y+R2H, "#888888", lw=1.3, rad=0.0)
arr(cx1, R1Y,        MX+MW/2,  R2Y+R2H, "#888888", lw=1.3, rad=0.0)
arr(cx1, R1Y,        RX+RW/2,  R2Y+R2H, C_PAIRED,  lw=1.5, rad=0.0)

# ══════════════════════════════════════════════════════════════════════
# ROW 3 — 分析模块（4 个等宽框）
# ══════════════════════════════════════════════════════════════════════
R3W  = 3.60
R3H  = 2.30
R3Y  = 3.30
GAP3 = (W - 4*R3W) / 5
xs3  = [GAP3 + i*(R3W + GAP3) for i in range(4)]
cxs3 = [x + R3W/2 for x in xs3]

modules = [
    (C_ANALYSIS[0],
     "UMAP Overview",
     "Cohort-level distribution",
     "Colored by compartment (CD20/CD68)\n"
     "and tissue type (DLBCL/Normal)",
     f"All {n_total} AOIs · post-Harmony embedding"),

    (C_ANALYSIS[1],
     "Differential Expression\n& Marker Comparison",
     "DLBCL vs Normal · CD68+ macrophages",
     "Scanpy Wilcoxon · BH FDR\n"
     "|LFC| ≥ 0.5 · padj ≤ 0.05\n"
     "Highlights: HMGB1 · HAVCR2 · MAFB",
     f"CD68+ AOIs  n = {n_cd68}"),

    (C_ANALYSIS[2],
     "Cross-Compartment\nPaired Correlations",
     "B-cell HMGB1 vs macrophage response",
     "Spearman r · OLS regression\n"
     "B-cell HMGB1 → macrophage HAVCR2\n"
     "B-cell HMGB1 → macrophage LA_TAM score",
     f"Paired ROIs  n = {n_paired}"),

    (C_ANALYSIS[3],
     "Spatial Co-evolution\nHeatmap & GSEA",
     "LA_TAM-rich vs LA_TAM-poor niches",
     "ssGSEA · Hallmark gene sets\n"
     "HMGB1 → HAVCR2 → MAFB → Target score\n"
     "Sorted by B-cell HMGB1 expression",
     f"Paired ROIs  n = {n_paired}"),
]

for i, (color, title, subtitle, methods, source) in enumerate(modules):
    x0 = xs3[i]
    cx = cxs3[i]

    # 框体
    rbox(x0, R3Y, R3W, R3H, color, "white", lw=1.5)

    # 顶部色条
    rbox(x0, R3Y+R3H-0.40, R3W, 0.40, color, color, lw=0, zorder=3)
    label(cx, R3Y+R3H-0.20, title,
          fs=12, fw="bold", color="white", zorder=4)

    # 副标题
    label(cx, R3Y+R3H-0.68, subtitle,
          fs=11, color=color, fw="bold")

    # 分析方法
    label(cx, R3Y+1, methods,
          fs=10, color="#444444")

    # 数据来源标签
    rbox(x0+0.18, R3Y+0.10, R3W-0.36, 0.36,
         color, "#F8F8F8", lw=1.0, zorder=3)
    label(cx, R3Y+0.28, source,
          fs=10, color=color, fw="bold", zorder=4)

# ── Arrows: Row2 → Row3（全部直线）───────────────────────────────────
arr(LX+LW/2,  R2Y, cxs3[0], R3Y+R3H, C_ANALYSIS[0], lw=1.2, rad=0.0)
arr(MX+MW/2,  R2Y, cxs3[1], R3Y+R3H, C_ANALYSIS[1], lw=1.2, rad=0.0)
arr(RX+RW/2,  R2Y, cxs3[2], R3Y+R3H, C_ANALYSIS[2], lw=1.2, rad=0.0)
arr(RX+RW/2,  R2Y, cxs3[3], R3Y+R3H, C_ANALYSIS[3], lw=1.2, rad=0.15)



plt.tight_layout(rect=[0, 0.08, 1, 1])
savefig(fig, "S5A_ROI_AOI_design_schematic")



# ═══════════════════════════════════════════════════════════════
# S5B — PCA / UMAP colored by Sample / Patient (batch check)
# ═══════════════════════════════════════════════════════════════
print("S5B: PCA/UMAP by sample (batch check)...")

# Detect PCA columns
pc1_col = [c for c in pca_df.columns if "pc1" in c.lower() or
           ("pc" in c.lower() and "1" in c)][0]
pc2_col = [c for c in pca_df.columns if "pc2" in c.lower() or
           ("pc" in c.lower() and "2" in c)][0]

# Detect slide/sample ID column (try several common names)
slide_candidates = [c for c in pca_df.columns
                    if any(k in c.lower() for k in
                           ["slide","patient","sample","batch","subject","case"])]
slide_col_pca = slide_candidates[0] if slide_candidates else pca_df.columns[0]

# Detect mask & tissue in pca_df
mask_pca   = [c for c in pca_df.columns if "mask" in c.lower()][0]   if any("mask" in c.lower() for c in pca_df.columns) else None
tissue_pca = [c for c in pca_df.columns if "tissue" in c.lower()][0] if any("tissue" in c.lower() for c in pca_df.columns) else None

# UMAP: use same slide_col if present, else merge from meta
umap_slide_candidates = [c for c in umap_df.columns
                          if any(k in c.lower() for k in
                                 ["slide","patient","sample","batch","subject","case"])]
slide_col_umap = umap_slide_candidates[0] if umap_slide_candidates else None

# Build a qualitative colormap for slides
def make_slide_colors(series):
    cats = sorted(series.dropna().unique())
    base_colors = plt.cm.tab20.colors + plt.cm.tab20b.colors
    return {cat: base_colors[i % len(base_colors)] for i, cat in enumerate(cats)}

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.patch.set_facecolor("#FFFFFF")

# ── Left: PCA colored by Slide/Sample ──
ax = axes[0]
ax.set_facecolor("#FFFFFF")
slide_vals_pca = pca_df[slide_col_pca].astype(str)
slide_cmap_pca = make_slide_colors(slide_vals_pca)
ax.scatter(pca_df[pc1_col], pca_df[pc2_col],
           c="#DDDDDD", s=10, alpha=0.25, linewidths=0, zorder=1, rasterized=True)
for slide in sorted(slide_vals_pca.unique()):
    idx = slide_vals_pca == slide
    ax.scatter(pca_df.loc[idx, pc1_col], pca_df.loc[idx, pc2_col],
               c=[slide_cmap_pca[slide]], s=18, alpha=0.75,
               linewidths=0, zorder=2, rasterized=True, label=slide)
ax.set_xlabel("PC1", fontsize=11); ax.set_ylabel("PC2", fontsize=11)
ax.set_title(f"PCA — colored by {slide_col_pca} (batch / sample clustering check)",
             fontsize=11, fontweight="bold", pad=8)
# Legend: only if ≤ 20 categories, else note
n_slides = slide_vals_pca.nunique()
if n_slides <= 20:
    ax.legend(fontsize=7, framealpha=0.85, edgecolor="#CCCCCC",
              loc="best", markerscale=1.2,
              ncol=2 if n_slides > 10 else 1,
              title=slide_col_pca, title_fontsize=7)
else:
    ax.text(0.97, 0.03, f"572 unique values (legend omitted for clarity)",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=12, color="#666666", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="#CCCCCC", alpha=0.85))
ax.tick_params(labelsize=9)
for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")

# ── Right: UMAP colored by Slide/Sample ──
ax2 = axes[1]
ax2.set_facecolor("#FFFFFF")
if slide_col_umap:
    slide_vals_umap = umap_df[slide_col_umap].astype(str)
    slide_cmap_umap = make_slide_colors(slide_vals_umap)
    ax2.scatter(umap_df[u1_col], umap_df[u2_col],
                c="#DDDDDD", s=10, alpha=0.25, linewidths=0, zorder=1, rasterized=True)
    for slide in sorted(slide_vals_umap.unique()):
        idx = slide_vals_umap == slide
        ax2.scatter(umap_df.loc[idx, u1_col], umap_df.loc[idx, u2_col],
                    c=[slide_cmap_umap[slide]], s=18, alpha=0.75,
                    linewidths=0, zorder=2, rasterized=True, label=slide)
    n_slides_u = slide_vals_umap.nunique()
    col_label  = slide_col_umap
    if n_slides_u <= 20:
        ax2.legend(fontsize=7, framealpha=0.85, edgecolor="#CCCCCC",
                   loc="best", markerscale=1.2,
                   ncol=2 if n_slides_u > 10 else 1,
                   title=col_label, title_fontsize=7)
    else:
        ax2.text(0.97, 0.03, f"{n_slides_u} unique values (legend omitted for clarity)",
                 transform=ax2.transAxes, ha="right", va="bottom",
                 fontsize=12, color="#666666", style="italic",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white",
                           ec="#CCCCCC", alpha=0.85))
else:
    # Fallback: color by Mask if slide column not in UMAP file
    mask_vals_u = umap_df[[c for c in umap_df.columns
                            if "mask" in c.lower()][0]]
    ax2.scatter(umap_df[u1_col], umap_df[u2_col],
                c="#DDDDDD", s=10, alpha=0.25, linewidths=0, zorder=1, rasterized=True)
    for grp in sorted(mask_vals_u.unique()):
        idx = mask_vals_u == grp
        ax2.scatter(umap_df.loc[idx, u1_col], umap_df.loc[idx, u2_col],
                    c=[MASK_COLORS.get(grp,"#888888")], s=18, alpha=0.75,
                    linewidths=0, zorder=2, rasterized=True, label=grp)
    ax2.legend(fontsize=9, framealpha=0.9, edgecolor="#CCCCCC", loc="best")
    col_label = "Mask"

ax2.set_xlabel("UMAP 1", fontsize=11); ax2.set_ylabel("UMAP 2", fontsize=11)
ax2.set_title(f"UMAP — colored by {col_label} (batch / sample clustering check)",
              fontsize=11, fontweight="bold", pad=8)
ax2.tick_params(labelsize=9)
for sp in ax2.spines.values(): sp.set_edgecolor("#CCCCCC")

fig.suptitle("S5B — PCA & UMAP: Sample / Batch Distribution Check"
             "GSE232853 · GeoMx DSP · No dominant sample-level clustering expected after Harmony",
             fontsize=11, fontweight="bold", y=1.02)
plt.tight_layout()
savefig(fig, "S5B_PCA_UMAP_by_sample")


# ============================================================
#  Volcano Plot — 本地 Windows 版（修正版）
#  pip install pandas matplotlib adjusttext numpy
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches
from adjustText import adjust_text

# ── 路径配置 ────────────────────────────────────────────────
DATA_DIR = translate(r"D:\bulk-download\GSE232853_v2")
CSV_PATH = os.path.join(DATA_DIR, "fig3_DEA_CD68_scanpy.csv")
OUT_PATH = os.path.join(DATA_DIR, "Fig3_volcano_CD68.png")

# ── 参数 ────────────────────────────────────────────────────
LFC_THR  = 0.5
PADJ_THR = 0.05
TOP_N    = 15
FORCE    = ["HMGB1", "HAVCR2", "MAFB"]

COL = {"Up":   "#C0392B",
       "Down": "#2980B9",
       "Sig":  "#F39C12",
       "NS":   "#CCCCCC"}

# ── 读取数据 ────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH, sep=",")
df.columns = df.columns.str.strip()

# ── 分类 ────────────────────────────────────────────────────
def classify(row):
    lfc, p = row["logfoldchanges"], row["pvals_adj"]
    if p <= PADJ_THR and lfc >= LFC_THR:
        return "Up"
    elif p <= PADJ_THR and lfc <= -LFC_THR:
        return "Down"
    elif p <= PADJ_THR:
        return "Sig"
    else:
        return "NS"

df["Category"]    = df.apply(classify, axis=1)
df["neg_log10_p"] = df["pvals_adj"].apply(
    lambda x: -np.log10(x) if x > 0 else 300
)

# ── 选取标注基因 ─────────────────────────────────────────────
up_genes   = (df[df["Category"] == "Up"]
              .nlargest(TOP_N, "neg_log10_p")["gene"].tolist())
down_genes = (df[df["Category"] == "Down"]
              .nlargest(TOP_N, "neg_log10_p")["gene"].tolist())
label_set  = set(up_genes + down_genes + FORCE)

# ── 绘图 ────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "font.family":  "Arial",
    "font.size":    10,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(7, 6))

for cat in ["NS", "Sig", "Down", "Up"]:
    sub  = df[df["Category"] == cat]
    size = 8 if cat == "NS" else 14
    alpha = 0.6 if cat == "NS" else 0.85
    ax.scatter(sub["logfoldchanges"], sub["neg_log10_p"],
               c=COL[cat], s=size, alpha=alpha,
               linewidths=0, rasterized=True, label=cat, zorder=2)

ax.axhline(-np.log10(PADJ_THR), color="#888888", lw=0.8, ls="--", zorder=1)
ax.axvline( LFC_THR,            color="#888888", lw=0.8, ls="--", zorder=1)
ax.axvline(-LFC_THR,            color="#888888", lw=0.8, ls="--", zorder=1)

# ── 标注 ────────────────────────────────────────────────────
texts = []
for _, row in df[df["gene"].isin(label_set)].iterrows():
    gene     = row["gene"]
    x, y     = row["logfoldchanges"], row["neg_log10_p"]
    is_force = gene in FORCE

    if is_force:
        ax.scatter(x, y, s=60, facecolors="none",
                   edgecolors="black", linewidths=1.2, zorder=5)

    label = f"▲{gene}" if is_force else gene
    t = ax.text(x, y, label,
                fontsize=10 if is_force else 9,
                fontweight="bold" if is_force else "normal",
                color="black", zorder=6)
    texts.append(t)

adjust_text(texts, ax=ax,
            arrowprops=dict(
                arrowstyle="-",
                color="#AAAAAA",
                lw=0.5,
                shrinkA=8,      # 箭头起点缩进，避免穿过标签
                shrinkB=3,      # 箭头终点缩进，避免穿过散点
            ),
            expand_points=(1.3, 1.5),
            expand_text=(1.2, 1.4),
            force_text=(0.4, 0.5))


# ── 图例与轴标签 ─────────────────────────────────────────────
legend_labels = {
    "Up":   f"Up (LFC≥{LFC_THR}, padj≤{PADJ_THR})",
    "Down": f"Down (LFC≤-{LFC_THR}, padj≤{PADJ_THR})",
    "Sig":  "Sig (padj≤0.05, |LFC|<0.5)",
    "NS":   "NS",
}
handles = [matplotlib.patches.Patch(facecolor=COL[k], label=v)
           for k, v in legend_labels.items()]
ax.legend(handles=handles, frameon=False, fontsize=10,
          loc="lower left", handlelength=1)

n_up   = (df["Category"] == "Up").sum()
n_down = (df["Category"] == "Down").sum()
ax.text(0.98, 0.02,
        f"Up: {n_up}  |  Down: {n_down}",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=11, color="#444444")

ax.set_xlabel("log$_2$ Fold Change (DLBCL vs Normal)", fontsize=11)
ax.set_ylabel("-log$_{10}$(adjusted p-value)",          fontsize=11)
ax.set_title("CD68$^+$ Macrophages — DLBCL vs Normal",
             fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
print(f"图已保存至：{OUT_PATH}")
plt.show()




import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import warnings
import os

warnings.filterwarnings("ignore")

# ══════════════════════════════════════════
# 路径配置
# ══════════════════════════════════════════
BASE_DIR   = translate(r"D:\bulk-download\GSE232853_v2")
EXPR_CSV   = os.path.join(BASE_DIR, "expression_lognorm_all_samples.csv")
PAIRED_CSV = os.path.join(BASE_DIR, "task5_paired_ROI_LA_TAM.csv")
DEA_CSV    = os.path.join(BASE_DIR, "task4_DEA_CD68_DLBCL_vs_Normal.csv")
OUTPUT_DIR = BASE_DIR

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ══════════════════════════════════════════
# 1. 载入表达矩阵
#    index 格式：'DLBCL A3 | 001 | CD20'
# ══════════════════════════════════════════
print("Loading expression matrix...")
expr_df = pd.read_csv(EXPR_CSV, index_col=0)
print(f"  expr_df shape: {expr_df.shape}  (samples × genes)")
print(f"  Index sample: {expr_df.index[0]!r}")

# ══════════════════════════════════════════
# 2. 解析 index → Slide_ID / ROI_Num / Mask
#    格式：'DLBCL A3 | 001 | CD20'
#    分隔符：' | '（注意两侧有空格）
# ══════════════════════════════════════════
def parse_sample_id(sample_id):
    """
    解析 'DLBCL A3 | 001 | CD20' 格式的 Sample_ID。
    返回 (slide_id, roi_num_str, mask)
    roi_num_str 保留前导零以便匹配，同时也转为 int 供 paired 匹配。
    """
    parts = [p.strip() for p in sample_id.split("|")]
    if len(parts) != 3:
        return None, None, None
    slide_id = parts[0].strip()
    roi_str  = parts[1].strip()      # '001'
    mask     = parts[2].strip()      # 'CD20' or 'CD68'
    return slide_id, roi_str, mask

# 构建 meta DataFrame，直接从 expr_df.index 解析
meta_records = []
for sid in expr_df.index:
    slide_id, roi_str, mask = parse_sample_id(sid)
    if slide_id is None:
        print(f"  [WARN] Cannot parse sample ID: {sid!r}, skipping")
        continue
    try:
        roi_num = int(roi_str)
    except ValueError:
        print(f"  [WARN] Cannot convert ROI_Num to int: {roi_str!r}, skipping")
        continue
    meta_records.append({
        "Sample_ID": sid,
        "Slide_ID":  slide_id,
        "ROI_Num":   roi_num,
        "Mask":      mask
    })

meta = pd.DataFrame(meta_records).set_index("Sample_ID")
print(f"  Parsed meta: {meta.shape}")
print(f"  Mask values: {meta['Mask'].unique().tolist()}")

# ══════════════════════════════════════════
# 3. 载入 paired 和 DEA
# ══════════════════════════════════════════
paired = pd.read_csv(PAIRED_CSV)
dea    = pd.read_csv(DEA_CSV)
print(f"  paired: {paired.shape},  dea: {dea.shape}")

# ══════════════════════════════════════════
# 4. 检测 DEA 列名（兼容多种命名）
# ══════════════════════════════════════════
logfc_col = next(
    (c for c in ["LogFC_DLBCL_vs_Normal", "logFC", "log2FoldChange", "LogFC"]
     if c in dea.columns),
    None
)
if logfc_col is None:
    raise ValueError(f"DEA 文件找不到 LogFC 列，现有列：{dea.columns.tolist()}")

gene_col = "Gene" if "Gene" in dea.columns else dea.columns[0]

# ══════════════════════════════════════════
# 5. 构建 HEATMAP_GENES
# ══════════════════════════════════════════
top_up = (dea[dea[logfc_col] >  0.5].nlargest(12, logfc_col)[gene_col].tolist())
top_dn = (dea[dea[logfc_col] < -0.5].nsmallest(8,  logfc_col)[gene_col].tolist())
latam  = ["HAVCR2", "MAFB", "APOE", "LGALS9C", "C1QA", "C1QB", "FCER1G", "HMGB1"]

HEATMAP_GENES = list(dict.fromkeys(latam + top_up + top_dn))
HEATMAP_GENES = [g for g in HEATMAP_GENES if g in expr_df.columns]
print(f"  Heatmap genes ({len(HEATMAP_GENES)}): {HEATMAP_GENES}")

if len(HEATMAP_GENES) == 0:
    print(f"  [DEBUG] expr_df columns (first 20): {expr_df.columns[:20].tolist()}")
    print(f"  [DEBUG] latam not found: {[g for g in latam if g not in expr_df.columns]}")
    raise ValueError(
        "HEATMAP_GENES 为空。\n"
        "请检查 [DEBUG] 输出，确认基因名格式是否一致（大小写、空格等）"
    )

# ══════════════════════════════════════════
# 6. 构建 _key 用于配对匹配
#    paired 里 ROI_Num 是整数，meta 里也转成了整数
# ══════════════════════════════════════════
meta["_key"]   = meta["Slide_ID"].astype(str) + "_" + meta["ROI_Num"].astype(str)
paired["_key"] = paired["Slide_ID"].astype(str) + "_" + paired["ROI_Num"].astype(str)

paired_sorted = paired.sort_values(
    ["LA_TAM_Group", "Tissue_Type", "LA_TAM_Score_CD68"],
    ascending=[False, True, False]
).reset_index(drop=True)

# ══════════════════════════════════════════
# 7. 构建配对矩阵
# ══════════════════════════════════════════
mat_cd20, mat_cd68, valid_rows = [], [], []
skipped_detail = []

for _, row in paired_sorted.iterrows():
    k = row["_key"]

    obs_cd20 = meta[(meta["_key"] == k) & (meta["Mask"] == "CD20")].index
    obs_cd68 = meta[(meta["_key"] == k) & (meta["Mask"] == "CD68")].index

    if len(obs_cd20) == 1 and len(obs_cd68) == 1:
        mat_cd20.append(expr_df.loc[obs_cd20[0], HEATMAP_GENES].values)
        mat_cd68.append(expr_df.loc[obs_cd68[0], HEATMAP_GENES].values)
        valid_rows.append(row)
    else:
        skipped_detail.append(
            f"  SKIP key={k!r}: CD20={len(obs_cd20)}, CD68={len(obs_cd68)}"
        )

print(f"  Valid pairs: {len(valid_rows)},  Skipped: {len(skipped_detail)}")
if skipped_detail:
    show = skipped_detail[:5]
    for s in show:
        print(s)
    if len(skipped_detail) > 5:
        print(f"  ... and {len(skipped_detail) - 5} more")

if len(valid_rows) == 0:
    # 诊断：打印 meta 和 paired 的 key 样例，帮助排查格式不一致
    print(f"  [DEBUG] meta _key samples:   {meta['_key'].head(5).tolist()}")
    print(f"  [DEBUG] paired _key samples: {paired['_key'].head(5).tolist()}")
    raise ValueError(
        "没有找到有效配对 ROI。\n"
        "请对照 [DEBUG] 输出，检查 Slide_ID / ROI_Num 格式是否一致。\n"
        "常见问题：paired 的 ROI_Num 是整数 1，而 expr index 里是 '001'"
    )

mat_cd20 = np.array(mat_cd20, dtype=float)
mat_cd68 = np.array(mat_cd68, dtype=float)
valid_df  = pd.DataFrame(valid_rows).reset_index(drop=True)

# ══════════════════════════════════════════
# 8. Z-score（pooled）+ 层次聚类
# ══════════════════════════════════════════
combined   = np.vstack([mat_cd20, mat_cd68])
gene_mean  = combined.mean(axis=0)
gene_std   = combined.std(axis=0) + 1e-8
mat_cd20_z = (mat_cd20 - gene_mean) / gene_std
mat_cd68_z = (mat_cd68 - gene_mean) / gene_std

link  = linkage(pdist(mat_cd68_z, metric="euclidean"), method="ward")
order = leaves_list(link)
mat_cd20_z = mat_cd20_z[order]
mat_cd68_z = mat_cd68_z[order]
valid_df   = valid_df.iloc[order].reset_index(drop=True)

# ══════════════════════════════════════════
# 9. Annotation bar encoding
# ══════════════════════════════════════════
def encode_group(row):
    if row["LA_TAM_Group"] == "LA_TAM_Rich" and row["Tissue_Type"] == "DLBCL":  return 0
    if row["LA_TAM_Group"] == "LA_TAM_Rich" and row["Tissue_Type"] == "Normal": return 1
    if row["LA_TAM_Group"] == "LA_TAM_Poor" and row["Tissue_Type"] == "DLBCL":  return 2
    return 3

bar_vals = valid_df.apply(encode_group, axis=1).values.reshape(1, -1).astype(float)
bar_cmap = ListedColormap(["#C0392B", "#E74C3C", "#2980B9", "#5DADE2"])

# ══════════════════════════════════════════
# 10. 绘图
# ══════════════════════════════════════════
print("Plotting...")
fig = plt.figure(figsize=(18, 9))
fig.patch.set_facecolor("#FFFFFF")

gs = gridspec.GridSpec(
    3, 3, figure=fig,
    height_ratios=[0.06, 1, 1],
    width_ratios=[0.03, 1, 0.03],
    hspace=0.04, wspace=0.03
)
ax_bar  = fig.add_subplot(gs[0, 1])
ax_cd20 = fig.add_subplot(gs[1, 1])
ax_cd68 = fig.add_subplot(gs[2, 1])
ax_cbar = fig.add_subplot(gs[1:, 2])

norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)

# Annotation bar
ax_bar.imshow(bar_vals, aspect="auto", cmap=bar_cmap,
              vmin=-0.5, vmax=3.5, interpolation="nearest")
ax_bar.set_xticks([]); ax_bar.set_yticks([])
ax_bar.set_ylabel("Group", fontsize=7, rotation=0, labelpad=30, va="center")

# CD20 heatmap
im = ax_cd20.imshow(mat_cd20_z.T, aspect="auto", cmap="RdBu_r",
                    norm=norm, interpolation="nearest")
ax_cd20.set_yticks(range(len(HEATMAP_GENES)))
ax_cd20.set_yticklabels(HEATMAP_GENES, fontsize=7.5)
ax_cd20.set_xticks([])
ax_cd20.set_ylabel("CD20+ ROIs", fontsize=9, fontweight="bold")

# CD68 heatmap
ax_cd68.imshow(mat_cd68_z.T, aspect="auto", cmap="RdBu_r",
               norm=norm, interpolation="nearest")
ax_cd68.set_yticks(range(len(HEATMAP_GENES)))
ax_cd68.set_yticklabels(HEATMAP_GENES, fontsize=7.5)
ax_cd68.set_xlabel(
    "Paired ROIs (ordered by CD68 hierarchical clustering)", fontsize=9
)
ax_cd68.set_ylabel("CD68+ ROIs", fontsize=9, fontweight="bold")

# Colorbar
plt.colorbar(im, cax=ax_cbar, label="Z-score")
ax_cbar.tick_params(labelsize=8)

# Legend
patches = [
    mpatches.Patch(color="#C0392B", label="LA_TAM_Rich · DLBCL"),
    mpatches.Patch(color="#2471A3", label="LA_TAM_Rich · Normal"),
    mpatches.Patch(color="#F1948A", label="LA_TAM_Poor · DLBCL"),
    mpatches.Patch(color="#AED6F1", label="LA_TAM_Poor · Normal"),
]
fig.legend(handles=patches, loc="upper center", fontsize=10,
           framealpha=0.92, edgecolor="#CCCCCC",
           bbox_to_anchor=(0.5, 0.93), ncol=4)

fig.suptitle(
    "Dual-Chamber Heatmap: CD20+ vs CD68+ Paired ROIs\n"
    "GSE232853 · LA-TAM Markers + Top DEGs · Z-scored Q3-Normalized Expression",
    fontsize=11, fontweight="bold", y=1.01
)

# ══════════════════════════════════════════
# 11. 保存
# ══════════════════════════════════════════
out_png = os.path.join(OUTPUT_DIR, "task5_dual_chamber_heatmap.png")
out_svg = os.path.join(OUTPUT_DIR, "task5_dual_chamber_heatmap.svg")
plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor="#FFFFFF")
plt.savefig(out_svg, bbox_inches="tight", facecolor="#FFFFFF")
plt.close()
print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
print("Done.")



# ═══════════════════════════════════════════════════════════════
# EXTRA EXPORTS — Results Section 5 Quantitative Tables
# Append after print("Done.") — reuses all already-loaded objects
# Output: D:\bulk-download\GSE232853_v2\figures_reproduced\revision_tables\
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("EXTRA EXPORT: Results Section 5 tables")
print("="*60)

REV_OUT = os.path.join(OUT, "revision_tables")
os.makedirs(REV_OUT, exist_ok=True)



# ────────────────────────────────────────────────────────────
# TABLE 2 — S9 Paired ROI Summary Stats
#   HMGB1 / HAVCR2 / MAFB: CD20 vs CD68 within paired ROIs
#   Groups: All / DLBCL / Normal
#   Values: n, median, q25, q75, median_diff, MWU p
# ────────────────────────────────────────────────────────────
print("\n[2/5] Exporting S9 paired ROI summary stats...")

# Auto-detect CD20/CD68 paired columns for each gene
def _find_paired_col(df, mask, gene):
    """Find column like CD20_HMGB1 or HMGB1_CD20 (case-insensitive)."""
    for c in df.columns:
        cl = c.lower()
        if mask.lower() in cl and gene.lower() in cl:
            return c
    return None

_s9_genes = ["HMGB1", "HAVCR2", "MAFB"]
_s9_rows  = []
_s9_missing = []

for _gene in _s9_genes:
    _col_cd20 = _find_paired_col(paired, "CD20", _gene)
    _col_cd68 = _find_paired_col(paired, "CD68", _gene)

    if _col_cd20 is None or _col_cd68 is None:
        _s9_missing.append(_gene)
        continue

    for _tname, _tmask in [("All",    np.ones(len(paired), dtype=bool)),
                           ("DLBCL",  _dlbcl),
                           ("Normal", _norm)]:
        _sub = paired[_tmask][[_col_cd20, _col_cd68, _tissue_c]].dropna(
            subset=[_col_cd20, _col_cd68])
        _n = len(_sub)
        if _n == 0:
            continue
        _a = _sub[_col_cd20].astype(float).values
        _b = _sub[_col_cd68].astype(float).values
        _diff = _b - _a
        try:
            _, _wp = mannwhitneyu(_a, _b, alternative="two-sided")
        except Exception:
            _wp = float("nan")
        _s9_rows.append({
            "gene":                    _gene,
            "tissue":                  _tname,
            "n_pairs":                 _n,
            "median_CD20":             round(float(np.median(_a)), 4),
            "q25_CD20":                round(float(np.percentile(_a, 25)), 4),
            "q75_CD20":                round(float(np.percentile(_a, 75)), 4),
            "median_CD68":             round(float(np.median(_b)), 4),
            "q25_CD68":                round(float(np.percentile(_b, 25)), 4),
            "q75_CD68":                round(float(np.percentile(_b, 75)), 4),
            "median_diff_CD68_minus_CD20": round(float(np.median(_diff)), 4),
            "MWU_p":                   float(_wp),
            "col_CD20":                _col_cd20,
            "col_CD68":                _col_cd68,
        })

if _s9_missing:
    print(f"  [WARNING] Could not find paired columns for: {_s9_missing}")
    print(f"  Available paired columns: {[c for c in paired.columns if 'CD20' in c or 'CD68' in c]}")

if _s9_rows:
    _s9_df = pd.DataFrame(_s9_rows)
    _s9_df.to_csv(os.path.join(OUT, "paired_roi_summary_stats.csv"), index=False)
    print(f"  ✓ paired_roi_summary_stats.csv  ({len(_s9_df)} rows)")
else:
    print("  [SKIP] No paired summary rows generated — check column names above")

# ────────────────────────────────────────────────────────────
# TABLE 3 — S10 GSEA Top Pathways
#   LA_TAM_Rich vs LA_TAM_Poor, Hallmark
#   All significant pathways (FDR < 0.25), sorted by NES
#   Values: Term, NES (2dp), FDR q-val, Direction
# ────────────────────────────────────────────────────────────
print("\n[3/5] Exporting S10 GSEA top pathways...")

_gsea_sig = gsea_df[gsea_df["FDR q-val"] < 0.25].copy()
_gsea_up  = _gsea_sig[_gsea_sig["NES"] > 0].nlargest(10, "NES")
_gsea_dn  = _gsea_sig[_gsea_sig["NES"] < 0].nsmallest(8,  "NES")
_gsea_out = pd.concat([_gsea_dn, _gsea_up]).sort_values("NES", ascending=False).reset_index(drop=True)

_gsea_out["Term_clean"]      = _gsea_out["Term"].apply(clean_term)
_gsea_out["NES_2dp"]         = _gsea_out["NES"].apply(lambda x: round(x, 2))
_gsea_out["FDR_formatted"]   = _gsea_out["FDR q-val"].apply(
    lambda q: "< 0.001" if q < 0.001 else f"{q:.3f}")
_gsea_out["Direction"]       = _gsea_out["NES"].apply(
    lambda x: "Enriched_LA_TAM_Rich" if x > 0 else "Enriched_LA_TAM_Poor")

_gsea_export_cols = ["Term", "Term_clean", "NES", "NES_2dp",
                     "FDR q-val", "FDR_formatted", "Direction"]
_gsea_export_cols = [c for c in _gsea_export_cols if c in _gsea_out.columns]
_gsea_out[_gsea_export_cols].to_csv(
    os.path.join(REV_OUT, "gsea_top_pathways_for_results.csv"), index=False)
print(f"  ✓ gsea_top_pathways_for_results.csv  ({len(_gsea_out)} pathways)")

# ────────────────────────────────────────────────────────────
# TABLE 4 — S11 Key DEA Genes (HAVCR2, MAFB, HMGB1)
#   CD68+ DLBCL vs Normal
#   Values: gene, log2FC (2dp), padj, Direction
# ────────────────────────────────────────────────────────────
print("\n[4/5] Exporting S11 key DEA genes...")

_lfc_col  = [c for c in dea.columns if any(k in c.lower() for k in
             ["logfold", "log2fc", "logfc", "lfc"])][0]
_pval_col = [c for c in dea.columns if any(k in c.lower() for k in
             ["padj", "pval_adj", "p_val_adj", "fdr", "pvals_adj"])][0]
_gene_col = [c for c in dea.columns if any(k in c.lower() for k in
             ["gene", "names"])][0]

_key_genes = ["HAVCR2", "MAFB", "HMGB1"]
_dea_sub = dea[dea[_gene_col].isin(_key_genes)].copy()

if len(_dea_sub) == 0:
    # Try case-insensitive match
    _dea_sub = dea[dea[_gene_col].str.upper().isin([g.upper() for g in _key_genes])].copy()

_dea_sub["log2FC_2dp"]    = _dea_sub[_lfc_col].apply(lambda x: round(float(x), 2))
_dea_sub["padj_formatted"] = _dea_sub[_pval_col].apply(
    lambda x: "< 0.001" if float(x) < 0.001 else f"{float(x):.3e}")
_dea_sub["Direction"] = _dea_sub[_lfc_col].apply(
    lambda x: "Up_in_DLBCL" if float(x) > 0 else "Down_in_DLBCL")

_dea_export_cols = [_gene_col, _lfc_col, "log2FC_2dp",
                    _pval_col, "padj_formatted", "Direction"]
_dea_export_cols = [c for c in _dea_export_cols if c in _dea_sub.columns]
_dea_sub[_dea_export_cols].to_csv(
    os.path.join(REV_OUT, "dea_key_genes_for_results.csv"), index=False)
print(f"  ✓ dea_key_genes_for_results.csv  ({len(_dea_sub)} genes found)")

if len(_dea_sub) < len(_key_genes):
    _found = _dea_sub[_gene_col].tolist()
    _missing_genes = [g for g in _key_genes if g not in _found]
    print(f"  [WARNING] Not found in DEA: {_missing_genes}")
    print(f"  DEA gene column sample: {dea[_gene_col].head(10).tolist()}")





# ═══════════════════════════════════════════════════════════════
# REBUILD: Fig2 & S5C — 使用正确的 Q3 标准化数据重绘
# 数据来源：
#   expression_raw_filtered.csv  (基因 × 样本，原始计数)
#   meta_filtered.csv            (样本元数据)
#   Q3_normalization_stats.csv   (每样本 Q3 值及 scale_factor)
# ═══════════════════════════════════════════════════════════════

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

BASE_DIR = translate(r"D:\bulk-download\GSE232853_v2")
OUT      = os.path.join(BASE_DIR, "figures_reproduced")
os.makedirs(OUT, exist_ok=True)

MASK_COLORS   = {"CD20": "#2166AC", "CD68": "#D6604D"}
TISSUE_COLORS = {"DLBCL": "#B2182B", "Normal": "#4393C3"}

def sig_stars(p):
    if p < 0.0001: return "****"
    elif p < 0.001: return "***"
    elif p < 0.01:  return "**"
    elif p < 0.05:  return "*"
    return "ns"

def savefig(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {name}")

# ─────────────────────────────────────────────────────────────
# Step 1: 载入原始计数矩阵与元数据
# ─────────────────────────────────────────────────────────────
print("Loading raw expression matrix...")

# expression_raw_filtered.csv: 行=基因, 列=样本
df_raw = pd.read_csv(
    os.path.join(BASE_DIR, "expression_raw_filtered.csv"),
    index_col=0
)
# 转置为 样本×基因
expr = df_raw.T.copy()
print(f"  Expression: {expr.shape[0]} samples × {expr.shape[1]} genes")

# meta_filtered.csv
meta = pd.read_csv(os.path.join(BASE_DIR, "meta_filtered.csv"))
meta.columns = [c.strip() for c in meta.columns]

# 检测 Sample_ID 列（兼容 meta 里可能存在的重复/多余列）
sid_col = [c for c in meta.columns if "sample_id" in c.lower()][0]
meta = meta.set_index(sid_col)
# 去除重复 index（meta 里有时会有 Sample_ID.1 等冗余列）
meta = meta[~meta.index.duplicated(keep="first")]

# Q3_normalization_stats.csv
q3_stats = pd.read_csv(os.path.join(BASE_DIR, "Q3_normalization_stats.csv"))
q3_stats.columns = [c.strip() for c in q3_stats.columns]
q3_stats = q3_stats.set_index("Sample_ID")

# ─────────────────────────────────────────────────────────────
# Step 2: Q3 标准化 + log1p
#   公式：X_q3 = (X / Q3_i) * median_Q3
#   若 Q3_normalization_stats.csv 已存有 Q3_value，直接使用；
#   否则从原始矩阵重新计算。
# ─────────────────────────────────────────────────────────────
print("Applying Q3 normalization...")

X_raw = expr.values.astype(np.float64)
sample_ids = expr.index.tolist()

if "Q3_value" in q3_stats.columns:
    # 优先使用已保存的 Q3 值（与原始流程完全一致）
    q3_per_sample = np.array([
        q3_stats.loc[sid, "Q3_value"] if sid in q3_stats.index else np.nan
        for sid in sample_ids
    ])
    # 对缺失样本回退到现场计算
    for i, sid in enumerate(sample_ids):
        if np.isnan(q3_per_sample[i]):
            nonzero = X_raw[i, X_raw[i, :] > 0]
            q3_per_sample[i] = np.percentile(nonzero, 75) if len(nonzero) > 0 else 1.0
            print(f"  [WARN] Q3 not found for {sid!r}, computed from raw: {q3_per_sample[i]:.2f}")
else:
    # 从原始矩阵重新计算
    q3_per_sample = np.zeros(X_raw.shape[0])
    for i in range(X_raw.shape[0]):
        nonzero = X_raw[i, X_raw[i, :] > 0]
        q3_per_sample[i] = np.percentile(nonzero, 75) if len(nonzero) > 0 else 1.0

median_q3 = np.median(q3_per_sample)
print(f"  Q3 range: {q3_per_sample.min():.2f} – {q3_per_sample.max():.2f}")
print(f"  Median Q3 (scale target): {median_q3:.2f}")

# 标准化
X_q3  = X_raw / q3_per_sample[:, np.newaxis] * median_q3
X_log = np.log1p(X_q3)
print(f"  Q3 log1p range: {X_log.min():.3f} – {X_log.max():.3f}")

# 重建为 DataFrame（样本×基因），列名与原始一致
expr_q3 = pd.DataFrame(X_log, index=sample_ids, columns=expr.columns)

# ─────────────────────────────────────────────────────────────
# Step 3: 构建 tgt_expr（与原始代码接口完全兼容）
#   格式：行=样本，列包含基因列 + Mask 列 + Tissue_Type 列
# ─────────────────────────────────────────────────────────────
print("Building tgt_expr compatible table...")

# 对齐 meta 与 expr_q3（取交集，保证顺序一致）
common_samples = [s for s in sample_ids if s in meta.index]
missing = len(sample_ids) - len(common_samples)
if missing > 0:
    print(f"  [WARN] {missing} samples not found in meta, dropped.")

expr_q3 = expr_q3.loc[common_samples]
meta_aligned = meta.loc[common_samples]

# 检测 Mask / Tissue_Type 列名（兼容大小写变体）
mask_col_meta   = [c for c in meta_aligned.columns if "mask"   in c.lower()][0]
tissue_col_meta = [c for c in meta_aligned.columns if "tissue" in c.lower()][0]

tgt_expr = expr_q3.copy()
tgt_expr["Mask"]        = meta_aligned[mask_col_meta].values
tgt_expr["Tissue_Type"] = meta_aligned[tissue_col_meta].values

# 目标三基因的列名（直接使用基因名，不加后缀）
GENES    = ["HMGB1", "HAVCR2", "MAFB"]
GENE_COL = {}
for g in GENES:
    if g in tgt_expr.columns:
        GENE_COL[g] = g
    else:
        # 兼容带后缀的列名（如 HMGB1_logCPM）
        candidates = [c for c in tgt_expr.columns if c.startswith(g)]
        if candidates:
            GENE_COL[g] = candidates[0]
            print(f"  [INFO] Gene {g} mapped to column: {GENE_COL[g]}")
        else:
            raise ValueError(f"Gene {g} not found in expression matrix. "
                             f"Check column names: {tgt_expr.columns[:10].tolist()}")

mask_c   = "Mask"
tissue_c = "Tissue_Type"

cd20_idx = tgt_expr[mask_c] == "CD20"
cd68_idx = tgt_expr[mask_c] == "CD68"

print(f"  CD20 ROIs: {cd20_idx.sum()},  CD68 ROIs: {cd68_idx.sum()}")
print(f"  DLBCL ROIs: {(tgt_expr[tissue_c]=='DLBCL').sum()},  "
      f"Normal ROIs: {(tgt_expr[tissue_c]=='Normal').sum()}")

# ═══════════════════════════════════════════════════════════════
# FIG 2 — Violin: HMGB1 / HAVCR2 / MAFB (CD20 vs CD68)
#          使用正确的 Q3 标准化数据重绘
# ═══════════════════════════════════════════════════════════════
print("\nFig2 (corrected): Violin HMGB1/HAVCR2/MAFB...")

stats_rows = []
for g in GENES:
    col = GENE_COL[g]
    a = tgt_expr.loc[cd20_idx, col].values.astype(float)
    b = tgt_expr.loc[cd68_idx, col].values.astype(float)
    stat, pval = mannwhitneyu(a, b, alternative="two-sided")
    stats_rows.append({
        "Gene": g, "col": col,
        "median_CD20": np.median(a),
        "median_CD68": np.median(b),
        "p_value": pval
    })

stats_df = pd.DataFrame(stats_rows)
_, padj, _, _ = multipletests(stats_df["p_value"], method="fdr_bh")
stats_df["padj_BH"] = padj

rng = np.random.default_rng(42)
fig, axes = plt.subplots(1, 3, figsize=(13, 6))
fig.patch.set_facecolor("#FFFFFF")

for ax, gene, row in zip(axes, GENES, stats_df.itertuples()):
    col = GENE_COL[gene]
    vals = {
        "CD20": tgt_expr.loc[cd20_idx, col].values.astype(float),
        "CD68": tgt_expr.loc[cd68_idx, col].values.astype(float)
    }
    vp = ax.violinplot(
        [vals["CD20"], vals["CD68"]], positions=[0, 1],
        widths=0.65, showmedians=False, showextrema=False
    )
    for body, key in zip(vp["bodies"], ["CD20", "CD68"]):
        body.set_facecolor(MASK_COLORS[key]); body.set_alpha(0.30)
        body.set_edgecolor(MASK_COLORS[key]); body.set_linewidth(0.8)

    bp = ax.boxplot(
        [vals["CD20"], vals["CD68"]], positions=[0, 1], widths=0.18,
        patch_artist=True,
        medianprops=dict(color="#222222", linewidth=2.0),
        whiskerprops=dict(color="#555555", linewidth=1.0),
        capprops=dict(color="#555555", linewidth=1.0),
        flierprops=dict(marker="", markersize=0), zorder=4
    )
    for patch, key in zip(bp["boxes"], ["CD20", "CD68"]):
        patch.set_facecolor(MASK_COLORS[key]); patch.set_alpha(0.55)
        patch.set_edgecolor(MASK_COLORS[key])

    for pos, (key, v) in zip([0, 1], vals.items()):
        jit = rng.uniform(-0.22, 0.22, size=len(v))
        ax.scatter(pos + jit, v, c=MASK_COLORS[key], s=5, alpha=0.35,
                   linewidths=0, zorder=3, rasterized=True)

    y_max = max(np.max(vals["CD20"]), np.max(vals["CD68"]))
    y_br  = y_max * 1.06
    s     = sig_stars(row.padj_BH)
    ax.plot([0, 0, 1, 1], [y_br, y_br + 0.12, y_br + 0.12, y_br],
            lw=1.2, c="#333333", zorder=5)
    ax.text(0.5, y_br + 0.14, s, ha="center", va="bottom",
            fontsize=15, fontweight="bold", color="#333333", zorder=5)
    ax.text(0.5, 0.9, f"padj = {row.padj_BH:.2e}",
            transform=ax.transAxes, ha="center", fontsize=10, color="#555555")

    n20 = cd20_idx.sum(); n68 = cd68_idx.sum()
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"CD20+\n(n={n20})", f"CD68+\n(n={n68})"], fontsize=14)
    ax.set_ylabel("Q3 log-norm expression", fontsize=9)
    ax.set_title(gene, fontsize=15, fontweight="bold", pad=6)
    ax.set_facecolor("#FFFFFF"); ax.tick_params(labelsize=8)
    for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")

handles = [mpatches.Patch(color=MASK_COLORS[k], alpha=0.7, label=k)
           for k in ["CD20", "CD68"]]
fig.legend(handles=handles, loc="upper right", fontsize=9,
           framealpha=0.92, edgecolor="#CCCCCC", bbox_to_anchor=(0.99, 0.85))
fig.suptitle(
    "Core Communication Genes: CD20+ vs CD68+ ROIs\n"
    "GSE232853 · Wilcoxon rank-sum · BH-corrected padj · Q3 log-norm (corrected)",
    fontsize=11, fontweight="bold", y=1.02
)
plt.tight_layout()
savefig(fig, "fig2_violin_HMGB1_HAVCR2_MAFB_Q3corrected")

# ─────────────────────────────────────────────────────────────
# Export: Fig2 绘图原始数据
# 保存每个样本在三个基因上的 Q3 log-norm 表达值 + Mask 标签
# ─────────────────────────────────────────────────────────────
print("Saving Fig2 plot data...")

fig2_cols = [GENE_COL["HMGB1"], GENE_COL["HAVCR2"], GENE_COL["MAFB"], mask_c]
fig2_data = tgt_expr[fig2_cols].copy()

# 统一列名，方便后续读取
fig2_data.columns = ["HMGB1", "HAVCR2", "MAFB", "Mask"]
fig2_data.index.name = "Sample_ID"
fig2_data.insert(fig2_data.columns.get_loc("Mask") + 1,
                 "Tissue_Type",
                 tgt_expr[tissue_c].values)


out_fig2 = os.path.join(OUT, "fig2_plot_data.csv")
fig2_data.to_csv(out_fig2)
print(f"  ✓ Fig2 plot data saved: {out_fig2}")
print(f"  Shape: {fig2_data.shape[0]} samples × {fig2_data.shape[1]} columns")
print(f"  CD20: {(fig2_data['Mask']=='CD20').sum()},  CD68: {(fig2_data['Mask']=='CD68').sum()}")

# ═══════════════════════════════════════════════════════════════
# Export Table S8 — AOI-level descriptive statistics and
# Mann–Whitney comparison for HMGB1 / HAVCR2 / MAFB
# ═══════════════════════════════════════════════════════════════
print("Exporting Table S8...")

# ── 直接复用前面 Step 3 已经建好的 GENE_COL 和 mask_c ──
# GENE_COL = {"HMGB1": "HMGB1", "HAVCR2": "HAVCR2", "MAFB": "MAFB"}  ← 已在上文定义
# mask_c = "Mask"  ← 已在上文定义

genes_s8  = ["HMGB1", "HAVCR2", "MAFB"]
rows_s8   = []
pvals_s8  = []
gene_p_map = {}

# ---- first pass: descriptive stats + raw p-values ----
for gene in genes_s8:
    col = GENE_COL[gene]   # 现在正确指向 "HMGB1" / "HAVCR2" / "MAFB"

    vals_cd20 = tgt_expr.loc[tgt_expr[mask_c] == "CD20", col].dropna().astype(float).values
    vals_cd68 = tgt_expr.loc[tgt_expr[mask_c] == "CD68", col].dropna().astype(float).values

    stat, pval = mannwhitneyu(vals_cd20, vals_cd68, alternative="two-sided")
    gene_p_map[gene] = pval
    pvals_s8.append(pval)

    for grp, vals in [("CD20", vals_cd20), ("CD68", vals_cd68)]:
        rows_s8.append({
            "Gene":    gene,
            "Group":   grp,
            "n":       len(vals),
            "mean":    float(np.mean(vals)),
            "median":  float(np.median(vals)),
            "std":     float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan,
            "q25":     float(np.percentile(vals, 25)),
            "q75":     float(np.percentile(vals, 75)),
            "p_value": pval
        })

# ---- BH correction ----
_, padj_vals, _, _ = multipletests(pvals_s8, method="fdr_bh")
gene_padj_map = dict(zip(genes_s8, padj_vals))
for row in rows_s8:
    row["padj_BH"] = gene_padj_map[row["Gene"]]

s8_long = pd.DataFrame(rows_s8)

# ---- wide format ----
wide_rows = []
for gene in genes_s8:
    sub = s8_long[s8_long["Gene"] == gene].set_index("Group")
    wide_rows.append({
        "Gene":        gene,
        "n_CD20":      int(sub.loc["CD20", "n"]),
        "mean_CD20":   float(sub.loc["CD20", "mean"]),
        "median_CD20": float(sub.loc["CD20", "median"]),
        "std_CD20":    float(sub.loc["CD20", "std"]),
        "q25_CD20":    float(sub.loc["CD20", "q25"]),
        "q75_CD20":    float(sub.loc["CD20", "q75"]),
        "n_CD68":      int(sub.loc["CD68", "n"]),
        "mean_CD68":   float(sub.loc["CD68", "mean"]),
        "median_CD68": float(sub.loc["CD68", "median"]),
        "std_CD68":    float(sub.loc["CD68", "std"]),
        "q25_CD68":    float(sub.loc["CD68", "q25"]),
        "q75_CD68":    float(sub.loc["CD68", "q75"]),
        "p_value":     float(sub.loc["CD20", "p_value"]),
        "padj_BH":     float(sub.loc["CD20", "padj_BH"]),
    })

s8_wide = pd.DataFrame(wide_rows)

# ---- save ----
out_s8_long = os.path.join(OUT, "Table_S8_AOI_expression_summary_long.csv")
out_s8_wide = os.path.join(OUT, "Table_S8_AOI_expression_summary_wide.csv")
s8_long.to_csv(out_s8_long, index=False)
s8_wide.to_csv(out_s8_wide, index=False)
print(f"  ✓ S8 long format saved: {out_s8_long}")
print(f"  ✓ S8 wide format saved: {out_s8_wide}")


# ═══════════════════════════════════════════════════════════════
# S5C — 2×2 Stratified Violin（正确 Q3 标准化版）
# Normal CD20 | Normal CD68 | DLBCL CD20 | DLBCL CD68
# ═══════════════════════════════════════════════════════════════
print("\nS5C (corrected): 2×2 stratified violin...")

tissue_tgt = "Tissue_Type"
mask_tgt   = "Mask"

GROUPS_S5C = [
    ("Normal", "CD20", "#AED6F1", "#2471A3"),
    ("Normal", "CD68", "#2471A3", "#2471A3"),
    ("DLBCL",  "CD20", "#F1948A", "#C0392B"),
    ("DLBCL",  "CD68", "#C0392B", "#C0392B"),
]
GROUP_LABELS = ["Normal CD20+", "Normal CD68+", "DLBCL CD20+", "DLBCL CD68+"]
POSITIONS    = [1, 2, 3.6, 4.6]

COMPARISONS = {
    "HMGB1":  [
        (0, 1, "Normal: CD20 vs CD68"),
        (2, 3, "DLBCL: CD20 vs CD68"),
        (0, 2, "CD20: Normal vs DLBCL"),
    ],
    "HAVCR2": [
        (0, 1, "Normal: CD20 vs CD68"),
        (2, 3, "DLBCL: CD20 vs CD68"),
        (1, 3, "CD68: Normal vs DLBCL"),
    ],
    "MAFB": [
        (0, 1, "Normal: CD20 vs CD68"),
        (2, 3, "DLBCL: CD20 vs CD68"),
        (1, 3, "CD68: Normal vs DLBCL"),
    ],
}

fig, axes = plt.subplots(1, 3, figsize=(17, 7))
fig.patch.set_facecolor("#FFFFFF")
rng_s5c = np.random.default_rng(42)

for ax, gene in zip(axes, GENES):
    ax.set_facecolor("#FFFFFF")
    col = GENE_COL[gene]

    # 提取四组数据
    group_vals = []
    for tissue_g, mask_g, *_ in GROUPS_S5C:
        idx = (
            (tgt_expr[tissue_tgt] == tissue_g) &
            (tgt_expr[mask_tgt]   == mask_g)
        )
        group_vals.append(tgt_expr.loc[idx, col].dropna().values)

    y_min_all = min(np.min(v) for v in group_vals)
    y_max_all = max(np.max(v) for v in group_vals)
    y_span    = y_max_all - y_min_all

    # Violin
    vp = ax.violinplot(
        group_vals, positions=POSITIONS,
        widths=0.60, showmedians=False, showextrema=False
    )
    for body, (_, _, fill_c, edge_c) in zip(vp["bodies"], GROUPS_S5C):
        body.set_facecolor(fill_c); body.set_edgecolor(edge_c)
        body.set_alpha(0.75); body.set_linewidth(1.2)

    # IQR bar + median dot + jitter
    for pos, vals, (_, _, fill_c, edge_c) in zip(POSITIONS, group_vals, GROUPS_S5C):
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.plot([pos, pos], [q1, q3], color=edge_c, lw=2.5, zorder=4)
        ax.scatter(pos, med, color="white", s=55, zorder=5,
                   edgecolors=edge_c, linewidths=1.8)
        jit = rng_s5c.uniform(-0.14, 0.14, size=len(vals))
        ax.scatter(pos + jit, vals, color=edge_c,
                   s=5, alpha=0.30, linewidths=0, zorder=3, rasterized=True)

    # 显著性 bracket
    bracket_fracs = [0.06, 0.16, 0.26]
    for (idx_a, idx_b, comp_lbl), y_frac in zip(COMPARISONS[gene], bracket_fracs):
        vals_a = group_vals[idx_a]
        vals_b = group_vals[idx_b]
        _, p   = mannwhitneyu(vals_a, vals_b, alternative="two-sided")

        p_txt  = f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.2e}"
        stars  = sig_stars(p)
        label_txt = p_txt if 0.05 < p < 1 else f"{stars}  {p_txt}"

        x_a  = POSITIONS[idx_a]; x_b = POSITIONS[idx_b]
        y_br = y_max_all + y_span * y_frac
        tip  = y_span * 0.022

        ax.plot([x_a, x_a, x_b, x_b],
                [y_br, y_br + tip, y_br + tip, y_br],
                lw=1.2, c="#333333", zorder=5)
        ax.text((x_a + x_b) / 2, y_br + tip + y_span * 0.008,
                label_txt, ha="center", va="bottom",
                fontsize=11, fontweight="bold", color="#333333")

    # 中间分隔线
    ax.axvline(x=2.8, color="#BBBBBB", lw=1.0, ls=":", zorder=1)

    ax.set_xticks(POSITIONS)
    ax.set_xticklabels(GROUP_LABELS, fontsize=9)
    ax.set_ylabel("Expression (Q3 log-norm)", fontsize=9)
    ax.set_title(gene, fontsize=15, fontweight="bold", pad=8)
    ax.tick_params(axis="y", labelsize=8)
    for sp in ax.spines.values(): sp.set_edgecolor("#CCCCCC")

    # n 标签
    for pos, vals in zip(POSITIONS, group_vals):
        ax.text(
            pos,
            ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05,
            f"n={len(vals)}", ha="center", va="top",
            fontsize=11, color="#555555"
        )

legend_handles = [
    mpatches.Patch(facecolor=fill_c, edgecolor=edge_c,
                   alpha=0.85, label=lbl, linewidth=1.2)
    for (_, _, fill_c, edge_c), lbl in zip(GROUPS_S5C, GROUP_LABELS)
]
fig.legend(
    handles=legend_handles, loc="upper center", ncol=4,
    fontsize=11, framealpha=0.92, edgecolor="#CCCCCC",
    bbox_to_anchor=(0.5, 1.02),
    title="Group (Tissue × Compartment)", title_fontsize=11
)
fig.suptitle(
    "S5C — HMGB1 / HAVCR2 / MAFB: Expression by Tissue × Compartment\n"
    "Normal vs DLBCL  ×  CD20+ vs CD68+  ·  GSE232853 · Q3 log-norm (corrected)",
    fontsize=11, fontweight="bold", y=1.08
)
plt.tight_layout()
savefig(fig, "S5C_2x2_stratified_violin_Q3corrected")

print("\nDone. Corrected figures saved to:", OUT)

# ─────────────────────────────────────────────────────────────
# Export: S5C 绘图原始数据
# 保存每个样本在三个基因上的 Q3 log-norm 表达值 + Mask + Tissue_Type
# ─────────────────────────────────────────────────────────────
print("Saving S5C plot data...")

s5c_cols = [GENE_COL["HMGB1"], GENE_COL["HAVCR2"], GENE_COL["MAFB"], mask_c, tissue_c]
s5c_data = tgt_expr[s5c_cols].copy()

# 统一列名
s5c_data.columns = ["HMGB1", "HAVCR2", "MAFB", "Mask", "Tissue_Type"]
s5c_data.index.name = "Sample_ID"

out_s5c = os.path.join(OUT, "S5C_plot_data.csv")
s5c_data.to_csv(out_s5c)
print(f"  ✓ S5C plot data saved: {out_s5c}")
print(f"  Shape: {s5c_data.shape[0]} samples × {s5c_data.shape[1]} columns")

# 打印四组样本量，方便核对
for tissue_g in ["Normal", "DLBCL"]:
    for mask_g in ["CD20", "CD68"]:
        n = ((s5c_data["Tissue_Type"] == tissue_g) & (s5c_data["Mask"] == mask_g)).sum()
        print(f"  {tissue_g} {mask_g}: n={n}")


print("\n[5/5] Exporting Fig2/S8 three-gene summary...")

_GENES_S8 = ["HMGB1", "HAVCR2", "MAFB"]
_GENE_COL_MAP = {
    "HMGB1":  "HMGB1_logCPM",
    "HAVCR2": "HAVCR2_logCPM",
    "MAFB":   "MAFB_logCPM",
}
_mask_c   = [c for c in tgt_expr.columns if "mask"   in c.lower()][0]
_cd20_idx = tgt_expr[_mask_c] == "CD20"
_cd68_idx = tgt_expr[_mask_c] == "CD68"

_s8_rows = []
for _g in _GENES_S8:
    _col = _GENE_COL_MAP[_g]
    if _col not in tgt_expr.columns:
        # Fallback: find column containing gene name
        _candidates = [c for c in tgt_expr.columns if _g.lower() in c.lower()]
        if not _candidates:
            print(f"  [WARNING] Column for {_g} not found. Skipping.")
            continue
        _col = _candidates[0]
    _a = tgt_expr.loc[_cd20_idx, _col].astype(float).values
    _b = tgt_expr.loc[_cd68_idx, _col].astype(float).values
    _, _pv = mannwhitneyu(_a, _b, alternative="two-sided")
    _s8_rows.append({
        "gene":           _g,
        "n_CD20":         len(_a),
        "median_CD20":    round(float(np.median(_a)), 4),
        "q25_CD20":       round(float(np.percentile(_a, 25)), 4),
        "q75_CD20":       round(float(np.percentile(_a, 75)), 4),
        "n_CD68":         len(_b),
        "median_CD68":    round(float(np.median(_b)), 4),
        "q25_CD68":       round(float(np.percentile(_b, 25)), 4),
        "q75_CD68":       round(float(np.percentile(_b, 75)), 4),
        "MWU_p_raw":      float(_pv),
        "expr_col_used":  _col,
    })

_s8_df = pd.DataFrame(_s8_rows)
if len(_s8_df) > 0:
    _, _padj_s8, _, _ = multipletests(_s8_df["MWU_p_raw"], method="fdr_bh")
    _s8_df["padj_BH"] = _padj_s8
    _s8_df["padj_formatted"] = _s8_df["padj_BH"].apply(
        lambda x: "< 0.001" if x < 0.001 else f"{x:.3e}")
    _s8_df["Direction"] = _s8_df.apply(
        lambda r: "Higher_in_CD20" if r["median_CD20"] > r["median_CD68"]
                  else "Higher_in_CD68", axis=1)
    _s8_df.to_csv(os.path.join(OUT, "fig2_three_gene_summary.csv"), index=False)
    print(f"  ✓ fig2_three_gene_summary.csv  ({len(_s8_df)} genes)")

    


# ─────────────────────────────────────────────────────────────
# Step 2.5: 保存 Q3 标准化后的表达矩阵
# 格式与原始 expression_raw_filtered.csv 一致：行=基因, 列=样本
# ─────────────────────────────────────────────────────────────
print("Saving Q3-normalized expression matrix...")

# expr_q3 当前是 样本×基因，转置回 基因×样本 再保存
expr_q3_save = pd.DataFrame(
    X_log,
    index=sample_ids,
    columns=expr.columns
).T  # 转置为 基因×样本

out_path = os.path.join(BASE_DIR, "expression_q3_lognorm.csv")
expr_q3_save.to_csv(out_path)
print(f"  ✓ Saved: {out_path}")
print(f"  Shape: {expr_q3_save.shape[0]} genes × {expr_q3_save.shape[1]} samples")

