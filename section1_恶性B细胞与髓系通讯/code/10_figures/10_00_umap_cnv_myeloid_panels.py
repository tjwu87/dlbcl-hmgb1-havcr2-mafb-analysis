# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# ============================================================
# GSE182434 本地绘图脚本
# 直接加载 Agent 保存的数据，复现全部 7 张图
# 修改任意 "# [PARAM]" 注释行下方的参数即可调整图形外观
# ============================================================

import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib

matplotlib.use("Agg")  # [PARAM] 非交互后端；改为 'TkAgg' 或注释掉此行可弹出窗口
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")

# ============================================================
# 0. 路径配置  ← 只需修改这里
# ============================================================
BASE_DIR = "."  # [PARAM] 项目根目录，改为你的实际路径
H5AD_PATH = translate(f"D:/bulk-download/GSE182434/adata_processed.h5ad")
MARKERS_CSV = translate(f"D:/bulk-download/GSE182434/comparison/cluster_markers.csv")
CNV_CSV = translate(f"D:/bulk-download/GSE182434/cnv/malignancy_classification.csv")
AUCELL_CSV = translate(f"D:/bulk-download/GSE182434/cnv/aucell_metabolic_scores.csv")

OUT_UMAP = translate(f"D:/bulk-download/GSE182434/umap")
OUT_ANNOT = translate(f"D:/bulk-download/GSE182434/annotation")
OUT_COMP = translate(f"D:/bulk-download/GSE182434/comparison")
OUT_CNV = translate(f"D:/bulk-download/GSE182434/cnv")

for d in [OUT_UMAP, OUT_ANNOT, OUT_COMP, OUT_CNV]:
    os.makedirs(d, exist_ok=True)

# [PARAM] 输出格式列表，可选 'png'、'svg'、'pdf'
SAVE_FORMATS = ["png", "svg"]

# [PARAM] 全局 DPI（分辨率），越高文件越大；屏幕预览用 100，出版用 300
DPI = 300


# ============================================================
# 1. 加载数据
# ============================================================
print("Loading adata_processed.h5ad ...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  adata: {adata.shape}")
print(f"  obs columns: {adata.obs.columns.tolist()}")
print(f"  obsm keys:   {list(adata.obsm.keys())}")
print(f"  layers:      {list(adata.layers.keys())}")

# UMAP 坐标
umap = adata.obsm["X_umap"]

# 加载 CSV 表格
markers_df = pd.read_csv(MARKERS_CSV)
cnv_meta = pd.read_csv(CNV_CSV, index_col=0)
aucell_df = pd.read_csv(AUCELL_CSV, index_col=0)

# 把 CNV/malignancy 信息合并回 adata.obs（仅 B 细胞有值）
adata.obs["cnv_score"] = cnv_meta["cnv_score"].reindex(adata.obs_names)
adata.obs["malignancy"] = cnv_meta["malignancy"].reindex(adata.obs_names)

print("Data loaded successfully.\n")


# ============================================================
# 公共调色板（所有图共享，统一修改颜色在这里）
# ============================================================

# [PARAM] Tissue 颜色：键=Tissue 名，值=十六进制颜色
tissue_colors = {
    "DLBCL": "#0279EE",
    "Tonsil": "#FF9400",
}

# [PARAM] Patient/batch 颜色
patient_colors = {
    "DLBCL002": "#1f77b4",
    "DLBCL007": "#ff7f0e",
    "DLBCL008": "#2ca02c",
    "DLBCL111": "#d62728",
    "T2": "#9467bd",
}

# [PARAM] GEO 原始粗粒度细胞类型颜色（用于 umap_overview Panel 4）
celltype_colors_geo = {
    "B cells": "#1f77b4",
    "T cells CD8": "#d62728",
    "T cells CD4": "#ff7f0e",
    "Tregs": "#9467bd",
    "TFH": "#8c564b",
    "Monocytes and Macrophages": "#e377c2",
    "NK cells": "#7f7f7f",
    "Plasma cells": "#bcbd22",
    "Others": "#17becf",
}

# [PARAM] 精细化细胞类型颜色（CellTypist 细分 B 细胞后，共 14 种）
celltype_palette = {
    "CD8 T cells": "#1f77b4",
    "CD4 T cells": "#aec7e8",
    "Tregs": "#9467bd",
    "TFH cells": "#c5b0d5",
    "NK cells": "#17becf",
    "Naive B cells": "#2ca02c",
    "GC B cells": "#98df8a",
    "Proliferative GC B cells": "#d62728",
    "Memory B cells": "#ff9896",
    "Age-associated B cells": "#8c564b",
    "B cells (other)": "#c49c94",
    "Monocytes/Macrophages": "#e377c2",
    "Plasma cells": "#ff7f0e",
    "pDC/Other": "#7f7f7f",
}

# [PARAM] 细胞类型排列顺序（影响堆叠柱状图、水平柱图、violin 图的顺序）
ct_order = [
    "CD8 T cells",
    "CD4 T cells",
    "Tregs",
    "TFH cells",
    "NK cells",
    "Naive B cells",
    "GC B cells",
    "Proliferative GC B cells",
    "Memory B cells",
    "Age-associated B cells",
    "B cells (other)",
    "Monocytes/Macrophages",
    "Plasma cells",
    "pDC/Other",
]

# [PARAM] 恶性分类颜色
group_palette = {
    "Normal B cell": "#2ca02c",
    "Normal B cell (DLBCL)": "#74b9ff",
    "Malignant B cell": "#d62728",
}
group_order = ["Normal B cell", "Normal B cell (DLBCL)", "Malignant B cell"]


# ============================================================
# 辅助函数：保存图片
# ============================================================
def save_fig(fig, path_no_ext):
    for fmt in SAVE_FORMATS:
        fig.savefig(f"{path_no_ext}.{fmt}", dpi=DPI, bbox_inches="tight", format=fmt)
    plt.close(fig)
    print(f"  ✓ Saved: {path_no_ext}.{{{','.join(SAVE_FORMATS)}}}")


# ============================================================
# 图 1: umap_overview — 4-panel UMAP
# ============================================================
print("=== Figure 1: UMAP Overview ===")

# [PARAM] 整体图尺寸 (宽, 高) 英寸
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# [PARAM] 总标题文字、字号、粗体
fig.suptitle(
    "GSE182434: DLBCL vs Tonsil scRNA-seq\n(Harmony batch-corrected)",
    fontsize=14,  # [PARAM] 标题字号
    fontweight="bold",
)

# --- Panel 1: Tissue ---
ax = axes[0, 0]
for tissue, color in tissue_colors.items():
    mask = adata.obs["Tissue"] == tissue
    ax.scatter(
        umap[mask, 0],
        umap[mask, 1],
        c=color,
        s=2,  # [PARAM] 点大小（越小越不遮挡，建议 1~5）
        alpha=0.4,  # [PARAM] 透明度 0~1（越小越透明）
        rasterized=True,
        label=tissue,
    )
ax.set_title("Tissue", fontsize=14)  # [PARAM] 子图标题字号
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(
    markerscale=4,  # [PARAM] 图例点放大倍数
    frameon=False,  # [PARAM] True=有边框，False=无边框
    fontsize=11,  # [PARAM] 图例字号
)
ax.set_aspect("equal")  # [PARAM] 'equal'=等比例坐标轴，'auto'=自适应

# --- Panel 2: Patient (batch) ---
ax = axes[0, 1]
for patient, color in patient_colors.items():
    mask = adata.obs["Patient"] == patient
    ax.scatter(
        umap[mask, 0],
        umap[mask, 1],
        c=color,
        s=2,
        alpha=0.4,
        rasterized=True,
        label=patient,
    )
ax.set_title("Patient (Batch)", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(markerscale=4, frameon=False, fontsize=11)
ax.set_aspect("equal")

# --- Panel 3: Leiden 0.6 clusters ---
ax = axes[1, 0]
clusters = adata.obs["leiden_0.6"].astype(str)
unique_clusters = sorted(clusters.unique(), key=lambda x: int(x))
# [PARAM] 聚类配色方案，可换为 'tab10'、'Set1'、'hsv' 等
cmap = plt.cm.get_cmap("tab20", len(unique_clusters))
for i, cl in enumerate(unique_clusters):
    mask = clusters == cl
    ax.scatter(
        umap[mask, 0],
        umap[mask, 1],
        c=[cmap(i)],
        s=2,
        alpha=0.5,
        rasterized=True,
        label=cl,
    )
ax.set_title("Leiden Clusters (res=0.6)", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(
    markerscale=4,
    frameon=False,
    fontsize=11,  # [PARAM] 聚类图例字号（聚类多时需调小）
    ncol=3,  # [PARAM] 图例列数
    title="Cluster",
)
ax.set_aspect("equal")

# --- Panel 4: Cell type (GEO 原始粗粒度标注) ---
ax = axes[1, 1]
for ct, color in celltype_colors_geo.items():
    mask = adata.obs["CellType"] == ct
    if mask.sum() > 0:
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=color,
            s=2,
            alpha=0.5,
            rasterized=True,
            label=ct,
        )
ax.set_title("Cell Type (GEO annotation)", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(markerscale=4, frameon=False, fontsize=9, title="Cell Type")
ax.set_aspect("equal")

plt.tight_layout()
save_fig(fig, f"{OUT_UMAP}/umap_overview2")


# ============================================================
# 图 2: umap_celltypes — 精细化细胞类型 + Tissue 2-panel UMAP
# ============================================================
print("=== Figure 2: Cell Type UMAP ===")

# [PARAM] 整体图尺寸
fig, axes = plt.subplots(1, 2, figsize=(9.8, 8.4),
                        gridspec_kw=dict(wspace=0.25))
fig.suptitle(
    "GSE182434: Cell Type Annotation\n(GEO labels + CellTypist B cell refinement)",
    fontsize=13,
    fontweight="bold",
)

# --- Panel 1: 精细化细胞类型 ---
ax = axes[0]
for ct, color in celltype_palette.items():
    mask = adata.obs["cell_type"] == ct
    if mask.sum() > 0:
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=color,
            s=2,
            alpha=0.5,
            rasterized=True,
            label=f"{ct} ({mask.sum()})",  # [PARAM] 图例格式，可去掉 (n) 改为只写 ct
        )
ax.set_title("Refined Cell Types", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(
    markerscale=5,
    frameon=False,
    fontsize=8.8,
    loc="upper center",
    bbox_to_anchor=(1.075, -0.11),  # 锚点右移：5 列图例以整幅画布为中心（左轴只占 44% 宽）
    ncol=5,
    columnspacing=0.5,
    handlelength=1.0,
    title="Cell Type (n)",
)
ax.set_aspect("equal")

# --- Panel 2: Tissue 覆盖 ---
ax = axes[1]
for tissue, color in tissue_colors.items():
    mask = adata.obs["Tissue"] == tissue
    ax.scatter(
        umap[mask, 0],
        umap[mask, 1],
        c=color,
        s=2,
        alpha=0.35,  # [PARAM] Tissue panel 透明度（稍低以显示重叠）
        rasterized=True,
        label=f"{tissue} (n={mask.sum()})",
    )
ax.set_title("Tissue (DLBCL vs Tonsil)", fontsize=14)
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.legend(markerscale=5, frameon=False, fontsize=12, title="Tissue", loc="lower right")
ax.set_aspect("equal")

plt.tight_layout()
save_fig(fig, f"{OUT_ANNOT}/umap_celltypes2")


# ============================================================
# 图 3: marker_genes — 12 个 marker 基因表达 UMAP（3×4）
# ============================================================
print("=== Figure 3: Marker Gene Expression ===")

# [PARAM] 要展示的 marker 基因字典：键=基因名，值=注释标签
markers = {
    "CD3D": "T cells",
    "CD8A": "CD8 T cells",
    "CD4": "CD4 T cells",
    "FOXP3": "Tregs",
    "CXCR5": "TFH",
    "MS4A1": "B cells (CD20)",
    "CD38": "Plasma/GC B",
    "MKI67": "Proliferating",
    "LYZ": "Monocytes/Mac",
    "GNLY": "NK cells",
    "IGHG1": "Plasma cells",
    "BCL6": "GC B cells",
}

X_log = adata.layers["log1p_norm"]
gene_names = adata.var_names.tolist()

# [PARAM] 整体图尺寸
fig, axes = plt.subplots(3, 4, figsize=(18, 13))
fig.suptitle("Marker Gene Expression (log-normalized)", fontsize=13, fontweight="bold")

for ax, (gene, label) in zip(axes.flatten(), markers.items()):
    if gene in gene_names:
        idx = gene_names.index(gene)
        expr = (
            np.array(X_log[:, idx].todense()).flatten()
            if sp.issparse(X_log)
            else X_log[:, idx]
        )

        # 按表达量排序，高表达点画在最上层
        order = np.argsort(expr)
        sc_plot = ax.scatter(
            umap[order, 0],
            umap[order, 1],
            c=expr[order],
            s=1.5,  # [PARAM] 基因表达图点大小
            alpha=0.7,  # [PARAM] 基因表达图透明度
            cmap="YlOrRd",  # [PARAM] 颜色映射，可换 'Blues'、'viridis'、'RdPu' 等
            rasterized=True,
            vmin=0,  # [PARAM] 颜色映射最小值（0=从无表达开始）
        )
        plt.colorbar(
            sc_plot,
            ax=ax,
            shrink=0.7,  # [PARAM] colorbar 相对高度 0~1
            pad=0.02,
        )
        ax.set_title(f"{gene}\n({label})", fontsize=9)
    else:
        ax.set_title(f"{gene} (not found)", fontsize=9)
        ax.text(
            0.5, 0.5, "Not in dataset", ha="center", va="center", transform=ax.transAxes
        )

    ax.set_xlabel("UMAP1", fontsize=7)
    ax.set_ylabel("UMAP2", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_aspect("equal")

plt.tight_layout()
save_fig(fig, f"{OUT_ANNOT}/marker_genes2")


# ============================================================
# 图 4 & 5: 细胞类型组成图（堆叠柱状图 + 平均比例水平柱图）
# ============================================================
print("=== Figure 4: Per-sample stacked bar ===")
print("=== Figure 5: Mean proportion comparison ===")

# --- 计算每样本各细胞类型比例 ---
obs = adata.obs[["Sample", "Tissue", "Patient", "cell_type"]].copy()
counts_df = obs.groupby(["Sample", "cell_type"]).size().reset_index(name="n")
totals = counts_df.groupby("Sample")["n"].sum().rename("total")
counts_df = counts_df.join(totals, on="Sample")
counts_df["proportion"] = counts_df["n"] / counts_df["total"]

wide = counts_df.pivot_table(
    index="Sample",
    columns="cell_type",
    values="proportion",
    aggfunc="first",
    fill_value=0,
)
meta_df = obs[["Sample", "Tissue", "Patient"]].drop_duplicates().set_index("Sample")
wide = wide.join(meta_df).reset_index()
wide = wide.sort_values(["Tissue", "Patient"]).reset_index(drop=True)
ct_cols = [c for c in ct_order if c in wide.columns]

n = len(wide)
n_dlbcl = (wide["Tissue"] == "DLBCL").sum()


def short_label(row):
    """X 轴样本标签格式，可按需修改"""
    s = row["Sample"]
    suffix = "NB" if s.endswith("NB") else "B"
    return f"{row['Patient']}\n({suffix})"


sample_labels = [short_label(row) for _, row in wide.iterrows()]

# ---- 图 4: 每样本堆叠柱状图 ----
# [PARAM] 整体图尺寸
fig1, ax = plt.subplots(figsize=(12, 5.5))
fig1.suptitle(
    "Per-Sample Cell Type Composition: DLBCL vs Tonsil", fontsize=13, fontweight="bold"
)

x = np.arange(n)
bottom = np.zeros(n)
for ct in ct_cols:
    vals = wide[ct].values
    ax.bar(
        x,
        vals,
        bottom=bottom,
        color=celltype_palette[ct],
        label=ct,
        width=0.72,  # [PARAM] 柱宽 0~1（越接近1越宽）
        linewidth=0,  # [PARAM] 柱边线宽度，0=无边线
    )
    bottom += vals

# [PARAM] 两组背景底色透明度
ax.axvspan(-0.5, n_dlbcl - 0.5, alpha=0.07, color="#0279EE", zorder=0)
ax.axvspan(n_dlbcl - 0.5, n - 0.5, alpha=0.07, color="#FF9400", zorder=0)
ax.axvline(
    x=n_dlbcl - 0.5,
    color="black",
    linewidth=1.5,
    linestyle="--",
    alpha=0.7,  # [PARAM] 分隔线样式
)

# [PARAM] 组标签位置和样式
ax.text(
    (n_dlbcl - 1) / 2,
    1.04,
    "DLBCL",
    ha="center",
    fontsize=14,
    fontweight="bold",
    color="#0279EE",
)
ax.text(
    n_dlbcl + (n - n_dlbcl - 1) / 2,
    1.04,
    "Tonsil",
    ha="center",
    fontsize=14,
    fontweight="bold",
    color="#FF9400",
)

ax.set_xticks(x)
ax.set_xticklabels(sample_labels, fontsize=11, ha="center")
ax.set_ylabel("Cell Type Proportion", fontsize=12)
ax.set_ylim(0, 1.10)  # [PARAM] Y 轴范围（上限留空间给组标签）
ax.set_xlim(-0.5, n - 0.5)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

handles = [mpatches.Patch(color=celltype_palette[ct], label=ct) for ct in ct_cols]
ax.legend(
    handles=handles,
    bbox_to_anchor=(1.01, 1),  # [PARAM] 图例位置
    loc="upper left",
    fontsize=11,  # [PARAM] 图例字号
    frameon=False,
    title="Cell Type",
    title_fontsize=12,
)

plt.tight_layout()
save_fig(fig1, f"{OUT_COMP}/celltype_composition_per_sample2")

# ---- 图 5: 平均比例水平柱图 ----
mean_by_tissue = wide.groupby("Tissue")[ct_cols].mean()
mean_pivot = mean_by_tissue.T.reindex(
    [c for c in ct_order if c in mean_by_tissue.columns]
)

# [PARAM] 整体图尺寸
fig2, ax2 = plt.subplots(figsize=(8, 7))
fig2.suptitle(
    "Mean Cell Type Proportion: DLBCL vs Tonsil", fontsize=13, fontweight="bold"
)

y_pos = np.arange(len(mean_pivot))
dlbcl_vals = mean_pivot.get("DLBCL", pd.Series(0, index=mean_pivot.index)).values
tonsil_vals = mean_pivot.get("Tonsil", pd.Series(0, index=mean_pivot.index)).values

bars_d = ax2.barh(
    y_pos - 0.2,
    dlbcl_vals,
    height=0.35,  # [PARAM] 水平柱高度
    color="#0279EE",  # [PARAM] DLBCL 柱颜色
    alpha=0.85,
    label="DLBCL",
)
bars_t = ax2.barh(
    y_pos + 0.2,
    tonsil_vals,
    height=0.35,
    color="#FF9400",  # [PARAM] Tonsil 柱颜色
    alpha=0.85,
    label="Tonsil",
)

# 数值标签（仅显示 > 0.01 的值）
for bar, val in zip(bars_d, dlbcl_vals):
    if val > 0.01:  # [PARAM] 显示数值标签的最小阈值
        ax2.text(
            val + 0.004,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.2f}",  # [PARAM] 数值格式，'.2f'=两位小数
            va="center",
            fontsize=7.5,
            color="#0055bb",
        )
for bar, val in zip(bars_t, tonsil_vals):
    if val > 0.01:
        ax2.text(
            val + 0.004,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.2f}",
            va="center",
            fontsize=7.5,
            color="#cc6600",
        )

ax2.set_yticks(y_pos)
ax2.set_yticklabels(mean_pivot.index, fontsize=10)
ax2.set_xlabel("Mean Proportion", fontsize=11)
ax2.legend(fontsize=10, frameon=False, loc="lower right")
ax2.invert_yaxis()
ax2.set_xlim(
    0, max(dlbcl_vals.max(), tonsil_vals.max()) * 1.3
)  # [PARAM] X 轴右侧留白倍数
for spine in ["top", "right"]:
    ax2.spines[spine].set_visible(False)

plt.tight_layout()
save_fig(fig2, f"{OUT_COMP}/celltype_mean_proportion2")


# ============================================================
# 图 6: metabolic_violin — 10 个代谢通路 violin 图（2×5）
# ============================================================
print("=== Figure 6: Metabolic Pathway Violin ===")

# [PARAM] 要展示的代谢通路（必须是 aucell_df 的列名）
key_pathways = [
    "Glycolysis_Gluconeogenesis",
    "Oxidative_Phosphorylation",
    "Glutathione_Metabolism",
    "Cysteine_Methionine_Metabolism",
    "Glutamate_Metabolism",
    "TCA_Cycle",
    "Pentose_Phosphate_Pathway",
    "HIF1_Signaling",
    "Hallmark_Glycolysis",
    "Hallmark_Oxidative_Phosphorylation",
]

# [PARAM] 通路显示名称（换行用 \n）
pathway_labels = {
    "Glycolysis_Gluconeogenesis": "Glycolysis /\nGluconeogenesis",
    "Oxidative_Phosphorylation": "Oxidative\nPhosphorylation",
    "Glutathione_Metabolism": "Glutathione\nMetabolism",
    "Cysteine_Methionine_Metabolism": "Cysteine &\nMethionine Metab.",
    "Glutamate_Metabolism": "Glutamate\nMetabolism",
    "TCA_Cycle": "TCA Cycle",
    "Pentose_Phosphate_Pathway": "Pentose Phosphate\nPathway",
    "HIF1_Signaling": "HIF-1 Signaling",
    "Hallmark_Glycolysis": "Hallmark\nGlycolysis",
    "Hallmark_Oxidative_Phosphorylation": "Hallmark\nOXPHOS",
}

# 把 malignancy 合并进 aucell_df
scores_df = aucell_df.copy()
# aucell_df 已含 malignancy 列（Agent 保存时已写入）


def add_bracket(ax, x1, x2, y_top, sig_text, color="black"):
    """在两组之间画显著性括号"""
    h = y_top * 0.035
    ax.plot(
        [x1, x1, x2, x2],
        [y_top - h, y_top, y_top, y_top - h],
        color=color,
        linewidth=1.1,
        clip_on=False,
    )
    ax.text(
        (x1 + x2) / 2,
        y_top + h * 0.2,
        sig_text,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        color=color,
    )


def sig_stars(p):
    """p 值转显著性星号"""
    # [PARAM] 显著性阈值，可自行调整
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


# [PARAM] 整体图尺寸（2行×5列，共10个子图）
fig, axes = plt.subplots(
    2,
    5,
    figsize=(22, 10),
    gridspec_kw={
        "hspace": 0.55,  # [PARAM] 行间距，越大行越分开
        "wspace": 0.38,  # [PARAM] 列间距
    },
)
axes_flat = axes.flatten()

fig.suptitle(
    "Metabolic Pathway Activity in DLBCL B Cells"
    "Malignant vs Normal B Cells (AUCell Scores, MSigDB/KEGG)",
    fontsize=14,
    fontweight="bold",
    y=1.02,  # [PARAM] 总标题垂直位置，>1 在图上方
)

for ax, pw in zip(axes_flat, key_pathways):
    sub = scores_df[["malignancy", pw]].rename(columns={pw: "score"})

    sns.violinplot(
        data=sub,
        x="malignancy",
        y="score",
        order=group_order,
        palette=group_palette,
        inner="box",  # [PARAM] violin 内部显示方式：'box'=箱线，'quartile'=四分位线，'point'=点，None=空
        linewidth=0.8,  # [PARAM] violin 轮廓线宽
        cut=0,  # [PARAM] violin 末端截断：0=截至数据范围，2=延伸2个带宽
        ax=ax,
        density_norm="width",  # [PARAM] violin 宽度归一化：'width'=等宽，'area'=等面积，'count'=按样本量
    )

    mal_s = sub[sub["malignancy"] == "Malignant B cell"]["score"].values
    norm_s = sub[sub["malignancy"] == "Normal B cell"]["score"].values
    ndlbcl = sub[sub["malignancy"] == "Normal B cell (DLBCL)"]["score"].values

    _, p1 = stats.mannwhitneyu(mal_s, norm_s, alternative="two-sided")
    _, p2 = stats.mannwhitneyu(mal_s, ndlbcl, alternative="two-sided")

    y_max = sub["score"].quantile(
        0.995
    )  # [PARAM] 括号起始高度基准分位数，调大可避免括号压点
    y_br1 = y_max * 1.10  # [PARAM] 第一条括号高度倍数（Malignant vs Normal Tonsil）
    y_br2 = y_max * 1.26  # [PARAM] 第二条括号高度倍数（Malignant vs Normal DLBCL）

    add_bracket(ax, 0, 2, y_br1, sig_stars(p1), color="black")
    add_bracket(ax, 1, 2, y_br2, sig_stars(p2), color="#444444")

    ax.set_ylim(0, y_br2 * 1.20)  # [PARAM] Y 轴上限倍数，留出括号空间
    ax.set_title(pathway_labels[pw], fontsize=14, fontweight="bold", pad=4)
    ax.set_xlabel("")
    ax.set_ylabel("AUCell Score", fontsize=10)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(
        ["Normal\n(Tonsil)", "Normal\n(DLBCL)", "Malignant"],
        fontsize=12,  # [PARAM] X 轴刻度字号
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

# 底部公共图例
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# 用一个透明的 patch 当 title 占位符
title_handle = mpatches.Patch(
    color='none',
    label='Cell Group:'    # 这就是你的 title 文字
)

fig.legend(
    handles=[title_handle] + handles,   # title 排第一个
    loc="lower center",
    ncol=len(handles) + 1,              # 列数 = 条目数 + 1（title 占一格）
    fontsize=14,
    frameon=False,
    bbox_to_anchor=(0.5, 0),
    title=None,                         # 关掉原来的 title
)


save_fig(fig, f"{OUT_CNV}/metabolic_violin2")


# ============================================================
# 图 7: metabolic_heatmap — 19 个代谢通路 heatmap（原始 + Z-score）
# ============================================================
print("=== Figure 7: Metabolic Pathway Heatmap ===")

import matplotlib.gridspec as gridspec
from scipy.stats import zscore as scipy_zscore

# [PARAM] 所有通路的可读名称映射
pw_labels = {
    "Glycolysis_Gluconeogenesis": "Glycolysis / Gluconeogenesis",
    "TCA_Cycle": "TCA Cycle",
    "Oxidative_Phosphorylation": "Oxidative Phosphorylation",
    "Glutathione_Metabolism": "Glutathione Metabolism",
    "Cysteine_Methionine_Metabolism": "Cysteine & Methionine Metabolism",
    "Glutamate_Metabolism": "Glutamate Metabolism",
    "Pentose_Phosphate_Pathway": "Pentose Phosphate Pathway",
    "Fatty_Acid_Degradation": "Fatty Acid Degradation",
    "Purine_Metabolism": "Purine Metabolism",
    "Central_Carbon_Metabolism_Cancer": "Central Carbon Metab. (Cancer)",
    "HIF1_Signaling": "HIF-1 Signaling",
    "mTOR_Signaling": "mTOR Signaling",
    "Hallmark_Hypoxia": "Hallmark: Hypoxia",
    "Hallmark_Cholesterol_Homeostasis": "Hallmark: Cholesterol Homeostasis",
    "Hallmark_mTORC1_Signaling": "Hallmark: mTORC1 Signaling",
    "Hallmark_Fatty_Acid_Metabolism": "Hallmark: Fatty Acid Metabolism",
    "Hallmark_Oxidative_Phosphorylation": "Hallmark: Oxidative Phosphorylation",
    "Hallmark_Glycolysis": "Hallmark: Glycolysis",
    "Hallmark_Reactive_Oxygen_Species_Pathway": "Hallmark: ROS Pathway",
}

# 取 aucell_df 中实际存在的通路列
pathway_cols_all = [c for c in pw_labels.keys() if c in aucell_df.columns]

# 计算每组均值
group_means = scores_df.groupby("malignancy")[pathway_cols_all].mean()
group_means = group_means.reindex([g for g in group_order if g in group_means.index])

# Z-score（对每个通路跨组标准化）
heatmap_data = group_means.copy()
heatmap_z = heatmap_data.apply(scipy_zscore, axis=0)

# 按 Malignant/Normal 倍数排序通路
fc_order = (
    group_means.loc["Malignant B cell"] / (group_means.loc["Normal B cell"] + 1e-10)
).sort_values(ascending=False)
sorted_pws = [c for c in fc_order.index]  # 原始列名顺序

# 重命名为可读标签
heatmap_data = heatmap_data[sorted_pws].rename(columns=pw_labels)
heatmap_z = heatmap_z[sorted_pws].rename(columns=pw_labels)

# X 轴标签（组名）
x_labels = ["Normal\n(Tonsil)", "Normal\n(DLBCL)", "Malignant"]

# [PARAM] 整体图尺寸
fig = plt.figure(figsize=(14, 9))
gs = gridspec.GridSpec(
    1,
    2,
    figure=fig,
    width_ratios=[1, 1],  # [PARAM] 左右子图宽度比
    wspace=0.35,  # [PARAM] 左右子图间距
)

# --- 左图: 原始 AUCell 均值 ---
ax1 = fig.add_subplot(gs[0])
sns.heatmap(
    heatmap_data.T,
    ax=ax1,
    cmap="YlOrRd",  # [PARAM] 颜色映射，可换 'Reds'、'plasma' 等
    annot=True,  # [PARAM] True=在格子里显示数值
    fmt=".3f",  # [PARAM] 数值格式，'.3f'=三位小数
    annot_kws={"size": 11},  # [PARAM] 数值字号
    linewidths=0.5,  # [PARAM] 格子间线宽，0=无线
    linecolor="white",
    cbar_kws={
        "label": "Mean AUCell Score",
        "shrink": 0.6,  # [PARAM] colorbar 相对高度
    },
    xticklabels=x_labels,
)
ax1.set_title("Mean AUCell Score per Cell Group", fontsize=13, fontweight="bold")
ax1.set_xlabel("")
ax1.set_ylabel("")
ax1.tick_params(axis="y", labelsize=12)
ax1.tick_params(axis="x", labelsize=12)

# --- 右图: Z-score 相对富集 ---
ax2 = fig.add_subplot(gs[1])
sns.heatmap(
    heatmap_z.T,
    ax=ax2,
    cmap="RdBu_r",  # [PARAM] 颜色映射，双向分歧色板；可换 'coolwarm'、'bwr'
    center=0,  # [PARAM] 颜色中心值（0=白色对应 z=0）
    vmin=-2,  # [PARAM] 颜色下限（z-score）
    vmax=2,  # [PARAM] 颜色上限（z-score）
    annot=True,
    fmt=".2f",  # [PARAM] Z-score 数值格式
    annot_kws={"size": 11},
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Z-score", "shrink": 0.6},
    xticklabels=x_labels,
)
ax2.set_title(
    "Z-scored AUCell Score(relative enrichment)", fontsize=13, fontweight="bold"
)
ax2.set_xlabel("")
ax2.set_ylabel("")
ax2.set_yticklabels([])   # 清空标签文字
# 或者连刻度线也一起去掉
ax2.tick_params(axis="y", left=False, labelleft=False)

ax2.tick_params(axis="x", labelsize=12)

fig.suptitle(
    "Metabolic Pathway Activity HeatmapDLBCL Malignant vs Normal B Cells (19 Pathways)",
    fontsize=13,
    fontweight="bold",
    y=1.01,  # [PARAM] 总标题垂直位置
)

save_fig(fig, f"{OUT_CNV}/metabolic_heatmap2")

print("All figures saved successfully!")

# ============================================================
# GSE182434_Bcell_inferCNV 本地绘图脚本
# 加载 Agent 保存的 DEG CSV，复现 Volcano plot
# 同时提供 Mono/Mac CellTypist 亚型 UMAP（如数据可用）
# 所有 "# [PARAM]" 注释标注可调参数
# ============================================================

# ============================================================
# GSE182434 Volcano Plot — 最终修正版
# 修复 1: 不转换 log2FC（已是真正 log2 单位）
# 修复 2: Y 轴用 -log10(pval)，显著性分类用 padj
# ============================================================

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # [PARAM] 改为 'TkAgg' 可弹窗预览
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")


# ============================================================
# 0. 路径配置
# ============================================================
BASE_RESULTS = translate("D:/bulk-download")

DEG_FULL = f"{BASE_RESULTS}/deg_malignant_vs_normal_bcell_full.csv"
DEG_SIG = f"{BASE_RESULTS}/deg_malignant_vs_normal_bcell_significant.csv"
OUT_DIR = f"{BASE_RESULTS}/GSE182434/bcell_deg"
os.makedirs(OUT_DIR, exist_ok=True)

SAVE_FORMATS = ["png", "svg"]  # [PARAM] 输出格式
DPI = 300  # [PARAM] 分辨率


def save_fig(fig, path_no_ext):
    for fmt in SAVE_FORMATS:
        fig.savefig(f"{path_no_ext}.{fmt}", dpi=DPI, bbox_inches="tight", format=fmt)
    plt.close(fig)
    print(f"  ✓ Saved: {path_no_ext}.*")


# ============================================================
# 1. 加载数据
# ============================================================
print("Loading DEG table ...")
deg = pd.read_csv(DEG_FULL)
print(f"  Shape: {deg.shape}")

# ── log2FC：最大值已超过 10，确认是真正 log2FC，不需要转换 ──────────────
print(f"  log2FC range: {deg['log2FC'].min():.3f} ~ {deg['log2FC'].max():.3f}")

# ── Y 轴：用原始 pval 做 -log10，padj 仅用于显著性分类 ─────────────────
# padj 列大量为 1.0（不显著基因被 scanpy 填充），直接 -log10 会全堆在 Y=0
# 改用 pval 列，保留真实的统计信号强度用于可视化
deg["pval_clip"] = deg["pval"].clip(lower=1e-300)
deg["neg_log10_pval"] = -np.log10(deg["pval_clip"])

print(
    f"  neg_log10_pval range: "
    f"{deg['neg_log10_pval'].min():.2f} ~ {deg['neg_log10_pval'].max():.2f}"
)

# ── 重新计算 direction（用 padj 做分类判断，log2FC 不变）──────────────
# [PARAM] 显著性截断值
PADJ_CUTOFF = 0.05  # [PARAM] FDR 阈值（用 padj 列判断，padj=0 的基因肯定显著）
LFC_CUTOFF = 1  # [PARAM] |log2FC| 阈值

# padj=0.0 的行是极显著基因（p 值下溢为 0），视为显著
sig_mask = deg["padj"] < PADJ_CUTOFF  # 0.0 < 0.05 → True，1.0 < 0.05 → False

deg["direction"] = np.where(
    sig_mask & (deg["log2FC"] >= LFC_CUTOFF),
    "Up in Malignant",
    np.where(
        sig_mask & (deg["log2FC"] <= -LFC_CUTOFF), "Up in Normal", "Not significant"
    ),
)

n_up_mal = (deg["direction"] == "Up in Malignant").sum()
n_up_norm = (deg["direction"] == "Up in Normal").sum()
n_ns = (deg["direction"] == "Not significant").sum()
print(f"  Up in Malignant: {n_up_mal}")
print(f"  Up in Normal:    {n_up_norm}")
print(f"  Not significant: {n_ns}")


# ============================================================
# 2. 火山图（X = log2FC，Y = -log10(pval)）
# ============================================================
print("\n=== Volcano Plot ===")

# [PARAM] 三类点颜色
volcano_colors = {
    "Up in Malignant": "#d62728",  # 红
    "Up in Normal": "#1f77b4",  # 蓝
    "Not significant": "#cccccc",  # 灰
}

# [PARAM] 标注 top N 基因（按 log2FC 绝对值排序）
N_LABEL_UP = 20  # [PARAM] 标注上调基因数
N_LABEL_DOWN = 15  # [PARAM] 标注下调基因数

deg_sig_only = deg[deg["direction"] != "Not significant"].copy()
label_up = deg_sig_only[deg_sig_only["direction"] == "Up in Malignant"].nlargest(
    N_LABEL_UP, "log2FC"
)
label_down = deg_sig_only[deg_sig_only["direction"] == "Up in Normal"].nsmallest(
    N_LABEL_DOWN, "log2FC"
)
label_genes = pd.concat([label_up, label_down])

# [PARAM] 图尺寸
fig, ax = plt.subplots(figsize=(11, 8))

# 分层绘制（灰色先画，彩色压在上层）
for direction in ["Not significant", "Up in Normal", "Up in Malignant"]:
    mask = deg["direction"] == direction
    ax.scatter(
        deg.loc[mask, "log2FC"],
        deg.loc[mask, "neg_log10_pval"],
        c=volcano_colors[direction],
        s=6 if direction == "Not significant" else 10,  # [PARAM] 点大小
        alpha=0.25 if direction == "Not significant" else 0.7,  # [PARAM] 透明度
        linewidths=0,
        rasterized=True,
        label=f"{direction} (n={mask.sum()})",
    )

# Y 轴截断线（对应 padj=0.05 时 pval 的大致位置，仅供参考）
# 因为 pval 和 padj 不是线性关系，这里用一条视觉参考线
# [PARAM] 手动设置参考线高度，None=不画
HLINE_Y = None  # 例：设为 10 表示在 -log10(pval)=10 处画水平线
if HLINE_Y:
    ax.axhline(
        HLINE_Y,
        color="grey",
        linewidth=1.0,
        linestyle="--",
        alpha=0.7,
        label=f"-log10(pval) = {HLINE_Y}",
    )

# 垂直截断线
ax.axvline(LFC_CUTOFF, color="grey", linewidth=1.0, linestyle="--", alpha=0.7)
ax.axvline(-LFC_CUTOFF, color="grey", linewidth=1.0, linestyle="--", alpha=0.7)
ax.axvline(0, color="black", linewidth=0.6, alpha=0.4)

# 基因名标注
for _, row in label_genes.iterrows():
    x_off = 5 if row["log2FC"] > 0 else -5  # [PARAM] 标注偏移方向
    ax.annotate(
        row["gene"],
        xy=(row["log2FC"], row["neg_log10_pval"]),
        xytext=(x_off, 2),  # [PARAM] 文字偏移量 (x, y) 像素
        textcoords="offset points",
        fontsize=6.5,  # [PARAM] 基因名字号
        color="#222222",
        arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.5),
    )

ax.set_xlabel("log$_2$ Fold Change  (Malignant / Normal B cells)", fontsize=12)
ax.set_ylabel("$-$log$_{10}$(p-value)", fontsize=12)
ax.set_title(
    "Differential Gene Expression\nMalignant vs Normal B Cells  (Wilcoxon, FDR < 0.05, |log₂FC| ≥ 0.5)",
    fontsize=13,
    fontweight="bold",
)

# [PARAM] 手动轴范围（None=自动）
X_LIM = None  # 例：(-15, 30)
Y_LIM = None  # 例：(0, 320)
if X_LIM:
    ax.set_xlim(X_LIM)
if Y_LIM:
    ax.set_ylim(Y_LIM)

ax.legend(markerscale=2, frameon=False, fontsize=9, loc="upper left")
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

save_fig(fig, f"{OUT_DIR}/volcano_malignant_vs_normal_bcell2")


# ============================================================
# 3. Top DEG 水平柱状图
# ============================================================
print("=== Top DEG Bar Chart ===")

N_TOP = 20  # [PARAM] 每个方向展示 top N 基因

top_up = deg_sig_only[deg_sig_only["direction"] == "Up in Malignant"].nlargest(
    N_TOP, "log2FC"
)
top_down = deg_sig_only[deg_sig_only["direction"] == "Up in Normal"].nsmallest(
    N_TOP, "log2FC"
)
combined = pd.concat([top_up, top_down]).sort_values("log2FC", ascending=True)

bar_colors = [
    volcano_colors["Up in Malignant"] if v > 0 else volcano_colors["Up in Normal"]
    for v in combined["log2FC"]
]

fig, ax = plt.subplots(figsize=(8, 10))  # [PARAM] 图尺寸

ax.barh(
    combined["gene"],
    combined["log2FC"],
    color=bar_colors,
    height=0.7,  # [PARAM] 柱高
    linewidth=0,
    alpha=0.85,
)
ax.axvline(0, color="black", linewidth=0.8)
ax.axvline(LFC_CUTOFF, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)
ax.axvline(-LFC_CUTOFF, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)

ax.set_xlabel("log$_2$ Fold Change", fontsize=11)
ax.set_title(
    f"Top {N_TOP} Up/Down-regulated Genes\nMalignant vs Normal B Cells",
    fontsize=12,
    fontweight="bold",
)
handles = [
    mpatches.Patch(color=volcano_colors["Up in Malignant"], label="Up in Malignant"),
    mpatches.Patch(color=volcano_colors["Up in Normal"], label="Up in Normal"),
]
ax.legend(handles=handles, frameon=False, fontsize=9, loc="lower right")
ax.tick_params(axis="y", labelsize=8.5)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

save_fig(fig, f"{OUT_DIR}/top_deg_barplot")

print("\nAll done!")

# ============================================================
# GSE182434 inferCNV 染色体热图 — 本地重建版
# 因为 X_cnv 矩阵未单独保存，需从 adata_processed.h5ad 重建
# 重建完成后缓存为 adata_cnv.h5ad，下次直接加载
# GSE182434 inferCNV 染色体热图 — 完整最终版 v6
# ============================================================
import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import h5py
import matplotlib
matplotlib.use('Agg')       # 改为 'TkAgg' 可弹窗预览
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import infercnvpy as cnvpy
warnings.filterwarnings('ignore')
sc.settings.autoshow = False


# ============================================================
# 0. 路径配置 ← 只需修改这里
# ============================================================
BASE         = translate("D:/bulk-download")
H5AD_MAIN    = translate(f"{BASE}/GSE182434/adata_processed.h5ad")
H5AD_CNV     = f"{BASE}/GSE182434/cnv/adata_cnv.h5ad"
CNV_META_CSV = f"{BASE}/GSE182434/cnv/malignancy_classification.csv"
OUT_DIR      = f"{BASE}/GSE182434/cnv"
os.makedirs(OUT_DIR, exist_ok=True)

SAVE_FORMATS = ['png', 'svg']
DPI          = 300
WINDOW_SIZE  = 200
STEP         = 10

def save_fig(fig, path_no_ext):
    for fmt in SAVE_FORMATS:
        fig.savefig(f"{path_no_ext}.{fmt}", dpi=DPI,
                    bbox_inches='tight', format=fmt)
    plt.close(fig)
    print(f"  ✓ Saved: {path_no_ext}.*")


# ============================================================
# 辅助函数
# ============================================================
def is_h5py_object(v):
    return 'h5py' in str(type(v).__module__)

def fix_chr_pos(uns):
    """在 deep_clean_uns 之前调用，把 chr_pos 转成纯 Python 类型"""
    if 'cnv' in uns and 'chr_pos' in uns.get('cnv', {}):
        uns['cnv']['chr_pos'] = {
            str(k): int(v)
            for k, v in uns['cnv']['chr_pos'].items()
        }
        print(f"  ✓ chr_pos fixed: {len(uns['cnv']['chr_pos'])} chromosomes")
    return uns

def deep_clean_uns(d, path='uns'):
    cleaned = {}
    for k, v in d.items():
        if is_h5py_object(v):
            print(f"    Removed {path}/{k} ({type(v).__name__})")
            continue
        if isinstance(v, dict):
            cleaned[k] = deep_clean_uns(v, f"{path}/{k}")
        else:
            cleaned[k] = v
    return cleaned

def compute_cnv_score(adata):
    X = adata.obsm['X_cnv']
    if sp.issparse(X):
        return np.array(np.abs(X).mean(axis=1)).flatten()
    return np.abs(X).mean(axis=1)

def assign_malignancy(adata, cnv_meta_csv, quantile=0.99):
    if os.path.exists(cnv_meta_csv):
        print("  Loading malignancy labels from CSV ...")
        meta = pd.read_csv(cnv_meta_csv, index_col=0)
        adata.obs['malignancy'] = meta['malignancy'].reindex(adata.obs_names)
        adata.obs['malignancy'] = adata.obs['malignancy'].fillna('Normal B cell')
    else:
        print("  Classifying malignancy by CNV score threshold ...")
        tonsil_scores = adata.obs.loc[
            adata.obs['cnv_group'] == 'Tonsil_ref', 'cnv_score']
        threshold = tonsil_scores.quantile(quantile)
        print(f"  Threshold (Tonsil {int(quantile*100)}th pct): {threshold:.4f}")
        dlbcl_mask = adata.obs['cnv_group'] == 'DLBCL_B'
        adata.obs['malignancy'] = 'Normal B cell'
        adata.obs.loc[
            dlbcl_mask & (adata.obs['cnv_score'] > threshold),
            'malignancy'] = 'Malignant B cell'
        adata.obs.loc[
            dlbcl_mask & (adata.obs['cnv_score'] <= threshold),
            'malignancy'] = 'Normal B cell (DLBCL)'
    print(f"  Malignancy: {adata.obs['malignancy'].value_counts().to_dict()}")


# ============================================================
# 1. 加载或重建 adata_cnv
# ============================================================
if os.path.exists(H5AD_CNV):
    print(f"Loading cached adata_cnv ...")
    adata_cnv = sc.read_h5ad(H5AD_CNV)
    print(f"  shape:     {adata_cnv.shape}")
    print(f"  obsm keys: {list(adata_cnv.obsm.keys())}")
    print(f"  uns keys:  {list(adata_cnv.uns.keys())}")

else:
    print("Cache not found. Running full pipeline ...")

    print("  Loading adata_processed.h5ad ...")
    adata = sc.read_h5ad(H5AD_MAIN)

    b_types = [
        'Naive B cells', 'GC B cells', 'Proliferative GC B cells',
        'Memory B cells', 'Age-associated B cells', 'B cells (other)', 'Plasma cells'
    ]
    adata_b = adata[adata.obs['cell_type'].isin(b_types)].copy()
    adata_b.X = adata_b.layers['counts'].copy()
    print(f"  B cells: {adata_b.shape[0]}")

    adata_b.obs['cnv_group'] = 'DLBCL_B'
    adata_b.obs.loc[adata_b.obs['Tissue'] == 'Tonsil', 'cnv_group'] = 'Tonsil_ref'
    print(f"  cnv_group: {adata_b.obs['cnv_group'].value_counts().to_dict()}")

    print("  Fetching gene positions from Ensembl biomart ...")
    annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name",
         "chromosome_name", "start_position", "end_position"],
        use_cache=True
    )
    valid_chroms = [str(i) for i in range(1, 23)] + ['X']
    annot_filt   = (annot[annot['chromosome_name'].isin(valid_chroms)]
                    .drop_duplicates('external_gene_name')
                    .set_index('external_gene_name'))
    adata_b.var['chromosome'] = annot_filt['chromosome_name'].reindex(adata_b.var_names)
    adata_b.var['start']      = annot_filt['start_position'].reindex(adata_b.var_names)
    adata_b.var['end']        = annot_filt['end_position'].reindex(adata_b.var_names)
    adata_b.var['chromosome'] = adata_b.var['chromosome'].apply(
        lambda x: f'chr{x}' if pd.notna(x) and not str(x).startswith('chr') else x
    )

    adata_cnv = adata_b[:, adata_b.var['chromosome'].notna()].copy()
    if sp.issparse(adata_cnv.X):
        adata_cnv.X = adata_cnv.X.astype(np.float32)
    print(f"  After filtering: {adata_cnv.shape[0]} cells × {adata_cnv.shape[1]} genes")

    print("  Running inferCNV (may take 10~30 min) ...")
    cnvpy.tl.infercnv(
        adata_cnv,
        reference_key='cnv_group',
        reference_cat=['Tonsil_ref'],  # ← 必须是 list！
        window_size=WINDOW_SIZE,
        step=STEP,
        lfc_clip=3.0,
        dynamic_threshold=1.5,
        chunksize=5000,
        n_jobs=1,
    )
    print(f"  ✓ inferCNV done.")

    # 固化 chr_pos 并清理
    adata_cnv.uns = fix_chr_pos(adata_cnv.uns)
    adata_cnv.uns = deep_clean_uns(adata_cnv.uns)

    # 计算 CNV 分数并分配恶性标签
    adata_cnv.obs['cnv_score'] = compute_cnv_score(adata_cnv)
    assign_malignancy(adata_cnv, CNV_META_CSV)

    print(f"  Saving cache ...")
    adata_cnv.write_h5ad(H5AD_CNV, compression='gzip')
    print("  ✓ Cache saved")


# ============================================================
# 2. 补全字段 + 强制 malignancy 为有序 Categorical（关键！）
# ============================================================
if 'cnv' not in adata_cnv.uns or 'chr_pos' not in adata_cnv.uns.get('cnv', {}):
    raise RuntimeError(
        f"chr_pos missing! Please delete {H5AD_CNV} and re-run."
    )
if 'cnv_score' not in adata_cnv.obs.columns:
    adata_cnv.obs['cnv_score'] = compute_cnv_score(adata_cnv)
if 'malignancy' not in adata_cnv.obs.columns:
    assign_malignancy(adata_cnv, CNV_META_CSV)

# ✅ 关键修复：将 malignancy 转为有序分类变量
DISPLAY_ORDER = ['Normal B cell', 'Normal B cell (DLBCL)', 'Malignant B cell']
present_cats = [cat for cat in DISPLAY_ORDER if cat in adata_cnv.obs['malignancy'].values]
adata_cnv.obs['malignancy'] = pd.Categorical(
    adata_cnv.obs['malignancy'],
    categories=present_cats,
    ordered=True
)

print(f"\nReady:")
print(f"  cells:      {adata_cnv.n_obs}")
print(f"  X_cnv:      {adata_cnv.obsm['X_cnv'].shape}")
print(f"  chr_pos:    {list(adata_cnv.uns['cnv']['chr_pos'].keys())}")
print(f"  malignancy: {adata_cnv.obs['malignancy'].value_counts().to_dict()}")
print(f"  categories: {adata_cnv.obs['malignancy'].cat.categories.tolist()}")


# ============================================================
# 图 1: CNV Score 分布图
# ============================================================
print("\n=== Figure 1: CNV Score Distribution ===")

dlbcl_scores  = adata_cnv.obs.loc[adata_cnv.obs['cnv_group'] == 'DLBCL_B',   'cnv_score']
tonsil_scores = adata_cnv.obs.loc[adata_cnv.obs['cnv_group'] == 'Tonsil_ref', 'cnv_score']
threshold_val = tonsil_scores.quantile(0.99)

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(tonsil_scores, bins=50, alpha=0.6, color='#2ca02c', label='Tonsil ref')
ax.hist(dlbcl_scores,  bins=50, alpha=0.6, color='#1f77b4', label='DLBCL B')
ax.axvline(threshold_val, color='red', linestyle='--', linewidth=1.2,
           label=f'Threshold (99th pct = {threshold_val:.4f})')
ax.set_xlabel('CNV Score', fontsize=12)
ax.set_ylabel('Cell Count', fontsize=12)
ax.set_title('CNV Score Distribution', fontsize=12, fontweight='bold')
ax.legend(frameon=False, fontsize=13)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
plt.tight_layout()
save_fig(fig, f"{OUT_DIR}/cnv_score_distribution")


# ============================================================
# 图 2: 染色体热图 — 正确顺序（无需手动排序！）
# ============================================================
print("\n=== Figure 2: Chromosome Heatmap ===")

GROUP_COLORS = {
    'Normal B cell':         '#2ca02c',
    'Normal B cell (DLBCL)': '#FF9400',
    'Malignant B cell':      '#d62728',
}

# ✅ 颜色顺序必须与 .cat.categories 一致
actual_categories = adata_cnv.obs['malignancy'].cat.categories.tolist()
adata_cnv.uns['malignancy_colors'] = [GROUP_COLORS[cat] for cat in actual_categories]

# 直接绘图（infercnvpy 会自动按 category 顺序排列行）
axes_dict = cnvpy.pl.chromosome_heatmap(
    adata_cnv,
    groupby='malignancy',
    figsize=(9.8, 6.2),
    show=False,
    cmap='bwr',
)


fig = list(axes_dict.values())[0].get_figure() if axes_dict else plt.gcf()

# ── 事后统一调整所有字体 ──────────────────────────────────────────────────
for ax in fig.axes:
    # x 轴（染色体号）
    ax.tick_params(axis='x', labelsize=10)
    for tick in ax.get_xticklabels():
        tick.set_fontsize(10)
        tick.set_rotation(45)
        tick.set_ha('right')

    # y 轴（分组标签）—— 竖排 + 放大
    ax.tick_params(axis='y', labelsize=13)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(10)
        tick.set_rotation(0)         # 横排，短标签不再互相叠印
        tick.set_va('center')
        # 组名过长会伸出画布左缘（QA 实测），统一缩短
        _short = {"Normal B cell (Tonsil)": "Normal (Tonsil)",
                  "Normal B cell (DLBCL)": "Normal (DLBCL)",
                  "Malignant B cell (DLBCL)": "Malignant",
                  "Malignant B cell": "Malignant",
                  "Normal B cell": "Normal"}.get(tick.get_text())
        if _short:
            tick.set_text(_short)
    # set_text 对固定 locator 无效，需整体替换一次
    for ax in fig.axes:
        _yl = ax.get_yticklabels()
        _txt = [t2.get_text() for t2 in _yl]
        if any("(Tonsil)" in s or "(DLBCL)" in s for s in _txt):
            ax.set_yticklabels([{"Normal B cell (Tonsil)": "Normal (Tonsil)",
                                 "Normal B cell (DLBCL)": "Normal (DLBCL)",
                                 "Malignant B cell (DLBCL)": "Malignant",
                                 "Malignant B cell": "Malignant",
                                 "Normal B cell": "Normal"}.get(s, s) for s in _txt],
                               fontsize=10)

    # colorbar 的刻度（colorbar 本质上也是一个 ax）
    if hasattr(ax, 'collections') and len(ax.collections) == 0:
        ax.tick_params(labelsize=10)

fig.suptitle(
    'Inferred Copy Number Variations Along Chromosomes\n'
    'DLBCL B Cells vs Tonsil Reference  |  inferCNVpy (window=200 genes)',
    fontsize=13, fontweight='bold', y=1.01
)

# 图例按 DISPLAY_ORDER 从左到右
legend_handles = [
    mpatches.Patch(color=GROUP_COLORS[g], label=g)
    for g in actual_categories
]
# 图例独占底部条带，避免压住热图本体（QA/作者反馈）
# 右侧 colorbar + 顶部 chr 标签都要留出空间；chr20-22 拥挤 -> 45 度右对齐
for _ax in fig.axes:
    _xt = _ax.get_xticklabels()
    if any(s.get_text().startswith("chr") for s in _xt):
        # 去掉 chr 前缀（纯数字更窄，chr20-22 不再互撞），并加轴标题
        _ax.set_xticklabels([s.get_text().replace("chr", "") for s in _xt], fontsize=10)
        for _t in _xt:
            _t.set_rotation(90)
            _t.set_ha("center")
        _ax.set_xlabel("Chromosome", fontsize=10)
fig.subplots_adjust(right=0.82, top=0.86, bottom=0.13)
fig.legend(
    handles=legend_handles,
    loc='lower center',
    fontsize=10,
    frameon=True,
    title='Cell Group',
    title_fontsize=12,
    bbox_to_anchor=(0.5, 0.005),
    ncol=len(actual_categories),    # 所有条目排成一行    
)

save_fig(fig, f"{OUT_DIR}/chromosome_heatmap")
# svglib 转 PDF 会丢失文字旋转（chr 标签变横排互撞），此面板直接存原生 PDF
fig.savefig(f"{OUT_DIR}/chromosome_heatmap.pdf")
print("All done!")

# ============================================================
# GSE182434 Mono/Mac 亚型分析 — 完整本地绘图版
# 包含新增的 mac_subtype UMAP 图
# ============================================================

# ============================================================
# GSE182434 — Monocytes / Macrophages 最终本地分析脚本
# ============================================================

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import warnings
warnings.filterwarnings("ignore")

# ----------------------------
# 路径设置（统一本地路径）
# ----------------------------
BASE_DIR = translate("D:/bulk-download/GSE182434")
FIG_DIR  = f"{BASE_DIR}/mono_mac"
MARKER_XLSX = translate("D:/bulk-download/Macro_mono_cellmarker.xlsx")

import os
os.makedirs(FIG_DIR, exist_ok=True)

# ----------------------------
# 1. 读取主 AnnData
# ----------------------------
adata = sc.read_h5ad(translate(f"{BASE_DIR}/adata_processed.h5ad"))

# 只保留 Mono / Mac
adata_mac = adata[adata.obs["cell_type"] == "Monocytes/Macrophages"].copy()
print(f"✅ Mono/Mac cells: {adata_mac.n_obs}")

# ============================================================
# ✅ CRITICAL FIX: 清理 h5py.Empty（必须在 HVG 前）
# ============================================================
def clean_uns_empty(adata):
    """递归删除 adata.uns 中所有 h5py.Empty 对象"""
    def _clean(obj):
        if isinstance(obj, dict):
            new = {}
            for k, v in obj.items():
                if isinstance(v, h5py._hl.base.Empty):
                    print(f"  Removed uns['{k}'] (h5py.Empty)")
                    continue
                new[k] = _clean(v)
            return new
        return obj
    adata.uns = _clean(adata.uns)

clean_uns_empty(adata_mac)

# ============================================================
# 2. 重新 PCA → neighbors → UMAP → Leiden
# ============================================================
sc.pp.highly_variable_genes(
    adata_mac,
    n_top_genes=2000,
    flavor="seurat",
    subset=False
)

sc.tl.pca(adata_mac, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata_mac, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata_mac)
sc.tl.leiden(adata_mac, resolution=0.4, key_added="leiden_mac")

print("Leiden clusters:")
print(adata_mac.obs["leiden_mac"].value_counts())

# ============================================================
# 3. UMAP：Leiden + LRP1
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(
    adata_mac,
    color="leiden_mac",
    ax=axes[0],
    show=False,
    frameon=False,
    legend_loc="on data",
    title="Mono/Mac Leiden clusters",
    size=40
)

sc.pl.umap(
    adata_mac,
    color="LRP1",
    ax=axes[1],
    show=False,
    cmap="RdYlBu_r",
    frameon=False,
    title="LRP1 expression",
    size=40,
    vmin=0
)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/mac_umap_leiden_lrp12.png", dpi=300)
plt.close()

# ============================================================
# 4. M1 / M2 marker violin
# ============================================================
m1_markers = ["CD86", "CD80", "IL1B", "TNF", "CXCL10", "HLA-DRA"]
m2_markers = ["CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CCL18"]

genes = [g for g in m1_markers + m2_markers if g in adata_mac.var_names]

fig, axes = plt.subplots(4, 3, figsize=(16, 18))
axes = axes.flatten()

for i, g in enumerate(genes):
    sc.pl.violin(
        adata_mac,
        keys=g,
        groupby="leiden_mac",
        ax=axes[i],
        show=False,
        stripplot=False
    )
    color = "#C0392B" if g in m1_markers else "#2471A3"
    label = "M1" if g in m1_markers else "M2"
    axes[i].set_title(f"{g} [{label}]", color=color, fontweight="bold")

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/mac_m1_m2_violin2.png", dpi=300)
plt.close()

# ============================================================
# 5. 读取巨噬细胞亚型 marker
# ============================================================
df_markers = pd.read_excel(MARKER_XLSX, sheet_name=0)

subtypes = ["Mono", "DC_1", "LA_TAM", "IFN_TAM", "DC_2"]
marker_dict = {
    ct: df_markers[df_markers["Celltype"] == ct]["gene"].tolist()
    for ct in subtypes
}

# ============================================================
# 6. 每个细胞做 marker score
# ============================================================
for ct, genes in marker_dict.items():
    genes = [g for g in genes if g in adata_mac.var_names]
    sc.tl.score_genes(
        adata_mac,
        gene_list=genes,
        score_name=f"score_{ct}",
        use_raw=False
    )

score_cols = [f"score_{ct}" for ct in subtypes]
score_mat = adata_mac.obs[score_cols].values
best_idx = np.argmax(score_mat, axis=1)

adata_mac.obs["mac_subtype"] = pd.Categorical(
    [subtypes[i] for i in best_idx],
    categories=subtypes
)

print("Subtype distribution:")
print(adata_mac.obs["mac_subtype"].value_counts())

# ============================================================
# 7. UMAP：Macrophage subtypes
# ============================================================
palette = {
    "Mono": "#1ABC9C",
    "DC_1": "#9B59B6",
    "LA_TAM": "#E67E22",
    "IFN_TAM": "#E74C3C",
    "DC_2": "#3498DB",
}

sc.pl.umap(
    adata_mac,
    color="mac_subtype",
    palette=palette,
    frameon=False,
    title="Macrophage subtypes",
    size=40,
    save=False,
    show=False
)

plt.savefig(f"{FIG_DIR}/mac_subtype_umap2.png", dpi=300, bbox_inches="tight")
plt.close()

# ============================================================
# 8. Dotplot（Top markers）
# ============================================================
top_n = 5
plot_genes = []

for ct in subtypes:
    top = (
        df_markers[df_markers["Celltype"] == ct]
        .sort_values("avg_log2FC", ascending=False)
        .head(top_n)["gene"]
        .tolist()
    )
    plot_genes.extend(top)

plot_genes = list(dict.fromkeys(plot_genes))
plot_genes = [g for g in plot_genes if g in adata_mac.var_names]

plot_genes = ['FCN1',
 'S100A9',
 'S100A8',
 'S100A4',
 'APOBEC3A',
 'LTB',
 'CLEC10A',
 'CD1C',
 'JAML',
 'CD1E',
 'APOC1',
 'CCL18',
 'APOE',
 'ACP5',
 'CTSD',
 'MT1H',
 'MT1G',
 'CCL8',
 'CCL2',
 'MT1X',
 'DNASE1L3',
 'RGCC',
 'CST3',
 'CLEC9A',
 'SNX3']

dp = sc.pl.dotplot(
    adata_mac,
    var_names=plot_genes,
    groupby="mac_subtype",
    standard_scale="var",
    cmap="YlOrRd",
    dot_max=0.8,
    show=False
)

# ✅ 关键修正：dp 本身就是 ax_dict
ax = dp["mainplot_ax"]

for tick in ax.get_yticklabels():
    ct = tick.get_text()
    if ct in palette:
        tick.set_color(palette[ct])
        tick.set_fontweight("bold")
plt.savefig(f"{FIG_DIR}/mac_subtype_dotplot2.png", dpi=300, bbox_inches="tight")
plt.close()


sc.pl.dotplot(
    adata_mac,
    var_names=plot_genes,
    groupby="mac_subtype",
    standard_scale="var",
    cmap="YlOrRd",
    dot_max=0.8,
    save=False,
    show=False
)

plt.savefig(f"{FIG_DIR}/mac_subtype_dotplot2.png", dpi=300, bbox_inches="tight")
plt.close()

# ============================================================
print("✅ Mono / Macrophage analysis finished successfully.")
print(f"📁 Figures saved to: {FIG_DIR}")

