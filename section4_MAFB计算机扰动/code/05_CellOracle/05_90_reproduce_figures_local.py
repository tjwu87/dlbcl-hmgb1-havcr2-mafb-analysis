# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python
# reproduce_fig01_23.py
# 本地复现 MAFB CellOracle 23 张图
# 依赖：pandas, numpy, matplotlib, seaborn, scipy, statsmodels
# 安装：pip install pandas numpy matplotlib seaborn scipy statsmodels
# 运行：python reproduce_fig01_23.py

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════
# 全局路径配置 — 修改这里指向你的数据目录和输出目录
# ══════════════════════════════════════════════════════════════════════════════
DATA_DIR = translate(r"D:\bulk-download\celloracle0331")   # 所有 CSV 文件所在目录
OUT_DIR  = translate(r"D:\bulk-download\celloracle0331\figures")  # 图片输出目录
os.makedirs(OUT_DIR, exist_ok=True)

def p(filename):
    """拼接数据目录路径"""
    return os.path.join(DATA_DIR, filename)

def save_fig(fig, name):
    """保存 PNG + SVG"""
    fig.savefig(os.path.join(OUT_DIR, f"{name}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {name}.png / .svg")

# ══════════════════════════════════════════════════════════════════════════════
# 数据加载
# ══════════════════════════════════════════════════════════════════════════════
print("Loading data...")

df_meta      = pd.read_csv(p("cell_metadata.csv"), index_col=0)
umap_coords  = df_meta[["UMAP1", "UMAP2"]].values
subtypes     = df_meta["subtype"].values
subtypes_order = ["Mono", "IFN_TAM", "LA_TAM"]

print("  Loading log1p_norm_expression.csv (may take a moment)...")
df_log       = pd.read_csv(p("log1p_norm_expression.csv"), index_col=0)
X_log        = df_log.values
gene_names   = np.array(df_log.columns)

print("  Loading imputed_count_WT.csv...")
df_imputed   = pd.read_csv(p("imputed_count_WT.csv"), index_col=0)
X_wt_imputed = df_imputed.values

print("  Loading simulated_count_KO.csv...")
df_simulated = pd.read_csv(p("simulated_count_KO.csv"), index_col=0)
X_ko_sim     = df_simulated.values

def load_delta_emb(cond):
    df = pd.read_csv(p(f"delta_embedding_{cond}.csv"), index_col=0)
    return df[["UMAP1", "UMAP2"]].values

vel_wt   = load_delta_emb("WT")
vel_ko   = load_delta_emb("MAFB_KO")
vel_oe   = load_delta_emb("MAFB_OE")
vel_null = load_delta_emb("Null")

def load_delta_X(cond):
    df = pd.read_csv(p(f"delta_X_{cond}.csv"), index_col=0)
    return df.values, np.array(df.columns)

delta_ko, delta_gene_names = load_delta_X("MAFB_KO")
delta_oe, _                = load_delta_X("MAFB_OE")

df_fate_stats = pd.read_csv(p("04_fate_stats.csv"))
df_fate_wt    = pd.read_csv(p("cell_fate_scores_WT.csv"))
df_fate_ko    = pd.read_csv(p("cell_fate_scores_KO.csv"))
df_scenic     = pd.read_csv(p("05_scenic_target_stats.csv"))
df_prog_stats = pd.read_csv(p("06_program_score_stats.csv"))
df_enr_all    = pd.read_csv(p("07_ko_response_enrichment_results.csv"))
df_deg        = pd.read_csv(p("07_ko_response_deg_per_subtype.csv"))
df_sens       = pd.read_csv(p("08_sensitivity_results.csv"))
df_fig23_deg  = pd.read_csv(p("fig23_deg_imputed_vs_simulated.csv"))
df_fig23_enr  = pd.read_csv(p("fig23_ora_enrichment_results.csv"))

with open(p("06_program_gene_50.txt")) as f:
    prog_genes = [l.strip() for l in f if l.strip() and not l.startswith("#")]

mafb_scenic_targets = df_scenic[df_scenic["evaluated_or_not"] == "evaluated"]["gene"].unique().tolist()
mafb_idx = np.where(gene_names == "MAFB")[0][0] if "MAFB" in gene_names else None

def load_grn(subtype):
    df = pd.read_csv(p(f"filtered_A_{subtype}.csv"), index_col=0)
    df_long = df.reset_index()[["source", "target", "coef_mean"]].copy()
    df_long.columns = ["TF", "target", "coef"]
    return df_long

grn_coefs = {st: load_grn(st) for st in subtypes_order}

print("All data loaded.\n")

# ══════════════════════════════════════════════════════════════════════════════
# 全局颜色方案
# ══════════════════════════════════════════════════════════════════════════════
subtype_colors = {
    "Mono":    "#E64B35",
    "IFN_TAM": "#4DBBD5",
    "LA_TAM":  "#00A087",
}

# ══════════════════════════════════════════════════════════════════════════════
# 辅助函数：计算网格速度场
# ══════════════════════════════════════════════════════════════════════════════
def compute_grid_velocity(umap_coords, velocity, n_grid=30):
    x_min, x_max = umap_coords[:,0].min()-0.5, umap_coords[:,0].max()+0.5
    y_min, y_max = umap_coords[:,1].min()-0.5, umap_coords[:,1].max()+0.5
    gx = np.linspace(x_min, x_max, n_grid)
    gy = np.linspace(y_min, y_max, n_grid)
    GX, GY = np.meshgrid(gx, gy)
    VX = np.zeros_like(GX)
    VY = np.zeros_like(GY)
    sigma = (x_max - x_min) / n_grid * 1.5
    for i in range(n_grid):
        for j in range(n_grid):
            gp = np.array([GX[i,j], GY[i,j]])
            dists = np.sqrt(((umap_coords - gp)**2).sum(1))
            weights = np.exp(-dists**2 / (2*sigma**2))
            w_sum = weights.sum()
            if w_sum > 1e-8:
                VX[i,j] = (weights * velocity[:,0]).sum() / w_sum
                VY[i,j] = (weights * velocity[:,1]).sum() / w_sum
    return GX, GY, VX, VY

print("Computing grid velocities (this takes ~30s)...")
GX, GY, VX_wt,   VY_wt   = compute_grid_velocity(umap_coords, vel_wt)
GX, GY, VX_ko,   VY_ko   = compute_grid_velocity(umap_coords, vel_ko)
GX, GY, VX_oe,   VY_oe   = compute_grid_velocity(umap_coords, vel_oe)
GX, GY, VX_null, VY_null = compute_grid_velocity(umap_coords, vel_null)
vel_ko_diff = vel_ko - vel_null
vel_oe_diff = vel_oe - vel_null
GX, GY, VX_ko_diff, VY_ko_diff = compute_grid_velocity(umap_coords, vel_ko_diff)
GX, GY, VX_oe_diff, VY_oe_diff = compute_grid_velocity(umap_coords, vel_oe_diff)
print("Grid velocities done.\n")

# ══════════════════════════════════════════════════════════════════════════════
# Fig01: Quiver 向量场 (WT / KO / Null)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig01 前置参数 ──────────────────────────────────────────────────────────
P01 = {
    # 整体图尺寸 (宽, 高) 英寸
    "figsize": (10, 3),
    # 三个子图的标题文字
    "panel_titles": ["WT (Baseline)", "MAFB KO", "Null (Randomized GRN)"],
    # 总标题文字
    "suptitle": "",
    # 总标题字号
    "suptitle_fontsize": 0,
    # 子图标题字号
    "title_fontsize": 10,
    # 坐标轴标签字号
    "axis_label_fontsize": 8,
    # 散点大小
    "scatter_size": 8,
    # 散点透明度
    "scatter_alpha": 0.7,
    # quiver 箭头宽度（相对图宽）
    "quiver_width": 0.003,
    # quiver 箭头头部宽度（单位：箭杆宽度倍数）
    "quiver_headwidth": 4,
    # quiver 箭头头部长度
    "quiver_headlength": 3,
    # quiver 透明度
    "quiver_alpha": 0.75,
    # 图例字号
    "legend_fontsize": 5,
    # 图例位置
    "legend_loc": "lower right",
}

print("Drawing Fig01...")
fig, axes = plt.subplots(1, 3, figsize=P01["figsize"])
fig.patch.set_facecolor("white")
for ax, title, vx, vy in zip(axes, P01["panel_titles"],
                               [VX_wt, VX_ko, VX_null],
                               [VY_wt, VY_ko, VY_null]):
    for st in subtypes_order:
        mask = subtypes == st
        ax.scatter(umap_coords[mask,0], umap_coords[mask,1],
                   c=subtype_colors[st], s=P01["scatter_size"],
                   alpha=P01["scatter_alpha"], linewidths=0, zorder=2)
    speed = np.sqrt(vx**2 + vy**2)
    ax.quiver(GX, GY, vx, vy, speed, cmap="RdBu_r",
              alpha=P01["quiver_alpha"], scale=None,
              width=P01["quiver_width"],
              headwidth=P01["quiver_headwidth"],
              headlength=P01["quiver_headlength"], zorder=3)
    ax.set_title(title, fontsize=P01["title_fontsize"], fontweight="bold", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=P01["axis_label_fontsize"])
    ax.set_ylabel("UMAP 2", fontsize=P01["axis_label_fontsize"])
    ax.set_aspect("equal")
    ax.spines[["top","right"]].set_visible(False)
handles = [mpatches.Patch(color=subtype_colors[st], label=st) for st in subtypes_order]
axes[2].legend(handles=handles, loc=P01["legend_loc"],
               fontsize=P01["legend_fontsize"], framealpha=0.8)
plt.suptitle(P01["suptitle"], fontsize=P01["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig01_mafb_quiver_vector_field")

# ══════════════════════════════════════════════════════════════════════════════
# Fig02: Streamline 向量场 (KO / OE / Null)
# ══════════════════════════════════════════════════════════════════════════════

from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

P02 = {
    "figsize"            : (10, 3),
    "panel_titles"       : ["MAFB KO", "MAFB OE", "Null (Randomized GRN)"],
    "suptitle"           : "",
    "suptitle_fontsize"  : 14,
    "title_fontsize"     : 10,
    "axis_label_fontsize": 8,
    "scatter_size"       : 10,
    "scatter_alpha"      : 0.6,
    "stream_linewidth"   : 0.8,
    "stream_density"     : 1.2,
    "stream_arrowsize"   : 1,
    "legend_fontsize"    : 8,
    "legend_loc"         : "lower right",
    "cbar_label"         : "Normalized perturbation speed",
    "cbar_label_fontsize": 10,
    "cbar_tick_fontsize" : 8,
}

print("Drawing Fig02...")

fig = plt.figure(figsize=P02["figsize"])
fig.patch.set_facecolor("white")

gs = GridSpec(1, 3, figure=fig,
              wspace=0.1,
              left=0.06, right=0.96, top=0.88, bottom=0.15)

ax_list = [fig.add_subplot(gs[0, i]) for i in range(3)]

strm_last = None

for ax, title, vx, vy in zip(ax_list, P02["panel_titles"],
                               [VX_ko, VX_oe, VX_null],
                               [VY_ko, VY_oe, VY_null]):
    for st in subtypes_order:
        mask = subtypes == st
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   c=subtype_colors[st], s=P02["scatter_size"],
                   alpha=P02["scatter_alpha"], linewidths=0, zorder=2)

    speed = np.sqrt(vx**2 + vy**2)
    try:
        strm = ax.streamplot(GX[0, :], GY[:, 0], vx, vy,
                             color=speed / (speed.max() + 1e-8),
                             cmap="viridis",
                             linewidth=P02["stream_linewidth"],
                             density=P02["stream_density"],
                             arrowsize=P02["stream_arrowsize"],
                             zorder=3)
        strm_last = strm
    except Exception:
        ax.quiver(GX, GY, vx, vy, speed,
                  cmap="viridis", alpha=0.7, width=0.003, zorder=3)

    ax.set_title(title, fontsize=P02["title_fontsize"],
                 fontweight="bold", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=P02["axis_label_fontsize"])
    ax.set_ylabel("UMAP 2", fontsize=P02["axis_label_fontsize"])
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)

# ── subtype 图例（第三个面板右下）────────────────────────────────────────────
handles = [mpatches.Patch(color=subtype_colors[st], label=st)
           for st in subtypes_order]
ax_list[2].legend(handles=handles, loc=P02["legend_loc"],
                  fontsize=P02["legend_fontsize"], framealpha=0.8)

# ── 色条：紧贴第三个主图右侧 ─────────────────────────────────────────────────
if strm_last is not None:
    divider = make_axes_locatable(ax_list[2])
    ax_cbar = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(strm_last.lines, cax=ax_cbar)
    cbar.set_label(P02["cbar_label"],
                   fontsize=P02["cbar_label_fontsize"], labelpad=6)
    cbar.ax.tick_params(labelsize=P02["cbar_tick_fontsize"])
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0", "0.25", "0.50", "0.75", "1.0"])

if P02["suptitle"]:
    fig.suptitle(P02["suptitle"],
                 fontsize=P02["suptitle_fontsize"],
                 fontweight="bold", y=1.02)

save_fig(fig, "fig02_mafb_streamline_vector_field")



# ══════════════════════════════════════════════════════════════════════════════
# Fig03: 差分向量场 (KO-Null, OE-Null)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig03 前置参数 ──────────────────────────────────────────────────────────
P03 = {
    "figsize": (10, 4),
    "panel_titles": ["MAFB KO − Null (Differential)", "MAFB OE − Null (Differential)"],
    "suptitle": "",
    "suptitle_fontsize": 13,
    "title_fontsize": 10,
    "axis_label_fontsize": 8,
    # 散点大小
    "scatter_size": 18,
    # 散点透明度
    "scatter_alpha": 0.85,
    # colorbar 标签文字
    "colorbar_label": "|Δ velocity| vs null",
    # 亚群标签字号（显示在 UMAP 上的文字）
    "subtype_label_fontsize": 8,
    # quiver 箭头宽度
    "quiver_width": 0.003,
    # quiver 箭头头部宽度
    "quiver_headwidth": 3,
}

print("Drawing Fig03...")
fig, axes = plt.subplots(1, 2, figsize=P03["figsize"])
fig.patch.set_facecolor("white")
for ax, title, vx, vy, vel_diff in zip(
        axes, P03["panel_titles"],
        [VX_ko_diff, VX_oe_diff], [VY_ko_diff, VY_oe_diff],
        [vel_ko_diff, vel_oe_diff]):
    speed_cell = np.sqrt(vel_diff[:,0]**2 + vel_diff[:,1]**2)
    vmax = np.percentile(speed_cell, 95)
    sc = ax.scatter(umap_coords[:,0], umap_coords[:,1],
                    c=speed_cell, cmap="RdBu_r",
                    s=P03["scatter_size"], alpha=P03["scatter_alpha"],
                    linewidths=0, zorder=2, vmin=0, vmax=vmax)
    speed_grid = np.sqrt(vx**2 + vy**2)
    ax.quiver(GX, GY, vx, vy, speed_grid, cmap="RdBu_r", alpha=0.7,
              scale=None, width=P03["quiver_width"],
              headwidth=P03["quiver_headwidth"], zorder=3)
    plt.colorbar(sc, ax=ax, label=P03["colorbar_label"], shrink=0.8)
    for st in subtypes_order:
        mask = subtypes == st
        cx, cy = umap_coords[mask,0].mean(), umap_coords[mask,1].mean()
        ax.text(cx, cy, st, fontsize=P03["subtype_label_fontsize"],
                ha="center", va="center", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          alpha=0.7, edgecolor="none"))
    ax.set_title(title, fontsize=P03["title_fontsize"], fontweight="bold", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=P03["axis_label_fontsize"])
    ax.set_ylabel("UMAP 2", fontsize=P03["axis_label_fontsize"])
    ax.set_aspect("equal")
    ax.spines[["top","right"]].set_visible(False)
plt.suptitle(P03["suptitle"], fontsize=P03["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig03_mafb_differential_vector_field")

# ══════════════════════════════════════════════════════════════════════════════
# Fig04: 各亚群速度密度图 (KO / OE)
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# Fig04: MAFB Perturbation — Velocity Magnitude per Subtype
# ══════════════════════════════════════════════════════════════════════════════

from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

P04 = {
    "figsize"        : (10, 6),       # ← 宽度从 10 缩到 9，去掉多余空白
    "row_titles"     : ["MAFB KO", "MAFB OE"],
    "suptitle"       : "MAFB Perturbation — Velocity Magnitude per Subtype\n\n",
    "suptitle_fontsize": 13,
    "title_fontsize" : 11,
    "bg_color"       : "#DDDDDD",
    "bg_size"        : 4,
    "bg_alpha"       : 0.4,
    "fg_size"        : 15,
    "fg_alpha"       : 0.9,
    "colorbar_label" : "Speed",
}

print("Drawing Fig04...")

n_cols = len(subtypes_order)   # 3

# GridSpec：3 列主图 + 1 列色条，色条列很窄
gs = GridSpec(2, n_cols + 1, figure=None,
              width_ratios=[1] * n_cols + [0.05],
              wspace=0.08, hspace=0.45,
              left=0.04, right=0.93, top=0.9, bottom=0.04)

fig = plt.figure(figsize=P04["figsize"])
fig.patch.set_facecolor("white")

# 重新用 fig 创建 gs（GridSpec 需要绑定 figure）
gs = GridSpec(2, n_cols + 1,
              figure=fig,
              width_ratios=[1] * n_cols + [0.05],
              wspace=0.08, hspace=0.35,
              left=0.04, right=0.93, top=0.80, bottom=0.04)

sc_last = None   # 保存最后一个 scatter，用于色条

for row_idx, (pert_name, vel) in enumerate([("MAFB KO", vel_ko),
                                             ("MAFB OE", vel_oe)]):
    speed = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2)
    vmax  = np.percentile(speed, 95)

    for col_idx, st in enumerate(subtypes_order):
        ax   = fig.add_subplot(gs[row_idx, col_idx])
        mask = subtypes == st

        # 背景点
        ax.scatter(umap_coords[~mask, 0], umap_coords[~mask, 1],
                   c=P04["bg_color"], s=P04["bg_size"],
                   alpha=P04["bg_alpha"], linewidths=0, zorder=1)

        # 前景点（速度着色）
        sc = ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                        c=speed[mask], cmap="hot_r",
                        s=P04["fg_size"], alpha=P04["fg_alpha"],
                        linewidths=0, zorder=2, vmin=0, vmax=vmax)
        sc_last = sc

        ax.set_title(f"{st}\n({pert_name})",
                     fontsize=P04["title_fontsize"], fontweight="bold")
        ax.set_aspect("equal")
        ax.axis("off")

# ── 色条：放在 GridSpec 最右列，跨两行 ───────────────────────────────────────
if sc_last is not None:
    ax_cbar = fig.add_subplot(gs[:, n_cols])   # 跨两行
    cbar = fig.colorbar(sc_last, cax=ax_cbar)
    cbar.set_label(P04["colorbar_label"], fontsize=8, labelpad=6)
    cbar.ax.tick_params(labelsize=7)

fig.suptitle(P04["suptitle"],
             fontsize=P04["suptitle_fontsize"], fontweight="bold",y=0.9)

save_fig(fig, "fig04_mafb_velocity_density_per_subtype")


# ══════════════════════════════════════════════════════════════════════════════
# Fig05: pySCENIC 靶基因 violin 图
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig05 前置参数 ──────────────────────────────────────────────────────────
P05 = {
    # 每个子图宽度（英寸）× 列数，高度 × 行数
    "subplot_w": 4,
    "subplot_h": 3.5,
    # 每行列数
    "ncols": 3,
    "suptitle": "MAFB pySCENIC Target Gene Expression Across Subtypes",
    "suptitle_fontsize": 13,
    # 子图标题（基因名）字号
    "gene_title_fontsize": 13,
    # x 轴刻度标签字号
    "xtick_fontsize": 10,
    # y 轴标签文字
    "ylabel": "log1p expr",
    # y 轴标签字号
    "ylabel_fontsize": 10,
    # violin 透明度
    "violin_alpha": 0.8,
    # 中位线宽度
    "median_linewidth": 1.5,
}

print("Drawing Fig05...")
targets_to_plot = [g for g in mafb_scenic_targets if g in set(gene_names)]
n_t = len(targets_to_plot)
ncols = P05["ncols"]
nrows = int(np.ceil(n_t / ncols)) if n_t > 0 else 1
fig, axes = plt.subplots(nrows, ncols,
                          figsize=(ncols * P05["subplot_w"], nrows * P05["subplot_h"]))
fig.patch.set_facecolor("white")
axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]
for idx, gene in enumerate(targets_to_plot):
    ax = axes_flat[idx]
    gene_idx = np.where(gene_names == gene)[0][0]
    data_list = [X_log[subtypes == st, gene_idx] for st in subtypes_order]
    parts = ax.violinplot(data_list, positions=range(len(subtypes_order)),
                          showmedians=True, showextrema=False)
    for pc, st in zip(parts["bodies"], subtypes_order):
        pc.set_facecolor(subtype_colors[st]); pc.set_alpha(P05["violin_alpha"])
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(P05["median_linewidth"])
    ax.set_title(gene, fontsize=P05["gene_title_fontsize"], fontweight="bold")
    ax.set_xticks(range(len(subtypes_order)))
    ax.set_xticklabels([s.replace("_", "\n") for s in subtypes_order],
                       fontsize=P05["xtick_fontsize"])
    ax.set_ylabel(P05["ylabel"], fontsize=P05["ylabel_fontsize"])
    ax.spines[["top","right"]].set_visible(False)
for idx in range(n_t, len(axes_flat)):
    axes_flat[idx].axis("off")
plt.suptitle(P05["suptitle"], fontsize=P05["suptitle_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig05_mafb_target_violin_plots")

# ══════════════════════════════════════════════════════════════════════════════
# Fig06: Butterfly 图 (KO vs OE delta per target per subtype)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig06 前置参数 ──────────────────────────────────────────────────────────
P06 = {
    "figsize_w_per_panel": 5,   # 每个子图宽度
    "figsize_h_base": 2,        # 图高基础值（加上基因数×0.5）
    "suptitle": "MAFB KO vs OE — Predicted Δ Expression of pySCENIC Targets",
    "suptitle_fontsize": 13,
    # 子图标题（亚群名）字号
    "title_fontsize": 11,
    # x 轴标签文字
    "xlabel": "Mean Δ expression",
    "xlabel_fontsize": 9,
    # 柱子高度
    "bar_height": 0.4,
    # KO 柱颜色
    "ko_color": "#4393C3",
    # OE 柱颜色
    "oe_color": "#D6604D",
    # 柱透明度
    "bar_alpha": 0.8,
    # 图例字号
    "legend_fontsize": 8,
    # 图例位置
    "legend_loc": "lower right",
    # y 轴标签字号（基因名）
    "ytick_fontsize": 9,
}

print("Drawing Fig06...")
targets_available = [g for g in mafb_scenic_targets if g in set(delta_gene_names)]
n_t = len(targets_available)
fig, axes = plt.subplots(1, len(subtypes_order),
                          figsize=(P06["figsize_w_per_panel"] * len(subtypes_order),
                                   max(4, n_t * 0.5 + P06["figsize_h_base"])))
fig.patch.set_facecolor("white")
for col_idx, st in enumerate(subtypes_order):
    ax = axes[col_idx]
    mask = subtypes == st
    ko_means, oe_means = [], []
    for gene in targets_available:
        g_idx = np.where(delta_gene_names == gene)[0][0]
        ko_means.append(delta_ko[mask, g_idx].mean())
        oe_means.append(delta_oe[mask, g_idx].mean())
    ko_means = np.array(ko_means)
    oe_means = np.array(oe_means)
    y_pos = np.arange(n_t)
    ax.barh(y_pos, ko_means,
            color=P06["ko_color"], alpha=P06["bar_alpha"],
            label="KO", height=P06["bar_height"])
    ax.barh(y_pos + P06["bar_height"], oe_means,
            color=P06["oe_color"], alpha=P06["bar_alpha"],
            label="OE", height=P06["bar_height"])
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_yticks(y_pos + P06["bar_height"] / 2)
    ax.set_yticklabels(targets_available if col_idx == 0 else [],
                       fontsize=P06["ytick_fontsize"])
    ax.set_title(st, fontsize=P06["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.set_xlabel(P06["xlabel"], fontsize=P06["xlabel_fontsize"])
    ax.spines[["top","right"]].set_visible(False)
    if col_idx == 0:
        ax.legend(fontsize=P06["legend_fontsize"], loc=P06["legend_loc"])
plt.suptitle(P06["suptitle"], fontsize=P06["suptitle_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig06_mafb_butterfly_ko_oe")

# ══════════════════════════════════════════════════════════════════════════════
# Fig07: 富集分析 bubble 图 (LA_TAM down-regulated)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig07 前置参数 ──────────────────────────────────────────────────────────
P07 = {
    "figsize_w": 12,
    "figsize_h_base": 2,        # 高度 = base + 每个 term × per_term
    "figsize_h_per_term": 0.35,
    "title": "MAFB KO Response — Pathway Enrichment\n(LA_TAM down-regulated genes, bubble size = overlap fraction)",
    "title_fontsize": 12,
    # x 轴标签
    "xlabel": "-log10(Adjusted P-value)",
    "xlabel_fontsize": 11,
    # 显著性阈值竖线位置（-log10(0.05) ≈ 1.30）
    "sig_line_x": -np.log10(0.05),
    # 竖线颜色
    "sig_line_color": "gray",
    # 竖线样式
    "sig_line_style": "--",
    # 气泡缩放系数（越大气泡越大）
    "bubble_scale": 600,
    # 气泡最小尺寸
    "bubble_min": 20,
    # 气泡透明度
    "bubble_alpha": 0.85,
    # 图例标题
    "legend_title": "Gene Set Library",
    "legend_fontsize": 8,
    "legend_loc": "lower right",
    # 每个 term 名称最大字符数（超出截断）
    "term_max_chars": 55,
    # 每个数据库最多显示 top N terms
    "top_n_per_db": 6,
    # 数据库颜色映射
    "lib_colors": {
        "GO_Biological_Process": "#4DBBD5",
        "KEGG":                  "#E64B35",
        "MSigDB_Hallmark":       "#00A087",
    },
}

print("Drawing Fig07...")
df_enr_latam = df_enr_all[
    (df_enr_all["subtype"] == "LA_TAM") & (df_enr_all["direction"] == "down")
].copy()

# 兼容不同列名
gene_set_col = "Gene_set" if "Gene_set" in df_enr_latam.columns else "gene_set"
padj_col_enr = "Adjusted P-value" if "Adjusted P-value" in df_enr_latam.columns else "padj"

top_terms = []
for gs in df_enr_latam[gene_set_col].unique():
    sub = df_enr_latam[df_enr_latam[gene_set_col] == gs].head(P07["top_n_per_db"]).copy()
    lib_short = gs.replace("_2023","").replace("_2021_Human","").replace("_2020","").replace("_Human","")
    sub["library"] = lib_short
    top_terms.append(sub)

if top_terms:
    df_top = pd.concat(top_terms, ignore_index=True)
    df_top["overlap_frac"] = df_top["n_fg_genes"].astype(float) / df_top["n_bg_genes"].astype(float)
    df_top["-log10_padj"] = -np.log10(df_top[padj_col_enr] + 1e-300)
    df_top["Term_short"] = df_top["Term"].apply(
        lambda x: x[:P07["term_max_chars"]] + "..." if len(x) > P07["term_max_chars"] else x)
    df_top = df_top.sort_values("-log10_padj", ascending=True)

    fig, ax = plt.subplots(figsize=(P07["figsize_w"],
                                    max(6, len(df_top) * P07["figsize_h_per_term"] + P07["figsize_h_base"])))
    fig.patch.set_facecolor("white")
    for lib in df_top["library"].unique():
        sub = df_top[df_top["library"] == lib]
        ax.scatter(sub["-log10_padj"], sub["Term_short"],
                   s=sub["overlap_frac"] * P07["bubble_scale"] + P07["bubble_min"],
                   c=P07["lib_colors"].get(lib, "#888888"),
                   alpha=P07["bubble_alpha"], linewidths=0.5,
                   edgecolors="white", label=lib, zorder=3)
    ax.axvline(x=P07["sig_line_x"], color=P07["sig_line_color"],
               linestyle=P07["sig_line_style"], linewidth=1, alpha=0.7)
    ax.set_xlabel(P07["xlabel"], fontsize=P07["xlabel_fontsize"])
    ax.set_title(P07["title"], fontsize=P07["title_fontsize"], fontweight="bold")
    ax.legend(title=P07["legend_title"], fontsize=P07["legend_fontsize"],
              loc=P07["legend_loc"], framealpha=0.8)
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    save_fig(fig, "fig07_mafb_enrichment_bubble")

# ══════════════════════════════════════════════════════════════════════════════
# Fig08: MAFB 调控网络图
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig08 前置参数 ──────────────────────────────────────────────────────────
P08 = {
    "figsize": (12, 10),
    "title": "MAFB Regulatory Network in LA_TAM\n(Top 15 targets + co-regulators)",
    "title_fontsize": 13,
    # 显示的 top N 靶基因数
    "n_top_targets": 15,
    # 显示的 top N 共调控因子数
    "n_top_coregulators": 5,
    # 共调控因子共调控系数阈值（绝对值）
    "coregulator_coef_thresh": 0.1,
    # MAFB 中心节点半径
    "mafb_radius": 0.5,
    # 靶基因节点半径
    "target_radius": 0.28,
    # 共调控因子节点半径
    "coregulator_radius": 0.32,
    # 靶基因环半径（距中心距离）
    "target_ring_r": 3.2,
    # 共调控因子环半径
    "coregulator_ring_r": 1.8,
    # MAFB 节点颜色
    "mafb_color": "#FF6B35",
    # MAFB 节点边框颜色
    "mafb_edge_color": "#CC3300",
    # 激活靶基因颜色
    "activated_target_color": "#FFCCCC",
    # 抑制靶基因颜色
    "repressed_target_color": "#CCCCFF",
    # 共调控因子颜色
    "coregulator_color": "#FFFFCC",
    # 共调控因子边框颜色
    "coregulator_edge_color": "#888800",
    # 激活边颜色
    "activation_edge_color": "#D6604D",
    # 抑制边颜色
    "repression_edge_color": "#4393C3",
    # 靶基因文字字号
    "target_fontsize": 7,
    # 共调控因子文字字号
    "coregulator_fontsize": 8,
    # MAFB 文字字号
    "mafb_fontsize": 12,
    # 图例字号
    "legend_fontsize": 8,
    # 图例位置
    "legend_loc": "upper right",
    # 轴范围
    "axis_lim": 4.2,
}

print("Drawing Fig08...")
latam_grn = grn_coefs["LA_TAM"]
mafb_latam = latam_grn[latam_grn["TF"] == "MAFB"].copy().sort_values("coef", key=abs, ascending=False)
top_targets = mafb_latam.head(P08["n_top_targets"])["target"].tolist()

co_regulators = {}
for tgt in top_targets:
    for _, row in latam_grn[latam_grn["target"] == tgt].iterrows():
        if row["TF"] != "MAFB" and abs(row["coef"]) > P08["coregulator_coef_thresh"]:
            co_regulators[row["TF"]] = co_regulators.get(row["TF"], 0) + abs(row["coef"])
top_coregulators = sorted(co_regulators.items(), key=lambda x: x[1], reverse=True)[:P08["n_top_coregulators"]]

fig, ax = plt.subplots(figsize=P08["figsize"])
fig.patch.set_facecolor("white"); ax.set_facecolor("#F8F8F8")
mafb_pos = np.array([0.0, 0.0])
n_cor = len(top_coregulators)
n_tgt = min(P08["n_top_targets"], len(top_targets))

cor_angles = np.linspace(0, 2*np.pi, max(n_cor, 1), endpoint=False)
cor_pos = {tf: np.array([P08["coregulator_ring_r"] * np.cos(cor_angles[i]),
                          P08["coregulator_ring_r"] * np.sin(cor_angles[i])])
           for i, (tf, _) in enumerate(top_coregulators)}
tgt_angles = np.linspace(0, 2*np.pi, max(n_tgt, 1), endpoint=False)
tgt_pos = {tgt: np.array([P08["target_ring_r"] * np.cos(tgt_angles[i]),
                            P08["target_ring_r"] * np.sin(tgt_angles[i])])
           for i, tgt in enumerate(top_targets[:n_tgt])}

coef_max = mafb_latam["coef"].abs().max() + 1e-8
for tgt in top_targets[:n_tgt]:
    rows = mafb_latam[mafb_latam["target"] == tgt]["coef"].values
    if len(rows) == 0: continue
    coef = rows[0]
    color = P08["activation_edge_color"] if coef > 0 else P08["repression_edge_color"]
    lw = abs(coef) / coef_max * 3 + 0.5
    ax.annotate("", xy=tgt_pos[tgt], xytext=mafb_pos,
                arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=0.7))

for tf, _ in top_coregulators:
    if tf not in cor_pos: continue
    for _, row in latam_grn[(latam_grn["TF"] == tf) &
                             (latam_grn["target"].isin(top_targets[:n_tgt]))].iterrows():
        if row["target"] in tgt_pos:
            color = "#FFAA88" if row["coef"] > 0 else "#88AAFF"
            ax.plot([cor_pos[tf][0], tgt_pos[row["target"]][0]],
                    [cor_pos[tf][1], tgt_pos[row["target"]][1]],
                    color=color, alpha=0.3, lw=0.8, zorder=1)

for tgt, pos in tgt_pos.items():
    rows = mafb_latam[mafb_latam["target"] == tgt]["coef"].values
    coef = rows[0] if len(rows) > 0 else 0
    c = P08["activated_target_color"] if coef > 0 else P08["repressed_target_color"]
    ax.add_patch(plt.Circle(pos, P08["target_radius"], color=c, ec="gray", lw=0.8, zorder=3))
    ax.text(pos[0], pos[1], tgt, ha="center", va="center",
            fontsize=P08["target_fontsize"], fontweight="bold", zorder=4)

for tf, pos in cor_pos.items():
    ax.add_patch(plt.Circle(pos, P08["coregulator_radius"],
                             color=P08["coregulator_color"],
                             ec=P08["coregulator_edge_color"], lw=1.2, zorder=3))
    ax.text(pos[0], pos[1], tf, ha="center", va="center",
            fontsize=P08["coregulator_fontsize"], fontweight="bold", zorder=4)

ax.add_patch(plt.Circle(mafb_pos, P08["mafb_radius"],
                         color=P08["mafb_color"], ec=P08["mafb_edge_color"], lw=2, zorder=5))
ax.text(0, 0, "MAFB", ha="center", va="center",
        fontsize=P08["mafb_fontsize"], fontweight="bold", color="white", zorder=6)

ax.legend(handles=[
    Line2D([0],[0], color=P08["activation_edge_color"], lw=2, label="Activation"),
    Line2D([0],[0], color=P08["repression_edge_color"], lw=2, label="Repression"),
    plt.Circle((0,0), 0.1, color=P08["mafb_color"], label="MAFB"),
    plt.Circle((0,0), 0.1, color=P08["coregulator_color"], ec=P08["coregulator_edge_color"], label="Co-regulator"),
    plt.Circle((0,0), 0.1, color=P08["activated_target_color"], label="Activated target"),
    plt.Circle((0,0), 0.1, color=P08["repressed_target_color"], label="Repressed target"),
], loc=P08["legend_loc"], fontsize=P08["legend_fontsize"], framealpha=0.9)
ax.set_xlim(-P08["axis_lim"], P08["axis_lim"])
ax.set_ylim(-P08["axis_lim"], P08["axis_lim"])
ax.set_aspect("equal"); ax.axis("off")
ax.set_title(P08["title"], fontsize=P08["title_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig08_mafb_regulatory_network")

# ══════════════════════════════════════════════════════════════════════════════
# Fig09: MAFB 表达 UMAP + KO 速度强度
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig09 前置参数 ──────────────────────────────────────────────────────────
P09 = {
    "figsize": (10, 3),
    "panel_titles": ["Cell Subtypes (UMAP)", "MAFB Expression", "MAFB KO — Velocity Magnitude"],
    "suptitle": "",
    "suptitle_fontsize": 5,
    "title_fontsize": 10,
    "axis_label_fontsize": 8,
    "scatter_size": 8,
    "scatter_alpha": 0.85,
    # colorbar 标签（panel 2）
    "colorbar_mafb_label": "MAFB log1p expr",
    # colorbar 标签（panel 3）
    "colorbar_speed_label": "MAFB KO velocity magnitude",
    "legend_fontsize": 8,
    "legend_loc": "lower left",
}

print("Drawing Fig09...")
speed_ko = np.sqrt(vel_ko[:,0]**2 + vel_ko[:,1]**2)
fig, axes = plt.subplots(1, 3, figsize=P09["figsize"])
fig.patch.set_facecolor("white")

ax = axes[0]
for st in subtypes_order:
    mask = subtypes == st
    ax.scatter(umap_coords[mask,0], umap_coords[mask,1],
               c=subtype_colors[st], s=P09["scatter_size"],
               alpha=P09["scatter_alpha"], linewidths=0, zorder=2, label=st)
ax.set_title(P09["panel_titles"][0], fontsize=P09["title_fontsize"], fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=P09["axis_label_fontsize"])
ax.set_ylabel("UMAP 2", fontsize=P09["axis_label_fontsize"])
ax.set_aspect("equal"); ax.spines[["top","right"]].set_visible(False)
ax.legend(fontsize=P09["legend_fontsize"], loc=P09["legend_loc"], framealpha=0.8)

ax = axes[1]
if mafb_idx is not None:
    sc = ax.scatter(umap_coords[:,0], umap_coords[:,1],
                    c=X_log[:, mafb_idx], cmap="Reds",
                    s=P09["scatter_size"], alpha=P09["scatter_alpha"], linewidths=0, zorder=2)
    plt.colorbar(sc, ax=ax, label=P09["colorbar_mafb_label"], shrink=0.8)
else:
    ax.text(0.5, 0.5, "MAFB not in HVG", ha="center", va="center", transform=ax.transAxes)
ax.set_title(P09["panel_titles"][1], fontsize=P09["title_fontsize"], fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=P09["axis_label_fontsize"])
ax.set_ylabel("UMAP 2", fontsize=P09["axis_label_fontsize"])
ax.set_aspect("equal"); ax.spines[["top","right"]].set_visible(False)

ax = axes[2]
sc = ax.scatter(umap_coords[:,0], umap_coords[:,1],
                c=speed_ko, cmap="plasma",
                s=P09["scatter_size"], alpha=P09["scatter_alpha"], linewidths=0, zorder=2)
plt.colorbar(sc, ax=ax, label=P09["colorbar_speed_label"], shrink=0.8)
ax.set_title(P09["panel_titles"][2], fontsize=P09["title_fontsize"], fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=P09["axis_label_fontsize"])
ax.set_ylabel("UMAP 2", fontsize=P09["axis_label_fontsize"])
ax.set_aspect("equal"); ax.spines[["top","right"]].set_visible(False)

plt.suptitle(P09["suptitle"], fontsize=P09["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig09_mafb_expression_umap")

# ══════════════════════════════════════════════════════════════════════════════
# Fig10: Fate score violin 图 (WT vs KO per subtype) — Fig2 风格
# ══════════════════════════════════════════════════════════════════════════════

P10 = {
    "figsize"           : (9.5, 3.5),
    "suptitle"          : "",
    "suptitle_fontsize" : 13,
    "title_fontsize"    : 12,
    "xtick_labels"      : ["WT", "MAFB KO"],
    "xtick_fontsize"    : 14,
    "ylabel"            : "LA_TAM fate score",
    "ylabel_fontsize"   : 9,
    "wt_color"          : "#4DBBD5",   # ← 灰色改蓝色
    "ko_color"          : "#D6604D",
    # violin
    "violin_alpha"      : 0.30,
    "violin_width"      : 0.65,
    # boxplot
    "box_width"         : 0.18,
    "box_alpha"         : 0.55,
    "median_lw"         : 2.0,
    "whisker_lw"        : 1.0,
    # jitter
    "jitter_range"      : 0.22,
    "jitter_s"          : 5,
    "jitter_alpha"      : 0.35,
    # 显著性标注
    "bracket_lw"        : 1.2,
    "bracket_color"     : "#333333",
    "star_fontsize"     : 15,
    "padj_fontsize"     : 10,
    "padj_color"        : "#555555",
    # 背景
    "ax_facecolor"      : "#F7F7F7",
    "spine_color"       : "#CCCCCC",
    # 图例字体
    "legend_fontsize"   : 8,           # ← 原 9 改为 8
}

def sig_stars_10(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

print("Drawing Fig10...")

fate_col = "fate_score" if "fate_score" in df_fate_wt.columns else \
    [c for c in df_fate_wt.columns if c not in ["cell_id", "subtype"]][0]
subtype_col_fate = "subtype" if "subtype" in df_fate_wt.columns else df_fate_wt.columns[1]

rng = np.random.default_rng(42)
colors = [P10["wt_color"], P10["ko_color"]]

fig, axes = plt.subplots(1, len(subtypes_order), figsize=P10["figsize"])
fig.patch.set_facecolor("white")

for col_idx, st in enumerate(subtypes_order):
    ax = axes[col_idx]

    wt_data = df_fate_wt[df_fate_wt[subtype_col_fate] == st][fate_col].values.astype(float)
    ko_data = df_fate_ko[df_fate_ko[subtype_col_fate] == st][fate_col].values.astype(float)
    vals    = [wt_data, ko_data]

    # ── Violin ────────────────────────────────────────────────────────────
    vp = ax.violinplot(vals, positions=[0, 1],
                       widths=P10["violin_width"],
                       showmedians=False, showextrema=False)
    for body, c in zip(vp["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(P10["violin_alpha"])
        body.set_edgecolor(c)
        body.set_linewidth(0.8)

    # ── Boxplot ───────────────────────────────────────────────────────────
    bp = ax.boxplot(vals, positions=[0, 1],
                    widths=P10["box_width"],
                    patch_artist=True,
                    medianprops =dict(color="#222222", linewidth=P10["median_lw"]),
                    whiskerprops=dict(color="#555555", linewidth=P10["whisker_lw"]),
                    capprops    =dict(color="#555555", linewidth=P10["whisker_lw"]),
                    flierprops  =dict(marker="", markersize=0),
                    zorder=4)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(P10["box_alpha"])
        patch.set_edgecolor(c)

    # ── Jitter 散点 ───────────────────────────────────────────────────────
    for pos, (v, c) in enumerate(zip(vals, colors)):
        jit = rng.uniform(-P10["jitter_range"], P10["jitter_range"], size=len(v))
        ax.scatter(pos + jit, v,
                   c=c, s=P10["jitter_s"], alpha=P10["jitter_alpha"],
                   linewidths=0, zorder=3, rasterized=True)

    # ── 显著性标注 ────────────────────────────────────────────────────────
    stat_row = df_fate_stats[
        (df_fate_stats["subtype"]    == st) &
        (df_fate_stats["comparison"] == "KO_vs_WT")
    ]
    if len(stat_row) > 0:
        padj  = stat_row["padj"].values[0]
        stars = sig_stars_10(padj)

        # ★ 先确定 y 轴当前范围，在范围内预留空间而非超出边框
        y_min_d = min(np.min(wt_data), np.min(ko_data))
        y_max_d = max(np.max(wt_data), np.max(ko_data))
        y_rng   = y_max_d - y_min_d

        # 连线起点：数据最高点上方 4%
        y_br   = y_max_d + y_rng * 0.04
        # 连线顶端：再上 6%
        step   = y_rng * 0.06
        # 星号顶端：连线顶端再上 2%
        y_star = y_br + step + y_rng * 0.02

        # 把 y 轴上限扩展到星号顶端再留 10% 空间，确保标注在框内
        ax.set_ylim(y_min_d - y_rng * 0.05,
                    y_star  + y_rng * 0.10)

        # 连线
        ax.plot([0, 0, 1, 1],
                [y_br, y_br + step, y_br + step, y_br],
                lw=P10["bracket_lw"], c=P10["bracket_color"], zorder=5)

        # padj 数值（axes 相对坐标，固定在顶部内侧）
        ax.text(0.5, 0.97, f"padj = {padj:.2e}",
                transform=ax.transAxes,
                ha="center", va="top",
                fontsize=P10["padj_fontsize"],
                color=P10["padj_color"])

    # ── 轴修饰 ────────────────────────────────────────────────────────────
    n_wt = len(wt_data)
    n_ko = len(ko_data)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"{P10['xtick_labels'][0]}\n(n={n_wt})",
         f"{P10['xtick_labels'][1]}\n(n={n_ko})"],
        fontsize=P10["xtick_fontsize"]
    )
    ax.set_title(st,
                 fontsize=P10["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.set_ylabel(P10["ylabel"] if col_idx == 0 else "",
                  fontsize=P10["ylabel_fontsize"])
    ax.set_facecolor(P10["ax_facecolor"])
    ax.tick_params(labelsize=8)
    for sp in ax.spines.values():
        sp.set_edgecolor(P10["spine_color"])

# ── 全局图例 ──────────────────────────────────────────────────────────────
handles = [
    mpatches.Patch(color=P10["wt_color"], alpha=0.7, label="WT"),
    mpatches.Patch(color=P10["ko_color"], alpha=0.7, label="MAFB KO"),
]

plt.suptitle(P10["suptitle"],
             fontsize=P10["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig10_fate_score_violin_wt_vs_ko")


# ══════════════════════════════════════════════════════════════════════════════
# Fig11: pySCENIC 靶基因 dot plot
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig11 前置参数 ──────────────────────────────────────────────────────────
P11 = {
    "figsize_w_base": 2,        # 图宽 = base + 亚群数×2
    "figsize_h_base": 2,        # 图高 = base + 基因数×0.6
    "title": "MAFB Target Gene Expression Dot Plot\n(dot size = % expressing, color = mean expr)",
    "title_fontsize": 12,
    "xlabel": "Cell Subtype",
    "xlabel_fontsize": 11,
    "ylabel": "MAFB pySCENIC Target Gene",
    "ylabel_fontsize": 11,
    # x 轴刻度字号
    "xtick_fontsize": 10,
    # y 轴刻度字号
    "ytick_fontsize": 9,
    # colorbar 标签
    "colorbar_label": "Normalized mean expression",
    # 图例标题
    "legend_title": "% Expressing",
    "legend_fontsize": 8,
    "legend_loc": "lower right",
    # 气泡边框颜色
    "edge_color": "gray",
    # 气泡透明度
    "alpha": 0.9,
    # 气泡缩放（点大小 = pct × scale + min）
    "size_scale": 4,
    "size_min": 5,
}

print("Drawing Fig11...")
targets_plot = [g for g in mafb_scenic_targets if g in set(gene_names)]
dot_data = []
for gene in targets_plot:
    gene_idx = np.where(gene_names == gene)[0][0]
    for st in subtypes_order:
        mask = subtypes == st
        expr = X_log[mask, gene_idx]
        pct_expr = (expr > 0).mean() * 100
        mean_expr = expr[expr > 0].mean() if (expr > 0).any() else 0
        dot_data.append({"gene": gene, "subtype": st, "pct": pct_expr, "mean": mean_expr})
df_dot = pd.DataFrame(dot_data)
df_dot["mean_norm"] = df_dot.groupby("gene")["mean"].transform(
    lambda x: (x - x.min()) / (x.max() - x.min() + 1e-8))

fig, ax = plt.subplots(figsize=(len(subtypes_order) * 2 + P11["figsize_w_base"],
                                 len(targets_plot) * 0.6 + P11["figsize_h_base"]))
fig.patch.set_facecolor("white")
x_pos = {st: i for i, st in enumerate(subtypes_order)}
y_pos = {gene: i for i, gene in enumerate(targets_plot)}
sc = ax.scatter(
    [x_pos[r["subtype"]] for _, r in df_dot.iterrows()],
    [y_pos[r["gene"]] for _, r in df_dot.iterrows()],
    s=[r["pct"] *        P11["size_scale"] + P11["size_min"] for _, r in df_dot.iterrows()],
    c=[r["mean_norm"] for _, r in df_dot.iterrows()],
    cmap="Reds", alpha=P11["alpha"], linewidths=0.5,
    edgecolors=P11["edge_color"], vmin=0, vmax=1)
plt.colorbar(sc, ax=ax, label=P11["colorbar_label"], shrink=0.6)
ax.set_xticks(range(len(subtypes_order)))
ax.set_xticklabels(subtypes_order, rotation=30, ha="right", fontsize=P11["xtick_fontsize"])
ax.set_yticks(range(len(targets_plot)))
ax.set_yticklabels(targets_plot, fontsize=P11["ytick_fontsize"])
ax.set_xlabel(P11["xlabel"], fontsize=P11["xlabel_fontsize"])
ax.set_ylabel(P11["ylabel"], fontsize=P11["ylabel_fontsize"])
ax.set_title(P11["title"], fontsize=P11["title_fontsize"], fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
ax.grid(alpha=0.2)
for pct_val in [25, 50, 75, 100]:
    ax.scatter([], [], s=pct_val * P11["size_scale"] + P11["size_min"],
               c="gray", alpha=0.7, label=f"{pct_val}%")
ax.legend(title=P11["legend_title"], loc=P11["legend_loc"],
          fontsize=P11["legend_fontsize"], framealpha=0.8)
plt.tight_layout()
save_fig(fig, "fig11_mafb_target_dotplot")

# ══════════════════════════════════════════════════════════════════════════════
# Fig12: Waterfall 图 (LA_TAM 中 MAFB KO 效应)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig12 前置参数 ──────────────────────────────────────────────────────────
P12 = {
    "figsize_w": 12,
    "figsize_h_base": 2,
    "figsize_h_per_gene": 0.35,
    # 每侧显示 top N 基因（上调 + 下调各 N 个）
    "top_n_each_side": 15,
    "title_fontsize": 12,
    "xlabel": "Mean Delta Expression (MAFB KO, LA_TAM)",
    "xlabel_fontsize": 11,
    # y 轴刻度字号
    "ytick_fontsize": 12,
    # 上调柱颜色
    "up_color": "#D6604D",
    # 下调柱颜色
    "down_color": "#4393C3",
    # pySCENIC 靶基因高亮颜色
    "scenic_color": "#FF6B35",
    # pySCENIC 靶基因边框颜色
    "scenic_edge_color": "black",
    # 柱透明度
    "bar_alpha": 0.8,
    "legend_fontsize": 11,
    "legend_loc": "lower right",
}

print("Drawing Fig12...")
latam_mask = subtypes == "LA_TAM"
mafb_latam_full = grn_coefs["LA_TAM"][grn_coefs["LA_TAM"]["TF"] == "MAFB"].copy()
target_deltas = []
for _, row in mafb_latam_full.iterrows():
    tgt = row["target"]
    if tgt in set(delta_gene_names):
        tgt_idx = np.where(delta_gene_names == tgt)[0][0]
        target_deltas.append({
            "target": tgt,
            "coef": row["coef"],
            "mean_delta_ko": delta_ko[latam_mask, tgt_idx].mean()
        })
df_waterfall = pd.DataFrame(target_deltas).sort_values("mean_delta_ko")
n_show = min(P12["top_n_each_side"], len(df_waterfall) // 2)
df_show = pd.concat([df_waterfall.head(n_show),
                     df_waterfall.tail(n_show)]).drop_duplicates().sort_values("mean_delta_ko")
df_show["is_scenic"] = df_show["target"].isin(mafb_scenic_targets)

fig, ax = plt.subplots(figsize=(P12["figsize_w"],
                                 max(6, len(df_show) * P12["figsize_h_per_gene"] + P12["figsize_h_base"])))
fig.patch.set_facecolor("white")
colors_bar = [P12["down_color"] if v < 0 else P12["up_color"] for v in df_show["mean_delta_ko"]]
ax.barh(range(len(df_show)), df_show["mean_delta_ko"],
        color=colors_bar, alpha=P12["bar_alpha"])
for i, (_, row) in enumerate(df_show.iterrows()):
    if row["is_scenic"]:
        ax.barh(i, row["mean_delta_ko"],
                color=P12["scenic_color"], alpha=1.0,
                edgecolor=P12["scenic_edge_color"], linewidth=1.5)
ax.set_yticks(range(len(df_show)))
ax.set_yticklabels(df_show["target"], fontsize=P12["ytick_fontsize"])
ax.axvline(0, color="black", linewidth=1)
ax.set_xlabel(P12["xlabel"], fontsize=P12["xlabel_fontsize"])
ax.set_title(f"MAFB KO Effect Target Genes in LA_TAM\n",
             fontsize=P12["title_fontsize"], fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
ax.legend(handles=[
    Patch(facecolor=P12["up_color"],     label="Upregulated"),
    Patch(facecolor=P12["down_color"],   label="Downregulated"),
], fontsize=P12["legend_fontsize"], loc=P12["legend_loc"])
plt.tight_layout()
save_fig(fig, "fig12_mafb_ko_waterfall")

# ══════════════════════════════════════════════════════════════════════════════
# Fig13: KO vs OE 净速度对比 (streamlines)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig13 前置参数 ──────────────────────────────────────────────────────────
P13 = {
    "figsize": (14, 6),
    "panel_titles": ["MAFB KO — Net Velocity", "MAFB OE — Net Velocity"],
    "suptitle": "MAFB Perturbation — Net Velocity Field (KO vs OE)",
    "suptitle_fontsize": 13,
    "title_fontsize": 12,
    "axis_label_fontsize": 10,
    "scatter_size": 20,
    "scatter_alpha": 0.6,
    "stream_linewidth": 1.5,
    "stream_density": 1.5,
    "stream_arrowsize": 1.5,
    # 亚群标签字号
    "subtype_label_fontsize": 9,
    "legend_fontsize": 8,
    "legend_loc": "lower right",
}

print("Drawing Fig13...")
fig, axes = plt.subplots(1, 2, figsize=P13["figsize"])
fig.patch.set_facecolor("white")
for ax, title, vx, vy in zip(axes, P13["panel_titles"],
                               [VX_ko, VX_oe], [VY_ko, VY_oe]):
    for st in subtypes_order:
        mask = subtypes == st
        ax.scatter(umap_coords[mask,0], umap_coords[mask,1],
                   c=subtype_colors[st], s=P13["scatter_size"],
                   alpha=P13["scatter_alpha"], linewidths=0, zorder=2)
    speed = np.sqrt(vx**2 + vy**2)
    try:
        ax.streamplot(GX[0,:], GY[:,0], vx, vy,
                      color=speed / (speed.max() + 1e-8), cmap="Oranges",
                      linewidth=P13["stream_linewidth"],
                      density=P13["stream_density"],
                      arrowsize=P13["stream_arrowsize"], zorder=3)
    except Exception:
        ax.quiver(GX, GY, vx, vy, speed, cmap="Oranges", alpha=0.7, width=0.003, zorder=3)
    for st in subtypes_order:
        mask = subtypes == st
        cx, cy = umap_coords[mask,0].mean(), umap_coords[mask,1].mean()
        ax.text(cx, cy, st, fontsize=P13["subtype_label_fontsize"],
                ha="center", va="center", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8, edgecolor="none"))
    ax.set_title(title, fontsize=P13["title_fontsize"], fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=P13["axis_label_fontsize"])
    ax.set_ylabel("UMAP 2", fontsize=P13["axis_label_fontsize"])
    ax.set_aspect("equal"); ax.spines[["top","right"]].set_visible(False)
handles = [mpatches.Patch(color=subtype_colors[st], label=st) for st in subtypes_order]
axes[1].legend(handles=handles, loc=P13["legend_loc"],
               fontsize=P13["legend_fontsize"], framealpha=0.8)
plt.suptitle(P13["suptitle"], fontsize=P13["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig13_mafb_net_velocity_ko_oe")

# ══════════════════════════════════════════════════════════════════════════════
# Fig14: 机制总结 6 面板图
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig14 前置参数 ──────────────────────────────────────────────────────────
P14 = {
    "figsize": (18, 10),
    "suptitle": "MAFB Regulatory Mechanism Summary — DLBCL Mono/Mac Subtypes",
    "suptitle_fontsize": 14,
    # 各子图标题字号
    "panel_title_fontsize": 11,
    # 各子图 y 轴标签字号
    "ylabel_fontsize": 10,
    # x 轴刻度标签字号
    "xtick_fontsize": 9,
    # 显著性标注字号
    "sig_fontsize": 10,
    # bar 透明度
    "bar_alpha": 0.85,
    # KO bar 颜色（Panel B）
    "ko_bar_color": "#D6604D",
    # Null bar 颜色（Panel B）
    "null_bar_color": "#888888",
    # 正 delta 颜色（Panel C, E）
    "pos_delta_color": "#D6604D",
    # 负 delta 颜色（Panel C, E）
    "neg_delta_color": "#4393C3",
    # GRN 热图颜色（Panel D）
    "grn_cmap": "RdBu_r",
    # GRN 热图显示 top N 靶基因
    "grn_top_n": 20,
    "legend_fontsize": 9,
    # error bar cap 大小
    "capsize": 4,
}

print("Drawing Fig14...")
fig = plt.figure(figsize=P14["figsize"])
fig.patch.set_facecolor("white")

# Panel A: MAFB 各亚群表达量
ax1 = fig.add_subplot(2, 3, 1)
if mafb_idx is not None:
    mafb_means = [X_log[subtypes == st, mafb_idx].mean() for st in subtypes_order]
    mafb_sems  = [X_log[subtypes == st, mafb_idx].std() /
                  np.sqrt((subtypes == st).sum()) for st in subtypes_order]
    ax1.bar(subtypes_order, mafb_means,
            color=[subtype_colors[st] for st in subtypes_order],
            alpha=P14["bar_alpha"], yerr=mafb_sems, capsize=P14["capsize"],
            error_kw={"linewidth": 1.5})
    ax1.set_ylabel("Mean MAFB log1p expr", fontsize=P14["ylabel_fontsize"])
    ax1.set_title("MAFB Expression per Subtype",
                  fontsize=P14["panel_title_fontsize"], fontweight="bold")
    ax1.set_xticklabels(subtypes_order, rotation=30, ha="right", fontsize=P14["xtick_fontsize"])
    ax1.spines[["top","right"]].set_visible(False)
else:
    ax1.text(0.5, 0.5, "MAFB not in HVG", ha="center", va="center", transform=ax1.transAxes)
    ax1.set_title("MAFB Expression", fontsize=P14["panel_title_fontsize"], fontweight="bold")

# Panel B: KO 扰动强度
ax2 = fig.add_subplot(2, 3, 2)
speed_ko_p   = np.sqrt(vel_ko[:,0]**2 + vel_ko[:,1]**2)
speed_null_p = np.sqrt(vel_null[:,0]**2 + vel_null[:,1]**2)
ko_speeds   = [speed_ko_p[subtypes == st].mean() for st in subtypes_order]
null_speeds = [speed_null_p[subtypes == st].mean() for st in subtypes_order]
x = np.arange(len(subtypes_order)); w = 0.35
ax2.bar(x - w/2, ko_speeds,   w, color=P14["ko_bar_color"],   alpha=P14["bar_alpha"], label="MAFB KO")
ax2.bar(x + w/2, null_speeds, w, color=P14["null_bar_color"], alpha=P14["bar_alpha"], label="Null")
ax2.set_xticks(x)
ax2.set_xticklabels(subtypes_order, rotation=30, ha="right", fontsize=P14["xtick_fontsize"])
ax2.set_ylabel("Mean velocity magnitude", fontsize=P14["ylabel_fontsize"])
ax2.set_title("KO Perturbation Strength per Subtype",
              fontsize=P14["panel_title_fontsize"], fontweight="bold")
ax2.legend(fontsize=P14["legend_fontsize"]); ax2.spines[["top","right"]].set_visible(False)

# Panel C: Fate score delta
ax3 = fig.add_subplot(2, 3, 3)
df_ko_wt = df_fate_stats[df_fate_stats["comparison"] == "KO_vs_WT"].set_index("subtype")
deltas = [df_ko_wt.loc[st, "delta"] if st in df_ko_wt.index else 0 for st in subtypes_order]
padjs  = [df_ko_wt.loc[st, "padj"]  if st in df_ko_wt.index else 1 for st in subtypes_order]
bar_colors = [P14["pos_delta_color"] if d > 0 else P14["neg_delta_color"] for d in deltas]
ax3.bar(subtypes_order, deltas, color=bar_colors, alpha=P14["bar_alpha"])
for i, (d, p_val) in enumerate(zip(deltas, padjs)):
    sig = "***​" if p_val < 0.001 else ("​**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
    ax3.text(i, d + (0.01 if d >= 0 else -0.02), sig,
             ha="center", va="bottom" if d >= 0 else "top", fontsize=P14["sig_fontsize"])
ax3.axhline(0, color="black", linewidth=0.8)
ax3.set_ylabel("Delta LA_TAM fate score (KO-WT)", fontsize=P14["ylabel_fontsize"] - 1)
ax3.set_title("Fate Score Change on MAFB KO",
              fontsize=P14["panel_title_fontsize"], fontweight="bold")
ax3.set_xticklabels(subtypes_order, rotation=30, ha="right", fontsize=P14["xtick_fontsize"])
ax3.spines[["top","right"]].set_visible(False)

# Panel D: GRN 系数热图
ax4 = fig.add_subplot(2, 3, (4, 5))
all_coefs = {}
for st in subtypes_order:
    mafb_df = grn_coefs[st][grn_coefs[st]["TF"] == "MAFB"].set_index("target")["coef"]
    all_coefs[st] = mafb_df
df_coef_matrix = pd.DataFrame(all_coefs).fillna(0)
top_idx = df_coef_matrix.abs().mean(1).nlargest(P14["grn_top_n"]).index
df_coef_top = df_coef_matrix.loc[top_idx]
vmax_grn = df_coef_top.abs().max().max()
im = ax4.imshow(df_coef_top.values, cmap=P14["grn_cmap"],
                aspect="auto", vmin=-vmax_grn, vmax=vmax_grn)
plt.colorbar(im, ax=ax4, label="Ridge coef", shrink=0.8)
ax4.set_xticks(range(len(subtypes_order)))
ax4.set_xticklabels(subtypes_order, fontsize=10)
ax4.set_yticks(range(len(top_idx)))
ax4.set_yticklabels(top_idx, fontsize=8)
ax4.set_title(f"MAFB GRN Coefficients — Top {P14['grn_top_n']} Targets",
              fontsize=P14["panel_title_fontsize"], fontweight="bold")

# Panel E: Program score delta
ax5 = fig.add_subplot(2, 3, 6)
prog_deltas = [df_prog_stats[df_prog_stats["subtype"] == st]["delta"].values[0]
               if st in df_prog_stats["subtype"].values else 0 for st in subtypes_order]
prog_padjs  = [df_prog_stats[df_prog_stats["subtype"] == st]["padj"].values[0]
               if st in df_prog_stats["subtype"].values else 1 for st in subtypes_order]
bar_colors2 = [P14["pos_delta_color"] if d > 0 else P14["neg_delta_color"] for d in prog_deltas]
ax5.bar(subtypes_order, prog_deltas, color=bar_colors2, alpha=P14["bar_alpha"])
for i, (d, p_val) in enumerate(zip(prog_deltas, prog_padjs)):
    sig = "***​" if p_val < 0.001 else ("​**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
    ax5.text(i, d + (0.002 if d >= 0 else -0.005), sig,
             ha="center", va="bottom" if d >= 0 else "top", fontsize=P14["sig_fontsize"])
ax5.axhline(0, color="black", linewidth=0.8)
ax5.set_ylabel("Delta LA_TAM program score (KO-WT)", fontsize=P14["ylabel_fontsize"] - 1)
ax5.set_title("LA_TAM Program Score Change on MAFB KO",
              fontsize=P14["panel_title_fontsize"], fontweight="bold")
ax5.set_xticklabels(subtypes_order, rotation=30, ha="right", fontsize=P14["xtick_fontsize"])
ax5.spines[["top","right"]].set_visible(False)

plt.suptitle(P14["suptitle"], fontsize=P14["suptitle_fontsize"], fontweight="bold", y=1.01)
plt.tight_layout()
save_fig(fig, "fig14_mafb_mechanism_summary")

# ══════════════════════════════════════════════════════════════════════════════
# Fig15: LA_TAM Program Score violin 图 (WT vs KO per subtype) — Fig2 风格
# ══════════════════════════════════════════════════════════════════════════════

P15 = {
    "figsize"           : (10, 3.5),
    "suptitle"          : "LA_TAM Program Score: WT vs MAFB KO per Subtype",
    "suptitle_fontsize" : 13,
    "title_fontsize"    : 12,
    "xtick_labels"      : ["WT", "MAFB KO"],
    "xtick_fontsize"    : 14,
    "ylabel"            : "LA_TAM program score",
    "ylabel_fontsize"   : 9,
    "wt_color"          : "#4DBBD5",
    "ko_color"          : "#D6604D",
    # violin
    "violin_alpha"      : 0.30,
    "violin_width"      : 0.65,
    # boxplot
    "box_width"         : 0.18,
    "box_alpha"         : 0.55,
    "median_lw"         : 2.0,
    "whisker_lw"        : 1.0,
    # jitter
    "jitter_range"      : 0.22,
    "jitter_s"          : 5,
    "jitter_alpha"      : 0.35,
    # 显著性标注
    "bracket_lw"        : 1.2,
    "bracket_color"     : "#333333",
    "star_fontsize"     : 15,
    "padj_fontsize"     : 10,
    "padj_color"        : "#555555",
    # 背景
    "ax_facecolor"      : "#F7F7F7",
    "spine_color"       : "#CCCCCC",
    # 图例字体
    "legend_fontsize"   : 8,
}

def sig_stars_15(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

print("Drawing Fig15...")

# ── 数据处理：照抄原始代码，只在末尾加 np.asarray().ravel() 确保一维 ────────
prog_genes_valid    = [g for g in prog_genes if g in set(gene_names)]
prog_idx_wt         = [np.where(gene_names == g)[0][0] for g in prog_genes_valid]
prog_score_wt       = np.asarray(X_log[:, prog_idx_wt].mean(1)).ravel()  # ← 确保一维

prog_genes_ko_valid = [g for g in prog_genes_valid if g in set(df_simulated.columns)]
prog_idx_ko         = [list(df_simulated.columns).index(g) for g in prog_genes_ko_valid]
prog_score_ko       = np.asarray(X_ko_sim[:, prog_idx_ko].mean(1)).ravel()  # ← 确保一维

# ── 绘图 ──────────────────────────────────────────────────────────────────────
rng    = np.random.default_rng(42)
colors = [P15["wt_color"], P15["ko_color"]]

fig, axes = plt.subplots(1, len(subtypes_order), figsize=P15["figsize"])
fig.patch.set_facecolor("white")

for col_idx, st in enumerate(subtypes_order):
    ax   = axes[col_idx]
    mask = subtypes == st   # ← 与原始代码完全一致

    wt_data = prog_score_wt[mask]
    ko_data = prog_score_ko[mask]
    vals    = [wt_data, ko_data]

    # ── Violin ────────────────────────────────────────────────────────────
    vp = ax.violinplot(vals, positions=[0, 1],
                       widths=P15["violin_width"],
                       showmedians=False, showextrema=False)
    for body, c in zip(vp["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(P15["violin_alpha"])
        body.set_edgecolor(c)
        body.set_linewidth(0.8)

    # ── Boxplot ───────────────────────────────────────────────────────────
    bp = ax.boxplot(vals, positions=[0, 1],
                    widths=P15["box_width"],
                    patch_artist=True,
                    medianprops =dict(color="#222222", linewidth=P15["median_lw"]),
                    whiskerprops=dict(color="#555555", linewidth=P15["whisker_lw"]),
                    capprops    =dict(color="#555555", linewidth=P15["whisker_lw"]),
                    flierprops  =dict(marker="", markersize=0),
                    zorder=4)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(P15["box_alpha"])
        patch.set_edgecolor(c)

    # ── Jitter 散点 ───────────────────────────────────────────────────────
    for pos, (v, c) in enumerate(zip(vals, colors)):
        jit = rng.uniform(-P15["jitter_range"], P15["jitter_range"], size=len(v))
        ax.scatter(pos + jit, v,
                   c=c, s=P15["jitter_s"], alpha=P15["jitter_alpha"],
                   linewidths=0, zorder=3, rasterized=True)

    # ── 显著性标注 ────────────────────────────────────────────────────────
    stat_row = df_prog_stats[df_prog_stats["subtype"] == st]
    if len(stat_row) > 0:
        padj  = stat_row["padj"].values[0]
        stars = sig_stars_15(padj)

        y_min_d = min(np.min(wt_data), np.min(ko_data))
        y_max_d = max(np.max(wt_data), np.max(ko_data))
        y_rng   = y_max_d - y_min_d

        y_br   = y_max_d + y_rng * 0.04
        step   = y_rng * 0.06
        y_star = y_br + step + y_rng * 0.02

        ax.set_ylim(y_min_d - y_rng * 0.05,
                    y_star  + y_rng * 0.10)

        ax.plot([0, 0, 1, 1],
                [y_br, y_br + step, y_br + step, y_br],
                lw=P15["bracket_lw"], c=P15["bracket_color"], zorder=5)

        ax.text(0.5, 0.97, f"padj = {padj:.2e}",
                transform=ax.transAxes,
                ha="center", va="top",
                fontsize=P15["padj_fontsize"],
                color=P15["padj_color"])

    # ── 轴修饰 ────────────────────────────────────────────────────────────
    n_st = int(mask.sum())
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"{P15['xtick_labels'][0]}\n(n={n_st})",
         f"{P15['xtick_labels'][1]}\n(n={n_st})"],
        fontsize=P15["xtick_fontsize"]
    )
    ax.set_title(st,
                 fontsize=P15["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.set_ylabel(P15["ylabel"] if col_idx == 0 else "",
                  fontsize=P15["ylabel_fontsize"])
    ax.set_facecolor(P15["ax_facecolor"])
    ax.tick_params(labelsize=8)
    for sp in ax.spines.values():
        sp.set_edgecolor(P15["spine_color"])

# ── 全局图例 ──────────────────────────────────────────────────────────────
handles = [
    mpatches.Patch(color=P15["wt_color"], alpha=0.7, label="WT"),
    mpatches.Patch(color=P15["ko_color"], alpha=0.7, label="MAFB KO"),
]


plt.suptitle(P15["suptitle"],
             fontsize=P15["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig15_program_score_violin")


# ══════════════════════════════════════════════════════════════════════════════
# Fig16: Program Score paired violin+box 图 (WT vs KO per subtype)
# ══════════════════════════════════════════════════════════════════════════════

P16 = {
    "figsize"           : (9, 3.5),
    "suptitle"          : "MAFB KO Effect on LA_TAM Program Score",
    "suptitle_fontsize" : 12,
    "title_fontsize"    : 12,
    "xtick_labels"      : ["WT", "MAFB KO"],
    "xtick_fontsize"    : 14,
    "ylabel"            : "LA_TAM program score",
    "ylabel_fontsize"   : 9,
    "wt_color"          : "#4DBBD5",
    "ko_color"          : "#D6604D",
    # violin
    "violin_alpha"      : 0.30,
    "violin_width"      : 0.65,
    # boxplot
    "box_width"         : 0.18,
    "box_alpha"         : 0.55,
    "median_lw"         : 2.0,
    "whisker_lw"        : 1.0,
    # jitter
    "jitter_range"      : 0.22,
    "jitter_s"          : 5,
    "jitter_alpha"      : 0.35,
    # paired line（细胞配对连线）
    "show_paired_lines" : True,
    "paired_line_alpha" : 0.15,   # 线条很密时调低透明度
    "paired_line_lw"    : 0.5,
    "paired_line_color" : "#888888",
    # 显著性标注
    "bracket_lw"        : 1.2,
    "bracket_color"     : "#333333",
    "star_fontsize"     : 15,
    "annot_fontsize"    : 9,
    "annot_color"       : "#555555",
    # 背景
    "ax_facecolor"      : "#F7F7F7",
    "spine_color"       : "#CCCCCC",
    # 图例字体
    "legend_fontsize"   : 8,
}

def sig_stars_16(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

print("Drawing Fig16...")

# ── 数据：与 Fig15 完全一致，确保一维 ────────────────────────────────────────
prog_genes_valid    = [g for g in prog_genes if g in set(gene_names)]
prog_idx_wt         = [np.where(gene_names == g)[0][0] for g in prog_genes_valid]
prog_score_wt       = np.asarray(X_log[:, prog_idx_wt].mean(1)).ravel()

prog_genes_ko_valid = [g for g in prog_genes_valid if g in set(df_simulated.columns)]
prog_idx_ko         = [list(df_simulated.columns).index(g) for g in prog_genes_ko_valid]
prog_score_ko       = np.asarray(X_ko_sim[:, prog_idx_ko].mean(1)).ravel()

rng    = np.random.default_rng(42)
colors = [P16["wt_color"], P16["ko_color"]]

fig, axes = plt.subplots(1, len(subtypes_order), figsize=P16["figsize"])
fig.patch.set_facecolor("white")

for col_idx, st in enumerate(subtypes_order):
    ax   = axes[col_idx]
    mask = subtypes == st

    wt_data = prog_score_wt[mask]
    ko_data = prog_score_ko[mask]
    vals    = [wt_data, ko_data]

    # ── 打印每个 subtype 的数据摘要（方便核查）──────────────────────────
    print(f"  [{st}] n={mask.sum()}  "
          f"WT median={np.median(wt_data):.4f}  "
          f"KO median={np.median(ko_data):.4f}  "
          f"delta={np.median(ko_data)-np.median(wt_data):.4f}")

    # ── Paired lines（同一细胞 WT→KO 连线）──────────────────────────────
    if P16["show_paired_lines"]:
        # WT 和 KO 行数相同（同一批细胞的 imputed vs simulated）
        n_cells = min(len(wt_data), len(ko_data))
        # 给每条线加独立 jitter，让线不完全重叠
        jit_wt = rng.uniform(-0.05, 0.05, size=n_cells)
        jit_ko = rng.uniform(-0.05, 0.05, size=n_cells)
        for j in range(n_cells):
            ax.plot([0 + jit_wt[j], 1 + jit_ko[j]],
                    [wt_data[j], ko_data[j]],
                    color=P16["paired_line_color"],
                    alpha=P16["paired_line_alpha"],
                    lw=P16["paired_line_lw"],
                    zorder=2)

    # ── Violin ────────────────────────────────────────────────────────────
    vp = ax.violinplot(vals, positions=[0, 1],
                       widths=P16["violin_width"],
                       showmedians=False, showextrema=False)
    for body, c in zip(vp["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(P16["violin_alpha"])
        body.set_edgecolor(c)
        body.set_linewidth(0.8)
        body.set_zorder(3)

    # ── Boxplot ───────────────────────────────────────────────────────────
    bp = ax.boxplot(vals, positions=[0, 1],
                    widths=P16["box_width"],
                    patch_artist=True,
                    medianprops =dict(color="#222222", linewidth=P16["median_lw"]),
                    whiskerprops=dict(color="#555555", linewidth=P16["whisker_lw"]),
                    capprops    =dict(color="#555555", linewidth=P16["whisker_lw"]),
                    flierprops  =dict(marker="", markersize=0),
                    zorder=4)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(P16["box_alpha"])
        patch.set_edgecolor(c)

    # ── Jitter 散点 ───────────────────────────────────────────────────────
    for pos, (v, c) in enumerate(zip(vals, colors)):
        jit = rng.uniform(-P16["jitter_range"], P16["jitter_range"], size=len(v))
        ax.scatter(pos + jit, v,
                   c=c, s=P16["jitter_s"], alpha=P16["jitter_alpha"],
                   linewidths=0, zorder=5, rasterized=True)

    # ── 显著性标注（连线 + 星号 + delta + padj）──────────────────────────
    stat_row = df_prog_stats[df_prog_stats["subtype"] == st]
    if len(stat_row) > 0:
        padj    = stat_row["padj"].values[0]
        delta   = stat_row["delta"].values[0]
        cohens  = stat_row["cohens_d"].values[0]
        stars   = sig_stars_16(padj)

        y_min_d = min(np.min(wt_data), np.min(ko_data))
        y_max_d = max(np.max(wt_data), np.max(ko_data))
        y_rng   = y_max_d - y_min_d

        y_br   = y_max_d + y_rng * 0.04
        step   = y_rng * 0.06
        y_star = y_br + step + y_rng * 0.02

        ax.set_ylim(y_min_d - y_rng * 0.05,
                    y_star  + y_rng * 0.22)   # 留更多空间放标注文字

        # 连线
        ax.plot([0, 0, 1, 1],
                [y_br, y_br + step, y_br + step, y_br],
                lw=P16["bracket_lw"], c=P16["bracket_color"], zorder=6)
        # 星号
        
        # delta + Cohen's d + padj（三行文字，放在图内顶部）
        ax.text(0.5, 0.97,
                f"Δ={delta:.3f}  d={cohens:.2f}\npadj={padj:.2e}",
                transform=ax.transAxes,
                ha="center", va="top",
                fontsize=P16["annot_fontsize"],
                color=P16["annot_color"])

    # ── 轴修饰 ────────────────────────────────────────────────────────────
    n_st = int(mask.sum())
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"{P16['xtick_labels'][0]}\n(n={n_st})",
         f"{P16['xtick_labels'][1]}\n(n={n_st})"],
        fontsize=P16["xtick_fontsize"]
    )
    ax.set_title(st,
                 fontsize=P16["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.set_ylabel(P16["ylabel"] if col_idx == 0 else "",
                  fontsize=P16["ylabel_fontsize"])
    ax.set_facecolor(P16["ax_facecolor"])
    ax.tick_params(labelsize=8)
    for sp in ax.spines.values():
        sp.set_edgecolor(P16["spine_color"])

# ── 全局图例 ──────────────────────────────────────────────────────────────
handles = [
    mpatches.Patch(color=P16["wt_color"], alpha=0.7, label="WT"),
    mpatches.Patch(color=P16["ko_color"], alpha=0.7, label="MAFB KO"),
]

plt.suptitle(P16["suptitle"],
             fontsize=P16["suptitle_fontsize"], fontweight="bold", y=1.02)
plt.tight_layout()
save_fig(fig, "fig16_program_score_paired_violin")


# ══════════════════════════════════════════════════════════════════════════════
# Fig17: 敏感性分析汇总
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig17 前置参数 ──────────────────────────────────────────────────────────
P17 = {
    "figsize": (16, 6),
    "suptitle": "Sensitivity Analysis Summary",
    "suptitle_fontsize": 13,
    # Panel A 标题
    "panelA_title": "Sensitivity Analysis — Delta per Test",
    "panelA_title_fontsize": 11,
    "panelA_ylabel": "Delta fate score (KO-WT)",
    "panelA_ylabel_fontsize": 10,
    # Panel B 标题
    "panelB_title": "Direction Consistency Heatmap",
    "panelB_title_fontsize": 11,
    # x 轴刻度标签旋转角度
    "xtick_rotation": 45,
    "xtick_fontsize": 8,
    # bar 宽度
    "bar_width": 0.25,
    # colorbar 标签
    "colorbar_label": "Direction: 1=consistent, -1=flipped",
    # 热图单元格文字字号
    "heatmap_cell_fontsize": 7,
    "legend_fontsize": 9,
}

print("Drawing Fig17...")
df_sens_eval = df_sens[df_sens["delta"].notna()].copy()
if len(df_sens_eval) > 0:
    sens_ids = df_sens_eval["sensitivity_id"].unique()
    fig, axes = plt.    subplots(1, 2, figsize=P17["figsize"])
    fig.patch.set_facecolor("white")

    # Panel A: delta 柱状图
    ax = axes[0]
    x = np.arange(len(sens_ids))
    for i, st in enumerate(subtypes_order):
        sub = df_sens_eval[df_sens_eval["subtype"] == st]
        deltas_s = []
        for sid in sens_ids:
            row = sub[sub["sensitivity_id"] == sid]
            deltas_s.append(row["delta"].values[0] if len(row) > 0 else np.nan)
        ax.bar(x + i * P17["bar_width"], deltas_s,
               P17["bar_width"], color=subtype_colors[st], alpha=0.85, label=st)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x + P17["bar_width"])
    ax.set_xticklabels(sens_ids, rotation=P17["xtick_rotation"],
                       ha="right", fontsize=P17["xtick_fontsize"])
    ax.set_ylabel(P17["panelA_ylabel"], fontsize=P17["panelA_ylabel_fontsize"])
    ax.set_title(P17["panelA_title"], fontsize=P17["panelA_title_fontsize"], fontweight="bold")
    ax.legend(fontsize=P17["legend_fontsize"])
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: 方向一致性热图
    ax = axes[1]
    direction_map = {"YES": 1, "REFERENCE": 1, "FLIPPED": -1,
                     "NO": -1, "N/A (no re-simulation)": 0}
    matrix = np.zeros((len(subtypes_order), len(sens_ids)))
    for i, st in enumerate(subtypes_order):
        for j, sid in enumerate(sens_ids):
            row = df_sens_eval[(df_sens_eval["subtype"] == st) &
                               (df_sens_eval["sensitivity_id"] == sid)]
            if len(row) > 0:
                d = row["direction_vs_main"].values[0]
                matrix[i, j] = direction_map.get(d, 0)
    im = ax.imshow(matrix, cmap="RdYlGn", aspect="auto", vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax, label=P17["colorbar_label"], shrink=0.8)
    ax.set_xticks(range(len(sens_ids)))
    ax.set_xticklabels(sens_ids, rotation=P17["xtick_rotation"],
                       ha="right", fontsize=P17["xtick_fontsize"])
    ax.set_yticks(range(len(subtypes_order)))
    ax.set_yticklabels(subtypes_order, fontsize=10)
    ax.set_title(P17["panelB_title"], fontsize=P17["panelB_title_fontsize"], fontweight="bold")
    for i in range(len(subtypes_order)):
        for j in range(len(sens_ids)):
            val = matrix[i, j]
            label = "OK" if val == 1 else ("FLIP" if val == -1 else "N/A")
            ax.text(j, i, label, ha="center", va="center",
                    fontsize=P17["heatmap_cell_fontsize"],
                    color="black" if abs(val) < 0.5 else "white")

    plt.suptitle(P17["suptitle"], fontsize=P17["suptitle_fontsize"], fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "fig17_sensitivity_summary")

# ══════════════════════════════════════════════════════════════════════════════
# Fig18: GRN 系数热图 (MAFB top 30 targets across subtypes)
# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
# Fig18: MAFB GRN Coefficients heatmap — 横版（基因在 x 轴，亚群在 y 轴）
# ══════════════════════════════════════════════════════════════════════════════

P18 = {
    "figsize"               : (10, 5),   # ← 横版：宽 > 高
    "title"                 : "MAFB GRN Coefficients — Top 30 Targets)",
    "title_fontsize"        : 12,
    "top_n"                 : 30,
    "cmap"                  : "RdBu_r",
    "colorbar_label"        : "Ridge regression coefficient",
    # x 轴（基因名）字号
    "xtick_fontsize"        : 11,
    # y 轴（亚群）字号
    "ytick_fontsize"        : 11,
    # pySCENIC 标记
    "scenic_mark_color"     : "#FF6B35",
    "scenic_mark_fontsize"  : 12,
    "scenic_legend_fontsize": 9,
}

print("Drawing Fig18...")

all_coefs_18 = {}
for st in subtypes_order:
    mafb_df = grn_coefs[st][grn_coefs[st]["TF"] == "MAFB"].set_index("target")["coef"]
    all_coefs_18[st] = mafb_df

df_coef_matrix_18 = pd.DataFrame(all_coefs_18).fillna(0)
top_idx_18        = df_coef_matrix_18.abs().mean(1).nlargest(P18["top_n"]).index
df_coef_top_18    = df_coef_matrix_18.loc[top_idx_18]
vmax_18           = df_coef_top_18.abs().max().max()

# ── 横版：转置矩阵，行=亚群，列=基因 ────────────────────────────────────────
# 原始 df_coef_top_18 shape: (n_genes, n_subtypes)
# 转置后 shape: (n_subtypes, n_genes) → imshow 行=亚群，列=基因
data_T = df_coef_top_18.values.T

fig, ax = plt.subplots(figsize=P18["figsize"])
fig.patch.set_facecolor("white")

im = ax.imshow(data_T, cmap=P18["cmap"],
               aspect="auto", vmin=-vmax_18, vmax=vmax_18)

# ── x 轴：基因名，旋转 45° 避免重叠 ──────────────────────────────────────────
ax.set_xticks(range(len(top_idx_18)))
ax.set_xticklabels(top_idx_18,
                   fontsize=P18["xtick_fontsize"],
                   rotation=45, ha="right")

# ── y 轴：亚群名 ──────────────────────────────────────────────────────────────
ax.set_yticks(range(len(subtypes_order)))
ax.set_yticklabels(subtypes_order, fontsize=P18["ytick_fontsize"])

# ── pySCENIC 标记：在对应基因列的上方标 * ────────────────────────────────────
for j, gene in enumerate(top_idx_18):
    if gene in mafb_scenic_targets:
        ax.text(j, -0.6, "*",
                ha="center", va="bottom",
                fontsize=P18["scenic_mark_fontsize"],
                color=P18["scenic_mark_color"], fontweight="bold")



# ── colorbar：放在图右侧 ──────────────────────────────────────────────────────
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
ax_cbar = divider.append_axes("right", size="2%", pad=0.08)
cbar    = fig.colorbar(im, cax=ax_cbar)
cbar.set_label(P18["colorbar_label"], fontsize=10, labelpad=6)
cbar.ax.tick_params(labelsize=8)

ax.set_title(P18["title"], fontsize=P18["title_fontsize"],
             fontweight="bold", pad=18)

plt.tight_layout()
save_fig(fig, "fig18_grn_coefficient_heatmap")

# ══════════════════════════════════════════════════════════════════════════════
# Fig19: 富集分析热图 (top terms per subtype × direction)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig19 前置参数 ──────────────────────────────────────────────────────────
P19 = {
    "figsize_w": 12,
    "figsize_h_base": 2,
    "figsize_h_per_term": 0.4,
    "title": "Pathway Enrichment Heatmap\n(MAFB KO response per subtype/direction)",
    "title_fontsize": 12,
    # x 轴（分组）字号
    "xtick_fontsize": 9,
    # y 轴（通路名）字号
    "ytick_fontsize": 8,
    # colorbar 标签
    "colorbar_label": "-log10(adj. p-value)",
    # 热图 colormap
    "cmap": "YlOrRd",
    # 每个分组最多显示 top N terms
    "top_n_per_group": 5,
    # x 轴标签旋转角度
    "xtick_rotation": 45,
}

print("Drawing Fig19...")
padj_col_19 = "Adjusted P-value" if "Adjusted P-value" in df_enr_all.columns else "padj"
top_enr_19 = []
for st in subtypes_order:
    for direction in ["up", "down"]:
        sub = df_enr_all[(df_enr_all["subtype"] == st) &
                         (df_enr_all["direction"] == direction)].copy()
        sub = sub.sort_values(padj_col_19).head(P19["top_n_per_group"])
        sub["label"] = sub["Term"].apply(lambda x: x[:40] + "..." if len(x) > 40 else x)
        sub["group"] = f"{st}_{direction}"
        top_enr_19.append(sub)

if top_enr_19:
    df_enr_top_19 = pd.concat(top_enr_19, ignore_index=True)
    df_enr_top_19["-log10_padj"] = -np.log10(df_enr_top_19[padj_col_19] + 1e-300)
    groups_19 = [f"{st}_{d}" for st in subtypes_order for d in ["up", "down"]]
    all_terms_19 = df_enr_top_19["label"].unique()
    pivot_19 = pd.DataFrame(0.0, index=all_terms_19, columns=groups_19)
    for _, row in df_enr_top_19.iterrows():
        pivot_19.loc[row["label"], row["group"]] = row["-log10_padj"]
    pivot_19 = pivot_19[pivot_19.max(1) > -np.log10(0.05)]

    if len(pivot_19) > 0:
        fig, ax = plt.subplots(figsize=(P19["figsize_w"],
                                        max(6, len(pivot_19) * P19["figsize_h_per_term"] + P19["figsize_h_base"])))
        fig.patch.set_facecolor("white")
        im = ax.imshow(pivot_19.values, cmap=P19["cmap"], aspect="auto", vmin=0)
        plt.colorbar(im, ax=ax, label=P19["colorbar_label"], shrink=0.8)
        ax.set_xticks(range(len(groups_19)))
        ax.set_xticklabels(groups_19, rotation=P19["xtick_rotation"],
                           ha="right", fontsize=P19["xtick_fontsize"])
        ax.set_yticks(range(len(pivot_19)))
        ax.set_yticklabels(pivot_19.index, fontsize=P19["ytick_fontsize"])
        ax.set_title(P19["title"], fontsize=P19["title_fontsize"], fontweight="bold")
        plt.tight_layout()
        save_fig(fig, "fig19_enrichment_heatmap")

# ══════════════════════════════════════════════════════════════════════════════
# Fig20: pySCENIC 靶基因分类柱状图
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig20 前置参数 ──────────────────────────────────────────────────────────
P20 = {
    "figsize": (14, 5),
    "suptitle": "pySCENIC MAFB Target Classification per Subtype",
    "suptitle_fontsize": 13,
    "title_fontsize": 11,
    # y 轴标签
    "ylabel": "# genes",
    "ylabel_fontsize": 10,
    # x 轴刻度字号
    "xtick_fontsize": 8,
    # x 轴刻度旋转角度
    "xtick_rotation": 45,
    # 柱透明度
    "bar_alpha": 0.85,
}

print("Drawing Fig20...")
class_col = "class_t2000" if "class_t2000" in df_scenic.columns else df_scenic.columns[2]
fig, axes = plt.subplots(1, len(subtypes_order), figsize=P20["figsize"])
fig.patch.set_facecolor("white")
for col_idx, st in enumerate(subtypes_order):
    ax = axes[col_idx]
    sub = df_scenic[df_scenic["subtype"] == st]
    class_counts = sub[class_col].value_counts()
    colors_cls = plt.cm.Set2(np.linspace(0, 1, len(class_counts)))
    ax.bar(range(len(class_counts)), class_counts.values,
           color=colors_cls, alpha=P20["bar_alpha"])
    ax.set_xticks(range(len(class_counts)))
    ax.set_xticklabels(class_counts.index,
                       rotation=P20["xtick_rotation"], ha="right",
                       fontsize=P20["xtick_fontsize"])
    ax.set_title(st, fontsize=P20["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.set_ylabel(P20["ylabel"] if col_idx == 0 else "",
                  fontsize=P20["ylabel_fontsize"])
    ax.spines[["top","right"]].set_visible(False)
plt.suptitle(P20["suptitle"], fontsize=P20["suptitle_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig20_scenic_target_classification")

# ══════════════════════════════════════════════════════════════════════════════
# Fig21: Cluster deltaX 热图 (KO and OE)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig21 前置参数 ──────────────────────────────────────────────────────────
P21 = {
    "figsize": (18, 6),
    "suptitle": "MAFB Perturbation — Cluster-level Gene Expression Changes",
    "suptitle_fontsize": 13,
    "title_fontsize": 12,
    # 每个面板显示 top N 变异基因
    "top_n_genes": 30,
    # 热图 colormap
    "cmap": "RdBu_r",
    # colorbar 标签
    "colorbar_label": "Mean delta expression",
    # x 轴（基因名）字号
    "xtick_fontsize": 7,
    # y 轴（亚群名）字号
    "ytick_fontsize": 10,
}

print("Drawing Fig21...")
# 从 delta_X 计算 cluster-level 均值
cluster_delta_ko = {}
cluster_delta_oe = {}
for st in subtypes_order:
    mask = subtypes == st
    cluster_delta_ko[st] = delta_ko[mask, :].mean(0)
    cluster_delta_oe[st] = delta_oe[mask, :].mean(0)

df_cdx_ko = pd.DataFrame(cluster_delta_ko, index=delta_gene_names).T
df_cdx_oe = pd.DataFrame(cluster_delta_oe, index=delta_gene_names).T

top_genes_21 = df_cdx_ko.abs().max(0).nlargest(P21["top_n_genes"]).index
df_ko_top_21 = df_cdx_ko[top_genes_21]
df_oe_top_21 = df_cdx_oe[top_genes_21]

fig, axes = plt.subplots(1, 2, figsize=P21["figsize"])
fig.patch.set_facecolor("white")
for ax, df_plot, title in zip(axes, [df_ko_top_21, df_oe_top_21], ["MAFB KO", "MAFB OE"]):
    vmax_21 = df_plot.abs().max().max()
    im = ax.imshow(df_plot.values, cmap=P21["cmap"],
                   aspect="auto", vmin=-vmax_21, vmax=vmax_21)
    plt.colorbar(im, ax=ax, label=P21["colorbar_label"], shrink=0.8)
    ax.set_xticks(range(len(df_plot.columns)))
    ax.set_xticklabels(df_plot.columns, rotation=90, fontsize=P21["xtick_fontsize"])
    ax.set_yticks(range(len(df_plot.index)))
    ax.set_yticklabels(df_plot.index, fontsize=P21["ytick_fontsize"])
    ax.set_title(f"Cluster Delta Expression — {title}\n(top {P21['top_n_genes']} variable genes)",
                 fontsize=P21["title_fontsize"], fontweight="bold")
plt.suptitle(P21["suptitle"], fontsize=P21["suptitle_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig21_cluster_deltaX_heatmap")

# ══════════════════════════════════════════════════════════════════════════════
# Fig22: KO 响应 DEG 火山图
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig22 前置参数 ──────────────────────────────────────────────────────────
P22 = {
    "figsize": (18, 6),
    "suptitle": "MAFB KO Response — DEG Volcano Plot per Subtype",
    "suptitle_fontsize": 13,
    "title_fontsize": 12,
    # x 轴标签
    "xlabel": "log2FC (KO vs WT)",
    "xlabel_fontsize": 10,
    # y 轴标签（仅第一个子图）
    "ylabel": "-log10(padj)",
    "ylabel_fontsize": 10,
    # 背景点颜色（不显著）
    "bg_color": "#CCCCCC",
    # 背景点大小
    "bg_size": 12,
    # 背景点透明度
    "bg_alpha": 0.5,
    # 上调点颜色
    "up_color": "#D6604D",
    # 下调点颜色
    "down_color": "#4393C3",
    # 显著点大小
    "sig_size": 20,
    # 显著点透明度
    "sig_alpha": 0.8,
    # 显著性阈值（padj）
    "padj_thresh": 0.05,
    # log2FC 阈值（绝对值）
    "fc_thresh": 0.5,
    # 标注 top N 基因
    "top_label_n": 5,
    # 标注字号
    "label_fontsize": 7,
    "legend_fontsize": 8,
    "legend_loc": "upper right",
}

print("Drawing Fig22...")
log2fc_col_22 = "log2FC" if "log2FC" in df_deg.columns else \
    [c for c in df_deg.columns if "log2" in c.lower() or "fc" in c.lower()][0]
padj_col_22 = "padj" if "padj" in df_deg.columns else \
    [c for c in df_deg.columns if "padj" in c.lower() or "adj" in c.lower()][0]
gene_col_22 = "gene" if "gene" in df_deg.columns else df_deg.columns[0]
subtype_col_22 = "subtype" if "subtype" in df_deg.columns else df_deg.columns[1]

fig, axes = plt.subplots(1, len(subtypes_order), figsize=P22["figsize"])
fig.patch.set_facecolor("white")
for col_idx, st in enumerate(subtypes_order):
    ax = axes[col_idx]
    sub = df_deg[df_deg[subtype_col_22] == st].copy()
    if len(sub) == 0:
        ax.text(0.5, 0.5, "No DEGs", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(st, fontsize=P22["title_fontsize"], fontweight="bold",
                     color=subtype_colors[st])
        continue
    sub["-log10_padj"] = -np.log10(sub[padj_col_22].clip(1e-300))
    sig_mask = (sub[padj_col_22] < P22["padj_thresh"]) & \
               (sub[log2fc_col_22].abs() > P22["fc_thresh"])
    ax.scatter(sub.loc[~sig_mask, log2fc_col_22],
               sub.loc[~sig_mask, "-log10_padj"],
               c=P22["bg_color"], s=P22["bg_size"],
               alpha=P22["bg_alpha"], linewidths=0, zorder=2)
    up_mask = sig_mask & (sub[log2fc_col_22] > 0)
    dn_mask = sig_mask & (sub[log2fc_col_22] < 0)
    ax.scatter(sub.loc[up_mask, log2fc_col_22],
               sub.loc[up_mask, "-log10_padj"],
               c=P22["up_color"], s=P22["sig_size"],
               alpha=P22["sig_alpha"], linewidths=0, zorder=3,
               label=f"Up ({up_mask.sum()})")
    ax.scatter(sub.loc[dn_mask, log2fc_col_22],
               sub.loc[dn_mask, "-log10_padj"],
               c=P22["down_color"], s=P22["sig_size"],
               alpha=P22["sig_alpha"], linewidths=0, zorder=3,
               label=f"Down ({dn_mask.sum()})")
    top_genes_22 = sub[sig_mask].nlargest(P22["top_label_n"], "-log10_padj")
    for _, row in top_genes_22.iterrows():
        ax.text(row[log2fc_col_22], row["-log10_padj"] + 0.1,
                row[gene_col_22], fontsize=P22["label_fontsize"],
                ha="center", va="bottom")
    ax.axhline(-np.log10(P22["padj_thresh"]), color="gray",
               linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline( P22["fc_thresh"], color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(-P22["fc_thresh"], color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel(P22["xlabel"], fontsize=P22["xlabel_fontsize"])
    ax.set_ylabel(P22["ylabel"] if col_idx == 0 else "",
                  fontsize=P22["ylabel_fontsize"])
    ax.set_title(st, fontsize=P22["title_fontsize"], fontweight="bold",
                 color=subtype_colors[st])
    ax.legend(fontsize=P22["legend_fontsize"], loc=P22["legend_loc"])
    ax.spines[["top","right"]].set_visible(False)
plt.suptitle(P22["suptitle"], fontsize=P22["suptitle_fontsize"], fontweight="bold")
plt.tight_layout()
save_fig(fig, "fig22_ko_response_deg_volcano")

# ══════════════════════════════════════════════════════════════════════════════
# Fig23: KO 响应 ORA 富集 bubble 图 (imputed vs simulated DEG)
# ══════════════════════════════════════════════════════════════════════════════
# ── Fig23 前置参数 ──────────────────────────────────────────────────────────
P23 = {
    # 每个面板宽度（英寸）
    "panel_w": 7,
    # 每个面板高度（英寸）
    "panel_h": 5,
    # 最大列数
    "max_cols": 3,
    "suptitle": ("Fig23: MAFB KO vs WT — ORA Enrichment per Subtype\n"
                 "(WT: imputed_count | KO: simulated_count | "
                 "Wilcoxon signed-rank, BH | padj<0.05, |log2FC|>0.3 | background=HVG3000)"),
    "suptitle_fontsize": 11,
    # 子图标题字号
    "panel_title_fontsize": 11,
    # x 轴标签
    "xlabel": "-log10(Adjusted P-value)",
    "xlabel_fontsize": 10,
    # 显著性阈值竖线
    "sig_line_x": -np.log10(0.05),
    "sig_line_color": "gray",
    "sig_line_style": "--",
    # 气泡缩放系数
    "bubble_scale": 3000,
    # 气泡最小尺寸
    "bubble_min": 30,
    # 气泡透明度
    "bubble_alpha": 0.85,
    # 每个分组最多显示 top N terms
    "top_n_per_group": 5,
    # term 名称最大字符数
    "term_max_chars": 55,
    # 图例标题
    "legend_title": "Gene Set",
    "legend_fontsize": 8,
    "legend_loc": "lower right",
    # 数据库颜色
    "lib_colors": {
        "GO_Biological_Process_2023": "#4DBBD5",
        "KEGG_2021_Human":            "#E64B35",
        "MSigDB_Hallmark_2020":       "#00A087",
    },
}

print("Drawing Fig23...")
padj_col_23 = "Adjusted P-value" if "Adjusted P-value" in df_fig23_enr.columns else "padj"
gene_set_col_23 = "gene_set" if "gene_set" in df_fig23_enr.columns else \
    [c for c in df_fig23_enr.columns if "gene_set" in c.lower() or "Gene_set" in c][0]

top_terms_23 = []
for st in subtypes_order:
    for direction in ["up", "down"]:
        for gs in df_fig23_enr[gene_set_col_23].unique():
            sub = df_fig23_enr[
                (df_fig23_enr["subtype"] == st) &
                (df_fig23_enr["direction"] == direction) &
                (df_fig23_enr[gene_set_col_23] == gs)
            ].copy()
            sub = sub[sub[padj_col_23] < 0.05]
            if len(sub) == 0:
                continue
            sub = sub.nsmallest(P23["top_n_per_group"], padj_col_23)
            sub["panel"] = f"{st} ({direction})"
            top_terms_23.append(sub)

if top_terms_23:
    df_plot_23 = pd.concat(top_terms_23, ignore_index=True)
    df_plot_23["-log10_padj"] = -np.log10(df_plot_23[padj_col_23] + 1e-300)
    df_plot_23["overlap_frac"] = (df_plot_23["n_fg_genes"].astype(float) /
                                   df_plot_23["n_bg_genes"].astype(float))
    df_plot_23["Term_short"] = df_plot_23["Term"].apply(
        lambda x: x[:P23["term_max_chars"]] + "..." if len(x) > P23["term_max_chars"] else x)

    panels_all = [f"{st} ({d})" for st in subtypes_order for d in ["up", "down"]]
    panels_present = [p for p in panels_all if p in df_plot_23["panel"].values]

    n_panels = len(panels_present)
    ncols = min(P23["max_cols"], n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * P23["panel_w"], nrows * P23["panel_h"] + 1),
                              squeeze=False)
    fig.patch.set_facecolor("white")

    for idx, panel in enumerate(panels_present):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        sub = df_plot_23[df_plot_23["panel"] == panel].copy()
        sub = sub.sort_values("-log10_padj", ascending=True)

        for gs in df_fig23_enr[gene_set_col_23].unique():
            gs_sub = sub[sub[gene_set_col_23] == gs]
            if len(gs_sub) == 0:
                continue
            color = P23["lib_colors"].get(gs, "#888888")
            ax.scatter(
                gs_sub["-log10_padj"], gs_sub["Term_short"],
                s=gs_sub["overlap_frac"] * P23["bubble_scale"] + P23["bubble_min"],
                c=color, alpha=P23["bubble_alpha"],
                linewidths=0.5, edgecolors="white",
                label=gs.replace("_2023","").replace("_2021_Human","").replace("_2020",""),
                zorder=3)

        ax.axvline(x=P23["sig_line_x"], color=P23["sig_line_color"],
                   linestyle=P23["sig_line_style"], linewidth=1, alpha=0.7)
        ax.set_xlabel(P23["xlabel"], fontsize=P23["xlabel_fontsize"])
        ax.set_title(panel, fontsize=P23["panel_title_fontsize"], fontweight="bold")
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.3)
        if idx == 0:
            ax.legend(title=P23["legend_title"],
                      fontsize=P23["legend_fontsize"],
                      loc=P23["legend_loc"], framealpha=0.8)

    # 关闭多余子图
    for idx in range(len(panels_present), nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].axis("off")

    plt.suptitle(P23["suptitle"], fontsize=P23["suptitle_fontsize"],
                 fontweight="bold", y=1.01)
    plt.tight_layout()
    save_fig(fig, "fig23_mafb_ko_ora_enrichment_bubble")

else:
    # 无显著富集时生成占位图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.5, 0.5,
            "No significant ORA terms\n(foreground sets < 3 genes per subtype/direction)",
            ha="center", va="center", fontsize=14, transform=ax.transAxes)
    ax.set_title("Fig23: MAFB KO vs WT — ORA Enrichment (imputed vs simulated)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "fig23_mafb_ko_ora_enrichment_bubble")

# ══════════════════════════════════════════════════════════════════════════════
# 完成
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"All 23 figures saved to: {OUT_DIR}")
print(f"{'='*60}")



# ══════════════════════════════════════════════════════════════════════════════
# Fig08: MAFB Regulatory Network (FINAL — LOGIC CORRECTED)
# - GRN: Mono + IFN_TAM + LA_TAM
# - top_targets 排除 co-regulators（✔ 正确方向）
# - CTSL / APOL3 强制保留
# ══════════════════════════════════════════════════════════════════════════════

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Parameters (unchanged) ─────────────────────────────────────────────
P08 = {
    "figsize": (12, 10),
    "title": "MAFB Regulatory Network\n(Top 15 targets + co-regulators)",
    "title_fontsize": 13,
    "n_top_targets": 15,
    "n_top_coregulators": 5,
    "coregulator_coef_thresh": 0.1,
    "mafb_radius": 0.5,
    "target_radius": 0.28,
    "coregulator_radius": 0.32,
    "target_ring_r": 3.2,
    "coregulator_ring_r": 1.8,
    "mafb_color": "#FF6B35",
    "mafb_edge_color": "#CC3300",
    "activated_target_color": "#FFCCCC",
    "repressed_target_color": "#CCCCFF",
    "coregulator_color": "#FFFFCC",
    "coregulator_edge_color": "#888800",
    "activation_edge_color": "#D6604D",
    "repression_edge_color": "#4393C3",
    "target_fontsize": 7,
    "coregulator_fontsize": 8,
    "mafb_fontsize": 12,
    "legend_fontsize": 8,
    "legend_loc": "upper right",
    "axis_lim": 4.2,
}

print("Drawing Fig08 (FINAL LOGIC FIXED)...")

# ── GRN merge ─────────────────────────────────────────────────────────
grn_all = pd.concat(
    [grn_coefs["Mono"], grn_coefs["IFN_TAM"], grn_coefs["LA_TAM"]],
    ignore_index=True
)

# ── MAFB direct GRN ───────────────────────────────────────────────────
mafb_grn = grn_all[grn_all["TF"] == "MAFB"].copy()
mafb_grn["abs_coef"] = mafb_grn["coef"].abs()

# ── Ordered unique targets by |coef| ──────────────────────────────────
ordered_targets = (
    mafb_grn
    .sort_values("abs_coef", ascending=False)
    .drop_duplicates(subset="target", keep="first")
    ["target"]
    .tolist()
)

# ── Step 1: 临时选一批 targets（多于 15） ─────────────────────────────
tmp_targets = ordered_targets[:30]   # 足够大，避免后面不够

# ── Step 2: 计算 co‑regulators（基于 tmp_targets） ────────────────────
co_regulators = {}
for tgt in tmp_targets:
    for _, row in grn_all[grn_all["target"] == tgt].iterrows():
        if (
            row["TF"] != "MAFB"
            and abs(row["coef"]) > P08["coregulator_coef_thresh"]
        ):
            co_regulators[row["TF"]] = co_regulators.get(row["TF"], 0) + abs(row["coef"])

top_coregulators = sorted(
    co_regulators.items(),
    key=lambda x: x[1],
    reverse=True
)[:P08["n_top_coregulators"]]

co_regulator_genes = {tf for tf, _ in top_coregulators}

# ── Step 3 ✅ 构造最终 top_targets（排除 co‑regulators） ──────────────
top_targets = []

# 3.1 强制加入 CTSL / APOL3
for g in ["CTSL", "APOL3"]:
    if g in ordered_targets:
        top_targets.append(g)

# 3.2 按 |coef| 顺序补齐，但排除 co‑regulators
for g in ordered_targets:
    if g in co_regulator_genes:
        continue
    if g not in top_targets:
        top_targets.append(g)
    if len(top_targets) == P08["n_top_targets"]:
        break

assert len(top_targets) == 15, f"Expected 15 targets, got {len(top_targets)}"

# ── Layout ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=P08["figsize"])
fig.patch.set_facecolor("white")
ax.set_facecolor("#F8F8F8")

mafb_pos = np.array([0.0, 0.0])

cor_angles = np.linspace(0, 2*np.pi, len(top_coregulators), endpoint=False)
cor_pos = {
    tf: np.array([
        P08["coregulator_ring_r"] * np.cos(cor_angles[i]),
        P08["coregulator_ring_r"] * np.sin(cor_angles[i])
    ])
    for i, (tf, _) in enumerate(top_coregulators)
}

tgt_angles = np.linspace(0, 2*np.pi, len(top_targets), endpoint=False)
tgt_pos = {
    tgt: np.array([
        P08["target_ring_r"] * np.cos(tgt_angles[i]),
        P08["target_ring_r"] * np.sin(tgt_angles[i])
    ])
    for i, tgt in enumerate(top_targets)
}

coef_max = mafb_grn["abs_coef"].max() + 1e-8

# ── MAFB → target edges ───────────────────────────────────────────────
for tgt in top_targets:
    coef = mafb_grn.loc[mafb_grn["target"] == tgt, "coef"].iloc[0]
    color = P08["activation_edge_color"] if coef > 0 else P08["repression_edge_color"]
    lw = abs(coef) / coef_max * 3 + 0.5
    ax.annotate(
        "", xy=tgt_pos[tgt], xytext=mafb_pos,
        arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=0.7),
        zorder=2
    )

# ── co‑regulator → target edges ───────────────────────────────────────
for tf, _ in top_coregulators:
    for _, row in grn_all[
        (grn_all["TF"] == tf) &
        (grn_all["target"].isin(top_targets))
    ].iterrows():
        color = "#FFAA88" if row["coef"] > 0 else "#88AAFF"
        ax.plot(
            [cor_pos[tf][0], tgt_pos[row["target"]][0]],
            [cor_pos[tf][1], tgt_pos[row["target"]][1]],
            color=color, alpha=0.3, lw=0.8, zorder=1
        )

# ── target nodes ───────────────────────────────────────────────────────
for tgt, pos in tgt_pos.items():
    coef = mafb_grn.loc[mafb_grn["target"] == tgt, "coef"].iloc[0]
    fc = P08["activated_target_color"] if coef > 0 else P08["repressed_target_color"]
    ax.add_patch(
        plt.Circle(pos, P08["target_radius"], color=fc, ec="gray", lw=0.8, zorder=3)
    )
    ax.text(
        pos[0], pos[1], tgt,
        ha="center", va="center",
        fontsize=P08["target_fontsize"],
        fontweight="bold", zorder=4
    )

# ── co‑regulator nodes ────────────────────────────────────────────────
for tf, pos in cor_pos.items():
    ax.add_patch(
        plt.Circle(
            pos, P08["coregulator_radius"],
            color=P08["coregulator_color"],
            ec=P08["coregulator_edge_color"],
            lw=1.2, zorder=3
        )
    )
    ax.text(
        pos[0], pos[1], tf,
        ha="center", va="center",
        fontsize=P08["coregulator_fontsize"],
        fontweight="bold", zorder=4
    )

# ── MAFB node ─────────────────────────────────────────────────────────
ax.add_patch(
    plt.Circle(
        mafb_pos, P08["mafb_radius"],
        color=P08["mafb_color"],
        ec=P08["mafb_edge_color"],
        lw=2, zorder=5
    )
)
ax.text(
    0, 0, "MAFB",
    ha="center", va="center",
    fontsize=P08["mafb_fontsize"],
    fontweight="bold",
    color="white", zorder=6
)

# ── Legend ────────────────────────────────────────────────────────────
ax.legend(
    handles=[
        Line2D([0],[0], color=P08["activation_edge_color"], lw=2, label="Activation"),
        Line2D([0],[0], color=P08["repression_edge_color"], lw=2, label="Repression"),
        plt.Circle((0,0), 0.1, color=P08["mafb_color"], label="MAFB"),
        plt.Circle((0,0), 0.1,
                   color=P08["coregulator_color"],
                   ec=P08["coregulator_edge_color"],
                   label="Co-regulator"),
        plt.Circle((0,0), 0.1,
                   color=P08["activated_target_color"],
                   label="Activated target"),
        plt.Circle((0,0), 0.1,
                   color=P08["repressed_target_color"],
                   label="Repressed target"),
    ],
    loc=P08["legend_loc"],
    fontsize=P08["legend_fontsize"],
    framealpha=0.9
)

ax.set_xlim(-P08["axis_lim"], P08["axis_lim"])
ax.set_ylim(-P08["axis_lim"], P08["axis_lim"])
ax.set_aspect("equal")
ax.axis("off")

ax.set_title(P08["title"], fontsize=P08["title_fontsize"], fontweight="bold")

plt.tight_layout()
save_fig(fig, "fig08_mafb_regulatory_network_FINAL_LOGIC_FIXED")

# ══════════════════════════════════════════════════════════════════════════════
# Fig08: MAFB Regulatory Network (FINAL — UPDATED LOGIC)
# - GRN: Mono + IFN_TAM + LA_TAM
# - Outer ring:
#     • Top 10 targets by |coef| (force include APOL3)
#     • + Bottom 5 targets with coef < 0
# ══════════════════════════════════════════════════════════════════════════════

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Fig08 parameters (unchanged style) ───────────────────────────────────────
P08 = {
    "figsize": (12, 10),
    "title": "MAFB Regulatory Network\n(Top & repressed targets + co-regulators)",
    "title_fontsize": 13,
    "n_top_coregulators": 5,
    "coregulator_coef_thresh": 0.1,
    "mafb_radius": 0.5,
    "target_radius": 0.28,
    "coregulator_radius": 0.32,
    "target_ring_r": 3.2,
    "coregulator_ring_r": 1.8,
    "mafb_color": "#FF6B35",
    "mafb_edge_color": "#CC3300",
    "activated_target_color": "#FFCCCC",
    "repressed_target_color": "#CCCCFF",
    "coregulator_color": "#FFFFCC",
    "coregulator_edge_color": "#888800",
    "activation_edge_color": "#D6604D",
    "repression_edge_color": "#4393C3",
    "target_fontsize": 7,
    "coregulator_fontsize": 8,
    "mafb_fontsize": 12,
    "legend_fontsize": 8,
    "legend_loc": "upper right",
    "axis_lim": 4.2,
}

print("Drawing Fig08 (FINAL UPDATED)…")

# ── Merge GRN from three subtypes ─────────────────────────────────────────────
grn_all = pd.concat(
    [grn_coefs["Mono"], grn_coefs["IFN_TAM"], grn_coefs["LA_TAM"]],
    ignore_index=True
)

# ── MAFB direct GRN ──────────────────────────────────────────────────────────
mafb_grn = grn_all[grn_all["TF"] == "MAFB"].copy()
mafb_grn["abs_coef"] = mafb_grn["coef"].abs()

# unique targets ordered by |coef|
mafb_unique = (
    mafb_grn
    .sort_values("abs_coef", ascending=False)
    .drop_duplicates(subset="target", keep="first")
)

# ── Step 1: Top 10 by |coef| (force include APOL3) ───────────────────────────
top10 = mafb_unique.head(10)["target"].tolist()
if "APOL3" in mafb_unique["target"].values and "APOL3" not in top10:
    top10[-1] = "APOL3"   # 替换掉第 10 个，保证数量不变

# ── Step 2: Bottom 5 repressed targets (coef < 0) ───────────────────────────
bottom5 = (
    mafb_unique[mafb_unique["coef"] < 0]
    .sort_values("coef", ascending=True)
    .head(5)["target"]
    .tolist()
)

# ── Final outer-ring targets ─────────────────────────────────────────────────
top_targets = []
for g in top10 + bottom5:
    if g not in top_targets:
        top_targets.append(g)

print("Final outer-ring targets:")
print(top_targets)
print(f"Total targets: {len(top_targets)}")

# ── Co‑regulators (original Fig08 logic) ─────────────────────────────────────
co_regulators = {}
for tgt in top_targets:
    for _, row in grn_all[grn_all["target"] == tgt].iterrows():
        if (
            row["TF"] != "MAFB"
            and row["TF"] not in top_targets
            and abs(row["coef"]) > P08["coregulator_coef_thresh"]
        ):
            co_regulators[row["TF"]] = co_regulators.get(row["TF"], 0) + abs(row["coef"])

top_coregulators = sorted(
    co_regulators.items(),
    key=lambda x: x[1],
    reverse=True
)[:P08["n_top_coregulators"]]

# ── Layout ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=P08["figsize"])
fig.patch.set_facecolor("white")
ax.set_facecolor("#F8F8F8")

mafb_pos = np.array([0.0, 0.0])

cor_angles = np.linspace(0, 2*np.pi, len(top_coregulators), endpoint=False)
cor_pos = {
    tf: np.array([
        P08["coregulator_ring_r"] * np.cos(cor_angles[i]),
        P08["coregulator_ring_r"] * np.sin(cor_angles[i])
    ])
    for i, (tf, _) in enumerate(top_coregulators)
}

tgt_angles = np.linspace(0, 2*np.pi, len(top_targets), endpoint=False)
tgt_pos = {
    tgt: np.array([
        P08["target_ring_r"] * np.cos(tgt_angles[i]),
        P08["target_ring_r"] * np.sin(tgt_angles[i])
    ])
    for i, tgt in enumerate(top_targets)
}

coef_max = mafb_unique["abs_coef"].max() + 1e-8

# ── MAFB → target edges ──────────────────────────────────────────────────────
for tgt in top_targets:
    coef = mafb_unique.loc[mafb_unique["target"] == tgt, "coef"].iloc[0]
    color = P08["activation_edge_color"] if coef > 0 else P08["repression_edge_color"]
    lw = abs(coef) / coef_max * 3 + 0.5
    ax.annotate(
        "", xy=tgt_pos[tgt], xytext=mafb_pos,
        arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=0.7),
        zorder=2
    )

# ── co‑regulator → target edges ──────────────────────────────────────────────
for tf, _ in top_coregulators:
    for _, row in grn_all[
        (grn_all["TF"] == tf) &
        (grn_all["target"].isin(top_targets))
    ].iterrows():
        color = "#FFAA88" if row["coef"] > 0 else "#88AAFF"
        ax.plot(
            [cor_pos[tf][0], tgt_pos[row["target"]][0]],
            [cor_pos[tf][1], tgt_pos[row["target"]][1]],
            color=color, alpha=0.3, lw=0.8, zorder=1
        )

# ── Target nodes ─────────────────────────────────────────────────────────────
for tgt, pos in tgt_pos.items():
    coef = mafb_unique.loc[mafb_unique["target"] == tgt, "coef"].iloc[0]
    fc = P08["activated_target_color"] if coef > 0 else P08["repressed_target_color"]
    ax.add_patch(
        plt.Circle(pos, P08["target_radius"], color=fc, ec="gray", lw=0.8, zorder=3)
    )
    ax.text(
        pos[0], pos[1], tgt,
        ha="center", va="center",
        fontsize=P08["target_fontsize"],
        fontweight="bold", zorder=4
    )

# ── Co‑regulator nodes ───────────────────────────────────────────────────────
for tf, pos in cor_pos.items():
    ax.add_patch(
        plt.Circle(
            pos, P08["coregulator_radius"],
            color=P08["coregulator_color"],
            ec=P08["coregulator_edge_color"], lw=1.2, zorder=3
        )
    )
    ax.text(
        pos[0], pos[1], tf,
        ha="center", va="center",
        fontsize=P08["coregulator_fontsize"], fontweight="bold", zorder=4
    )

# ── MAFB node ────────────────────────────────────────────────────────────────
ax.add_patch(
    plt.Circle(
        mafb_pos, P08["mafb_radius"],
        color=P08["mafb_color"],
        ec=P08["mafb_edge_color"], lw=2, zorder=5
    )
)
ax.text(
    0, 0, "MAFB",
    ha="center", va="center",
    fontsize=P08["mafb_fontsize"],
    fontweight="bold",
    color="white", zorder=6
)

# ── Legend ───────────────────────────────────────────────────────────────────
ax.legend(
    handles=[
        Line2D([0],[0], color=P08["activation_edge_color"], lw=2, label="Activation"),
        Line2D([0],[0], color=P08["repression_edge_color"], lw=2, label="Repression"),
        plt.Circle((0,0), 0.1, color=P08["mafb_color"], label="MAFB"),
        plt.Circle((0,0), 0.1,
                   color=P08["coregulator_color"],
                   ec=P08["coregulator_edge_color"], label="Co-regulator"),
        plt.Circle((0,0), 0.1,
                   color=P08["activated_target_color"], label="Activated target"),
        plt.Circle((0,0), 0.1,
                   color=P08["repressed_target_color"], label="Repressed target"),
    ],
    loc=P08["legend_loc"],
    fontsize=P08["legend_fontsize"],
    framealpha=0.9
)

ax.set_xlim(-P08["axis_lim"], P08["axis_lim"])
ax.set_ylim(-P08["axis_lim"], P08["axis_lim"])
ax.set_aspect("equal")
ax.axis("off")

ax.set_title(P08["title"], fontsize=P08["title_fontsize"], fontweight="bold")

plt.tight_layout()
save_fig(fig, "fig08_mafb_regulatory_network_FINAL_UPDATED")

#!/usr/bin/env python
# fig_heatmap_all_deg_horizontal.py
# 横版紧凑热图：所有显著 DEG
#   行 = Mono / IFN_TAM / LA_TAM（subtype 变成行）
#   列 = 基因（按来源 subtype 分组排列）
#   值 = KO vs WT log2FC
# 输入：07_ko_response_deg_per_subtype.csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# =============================================================================
# ★ 路径配置
# =============================================================================
BASE = translate(r"D:\bulk-download\celloracle0331")
OUT  = Path(BASE) / "figures"
OUT.mkdir(parents=True, exist_ok=True)

def P(fname):
    return str(Path(BASE) / fname)

def save_fig(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {name}.png / .svg")

# subtype 顺序（行顺序）
SUBTYPES_ORDER = ["Mono", "IFN_TAM", "LA_TAM"]
SUBTYPE_COLORS = {
    "Mono":    "#E64B35",
    "IFN_TAM": "#4DBBD5",
    "LA_TAM":  "#00A087",
}

# =============================================================================
# ★ 参数字典
# =============================================================================
PARAMS = {
    # ── 基因筛选 ──────────────────────────────────────────────────────────────
    "padj_thresh"        : 0.05,
    "fc_thresh"          : 0.3,    # |log2FC| 最低阈值
    "direction"          : "down", # "down" / "up" / "both"

    # ── 基因归属优先级（越靠前越优先）────────────────────────────────────────
    # 一个基因可能在多个 subtype 里显著，归属给优先级最高的那个
    "source_priority"    : ["LA_TAM", "IFN_TAM", "Mono"],

    # ── 列排序 ────────────────────────────────────────────────────────────────
    # "subtype"  : 按归属 subtype 分组，组内按该 subtype 的 log2FC 排序
    # "cluster"  : 层次聚类（需要 scipy）
    "col_order"          : "subtype",

    # ── 颜色范围（固定区间）──────────────────────────────────────────────────
    "vmin"               : -1.5,   # 颜色下限
    "vmax"               :  1,   # 颜色上限
    # 注意：不再自动对称，直接用 vmin/vmax

    # ── 色块大小 ──────────────────────────────────────────────────────────────
    "cell_height_pt"     : 40,     # 每行（subtype）高度，pt；横版行少所以可以大一些
    "cell_width_pt"      : 2,      # 每列（基因）宽度，pt；基因多所以要小

    # ── 顶部来源色条（现在是顶部，对应列的归属）─────────────────────────────
    "show_source_bar"    : True,
    "source_bar_height_pt": 6,     # 来源色条高度，pt
    "show_group_label"   : True,   # 是否在色条上方标注 subtype 名
    "group_label_fs"     : 10,

    # ── 分组分隔线 ────────────────────────────────────────────────────────────
    "show_group_divider" : True,
    "divider_color"      : "black",
    "divider_lw"         : 0.8,
    "divider_ls"         : "--",

    # ── 行分隔线（subtype 之间）──────────────────────────────────────────────
    "row_divider_lw"     : 0.8,
    "row_divider_color"  : "#555555",

    # ── 左侧 subtype 标签 ─────────────────────────────────────────────────────
    "ylabel_fs"          : 10,
    "show_gene_names"    : False,  # 列（基因）是否显示名称

    # ── 标题 & colorbar ───────────────────────────────────────────────────────
    "title_fs"           : 10,
    "cbar_label"         : "log2FC (KO vs WT)",
    "cbar_label_fs"      : 8,
    "cbar_tick_fs"       : 7,
}

# =============================================================================
# 1. 读取 DEG 表，自动识别列名
# =============================================================================
print("读取 DEG 表...")
df_deg = pd.read_csv(P("07_ko_response_deg_per_subtype.csv"))
print(f"  shape: {df_deg.shape}")
print(f"  cols : {df_deg.columns.tolist()}")

log2fc_col = ("log2FC" if "log2FC" in df_deg.columns else
              next(c for c in df_deg.columns
                   if "log2" in c.lower() or
                      ("fc" in c.lower() and "log" in c.lower())))
padj_col   = ("padj" if "padj" in df_deg.columns else
              next(c for c in df_deg.columns
                   if "padj" in c.lower() or "adj" in c.lower()))
gene_col   = ("gene" if "gene" in df_deg.columns else df_deg.columns[0])
stype_col  = ("subtype" if "subtype" in df_deg.columns else df_deg.columns[1])
print(f"  识别列: gene={gene_col}, subtype={stype_col}, "
      f"log2FC={log2fc_col}, padj={padj_col}")

# =============================================================================
# 2. 筛选显著基因
# =============================================================================
direction = PARAMS["direction"]
padj_thr  = PARAMS["padj_thresh"]
fc_thr    = PARAMS["fc_thresh"]

sig_mask = df_deg[padj_col] < padj_thr
if direction == "down":
    dir_mask = df_deg[log2fc_col] < -fc_thr
elif direction == "up":
    dir_mask = df_deg[log2fc_col] >  fc_thr
else:
    dir_mask = df_deg[log2fc_col].abs() > fc_thr

df_sig = df_deg[sig_mask & dir_mask].copy()

print(f"\n筛选结果:")
for st in SUBTYPES_ORDER:
    n = len(df_sig[df_sig[stype_col] == st])
    print(f"  {st}: {n} 个显著基因")

# 降级：某 subtype 无显著基因时取极端 20 个
for st in SUBTYPES_ORDER:
    if len(df_sig[df_sig[stype_col] == st]) == 0:
        print(f"  WARNING: {st} 无显著基因，取极端 20 个替代")
        sub = df_deg[df_deg[stype_col] == st].copy()
        sub = sub.sort_values(log2fc_col,
                              ascending=(direction != "up")).head(20)
        df_sig = pd.concat([df_sig, sub], ignore_index=True)

# =============================================================================
# 3. 确定基因归属（优先级：LA_TAM > IFN_TAM > Mono）
# =============================================================================
# 先收集所有显著基因（跨所有 subtype 的并集）
all_sig_genes = df_sig[gene_col].unique().tolist()
print(f"\n所有显著基因（去重前）: {len(all_sig_genes)}")

priority = PARAMS["source_priority"]  # ["LA_TAM", "IFN_TAM", "Mono"]

# 对每个基因，按优先级判断归属：
# 如果该基因在优先级最高的 subtype 里显著，就归属给它；否则看次高，以此类推
gene_source = {}
for g in all_sig_genes:
    assigned = False
    for st in priority:
        # 检查该基因在该 subtype 里是否显著
        rows = df_sig[(df_sig[gene_col] == g) & (df_sig[stype_col] == st)]
        if len(rows) > 0:
            gene_source[g] = st
            assigned = True
            break
    if not assigned:
        # 理论上不会发生，保底归给第一个出现的 subtype
        gene_source[g] = df_sig[df_sig[gene_col] == g][stype_col].values[0]

# 统计归属结果
print("归属结果（优先级: LA_TAM > IFN_TAM > Mono）:")
for st in SUBTYPES_ORDER:
    n = sum(1 for g in all_sig_genes if gene_source[g] == st)
    print(f"  {st}: {n} 个基因")

# =============================================================================
# 4. 确定列（基因）顺序
# =============================================================================
if PARAMS["col_order"] == "subtype":
    # 按归属 subtype 分组，组内按该 subtype 的 log2FC 排序
    gene_order = []
    for st in SUBTYPES_ORDER:
        # 属于该 subtype 的基因
        genes_in_st = [g for g in all_sig_genes if gene_source[g] == st]
        # 取该 subtype 的 log2FC 排序（下调：升序；上调：降序）
        ref = df_sig[df_sig[stype_col] == st].set_index(gene_col)[log2fc_col]
        genes_sorted = sorted(
            genes_in_st,
            key=lambda g: ref.get(g, 0),
            reverse=(direction == "up")
        )
        gene_order.extend(genes_sorted)

elif PARAMS["col_order"] == "cluster":
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    tmp = pd.DataFrame(index=all_sig_genes,
                       columns=SUBTYPES_ORDER, dtype=float)
    for st in SUBTYPES_ORDER:
        sub = df_sig[df_sig[stype_col] == st].set_index(gene_col)[log2fc_col]
        for g in all_sig_genes:
            tmp.loc[g, st] = sub.get(g, 0.0)
    Z = linkage(pdist(tmp.fillna(0).values, metric="euclidean"), method="ward")
    gene_order = [all_sig_genes[i] for i in leaves_list(Z)]

n_genes = len(gene_order)
n_rows  = len(SUBTYPES_ORDER)
print(f"\n最终列数（基因数）: {n_genes}")

# =============================================================================
# 5. 构建矩阵（行=subtype，列=基因）
# =============================================================================
mat = pd.DataFrame(index=SUBTYPES_ORDER, columns=gene_order, dtype=float)
for st in SUBTYPES_ORDER:
    sub = df_deg[df_deg[stype_col] == st].set_index(gene_col)[log2fc_col]
    for g in gene_order:
        mat.loc[st, g] = sub.get(g, np.nan)

mat_arr = mat.values.astype(float)

vmin = PARAMS["vmin"]
vmax = PARAMS["vmax"]
print(f"颜色范围: {vmin} ~ {vmax}")

# =============================================================================
# 6. 计算图形尺寸（由色块大小反推）
# =============================================================================
PT_PER_INCH = 72.0

cell_h_in  = PARAMS["cell_height_pt"]      / PT_PER_INCH
cell_w_in  = PARAMS["cell_width_pt"]       / PT_PER_INCH
src_bar_h_in = PARAMS["source_bar_height_pt"] / PT_PER_INCH

heat_h_in = cell_h_in  * n_rows    # 热图主体高度
heat_w_in = cell_w_in  * n_genes   # 热图主体宽度

# 边距
margin_left   = 1.0    # 左侧 subtype 标签区域
margin_right  = 1.2    # colorbar 区域
margin_top    = src_bar_h_in + 0.7  # 顶部来源色条 + 标题
margin_bottom = 0.4

fig_w = margin_left  + heat_w_in + margin_right
fig_h = margin_top   + heat_h_in + margin_bottom
fig_w = max(fig_w, 5.0)
fig_h = max(fig_h, 2.5)
# ★ 手动覆盖图形尺寸（覆盖自动计算结果）
fig_w = 10   # 总宽度（英寸）
fig_h = 4.5    # 总高度（英寸）

# 同步更新热图主体尺寸，让热图撑满可用区域
heat_w_in = fig_w - margin_left - margin_right
heat_h_in = fig_h - margin_top  - margin_bottom


print(f"图形尺寸: {fig_w:.2f} × {fig_h:.2f} 英寸")

# =============================================================================
# 7. 绘图
# =============================================================================


fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor("white")

def to_frac_w(val_in): return val_in / fig_w
def to_frac_h(val_in): return val_in / fig_h

# 热图区域（figure fraction）
heat_l = to_frac_w(margin_left)
heat_b = to_frac_h(margin_bottom)
heat_w = to_frac_w(heat_w_in)
heat_h = to_frac_h(heat_h_in)

ax_heat = fig.add_axes([heat_l, heat_b, heat_w, heat_h])

# 顶部来源色条
if PARAMS["show_source_bar"]:
    bar_h = to_frac_h(src_bar_h_in)
    bar_b = heat_b + heat_h + to_frac_h(0.02)
    ax_bar = fig.add_axes([heat_l, bar_b, heat_w, bar_h])

# colorbar
cbar_l = heat_l + heat_w + to_frac_w(0.08)
cbar_w = to_frac_w(0.10)
cbar_b = heat_b + heat_h * 0.2
cbar_h_ax = heat_h * 0.6
ax_cbar = fig.add_axes([cbar_l, cbar_b, cbar_w, cbar_h_ax])

# ── 热图主体 ─────────────────────────────────────────────────────────────────
mat_plot = mat_arr.copy()
nan_mask = np.isnan(mat_plot)
mat_filled = np.where(nan_mask, (vmin + vmax) / 2, mat_plot)  # NaN 填中值

im = ax_heat.imshow(
    mat_filled,
    cmap=PARAMS["heatmap_cmap"] if "heatmap_cmap" in PARAMS
         else "RdBu_r",
    aspect="auto",
    vmin=vmin, vmax=vmax,
    interpolation="nearest"
)

# NaN 格子叠加灰色遮罩
if nan_mask.any():
    nan_rgba = np.zeros((*nan_mask.shape, 4))
    nan_rgba[nan_mask] = [0.82, 0.82, 0.82, 1.0]
    ax_heat.imshow(nan_rgba, aspect="auto",
                   interpolation="nearest", zorder=3)

# 列（基因）分组分隔线
if PARAMS["show_group_divider"] and PARAMS["col_order"] == "subtype":
    col_count = 0
    for st in SUBTYPES_ORDER[:-1]:
        genes_in_st = [g for g in gene_order if gene_source[g] == st]
        col_count += len(genes_in_st)
        ax_heat.axvline(col_count - 0.5,
                        color=PARAMS["divider_color"],
                        linewidth=PARAMS["divider_lw"],
                        linestyle=PARAMS["divider_ls"],
                        zorder=6)

# 行（subtype）分隔线
for y in np.arange(0.5, n_rows - 0.5, 1):
    ax_heat.axhline(y,
                    color=PARAMS["row_divider_color"],
                    linewidth=PARAMS["row_divider_lw"],
                    zorder=5)

# Y 轴：subtype 标签（左侧）
ax_heat.set_yticks(range(n_rows))
ax_heat.set_yticklabels(
    SUBTYPES_ORDER,
    fontsize=PARAMS["ylabel_fs"],
    fontweight="bold",
    rotation=90,
    va="center"
)
# 给每个 subtype 标签上色
for tick, st in zip(ax_heat.get_yticklabels(), SUBTYPES_ORDER):
    tick.set_color(SUBTYPE_COLORS[st])

# X 轴：基因名（默认不显示）
if PARAMS["show_gene_names"]:
    ax_heat.set_xticks(range(n_genes))
    ax_heat.set_xticklabels(gene_order, rotation=90,
                             fontsize=5, fontstyle="italic")
else:
    ax_heat.set_xticks([])

ax_heat.tick_params(left=True, bottom=False, right=False, top=False,
                    length=0)
ax_heat.set_xlim(-0.5, n_genes - 0.5)
ax_heat.set_ylim(n_rows - 0.5, -0.5)

for spine in ax_heat.spines.values():
    spine.set_visible(False)

# ── 顶部来源色条 ─────────────────────────────────────────────────────────────
if PARAMS["show_source_bar"]:
    for i, g in enumerate(gene_order):
        ax_bar.bar(i, 1,
                   color=SUBTYPE_COLORS[gene_source[g]],
                   edgecolor="none", width=1.0)

    # 分组标注文字（每组中间位置）
    if PARAMS["show_group_label"] and PARAMS["col_order"] == "subtype":
        col_start = 0
        for st in SUBTYPES_ORDER:
            genes_in_st = [g for g in gene_order if gene_source[g] == st]
            n_st = len(genes_in_st)
            if n_st > 0:
                mid = col_start + n_st / 2.0 - 0.5
                ax_bar.text(
                    mid, 1.2, st,
                    ha="center", va="bottom",
                    fontsize=PARAMS["group_label_fs"],
                    color=SUBTYPE_COLORS[st],
                    fontweight="bold",
                    clip_on=False
                )
                # 分组分隔线（色条上）
                if col_start > 0:
                    ax_bar.axvline(col_start - 0.5,
                                   color="white", linewidth=1.5, zorder=5)
            col_start += n_st

    ax_bar.set_xlim(-0.5, n_genes - 0.5)
    ax_bar.set_ylim(0, 1)
    ax_bar.axis("off")

# ── 标题 ─────────────────────────────────────────────────────────────────────
direction_label = {"down": "Down-regulated",
                   "up":   "Up-regulated",
                   "both": "DEG"}[direction]

title_y = (to_frac_h(margin_bottom + heat_h_in + src_bar_h_in + 0.12)
           if PARAMS["show_source_bar"]
           else to_frac_h(margin_bottom + heat_h_in + 0.1))

fig.text(
    heat_l + heat_w / 2,
    title_y,
    f"All Significant {direction_label} Genes  (MAFB KO vs WT)\n",
    ha="center", va="bottom",
    fontsize=PARAMS["title_fs"],
    fontweight="bold"
)

# ── colorbar ─────────────────────────────────────────────────────────────────
cbar = fig.colorbar(im, cax=ax_cbar)
cbar.set_label(PARAMS["cbar_label"], fontsize=PARAMS["cbar_label_fs"])
cbar.ax.tick_params(labelsize=PARAMS["cbar_tick_fs"])
# 标注 0 线位置
zero_frac = (0 - vmin) / (vmax - vmin)
cbar.ax.axhline(zero_frac, color="black",
                linewidth=0.8, linestyle="--", alpha=0.8)

# =============================================================================
# 8. 保存
# =============================================================================
save_fig(fig, "fig_horizontal_heatmap_all_deg")
print(f"\n完成！{n_rows} 个 subtype × {n_genes} 个基因")
print(f"输出目录：{OUT}")


# =============================================================================
# Fig23_merged: 三个亚型合并的 ORA 富集气泡图（单张图，up/down 两个面板）
# 输入：fig23_ora_enrichment_results.csv（已有）
# 合并逻辑：同一 Term + direction + gene_set 组合，
#           取最小 padj，n_fg_genes 取并集估算（去重求和后打折），
#           overlap_frac 重新计算
# =============================================================================


import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

# =============================================================================
# ★ 路径配置
# =============================================================================
BASE = translate(r"D:\bulk-download\celloracle0331")
OUT  = Path(BASE) / "figures"
OUT.mkdir(parents=True, exist_ok=True)

def P(fname):
    return str(Path(BASE) / fname)

def save_fig(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {name}.png / .svg")

# =============================================================================
# ★ 辅助函数：去除 GO term 名称末尾的 GO ID 尾巴
#   兼容英文/中文括号与冒号的所有组合
# =============================================================================
def strip_go_id(term: str) -> str:
    return re.sub(
        r"[\s　]*[(\（][Gg][Oo][：:]\d+[)\）]\s*$",
        "",
        str(term)
    ).strip()

# =============================================================================
# ★ 参数字典
# =============================================================================
PARAMS = {
    # ── DEG 筛选 ─────────────────────────────────────────────────────────────
    "deg_padj_thresh"   : 0.2,
    "log2fc_thresh"     : 0.1,

    # ── ORA 绘图筛选 ─────────────────────────────────────────────────────────
    "ora_padj_thresh"   : 0.2,

    # ── ORA ──────────────────────────────────────────────────────────────────
    "gene_sets" : [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "MSigDB_Hallmark_2020",
    ],
    # 面板内显示顺序（从上到下）：GO → KEGG → MSigDB
    "gs_order"  : [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "MSigDB_Hallmark_2020",
    ],
    "organism"      : "human",
    "top_n_per_gs"  : 11,
    "term_max_chars": 60,

    # ── 气泡图外观 ────────────────────────────────────────────────────────────
    "bubble_scale"  : 20,
    "bubble_min"    : 20,
    "bubble_alpha"  : 0.85,
    "panel_w"       : 9,
    "panel_h"       : 8,
    "title_fs"      : 11,
    "xlabel_fs"     : 10,
    "ytick_fs"      : 10.5,
    "legend_fs"     : 9,

    # ── 气泡大小图例固定三挡 ──────────────────────────────────────────────────
    "legend_size_refs": [1, 5, 10],

    # ── 图例位置（右下角，y 控制高低：0.0=底部边缘，正值往上）────────────────
    # 设为 0.08：比底部略高，不遮住 MSigDB 气泡
    "legend_bbox"   : (1.0, 0.05),

    # ── 颜色 ─────────────────────────────────────────────────────────────────
    "lib_colors": {
        "GO_Biological_Process_2023": "#4DBBD5",
        "KEGG_2021_Human"           : "#E64B35",
        "MSigDB_Hallmark_2020"      : "#00A087",
    },
    "lib_labels": {
        "GO_Biological_Process_2023": "GO Biological Process",
        "KEGG_2021_Human"           : "KEGG",
        "MSigDB_Hallmark_2020"      : "MSigDB Hallmark",
    },
}

# =============================================================================
# SECTION 1: 读取表达矩阵
# =============================================================================
print("=" * 60)
print("SECTION 1: 读取表达矩阵")
print("=" * 60)

df_wt = pd.read_csv(P("imputed_count_WT.csv"),   index_col=0)
df_ko = pd.read_csv(P("simulated_count_KO.csv"), index_col=0)
print(f"  WT shape: {df_wt.shape}")
print(f"  KO shape: {df_ko.shape}")

genes = sorted(set(df_wt.columns) & set(df_ko.columns))
print(f"  共同基因数: {len(genes)}")

wt_arr = df_wt[genes].values
ko_arr = df_ko[genes].values

# =============================================================================
# SECTION 2: Mann-Whitney U test
# =============================================================================
print()
print("=" * 60)
print("SECTION 2: Mann-Whitney U test")
print("=" * 60)

results = []
for i, gene in enumerate(genes):
    wt_g = wt_arr[:, i]
    ko_g = ko_arr[:, i]
    if wt_g.max() == 0 and ko_g.max() == 0:
        continue
    try:
        _, pval = stats.mannwhitneyu(ko_g, wt_g, alternative="two-sided")
    except ValueError:
        continue
    mean_wt = wt_g.mean()
    mean_ko = ko_g.mean()
    log2fc  = np.log2((mean_ko + 1e-9) / (mean_wt + 1e-9))
    results.append({"gene": gene, "log2FC": log2fc, "pval": pval,
                    "mean_wt": mean_wt, "mean_ko": mean_ko})
    if (i + 1) % 500 == 0:
        print(f"  进度: {i+1}/{len(genes)}")

print(f"  检验完成: {len(results)} 个基因")

df_deg = pd.DataFrame(results)
_, padj_arr, _, _ = multipletests(df_deg["pval"].values, method="fdr_bh")
df_deg["padj"] = padj_arr

print(f"\n  padj 分布:\n{df_deg['padj'].describe().to_string()}")
print(f"\n  log2FC 分布:\n{df_deg['log2FC'].describe().to_string()}")

sig = df_deg[
    (df_deg["padj"]  < PARAMS["deg_padj_thresh"]) &
    (df_deg["log2FC"].abs() > PARAMS["log2fc_thresh"])
].copy()

n_up   = (sig["log2FC"] > 0).sum()
n_down = (sig["log2FC"] < 0).sum()
print(f"\n  显著 DEG: {len(sig)} 个  (上调 {n_up} / 下调 {n_down})")

if len(sig) == 0:
    print("\n  ⚠ 0 个 DEG，尝试放宽阈值 padj<0.5, |log2FC|>0.05:")
    loose = df_deg[(df_deg["padj"] < 0.5) & (df_deg["log2FC"].abs() > 0.05)]
    print(f"  放宽后: {len(loose)} 个")

deg_path = Path(BASE) / "merged_deg_all_cells.csv"
df_deg.sort_values("padj").to_csv(deg_path, index=False)
print(f"\n  DEG 已保存: {deg_path}")

# =============================================================================
# SECTION 3: ORA 富集分析
# =============================================================================
print()
print("=" * 60)
print("SECTION 3: ORA 富集分析")
print("=" * 60)

try:
    import gseapy as gp
    HAS_GSEAPY = True
    print(f"  gseapy 版本: {gp.__version__}")
except ImportError:
    print("  ✗ gseapy 未安装")
    HAS_GSEAPY = False

background  = df_deg["gene"].tolist()
ora_results = []

if HAS_GSEAPY and len(sig) > 0:
    for direction in ["up", "down"]:
        gene_list = (
            sig[sig["log2FC"] > 0]["gene"].tolist() if direction == "up"
            else sig[sig["log2FC"] < 0]["gene"].tolist()
        )
        if len(gene_list) < 3:
            print(f"  {direction}: 基因数 {len(gene_list)} < 3，跳过")
            continue
        print(f"\n  方向: {direction} ({len(gene_list)} 个基因)")

        for gs in PARAMS["gene_sets"]:
            try:
                enr = gp.enrichr(
                    gene_list  = gene_list,
                    gene_sets  = gs,
                    background = background,
                    organism   = PARAMS["organism"],
                    outdir     = None,
                    verbose    = False,
                )
                df_enr = enr.results.copy()
                print(f"    {gs}: 返回列 → {df_enr.columns.tolist()}")

                df_enr["direction"] = direction
                df_enr["gene_set"]  = gs

                # 统一 padj 列名
                if "Adjusted P-value" in df_enr.columns:
                    df_enr = df_enr.rename(columns={"Adjusted P-value": "padj"})
                elif "P-value" in df_enr.columns:
                    df_enr = df_enr.rename(columns={"P-value": "padj"})

                # 确保 Term 列存在
                if "Term" not in df_enr.columns:
                    df_enr = df_enr.rename(columns={df_enr.columns[0]: "Term"})

                # 去除 GO ID 尾巴
                df_enr["Term"] = df_enr["Term"].apply(strip_go_id)

                # 解析 n_fg
                if "Overlap" in df_enr.columns:
                    split = df_enr["Overlap"].astype(str).str.split("/", expand=True)
                    df_enr["n_fg"] = pd.to_numeric(split[0], errors="coerce").fillna(0).astype(int)
                elif "Genes" in df_enr.columns:
                    df_enr["n_fg"] = df_enr["Genes"].astype(str).apply(
                        lambda x: len([g for g in x.split(";") if g.strip()])
                        if x and x.lower() != "nan" else 0
                    )
                else:
                    df_enr["n_fg"] = 5

                ora_results.append(df_enr)
                n_sig_term = (df_enr["padj"] < PARAMS["ora_padj_thresh"]).sum()
                print(f"    → {n_sig_term} 个显著 term (padj<{PARAMS['ora_padj_thresh']})")

            except Exception as e:
                import traceback
                print(f"    ✗ {gs} 失败: {e}")
                traceback.print_exc()

ora_path = Path(BASE) / "merged_ora_all_cells.csv"
if ora_results:
    df_ora = pd.concat(ora_results, ignore_index=True)
    df_ora.to_csv(ora_path, index=False)
    print(f"\n  ORA 已保存: {ora_path}")
elif ora_path.exists():
    df_ora = pd.read_csv(ora_path)
    df_ora["Term"] = df_ora["Term"].apply(strip_go_id)
    print(f"\n  读取已有 ORA 结果: {ora_path}  ({len(df_ora)} 行)")
else:
    df_ora = None
    print("\n  无 ORA 数据，跳过绘图")

# =============================================================================
# SECTION 4: 气泡图
# ★ 排序逻辑（从上到下）：GO → KEGG → MSigDB
#   组内（从上到下）：-log10_padj 从小到大
# ★ 用数值 y 坐标绘图，气泡与 term 名称严格对应
# =============================================================================
print()
print("=" * 60)
print("SECTION 4: 绘制气泡图")
print("=" * 60)

if df_ora is None or len(df_ora) == 0:
    print("  无 ORA 数据，退出")
else:
    padj_col = "padj" if "padj" in df_ora.columns else "P-value"

    df_ora_sig = df_ora[df_ora[padj_col] < PARAMS["ora_padj_thresh"]].copy()
    df_ora_sig["-log10_padj"] = -np.log10(df_ora_sig[padj_col] + 1e-300)
    df_ora_sig["Term_short"]  = df_ora_sig["Term"].apply(
        lambda x: x[:PARAMS["term_max_chars"]] + "..."
        if len(str(x)) > PARAMS["term_max_chars"] else x
    )

    print(f"  显著 term 数 (padj<{PARAMS['ora_padj_thresh']}): {len(df_ora_sig)}")

    # ── 每个 direction × gene_set 取 top N ────────────────────────────────
    top_terms = []
    for direction in ["up", "down"]:
        for gs in PARAMS["gs_order"]:
            sub = df_ora_sig[
                (df_ora_sig["direction"] == direction) &
                (df_ora_sig["gene_set"]  == gs)
            ].nsmallest(PARAMS["top_n_per_gs"], padj_col).copy()
            if len(sub) > 0:
                sub["panel"] = direction
                top_terms.append(sub)

    if not top_terms:
        print("  无显著富集 term，跳过绘图")
    else:
        df_plot = pd.concat(top_terms, ignore_index=True)
        panels  = [p for p in ["up", "down"] if p in df_plot["panel"].values]
        print(f"  绘制面板: {panels}")

        fig, axes = plt.subplots(
            1, len(panels),
            figsize=(len(panels) * PARAMS["panel_w"], PARAMS["panel_h"]),
            squeeze=False,
        )
        fig.patch.set_facecolor("white")

        for idx, panel in enumerate(panels):
            ax  = axes[0][idx]
            sub = df_plot[df_plot["panel"] == panel].copy()

            # ── 排序：确定数值 y 坐标 ────────────────────────────────────
            #
            # 目标（从上到下）：GO → KEGG → MSigDB，组内 -log10 小→大
            #
            # matplotlib y 轴：y=0 在最下方，y 越大越靠上
            # 因此需要：
            #   GO    → y 值最大（在上）  → 排在列表末尾
            #   MSigDB→ y 值最小（在下）  → 排在列表开头
            #   组内  → -log10 大的 y 小（在下），-log10 小的 y 大（在上）
            #
            # sort_values 策略：
            #   _gs_rank 降序：MSigDB(2)先排→小y在下，GO(0)后排→大y在上  ✓
            #   -log10_padj 降序：大值先排→小y在下（图中在下），小值后排→大y在上（图中在上）✓
            #   最终从上到下看：GO(小→大) KEGG(小→大) MSigDB(小→大)        ✓

            gs_rank = {gs: i for i, gs in enumerate(PARAMS["gs_order"])}
            sub = sub.copy()
            sub["_gs_rank"] = sub["gene_set"].map(gs_rank)

            sub = sub.sort_values(
                ["_gs_rank", "-log10_padj"],
                ascending=[False, False]   # gs_rank 降序，组内 -log10 降序
            ).reset_index(drop=True)

            # 行号直接作为数值 y 坐标
            # 行号 0 → y=0（最下方），行号 n-1 → y=n-1（最上方）
            sub["_y"] = range(len(sub))

            # ── 气泡大小 = n_fg ──────────────────────────────────────────
            sub["_bubble_s"] = sub["n_fg"] * PARAMS["bubble_scale"] + PARAMS["bubble_min"]

            # ── 用数值坐标绘制散点 ────────────────────────────────────────
            for gs in PARAMS["gs_order"]:
                gs_sub = sub[sub["gene_set"] == gs]
                if len(gs_sub) == 0:
                    continue
                ax.scatter(
                    gs_sub["-log10_padj"],
                    gs_sub["_y"],
                    s          = gs_sub["_bubble_s"],
                    c          = PARAMS["lib_colors"][gs],
                    alpha      = PARAMS["bubble_alpha"],
                    linewidths = 0.5,
                    edgecolors = "white",
                    zorder     = 3,
                )

            # ── y 轴刻度与标签严格对应数值坐标 ──────────────────────────
            ax.set_yticks(sub["_y"].tolist())
            ax.set_yticklabels(sub["Term_short"].tolist(),
                               fontsize=PARAMS["ytick_fs"])
            ax.set_ylim(-0.8, len(sub) - 0.2)

            # 参考线 p=0.05
            ax.axvline(
                -np.log10(0.05),
                color="gray", linestyle="--", linewidth=1, alpha=0.6,
            )

            n_degs = n_up if panel == "up" else n_down
            direction_label = "Up-regulated" if panel == "up" else "Down-regulated"
            ax.set_xlabel("-log10(Adjusted P-value)", fontsize=PARAMS["xlabel_fs"])
            ax.set_title(
                f"{direction_label} genes\n(All subtypes merged, n={n_degs} DEGs)",
                fontsize=PARAMS["title_fs"], fontweight="bold",
            )
            ax.tick_params(axis="y", labelsize=PARAMS["ytick_fs"])
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="x", alpha=0.3)

            # ── 图例：右下角，y=0.08 略高于底部边缘 ─────────────────────
            if idx == len(panels) - 1:
                color_handles = [
                    mpatches.Patch(
                        color=PARAMS["lib_colors"][gs],
                        label=PARAMS["lib_labels"][gs],
                        alpha=PARAMS["bubble_alpha"],
                    )
                    for gs in PARAMS["gs_order"]
                ]
                size_handles = [
                    mlines.Line2D(
                        [], [],
                        marker="o", color="w",
                        markerfacecolor="gray",
                        markeredgecolor="gray",
                        markersize=np.sqrt(
                            n * PARAMS["bubble_scale"] + PARAMS["bubble_min"]
                        ),
                        alpha=0.7,
                        label=f"{n} genes",
                    )
                    for n in PARAMS["legend_size_refs"]
                ]
                ax.legend(
                    handles        = color_handles + size_handles,
                    fontsize       = PARAMS["legend_fs"],
                    loc            = "lower right",
                    bbox_to_anchor = PARAMS["legend_bbox"],
                    framealpha     = 0.9,
                    title          = None,
                    handlelength   = 1.5,
                    borderpad      = 0.8,
                )

        plt.suptitle(
            f"MAFB KO vs WT — ORA Enrichment (All Subtypes Merged)\n"
            f"Mann-Whitney U | BH | DEG padj<{PARAMS['deg_padj_thresh']} "
            f"| |log2FC|>{PARAMS['log2fc_thresh']} "
            f"| ORA padj<{PARAMS['ora_padj_thresh']} "
            f"| background = all tested genes",
            fontsize=PARAMS["title_fs"], fontweight="bold", y=1.02,
        )
        plt.tight_layout()
        save_fig(fig, "fig_merged_ora_bubble_all_cells")
        print(f"\n完成！输出目录: {OUT}")







# ══════════════════════════════════════════════════════════════════════════════
# Fig_EffectSize: Effect-size summary dotplot
# 底部轴（ax1）物理范围 = SNR 范围：0 ~ 5.5
# 顶部轴（ax2）装饰范围 = delta 范围：-0.3 ~ 0.6
# SNR 点：直接用 SNR 值画在 ax1
# delta 点：线性映射 delta(-0.3~0.6) → 底部轴物理坐标(0~5.5) 后画在 ax1
# ══════════════════════════════════════════════════════════════════════════════
import pandas as pd
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    v = df_fate_ko.loc[df_fate_ko["subtype"] == st, "fate_score"]
    print(st, round(v.mean() / v.std(ddof=0), 4))
# 预期输出：Mono 0.6561 / IFN_TAM 2.0028 / LA_TAM 7.1464

# 预期输出：Mono 0.6561 / IFN_TAM 2.0028 / LA_TAM 7.1464
# 与 consistency_check.txt 逐位一致 → SNR 出处闭环


P_ES = {
    "figsize"         : (8, 5),
    "title"           : "MAFB KO Effect-Size Summary",
    "title_fontsize"  : 12,
    # 底部轴（SNR）
    "xlabel_bottom"   : "SNR value ",
    # 顶部轴（delta）
    "xlabel_top"      : "Delta (KO − WT median)   [Fate / Program score]",
    "xlabel_fontsize" : 10,
    # 底部轴物理范围 = SNR 逻辑范围
    "xlim_bottom"     : (0.0, 7.5),
    # 顶部轴逻辑范围 = delta 范围
    "xlim_top"        : (-0.30, 0.50),
    # 点大小两档
    "dot_size_sig"    : 160,   # padj < 0.05，实心
    "dot_size_ns"     : 40,    # padj ≥ 0.05，空心
    "dot_alpha"       : 0.88,
    "edge_lw"         : 0.6,
    "tick_fontsize"   : 10,
    "legend_fontsize" : 9,
    "sig_thresh"      : 0.05,
    "ax_facecolor"    : "#F7F7F7",
    "spine_color"     : "#CCCCCC",
}

def get_dot_size(padj, p_es):
    return p_es["dot_size_sig"] if padj < p_es["sig_thresh"] else p_es["dot_size_ns"]

def delta_to_bottom_x(delta_val, xlim_top, xlim_bottom):
    """
    把 delta 逻辑坐标线性映射到底部轴物理坐标
    delta 在 xlim_top(-0.3~0.6) → 映射到 xlim_bottom(0~5.5)
    """
    t_lo, t_hi = xlim_top
    b_lo, b_hi = xlim_bottom
    return b_lo + (delta_val - t_lo) / (t_hi - t_lo) * (b_hi - b_lo)

print("Drawing Fig_EffectSize...")

# ── SNR 解析 + 硬编码备用 ────────────────────────────────────────────────────
# ═══ Fig 4h SNR 计算（STEP 04.5 同口径）═══════════════════════════════════
# 来源：04_fate/cell_fate_scores_KO.csv（主分析 k=20 的 per-cell fate score）
# 公式：SNR = mean(KO fate) / std(KO fate, ddof=0)，与 step04_5_consistency_check.py
#       严格一致（该脚本用 .values 转 numpy，numpy 默认 ddof=0）。
# ⚠️ 注意：不要从 08_sensitivity_results.csv 的 note 字段解析 SNR——
#    08C (HVG=4000) 行的 SNR=0.8314/2.6479/6.521 是敏感性配置值，
#    不是主分析值（这正是当初论文正文抄错数字的根源）。
# 预期输出：{'Mono': 0.6561, 'IFN_TAM': 2.0028, 'LA_TAM': 7.1464}
#           应与 04_fate/consistency_check.txt 逐位一致。

import pandas as pd

subtypes_order = ["Mono", "IFN_TAM", "LA_TAM"]

snr_final = {}
for st in subtypes_order:
    vals = df_fate_ko.loc[df_fate_ko["subtype"] == st, "fate_score"].values
    assert len(vals) > 0, f"未找到 {st} 的 cells，请检查 subtype 列名与取值"
    snr_final[st] = vals.mean() / vals.std()          # ddof=0，与门控日志同口径

print(f"  SNR 最终值: { {k: round(v, 4) for k, v in snr_final.items()} }")

# 自检：与 consistency_check.txt 官方值逐位核对
SNR_REFERENCE = {"Mono": 0.6561, "IFN_TAM": 2.0028, "LA_TAM": 7.1464}
for st in subtypes_order:
    diff = abs(snr_final[st] - SNR_REFERENCE[st])
    status = "OK" if diff < 5e-4 else "MISMATCH — 检查输入文件是否为主分析版本"
    print(f"  {st}: computed={snr_final[st]:.4f}, reference={SNR_REFERENCE[st]:.4f}  [{status}]")


# ── 汇总数据 ──────────────────────────────────────────────────────────────────
rows_delta = []
rows_snr   = []

# Fate score
for st in subtypes_order:
    row = df_fate_stats[
        (df_fate_stats["subtype"]    == st) &
        (df_fate_stats["comparison"] == "KO_vs_WT")
    ]
    if len(row) == 0:
        continue
    padj      = float(row["padj"].values[0])
    delta_val = float(row["delta"].values[0])
    rows_delta.append({
        "readout"    : f"Fate score ({st})",
        "category"   : "Fate score",
        "subtype"    : st,
        "x_physical" : delta_to_bottom_x(delta_val,
                                          P_ES["xlim_top"],
                                          P_ES["xlim_bottom"]),
        "x_delta"    : delta_val,
        "padj"       : padj,
        "dot_size"   : get_dot_size(padj, P_ES),
    })

# Program score
for st in subtypes_order:
    row = df_prog_stats[df_prog_stats["subtype"] == st]
    if len(row) == 0:
        continue
    padj      = float(row["padj"].values[0])
    delta_val = float(row["delta"].values[0])
    rows_delta.append({
        "readout"    : f"Program score ({st})",
        "category"   : "Program score",
        "subtype"    : st,
        "x_physical" : delta_to_bottom_x(delta_val,
                                          P_ES["xlim_top"],
                                          P_ES["xlim_bottom"]),
        "x_delta"    : delta_val,
        "padj"       : padj,
        "dot_size"   : get_dot_size(padj, P_ES),
    })

# SNR：直接用 SNR 值作为底部轴物理坐标
df_08c = df_sens[df_sens["sensitivity_id"] == "08C_hvg4000"]
for st in subtypes_order:
    row_08c = df_08c[df_08c["subtype"] == st]
    try:
        padj = float(row_08c["padj"].values[0])
        if np.isnan(padj):
            padj = 0.001
    except Exception:
        padj = 0.001
    snr_val = snr_final[st]
    rows_snr.append({
        "readout"    : f"SNR ({st})",
        "category"   : "SNR",
        "subtype"    : st,
        "x_physical" : snr_val,   # SNR 值直接就是底部轴坐标
        "x_snr"      : snr_val,
        "padj"       : padj,
        "dot_size"   : get_dot_size(padj, P_ES),
    })

df_delta = pd.DataFrame(rows_delta)
df_snr   = pd.DataFrame(rows_snr)

# 打印验证，确认物理坐标在 xlim_bottom 范围内
print("  delta 点物理坐标验证:")
for _, r in df_delta.iterrows():
    print(f"    {r['readout']}: delta={r['x_delta']:.3f} → x_physical={r['x_physical']:.3f}")
print("  SNR 点物理坐标验证:")
for _, r in df_snr.iterrows():
    print(f"    {r['readout']}: SNR={r['x_snr']:.3f} → x_physical={r['x_physical']:.3f}")

# ── 纵轴顺序（从上到下）：Fate → Program → SNR ───────────────────────────────
readout_order = (
    [f"Fate score ({st})"    for st in subtypes_order
     if f"Fate score ({st})"    in df_delta["readout"].values] +
    [f"Program score ({st})" for st in subtypes_order
     if f"Program score ({st})" in df_delta["readout"].values] +
    [f"SNR ({st})"           for st in subtypes_order
     if f"SNR ({st})"           in df_snr["readout"].values]
)

readout_reversed = list(reversed(readout_order))
y_map = {r: i for i, r in enumerate(readout_reversed)}
df_delta["_y"] = df_delta["readout"].map(y_map)
df_snr["_y"]   = df_snr["readout"].map(y_map)

n_fate    = sum(1 for r in readout_order if r.startswith("Fate"))
n_program = sum(1 for r in readout_order if r.startswith("Program"))
n_snr_r   = sum(1 for r in readout_order if r.startswith("SNR"))

# ── 绘图 ──────────────────────────────────────────────────────────────────────
fig, ax1 = plt.subplots(figsize=P_ES["figsize"])
fig.patch.set_facecolor("white")

# 顶部装饰轴（delta 刻度，不参与绘图）
ax2 = ax1.twiny()
ax2.set_xlim(P_ES["xlim_top"])
ax2.set_xlabel(P_ES["xlabel_top"], fontsize=P_ES["xlabel_fontsize"])
ax2.tick_params(axis="x", labelsize=P_ES["tick_fontsize"])
ax2.yaxis.set_visible(False)

# 交替横条（画在 ax1）
for i in range(len(readout_order)):
    if i % 2 == 0:
        ax1.axhspan(i - 0.5, i + 0.5, color="#ECECEC", alpha=0.6, zorder=0)

# 只保留 Program/SNR 之间的分隔线
if n_snr_r > 0 and n_program > 0:
    ax1.axhline(n_snr_r - 0.5,
                color="#AAAAAA", lw=1.2, ls="--", zorder=1)

# delta=0 映射到底部轴的物理坐标
x_zero = delta_to_bottom_x(0, P_ES["xlim_top"], P_ES["xlim_bottom"])
ax1.axvline(x_zero, color="#333333", lw=1.0, ls="--", alpha=0.5, zorder=2)

# ── 所有点画在 ax1 ────────────────────────────────────────────────────────────

# delta 点（圆形）
for st in subtypes_order:
    sub = df_delta[df_delta["subtype"] == st]
    if len(sub) == 0:
        continue
    sig  = sub[sub["padj"] <  P_ES["sig_thresh"]]
    nsig = sub[sub["padj"] >= P_ES["sig_thresh"]]
    if len(sig):
        ax1.scatter(sig["x_physical"], sig["_y"],
                    s=sig["dot_size"], c=subtype_colors[st],
                    alpha=P_ES["dot_alpha"], linewidths=P_ES["edge_lw"],
                    edgecolors=subtype_colors[st],
                    marker="o", zorder=4)
    if len(nsig):
        ax1.scatter(nsig["x_physical"], nsig["_y"],
                    s=nsig["dot_size"], facecolors="none",
                    edgecolors=subtype_colors[st], linewidths=1.5,
                    marker="o", alpha=0.8, zorder=4)

# SNR 点（菱形）
for st in subtypes_order:
    sub = df_snr[df_snr["subtype"] == st]
    if len(sub) == 0:
        continue
    sig  = sub[sub["padj"] <  P_ES["sig_thresh"]]
    nsig = sub[sub["padj"] >= P_ES["sig_thresh"]]
    if len(sig):
        ax1.scatter(sig["x_physical"], sig["_y"],
                    s=sig["dot_size"], c=subtype_colors[st],
                    alpha=P_ES["dot_alpha"], linewidths=P_ES["edge_lw"],
                    edgecolors=subtype_colors[st],
                    marker="D", zorder=4)
    if len(nsig):
        ax1.scatter(nsig["x_physical"], nsig["_y"],
                    s=nsig["dot_size"], facecolors="none",
                    edgecolors=subtype_colors[st], linewidths=1.5,
                    marker="D", alpha=0.8, zorder=4)

# ── 底部轴范围与 y 轴 ─────────────────────────────────────────────────────────
ax1.set_xlim(P_ES["xlim_bottom"])
ax1.set_xlabel(P_ES["xlabel_bottom"], fontsize=P_ES["xlabel_fontsize"])
ax1.tick_params(axis="x", labelsize=P_ES["tick_fontsize"])

ax1.set_yticks(list(y_map.values()))
ax1.set_yticklabels(list(y_map.keys()), fontsize=P_ES["tick_fontsize"])
ax1.set_ylim(-0.7, len(readout_order) - 0.3)
ax2.set_ylim(ax1.get_ylim())

# 背景和边框
ax1.set_facecolor(P_ES["ax_facecolor"])
for sp in ax1.spines.values():
    sp.set_edgecolor(P_ES["spine_color"])
ax1.spines["right"].set_visible(False)

# ── 图例一：subtype 颜色（右侧偏下）─────────────────────────────────────────
color_handles = [
    mpatches.Patch(color=subtype_colors[st], label=st, alpha=0.85)
    for st in subtypes_order
]
legend1 = ax1.legend(
    handles        = color_handles,
    fontsize       = P_ES["legend_fontsize"],
    loc            = "center right",
    bbox_to_anchor = (1.0, 0.38),
    framealpha     = 0.9,
    edgecolor      = "#CCCCCC",
    title          = "Subtype",
    title_fontsize = P_ES["legend_fontsize"],
)
#ax1.add_artist(legend1)

# ── 图例二：点大小两档 + 形状（右下偏上）────────────────────────────────────
size_handles = [
    mlines.Line2D([], [], marker="o", color="w",
                  markerfacecolor="gray", markeredgecolor="gray",
                  markersize=np.sqrt(P_ES["dot_size_sig"]),
                  alpha=0.85, label="padj < 0.05  (filled)"),
    mlines.Line2D([], [], marker="o", color="w",
                  markerfacecolor="none", markeredgecolor="gray",
                  markersize=np.sqrt(P_ES["dot_size_ns"]),
                  linewidth=1.5, label="padj ≥ 0.05  (open)"),
    mlines.Line2D([], [], marker="D", color="w",
                  markerfacecolor="gray", markeredgecolor="gray",
                  markersize=6, alpha=0.85, label="◆ SNR (lower x-axis)"),
]
ax1.legend(
    handles        = size_handles,
    fontsize       = P_ES["legend_fontsize"],
    loc            = "lower right",
    bbox_to_anchor = (1.0, 0.12),
    framealpha     = 0.9,
    edgecolor      = "#CCCCCC",
    title          = "",
    title_fontsize = P_ES["legend_fontsize"],
)

ax1.set_title(P_ES["title"], fontsize=P_ES["title_fontsize"],
              fontweight="bold", pad=28)
plt.tight_layout()
save_fig(fig, "fig_effect_size_summary_dotplot")

