# -*- coding: utf-8 -*-
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))  # ★ 仓库根目录动态推导
from config.paths import translate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.patheffects as pe
import os
from config.paths import translate as _t
VOLCANO_CSV = _t(r'D:/bulk-download/GSE182434/cellcomm/volcano_data_malignant_vs_normal_monomac.csv')
OUTPUT_DIR = str(_Path(__file__).resolve().parents[1] / 'results' / 'uav_panels' / 'fig1f_v2')  # ★ 动态推导
DPI = 300

# Fig6 前置参数：Volcano Plot — Malignant vs Normal B → Monocytes/Macrophages
# ══════════════════════════════════════════════════════════════════════════════

P06 = {

    # ── 数据路径 ──────────────────────────────────────────────────────────────
    "volcano_csv"             : VOLCANO_CSV,

    # ── Y 轴截断与抖动 ────────────────────────────────────────────────────────
    "y_jitter_cap"            : 9.9,     # neg_log10_pval >= 此值则截断
    "y_jitter_base"           : 10.0,    # 截断后基准位置
    "y_jitter_range"          : 2.0,     # 抖动范围 uniform(0, range)
    "random_seed"             : 42,

    # ── 颜色映射（4 类别）────────────────────────────────────────────────────
    "palette" : {
        'Malignant-enriched':    '#d62728',
        'Normal-enriched':       '#2ca02c',
        'Significant (small Δ)': '#ff7f0e',
        'Not significant':       '#e0e0e0',
    },

    # ── 点大小（面积单位）────────────────────────────────────────────────────
    "size_map" : {
        'Malignant-enriched':    70,
        'Normal-enriched':       70,
        'Significant (small Δ)': 50,
        'Not significant':       30,
    },

    # ── 点透明度 ──────────────────────────────────────────────────────────────
    "alpha_map" : {
        'Malignant-enriched':    0.85,
        'Normal-enriched':       0.85,
        'Significant (small Δ)': 0.75,
        'Not significant':       0.40,
    },

    # ── 点边框 ────────────────────────────────────────────────────────────────
    "edge_color_ns"           : 'white',
    "edge_color_sig"          : '#333333',
    "edge_lw_ns"              : 0.3,
    "edge_lw_sig"             : 0.5,

    # ── 图像尺寸 ──────────────────────────────────────────────────────────────
    "figsize"                 : (9.8, 6.96),

    # ── 阈值参考线 ────────────────────────────────────────────────────────────
    "pval_threshold"          : 0.05,    # 横线位置 -log10(0.05)
    "delta_threshold"         : 0.1,     # 竖线位置 ±0.1
    "hline_color"             : '#555555',
    "vline_color"             : '#555555',
    "hline_lw"                : 1.0,
    "vline_lw"                : 1.0,
    "hline_alpha"             : 0.7,
    "vline_alpha"             : 0.6,
    "zero_line_color"         : '#aaaaaa',
    "zero_line_lw"            : 0.8,
    "zero_line_alpha"         : 0.5,

    # ── 标注点选择 ────────────────────────────────────────────────────────────
    "n_label_mal_top"         : 4,       # Malignant-enriched，Delta 最大前 N 个
    "n_label_nor_top"         : 4,       # Normal-enriched，Delta 最小前 N 个
    "n_label_mid_sig"         : 3,       # 中部显著点（p 值最大，排除截断点）
    "label_offset_x"          : 0.0,
    "label_offset_y"          : 0.0,
    "label_fontsize"          : 9,
    "label_stroke_lw"         : 2.0,
    "label_stroke_color"      : 'white',
    "label_arrow_color"       : 'gray',
    "label_arrow_lw"          : 0.8,

    # ── 统计信息文字 ──────────────────────────────────────────────────────────
    "stats_text_x"            : 0.63,
    "stats_text_y"            : 0.02,
    "stats_text_fontsize"     : 12,
    "stats_text_color"        : '#666666',

    # ── 轴范围 ────────────────────────────────────────────────────────────────
    "ymax"                    : 13,
    "ymin"                    : -0.3,
    "xmargin"                 : 0.08,

    # ── 方向提示文字 ──────────────────────────────────────────────────────────
    "dir_text_y_ratio"        : 0.7,
    "dir_text_fontsize"       : 15,
    "mal_dir_text"            : 'Malignant-enriched →',
    "nor_dir_text"            : '← Normal-enriched',
    "pval_label_text"         : 'p = 0.05',
    "pval_label_fontsize"     : 12,

    # ── 轴标签 ────────────────────────────────────────────────────────────────
    "xlabel"                  : 'ΔMagnitude  (Malignant B − Normal B interaction strength)',
    "xlabel_fontsize"         : 15,
    "xlabel_labelpad"         : 8,
    "ylabel"                  : '−log10(p-value)  [CellPhoneDB permutation test]',
    "ylabel_fontsize"         : 15,
    "ylabel_labelpad"         : 8,

    # ── 图标题 ────────────────────────────────────────────────────────────────
    "title"                   : (''),
    "title_fontsize"          : 12,
    "title_pad"               : 14,

    # ── 图例 ──────────────────────────────────────────────────────────────────
    "legend_loc"              : 'lower right',  # 虚线上方 5mm、右对齐
    "legend_fontsize"         : 10,
    "legend_title"            : 'Category',
    "legend_title_fontsize"   : 13,
    "legend_bbox"             : (1.0, 0.156),  # p=0.05 线 yfrac=0.12 + 5mm
    "legend_framealpha"       : 0.9,

    # ── 高显著性区域背景 ──────────────────────────────────────────────────────
    "highlight_ymin"          : 9.9,
    "highlight_alpha"         : 0.04,
    "highlight_color"         : 'gray',
    "capped_label_text"       : 'p ≤ 1e−10\n(capped, jittered)',
    "capped_label_fontsize"   : 12,
    "capped_label_color"      : '#888888',
    "capped_label_style"      : 'italic',
    "capped_label_y"          : 9,      # get_yaxis_transform 坐标下的 y 位置

    # ── Seaborn 主题 ──────────────────────────────────────────────────────────
    "sns_style"               : 'ticks',
    "sns_font_scale"          : 1.0,

    # ── 输出文件名 ────────────────────────────────────────────────────────────
    "output_name"             : 'volcano_malignant_vs_normal_monomac',
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig6 绘图
# ══════════════════════════════════════════════════════════════════════════════

print("\n[Fig 6] Drawing volcano plot (Malignant vs Normal → Mono/Mac)...")

# ── 读取数据 ──────────────────────────────────────────────────────────────────
df_vol = pd.read_csv(P06["volcano_csv"])
print(f"  Volcano data loaded: {df_vol.shape}")

# ── Y 轴截断与抖动 ────────────────────────────────────────────────────────────
np.random.seed(P06["random_seed"])
cap_mask = df_vol['neg_log10_pval'] >= P06["y_jitter_cap"]
df_vol['y_plot'] = df_vol['neg_log10_pval'].copy()
df_vol.loc[cap_mask, 'y_plot'] = (
    P06["y_jitter_base"]
    + np.random.uniform(0, P06["y_jitter_range"], cap_mask.sum())
)

# ── Seaborn 主题 ──────────────────────────────────────────────────────────────
sns.set_theme(style=P06["sns_style"], font_scale=P06["sns_font_scale"])

fig6, ax6 = plt.subplots(figsize=P06["figsize"])

# ── 按类别分层绘制 ────────────────────────────────────────────────────────────
for cat in ['Not significant', 'Significant (small Δ)', 'Normal-enriched', 'Malignant-enriched']:
    sub_vol    = df_vol[df_vol['category'] == cat]
    edge_color = P06["edge_color_ns"] if cat == 'Not significant' else P06["edge_color_sig"]
    edge_lw    = P06["edge_lw_ns"]    if cat == 'Not significant' else P06["edge_lw_sig"]
    ax6.scatter(
        sub_vol['Delta_Magnitude'], sub_vol['y_plot'],
        s          = P06["size_map"][cat],
        c          = P06["palette"][cat],
        alpha      = P06["alpha_map"][cat],
        edgecolors = edge_color,
        linewidths = edge_lw,
        label      = f"{cat} (n={len(sub_vol)})",
        zorder     = 3,
    )

# ── 阈值参考线 ────────────────────────────────────────────────────────────────
pval_line = -np.log10(P06["pval_threshold"])
ax6.axhline(pval_line,
            color=P06["hline_color"], linewidth=P06["hline_lw"],
            linestyle='--', alpha=P06["hline_alpha"], zorder=2)
ax6.axvline( P06["delta_threshold"],
            color=P06["vline_color"], linewidth=P06["vline_lw"],
            linestyle='--', alpha=P06["vline_alpha"], zorder=2)
ax6.axvline(-P06["delta_threshold"],
            color=P06["vline_color"], linewidth=P06["vline_lw"],
            linestyle='--', alpha=P06["vline_alpha"], zorder=2)
ax6.axvline(0,
            color=P06["zero_line_color"], linewidth=P06["zero_line_lw"],
            linestyle='-', alpha=P06["zero_line_alpha"], zorder=1)

# ── 标注关键 LR pair ──────────────────────────────────────────────────────────
top_mal_vol = df_vol[df_vol['category'] == 'Malignant-enriched'].nlargest(
    P06["n_label_mal_top"], 'Delta_Magnitude')
top_nor_vol = df_vol[df_vol['category'] == 'Normal-enriched'].nsmallest(
    P06["n_label_nor_top"], 'Delta_Magnitude')
mid_sig_vol = (df_vol[(df_vol['category'] == 'Malignant-enriched') &
                       (df_vol['neg_log10_pval'] < P06["y_jitter_cap"])]
               .nlargest(P06["n_label_mid_sig"], 'neg_log10_pval'))

to_label_vol = pd.concat([top_mal_vol, top_nor_vol, mid_sig_vol]).drop_duplicates('LR_pair')

texts6 = []
for _, row in to_label_vol.iterrows():
    color    = P06["palette"].get(row['category'], '#333333')
    offset_x = P06["label_offset_x"] if row['Delta_Magnitude'] > 0 else -P06["label_offset_x"]
    t = ax6.annotate(
        row['LR_pair'],
        xy     = (row['Delta_Magnitude'], row['y_plot']),
        xytext = (row['Delta_Magnitude'] + offset_x,
                  row['y_plot'] + P06["label_offset_y"]),
        fontsize     = P06["label_fontsize"],
        color        = color,
        fontweight   = 'bold',
        path_effects = [pe.withStroke(linewidth=P06["label_stroke_lw"],
                                      foreground=P06["label_stroke_color"])],
        arrowprops   = dict(arrowstyle='->', color=P06["label_arrow_color"],
                            lw=P06["label_arrow_lw"], shrinkA=3, shrinkB=3),
        zorder = 5,
    )
    texts6.append(t)

# 如需启用 adjustText 防重叠，取消下面注释：
from adjustText import adjust_text
adjust_text(
    texts6,
    x=df_vol['Delta_Magnitude'].values,
    y=df_vol['y_plot'].values,
    arrowprops=dict(arrowstyle='->', color='gray', lw=0.5),
    expand_points=(1.5, 1.5),   # 标注与数据点的排斥范围
    expand_text=(1.3, 1.3),     # 标注之间的排斥范围
    force_points=0.5,           # 对数据点的推力
    force_text=0.8,             # 对文字的推力
    max_move=3.0,               # 最大移动距离（数据坐标单位）
    ax=ax6,
)


# ── 统计信息文字 ──────────────────────────────────────────────────────────────
stats_text = (
    f"Total: {len(df_vol)}\n"
    f"Malignant-enriched: {len(df_vol[df_vol['category']=='Malignant-enriched'])}\n"
    f"Normal-enriched: {len(df_vol[df_vol['category']=='Normal-enriched'])}"
)
ax6.text(P06["stats_text_x"], P06["stats_text_y"], stats_text,
         transform=ax6.transAxes, fontsize=P06["stats_text_fontsize"],
         ha='right', va='bottom', color=P06["stats_text_color"])

# ── 轴范围 ────────────────────────────────────────────────────────────────────
ax6.set_ylim(P06["ymin"], P06["ymax"])
ax6.set_xlim(
    df_vol['Delta_Magnitude'].min() - P06["xmargin"],
    df_vol['Delta_Magnitude'].max() + P06["xmargin"],
)

# ── 方向提示文字 ──────────────────────────────────────────────────────────────
dir_y = P06["ymax"] * P06["dir_text_y_ratio"]
ax6.text( 0.2, dir_y, P06["mal_dir_text"],
          fontsize=P06["dir_text_fontsize"],
          color=P06["palette"]['Malignant-enriched'],
          fontweight='bold', ha='left', va='top')
ax6.text(-0.2, dir_y, P06["nor_dir_text"],
          fontsize=P06["dir_text_fontsize"],
          color=P06["palette"]['Normal-enriched'],
          fontweight='bold', ha='right', va='top')
ax6.text(df_vol['Delta_Magnitude'].max() + P06["xmargin"],
         pval_line + 0.1,
         P06["pval_label_text"],
         fontsize=P06["pval_label_fontsize"],
         color=P06["hline_color"], va='bottom', ha='right')

# ── 轴标签 & 标题 ─────────────────────────────────────────────────────────────
ax6.set_xlabel(P06["xlabel"],
               fontsize=P06["xlabel_fontsize"], labelpad=P06["xlabel_labelpad"])
ax6.set_ylabel(P06["ylabel"],
               fontsize=P06["ylabel_fontsize"], labelpad=P06["ylabel_labelpad"])
ax6.set_title(P06["title"],
              fontsize=P06["title_fontsize"], fontweight='bold', pad=P06["title_pad"])

# ── 图例 ──────────────────────────────────────────────────────────────────────
ax6.legend(loc=P06["legend_loc"], fontsize=P06["legend_fontsize"],
           frameon=True, framealpha=P06["legend_framealpha"],
           title=P06["legend_title"], title_fontsize=P06["legend_title_fontsize"],
           bbox_to_anchor=P06["legend_bbox"])

# ── 高显著性区域背景 ──────────────────────────────────────────────────────────
ax6.axhspan(P06["highlight_ymin"], P06["ymax"],
            alpha=P06["highlight_alpha"], color=P06["highlight_color"], zorder=0)
ax6.text(0.88, P06["capped_label_y"], P06["capped_label_text"],
         fontsize=P06["capped_label_fontsize"],
         color=P06["capped_label_color"],
         va='center', ha='left',
         style=P06["capped_label_style"],
         transform=ax6.get_yaxis_transform())

sns.despine(ax=ax6)
plt.tight_layout()

for fmt in ['png', 'svg']:
    fig6.savefig(
        os.path.join(OUTPUT_DIR, f'{P06["output_name"]}.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )

plt.close(fig6)
print("  ✓ volcano_malignant_vs_normal_monomac saved")


# ================================