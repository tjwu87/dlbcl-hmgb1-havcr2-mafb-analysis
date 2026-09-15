# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# =============================================================================
#  LIANA 细胞通讯可视化 — 本地复现脚本
#  数据根目录: D:/bulk-download/GSE182434/cellcomm/
#  所有图 DPI = 300
# =============================================================================

import matplotlib
matplotlib.use('Agg')   # 非交互式后端，适合本地保存图片；若想弹窗预览可改为 'TkAgg'

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch

import seaborn as sns
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
#  ① 路径配置 — 修改这里即可切换数据目录和输出目录
# =============================================================================
DATA_DIR   = translate(r'D:/bulk-download/GSE182434/cellcomm')   # CSV 数据所在目录
OUTPUT_DIR = translate(r'D:/bulk-download/GSE182434/cellcomm')   # 图片输出目录（可改为其他路径）
os.makedirs(OUTPUT_DIR, exist_ok=True)

VOLCANO_CSV = translate(r'D:/bulk-download/GSE182434/cellcomm/volcano_data_malignant_vs_normal_monomac.csv')
# 火山图数据单独放置，路径与 Agent 本地保存路径一致；如有不同请修改上面这行

DPI = 300   # ← 全局 DPI，所有图均使用此值；Agent 原始值为 150

# =============================================================================
#  ② 加载 CSV 数据
# =============================================================================
print("Loading data...")

# 完整 LIANA 结果（含所有细胞类型对）
liana_res = pd.read_csv(os.path.join(DATA_DIR, 'liana_results_full.csv'))

# 显著交互（magnitude_rank < 0.05）
sig = pd.read_csv(os.path.join(DATA_DIR, 'liana_results_significant.csv'))

# 交互权重矩阵（source × target，值为 Σ lr_means）
weight_matrix = pd.read_csv(os.path.join(DATA_DIR, 'interaction_weight_matrix.csv'), index_col=0)

print(f"liana_res shape      : {liana_res.shape}")
print(f"sig shape            : {sig.shape}")
print(f"weight_matrix shape  : {weight_matrix.shape}")

# =============================================================================
#  ③ 公共配色与短标签
# =============================================================================

# 每个细胞类型对应的十六进制颜色
# ← 修改这里的颜色值可以改变所有图中对应细胞类型的颜色
ct_palette = {
    'CD8 T cells':                '#1f77b4',
    'CD4 T cells':                '#aec7e8',
    'Tregs':                      '#17becf',
    'TFH cells':                  '#9edae5',
    'NK cells':                   '#2ca02c',
    'GC B cells':                 '#98df8a',
    'Proliferative GC B cells':   '#d62728',
    'Naive B cells':              '#ff9896',
    'Memory B cells':             '#9467bd',
    'Age-associated B cells':     '#c5b0d5',
    'B cells (other)':            '#8c564b',
    'Plasma cells':               '#ff7f0e',
    'Monocytes/Macrophages':      '#e377c2',
    'pDC/Other':                  '#7f7f7f',
}

# 图中显示的简短标签；修改右侧字符串可改变图中文字
ct_short = {
    'CD8 T cells': 'CD8 T', 'CD4 T cells': 'CD4 T',
    'Tregs': 'Tregs', 'TFH cells': 'TFH', 'NK cells': 'NK',
    'GC B cells': 'GC B', 'Proliferative GC B cells': 'Prolif. GC B',
    'Naive B cells': 'Naive B', 'Memory B cells': 'Memory B',
    'Age-associated B cells': 'Age-assoc. B', 'B cells (other)': 'B (other)',
    'Plasma cells': 'Plasma', 'Monocytes/Macrophages': 'Mono/Mac',
    'pDC/Other': 'pDC/Other',
}

# Malignant vs Normal 专用配色
ct_palette_mal = {
    'Malignant B cell':       '#d62728',
    'Normal B cell':          '#2ca02c',
    'CD8 T cells':            '#1f77b4',
    'CD4 T cells':            '#6baed6',
    'Tregs':                  '#17becf',
    'TFH cells':              '#31a354',
    'NK cells':               '#a1d99b',
    'Monocytes/Macrophages':  '#e377c2',
    'pDC/Other':              '#636363',
}
ct_short_mal = {
    'Malignant B cell': 'Malignant B', 'Normal B cell': 'Normal B',
    'CD8 T cells': 'CD8 T', 'CD4 T cells': 'CD4 T',
    'Tregs': 'Tregs', 'TFH cells': 'TFH', 'NK cells': 'NK',
    'Monocytes/Macrophages': 'Mono/Mac', 'pDC/Other': 'pDC/Other',
}

# =============================================================================
#  图1：Interaction Weight Heatmap（交互权重热图）
# =============================================================================
print("\n[Fig 1] Drawing interaction weight heatmap...")

fig1, ax1 = plt.subplots(
    figsize=(12, 10)   # ← (宽, 高) 英寸；增大可让标签更宽松
)

wm = weight_matrix.copy()
wm_log = np.log1p(wm)   # log1p 变换，压缩大值以便颜色分布更均匀

# 按行/列总信号强度降序排列，最活跃的细胞类型排在左上角
row_order = wm_log.sum(axis=1).sort_values(ascending=False).index
col_order = wm_log.sum(axis=0).sort_values(ascending=False).index
wm_plot = wm_log.loc[row_order, col_order]

# 替换为短标签
wm_plot.index   = [ct_short.get(i, i) for i in wm_plot.index]
wm_plot.columns = [ct_short.get(c, c) for c in wm_plot.columns]

sns.heatmap(
    wm_plot, ax=ax1,
    cmap='YlOrRd',          # ← 颜色方案；可换 'Blues', 'viridis', 'RdYlBu_r' 等
    linewidths=0.4,          # ← 格子间分割线宽度；0 表示无分割线
    linecolor='white',       # ← 分割线颜色
    cbar_kws={
        'label': 'Σ LR means (log1p)',   # ← colorbar 标签文字
        'shrink': 0.7,                   # ← colorbar 相对高度缩放比例
    },
    annot=True,              # ← True = 在格子内显示数值；False = 不显示
    fmt='.0f',               # ← 数值格式；'.1f' 保留1位小数，'.0f' 取整
    annot_kws={'size': 7.5}, # ← 格子内数字字体大小
)
ax1.set_xlabel('Target (Receiver)', fontsize=11, labelpad=8)   # ← X 轴标签及字号
ax1.set_ylabel('Source (Sender)',   fontsize=11, labelpad=8)   # ← Y 轴标签及字号
ax1.set_title(
    'Cell-Cell Communication Interaction Weight Matrix\n'
    'DLBCL — LIANA rank_aggregate (magnitude_rank < 0.05)',
    fontsize=12,             # ← 标题字号
    fontweight='bold',
)
ax1.tick_params(axis='x', rotation=45, labelsize=9)   # ← X 轴刻度旋转角度 & 字号
ax1.tick_params(axis='y', rotation=0,  labelsize=9)   # ← Y 轴刻度字号

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig1.savefig(
        os.path.join(OUTPUT_DIR, f'interaction_heatmap.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )
plt.close(fig1)
print("  ✓ interaction_heatmap saved")

# =============================================================================
#  图2：Circle Plot — 全细胞类型通讯网络
# =============================================================================
print("\n[Fig 2] Drawing circle communication plot (all cell types)...")

def draw_circle_comm(liana_df, palette, short_labels, title,
                     top_n=200,          # ← 取交互强度 Top N 条记录绘图；增大显示更多边
                     figsize=(13, 13),   # ← 图像尺寸（宽, 高）英寸
                     save_path=None):
    """
    自定义圆形通讯网络图。
    节点 = 细胞类型，边 = 通讯强度（越粗越强），节点大小 = 总通讯量。
    """
    df = liana_df.sort_values('magnitude_rank').head(top_n).copy()
    df['weight'] = 1 - df['magnitude_rank']   # 转化为"强度"，越大越强

    # 聚合：每对 source→target 的总权重
    pair_df = (df.groupby(['source', 'target'])
               .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
               .reset_index())

    active_cts = [c for c in palette
                  if c in set(pair_df['source']) | set(pair_df['target'])]
    n_ct = len(active_cts)

    # 每个节点的总通讯量（出 + 入）
    node_tot = {ct: (pair_df[pair_df['source']==ct]['total_weight'].sum() +
                     pair_df[pair_df['target']==ct]['total_weight'].sum())
                for ct in active_cts}

    # 节点均匀分布在圆上，从顶部顺时针排列
    angles = {ct: np.pi/2 - 2*np.pi*i/n_ct for i, ct in enumerate(active_cts)}
    pos    = {ct: (np.cos(angles[ct]), np.sin(angles[ct])) for ct in active_cts}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    w_vals = pair_df['total_weight'].values
    w_min, w_max = w_vals.min(), w_vals.max()

    # ── 绘制边（箭头）──────────────────────────────────────────────────────
    for _, row in pair_df.sort_values('total_weight').iterrows():
        src, tgt = row['source'], row['target']
        if src not in pos or tgt not in pos or src == tgt:
            continue
        w_norm = (row['total_weight'] - w_min) / (w_max - w_min + 1e-10)

        lw    = (0.4 + 5.6 * w_norm) * 3.0   # ← 边宽系数；3.0 = Agent 最终版本 ×3 加粗
        rad   = 0.2 + 0.15 * w_norm           # ← 弧度弯曲程度；越大弧越弯
        color = palette.get(src, '#888888')
        hw    = 0.015 + 0.025 * w_norm        # ← 箭头头部宽度
        hl    = 0.025 + 0.030 * w_norm        # ← 箭头头部长度

        ax.annotate(
            '', xy=pos[tgt], xytext=pos[src],
            arrowprops=dict(
                arrowstyle=f'-|>,head_width={hw:.3f},head_length={hl:.3f}',
                connectionstyle=f'arc3,rad={rad:.2f}',
                color=color,
                lw=lw,
                alpha=1.0,    # ← 边透明度；1.0 = 完全不透明
                shrinkA=8,   # ← 关键
                shrinkB=8    # ← 关键   
            ),
            zorder=2,
        )

    # ── 绘制节点（圆形）────────────────────────────────────────────────────
    ns_max = max(node_tot.values()) if node_tot else 1
    for ct in active_cts:
        x, y = pos[ct]
        ns_norm = node_tot.get(ct, 0) / ns_max
        r = (0.055 + 0.095 * ns_norm) * 0.75   # ← 0.75 = 节点半径缩放系数；1.0 为原始大小
        circle = plt.Circle((x, y), r,
                             color=palette.get(ct, '#888888'),
                             zorder=5,
                             linewidth=1.8,      # ← 节点白色描边宽度
                             edgecolor='white')
        ax.add_patch(circle)

    # ── 绘制标签 ────────────────────────────────────────────────────────────
    label_r = 1.20   # ← 标签距圆心的距离；增大可避免与节点重叠
    for ct in active_cts:
        ang = angles[ct]
        deg = np.degrees(ang) % 360
        ha  = 'right' if 80 < deg < 280 else 'left'   # ← 左半圆右对齐，右半圆左对齐
        lx  = (label_r + 0.08) * np.cos(ang)
        ly  = (label_r + 0.08) * np.sin(ang)

        ax.text(lx, ly, short_labels.get(ct, ct),
                ha=ha, va='center',
                fontsize=8.5,        # ← 标签字号
                fontweight='bold',
                color=palette.get(ct, '#333333'),
                path_effects=[pe.withStroke(linewidth=3.0, foreground='white')],
                zorder=6)

        nx, ny = pos[ct]
        ax.plot([nx * 1.07, label_r * np.cos(ang)],
                [ny * 1.07, label_r * np.sin(ang)],
                color=palette.get(ct, '#888888'),
                linewidth=0.6, alpha=0.5, zorder=3)   # ← 连接线宽度 & 透明度

    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.set_title(title,
                 fontsize=12,       # ← 标题字号
                 fontweight='bold',
                 pad=14)

    legend_els = [
        mpatches.Patch(color='#555555', label='Strong interaction (thick edge)'),
        mpatches.Patch(color='#aaaaaa', label='Weak interaction (thin edge)'),
        plt.scatter([], [], s=100, c='#555555', label='High comm. (large node)'),
        plt.scatter([], [], s=30,  c='#555555', label='Low comm. (small node)'),
    ]
    ax.legend(handles=legend_els,
              loc='lower left', fontsize=7.5,
              frameon=True, framealpha=0.9,
              bbox_to_anchor=(-0.08, -0.08),
              title='Edge / Node encoding', title_fontsize=8)

    if save_path:
        for fmt in ['png', 'svg']:
            fig.savefig(f'{save_path}.{fmt}',
                        dpi=DPI, bbox_inches='tight', format=fmt)
        plt.close(fig)
        print(f"  ✓ Saved: {save_path}")
    return fig


# 绘制全细胞类型圆形图
draw_circle_comm(
    liana_res, ct_palette, ct_short,
    title='Cell-Cell Communication Network\nDLBCL — LIANA rank_aggregate (Top 200 interactions)',
    top_n=200,          # ← 展示最强的前 200 条交互；增大可显示更多弱交互
    figsize=(13, 13),
    save_path=os.path.join(OUTPUT_DIR, 'circle_plot_v2')
)


# =============================================================================
#  图4：Dotplot — B 细胞亚型 → 免疫伙伴（Top LR Pairs）
# =============================================================================
print("\n[Fig 4] Drawing dotplot (B cells → immune partners)...")

# ── 筛选 source/target ──────────────────────────────────────────────────────
b_sources = ['GC B cells', 'Proliferative GC B cells', 'Memory B cells', 'Plasma cells']
# ← 修改上面列表可增减作为"信号发出方"的 B 细胞亚型

key_targets = ['CD4 T cells', 'CD8 T cells', 'Tregs', 'TFH cells', 'Monocytes/Macrophages']
# ← 修改上面列表可增减"信号接收方"细胞类型

sub = liana_res[
    liana_res['source'].isin(b_sources) &
    liana_res['target'].isin(key_targets)
].copy()

top_n_lr = 20   # ← 展示 Top N 个 LR pair；增大可显示更多行
top_lr = (sub.groupby(['ligand_complex', 'receptor_complex'])['magnitude_rank']
          .mean().sort_values().head(top_n_lr).index)

sub_top = sub[
    sub.apply(lambda r: (r['ligand_complex'], r['receptor_complex']) in top_lr, axis=1)
].copy()

src_short = {'GC B cells': 'GC B', 'Proliferative GC B cells': 'Prolif.GC B',
             'Memory B cells': 'Memory B', 'Plasma cells': 'Plasma'}
tgt_short = {'CD4 T cells': 'CD4 T', 'CD8 T cells': 'CD8 T',
             'Tregs': 'Tregs', 'TFH cells': 'TFH', 'Monocytes/Macrophages': 'Mono/Mac'}

sub_top['src_s'] = sub_top['source'].map(src_short)
sub_top['tgt_s'] = sub_top['target'].map(tgt_short)
sub_top['st']    = sub_top['src_s'] + '\n→ ' + sub_top['tgt_s']
sub_top['LR']    = sub_top['ligand_complex'] + ' → ' + sub_top['receptor_complex']

# ── 计算用于颜色和大小的分数 ────────────────────────────────────────────────
# spec_weight: 越高越细胞类型特异；log10 变换提升颜色对比度
sub_top['log_spec'] = np.log10(sub_top['spec_weight'] + 1e-6)
# magnitude_rank: 越低越强 → 转为 -log10 使"大点 = 强交互"
sub_top['neg_log_mag'] = -np.log10(sub_top['magnitude_rank'].clip(1e-10))

lr_order = [f"{l} → {r}" for l, r in top_lr]
st_order = sorted(sub_top['st'].unique())
lr_idx   = {lr: i for i, lr in enumerate(lr_order)}
st_idx   = {st: i for i, st in enumerate(st_order)}

# 颜色范围：去掉极端 5%/95% 分位数，避免颜色被极值拉偏
vmin_spec = sub_top['log_spec'].quantile(0.05)
vmax_spec = sub_top['log_spec'].quantile(0.95)
cmap_dot  = plt.cm.RdYlGn   # ← 颜色方案；红→黄→绿，绿=高特异性
norm_spec = mcolors.Normalize(vmin=vmin_spec, vmax=vmax_spec)

s_min = sub_top['neg_log_mag'].min()
s_max = sub_top['neg_log_mag'].max()

fig4, ax4 = plt.subplots(
    figsize=(16, 10)   # ← 图像尺寸；列多时可适当加宽
)

for _, row in sub_top.iterrows():
    lr = row['LR']
    st = row['st']
    if lr not in lr_idx or st not in st_idx:
        continue
    x = st_idx[st]
    y = lr_idx[lr]
    s = 30 + 370 * (row['neg_log_mag'] - s_min) / (s_max - s_min + 1e-10)
    # ← 点大小范围：最小 30，最大 400（30+370）；修改这两个数字调整大小范围
    rgba = list(cmap_dot(norm_spec(row['log_spec'])))
    ax4.scatter(x, y, s=s, color=rgba, alpha=0.88,   # ← alpha: 点的透明度
                linewidths=0.5, edgecolors='#444444', zorder=3)

# 背景网格线
for i in range(len(lr_order)):
    ax4.axhline(i, color='#ebebeb', linewidth=0.6, zorder=1)
for j in range(len(st_order)):
    ax4.axvline(j, color='#ebebeb', linewidth=0.6, zorder=1)

ax4.set_xticks(range(len(st_order)))
ax4.set_xticklabels(st_order, fontsize=8.5, rotation=0, ha='center', va='top')
ax4.set_yticks(range(len(lr_order)))
ax4.set_yticklabels(lr_order, fontsize=8.5)
ax4.set_xlim(-0.6, len(st_order) - 0.4)
ax4.set_ylim(-0.6, len(lr_order) - 0.4)
ax4.set_xlabel('Source → Target', fontsize=10, labelpad=12)
ax4.set_ylabel('Ligand → Receptor', fontsize=10)
ax4.set_title(
    'Top Ligand-Receptor Interactions: B Cells → Immune Partners\n'
    'DLBCL — LIANA rank_aggregate  |  Color: Specificity weight  |  Size: Magnitude',
    fontsize=11, fontweight='bold', pad=12
)

# Colorbar
sm4 = ScalarMappable(cmap=cmap_dot, norm=norm_spec)
sm4.set_array([])
cbar4 = plt.colorbar(sm4, ax=ax4, shrink=0.4, pad=0.02, aspect=18)
cbar4.set_label('Specificity weight\n(log₁₀)', fontsize=9)
cbar4.ax.tick_params(labelsize=8)

# 大小图例
for pct, label in [(0.1, 'Low'), (0.5, 'Medium'), (1.0, 'High')]:
    v  = s_min + pct * (s_max - s_min)
    sz = 30 + 370 * (v - s_min) / (s_max - s_min + 1e-10)
    ax4.scatter([], [], s=sz, color='#888888', alpha=0.85,
                edgecolors='#444444', linewidths=0.5, label=label)
ax4.legend(title='Magnitude', title_fontsize=9, fontsize=8.5,
           frameon=True, framealpha=0.9, loc='lower right',
           bbox_to_anchor=(1.22, 0.0))

for spine in ['top', 'right']:
    ax4.spines[spine].set_visible(False)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig4.savefig(
        os.path.join(OUTPUT_DIR, f'dotplot_top_lr.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )
plt.close(fig4)
print("  ✓ dotplot_top_lr saved")


# ══════════════════════════════════════════════════════════════════════════════
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
    "legend_loc"              : 'upper left',
    "legend_fontsize"         : 10,
    "legend_title"            : 'Category',
    "legend_title_fontsize"   : 13,
    "legend_bbox"             : (0.76, 0.32),
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


# =============================================================================
#  图3 + 图5 + 图7：Malignant vs Normal B cell 通讯可视化
#  数据文件：liana_results_13subtypes_full.csv
#  注意：该文件中细胞类型名称为 'Malignant B' 和 'Normal B'
#        （不含 "cell"，与之前 liana_results_full.csv 中的命名不同）
# =============================================================================

import matplotlib
matplotlib.use('Agg')   # 非交互后端；本地弹窗预览改为 'TkAgg'

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import seaborn as sns
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
#  ① 路径 & 全局参数
# =============================================================================
DATA_DIR   = translate(r'D:/bulk-download/GSE182434/cellcomm')
OUTPUT_DIR = translate(r'D:/bulk-download/GSE182434/cellcomm')
os.makedirs(OUTPUT_DIR, exist_ok=True)

DPI = 300   # ← 全局 DPI；Agent 原始值为 150

# =============================================================================
#  ② 加载数据
#  注意：该 CSV 中 Malignant/Normal B cell 的名称为 'Malignant B' / 'Normal B'
# =============================================================================
liana_res_mal = pd.read_csv(
    os.path.join(DATA_DIR, 'liana_results_13subtypes_full.csv')
)
print(f"Loaded: {liana_res_mal.shape}")
print("All sources:", sorted(liana_res_mal['source'].unique()))

# 确认关键细胞类型存在
for name in ['Malignant B', 'Normal B']:
    assert name in liana_res_mal['source'].unique() or name in liana_res_mal['target'].unique(), \
        f"'{name}' not found in data!"
print("✓ Malignant B and Normal B confirmed in data.")

# =============================================================================
#  ③ 配色 & 短标签（Malignant vs Normal 专用）
# =============================================================================
ct_palette_mal = {
    'Malignant B':            '#d62728',   # ← 恶性 B 细胞颜色
    'Normal B':               '#2ca02c',   # ← 正常 B 细胞颜色
    'CD8 T':                  '#1f77b4',
    'CD4 T':                  '#6baed6',
    'Tregs':                  '#17becf',
    'TFH':                    '#31a354',
    'NK':                     '#a1d99b',
    'Mono/Mac':               '#e377c2',
    'pDC':                    '#636363',
}

# 细胞类型全名 → 短名映射（用于 Circle Plot 节点标签）
# 注意：CSV 中的名称已经较短，这里做直通映射，需要与 source/target 列一致
ct_short_mal = {
    'Malignant B':            'Malignant B',
    'Normal B':               'Normal B',
    'CD8 T':                  'CD8 T',
    'CD4 T':                  'CD4 T',
    'Tregs':                  'Tregs',
    'TFH':                    'TFH',
    'NK':                     'NK',
    'Mono':                   'Mono/Mac',
    'pDC':                    'pDC',
    'LA_TAM':                 'LA-TAM',
    'IFN_TAM':                'IFN-TAM',
    'DC_1':                   'DC1',
    'DC_2':                   'DC2',
}

# Circle Plot 节点颜色（按 CSV 中的实际细胞类型名称）
ct_palette_circle = {
    'Malignant B':  '#d62728',
    'Normal B':     '#2ca02c',
    'CD8 T':        '#1f77b4',
    'CD4 T':        '#6baed6',
    'Tregs':        '#17becf',
    'TFH':          '#31a354',
    'NK':           '#a1d99b',
    'Mono':         '#e377c2',
    'pDC':          '#636363',
    'LA_TAM':       '#8c564b',
    'IFN_TAM':      '#ff7f0e',
    'DC_1':         '#9467bd',
    'DC_2':         '#c5b0d5',
}

# =============================================================================
#  图3：Circle Plot — Malignant vs Normal B cell 通讯网络
# =============================================================================
print("\n[Fig 3] Drawing circle communication plot (Malignant vs Normal)...")

def draw_circle_comm(liana_df, palette, short_labels, title,
                     top_n=200,          # ← 取交互强度 Top N 条；增大显示更多边
                     figsize=(13, 13),   # ← 图像尺寸（宽, 高）英寸
                     save_path=None):
    df = liana_df.sort_values('magnitude_rank').head(top_n).copy()
    df['weight'] = 1 - df['magnitude_rank']

    pair_df = (df.groupby(['source', 'target'])
               .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
               .reset_index())

    active_cts = [c for c in palette
                  if c in set(pair_df['source']) | set(pair_df['target'])]
    n_ct = len(active_cts)
    if n_ct == 0:
        print("  WARNING: no active cell types found, skipping.")
        return

    node_tot = {ct: (pair_df[pair_df['source']==ct]['total_weight'].sum() +
                     pair_df[pair_df['target']==ct]['total_weight'].sum())
                for ct in active_cts}

    angles = {ct: np.pi/2 - 2*np.pi*i/n_ct for i, ct in enumerate(active_cts)}
    pos    = {ct: (np.cos(angles[ct]), np.sin(angles[ct])) for ct in active_cts}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    w_vals = pair_df['total_weight'].values
    w_min, w_max = w_vals.min(), w_vals.max()

    # ── 绘制边（箭头）──────────────────────────────────────────────────────
    for _, row in pair_df.sort_values('total_weight').iterrows():
        src, tgt = row['source'], row['target']
        if src not in pos or tgt not in pos or src == tgt:
            continue
        w_norm = (row['total_weight'] - w_min) / (w_max - w_min + 1e-10)

        lw    = (0.4 + 5.6 * w_norm) * 3.0   # ← 边宽系数；3.0 = ×3 加粗
        rad   = 0.2 + 0.15 * w_norm           # ← 弧度弯曲程度
        color = palette.get(src, '#888888')
        hw    = 0.015 + 0.025 * w_norm        # ← 箭头头部宽度
        hl    = 0.025 + 0.030 * w_norm        # ← 箭头头部长度

        ax.annotate(
            '', xy=pos[tgt], xytext=pos[src],
            arrowprops=dict(
                arrowstyle=f'-|>,head_width={hw:.3f},head_length={hl:.3f}',
                connectionstyle=f'arc3,rad={rad:.2f}',
                color=color,
                lw=lw,
                alpha=1.0,   # ← 边透明度；1.0 = 完全不透明
            ),
            zorder=2,
        )

    # ── 绘制节点 ────────────────────────────────────────────────────────────
    ns_max = max(node_tot.values()) if node_tot else 1
    for ct in active_cts:
        x, y = pos[ct]
        ns_norm = node_tot.get(ct, 0) / ns_max
        r = (0.055 + 0.095 * ns_norm) * 0.75   # ← 0.75 = 节点半径缩放系数
        circle = plt.Circle((x, y), r,
                             color=palette.get(ct, '#888888'),
                             zorder=5,
                             linewidth=1.8,
                             edgecolor='white')
        ax.add_patch(circle)

    # ── 绘制标签 ────────────────────────────────────────────────────────────
    label_r = 1.20   # ← 标签距圆心距离
    for ct in active_cts:
        ang = angles[ct]
        deg = np.degrees(ang) % 360
        ha  = 'right' if 80 < deg < 280 else 'left'
        lx  = (label_r + 0.08) * np.cos(ang)
        ly  = (label_r + 0.08) * np.sin(ang)

        ax.text(lx, ly, short_labels.get(ct, ct),
                ha=ha, va='center',
                fontsize=8.5,       # ← 标签字号
                fontweight='bold',
                color=palette.get(ct, '#333333'),
                path_effects=[pe.withStroke(linewidth=3.0, foreground='white')],
                zorder=6)

        nx, ny = pos[ct]
        ax.plot([nx * 1.07, label_r * np.cos(ang)],
                [ny * 1.07, label_r * np.sin(ang)],
                color=palette.get(ct, '#888888'),
                linewidth=0.6, alpha=0.5, zorder=3)

    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=14)

    legend_els = [
        mpatches.Patch(color='#555555', label='Strong interaction (thick edge)'),
        mpatches.Patch(color='#aaaaaa', label='Weak interaction (thin edge)'),
        plt.scatter([], [], s=100, c='#555555', label='High comm. (large node)'),
        plt.scatter([], [], s=30,  c='#555555', label='Low comm. (small node)'),
    ]
    ax.legend(handles=legend_els, loc='lower left', fontsize=7.5,
              frameon=True, framealpha=0.9,
              bbox_to_anchor=(-0.08, -0.08),
              title='Edge / Node encoding', title_fontsize=8)

    if save_path:
        for fmt in ['png', 'svg']:
            fig.savefig(f'{save_path}.{fmt}',
                        dpi=DPI, bbox_inches='tight', format=fmt)
        plt.close(fig)
        print(f"  ✓ Saved: {save_path}")
    return fig


draw_circle_comm(
    liana_res_mal, ct_palette_circle, ct_short_mal,
    title='Cell-Cell Communication Network (Malignant vs Normal B cells)\nDLBCL — LIANA rank_aggregate (Top 200 interactions)',
    top_n=200,          # ← 展示最强的前 200 条交互
    figsize=(13, 13),
    save_path=os.path.join(OUTPUT_DIR, 'circle_plot_malignant_vs_normal_v2')
)

# =============================================================================
#  图5：Dotplot — Malignant/Normal B cell × Immune Partners（透视热图风格）
# =============================================================================
# ══════════════════════════════════════════════════════════════════════════════
# Fig5 前置参数：Dotplot — Malignant vs Normal B cell interactions
# ══════════════════════════════════════════════════════════════════════════════

P05 = {

    # ── 数据筛选 ──────────────────────────────────────────────────────────────
    # CSV 中 Malignant/Normal B 细胞的实际名称
    "b_cell_types"      : ['Malignant B', 'Normal B'],
    # 展示的 source→target 方向列表（不在数据中的方向自动过滤）
    "directions"        : [
        'Malignant B → CD8 T',
        'Malignant B → CD4 T',
        'Malignant B → Mono',
        'Malignant B → NK',
        'Malignant B → Tregs',
        'Malignant B → TFH',
        'CD8 T → Malignant B',
        'CD4 T → Malignant B',
        'Mono → Malignant B',
        'NK → Malignant B',
        'Normal B → CD8 T',
        'Normal B → CD4 T',
        'Normal B → Mono',
        'CD8 T → Normal B',
        'CD4 T → Normal B',
        'Mono → Normal B',
    ],
    # 展示 Top N 个 LR pair（增大可显示更多行）
    "top_n_lr"          : 30,

    # ── 图像尺寸（自动根据行列数计算，此为每格基准尺寸）────────────────────
    "col_width"         : 1.1,    # 每列宽度（英寸）
    "col_padding"       : 4.0,    # 列方向额外留白（英寸，给 y 轴标签）
    "row_height"        : 0.52,   # 每行高度（英寸）
    "row_padding"       : 3.0,    # 行方向额外留白（英寸，给标题/图例）

    # ── 点大小 ────────────────────────────────────────────────────────────────
    "dot_size_base"     : 60,     # 最小点面积（specificity=0 时）
    "dot_size_range"    : 300,    # 点面积变化幅度（specificity 从 0→1 增加的量）
    "dot_size_scale"    : 3.0,    # 整体缩放系数（同比放大/缩小所有点）

    # ── 点边框 ────────────────────────────────────────────────────────────────
    "dot_edge_color"    : '#333333',  # 点边框颜色
    "dot_edge_lw"       : 0.5,        # 点边框线宽

    # ── 颜色映射 ──────────────────────────────────────────────────────────────
    "cmap"              : 'plasma',   # 颜色方案（可换 'viridis'/'magma'/'RdYlGn'）
    "gamma"             : 0.25,       # PowerNorm gamma；< 1 拉伸低端颜色差异
    # score 分位数裁剪（vmin 取第 5 百分位，避免极端值影响颜色映射）
    "score_vmin_pct"    : 5,

    # ── 网格线 ────────────────────────────────────────────────────────────────
    "grid_color"        : '#eeeeee',  # 网格线颜色
    "grid_lw"           : 1,        # 网格线宽度

    # ── 分隔线（Malignant / Normal B 列组之间）───────────────────────────────
    "sep_line_color"    : '#333333',  # 分隔线颜色
    "sep_line_lw"       : 1.5,        # 分隔线宽度
    "sep_line_style"    : '--',        # 分隔线样式
    "sep_line_alpha"    : 0.7,         # 分隔线透明度

    # ── 分组标题（列组上方文字）──────────────────────────────────────────────
    "group_title_fontsize"  : 12,
    "group_title_y_offset"  : 0.15,   # 标题距热图顶部的偏移（数据坐标单位）
    "mal_title_color"       : '#d62728',   # Malignant B 组标题颜色
    "nor_title_color"       : '#2ca02c',   # Normal B 组标题颜色
    "mal_title_text"        : 'Malignant B interactions',
    "nor_title_text"        : 'Normal B interactions',

    # ── colorbar ──────────────────────────────────────────────────────────────
    "cbar_shrink"       : 0.45,    # colorbar 相对高度（0~1）
    "cbar_pad"          : 0.01,    # colorbar 与主图间距
    "cbar_aspect"       : 20,      # colorbar 长宽比
    "cbar_n_ticks"      : 7,       # colorbar 刻度数量
    "cbar_label"        : 'Interaction score\n(1 − magnitude_rank)',
    "cbar_label_fontsize"   : 11,
    "cbar_tick_fontsize"    : 10,

    # ── 点大小图例（Specificity）─────────────────────────────────────────────
    # 每项：(specificity 示例值, 标签文字)
    "legend_spec_examples"  : [(0.2, 'Low'), (0.55, 'Medium'), (0.8, 'High')],
    "legend_dot_color"      : '#aaaaaa',   # 图例点填充色
    "legend_fontsize"       : 11,
    "legend_title_fontsize" : 12,
    "legend_title"          : 'Specificity\n(dot size)',
    "legend_loc"            : 'upper left',
    "legend_bbox"           : (1.10, 1.0),  # 图例锚点位置
    "legend_framealpha"     : 0.9,

    # ── 轴刻度标签 ────────────────────────────────────────────────────────────
    "xtick_fontsize"    : 12,
    "xtick_rotation"    : 40,      # x 轴标签旋转角度
    "ytick_fontsize"    : 12,

    # ── 轴标签 ────────────────────────────────────────────────────────────────
    "xlabel"            : 'Interaction direction (Source → Target)',
    "xlabel_fontsize"   : 12,
    "xlabel_labelpad"   : 8,
    "ylabel"            : 'Ligand → Receptor pair',
    "ylabel_fontsize"   : 12,
    "ylabel_labelpad"   : 8,

    # ── 图标题 ────────────────────────────────────────────────────────────────
    "title"             : (''),
    "title_fontsize"    : 11,
    "title_pad"         : 22,

    # ── 输出文件名 ────────────────────────────────────────────────────────────
    "output_name"       : 'dotplot_malignant_vs_normal_v2',
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig5 绘图
# ══════════════════════════════════════════════════════════════════════════════

print("\n[Fig 5] Drawing dotplot (Malignant vs Normal B cells, pivot style)...")

# ── 筛选数据 ──────────────────────────────────────────────────────────────────
df_b = liana_res_mal[
    liana_res_mal['source'].isin(P05["b_cell_types"]) |
    liana_res_mal['target'].isin(P05["b_cell_types"])
].copy()

df_b['lr_pair']     = df_b['ligand_complex'] + ' → ' + df_b['receptor_complex']
df_b['pair_label']  = df_b['source'] + ' → ' + df_b['target']
df_b['score']       = 1 - df_b['magnitude_rank']
df_b['specificity'] = 1 - df_b['specificity_rank']

directions = [d for d in P05["directions"] if d in df_b['pair_label'].unique()]
print(f"  Available directions: {directions}")

df_filt    = df_b[df_b['pair_label'].isin(directions)].copy()
lr_mean    = df_filt.groupby('lr_pair')['score'].mean().nlargest(P05["top_n_lr"])
top_lr_mal = lr_mean.index.tolist()
df_filt    = df_filt[df_filt['lr_pair'].isin(top_lr_mal)]

score_mat = df_filt.pivot_table(index='lr_pair', columns='pair_label',
                                values='score', aggfunc='mean')
spec_mat  = df_filt.pivot_table(index='lr_pair', columns='pair_label',
                                values='specificity', aggfunc='mean')

score_mat = score_mat.reindex(columns=directions).dropna(how='all')
spec_mat  = spec_mat.reindex(columns=directions).reindex(score_mat.index)

row_order5 = score_mat.fillna(0).sum(axis=1).sort_values(ascending=True).index
score_mat  = score_mat.loc[row_order5]
spec_mat   = spec_mat.loc[row_order5]

n_rows5, n_cols5 = score_mat.shape
print(f"  Dotplot grid: {n_rows5} LR pairs × {n_cols5} directions")

# ── 颜色映射 ──────────────────────────────────────────────────────────────────
vals5  = score_mat.values.flatten()
vals5  = vals5[~np.isnan(vals5)]
vmin5  = float(np.percentile(vals5, P05["score_vmin_pct"]))
vmax5  = float(vals5.max())
norm5  = mcolors.PowerNorm(gamma=P05["gamma"], vmin=vmin5, vmax=vmax5)
cmap5  = plt.cm.get_cmap(P05["cmap"])

# ── 绘图 ──────────────────────────────────────────────────────────────────────
fig5, ax5 = plt.subplots(
    figsize=(
        n_cols5 * P05["col_width"] + P05["col_padding"],
        n_rows5 * P05["row_height"] + P05["row_padding"]
    )
)

for i, lr in enumerate(score_mat.index):
    for j, direction in enumerate(score_mat.columns):
        sv  = score_mat.loc[lr, direction]
        spv = spec_mat.loc[lr, direction]
        if pd.isna(sv):
            continue
        color = cmap5(norm5(np.clip(sv, vmin5, vmax5)))
        size  = (
            P05["dot_size_base"]
            + P05["dot_size_range"] * (spv if not pd.isna(spv) else 0)
        ) * P05["dot_size_scale"]
        ax5.scatter(j, i, s=size, color=color,
                    edgecolors=P05["dot_edge_color"],
                    linewidths=P05["dot_edge_lw"],
                    zorder=3)

# ── 轴刻度 ────────────────────────────────────────────────────────────────────
ax5.set_xticks(range(n_cols5))
ax5.set_xticklabels(score_mat.columns,
                    fontsize=P05["xtick_fontsize"],
                    rotation=P05["xtick_rotation"], ha='right')
ax5.set_yticks(range(n_rows5))
ax5.set_yticklabels(score_mat.index, fontsize=P05["ytick_fontsize"])
ax5.set_xlim(-0.5, n_cols5 - 0.5)
ax5.set_ylim(-0.5, n_rows5 - 0.5)

# ── 网格线 ────────────────────────────────────────────────────────────────────
ax5.grid(True, color=P05["grid_color"],
         linewidth=P05["grid_lw"], zorder=0)
ax5.set_axisbelow(True)

# ── 分隔线 & 分组标题 ─────────────────────────────────────────────────────────
mal_cols5 = [j for j, d in enumerate(score_mat.columns) if 'Malignant' in d]
nor_cols5  = [j for j, d in enumerate(score_mat.columns) if 'Normal'   in d]

if mal_cols5 and nor_cols5:
    sep5 = (max(mal_cols5) + min(nor_cols5)) / 2
    ax5.axvline(sep5,
                color=P05["sep_line_color"],
                linewidth=P05["sep_line_lw"],
                linestyle=P05["sep_line_style"],
                alpha=P05["sep_line_alpha"])

if mal_cols5:
    ax5.text(np.mean(mal_cols5),
             n_rows5 + P05["group_title_y_offset"],
             P05["mal_title_text"],
             ha='center', va='bottom',
             fontsize=P05["group_title_fontsize"],
             color=P05["mal_title_color"],
             fontweight='bold')

if nor_cols5:
    ax5.text(np.mean(nor_cols5),
             n_rows5 + P05["group_title_y_offset"],
             P05["nor_title_text"],
             ha='center', va='bottom',
             fontsize=P05["group_title_fontsize"],
             color=P05["nor_title_color"],
             fontweight='bold')

# ── colorbar ──────────────────────────────────────────────────────────────────
sm5 = plt.cm.ScalarMappable(cmap=cmap5, norm=norm5)
sm5.set_array([])
cbar5 = fig5.colorbar(sm5, ax=ax5,
                      shrink=P05["cbar_shrink"],
                      pad=P05["cbar_pad"],
                      aspect=P05["cbar_aspect"])
cbar5.set_label(P05["cbar_label"], fontsize=P05["cbar_label_fontsize"])
tick_scores5 = np.linspace(vmin5, vmax5, P05["cbar_n_ticks"])
cbar5.set_ticks([norm5(v) for v in tick_scores5])
cbar5.set_ticklabels([f'{v:.2f}' for v in tick_scores5],
                     fontsize=P05["cbar_tick_fontsize"])

# ── 点大小图例（Specificity）─────────────────────────────────────────────────
for spec_ex, label in P05["legend_spec_examples"]:
    ax5.scatter([], [],
                s=(P05["dot_size_base"] + P05["dot_size_range"] * spec_ex)
                  * P05["dot_size_scale"],
                color=P05["legend_dot_color"],
                edgecolors=P05["dot_edge_color"],
                linewidths=P05["dot_edge_lw"],
                label=label)

ax5.legend(title=P05["legend_title"],
           fontsize=P05["legend_fontsize"],
           title_fontsize=P05["legend_title_fontsize"],
           loc=P05["legend_loc"],
           bbox_to_anchor=P05["legend_bbox"],
           frameon=True,
           framealpha=P05["legend_framealpha"])

# ── 标题 & 轴标签 ─────────────────────────────────────────────────────────────
ax5.set_title(P05["title"],
              fontsize=P05["title_fontsize"],
              fontweight='bold',
              pad=P05["title_pad"])
ax5.set_xlabel(P05["xlabel"],
               fontsize=P05["xlabel_fontsize"],
               labelpad=P05["xlabel_labelpad"])
ax5.set_ylabel(P05["ylabel"],
               fontsize=P05["ylabel_fontsize"],
               labelpad=P05["ylabel_labelpad"])

plt.tight_layout()

for fmt in ['png', 'svg']:
    fig5.savefig(
        os.path.join(OUTPUT_DIR, f'{P05["output_name"]}.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )

plt.close(fig5)
print("  ✓ dotplot_malignant_vs_normal_v2 saved")

# =============================================================================
#  图7：Dotplot — Malignant/Normal B cell → Immune Partners（散点矩阵风格）
# =============================================================================
# ══════════════════════════════════════════════════════════════════════════════
# Fig7 前置参数：Dotplot — Malignant vs Normal B → Immune Partners
# ══════════════════════════════════════════════════════════════════════════════
DATA_DIR_7   = translate(r'D:/bulk-download/GSE182434')
liana_res_mal_7 = pd.read_csv(
    os.path.join(DATA_DIR_7, 'liana_res_mal.csv')
)
P07 = {

    # ── 数据筛选 ──────────────────────────────────────────────────────────────
    # 信号发出方（source）细胞类型
    "b_sources"             : ['Malignant B cell', 'Normal B cell'],
    # 信号接收方（target）细胞类型（需与 CSV 中 target 列一致）
    "key_targets"           : ['CD4 T cells', 'CD8 T cells', 'Tregs', 'TFH cells', 'NK cells', 'Monocytes/Macrophages'],
    # 展示 Top N 个 LR pair（增大可显示更多行）
    "top_n_lr"              : 20,
    # source 短标签映射（key=CSV 中名称，value=显示名称）
    "src_short"             : {'Malignant B cell': 'Malignant B', 'Normal B cell': 'Normal B'},
    # target 短标签映射
    "tgt_short"             : {
        'CD4 T cells': 'CD4 T', 'CD8 T cells': 'CD8 T',
        'Tregs': 'Tregs', 'TFH cells': 'TFH',
        'NK cells': 'NK', 'Monocytes/Macrophages': 'Mono/Mac',
    },

    # ── 图像尺寸（每格基准尺寸）──────────────────────────────────────────────
    "col_width"             : 0.56,   # 每列宽度（英寸）
    "col_padding"           : 3.5,    # 列方向额外留白（英寸，给 y 轴标签）
    "row_height"            : 0.3,   # 每行高度（英寸）
    "row_padding"           : 3.5,    # 行方向额外留白（英寸，给标题/图例）

    # ── 点大小（按 -log10(magnitude_rank) 线性缩放）──────────────────────────
    "dot_size_min"          : 60,     # 最小点面积
    "dot_size_range"        : 350,    # 点面积变化幅度（最大 = min + range = 400）
    "dot_alpha"             : 0.88,   # 点透明度

    # ── 点边框 ────────────────────────────────────────────────────────────────
    "dot_edge_color"        : '#555555',  # 点边框颜色
    "dot_edge_lw"           : 0.4,        # 点边框线宽

    # ── 颜色映射（颜色编码 spec_weight）──────────────────────────────────────
    "color_col"             : 'spec_weight',  # 颜色编码列名
    "cmap"                  : 'RdYlGn',       # 颜色方案（红=低特异，绿=高特异）
    "color_vmin_pct"        : 5,              # 颜色下限分位数
    "color_vmax_pct"        : 95,             # 颜色上限分位数

    # ── 网格线 ────────────────────────────────────────────────────────────────
    "hgrid_color"           : '#e0e0e0',  # 横向网格线颜色
    "hgrid_lw"              : 1,        # 横向网格线宽度
    "vgrid_color"           : '#e0e0e0',  # 纵向网格线颜色
    "vgrid_lw"              : 1,        # 纵向网格线宽度

    # ── 分隔线（Malignant / Normal B 列组之间）───────────────────────────────
    "sep_line_color"        : '#333333',  # 分隔线颜色
    "sep_line_lw"           : 1.5,        # 分隔线宽度
    "sep_line_style"        : '--',        # 分隔线样式
    "sep_line_alpha"        : 0.7,         # 分隔线透明度

    # ── 分组标题（列组上方文字）──────────────────────────────────────────────
    "group_title_fontsize"  : 15,
    "group_title_y_offset"  : -0.3,       # 距热图顶部偏移（数据坐标单位）
    "mal_title_color"       : '#d62728',  # Malignant B 组标题颜色
    "nor_title_color"       : '#2ca02c',  # Normal B 组标题颜色
    "mal_title_text"        : 'Malignant B cell',
    "nor_title_text"        : 'Normal B cell',

    # ── colorbar ──────────────────────────────────────────────────────────────
    "cbar_shrink"           : 0.45,   # colorbar 相对高度（0~1）
    "cbar_pad"              : 0.02,   # colorbar 与主图间距
    "cbar_aspect"           : 20,     # colorbar 长宽比
    "cbar_label"            : 'Specificity score\n(spec_weight)',
    "cbar_label_fontsize"   : 12,
    "cbar_tick_fontsize"    : 11,

    # ── 点大小图例（Magnitude）────────────────────────────────────────────────
    # 每项：(magnitude 百分比位置 0~1, 标签文字)
    "legend_mag_examples"   : [(0.1, 'Low'), (0.5, 'Medium'), (1.0, 'High')],
    "legend_dot_color"      : '#888888',  # 图例点填充色
    "legend_dot_alpha"      : 0.8,
    "legend_fontsize"       : 11,
    "legend_title_fontsize" : 12,
    "legend_title"          : 'Magnitude',
    "legend_loc"            : 'lower right',
    "legend_bbox"           : (1.0, 0.0),   # 图例锚点位置
    "legend_framealpha"     : 0.9,

    # ── 轴刻度标签 ────────────────────────────────────────────────────────────
    "xtick_fontsize"        : 11,
    "xtick_rotation"        : 45,     # x 轴标签旋转角度
    "ytick_fontsize"        : 15,

    # ── 轴标签 ────────────────────────────────────────────────────────────────
    "xlabel"                : 'Source → Target',
    "xlabel_fontsize"       : 14,
    "xlabel_labelpad"       : 10,
    "ylabel"                : 'Ligand → Receptor',
    "ylabel_fontsize"       : 14,
    "ylabel_labelpad"       : 8,

    # ── 图标题 ────────────────────────────────────────────────────────────────
    "title"                 : (''),
    "title_fontsize"        : 11,
    "title_pad"             : 18,

    # ── 输出文件名 ────────────────────────────────────────────────────────────
    "output_name"           : 'dotplot_malignant_normal_toplr',
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig7 绘图
# ══════════════════════════════════════════════════════════════════════════════

print("\n[Fig 7] Drawing dotplot (Malignant vs Normal B cells, scatter matrix)...")

# ── 筛选数据 ──────────────────────────────────────────────────────────────────
sub7 = liana_res_mal_7[
    liana_res_mal_7['source'].isin(P07["b_sources"]) &
    liana_res_mal_7['target'].isin(P07["key_targets"])
].copy()

top_lr7 = (sub7.groupby(['ligand_complex', 'receptor_complex'])['magnitude_rank']
           .mean()
           .sort_values()
           .head(P07["top_n_lr"])
           .index)

sub_top7 = sub7[
    sub7.apply(lambda r: (r['ligand_complex'], r['receptor_complex']) in top_lr7, axis=1)
].copy()

sub_top7['LR_pair']   = sub_top7['ligand_complex'] + ' → ' + sub_top7['receptor_complex']
sub_top7['src_short'] = sub_top7['source'].map(P07["src_short"])
sub_top7['tgt_short'] = sub_top7['target'].map(P07["tgt_short"])
sub_top7['st_label']  = sub_top7['src_short'] + '\n→ ' + sub_top7['tgt_short']
sub_top7['neg_log_mag'] = -np.log10(sub_top7['magnitude_rank'].clip(1e-10))

print("  spec_weight stats:")
print(sub_top7[P07["color_col"]].describe())

# ── 坐标轴顺序：Malignant 列在左，Normal 列在右 ───────────────────────────────
lr_order7   = [f"{l} → {r}" for l, r in top_lr7]
mal_labels7 = sorted([s for s in sub_top7['st_label'].unique() if 'Malignant' in s])
nor_labels7 = sorted([s for s in sub_top7['st_label'].unique() if 'Normal'    in s])
st_order7   = mal_labels7 + nor_labels7

lr_idx7 = {lr: i for i, lr in enumerate(lr_order7)}
st_idx7 = {st: i for i, st in enumerate(st_order7)}

# ── 颜色映射 ──────────────────────────────────────────────────────────────────
cv7    = sub_top7[P07["color_col"]].values
vmin_c = np.percentile(cv7, P07["color_vmin_pct"])
vmax_c = np.percentile(cv7, P07["color_vmax_pct"])
cmap7  = plt.cm.get_cmap(P07["cmap"])
norm_c = Normalize(vmin=vmin_c, vmax=vmax_c)

# ── 点大小缩放函数 ────────────────────────────────────────────────────────────
mv7          = sub_top7['neg_log_mag'].values
mag_min7     = mv7.min()
mag_max7     = mv7.max()

def scale_size7(v):
    return (P07["dot_size_min"]
            + P07["dot_size_range"]
            * (v - mag_min7) / (mag_max7 - mag_min7 + 1e-10))

# ── 绘图 ──────────────────────────────────────────────────────────────────────
fig7, ax7 = plt.subplots(
    figsize=(
        len(st_order7) * P07["col_width"] + P07["col_padding"],
        len(lr_order7) * P07["row_height"] + P07["row_padding"]
    )
)

for _, row in sub_top7.iterrows():
    lr = row['LR_pair']
    st = row['st_label']
    if lr not in lr_idx7 or st not in st_idx7:
        continue
    ax7.scatter(st_idx7[st], lr_idx7[lr],
                s=scale_size7(row['neg_log_mag']),
                c=[cmap7(norm_c(row[P07["color_col"]]))],
                alpha=P07["dot_alpha"],
                linewidths=P07["dot_edge_lw"],
                edgecolors=P07["dot_edge_color"],
                zorder=3)

# ── 网格线 ────────────────────────────────────────────────────────────────────
for i in range(len(lr_order7)):
    ax7.axhline(i, color=P07["hgrid_color"],
                linewidth=P07["hgrid_lw"], zorder=1)
for j in range(len(st_order7)):
    ax7.axvline(j, color=P07["vgrid_color"],
                linewidth=P07["vgrid_lw"], zorder=1)

# ── 分隔线 & 分组标题 ─────────────────────────────────────────────────────────
if mal_labels7 and nor_labels7:
    ax7.axvline(len(mal_labels7) - 0.5,
                color=P07["sep_line_color"],
                linewidth=P07["sep_line_lw"],
                linestyle=P07["sep_line_style"],
                alpha=P07["sep_line_alpha"],
                zorder=4)
    ax7.text(np.mean(range(len(mal_labels7))),
             len(lr_order7) + P07["group_title_y_offset"],
             P07["mal_title_text"],
             ha='center', va='bottom',
             fontsize=P07["group_title_fontsize"],
             color=P07["mal_title_color"],
             fontweight='bold')
    ax7.text(len(mal_labels7) + np.mean(range(len(nor_labels7))),
             len(lr_order7) + P07["group_title_y_offset"],
             P07["nor_title_text"],
             ha='center', va='bottom',
             fontsize=P07["group_title_fontsize"],
             color=P07["nor_title_color"],
             fontweight='bold')

# ── 轴刻度 ────────────────────────────────────────────────────────────────────
ax7.set_xticks(range(len(st_order7)))
ax7.set_xticklabels(st_order7,
                    fontsize=P07["xtick_fontsize"],
                    rotation=P07["xtick_rotation"], ha='right')
ax7.set_yticks(range(len(lr_order7)))
ax7.set_yticklabels(lr_order7, fontsize=P07["ytick_fontsize"])
ax7.set_xlim(-0.6, len(st_order7) - 0.4)
ax7.set_ylim(-0.6, len(lr_order7) - 0.4)

# ── 轴标签 & 标题 ─────────────────────────────────────────────────────────────
ax7.set_xlabel(P07["xlabel"],
               fontsize=P07["xlabel_fontsize"],
               labelpad=P07["xlabel_labelpad"])
ax7.set_ylabel(P07["ylabel"],
               fontsize=P07["ylabel_fontsize"],
               labelpad=P07["ylabel_labelpad"])
ax7.set_title(P07["title"],
              fontsize=P07["title_fontsize"],
              fontweight='bold',
              pad=P07["title_pad"])

# ── colorbar ──────────────────────────────────────────────────────────────────
sm7 = ScalarMappable(cmap=cmap7, norm=norm_c)
sm7.set_array([])
cbar7 = plt.colorbar(sm7, ax=ax7,
                     shrink=P07["cbar_shrink"],
                     pad=P07["cbar_pad"],
                     aspect=P07["cbar_aspect"])
cbar7.set_label(P07["cbar_label"], fontsize=P07["cbar_label_fontsize"])
cbar7.ax.tick_params(labelsize=P07["cbar_tick_fontsize"])

# ── 点大小图例（Magnitude）────────────────────────────────────────────────────
for pct, label in P07["legend_mag_examples"]:
    v = mag_min7 + pct * (mag_max7 - mag_min7)
    ax7.scatter([], [], s=scale_size7(v),
                c=P07["legend_dot_color"],
                alpha=P07["legend_dot_alpha"],
                edgecolors=P07["dot_edge_color"],
                linewidths=P07["dot_edge_lw"],
                label=label)

ax7.legend(title=P07["legend_title"],
           title_fontsize=P07["legend_title_fontsize"],
           fontsize=P07["legend_fontsize"],
           frameon=True,
           framealpha=P07["legend_framealpha"],
           loc=P07["legend_loc"],
           bbox_to_anchor=P07["legend_bbox"])

plt.tight_layout()

for fmt in ['png', 'svg']:
    fig7.savefig(
        os.path.join(OUTPUT_DIR, f'{P07["output_name"]}.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )

plt.close(fig7)
print("  ✓ dotplot_malignant_normal_toplr saved")
print("\n========== All figures saved to:", OUTPUT_DIR, "==========")


# =============================================================================
#  LIANA2 — 13 细胞亚型通讯可视化（本地复现）
#  图1: 全局 Circle Plot（13 subtypes）
#  图2: Interaction Weight Heatmap（13 subtypes）
#  图3-8: Focused Circle Plot × 6（Malignant B / LA_TAM / CD8 T，各 sender/receiver）
#
#  数据文件: liana_results_13subtypes_full.csv
#  路径规则: /mnt/results/GSE182434/cellcomm/ → D:/bulk-download/GSE182434/cellcomm/
# =============================================================================

import matplotlib
matplotlib.use('Agg')   # 非交互后端；本地弹窗预览改为 'TkAgg'

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import seaborn as sns
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
#  ① 路径 & 全局参数
# =============================================================================
DATA_DIR   = translate(r'D:/bulk-download/GSE182434/cellcomm')   # CSV 数据目录
OUTPUT_DIR = translate(r'D:/bulk-download/GSE182434/cellcomm')   # 图片输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

DPI = 300   # ← 全局 DPI；Agent 原始值为 150

# =============================================================================
#  ② 加载数据
# =============================================================================
liana_res = pd.read_csv(
    os.path.join(DATA_DIR, 'liana_results_13subtypes_full.csv')
)
print(f"Loaded: {liana_res.shape}")
print("Sources:", sorted(liana_res['source'].unique()))

# =============================================================================
#  ③ 公共配色（13 亚型）
#  ← 修改右侧十六进制颜色值可改变所有图中对应细胞类型的颜色
# =============================================================================
ct_palette = {
    'Malignant B':  '#d62728',
    'Normal B':     '#2ca02c',
    'CD8 T':        '#1f77b4',
    'CD4 T':        '#6baed6',
    'Tregs':        '#17becf',
    'TFH':          '#31a354',
    'NK':           '#a1d99b',
    'Mono':         '#3498DB',
    'DC_1':         '#E74C3C',
    'LA_TAM':       '#2ECC71',
    'IFN_TAM':      '#E67E22',
    'DC_2':         '#9B59B6',
    'pDC':          '#7f7f7f',
}
ct_short = {ct: ct for ct in ct_palette}   # 标签已足够短，直通映射

# =============================================================================
#  图1：全局 Circle Plot — 13 细胞亚型通讯网络
# =============================================================================
print("\n[Fig 1] Drawing global circle plot (13 subtypes)...")

def draw_circle_comm_v3(liana_df, palette, short_labels, title,
                        top_n=200,          # ← 取交互强度 Top N 条；增大显示更多边
                        figsize=(13, 13),   # ← 图像尺寸（宽, 高）英寸
                        save_path=None):
    df = liana_df.sort_values('magnitude_rank').head(top_n).copy()
    df['weight'] = 1 - df['magnitude_rank']

    pair_df = (df.groupby(['source', 'target'])
               .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
               .reset_index())

    active_cts = [c for c in palette
                  if c in set(pair_df['source']) | set(pair_df['target'])]
    n_ct = len(active_cts)

    node_tot = {ct: (pair_df[pair_df['source']==ct]['total_weight'].sum() +
                     pair_df[pair_df['target']==ct]['total_weight'].sum())
                for ct in active_cts}

    # 节点均匀分布在圆上，从顶部顺时针排列
    angles = {ct: np.pi/2 - 2*np.pi*i/n_ct for i, ct in enumerate(active_cts)}
    pos    = {ct: (np.cos(angles[ct]), np.sin(angles[ct])) for ct in active_cts}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    w_vals = pair_df['total_weight'].values
    w_min, w_max = w_vals.min(), w_vals.max()

    # ── 绘制边（箭头）──────────────────────────────────────────────────────
    for _, row in pair_df.sort_values('total_weight').iterrows():
        src, tgt = row['source'], row['target']
        if src not in pos or tgt not in pos or src == tgt:
            continue
        w_norm = (row['total_weight'] - w_min) / (w_max - w_min + 1e-10)

        lw    = (0.4 + 5.6 * w_norm) * 3.0   # ← 边宽系数；3.0 = ×3 加粗
        rad   = 0.2 + 0.15 * w_norm           # ← 弧度弯曲程度
        color = palette.get(src, '#888888')
        hw    = 0.015 + 0.025 * w_norm        # ← 箭头头部宽度
        hl    = 0.025 + 0.030 * w_norm        # ← 箭头头部长度

        ax.annotate(
            '', xy=pos[tgt], xytext=pos[src],
            arrowprops=dict(
                arrowstyle=f'-|>,head_width={hw:.3f},head_length={hl:.3f}',
                connectionstyle=f'arc3,rad={rad:.2f}',
                color=color, lw=lw, alpha=1.0,   # ← alpha: 边透明度
            ),
            zorder=2,
        )

    # ── 绘制节点 ────────────────────────────────────────────────────────────
    ns_max = max(node_tot.values()) if node_tot else 1
    for ct in active_cts:
        x, y = pos[ct]
        ns_norm = node_tot.get(ct, 0) / ns_max
        r = (0.055 + 0.095 * ns_norm) * 0.75   # ← 0.75 = 节点半径缩放系数
        circle = plt.Circle((x, y), r,
                             color=palette.get(ct, '#888888'),
                             zorder=5, linewidth=1.8, edgecolor='white')
        ax.add_patch(circle)

    # ── 绘制标签 ────────────────────────────────────────────────────────────
    label_r = 1.22   # ← 标签距圆心距离；增大可避免与节点重叠
    for ct in active_cts:
        ang = angles[ct]
        deg = np.degrees(ang) % 360
        ha  = 'right' if 80 < deg < 280 else 'left'
        lx  = (label_r + 0.08) * np.cos(ang)
        ly  = (label_r + 0.08) * np.sin(ang)

        ax.text(lx, ly, short_labels.get(ct, ct),
                ha=ha, va='center',
                fontsize=20,          # ← 标签字号
                fontweight='bold',
                color=palette.get(ct, '#333333'),
                path_effects=[pe.withStroke(linewidth=3.0, foreground='white')],
                zorder=6)
        nx, ny = pos[ct]
        ax.plot([nx * 1.07, label_r * np.cos(ang)],
                [ny * 1.07, label_r * np.sin(ang)],
                color=palette.get(ct, '#888888'),
                linewidth=0.6, alpha=0.5, zorder=3)

    ax.set_xlim(-1.70, 1.70)
    ax.set_ylim(-1.70, 1.70)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=14)

    legend_els = [
        mpatches.Patch(color='#555555', label='Strong interaction (thick edge)'),
        mpatches.Patch(color='#aaaaaa', label='Weak interaction (thin edge)'),
        plt.scatter([], [], s=100, c='#555555', label='High comm. (large node)'),
        plt.scatter([], [], s=30,  c='#555555', label='Low comm. (small node)'),
    ]
    ax.legend(handles=legend_els, loc='lower left', fontsize=11,
              frameon=True, framealpha=0.9,
              bbox_to_anchor=(-0.08, -0.08),
              title='Edge / Node encoding', title_fontsize=12)

    if save_path:
        for fmt in ['png', 'svg']:
            fig.savefig(f'{save_path}.{fmt}', dpi=DPI, bbox_inches='tight', format=fmt)
        plt.close(fig)
        print(f"  ✓ Saved: {os.path.basename(save_path)}")
    return fig


draw_circle_comm_v3(
    liana_res, ct_palette, ct_short,
    title='Cell-Cell Communication Network\nDLBCL — 13 Cell Subtypes (LIANA rank_aggregate, Top 200)',
    top_n=200,          # ← 展示最强的前 200 条交互；增大可显示更多弱交互
    figsize=(13, 13),
    save_path=os.path.join(OUTPUT_DIR, 'circle_plot_13subtypes')
)

# =============================================================================
#  图2：Interaction Weight Heatmap — 13 细胞亚型
# =============================================================================
print("\n[Fig 2] Drawing interaction weight heatmap (13 subtypes)...")

top_n_heat = 300   # ← 用于构建热图的 Top N 交互；增大纳入更多弱交互
df_heat = liana_res.sort_values('magnitude_rank').head(top_n_heat).copy()
df_heat['weight'] = 1 - df_heat['magnitude_rank']

pair_agg_heat = (df_heat.groupby(['source', 'target'])
                 .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
                 .reset_index())

# 细胞类型排列顺序：B细胞在前，T/NK 居中，髓系在后
# ← 修改此列表可调整热图行/列的顺序
order_heat = [
    'Malignant B', 'Normal B',
    'CD8 T', 'CD4 T', 'Tregs', 'TFH', 'NK',
    'Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2', 'pDC',
]
order_heat = [c for c in order_heat
              if c in set(pair_agg_heat['source']) | set(pair_agg_heat['target'])]

mat = pd.DataFrame(0.0, index=order_heat, columns=order_heat)
for _, row in pair_agg_heat.iterrows():
    if row['source'] in order_heat and row['target'] in order_heat:
        mat.loc[row['source'], row['target']] = row['total_weight']

fig2, ax2 = plt.subplots(
    figsize=(12, 10)   # ← 图像尺寸；列多时可适当加宽
)

# 对角线 mask（自环不显示）
mask_diag = np.eye(len(order_heat), dtype=bool)

sns.heatmap(
    mat, mask=mask_diag,
    cmap='YlOrRd',          # ← 颜色方案；可换 'Blues', 'RdPu', 'viridis' 等
    linewidths=0.5,          # ← 格子分割线宽度
    linecolor='#dddddd',     # ← 分割线颜色
    annot=True,              # ← True = 格子内显示数值；False = 不显示
    fmt='.2f',               # ← 数值格式；'.1f' 保留1位小数
    annot_kws={'size': 7.5}, # ← 格子内数字字号
    ax=ax2,
    cbar_kws={
        'label': 'Aggregated interaction weight\n(sum of 1 − magnitude_rank)',
        'shrink': 0.7,   # ← colorbar 相对高度缩放
    },
    square=True,
)

ax2.set_title(
    'Cell-Cell Communication Strength\n'
    'DLBCL — 13 Cell Subtypes (LIANA rank_aggregate, Top 300)',
    fontsize=12, fontweight='bold', pad=14
)
ax2.set_xlabel('Target cell type', fontsize=10, labelpad=8)
ax2.set_ylabel('Source cell type', fontsize=10, labelpad=8)
ax2.tick_params(axis='x', rotation=40, labelsize=9)   # ← X 轴刻度旋转角度 & 字号
ax2.tick_params(axis='y', rotation=0,  labelsize=9)

# 高亮 Malignant/Normal B 的行列标签颜色
b_colors = {'Malignant B': '#d62728', 'Normal B': '#2ca02c'}
for tick in ax2.get_xticklabels():
    if tick.get_text() in b_colors:
        tick.set_color(b_colors[tick.get_text()])
        tick.set_fontweight('bold')
for tick in ax2.get_yticklabels():
    if tick.get_text() in b_colors:
        tick.set_color(b_colors[tick.get_text()])
        tick.set_fontweight('bold')

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig2.savefig(
        os.path.join(OUTPUT_DIR, f'heatmap_13subtypes.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )
plt.close(fig2)
print("  ✓ heatmap_13subtypes saved")

# =============================================================================
#  图3-8：Focused Circle Plot × 6
#  每张图聚焦一个细胞类型，展示其作为 sender 或 receiver 时的通讯格局
#  节点均匀分布在完整圆上，focal 节点在顶部并加粗高亮
# =============================================================================
print("\n[Fig 3-8] Drawing 6 focused circle plots...")

def label_align(ang_rad):
    """根据角度返回 (ha, va, name_dy, stat_dy) 标签对齐参数"""
    deg = np.degrees(ang_rad) % 360
    if deg < 20 or deg > 340:
        return 'center', 'bottom', +0.00, -0.16
    elif 20 <= deg <= 80:
        return 'left',   'center', +0.08, -0.10
    elif 80 < deg <= 100:
        return 'left',   'center', +0.00, -0.14
    elif 100 < deg <= 160:
        return 'left',   'center', -0.06, -0.14
    elif 160 < deg <= 200:
        return 'center', 'top',    +0.00, +0.14
    elif 200 < deg <= 260:
        return 'right',  'center', -0.06, -0.14
    elif 260 < deg <= 280:
        return 'right',  'center', +0.00, -0.14
    else:
        return 'right',  'center', +0.08, -0.10


def draw_focused_circle_full(liana_df, focal_ct, role, palette,
                              top_n=100,          # ← 取 Top N 条交互；增大显示更多伙伴
                              figsize=(13, 13),   # ← 图像尺寸
                              save_path=None):
    """
    聚焦式圆形通讯图。
    focal 节点置于顶部，所有伙伴节点均匀分布在完整圆上。
    role='receiver'：箭头指向 focal；role='sender'：箭头从 focal 出发。
    """
    df = liana_df.copy()

    if role == 'receiver':
        df = df[df['target'] == focal_ct].sort_values('magnitude_rank').head(top_n).copy()
        partner_col = 'source'
        arrow_note  = f'Arrows: other cells  \u2192  {focal_ct}'
        role_str    = 'Receiver'
    else:
        df = df[df['source'] == focal_ct].sort_values('magnitude_rank').head(top_n).copy()
        partner_col = 'target'
        arrow_note  = f'Arrows: {focal_ct}  \u2192  other cells'
        role_str    = 'Sender'

    df['weight'] = 1 - df['magnitude_rank']

    pair_agg = (df.groupby(partner_col)
                  .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
                  .reset_index()
                  .rename(columns={partner_col: 'partner'}))
    pair_agg = pair_agg[pair_agg['partner'] != focal_ct].sort_values('total_weight', ascending=False)

    if pair_agg.empty:
        print(f"  WARNING: No interactions for {focal_ct} as {role}, skipping.")
        return

    partners   = pair_agg['partner'].tolist()
    n_partners = len(partners)
    n_total    = n_partners + 1   # focal + partners

    # focal 在顶部（90°），其余均匀分布
    all_nodes  = [focal_ct] + partners
    all_angles = [np.pi/2 - 2*np.pi*i/n_total for i in range(n_total)]
    node_pos   = {node: (np.cos(a), np.sin(a)) for node, a in zip(all_nodes, all_angles)}

    focal_pos   = node_pos[focal_ct]
    partner_pos = {p: node_pos[p] for p in partners}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-1.95, 1.95)
    ax.set_ylim(-1.95, 1.95)

    # ── 绘制边（箭头）──────────────────────────────────────────────────────
    w_vals = pair_agg['total_weight'].values
    w_min, w_max = w_vals.min(), w_vals.max()

    for _, row in pair_agg.sort_values('total_weight').iterrows():
        partner = row['partner']
        w_norm  = (row['total_weight'] - w_min) / (w_max - w_min + 1e-10)

        lw  = (0.6 + 6.4 * w_norm) * 2.0   # ← 边宽；2.0 = 全局缩放系数
        hw  = 0.018 + 0.030 * w_norm        # ← 箭头头部宽度
        hl  = 0.028 + 0.035 * w_norm        # ← 箭头头部长度
        rad = 0.30                           # ← 弧度弯曲程度；增大弧更弯

        if role == 'receiver':
            src_pos = partner_pos[partner]
            tgt_pos = focal_pos
            color   = palette.get(partner, '#888888')
        else:
            src_pos = focal_pos
            tgt_pos = partner_pos[partner]
            color   = palette.get(focal_ct, '#888888')

        ax.annotate(
            '', xy=tgt_pos, xytext=src_pos,
            arrowprops=dict(
                arrowstyle=f'-|>,head_width={hw:.3f},head_length={hl:.3f}',
                connectionstyle=f'arc3,rad={rad:.2f}',
                color=color, lw=lw, alpha=0.85,   # ← 边透明度
            ),
            zorder=2,
        )

    # ── 绘制节点 ────────────────────────────────────────────────────────────
    w_by_partner = dict(zip(pair_agg['partner'], pair_agg['total_weight']))
    n_by_partner = dict(zip(pair_agg['partner'], pair_agg['n'].astype(int)))
    w_max_node   = max(w_by_partner.values()) if w_by_partner else 1

    # 伙伴节点
    for partner, pos in partner_pos.items():
        w = w_by_partner.get(partner, 0)
        r = 0.045 + 0.075 * (w / w_max_node)   # ← 节点半径范围 0.045~0.12
        ax.add_patch(plt.Circle(pos, r,
                                color=palette.get(partner, '#888888'),
                                zorder=5, linewidth=1.5, edgecolor='white'))

    # focal 节点（更大 + 双圈高亮）
    ax.add_patch(plt.Circle(focal_pos, 0.115,
                            color=palette.get(focal_ct, '#888888'),
                            zorder=6, linewidth=3.0, edgecolor='white'))
    ax.add_patch(plt.Circle(focal_pos, 0.130,   # ← 外圈半径；增大可让高亮环更明显
                            fill=False,
                            edgecolor=palette.get(focal_ct, '#888888'),
                            linewidth=2.0, zorder=6, alpha=0.45))

    # ── 绘制标签 ────────────────────────────────────────────────────────────
    label_r = 1.32   # ← 标签圆环半径；增大可避免标签与节点重叠

    for node, ang in zip(all_nodes, all_angles):
        is_focal = (node == focal_ct)
        ha, va, name_dy, stat_dy = label_align(ang)

        lx = label_r * np.cos(ang)
        ly = label_r * np.sin(ang)

        # 节点名称
        ax.text(lx, ly + name_dy, node,
                ha=ha, va='center',
                fontsize=20 if is_focal else 17,   # ← focal 标签字号 / 伙伴标签字号
                fontweight='bold',
                color=palette.get(node, '#333333'),
                path_effects=[pe.withStroke(linewidth=3, foreground='white')],
                zorder=8)

        # 伙伴统计信息（n = LR pair 数，w = 总权重）
        if not is_focal:
            w    = w_by_partner.get(node, 0)
            n_lr = n_by_partner.get(node, 0)
            ax.text(lx, ly + stat_dy,
                    f'n={n_lr}  w={w:.1f}',
                    ha=ha, va='center',
                    fontsize=14,        # ← 统计信息字号
                    color='#555555',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')],
                    zorder=8)

        # 连接线（节点边缘 → 标签）
        nx, ny  = node_pos[node]
        node_r  = 0.115 if is_focal else (0.045 + 0.075 * (w_by_partner.get(node, 0) / w_max_node))
        ax.plot([nx + node_r * np.cos(ang), label_r * 0.94 * np.cos(ang)],
                [ny + node_r * np.sin(ang), label_r * 0.94 * np.sin(ang)],
                color=palette.get(node, '#888888'),
                linewidth=0.6, alpha=0.4, zorder=3)

    # ── 中心箭头方向注释 ────────────────────────────────────────────────────
    ax.text(0, 0, arrow_note,
            ha='center', va='center',
            fontsize=14, style='italic', color='#444444',
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white')],
            zorder=9)

    # ── 标题 ────────────────────────────────────────────────────────────────
    ax.set_title(
        f'{focal_ct} as {role_str}\n'
        f'DLBCL | LIANA rank_aggregate | Top {top_n} interactions',
        fontsize=12,       # ← 标题字号
        fontweight='bold',
        pad=10
    )

    # ── 图例 ────────────────────────────────────────────────────────────────
    legend_els = [
        mpatches.Patch(color='#555555', label='Strong interaction (thick edge)'),
        mpatches.Patch(color='#aaaaaa', label='Weak interaction (thin edge)'),
        plt.scatter([], [], s=130, c='#555555', label='High total weight (large node)'),
        plt.scatter([], [], s=25,  c='#555555', label='Low total weight (small node)'),
    ]
    ax.legend(handles=legend_els, loc='lower left', fontsize=8,
              frameon=True, framealpha=0.92,
              bbox_to_anchor=(-0.04, -0.04),
              title='Edge / Node encoding', title_fontsize=8.5)

    if save_path:
        for fmt in ['png', 'svg']:
            fig.savefig(f'{save_path}.{fmt}', dpi=DPI, bbox_inches='tight', format=fmt)
        plt.close(fig)
        print(f"  ✓ Saved: {os.path.basename(save_path)}")
    return fig


# ── 生成 8 张聚焦圆形图 ────────────────────────────────────────────────────────
# ← 修改 configs 列表可增减绘图目标；每项格式：(focal_ct, role, 输出文件名)
configs = [
    ('Malignant B', 'receiver', 'circle_MalignantB_receiver'),
    ('Malignant B', 'sender',   'circle_MalignantB_sender'),
    ('IFN_TAM',      'receiver', 'circle_IFNTAM_receiver'),
    ('IFN_TAM',      'sender',   'circle_IFNTAM_sender'),    
    ('LA_TAM',      'receiver', 'circle_LATAM_receiver'),
    ('LA_TAM',      'sender',   'circle_LATAM_sender'),
    ('CD8 T',       'receiver', 'circle_CD8T_receiver'),
    ('CD8 T',       'sender',   'circle_CD8T_sender'),
]

for focal_ct, role, fname in configs:
    draw_focused_circle_full(
        liana_res, focal_ct=focal_ct, role=role, palette=ct_palette,
        top_n=100,          # ← 每张图展示 Top N 条交互；增大可纳入更多伙伴节点
        figsize=(13, 13),
        save_path=os.path.join(OUTPUT_DIR, fname)
    )

print(f"\n========== All 8 figures saved to: {OUTPUT_DIR} ==========")
# =============================================================================
#  热图：7 种细胞，复用 focused circle 的 w 值
# =============================================================================

top_n = 100

order_heat = ['Malignant B', 'Normal B', 'Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']

mat = pd.DataFrame(0.0, index=order_heat, columns=order_heat)

# 每个 focal 只作为 sender，填自己那一行
for focal_ct in order_heat:
    df_s = (liana_res[liana_res['source'] == focal_ct]
            .sort_values('magnitude_rank').head(top_n).copy())
    df_s['weight'] = 1 - df_s['magnitude_rank']
    for tgt, w in df_s.groupby('target')['weight'].sum().items():
        if tgt in order_heat:
            mat.loc[focal_ct, tgt] = w   # 行=source，列=target

np.fill_diagonal(mat.values, np.nan)


label_colors = {
    'Malignant B': '#d62728', 'Normal B': '#2ca02c',
    'Mono': '#3498DB', 'DC_1': '#E74C3C',
    'LA_TAM': '#2ECC71', 'IFN_TAM': '#E67E22', 'DC_2': '#9B59B6',
}

fig, ax = plt.subplots(figsize=(7, 6))

sns.heatmap(
    mat,
    ax=ax,
    mask=np.eye(len(order_heat), dtype=bool),
    cmap='YlOrRd',
    linewidths=0.5,
    linecolor='white',
    annot=True,
    fmt='.1f',
    annot_kws={'size': 10},
    cbar_kws={
        'label': 'Total interaction weight\n(sum of 1 − magnitude_rank, Top 100 per focal)',
        'shrink': 0.75,
    },
    square=True,
)

for tick in ax.get_xticklabels():
    lbl = tick.get_text()
    if lbl in label_colors:
        tick.set_color(label_colors[lbl])
        tick.set_fontweight('bold')
for tick in ax.get_yticklabels():
    lbl = tick.get_text()
    if lbl in label_colors:
        tick.set_color(label_colors[lbl])
        tick.set_fontweight('bold')

ax.tick_params(axis='x', rotation=40, labelsize=10)
ax.tick_params(axis='y', rotation=0,  labelsize=10)
ax.set_xlabel('Target (Receiver)', fontsize=11, labelpad=8)
ax.set_ylabel('Source (Sender)',   fontsize=11, labelpad=8)
ax.set_title(
    'Subtype-focused summarized communication strength',
    fontsize=12, fontweight='bold', pad=14
)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig.savefig(
        os.path.join(OUTPUT_DIR, f'heatmap_7ct_focused_w.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt
    )
plt.close(fig)
print("  ✓ heatmap_7ct_focused_w saved")

#########MAFB与23个靶基因
# =============================================================================
#  Radial Layout Network  —  MAFB Regulon (RcisTarget)
#  中心节点: MAFB (TF)
#  环绕节点: 23 个 target genes
#  颜色风格与 volcano_receptor_pseudotime_corr 保持一致
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import networkx as nx

# ── 路径 ─────────────────────────────────────────────────────────────────────
DATA_DIR = translate(r'D:\bulk-download\GSE182434')
OUT_DIR  = DATA_DIR
DPI      = 150

# ── MAFB regulon targets（来自 RcisTarget CSV）───────────────────────────────
TF = 'MAFB'
TARGETS = [
    'CEP192', 'CHCHD10', 'NAGK', 'PPRC1', 'BASP1', 'MTR', 'TIMM17A',
    'AAK1', 'MGAT1', 'USP12', 'FAM20A', 'APOL3', 'GSDMD', 'ELMO1',
    'PICALM', 'RASSF4', 'IFIT3', 'NAE1', 'NAP1L4', 'PPP2R5A',
    'MRGBP', 'PPP1R18', 'CTSL',
]

# ── 颜色方案（与 volcano 图保持一致）─────────────────────────────────────────
# ── 调参说明 ─────────────────────────────────────────────────────────────────
# COL_TF      : 中心 TF 节点填充色（沿用 volcano 高亮深色 #2B2D42）
# COL_TF_EDGE : TF 节点描边色（沿用 volcano 高亮金黄 #FFB703）
# COL_TARGET  : target 节点填充色（沿用 volcano 显著正相关红色 #E63946）
# COL_EDGE    : 连边颜色
# COL_BG      : 画布背景（与 volcano 一致 #FAFAFA）
# ─────────────────────────────────────────────────────────────────────────────
COL_TF       = '#2B2D42'
COL_TF_EDGE  = '#FFB703'
COL_TARGET   = '#E63946'
COL_EDGE     = '#AAAAAA'
COL_BG       = '#FAFAFA'

# ── 构建图 ───────────────────────────────────────────────────────────────────
G = nx.DiGraph()
G.add_node(TF)
for t in TARGETS:
    G.add_node(t)
    G.add_edge(TF, t)

# ── Radial 布局：手动计算极坐标 ──────────────────────────────────────────────
# ── 调参说明 ─────────────────────────────────────────────────────────────────
# RADIUS      : 环状 target 节点距中心的半径
# START_ANGLE : 第一个 target 的起始角度（弧度），0 = 右侧，π/2 = 上方
# ─────────────────────────────────────────────────────────────────────────────
RADIUS      = 3.2
START_ANGLE = np.pi / 2   # 从顶部开始顺时针排列

n = len(TARGETS)
angles = [START_ANGLE - 2 * np.pi * i / n for i in range(n)]

pos = {TF: np.array([0.0, 0.0])}
for gene, angle in zip(TARGETS, angles):
    pos[gene] = np.array([RADIUS * np.cos(angle),
                           RADIUS * np.sin(angle)])

# ── 画布 ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 10))
fig.patch.set_facecolor('white')
ax.set_facecolor(COL_BG)
ax.set_aspect('equal')
ax.axis('off')

# ── 绘制边（带箭头）─────────────────────────────────────────────────────────
# ── 调参说明 ─────────────────────────────────────────────────────────────────
# connectionstyle : 'arc3,rad=0.08' 让边略微弯曲，避免全部压在中心点上
# arrowstyle      : '->' 单箭头；'-|>' 实心三角箭头
# node_size 传给 nx.draw_networkx_edges 的 node_size 参数用于控制箭头起止点缩进
# ─────────────────────────────────────────────────────────────────────────────
nx.draw_networkx_edges(
    G, pos, ax=ax,
    edgelist=[(TF, t) for t in TARGETS],
    edge_color=COL_EDGE,
    arrows=True,
    arrowstyle='-|>',
    arrowsize=14,
    width=1.2,
    alpha=0.7,
    node_size=[2800, *([900] * n)],      # 与下面 node_size 对应
    connectionstyle='arc3,rad=0.0',
    min_source_margin=28,
    min_target_margin=18,
)

# ── 绘制 target 节点 ─────────────────────────────────────────────────────────
nx.draw_networkx_nodes(
    G, pos, ax=ax,
    nodelist=TARGETS,
    node_color=COL_TARGET,
    node_size=900,
    alpha=0.85,
    linewidths=1.2,
    edgecolors='white',
)

# ── 绘制 TF 中心节点 ─────────────────────────────────────────────────────────
nx.draw_networkx_nodes(
    G, pos, ax=ax,
    nodelist=[TF],
    node_color=COL_TF,
    node_size=2800,
    alpha=1.0,
    linewidths=3.0,
    edgecolors=COL_TF_EDGE,
)

# ── TF 标签 ──────────────────────────────────────────────────────────────────
ax.text(0, 0, TF,
        fontsize=16, fontweight='bold', color='white',
        ha='center', va='center', zorder=10,
        path_effects=[pe.withStroke(linewidth=1.5, foreground=COL_TF)])

# ── target 标签：沿径向自动偏移，避免与节点重叠 ──────────────────────────────
# ── 调参说明 ─────────────────────────────────────────────────────────────────
# LABEL_OFFSET : 标签相对节点中心额外向外偏移的距离（数据坐标单位）
# 字体大小可通过 fontsize 调整；ha/va 根据节点所在象限自动选择
# ─────────────────────────────────────────────────────────────────────────────
LABEL_OFFSET = 0.55

for gene, angle in zip(TARGETS, angles):
    px, py = pos[gene]
    lx = px + LABEL_OFFSET * np.cos(angle)
    ly = py + LABEL_OFFSET * np.sin(angle)

    # 根据角度选择水平对齐方向
    if np.cos(angle) > 0.15:
        ha = 'left'
    elif np.cos(angle) < -0.15:
        ha = 'right'
    else:
        ha = 'center'

    ax.text(lx, ly, gene,
            fontsize=8.5, color='#2B2D42', ha=ha, va='center',
            fontweight='semibold',
            path_effects=[pe.withStroke(linewidth=2.0, foreground='white')])

# ── 标题 ─────────────────────────────────────────────────────────────────────
ax.set_title(f'MAFB Regulon  —  RcisTarget Predicted Targets\n'
             f'({n} target genes, MAFB(+) regulon)',
             fontsize=13, fontweight='bold', pad=16,
             color='#2B2D42')

plt.tight_layout()

# ── 保存 ─────────────────────────────────────────────────────────────────────
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'network_MAFB_regulon_radial.{fmt}')
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'Saved: {out_path}')

plt.close()
print('Done.')
# =============================================================================
#  Radial Layout Network  —  MAFB Regulon (RcisTarget)
#  中心节点: MAFB (TF)
#  环绕节点: 23 个 target genes（每个节点独立颜色，标签字号 x2）
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.cm as cm
import networkx as nx

# ── 路径 ──────────────────────────────────────────────────────────────────────
DATA_DIR = translate(r'D:\bulk-download\GSE182434')
OUT_DIR  = DATA_DIR
DPI      = 150

# ── MAFB regulon targets ──────────────────────────────────────────────────────
TF = 'MAFB'
TARGETS = [
    'CEP192', 'CHCHD10', 'NAGK', 'PPRC1', 'BASP1', 'MTR', 'TIMM17A',
    'AAK1', 'MGAT1', 'USP12', 'FAM20A', 'APOL3', 'GSDMD', 'ELMO1',
    'PICALM', 'RASSF4', 'IFIT3', 'NAE1', 'NAP1L4', 'PPP2R5A',
    'MRGBP', 'PPP1R18', 'CTSL',
]
n = len(TARGETS)   # 23

# ── 23 种区分度高的颜色：tab20(20) + Set1 前3色 ───────────────────────────────
# [PARAM] 如需换色板，替换下方两行即可
_c20   = [cm.tab20(i)  for i in range(20)]
_c_ext = [cm.Set1(i)   for i in range(3)]    # 补足到 23
TARGET_COLORS = _c20 + _c_ext                # list of 23 RGBA tuples

# ── 固定颜色（与 volcano 图一致）─────────────────────────────────────────────
COL_TF       = '#2B2D42'
COL_TF_EDGE  = '#FFB703'
COL_EDGE     = '#AAAAAA'
COL_BG       = '#FAFAFA'

# ── 构建图 ────────────────────────────────────────────────────────────────────
G = nx.DiGraph()
G.add_node(TF)
for t in TARGETS:
    G.add_node(t)
    G.add_edge(TF, t)

# ── Radial 布局 ───────────────────────────────────────────────────────────────
# [PARAM] RADIUS 调大到 4.2，为 17pt 标签留出足够空间
RADIUS      = 4.2
START_ANGLE = np.pi / 2

angles = [START_ANGLE - 2 * np.pi * i / n for i in range(n)]
pos = {TF: np.array([0.0, 0.0])}
for gene, angle in zip(TARGETS, angles):
    pos[gene] = np.array([RADIUS * np.cos(angle),
                           RADIUS * np.sin(angle)])

# ── 画布：figsize 放大以容纳更大字号 ─────────────────────────────────────────
# [PARAM] figsize 从 (10,10) 改为 (14,14)
fig, ax = plt.subplots(figsize=(14, 14))
fig.patch.set_facecolor('white')
ax.set_facecolor(COL_BG)
ax.set_aspect('equal')
ax.axis('off')

# ── 绘制边 ────────────────────────────────────────────────────────────────────
nx.draw_networkx_edges(
    G, pos, ax=ax,
    edgelist=[(TF, t) for t in TARGETS],
    edge_color=COL_EDGE,
    arrows=True,
    arrowstyle='-|>',
    arrowsize=14,
    width=1.2,
    alpha=0.7,
    node_size=[2800, *([900] * n)],
    connectionstyle='arc3,rad=0.0',
    min_source_margin=28,
    min_target_margin=18,
)

# ── 绘制 target 节点（每个节点独立颜色）──────────────────────────────────────
for i, gene in enumerate(TARGETS):
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        nodelist=[gene],
        node_color=[TARGET_COLORS[i]],   # 单独传入该基因的颜色
        node_size=900,
        alpha=0.90,
        linewidths=1.2,
        edgecolors='white',
    )

# ── 绘制 TF 中心节点 ──────────────────────────────────────────────────────────
nx.draw_networkx_nodes(
    G, pos, ax=ax,
    nodelist=[TF],
    node_color=COL_TF,
    node_size=2800,
    alpha=1.0,
    linewidths=3.0,
    edgecolors=COL_TF_EDGE,
)

# ── TF 标签 ───────────────────────────────────────────────────────────────────
ax.text(0, 0, TF,
        fontsize=16, fontweight='bold', color='white',
        ha='center', va='center', zorder=10,
        path_effects=[pe.withStroke(linewidth=1.5, foreground=COL_TF)])

# ── target 标签（字号翻倍：17pt，偏移量同步增大）────────────────────────────
# [PARAM] LABEL_OFFSET 从 0.55 → 0.85；fontsize 从 8.5 → 17
LABEL_OFFSET = 0.85

for i, (gene, angle) in enumerate(zip(TARGETS, angles)):
    px, py = pos[gene]
    lx = px + LABEL_OFFSET * np.cos(angle)
    ly = py + LABEL_OFFSET * np.sin(angle)

    if np.cos(angle) > 0.15:
        ha = 'left'
    elif np.cos(angle) < -0.15:
        ha = 'right'
    else:
        ha = 'center'

    ax.text(lx, ly, gene,
            fontsize=17,                        # [PARAM] 原 8.5 → 翻倍 17
            color=TARGET_COLORS[i],             # 标签颜色与节点颜色一致
            ha=ha, va='center',
            fontweight='semibold',
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])

# ── 标题 ──────────────────────────────────────────────────────────────────────
ax.set_title(f'MAFB Regulon  —  RcisTarget Predicted Targets\n'
             f'({n} target genes, MAFB(+) regulon)',
             fontsize=13, fontweight='bold', pad=16,
             color='#2B2D42')

plt.tight_layout()

# ── 保存 ──────────────────────────────────────────────────────────────────────
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'network_MAFB_regulon_radial.{fmt}')
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'Saved: {out_path}')
plt.close()
print('Done.')
# =============================================================================
#  Organic Tree Network  —  MAFB Regulon (RcisTarget)
#  主干: MAFB (根节点，底部居中)
#  分支: 23 个 target genes（叶节点，上方扇形分布 + 随机扰动）
#  连线: 三次贝塞尔曲线，控制点随机偏移，形成有机感
# =============================================================================

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.cm as cm
from matplotlib.path import Path
import matplotlib.colors as mcolors

# ── 路径 ──────────────────────────────────────────────────────────────────────
DATA_DIR = translate(r'D:\bulk-download\GSE182434')
OUT_DIR  = DATA_DIR
DPI      = 180

# ── MAFB regulon targets ──────────────────────────────────────────────────────
TF = 'MAFB'
TARGETS = [
    'CEP192', 'CHCHD10', 'NAGK',   'PPRC1',  'BASP1',
    'MTR',    'TIMM17A', 'AAK1',   'MGAT1',  'USP12',
    'FAM20A', 'APOL3',   'GSDMD',  'ELMO1',  'PICALM',
    'RASSF4', 'IFIT3',   'NAE1',   'NAP1L4', 'PPP2R5A',
    'MRGBP',  'PPP1R18', 'CTSL',
]
n = len(TARGETS)   # 23

# ── 颜色 ──────────────────────────────────────────────────────────────────────
_c20 = [cm.tab20(i) for i in range(20)]
_cx  = [cm.Set1(i)  for i in range(3)]
TARGET_COLORS = _c20 + _cx   # 23 RGBA

COL_TF      = '#2B2D42'
COL_TF_RING = '#FFB703'
COL_BG      = '#FAFAFA'

# ── 随机种子（固定，保证可复现）──────────────────────────────────────────────
RNG = np.random.default_rng(seed=8)

# ── 根节点位置 ────────────────────────────────────────────────────────────────
ROOT = np.array([0.0, 0.0])

# ── 叶节点位置：扇形均匀分布 + 独立径向 & 角向扰动 ──────────────────────────
# [PARAM] 叶节点平均半径、扇形角度范围、扰动幅度
RADIUS_MEAN  = 6.5     # 平均到根的距离
RADIUS_JITTER = 2    # 径向随机扰动 ±
ANGLE_SPAN   = 200     # 扇形总角度（度），以正上方为中心
ANGLE_JITTER = 9.0     # 角向随机扰动（度）

base_angles_deg = np.linspace(
    90 - ANGLE_SPAN / 2,
    90 + ANGLE_SPAN / 2,
    n
)
jitter_angles = RNG.uniform(-ANGLE_JITTER, ANGLE_JITTER, n)
jitter_r      = RNG.uniform(-RADIUS_JITTER, RADIUS_JITTER, n)

leaf_pos = {}
for i, gene in enumerate(TARGETS):
    angle_rad = np.deg2rad(base_angles_deg[i] + jitter_angles[i])
    r = RADIUS_MEAN + jitter_r[i]
    leaf_pos[gene] = np.array([r * np.cos(angle_rad),
                                r * np.sin(angle_rad)])

# ── 画布 ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 14))
fig.patch.set_facecolor('white')
ax.set_facecolor(COL_BG)
ax.set_aspect('equal')
ax.axis('off')

# ── 绘制贝塞尔曲线（MAFB → target）──────────────────────────────────────────
def draw_bezier(ax, p0, p3, color, lw=1.8, alpha=0.65):
    """
    三次贝塞尔曲线，两个控制点在 p0→p3 方向上随机偏移。
    p0: 起点 (根节点)
    p3: 终点 (叶节点)
    """
    vec   = p3 - p0
    perp  = np.array([-vec[1], vec[0]])          # 垂直方向
    perp_norm = perp / (np.linalg.norm(perp) + 1e-9)

    # 控制点1：从 p0 出发，沿主方向 30~55%，横向随机偏移
    t1    = RNG.uniform(0.25, 0.45)
    off1  = RNG.uniform(-1.8, 1.8)
    p1    = p0 + t1 * vec + off1 * perp_norm

    # 控制点2：靠近 p3，横向偏移方向相反（形成 S 形弯曲）
    t2    = RNG.uniform(0.55, 0.75)
    off2  = RNG.uniform(-1.8, 1.8)
    p2    = p0 + t2 * vec + off2 * perp_norm

    verts = [p0, p1, p2, p3]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path  = Path(verts, codes)
    patch = mpatches.FancyArrowPatch.__new__(mpatches.FancyArrowPatch)
    # 直接用 PathPatch 画曲线
    pp = mpatches.PathPatch(
        path,
        facecolor='none',
        edgecolor=color,
        linewidth=lw,
        alpha=alpha,
        zorder=1,
        capstyle='round',
    )
    ax.add_patch(pp)

    # 在终点附近加小箭头（用 annotate 实现方向感）
    # 取曲线末端切线方向（p2 → p3）
    tang = p3 - p2
    tang /= np.linalg.norm(tang) + 1e-9
    arrow_start = p3 - tang * 0.45
    ax.annotate(
        '', xy=p3, xytext=arrow_start,
        arrowprops=dict(
            arrowstyle='-|>',
            color=color,
            lw=0.8,
            mutation_scale=10,
        ),
        zorder=2,
    )

for i, gene in enumerate(TARGETS):
    draw_bezier(ax, ROOT, leaf_pos[gene],
                color=TARGET_COLORS[i],
                lw=1.8, alpha=0.70)

# ── 绘制叶节点 ────────────────────────────────────────────────────────────────
for i, gene in enumerate(TARGETS):
    px, py = leaf_pos[gene]
    circle = plt.Circle(
        (px, py), radius=0.38,
        color=TARGET_COLORS[i],
        alpha=0.92, zorder=4,
        linewidth=1.5,
        edgecolor='white',
    )
    ax.add_patch(circle)

# ── 绘制根节点（MAFB）────────────────────────────────────────────────────────
# 外圈光晕
glow = plt.Circle((0, 0), radius=0.85,
                   color=COL_TF_RING, alpha=0.25, zorder=3)
ax.add_patch(glow)
# 主节点
root_circle = plt.Circle((0, 0), radius=0.65,
                           color=COL_TF, alpha=1.0, zorder=5,
                           linewidth=3.0, edgecolor=COL_TF_RING)
ax.add_patch(root_circle)
ax.text(0, 0, TF,
        fontsize=17, fontweight='bold', color='white',
        ha='center', va='center', zorder=6,
        path_effects=[pe.withStroke(linewidth=2, foreground=COL_TF)])

# ── 叶节点标签 ────────────────────────────────────────────────────────────────
# [PARAM] fontsize=17（原 8.5 翻倍），偏移量 0.75
LABEL_OFFSET = 0.75

for i, gene in enumerate(TARGETS):
    px, py = leaf_pos[gene]

    # 标签方向：从根节点指向叶节点
    angle_rad = np.arctan2(py - ROOT[1], px - ROOT[0])
    lx = px + LABEL_OFFSET * np.cos(angle_rad)
    ly = py + LABEL_OFFSET * np.sin(angle_rad)

    # 水平对齐方式
    if np.cos(angle_rad) > 0.2:
        ha = 'left'
    elif np.cos(angle_rad) < -0.2:
        ha = 'right'
    else:
        ha = 'center'

    ax.text(lx, ly, gene,
            fontsize=17,
            color=TARGET_COLORS[i],
            ha=ha, va='center',
            fontweight='semibold',
            zorder=7,
            path_effects=[pe.withStroke(linewidth=2.8, foreground='white')])

# ── 标题 ──────────────────────────────────────────────────────────────────────
ax.set_title(
    f'MAFB Regulon  —  Organic Tree Layout\n'
    f'RcisTarget predicted targets  ({n} genes)',
    fontsize=14, fontweight='bold', pad=18, color='#2B2D42',
)

# ── 自动调整视图范围（留边距）────────────────────────────────────────────────
all_x = [ROOT[0]] + [leaf_pos[g][0] for g in TARGETS]
all_y = [ROOT[1]] + [leaf_pos[g][1] for g in TARGETS]
margin = 2.2
ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
ax.set_ylim(min(all_y) - margin, max(all_y) + margin)

plt.tight_layout()

# ── 保存 ──────────────────────────────────────────────────────────────────────
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'network_MAFB_tree_organic.{fmt}')
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'Saved: {out_path}')
plt.close()
print('Done.')
