# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# =============================================================================
#  导出 liana_res_mal.csv
#  前置条件：
#    - adata_processed.h5ad  在 D:\bulk-download\GSE182434\
#    - malignancy_classification.csv  在同目录
# =============================================================================

import scanpy as sc
import liana
import pandas as pd
import os

DATA_DIR = translate(r'D:\bulk-download\GSE182434')

# ── 1. 读取 adata ─────────────────────────────────────────────────────────────
adata = sc.read_h5ad(translate(os.path.join(DATA_DIR, 'adata_processed.h5ad')))
print(f"adata: {adata.shape}")
print(f"Tissue values: {adata.obs['Tissue'].unique().tolist()}")

adata_dlbcl = adata[adata.obs['Tissue'] == 'DLBCL'].copy()
adata_dlbcl.X = adata_dlbcl.layers['counts'].copy()
print(f"DLBCL subset: {adata_dlbcl.shape}")

# ── 2. 重新打 Malignant / Normal B 标签 ──────────────────────────────────────
mal_df = pd.read_csv(
    os.path.join(DATA_DIR, 'malignancy_classification.csv'), index_col=0
)
print(f"\nMalignancy labels:\n{mal_df['malignancy'].value_counts()}")

adata_mal = adata_dlbcl.copy()
adata_mal.obs['cell_type_mal'] = adata_mal.obs['cell_type'].astype(str)

for idx in adata_mal.obs_names:
    if idx in mal_df.index:
        status = mal_df.loc[idx, 'malignancy']
        if status == 'Malignant B cell':
            adata_mal.obs.loc[idx, 'cell_type_mal'] = 'Malignant B cell'
        elif 'Normal B cell' in str(status):
            adata_mal.obs.loc[idx, 'cell_type_mal'] = 'Normal B cell'

print(f"\nNew cell_type_mal distribution:")
print(adata_mal.obs['cell_type_mal'].value_counts())

# ── 3. 跑 LIANA rank_aggregate ────────────────────────────────────────────────
print("\nRunning LIANA rank_aggregate (Malignant/Normal B labels)...")
liana.mt.rank_aggregate(
    adata_mal,
    groupby='cell_type_mal',
    expr_prop=0.1,
    min_cells=5,
    use_raw=False,
    n_perms=100,
    seed=42,
    verbose=True,
)

liana_res_mal = adata_mal.uns['liana_res'].copy()
print(f"\nDone: {liana_res_mal.shape[0]:,} interactions")
print(f"Sources: {sorted(liana_res_mal['source'].unique())}")
print(f"Targets: {sorted(liana_res_mal['target'].unique())}")

# ── 4. 导出 ──────────────────────────────────────────────────────────────────
out_path = os.path.join(DATA_DIR, 'liana_res_mal.csv')
liana_res_mal.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
# =============================================================================
#  BCELL_LIANA — 本地版
#  所需文件：
#    1. liana_res_all.csv          ← 全细胞类型 LIANA 结果
#    2. liana_res_mal.csv          ← 恶性/正常 B 细胞 LIANA 结果
#    3. interaction_weight_matrix.csv  ← 已有
#    4. malignancy_classification.csv  ← 已有（仅用于参考，绘图不直接用）
#  输出目录：OUT_DIR（自动创建）
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ── 路径配置（按需修改）──────────────────────────────────────────────────────
DATA_DIR = translate(r'D:\bulk-download\GSE182434')   # [PARAM] 数据目录
OUT_DIR  = os.path.join(DATA_DIR, 'cellcomm')
os.makedirs(OUT_DIR, exist_ok=True)

LIANA_ALL_CSV = os.path.join(DATA_DIR, 'liana_res_all.csv')       # [PARAM]
LIANA_MAL_CSV = os.path.join(DATA_DIR, 'liana_res_mal.csv')       # [PARAM]
DPI = 300

# ── 读取数据 ──────────────────────────────────────────────────────────────────
liana_res     = pd.read_csv(LIANA_ALL_CSV)
liana_res_mal = pd.read_csv(LIANA_MAL_CSV)

print(f"liana_res (all):     {liana_res.shape}")
print(f"liana_res (mal/nor): {liana_res_mal.shape}")
print(f"Columns: {liana_res.columns.tolist()}")

# ── 颜色配置 ──────────────────────────────────────────────────────────────────
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
ct_short = {
    'CD8 T cells': 'CD8 T', 'CD4 T cells': 'CD4 T',
    'Tregs': 'Tregs', 'TFH cells': 'TFH', 'NK cells': 'NK',
    'GC B cells': 'GC B', 'Proliferative GC B cells': 'Prolif. GC B',
    'Naive B cells': 'Naive B', 'Memory B cells': 'Memory B',
    'Age-associated B cells': 'Age-assoc. B', 'B cells (other)': 'B (other)',
    'Plasma cells': 'Plasma', 'Monocytes/Macrophages': 'Mono/Mac',
    'pDC/Other': 'pDC/Other',
}

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
#  函数：circle communication plot（最终版 v3）
# =============================================================================
def draw_circle_comm_v3(liana_df, palette, short_labels, title,
                        top_n=200, figsize=(9.8, 9.8), save_path=None):
    df = liana_df.sort_values('magnitude_rank').head(top_n).copy()
    df['weight'] = 1 - df['magnitude_rank']

    pair_df = (df.groupby(['source', 'target'])
               .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
               .reset_index())

    active_cts = [c for c in palette
                  if c in set(pair_df['source']) | set(pair_df['target'])]
    n_ct = len(active_cts)

    node_tot = {ct: (pair_df[pair_df['source'] == ct]['total_weight'].sum() +
                     pair_df[pair_df['target'] == ct]['total_weight'].sum())
                for ct in active_cts}

    angles = {ct: np.pi / 2 - 2 * np.pi * i / n_ct
              for i, ct in enumerate(active_cts)}
    pos = {ct: (np.cos(angles[ct]), np.sin(angles[ct])) for ct in active_cts}

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    # ── 边 ────────────────────────────────────────────────────────────────────
    w_vals = pair_df['total_weight'].values
    w_min, w_max = w_vals.min(), w_vals.max()

    for _, row in pair_df.sort_values('total_weight').iterrows():
        src, tgt = row['source'], row['target']
        if src not in pos or tgt not in pos or src == tgt:
            continue
        w_norm = (row['total_weight'] - w_min) / (w_max - w_min + 1e-10)
        lw    = (0.4 + 5.6 * w_norm) * 3.0
        rad   = 0.2 + 0.15 * w_norm
        color = palette.get(src, '#888888')
        hw    = 0.015 + 0.025 * w_norm
        hl    = 0.025 + 0.030 * w_norm

        ax.annotate(
            '', xy=pos[tgt], xytext=pos[src],
            arrowprops=dict(
                arrowstyle=f'-|>,head_width={hw:.3f},head_length={hl:.3f}',
                connectionstyle=f'arc3,rad={rad:.2f}',
                color=color, lw=lw, alpha=1.0,
            ),
            zorder=2,
        )

    # ── 节点 ──────────────────────────────────────────────────────────────────
    ns_max = max(node_tot.values()) if node_tot else 1
    for ct in active_cts:
        x, y = pos[ct]
        ns_norm = node_tot.get(ct, 0) / ns_max
        r = (0.055 + 0.095 * ns_norm) * 0.75
        circle = plt.Circle((x, y), r,
                             color=palette.get(ct, '#888888'),
                             zorder=5, linewidth=1.8, edgecolor='white')
        ax.add_patch(circle)

    # ── 标签 ──────────────────────────────────────────────────────────────────
    label_r = 1.20
    for ct in active_cts:
        ang = angles[ct]
        deg = np.degrees(ang) % 360
        ha  = 'right' if 80 < deg < 280 else 'left'
        lx  = (label_r + 0.08) * np.cos(ang)
        ly  = (label_r + 0.08) * np.sin(ang)

        ax.text(lx, ly, short_labels.get(ct, ct),
                ha=ha, va='center', fontsize=15, fontweight='bold',
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
    ax.legend(handles=legend_els, loc='lower left', fontsize=12,
              frameon=True, framealpha=0.9,
              bbox_to_anchor=(-0.08, -0.08),
              title='Edge / Node encoding', title_fontsize=12)

    if save_path:
        for fmt in ['png', 'svg']:
            fig.savefig(f'{save_path}.{fmt}', dpi=DPI,
                        bbox_inches='tight', format=fmt)
        plt.close(fig)
        print(f"Saved: {save_path}")
    return fig


# =============================================================================
#  Fig 1：全细胞类型 circle plot
# =============================================================================
draw_circle_comm_v3(
    liana_res, ct_palette, ct_short,
    title='Cell-Cell Communication Network\nDLBCL — LIANA rank_aggregate (Top 200 interactions)',
    top_n=200, figsize=(9.8, 9.8),
    save_path=os.path.join(OUT_DIR, 'circle_plot_all')
)

# =============================================================================
#  Fig 2：Malignant vs Normal B cell circle plot
# =============================================================================
draw_circle_comm_v3(
    liana_res_mal, ct_palette_mal, ct_short_mal,
    title='Cell-Cell Communication Network (Malignant vs Normal B cells)\nDLBCL — LIANA rank_aggregate (Top 200 interactions)',
    top_n=200, figsize=(9.8, 9.8),
    save_path=os.path.join(OUT_DIR, 'Malignant_vs_Normal_B_cell_circle_plot')
)


# =============================================================================
#  Fig 3：Communication strength heatmap（Malignant/Normal B cell）
# =============================================================================
top_n_hm = 300
df_hm = liana_res_mal.sort_values('magnitude_rank').head(top_n_hm).copy()
df_hm['weight'] = 1 - df_hm['magnitude_rank']

pair_agg = (df_hm.groupby(['source', 'target'])
            .agg(total_weight=('weight', 'sum'), n=('weight', 'count'))
            .reset_index())

all_cts = sorted(set(pair_agg['source'].tolist() + pair_agg['target'].tolist()))
order = ['Malignant B cell', 'Normal B cell',
         'CD4 T cells', 'CD8 T cells', 'Tregs', 'TFH cells',
         'NK cells', 'Monocytes/Macrophages', 'pDC/Other']
order = [c for c in order if c in all_cts]

# ── 重命名显示标签 ────────────────────────────────────────────────────────────
label_map = {'Monocytes/Macrophages': 'Mono/Mac'}   # [PARAM] 按需增减
order_display = [label_map.get(c, c) for c in order]

mat = pd.DataFrame(0.0, index=order_display, columns=order_display)
for _, row in pair_agg.iterrows():
    src = label_map.get(row['source'], row['source'])
    tgt = label_map.get(row['target'], row['target'])
    if src in order_display and tgt in order_display:
        mat.loc[src, tgt] = row['total_weight']


fig, ax = plt.subplots(figsize=(9.8, 8.17))
mask = np.eye(len(order), dtype=bool)
sns.heatmap(
    mat, mask=mask,
    cmap='YlOrRd', linewidths=0.5, linecolor='#dddddd',
    annot=True, fmt='.2f', annot_kws={'size': 14},
    ax=ax,
    cbar_kws={'label': 'Aggregated interaction weight\n(sum of 1 − magnitude_rank)',
              'shrink': 0.7},
    square=True,
)
ax.set_title('Cell-Cell Communication Strength\nMalignant vs Normal B cells — LIANA rank_aggregate (Top 300)',
             fontsize=12, fontweight='bold', pad=14)
ax.set_xlabel('Target cell type', fontsize=14, labelpad=8)
ax.set_ylabel('Source cell type', fontsize=14, labelpad=8)
ax.tick_params(axis='x', rotation=40, labelsize=14)
ax.tick_params(axis='y', rotation=0,  labelsize=14)

for tick in ax.get_xticklabels() + ax.get_yticklabels():
    if 'Malignant' in tick.get_text():
        tick.set_color('#d62728'); tick.set_fontweight('bold')
    elif 'Normal B' in tick.get_text():
        tick.set_color('#2ca02c'); tick.set_fontweight('bold')

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig.savefig(os.path.join(OUT_DIR, f'heatmap_malignant_vs_normal.{fmt}'),
                dpi=DPI, bbox_inches='tight', format=fmt)
plt.close(fig)
print("Saved: heatmap_malignant_vs_normal")


# =============================================================================
#  Fig 4：LR pair dotplot（Malignant vs Normal B cell）
# =============================================================================
b_cells = ['Malignant B cell', 'Normal B cell']

df_b = liana_res_mal[
    liana_res_mal['source'].isin(b_cells) | liana_res_mal['target'].isin(b_cells)
].copy()
df_b['lr_pair']     = df_b['ligand_complex'] + ' → ' + df_b['receptor_complex']
df_b['pair_label']  = df_b['source'] + ' → ' + df_b['target']
df_b['score']       = 1 - df_b['magnitude_rank']
df_b['specificity'] = 1 - df_b['specificity_rank']

directions = [
    'Malignant B cell → CD8 T cells', 'Malignant B cell → CD4 T cells',
    'Malignant B cell → Monocytes/Macrophages', 'Malignant B cell → NK cells',
    'Malignant B cell → Tregs', 'Malignant B cell → TFH cells',
    'CD8 T cells → Malignant B cell', 'CD4 T cells → Malignant B cell',
    'Monocytes/Macrophages → Malignant B cell', 'NK cells → Malignant B cell',
    'Normal B cell → CD8 T cells', 'Normal B cell → CD4 T cells',
    'Normal B cell → Monocytes/Macrophages',
    'CD8 T cells → Normal B cell', 'CD4 T cells → Normal B cell',
    'Monocytes/Macrophages → Normal B cell',
]
directions = [d for d in directions if d in df_b['pair_label'].unique()]

df_filt = df_b[df_b['pair_label'].isin(directions)].copy()
lr_mean = df_filt.groupby('lr_pair')['score'].mean().nlargest(30)
top_lr  = lr_mean.index.tolist()
df_filt = df_filt[df_filt['lr_pair'].isin(top_lr)]

score_mat = df_filt.pivot_table(index='lr_pair', columns='pair_label',
                                values='score', aggfunc='mean')
spec_mat  = df_filt.pivot_table(index='lr_pair', columns='pair_label',
                                values='specificity', aggfunc='mean')

score_mat = score_mat.reindex(columns=directions).dropna(how='all')
spec_mat  = spec_mat.reindex(columns=directions).reindex(score_mat.index)

row_order = score_mat.fillna(0).sum(axis=1).sort_values(ascending=True).index
score_mat = score_mat.loc[row_order]
spec_mat  = spec_mat.loc[row_order]

n_rows, n_cols = score_mat.shape
print(f"Dotplot: {n_rows} LR pairs × {n_cols} directions")

vals = score_mat.values.flatten()
vals = vals[~np.isnan(vals)]
vmin = float(np.percentile(vals, 5))
vmax = float(vals.max())
norm = mcolors.PowerNorm(gamma=0.25, vmin=vmin, vmax=vmax)
cmap = plt.cm.plasma

fig, ax = plt.subplots(figsize=(n_cols * 1.1 + 4.0, n_rows * 0.52 + 3.0))

for i, lr in enumerate(score_mat.index):
    for j, direction in enumerate(score_mat.columns):
        sv  = score_mat.loc[lr, direction]
        spv = spec_mat.loc[lr, direction]
        if pd.isna(sv):
            continue
        color = cmap(norm(np.clip(sv, vmin, vmax)))
        size  = (30 + 300 * (spv if not pd.isna(spv) else 0)) * 3.0
        ax.scatter(j, i, s=size, color=color,
                   edgecolors='#333333', linewidths=0.5, zorder=3)

ax.set_xticks(range(n_cols))
ax.set_xticklabels(score_mat.columns, fontsize=9, rotation=40, ha='right')
ax.set_yticks(range(n_rows))
ax.set_yticklabels(score_mat.index, fontsize=9)
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_rows - 0.5)
ax.grid(True, color='#eeeeee', linewidth=0.7, zorder=0)
ax.set_axisbelow(True)

mal_cols = [j for j, d in enumerate(score_mat.columns) if 'Malignant' in d]
nor_cols  = [j for j, d in enumerate(score_mat.columns) if 'Normal B' in d]
if mal_cols and nor_cols:
    sep = (max(mal_cols) + min(nor_cols)) / 2
    ax.axvline(sep, color='#333333', linewidth=1.5, linestyle='--', alpha=0.7)
if mal_cols:
    ax.text(np.mean(mal_cols), n_rows + 0.15, 'Malignant B cell interactions',
            ha='center', va='bottom', fontsize=9.5, color='#d62728', fontweight='bold')
if nor_cols:
    ax.text(np.mean(nor_cols), n_rows + 0.15, 'Normal B cell interactions',
            ha='center', va='bottom', fontsize=9.5, color='#2ca02c', fontweight='bold')

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.01, aspect=20)
cbar.set_label('Interaction score\n(1 − magnitude_rank)', fontsize=9)
tick_scores = np.linspace(vmin, vmax, 7)
cbar.set_ticks([norm(v) for v in tick_scores])
cbar.set_ticklabels([f'{v:.2f}' for v in tick_scores], fontsize=8)

for spec_ex, label in [(0.2, 'Low'), (0.55, 'Medium'), (1.0, 'High')]:
    ax.scatter([], [], s=(30 + 300 * spec_ex) * 3.0,
               color='#aaaaaa', edgecolors='#333333', linewidths=0.5, label=label)
ax.legend(title='Specificity\n(dot size)', fontsize=8.5, title_fontsize=9,
          loc='upper left', bbox_to_anchor=(1.10, 1.0),
          frameon=True, framealpha=0.9)

ax.set_title('Top Ligand–Receptor Pairs: Malignant vs Normal B cell Interactions\nDLBCL — LIANA rank_aggregate',
             fontsize=11, fontweight='bold', pad=22)
ax.set_xlabel('Interaction direction (Source → Target)', fontsize=10, labelpad=8)
ax.set_ylabel('Ligand → Receptor pair', fontsize=10, labelpad=8)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig.savefig(os.path.join(OUT_DIR, f'dotplot_malignant_vs_normal.{fmt}'),
                dpi=DPI, bbox_inches='tight', format=fmt)
plt.close(fig)
print("Saved: dotplot_malignant_vs_normal")

print("\nAll done.")

# =============================================================================
#  B Cell Malignancy UMAP + Composition Figure — 本地版
#  所需文件：
#    - adata_processed.h5ad
#    - malignancy_classification.csv
#  输出：malignancy_bcell_umap_composition.png / .svg
# =============================================================================

import os
import scanpy as sc
import pandas as pd
import numpy as np
import harmonypy as hm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ── 路径配置 ──────────────────────────────────────────────────────────────────
DATA_DIR = translate(r'D:\bulk-download\GSE182434')
OUT_DIR  = os.path.join(DATA_DIR, 'cnv')
os.makedirs(OUT_DIR, exist_ok=True)
DPI = 300

# ── 1. 读取数据 ───────────────────────────────────────────────────────────────
adata = sc.read_h5ad(translate(os.path.join(DATA_DIR, 'adata_processed.h5ad')))
print(f"adata: {adata.shape}")

mal_df = pd.read_csv(
    os.path.join(DATA_DIR, 'malignancy_classification.csv'), index_col=0
)
print(f"Malignancy CSV: {mal_df.shape}")
print(mal_df['malignancy'].value_counts())

# ── 2. 提取 B 细胞，写入 malignancy 标签 ─────────────────────────────────────
b_types = [
    'Naive B cells', 'GC B cells', 'Proliferative GC B cells',
    'Memory B cells', 'Age-associated B cells', 'B cells (other)', 'Plasma cells'
]

adata_b = adata[adata.obs['cell_type'].isin(b_types)].copy()
print(f"B cells: {adata_b.shape}")

# 写入 malignancy（只更新有记录的 cell）
adata_b.obs['malignancy'] = 'Normal B cell'   # 默认值
common_idx = adata_b.obs_names.intersection(mal_df.index)
adata_b.obs.loc[common_idx, 'malignancy'] = mal_df.loc[common_idx, 'malignancy'].values

print(f"\nMalignancy distribution in B cells:")
print(adata_b.obs['malignancy'].value_counts())

# ── 3. Harmony + UMAP（B 细胞专属降维）──────────────────────────────────────
print("\nRunning Harmony + UMAP for B cells...")

adata_b.X = adata_b.layers['log1p_norm'].copy()

# ── 清除残留的 log1p uns 元数据，避免 highly_variable_genes 报错 ──────────────
if 'log1p' in adata_b.uns:
    del adata_b.uns['log1p']

sc.pp.highly_variable_genes(adata_b, n_top_genes=2000,
                             batch_key='batch', subset=True)

sc.tl.pca(adata_b, n_comps=50)

ho = hm.run_harmony(adata_b.obsm['X_pca'], adata_b.obs,
                    'batch', max_iter_harmony=20)
print(f"Z_corr shape: {ho.Z_corr.shape}")

# harmonypy 版本差异：Z_corr 可能是 (n_pcs, n_cells) 或 (n_cells, n_pcs)
# 统一处理：确保最终是 (n_cells, n_pcs)
if ho.Z_corr.shape[0] == adata_b.n_obs:
    adata_b.obsm['X_pca_harmony'] = ho.Z_corr
else:
    adata_b.obsm['X_pca_harmony'] = ho.Z_corr.T

print(f"Harmony embedding: {adata_b.obsm['X_pca_harmony'].shape}")

sc.pp.neighbors(adata_b, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_b, min_dist=0.3)
print(f"UMAP computed: {adata_b.obsm['X_umap'].shape}")

# ── 4. 颜色 / 标签配置 ────────────────────────────────────────────────────────
mal_palette = {
    'Normal B cell':         '#2ca02c',
    'Normal B cell (DLBCL)': '#74b9ff',
    'Malignant B cell':      '#d62728',
}
mal_order  = ['Normal B cell', 'Normal B cell (DLBCL)', 'Malignant B cell']
mal_labels = ['Normal (Tonsil)', 'Normal (DLBCL)', 'Malignant']

celltype_order = [
    'GC B cells', 'Proliferative GC B cells', 'Age-associated B cells',
    'B cells (other)', 'Plasma cells', 'Memory B cells', 'Naive B cells',
]
ct_short = {
    'Naive B cells':            'Naive B',
    'GC B cells':               'GC B',
    'Proliferative GC B cells': 'Prolif. GC B',
    'Memory B cells':           'Memory B',
    'Age-associated B cells':   'Age-assoc. B',
    'B cells (other)':          'B (other)',
    'Plasma cells':             'Plasma',
}
ct_short_inv = {v: k for k, v in ct_short.items()}

# 根据实际数据里的 Patient 值动态确定顺序
actual_patients = sorted(adata_b.obs['Patient'].unique().tolist())
print(f"Patients in data: {actual_patients}")

patient_order = [p for p in ['DLBCL002', 'DLBCL007', 'DLBCL008', 'DLBCL111', 'T2']
                 if p in actual_patients]
# 补上数据里有但预设列表里没有的 patient
for p in actual_patients:
    if p not in patient_order:
        patient_order.append(p)

patient_labels = {p: p for p in patient_order}
if 'T2' in patient_order:
    patient_labels['T2'] = 'T2 (Tonsil)'

# ── 5. 构建 UMAP DataFrame ────────────────────────────────────────────────────
umap_df = pd.DataFrame(
    adata_b.obsm['X_umap'],
    index=adata_b.obs_names,
    columns=['UMAP1', 'UMAP2']
)
umap_df['malignancy'] = adata_b.obs['malignancy'].values
umap_df['cell_type']  = adata_b.obs['cell_type'].values
umap_df['Patient']    = adata_b.obs['Patient'].values

# ── 6. 绘图 ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(21, 7))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.42,
                        width_ratios=[1.1, 1, 1])

# ── Panel 1：B 细胞 UMAP，按 malignancy 着色 ──────────────────────────────────
ax1 = fig.add_subplot(gs[0])

zorder_map = {'Normal B cell': 2, 'Normal B cell (DLBCL)': 3, 'Malignant B cell': 4}
for grp, lbl in zip(mal_order, mal_labels):
    sub = umap_df[umap_df['malignancy'] == grp]
    ax1.scatter(sub['UMAP1'], sub['UMAP2'],
                c=mal_palette[grp], s=7, alpha=0.65,
                linewidths=0, zorder=zorder_map[grp],
                rasterized=True,
                label=f"{lbl} (n={len(sub):,})")

ax1.set_xlabel('UMAP 1', fontsize=12)
ax1.set_ylabel('UMAP 2', fontsize=12)
ax1.set_title('B Cell UMAP — Malignancy Status\n'
              '(Harmony-integrated, B cells only)',
              fontsize=15, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.legend(fontsize=10, frameon=False, loc='lower left',
           handletextpad=0.4, labelspacing=0.35, markerscale=1.8)

# ── Panel 2：堆叠条形图 — 每个 B 细胞亚型的恶性比例 ──────────────────────────
ax2 = fig.add_subplot(gs[1])

ct_counts = (adata_b.obs.groupby(['cell_type', 'malignancy'])
             .size().unstack(fill_value=0))
for col in mal_order:
    if col not in ct_counts.columns:
        ct_counts[col] = 0
ct_counts = ct_counts[mal_order]
ct_pct = ct_counts.div(ct_counts.sum(axis=1), axis=0) * 100

# 按 celltype_order 排序（只保留数据里有的）
present_ct = [c for c in celltype_order if c in ct_pct.index]
ct_pct = ct_pct.reindex(present_ct)
ct_pct.index = [ct_short[i] for i in ct_pct.index]

bar_colors = [mal_palette[m] for m in mal_order]
ct_pct.plot(kind='barh', stacked=True, color=bar_colors,
            ax=ax2, width=0.65, edgecolor='white', linewidth=0.5, legend=False)

for i, (idx, row) in enumerate(ct_pct.iterrows()):
    mal_pct = row['Malignant B cell']
    x_start = row['Normal B cell'] + row['Normal B cell (DLBCL)']
    if mal_pct > 5:
        ax2.text(x_start + mal_pct / 2, i, f'{mal_pct:.0f}%',
                 ha='center', va='center', fontsize=13,
                 fontweight='bold', color='white')
    # 右侧 n 标注
    orig_name = ct_short_inv.get(idx, idx)
    n = ct_counts.loc[orig_name].sum() if orig_name in ct_counts.index else 0
    ax2.text(101, i, f'n={n}', va='center', ha='left',
             fontsize=11, color='#555555')

ax2.set_xlabel('Proportion (%)', fontsize=12)
ax2.set_ylabel('')
ax2.set_title('Malignancy Composition\nper B Cell Subtype',
              fontsize=15, fontweight='bold')
ax2.set_xlim(0, 100)
ax2.tick_params(axis='y', labelsize=15)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# ── Panel 3：堆叠条形图 — 每个 Patient 的恶性比例 ────────────────────────────
ax3 = fig.add_subplot(gs[2])

pt_counts = (adata_b.obs.groupby(['Patient', 'malignancy'])
             .size().unstack(fill_value=0))
for col in mal_order:
    if col not in pt_counts.columns:
        pt_counts[col] = 0
pt_counts = pt_counts[mal_order].reindex(patient_order)
pt_pct = pt_counts.div(pt_counts.sum(axis=1), axis=0) * 100
pt_pct.index = [patient_labels[p] for p in patient_order]

pt_pct.plot(kind='barh', stacked=True, color=bar_colors,
            ax=ax3, width=0.65, edgecolor='white', linewidth=0.5, legend=False)

for i, (idx, row) in enumerate(pt_pct.iterrows()):
    mal_pct   = row['Malignant B cell']
    x_start   = row['Normal B cell'] + row['Normal B cell (DLBCL)']
    ndlbcl_pct = row['Normal B cell (DLBCL)']

    if mal_pct > 3:
        ax3.text(x_start + mal_pct / 2, i, f'{mal_pct:.0f}%',
                 ha='center', va='center', fontsize=13,
                 fontweight='bold', color='white')
    if ndlbcl_pct > 5:
        ax3.text(row['Normal B cell'] + ndlbcl_pct / 2, i,
                 f'{ndlbcl_pct:.0f}%',
                 ha='center', va='center', fontsize=13, color='#333333')

    # 右侧 n 标注
    pt_name = patient_order[i]
    n = pt_counts.loc[pt_name].sum()
    ax3.text(101, i, f'n={n}', va='center', ha='left',
             fontsize=11, color='#555555')

ax3.set_xlabel('Proportion (%)', fontsize=12)
ax3.set_ylabel('')
ax3.set_title('Malignancy Composition\nper Patient',
              fontsize=15, fontweight='bold')
ax3.set_xlim(0, 100)
ax3.tick_params(axis='y', labelsize=15)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# ── 共享图例 ──────────────────────────────────────────────────────────────────
handles = [mpatches.Patch(color=mal_palette[m], label=l)
           for m, l in zip(mal_order, mal_labels)]
title_handle = mpatches.Patch(    
    color='none',    
    label='Cell Group:'    # 这就是你的 title 文字
    )
fig.legend(handles=[title_handle] + handles, loc='lower center', ncol=4, fontsize=14,
           frameon=False, bbox_to_anchor=(0.5, -0.04),
           title= None)

fig.suptitle(
    'Malignant vs Normal B Cells: UMAP Distribution and Composition\n'
    'GSE182434 DLBCL scRNA-seq  |  inferCNVpy (Tonsil reference)',
    fontsize=15, fontweight='bold', y=1.02
)

# ── 保存 ─────────────────────────────────────────────────────────────────────
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'malignancy_bcell_umap_composition.{fmt}')
    fig.savefig(out_path, dpi=DPI, bbox_inches='tight', format=fmt)
    print(f"Saved: {out_path}")

plt.close(fig)
print("Done.")
