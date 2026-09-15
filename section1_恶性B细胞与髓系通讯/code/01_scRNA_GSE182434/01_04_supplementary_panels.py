# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import os, re, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy.stats import zscore, spearmanr, mannwhitneyu
from scipy.ndimage import gaussian_filter1d
from statsmodels.stats.multitest import multipletests
from collections import Counter
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns
import scanpy as sc
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── 路径配置 ──────────────────────────────────────────────────────────────────
BASE_DIR     = translate(r'D:/bulk-download')
DATA_DIR     = os.path.join(BASE_DIR, 'GSE182434')
CELLCOMM_DIR = os.path.join(DATA_DIR, 'cellcomm')
TRAJ_DIR     = os.path.join(DATA_DIR, 'trajectory')
DEG_DIR      = os.path.join(DATA_DIR, 'IFN_TAM_vs_LA_TAM_DEG')
SCENIC_DIR   = os.path.join(DATA_DIR, 'scenic')
H5AD_PATH    = os.path.join(DATA_DIR, 'adata_processed.h5ad')
DPI          = 300

for d in [TRAJ_DIR, DEG_DIR]:
    os.makedirs(d, exist_ok=True)

# ── 颜色常量 ──────────────────────────────────────────────────────────────────
SUBTYPE_PALETTE = {
    'Mono': '#3498DB', 'DC_1': '#E74C3C', 'LA_TAM': '#2ECC71',
    'IFN_TAM': '#E67E22', 'DC_2': '#9B59B6',
}
SUBTYPES_ORDER    = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
PATH_COLORS       = {'Mono': '#E9C46A', 'IFN_TAM': '#E63946', 'LA_TAM': '#F4A261'}
SUBTYPE_COLORS_SC = {'Mono': '#4878CF', 'IFN_TAM': '#E8601C', 'LA_TAM': '#7BAE7F'}
TF_PALETTE = ['#E8601C', '#F6C141', '#1965B0', '#7BAFDE', '#4EB265',
              '#90C987', '#CAE0AB', '#DC050C', '#882E72', '#B178A6']

MARKER_GENES = {
    'Mono':    ['FCN1', 'S100A9', 'S100A8', 'S100A4', 'APOBEC3A'],
    'DC_1':    ['LTB', 'CLEC10A', 'CD1C', 'JAML', 'CD1E'],
    'LA_TAM':  ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD'],
    'IFN_TAM': ['MT1H', 'MT1G', 'CCL8', 'CCL2', 'MT1X'],
    'DC_2':    ['DNASE1L3', 'RGCC', 'CST3', 'CLEC9A', 'SNX3'],
}

# ── 工具函数 ──────────────────────────────────────────────────────────────────
def save_fig(fig, fname, out_dir=TRAJ_DIR):
    for fmt in ['png', 'svg']:
        fig.savefig(os.path.join(out_dir, f'{fname}.{fmt}'),
                    dpi=DPI, bbox_inches='tight', format=fmt)
    plt.close(fig)
    print(f"  ✓ {fname} saved")

def need_redraw(fname, out_dir=TRAJ_DIR):
    png_ok = os.path.exists(os.path.join(out_dir, f'{fname}.png'))
    svg_ok = os.path.exists(os.path.join(out_dir, f'{fname}.svg'))
    if png_ok and svg_ok:
        print(f"  SKIP {fname} (already exists)")
        return False
    return True

# =============================================================================
#  STEP 1 → 图 2  完整代码
#  100% 按照 trace.txt 原始逻辑还原
# =============================================================================
import os, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scanpy as sc
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── 路径配置 ──────────────────────────────────────────────────────────────────
BASE_DIR  = translate(r'D:/bulk-download')
DATA_DIR  = os.path.join(BASE_DIR, 'GSE182434')
TRAJ_DIR  = os.path.join(DATA_DIR, 'trajectory')
H5AD_PATH = translate(os.path.join(DATA_DIR, 'adata_processed.h5ad'))
DPI       = 150
os.makedirs(TRAJ_DIR, exist_ok=True)

SUBTYPE_PALETTE = {
    'Mono':    '#3498DB',
    'DC_1':    '#E74C3C',
    'LA_TAM':  '#2ECC71',
    'IFN_TAM': '#E67E22',
    'DC_2':    '#9B59B6',
}
SUBTYPES_ORDER = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']

# =============================================================================
print("=" * 60)
print("[STEP 1] Loading h5ad...")
print("=" * 60)

adata = sc.read_h5ad(H5AD_PATH)
adata_mac = adata[adata.obs['cell_type'] == 'Monocytes/Macrophages'].copy()
print(f"Mono/Mac cells: {adata_mac.shape[0]}")

# ── 亚型评分（与 trace.txt Cell 733 完全一致）────────────────────────────────
subtype_markers = {
    'Mono':    ['FCN1', 'S100A9', 'S100A8', 'S100A4', 'APOBEC3A'],
    'DC_1':    ['LTB', 'CLEC10A', 'CD1C', 'JAML', 'CD1E'],
    'LA_TAM':  ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD'],
    'IFN_TAM': ['MT1H', 'MT1G', 'CCL8', 'CCL2', 'MT1X'],
    'DC_2':    ['DNASE1L3', 'RGCC', 'CST3', 'CLEC9A', 'SNX3'],
}
subtypes = list(subtype_markers.keys())

for ct, genes in subtype_markers.items():
    present_genes = [g for g in genes if g in adata_mac.var_names]
    sc.tl.score_genes(adata_mac, gene_list=present_genes,
                      score_name=f'score_{ct}', use_raw=False)

score_cols = [f'score_{ct}' for ct in subtypes]
score_mat  = adata_mac.obs[score_cols].values
best_idx   = np.argmax(score_mat, axis=1)
adata_mac.obs['mac_subtype'] = pd.Categorical(
    [subtypes[i] for i in best_idx], categories=subtypes)

print("\nSubtype distribution:")
print(adata_mac.obs['mac_subtype'].value_counts())

# ── 关键修复：显式赋值 .X，清除 log1p 元数据 ─────────────────────────────────
# trace.txt 原始环境 .X 已是 log-norm，本地 subset 后需要显式赋值
_layer = adata_mac.layers['log1p_norm']
adata_mac.X = (_layer.tocsr().copy() if sp.issparse(_layer)
               else sp.csr_matrix(_layer.copy()))
print(f"  .X assigned: shape={adata_mac.X.shape}, max={adata_mac.X.max():.2f}")

if 'log1p' in adata_mac.uns:
    del adata_mac.uns['log1p']

# ── HVG + scale + PCA（与 trace.txt 完全一致）────────────────────────────────
sc.pp.highly_variable_genes(adata_mac, n_top_genes=2000, flavor='seurat')
print(f"HVGs: {adata_mac.var['highly_variable'].sum()}")

sc.pp.scale(adata_mac, max_value=10)

sc.tl.pca(adata_mac, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata_mac, n_neighbors=15, n_pcs=20)
sc.tl.diffmap(adata_mac, n_comps=15)
print("Diffusion map computed")

# ── 标准 UMAP（trace.txt 有这步，用于对比）───────────────────────────────────
sc.tl.umap(adata_mac, min_dist=0.3)
print("Standard UMAP computed")

# ── Leiden ────────────────────────────────────────────────────────────────────
sc.tl.leiden(adata_mac, resolution=0.5, key_added='leiden_mac')
print(f"Leiden clusters: {adata_mac.obs['leiden_mac'].nunique()}")

# ── PAGA ─────────────────────────────────────────────────────────────────────
# 强制 categories 与 SUBTYPES_ORDER 对齐，保证 PAGA 矩阵行列顺序正确
adata_mac.obs['mac_subtype'] = pd.Categorical(
    adata_mac.obs['mac_subtype'],
    categories=SUBTYPES_ORDER, ordered=False)

sc.tl.paga(adata_mac, groups='mac_subtype')
print("PAGA connectivity matrix:")
conn = adata_mac.uns['paga']['connectivities'].toarray()
conn_df = pd.DataFrame(conn, index=SUBTYPES_ORDER, columns=SUBTYPES_ORDER)
print(conn_df.round(3))

# ── DPT（完全按 trace.txt 原始写法：np.where → 直接整数位置）────────────────
mono_mask     = adata_mac.obs['mac_subtype'] == 'Mono'
mono_idx      = np.where(mono_mask)[0]                      # 整数位置数组
dc1_scores    = adata_mac.obs['score_Mono'].values
root_cell_idx = mono_idx[np.argmax(dc1_scores[mono_idx])]  # 直接整数位置
print(f"\nRoot cell index: {root_cell_idx} (Mono, highest score_Mono)")
adata_mac.uns['iroot'] = root_cell_idx                      # 原始写法
sc.tl.dpt(adata_mac, n_dcs=10)
print(f"DPT pseudotime range: "
      f"{adata_mac.obs['dpt_pseudotime'].min():.3f} – "
      f"{adata_mac.obs['dpt_pseudotime'].max():.3f}")
print("\nMean pseudotime per subtype:")
# 跑完 sc.tl.dpt 之后执行
print("neighbors params:", adata_mac.uns['neighbors']['params'])
print("\nMean pseudotime per subtype:")
print(adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean().sort_values().round(3))

pt_means = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean().sort_values()
print(pt_means.round(3))

# ── PAGA-initialized UMAP（与 trace.txt 完全一致）────────────────────────────
with plt.ioff():
    sc.pl.paga(adata_mac, show=False)   # 触发 paga pos 写入

# 若 pos 缺失则用 networkx 补算
if 'pos' not in adata_mac.uns.get('paga', {}):
    import networkx as nx
    G = nx.from_numpy_array(
        adata_mac.uns['paga']['connectivities'].toarray())
    pos_dict = nx.spring_layout(G, seed=42)
    adata_mac.uns['paga']['pos'] = np.array(
        [pos_dict[i] for i in range(len(SUBTYPES_ORDER))])
    print("  paga pos computed via networkx")

sc.tl.umap(adata_mac, init_pos='paga', min_dist=0.5, random_state=35)
adata_mac.obsm['X_umap_paga'] = adata_mac.obsm['X_umap'].copy()
print("PAGA-initialized UMAP computed")

# ── 路径子集 ──────────────────────────────────────────────────────────────────
path_mask  = adata_mac.obs['mac_subtype'].isin(['Mono', 'IFN_TAM', 'LA_TAM'])
adata_path = adata_mac[path_mask].copy()
adata_path = adata_path[adata_path.obs['dpt_pseudotime'].argsort()].copy()
pt_sorted      = adata_path.obs['dpt_pseudotime'].values
subtype_sorted = adata_path.obs['mac_subtype'].values
print(f"Path cells (Mono+IFN_TAM+LA_TAM): {adata_path.n_obs}")

# ── 共用绘图变量 ──────────────────────────────────────────────────────────────
umap_coords  = adata_mac.obsm['X_umap_paga']
subtypes_arr = adata_mac.obs['mac_subtype'].values
pt_arr       = adata_mac.obs['dpt_pseudotime'].values
conn         = adata_mac.uns['paga']['connectivities'].toarray()
paga_pos     = adata_mac.uns['paga']['pos']
pt_means     = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()
centroids    = {st: umap_coords[subtypes_arr == st].mean(axis=0)
                for st in SUBTYPES_ORDER}
n_cells      = {st: int((subtypes_arr == st).sum()) for st in SUBTYPES_ORDER}
node_sizes_scaled = (np.array([n_cells[st] for st in SUBTYPES_ORDER]) /
                     max(n_cells.values()) * 1200 + 300)

# 图 1 Panel C：固定轨迹箭头（trace.txt 原始硬编码）
trajectory_edges = [
    ('Mono',    'IFN_TAM'),
    ('IFN_TAM', 'LA_TAM'),
    ('Mono',    'DC_2'),
    ('DC_2',    'DC_1'),
]

# 图 2：动态有向边（trace.txt 原始：conn > 0.25 + PT 方向）
directed_edges = []
for i, st_i in enumerate(SUBTYPES_ORDER):
    for j, st_j in enumerate(SUBTYPES_ORDER):
        if j <= i: continue
        w = conn[i, j]
        if w < 0.25: continue
        if pt_means[st_i] < pt_means[st_j]:
            directed_edges.append((st_i, st_j, w))
        else:
            directed_edges.append((st_j, st_i, w))

print("\nDirected edges (图2):")
for src, tgt, w in sorted(directed_edges, key=lambda x: -x[2]):
    print(f"  {src} → {tgt}  (w={w:.3f}, ΔPT={pt_means[tgt]-pt_means[src]:.3f})")

print("STEP 1 complete.")
# =============================================================================
#  Supplementary Figures A–D  完整独立脚本
#  运行前提：已执行主脚本 STEP 1，以下变量已存在于内存中：
#    adata_mac, adata_path, pt_arr, subtypes_arr, pt_sorted, subtype_sorted
#    conn, SUBTYPE_PALETTE, SUBTYPES_ORDER, TRAJ_DIR, DPI
#    save_fig, need_redraw
#  若单独运行本脚本，请先取消注释 "── 独立运行时加载数据 ──" 段落
# =============================================================================

import os, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
warnings.filterwarnings('ignore')

# =============================================================================
#  ── 独立运行时加载数据（若已在主脚本环境中运行，注释掉此段）──────────────────
# =============================================================================
# BASE_DIR  = translate(r'D:/bulk-download')
# DATA_DIR  = os.path.join(BASE_DIR, 'GSE182434')
# TRAJ_DIR  = os.path.join(DATA_DIR, 'trajectory')
# CELLCOMM_DIR = os.path.join(DATA_DIR, 'cellcomm')
# H5AD_PATH = translate(os.path.join(DATA_DIR, 'adata_processed.h5ad'))
# DPI = 300
# os.makedirs(TRAJ_DIR, exist_ok=True)
#
# SUBTYPE_PALETTE = {
#     'Mono': '#3498DB', 'DC_1': '#E74C3C', 'LA_TAM': '#2ECC71',
#     'IFN_TAM': '#E67E22', 'DC_2': '#9B59B6',
# }
# SUBTYPES_ORDER = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
#
# # 重新加载 adata_mac（需已完成主脚本 STEP1 的预处理）
# adata = sc.read_h5ad(H5AD_PATH)
# # ... 此处省略，建议直接在主脚本环境中运行本脚本
#
# def save_fig(fig, fname, out_dir=TRAJ_DIR):
#     for fmt in ['png', 'svg']:
#         fig.savefig(os.path.join(out_dir, f'{fname}.{fmt}'),
#                     dpi=DPI, bbox_inches='tight', format=fmt)
#     plt.close(fig)
#     print(f"  ✓ {fname} saved")
#
# def need_redraw(fname, out_dir=TRAJ_DIR):
#     png_ok = os.path.exists(os.path.join(out_dir, f'{fname}.png'))
#     svg_ok = os.path.exists(os.path.join(out_dir, f'{fname}.svg'))
#     if png_ok and svg_ok:
#         print(f"  SKIP {fname} (already exists)")
#         return False
#     return True

# =============================================================================
#  Supplementary A：各亚群 DPT Pseudotime 分布
#  violin + boxplot + 均值点 + 相邻对显著性检验
# =============================================================================
if need_redraw('supp_A_pseudotime_distribution'):
    print("\n[Supp A] Pseudotime distribution by subtype...")

    # 按伪时间均值升序排列亚群
    pt_means_all = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()
    SUPP_ORDER_A = pt_means_all.sort_values().index.tolist()
    print(f"  Subtype order (by mean PT): {SUPP_ORDER_A}")

    pt_data   = [adata_mac.obs.loc[adata_mac.obs['mac_subtype'] == st,
                                   'dpt_pseudotime'].values
                 for st in SUPP_ORDER_A]
    n_cells_s = [len(d) for d in pt_data]
    colors_s  = [SUBTYPE_PALETTE[st] for st in SUPP_ORDER_A]

    # ── 相邻对 Mann-Whitney U 检验（BH 校正）─────────────────────────────────
    pairs    = list(zip(SUPP_ORDER_A[:-1], SUPP_ORDER_A[1:]))
    pvals_mw = []
    for st1, st2 in pairs:
        d1 = adata_mac.obs.loc[adata_mac.obs['mac_subtype'] == st1,
                                'dpt_pseudotime'].values
        d2 = adata_mac.obs.loc[adata_mac.obs['mac_subtype'] == st2,
                                'dpt_pseudotime'].values
        _, p = mannwhitneyu(d1, d2, alternative='two-sided')
        pvals_mw.append(p)
    _, padj_mw, _, _ = multipletests(pvals_mw, method='fdr_bh')

    def pval_star(p):
        if p < 0.001: return '***'
        if p < 0.01:  return '**'
        if p < 0.05:  return '*'
        return 'ns'

    # ── 绘图 ─────────────────────────────────────────────────────────────────
    fig_sa, ax_sa = plt.subplots(figsize=(11, 6.5))
    fig_sa.patch.set_facecolor('white')
    ax_sa.set_facecolor('white')

    # violin
    parts = ax_sa.violinplot(
        pt_data, positions=range(len(SUPP_ORDER_A)),
        showmedians=False, showextrema=False, widths=0.72)
    for pc, col in zip(parts['bodies'], colors_s):
        pc.set_facecolor(col)
        pc.set_alpha(0.50)
        pc.set_edgecolor('white')
        pc.set_linewidth(0)

    # boxplot（叠加）
    bp = ax_sa.boxplot(
        pt_data, positions=range(len(SUPP_ORDER_A)),
        widths=0.22, patch_artist=True,
        medianprops=dict(color='white', linewidth=2.5),
        whiskerprops=dict(color='#666666', linewidth=1.2),
        capprops=dict(color='#666666', linewidth=1.2),
        flierprops=dict(marker='o', markersize=2.0,
                        markerfacecolor='#AAAAAA',
                        markeredgewidth=0, alpha=0.35),
        zorder=3)
    for patch, col in zip(bp['boxes'], colors_s):
        patch.set_facecolor(col)
        patch.set_alpha(0.90)
        patch.set_edgecolor('white')
        patch.set_linewidth(1.0)

    # 均值白点
    for i, d in enumerate(pt_data):
        ax_sa.scatter(i, np.mean(d), color='white', s=50, zorder=6,
                      edgecolors='#333333', linewidths=1.5)

    # 各细胞散点（jitter，稀疏采样避免过密）
    rng = np.random.default_rng(42)
    for i, (d, col) in enumerate(zip(pt_data, colors_s)):
        n_show = min(len(d), 300)
        idx    = rng.choice(len(d), n_show, replace=False)
        jitter = rng.uniform(-0.12, 0.12, n_show)
        ax_sa.scatter(i + jitter, d[idx],
                      color=col, s=5, alpha=0.25,
                      linewidths=0, zorder=2)

    # 相邻对显著性括号
    y_top  = max(d.max() for d in pt_data)
    y_step = (y_top - ax_sa.get_ylim()[0]) * 0.055
    for k, ((st1, st2), p_adj) in enumerate(zip(pairs, padj_mw)):
        i1 = SUPP_ORDER_A.index(st1)
        i2 = SUPP_ORDER_A.index(st2)
        y_line = y_top + y_step * (k + 1.2)
        ax_sa.plot([i1, i1, i2, i2],
                   [y_line - y_step * 0.15, y_line,
                    y_line, y_line - y_step * 0.15],
                   color='#555555', lw=1.1)
        star = pval_star(p_adj)
        ax_sa.text((i1 + i2) / 2, y_line + y_step * 0.05,
                   star, ha='center', va='bottom',
                   fontsize=11 if star != 'ns' else 9,
                   color='#333333',
                   fontweight='bold' if star != 'ns' else 'normal')

    # 坐标轴
    ax_sa.set_xticks(range(len(SUPP_ORDER_A)))
    ax_sa.set_xticklabels(
        [f'{st}\n(n={n:,})' for st, n in zip(SUPP_ORDER_A, n_cells_s)],
        fontsize=11)
    for tick, col in zip(ax_sa.get_xticklabels(), colors_s):
        tick.set_color(col)
        tick.set_fontweight('bold')

    ax_sa.set_ylabel('DPT Pseudotime', fontsize=12)
    ax_sa.set_xlim(-0.6, len(SUPP_ORDER_A) - 0.4)
    ax_sa.spines[['top', 'right']].set_visible(False)

    ax_sa.set_title(
        'Supplementary A — DPT Pseudotime Distribution by Subtype\n'
        '(violin + box; white dot = mean; brackets = Mann-Whitney U, BH-corrected)',
        fontsize=12, fontweight='bold', pad=10)

    # 图例
    legend_sa = [mpatches.Patch(color=SUBTYPE_PALETTE[st], label=st)
                 for st in SUPP_ORDER_A]
    ax_sa.legend(handles=legend_sa, fontsize=9, loc='upper left',
                 framealpha=0.9, edgecolor='none',
                 title='Subtype', title_fontsize=9.5)

    plt.tight_layout()
    save_fig(fig_sa, 'supp_A_pseudotime_distribution')


# =============================================================================
#  Supplementary B：PAGA Connectivity Matrix Heatmap
#  配色参考 liana - 副本.txt：YlOrRd + 白色网格线 + sns.heatmap 风格
# =============================================================================

if need_redraw('supp_B_paga_connectivity_heatmap'):
    print("\n[Supp B] PAGA connectivity matrix heatmap...")

    import seaborn as sns
    import matplotlib.cm as mcm
    import matplotlib.colorbar as mcolorbar

    conn_full = adata_mac.uns['paga']['connectivities'].toarray()
    conn_df_b = pd.DataFrame(conn_full,
                              index=SUBTYPES_ORDER,
                              columns=SUBTYPES_ORDER)
    # 对称化，对角线设 NaN
    conn_sym = (conn_df_b + conn_df_b.T) / 2
    conn_plot = conn_sym.copy()
    np.fill_diagonal(conn_plot.values, np.nan)

    # ── 图布局 ────────────────────────────────────────────────────────────────
    fig_sb, ax_sb = plt.subplots(figsize=(7.5, 6.5))
    fig_sb.patch.set_facecolor('white')

    # 用于 sns.heatmap 的掩码（对角线不显示）
    mask_diag = np.eye(len(SUBTYPES_ORDER), dtype=bool)

    # ── sns.heatmap（与 liana 图1 完全一致的风格）────────────────────────────
    sns.heatmap(
        conn_plot,
        ax=ax_sb,
        cmap='YlOrRd',           # ← 与 liana interaction_heatmap 一致
        linewidths=0.4,           # ← 格子间分割线宽度
        linecolor='white',        # ← 分割线颜色
        mask=mask_diag,           # ← 对角线不填色
        annot=False,              # 先不用 annot，手动写数值（格式更可控）
        vmin=0,
        vmax=conn_plot.values[~np.isnan(conn_plot.values)].max(),
        cbar_kws={
            'label': 'PAGA Connectivity Weight',
            'shrink': 0.7,
        },
    )

    # ── 手动数值标注（对角线显示"—"，其余显示3位小数）────────────────────────
    n_st   = len(SUBTYPES_ORDER)
    vmax_b = conn_plot.values[~np.isnan(conn_plot.values)].max()

    for i in range(n_st):
        for j in range(n_st):
            if i == j:
                ax_sb.text(j + 0.5, i + 0.5, '—',
                           ha='center', va='center',
                           fontsize=11, color='#BBBBBB')
            else:
                val = conn_plot.values[i, j]
                if np.isnan(val):
                    continue
                text_col = 'white' if val > vmax_b * 0.60 else '#333333'
                ax_sb.text(j + 0.5, i + 0.5, f'{val:.3f}',
                           ha='center', va='center',
                           fontsize=11, color=text_col,
                           fontweight='bold')

    # ── 轴标签（与 liana 图1 风格一致：旋转 + 带颜色）────────────────────────
    ax_sb.set_xticklabels(
        SUBTYPES_ORDER,
        fontsize=11.5, rotation=45, ha='right')
    ax_sb.set_yticklabels(
        SUBTYPES_ORDER,
        fontsize=11.5, rotation=0)

    # 用各亚群颜色标注轴标签（与主图一致）
    for tick, st in zip(ax_sb.get_xticklabels(), SUBTYPES_ORDER):
        tick.set_color(SUBTYPE_PALETTE[st])
        tick.set_fontweight('bold')
    for tick, st in zip(ax_sb.get_yticklabels(), SUBTYPES_ORDER):
        tick.set_color(SUBTYPE_PALETTE[st])
        tick.set_fontweight('bold')

    # ── colorbar 字号（与 liana 图1 一致）───────────────────────────────────
    cbar_sb = ax_sb.collections[0].colorbar
    cbar_sb.ax.tick_params(labelsize=9)
    cbar_sb.set_label('PAGA Connectivity Weight', fontsize=10)

    ax_sb.set_xlabel('Target Subtype', fontsize=11, labelpad=8)
    ax_sb.set_ylabel('Source Subtype', fontsize=11, labelpad=8)
    ax_sb.set_title(
        'Supplementary B — PAGA Connectivity Matrix\n'
        '(edge weight between subtypes; "—" = self-connection, omitted)',
        fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()
    save_fig(fig_sb, 'supp_B_paga_connectivity_heatmap')
# =============================================================================
#  Supplementary C：各亚群组织来源 / 患者来源组成
#  配色参考 PAOTU.TXT：
#    sample_id  → patient_colors（DLBCL002/007/008/111 + Tonsil）
#    tissue_type → tissue_colors（DLBCL=#0279EE, Tonsil/Normal=#FF9400）
#  stacked barplot 风格与 PAOTU 图4 完全一致
# =============================================================================

if need_redraw('supp_C_sample_composition'):
    print("\n[Supp C] Sample / tissue composition by subtype...")

    # ── 与 PAOTU.TXT 完全一致的配色 ──────────────────────────────────────────
    # 患者/样本颜色（从 barcode 后缀推断）
    patient_colors_sc = {
        'DLBCL002': '#1f77b4',
        'DLBCL007': '#ff7f0e',
        'DLBCL008': '#2ca02c',
        'DLBCL111': '#d62728',
        'T2':       '#9467bd',   # Tonsil 样本
        'Unknown':  '#7f7f7f',
    }

    # 组织类型颜色（与 PAOTU tissue_colors 完全一致）
    tissue_colors_sc = {
        'DLBCL':  '#0279EE',
        'Tonsil': '#FF9400',
        'Normal': '#FF9400',    # 部分数据集用 Normal 代替 Tonsil
    }

    obs_cols = adata_mac.obs.columns.tolist()

    # ── 自动探测列名 ──────────────────────────────────────────────────────────
    sample_col = next(
        (c for c in ['sample_id', 'patient_id', 'sample', 'donor_id',
                      'orig.ident', 'orig_ident', 'batch', 'Sample', 'Patient']
         if c in obs_cols), None)
    tissue_col = next(
        (c for c in ['tissue_type', 'tissue', 'disease', 'condition',
                      'group', 'sample_type', 'Tissue']
         if c in obs_cols), None)

    # ── 若 obs 里没有现成列，从 barcode 后缀自动生成 ─────────────────────────
    if sample_col is None:
        print("  sample_col not found, extracting from barcode suffix...")
        adata_mac.obs['_sample_id'] = (
            adata_mac.obs_names.str.extract(r'_([A-Za-z0-9]+)$')[0]
            .fillna('Unknown'))
        sample_col = '_sample_id'

    if tissue_col is None:
        print("  tissue_col not found, inferring from sample_id...")
        def infer_tissue(s):
            s = str(s)
            if 'NB' in s or 'T2' in s or 'tonsil' in s.lower():
                return 'Tonsil'
            return 'DLBCL'
        adata_mac.obs['_tissue_type'] = (
            adata_mac.obs[sample_col].apply(infer_tissue))
        tissue_col = '_tissue_type'

    print(f"  sample_col = {sample_col}")
    print(f"  tissue_col = {tissue_col}")
    print(f"  Sample values: "
          f"{adata_mac.obs[sample_col].value_counts().to_dict()}")
    print(f"  Tissue values: "
          f"{adata_mac.obs[tissue_col].value_counts().to_dict()}")

    # ── 准备数据 ──────────────────────────────────────────────────────────────
    obs_sc = adata_mac.obs[['mac_subtype', sample_col, tissue_col]].copy()

    # ── 图布局：1×2 ───────────────────────────────────────────────────────────
    fig_sc, axes_sc = plt.subplots(1, 2, figsize=(18, 6.5))
    fig_sc.patch.set_facecolor('white')

    # ════════════════════════════════════════════════════════════════════════
    #  Panel 1：Stacked barplot by sample_id（参考 PAOTU 图4 风格）
    # ════════════════════════════════════════════════════════════════════════
    ax_sc1 = axes_sc[0]
    ax_sc1.set_facecolor('white')

    # 计算各亚群中每个 sample 的细胞比例
    ct_s = (obs_sc.groupby(['mac_subtype', sample_col])
            .size().unstack(fill_value=0))
    ct_s_pct = ct_s.div(ct_s.sum(axis=1), axis=0) * 100

    # 保证行顺序与 SUBTYPES_ORDER 一致
    for st in SUBTYPES_ORDER:
        if st not in ct_s_pct.index:
            ct_s_pct.loc[st] = 0.0
    ct_s_pct = ct_s_pct.loc[SUBTYPES_ORDER]

    # 样本配色：优先从 patient_colors_sc 匹配，否则用 tab20
    sample_list = ct_s_pct.columns.tolist()
    samp_color_map = {}
    used_tab20 = []
    tab20_cmap = plt.cm.get_cmap('tab20', max(len(sample_list), 2))
    tab20_idx  = 0
    for samp in sample_list:
        matched = next(
            (k for k in patient_colors_sc
             if k.upper() in str(samp).upper() or
                str(samp).upper() in k.upper()),
            None)
        if matched:
            samp_color_map[samp] = patient_colors_sc[matched]
        else:
            samp_color_map[samp] = tab20_cmap(tab20_idx)
            tab20_idx += 1

    # 绘制堆叠柱状图（与 PAOTU 图4 完全一致的方式）
    x_pos1  = np.arange(len(SUBTYPES_ORDER))
    bottom1 = np.zeros(len(SUBTYPES_ORDER))

    for samp in sample_list:
        vals = ct_s_pct[samp].values
        ax_sc1.bar(x_pos1, vals, bottom=bottom1,
                   color=samp_color_map[samp],
                   label=samp,
                   width=0.72,       # ← 与 PAOTU 图4 一致
                   linewidth=0,      # ← 无边线
                   edgecolor='white')
        # 比例 ≥ 7% 才标注（与 PAOTU 风格一致）
        for xi, (b, v) in enumerate(zip(bottom1, vals)):
            if v >= 7:
                ax_sc1.text(xi, b + v / 2, f'{v:.0f}%',
                            ha='center', va='center',
                            fontsize=8, color='white',
                            fontweight='bold')
        bottom1 += vals

    # X 轴标签（亚群名 + 细胞数，颜色与亚群一致）
    ax_sc1.set_xticks(x_pos1)
    ax_sc1.set_xticklabels(
        [f'{st}\n(n={int((subtypes_arr == st).sum()):,})'
         for st in SUBTYPES_ORDER],
        fontsize=10.5)
    for tick, st in zip(ax_sc1.get_xticklabels(), SUBTYPES_ORDER):
        tick.set_color(SUBTYPE_PALETTE[st])
        tick.set_fontweight('bold')

    ax_sc1.set_ylabel('Cell Proportion (%)', fontsize=11)
    ax_sc1.set_ylim(0, 110)
    ax_sc1.set_xlim(-0.5, len(SUBTYPES_ORDER) - 0.5)
    for spine in ['top', 'right']:
        ax_sc1.spines[spine].set_visible(False)

    # 图例（参考 PAOTU 图4：放在右侧）
    n_samps = len(sample_list)
    ax_sc1.legend(
        title=sample_col.replace('_', ' ').title(),
        fontsize=8.5, title_fontsize=9.5,
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),   # ← 与 PAOTU 图4 一致
        frameon=False,
        ncol=max(1, n_samps // 14))

    ax_sc1.set_title(
        f'Composition by {sample_col.replace("_", " ").title()}',
        fontsize=11, fontweight='bold', pad=8)

    # ════════════════════════════════════════════════════════════════════════
    #  Panel 2：Stacked barplot by tissue_type（参考 PAOTU tissue_colors）
    # ════════════════════════════════════════════════════════════════════════
    ax_sc2 = axes_sc[1]
    ax_sc2.set_facecolor('white')

    ct_t = (obs_sc.groupby(['mac_subtype', tissue_col])
            .size().unstack(fill_value=0))
    ct_t_pct = ct_t.div(ct_t.sum(axis=1), axis=0) * 100

    for st in SUBTYPES_ORDER:
        if st not in ct_t_pct.index:
            ct_t_pct.loc[st] = 0.0
    ct_t_pct = ct_t_pct.loc[SUBTYPES_ORDER]

    tissue_list = ct_t_pct.columns.tolist()

    # 组织类型配色：优先匹配 tissue_colors_sc，否则用备用色
    fallback_tissue = ['#8E44AD', '#16A085', '#C0392B', '#2980B9']
    tissue_color_map = {}
    fb_idx = 0
    for tis in tissue_list:
        matched = next(
            (k for k in tissue_colors_sc
             if k.upper() in str(tis).upper() or
                str(tis).upper() in k.upper()),
            None)
        if matched:
            tissue_color_map[tis] = tissue_colors_sc[matched]
        else:
            tissue_color_map[tis] = fallback_tissue[fb_idx % len(fallback_tissue)]
            fb_idx += 1

    x_pos2  = np.arange(len(SUBTYPES_ORDER))
    bottom2 = np.zeros(len(SUBTYPES_ORDER))

    for tis in tissue_list:
        vals = ct_t_pct[tis].values
        ax_sc2.bar(x_pos2, vals, bottom=bottom2,
                   color=tissue_color_map[tis],
                   label=tis,
                   width=0.72,
                   linewidth=0,
                   edgecolor='white')
        for xi, (b, v) in enumerate(zip(bottom2, vals)):
            if v >= 7:
                ax_sc2.text(xi, b + v / 2, f'{v:.0f}%',
                            ha='center', va='center',
                            fontsize=9, color='white',
                            fontweight='bold')
        bottom2 += vals

    # 分组背景色 + 分隔线（参考 PAOTU 图4 的 axvspan 风格）
    # 找 DLBCL 和 Tonsil 各自占的亚群（按 tissue 多数决）
    dominant_tissue = ct_t_pct.idxmax(axis=1)
    dlbcl_sts  = [i for i, st in enumerate(SUBTYPES_ORDER)
                  if 'DLBCL' in str(dominant_tissue.get(st, ''))]
    tonsil_sts = [i for i, st in enumerate(SUBTYPES_ORDER)
                  if 'Tonsil' in str(dominant_tissue.get(st, '')) or
                     'Normal' in str(dominant_tissue.get(st, ''))]

    # 若所有亚群都有 DLBCL/Tonsil 混合，则不画背景色，只画图例
    if dlbcl_sts:
        ax_sc2.axvspan(
            min(dlbcl_sts) - 0.5, max(dlbcl_sts) + 0.5,
            alpha=0.06, color='#0279EE', zorder=0)
    if tonsil_sts:
        ax_sc2.axvspan(
            min(tonsil_sts) - 0.5, max(tonsil_sts) + 0.5,
            alpha=0.06, color='#FF9400', zorder=0)

    ax_sc2.set_xticks(x_pos2)
    ax_sc2.set_xticklabels(
        [f'{st}\n(n={int((subtypes_arr == st).sum()):,})'
         for st in SUBTYPES_ORDER],
        fontsize=10.5)
    for tick, st in zip(ax_sc2.get_xticklabels(), SUBTYPES_ORDER):
        tick.set_color(SUBTYPE_PALETTE[st])
        tick.set_fontweight('bold')

    ax_sc2.set_ylabel('Cell Proportion (%)', fontsize=11)
    ax_sc2.set_ylim(0, 110)
    ax_sc2.set_xlim(-0.5, len(SUBTYPES_ORDER) - 0.5)
    for spine in ['top', 'right']:
        ax_sc2.spines[spine].set_visible(False)

    # 图例（与 PAOTU 图4 一致）
    handles_t = [mpatches.Patch(color=tissue_color_map[tis], label=tis)
                 for tis in tissue_list]
    ax_sc2.legend(
        handles=handles_t,
        title=tissue_col.replace('_', ' ').title(),
        fontsize=9.5, title_fontsize=10,
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),
        frameon=False)

    ax_sc2.set_title(
        f'Composition by {tissue_col.replace("_", " ").title()}',
        fontsize=11, fontweight='bold', pad=8)

    # ── 全图标题（参考 PAOTU 图4 suptitle 风格）─────────────────────────────
    fig_sc.suptitle(
        'Supplementary C — Cell Subtype Composition\n'
        'by Patient/Sample and Tissue Type (DLBCL vs Tonsil)',
        fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()
    save_fig(fig_sc, 'supp_C_sample_composition')


# =============================================================================
#  Supplementary D（重写）：Malignant B ↔ Myeloid Communication
#  数据来源：volcano_data_malignant_vs_normal_monomac.csv
#  列结构：LR_pair, ligand_complex, receptor_complex,
#          strength_mal, strength_nor, Delta_Magnitude,
#          pval, neg_log10_pval, category
#
#  Panel 1：Malignant-enriched vs Normal-enriched LR pair 数量 barplot
#  Panel 2：Top LR pairs strength_mal dot/barplot（按 Delta_Magnitude 排序）
#  Panel 3：Receptor gene 在各 myeloid subtype 的平均表达 heatmap
# =============================================================================

import os, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sps
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
warnings.filterwarnings('ignore')

# ── 路径（沿用主脚本变量）────────────────────────────────────────────────────
CELLCOMM_DIR = os.path.join(os.path.dirname(TRAJ_DIR), 'cellcomm')
long_path_v1 = os.path.join(CELLCOMM_DIR,
                             'volcano_data_malignant_vs_normal_monomac.csv')

if need_redraw('supp_D_cellcomm_strength'):
    print("\n[Supp D] Malignant B ↔ Myeloid communication strength...")

    if not os.path.exists(long_path_v1):
        print(f"  SKIP: {long_path_v1} not found")
    else:
        lr_df = pd.read_csv(long_path_v1)
        print(f"  Loaded: {lr_df.shape[0]} LR pairs")
        print(f"  Columns: {lr_df.columns.tolist()}")
        print(f"  Category counts:\n{lr_df['category'].value_counts()}")

        # ── 基础清洗 ──────────────────────────────────────────────────────────
        lr_df['strength_mal']      = pd.to_numeric(lr_df['strength_mal'],
                                                    errors='coerce').fillna(0)
        lr_df['strength_nor']      = pd.to_numeric(lr_df['strength_nor'],
                                                    errors='coerce').fillna(0)
        lr_df['Delta_Magnitude']   = pd.to_numeric(lr_df['Delta_Magnitude'],
                                                    errors='coerce').fillna(0)
        lr_df['neg_log10_pval']    = pd.to_numeric(lr_df['neg_log10_pval'],
                                                    errors='coerce').fillna(0)
        lr_df['pval']              = pd.to_numeric(lr_df['pval'],
                                                    errors='coerce').fillna(1)

        # 分组
        mal_enr = lr_df[lr_df['category'] == 'Malignant-enriched'].copy()
        nor_enr = lr_df[lr_df['category'] == 'Normal-enriched'].copy()
        ns_df   = lr_df[~lr_df['category'].isin(
                        ['Malignant-enriched', 'Normal-enriched'])].copy()

        print(f"  Malignant-enriched: {len(mal_enr)}, "
              f"Normal-enriched: {len(nor_enr)}, NS: {len(ns_df)}")

        # ── 提取 receptor 基因，查询各亚型平均表达 ───────────────────────────
        def extract_genes(s):
            return [p.strip() for p in str(s).split('_')
                    if p.strip() and p.strip() not in ('nan', '')]

        # 只用 Malignant-enriched 的 receptor genes
        rec_genes_all = []
        for rec in mal_enr['receptor_complex']:
            rec_genes_all.extend(extract_genes(rec))
        rec_genes = list(dict.fromkeys(rec_genes_all))   # 去重保序

        lig_genes_all = []
        for lig in mal_enr['ligand_complex']:
            lig_genes_all.extend(extract_genes(lig))
        lig_genes = list(dict.fromkeys(lig_genes_all))

        avail_rec = [g for g in rec_genes if g in adata_mac.var_names]
        avail_lig = [g for g in lig_genes if g in adata_mac.var_names]
        print(f"  Receptor genes in adata: {len(avail_rec)} / {len(rec_genes)}")
        print(f"  Ligand genes in adata:   {len(avail_lig)} / {len(lig_genes)}")

        # 各亚型平均表达矩阵（receptor）
        def mean_expr_by_subtype(genes, adata, subtypes_order):
            """返回 DataFrame: index=gene, columns=subtype"""
            genes_ok = [g for g in genes if g in adata.var_names]
            if not genes_ok:
                return pd.DataFrame()
            X = adata[:, genes_ok].layers['log1p_norm']
            if sps.issparse(X):
                X = X.toarray()
            df = pd.DataFrame(X, index=adata.obs_names,
                               columns=genes_ok)
            df['subtype'] = adata.obs['mac_subtype'].values
            mean_df = df.groupby('subtype')[genes_ok].mean()
            # 保证列顺序
            mean_df = mean_df.reindex(
                [s for s in subtypes_order if s in mean_df.index])
            return mean_df.T   # index=gene, columns=subtype

        rec_mean = mean_expr_by_subtype(avail_rec[:40], adata_mac, SUBTYPES_ORDER)
        lig_mean = mean_expr_by_subtype(avail_lig[:20], adata_mac, SUBTYPES_ORDER)

        # ── 图布局 ────────────────────────────────────────────────────────────
        fig_sd = plt.figure(figsize=(22, 8))
        fig_sd.patch.set_facecolor('white')
        gs_sd = fig_sd.add_gridspec(
            1, 3, wspace=0.42, width_ratios=[0.85, 1.0, 1.4])

        # ════════════════════════════════════════════════════════════════════
        #  Panel 1：Category summary barplot + Top LR pairs dot plot
        # ════════════════════════════════════════════════════════════════════
        ax_sd1 = fig_sd.add_subplot(gs_sd[0])
        ax_sd1.set_facecolor('white')

        # 上半：category 数量 bar
        gs_inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=gs_sd[0],
            height_ratios=[0.38, 0.62], hspace=0.45)
        ax_sd1a = fig_sd.add_subplot(gs_inner[0])   # 数量 bar
        ax_sd1b = fig_sd.add_subplot(gs_inner[1])   # top LR strength scatter
        ax_sd1.set_visible(False)

        # ── Panel 1a：category 数量 ──────────────────────────────────────────
        ax_sd1a.set_facecolor('white')
        cat_counts = lr_df['category'].value_counts()
        cat_order  = ['Malignant-enriched', 'Normal-enriched'] + \
                     [c for c in cat_counts.index
                      if c not in ('Malignant-enriched', 'Normal-enriched')]
        cat_order  = [c for c in cat_order if c in cat_counts.index]
        cat_colors = {'Malignant-enriched': '#D73027',
                      'Normal-enriched':    '#4575B4'}
        bar_colors_cat = [cat_colors.get(c, '#AAAAAA') for c in cat_order]
        bar_vals_cat   = [cat_counts[c] for c in cat_order]

        bars_cat = ax_sd1a.bar(range(len(cat_order)), bar_vals_cat,
                                color=bar_colors_cat, edgecolor='white',
                                linewidth=0.5, width=0.6)
        for bar, v in zip(bars_cat, bar_vals_cat):
            ax_sd1a.text(bar.get_x() + bar.get_width() / 2,
                         bar.get_height() + max(bar_vals_cat) * 0.02,
                         str(v), ha='center', va='bottom',
                         fontsize=10, fontweight='bold', color='#333333')

        ax_sd1a.set_xticks(range(len(cat_order)))
        def wrap_label(s, max_len=12):
            """超过 max_len 字符则在空格或连字符处换行"""
            import textwrap
            # 先把连字符替换为 "连字符+空格" 以便 textwrap 识别断点
            s_proc = s.replace('-', '- ')
            wrapped = textwrap.fill(s_proc, width=max_len)
            # 还原多余空格（连字符后不加空格）
            wrapped = wrapped.replace('- \n', '-\n').replace('- ', '-')
            return wrapped

        ax_sd1a.set_xticklabels(
            [wrap_label(c) for c in cat_order],
            fontsize=8.5, rotation=0, ha='center')

        for tick, col in zip(ax_sd1a.get_xticklabels(), bar_colors_cat):
            tick.set_color(col); tick.set_fontweight('bold')
        ax_sd1a.set_ylabel('# LR Pairs', fontsize=9)
        ax_sd1a.set_title('LR Pair Category', fontsize=10,
                           fontweight='bold', pad=6)
        ax_sd1a.spines[['top', 'right']].set_visible(False)

        # ── Panel 1b：Top Malignant-enriched LR pairs（strength scatter）────
        ax_sd1b.set_facecolor('white')
        top_mal = mal_enr.nlargest(
            min(20, len(mal_enr)), 'Delta_Magnitude').copy()
        top_mal = top_mal.sort_values('Delta_Magnitude', ascending=True)
        y_pos_b = np.arange(len(top_mal))

        # 颜色映射：neg_log10_pval
        norm_b  = plt.Normalize(vmin=top_mal['neg_log10_pval'].min(),
                                 vmax=top_mal['neg_log10_pval'].max())
        cmap_b  = plt.cm.RdYlBu
        dot_cols = [cmap_b(norm_b(v)) for v in top_mal['neg_log10_pval']]
        dot_sizes = (top_mal['strength_mal'] /
                     (top_mal['strength_mal'].max() + 1e-9) * 200 + 30)

        ax_sd1b.scatter(top_mal['Delta_Magnitude'].values, y_pos_b,
                         c=dot_cols, s=dot_sizes.values,
                         edgecolors='white', linewidths=0.8,
                         zorder=3)
        ax_sd1b.axvline(0, color='#AAAAAA', lw=0.8, ls='--', zorder=1)

        ax_sd1b.set_yticks(y_pos_b)
        ax_sd1b.set_yticklabels(top_mal['LR_pair'].values, fontsize=7.5)
        ax_sd1b.set_xlabel('ΔMagnitude\n(Malignant − Normal)', fontsize=8.5)
        ax_sd1b.set_title(
            f'Top {len(top_mal)} Malignant-enriched\nLR Pairs',
            fontsize=10, fontweight='bold', pad=6)
        ax_sd1b.spines[['top', 'right']].set_visible(False)
        ax_sd1b.grid(axis='x', color='#EEEEEE', lw=0.7, zorder=0)

        sm_b = plt.cm.ScalarMappable(cmap=cmap_b, norm=norm_b)
        sm_b.set_array([])
        cb_b = plt.colorbar(sm_b, ax=ax_sd1b, shrink=0.55,
                             pad=0.02, aspect=18)
        cb_b.set_label('−log₁₀(p)', fontsize=7.5)
        cb_b.ax.tick_params(labelsize=7)

        # ════════════════════════════════════════════════════════════════════
        #  Panel 2：Malignant-enriched vs Normal-enriched
        #           strength_mal / strength_nor scatter（bubble plot）
        # ════════════════════════════════════════════════════════════════════
        ax_sd2 = fig_sd.add_subplot(gs_sd[1])
        ax_sd2.set_facecolor('white')

        color_map_cat = {
            'Malignant-enriched': '#D73027',
            'Normal-enriched':    '#4575B4',
        }
        size_map_cat = {
            'Malignant-enriched': 45,
            'Normal-enriched':    45,
        }
        alpha_map = {
            'Malignant-enriched': 0.75,
            'Normal-enriched':    0.75,
        }

        # NS 先画（灰色底层）
        ax_sd2.scatter(ns_df['strength_nor'], ns_df['strength_mal'],
                       c='#CCCCCC', s=18, alpha=0.30,
                       linewidths=0, zorder=1,
                       label=f'NS (n={len(ns_df)})')

        for cat, sub_cat in [('Normal-enriched',    nor_enr),
                              ('Malignant-enriched', mal_enr)]:
            if len(sub_cat) == 0:
                continue
            ax_sd2.scatter(
                sub_cat['strength_nor'], sub_cat['strength_mal'],
                c=color_map_cat[cat],
                s=size_map_cat[cat],
                alpha=alpha_map[cat],
                edgecolors='white', linewidths=0.5,
                zorder=3, label=f'{cat} (n={len(sub_cat)})')

        # 对角线 y = x
        xy_max = max(lr_df['strength_mal'].max(),
                     lr_df['strength_nor'].max()) * 1.05
        ax_sd2.plot([0, xy_max], [0, xy_max],
                    color='#888888', lw=0.9, ls='--',
                    zorder=2, label='y = x')

        # 标注 top 5 Malignant-enriched
        top5 = mal_enr.nlargest(5, 'Delta_Magnitude')
        for _, row in top5.iterrows():
            ax_sd2.annotate(
                row['LR_pair'],
                xy=(row['strength_nor'], row['strength_mal']),
                xytext=(8, 4), textcoords='offset points',
                fontsize=7.5, color='#D73027', fontweight='bold',
                arrowprops=dict(arrowstyle='-', color='#D73027',
                                lw=0.6))

        ax_sd2.set_xlabel('Interaction Strength — Normal/Tonsil',
                           fontsize=10)
        ax_sd2.set_ylabel('Interaction Strength — Malignant (DLBCL)',
                           fontsize=10)
        ax_sd2.set_xlim(-0.005 * xy_max, xy_max)
        ax_sd2.set_ylim(-0.005 * xy_max, xy_max)
        ax_sd2.legend(fontsize=8.5, framealpha=0.9,
                       edgecolor='none', loc='upper left')
        ax_sd2.spines[['top', 'right']].set_visible(False)
        ax_sd2.set_title(
            'Malignant vs Normal Interaction Strength\n'
            '(each dot = one LR pair; above diagonal = Malignant-enriched)',
            fontsize=10, fontweight='bold', pad=8)

        # ════════════════════════════════════════════════════════════════════
        #  Panel 3：Receptor gene mean expression heatmap per myeloid subtype
        # ════════════════════════════════════════════════════════════════════
        ax_sd3 = fig_sd.add_subplot(gs_sd[2])
        ax_sd3.set_facecolor('white')

        if rec_mean.shape[0] > 0 and rec_mean.shape[1] > 0:
            from scipy.stats import zscore as scipy_zscore

            # 最多展示 top 30 基因（按跨亚型最大均值排序）
            rec_mean_show = rec_mean.copy()
            rec_mean_show['max_val'] = rec_mean_show.max(axis=1)
            rec_mean_show = (rec_mean_show
                             .sort_values('max_val', ascending=False)
                             .drop(columns='max_val')
                             .head(30))

            # Z-score（按行）
            rec_z = rec_mean_show.apply(
                lambda row: scipy_zscore(row)
                if row.std() > 0 else row,
                axis=1)
            rec_z = pd.DataFrame(rec_z.tolist(),
                                  index=rec_mean_show.index,
                                  columns=rec_mean_show.columns)

            cmap_sd3 = LinearSegmentedColormap.from_list(
                'rec_expr', ['#2166AC', '#F7F7F7', '#D6604D'])

            gs_inner3 = gridspec.GridSpecFromSubplotSpec(
                1, 2, subplot_spec=gs_sd[2],
                width_ratios=[1, 0.04], wspace=0.03)
            ax_hm3  = fig_sd.add_subplot(gs_inner3[0])
            ax_cb3  = fig_sd.add_subplot(gs_inner3[1])
            ax_sd3.set_visible(False)
            ax_hm3.set_facecolor('white')

            im_sd3 = ax_hm3.imshow(
                rec_z.values, aspect='auto',
                cmap=cmap_sd3, vmin=-2, vmax=2,
                interpolation='nearest')

            # 列（亚型）颜色标签
            ax_hm3.set_xticks(range(len(rec_z.columns)))
            ax_hm3.set_xticklabels(rec_z.columns,
                                    rotation=35, ha='right', fontsize=10)
            for tick, st in zip(ax_hm3.get_xticklabels(),
                                 rec_z.columns):
                tick.set_color(SUBTYPE_PALETTE.get(st, '#333333'))
                tick.set_fontweight('bold')

            # 行（基因）标签
            ax_hm3.set_yticks(range(len(rec_z.index)))
            ax_hm3.set_yticklabels(rec_z.index, fontsize=8)

            # 基因标签颜色：是否在 top Malignant-enriched receptor 里
            top_rec_genes = set()
            for rec in mal_enr.nlargest(
                    min(15, len(mal_enr)),
                    'Delta_Magnitude')['receptor_complex']:
                top_rec_genes.update(extract_genes(rec))

            for tick in ax_hm3.get_yticklabels():
                gene = tick.get_text()
                if gene in top_rec_genes:
                    tick.set_color('#D73027')
                    tick.set_fontweight('bold')
                else:
                    tick.set_color('#555555')

            # 网格分隔线
            for i in range(len(rec_z.index) + 1):
                ax_hm3.axhline(i - 0.5, color='white', lw=0.4)
            for j in range(len(rec_z.columns) + 1):
                ax_hm3.axvline(j - 0.5, color='white', lw=0.8)

            plt.colorbar(im_sd3, cax=ax_cb3, label='Z-score')
            ax_cb3.tick_params(labelsize=8)

            ax_hm3.set_title(
                'Malignant-enriched Receptor Gene\n'
                'Expression per Myeloid Subtype\n'
                '(red gene = top ΔMagnitude receptor; Z-score)',
                fontsize=10, fontweight='bold', pad=8)

            # 图例：红色 = top receptor
            legend_p3 = [
                mpatches.Patch(color='#D73027',
                               label='Top ΔMag receptor gene'),
                mpatches.Patch(color='#555555',
                               label='Other receptor gene'),
            ]
            ax_hm3.legend(handles=legend_p3,
                           fontsize=7.5, loc='lower right',
                           bbox_to_anchor=(1.0, -0.18),
                           framealpha=0.9, edgecolor='none',
                           ncol=2)
        else:
            ax_sd3.text(0.5, 0.5,
                        'No receptor genes found in adata_mac',
                        ha='center', va='center',
                        fontsize=10, color='#888888',
                        transform=ax_sd3.transAxes)

        # ── 全图标题 ──────────────────────────────────────────────────────────
        fig_sd.suptitle(
            'Supplementary D — Malignant B ↔ Myeloid Subtype Communication\n'
            '(Malignant-enriched vs Normal-enriched LR pairs; '
            'receptor expression by subtype)',
            fontsize=13, fontweight='bold', y=1.02)

        plt.tight_layout()
        save_fig(fig_sd, 'supp_D_cellcomm_strength')

# =============================================================================
print()
print("=" * 60)
print("Supp D DONE")
for fmt in ['png', 'svg']:
    p = os.path.join(TRAJ_DIR, f'supp_D_cellcomm_strength.{fmt}')
    status = '✓' if os.path.exists(p) else '✗ MISSING'
    print(f"  [{status}] {p}")
print("=" * 60)
