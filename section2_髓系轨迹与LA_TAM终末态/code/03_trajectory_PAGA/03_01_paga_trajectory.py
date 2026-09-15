# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# =============================================================================
#  单细胞 Monocyte/Macrophage 分化轨迹分析 — 完整复现脚本（最终版）
#  生成 11 张图：
#    Fig 1  monomac_paga_trajectory
#    Fig 2  monomac_trajectory_directed
#    Fig 3  monomac_gene_dynamics
#    Fig 4  scenic_tf_dynamics
#    Fig 5  scenic_tf_heatmap
#    Fig 6  receptor_pseudotime_heatmap_path
#    Fig 7  receptor_pseudotime_lineplots_path
#    Fig 8  volcano_IFN_TAM_vs_LA_TAM          → DEG_DIR
#    Fig 9  enrichment_barplot_IFN_TAM_vs_LA_TAM → DEG_DIR
#    Fig 10 cellrank_gene_relay_heatmap
#    Fig 11 cellrank_fate_trend_plots
# =============================================================================

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
H5AD_PATH    = translate(os.path.join(DATA_DIR, 'adata_processed.h5ad'))
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
#  图 1：PAGA 拓扑图 + PAGA-UMAP 亚型 + PAGA-UMAP 伪时间 (1×3)
# =============================================================================
print("\n[Fig 1] PAGA topology + UMAP panels...")

fig1 = plt.figure(figsize=(20, 6.5))
fig1.patch.set_facecolor('white')
gs1  = fig1.add_gridspec(1, 3, wspace=0.35)

# ── Panel A：PAGA 拓扑图 ──────────────────────────────────────────────────────
ax_paga = fig1.add_subplot(gs1[0])

threshold_edge = 0.15
for i, st_i in enumerate(SUBTYPES_ORDER):
    for j, st_j in enumerate(SUBTYPES_ORDER):
        if j <= i: continue
        w = conn[i, j]
        if w < threshold_edge: continue
        x0, y0 = paga_pos[i]; x1, y1 = paga_pos[j]
        ax_paga.plot([x0, x1], [y0, y1], '-', color='#888888',
                     lw=w * 8, alpha=min(0.9, 0.3 + w * 0.7), zorder=1)
        mx, my = (x0+x1)/2, (y0+y1)/2
        ax_paga.text(mx, my, f'{w:.2f}', fontsize=7,
                     ha='center', va='center', color='#555',
                     bbox=dict(boxstyle='round,pad=0.1', fc='white',
                               alpha=0.7, lw=0))

for i, st in enumerate(SUBTYPES_ORDER):
    x, y = paga_pos[i]
    ax_paga.scatter(x, y, s=node_sizes_scaled[i], c=SUBTYPE_PALETTE[st],
                    zorder=3, edgecolors='white', linewidths=2)
    ax_paga.text(x, y+0.08, st, ha='center', va='bottom',
                 fontsize=10, fontweight='bold', color=SUBTYPE_PALETTE[st],
                 bbox=dict(boxstyle='round,pad=0.2', fc='white',
                           alpha=0.8, lw=0))
    ax_paga.text(x, y-0.08, f'n={n_cells[st]}', ha='center', va='top',
                 fontsize=8, color='#555')

ax_paga.set_title('PAGA Topology Graph\n(edge weight = connectivity)',
                  fontsize=12, fontweight='bold')
ax_paga.set_xlim(paga_pos[:,0].min()-0.3, paga_pos[:,0].max()+0.3)
ax_paga.set_ylim(paga_pos[:,1].min()-0.3, paga_pos[:,1].max()+0.3)
ax_paga.axis('off')

# ── Panel B：PAGA-UMAP 亚型 ───────────────────────────────────────────────────
ax_sub = fig1.add_subplot(gs1[1])

for st in SUBTYPES_ORDER:
    mask = subtypes_arr == st
    ax_sub.scatter(umap_coords[mask,0], umap_coords[mask,1],
                   c=SUBTYPE_PALETTE[st], s=18, alpha=0.75,
                   linewidths=0, label=st, rasterized=True)

for st in SUBTYPES_ORDER:
    mask = subtypes_arr == st
    cx = umap_coords[mask,0].mean()
    cy = umap_coords[mask,1].mean()
    ax_sub.text(cx, cy, st, fontsize=8.5, fontweight='bold',
                ha='center', va='center', color='white',
                bbox=dict(boxstyle='round,pad=0.25',
                          fc=SUBTYPE_PALETTE[st], alpha=0.85, lw=0))

ax_sub.set_title('PAGA-initialized UMAP\n(Cell Subtypes)',
                 fontsize=12, fontweight='bold')
ax_sub.set_xlabel('UMAP1', fontsize=10)
ax_sub.set_ylabel('UMAP2', fontsize=10)
ax_sub.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_sub.spines.values(): sp.set_visible(False)
legend_handles = [mpatches.Patch(color=SUBTYPE_PALETTE[st], label=st)
                  for st in SUBTYPES_ORDER]
ax_sub.legend(handles=legend_handles, fontsize=8, loc='lower right',
              framealpha=0.85, edgecolor='none',
              title='Subtype', title_fontsize=8.5)

# ── Panel C：PAGA-UMAP 伪时间（固定轨迹箭头）────────────────────────────────
ax_pt = fig1.add_subplot(gs1[2])

sc_pt = ax_pt.scatter(umap_coords[:,0], umap_coords[:,1],
                      c=pt_arr, cmap='viridis', s=18, alpha=0.85,
                      linewidths=0, rasterized=True)

# trace.txt 原始固定4条箭头
for src, tgt in trajectory_edges:
    x0, y0 = centroids[src]; x1, y1 = centroids[tgt]
    ax_pt.annotate('', xy=(x1,y1), xytext=(x0,y0),
                   arrowprops=dict(arrowstyle='->', color='white',
                                   lw=2.5, mutation_scale=18))
    ax_pt.annotate('', xy=(x1,y1), xytext=(x0,y0),
                   arrowprops=dict(arrowstyle='->', color='#333333',
                                   lw=1.5, mutation_scale=16))

for st in SUBTYPES_ORDER:
    cx, cy = centroids[st]
    ax_pt.text(cx, cy, st, fontsize=8, fontweight='bold',
               ha='center', va='center', color='white',
               bbox=dict(boxstyle='round,pad=0.2',
                         fc=SUBTYPE_PALETTE[st], alpha=0.85, lw=0))

cbar1 = fig1.colorbar(sc_pt, ax=ax_pt, shrink=0.7, pad=0.02)
cbar1.set_label('Pseudotime\n(DPT)', fontsize=9)
cbar1.ax.tick_params(labelsize=8)
ax_pt.set_title('PAGA-initialized UMAP\n(DPT Pseudotime)',
                fontsize=12, fontweight='bold')
ax_pt.set_xlabel('UMAP1', fontsize=10)
ax_pt.set_ylabel('UMAP2', fontsize=10)
ax_pt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_pt.spines.values(): sp.set_visible(False)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig1.savefig(os.path.join(TRAJ_DIR, f'monomac_paga_trajectory.{fmt}'),
                 dpi=DPI, bbox_inches='tight', format=fmt)
plt.close(fig1)
print("  ✓ monomac_paga_trajectory saved")

# =============================================================================
#  图 2：有向轨迹 UMAP + 伪时间梯度 (1×2)
# =============================================================================
print("\n[Fig 2] Directed trajectory UMAP...")

fig2, axes2 = plt.subplots(1, 2, figsize=(15, 6.5))
fig2.patch.set_facecolor('white')

# ── Panel A：有向箭头 ─────────────────────────────────────────────────────────
ax2a = axes2[0]

for st in SUBTYPES_ORDER:
    mask = subtypes_arr == st
    ax2a.scatter(umap_coords[mask,0], umap_coords[mask,1],
                 c=SUBTYPE_PALETTE[st], s=20, alpha=0.65,
                 linewidths=0, rasterized=True)

for src, tgt, w in directed_edges:
    x0,y0 = centroids[src]; x1,y1 = centroids[tgt]
    dx,dy = x1-x0, y1-y0
    shrink = 0.22
    xs,ys = x0+dx*shrink, y0+dy*shrink
    xe,ye = x1-dx*shrink, y1-dy*shrink
    lw2 = 1.5 + w*4
    ax2a.annotate('', xy=(xe,ye), xytext=(xs,ys),
                  arrowprops=dict(arrowstyle='->', color='white',
                                  lw=lw2+1.5, mutation_scale=20))
    ax2a.annotate('', xy=(xe,ye), xytext=(xs,ys),
                  arrowprops=dict(arrowstyle='->', color='#333333',
                                  lw=lw2, mutation_scale=18))
    mx,my = (xs+xe)/2, (ys+ye)/2
    ax2a.text(mx, my, f'{w:.2f}', fontsize=7.5, ha='center', va='center',
              color='#333',
              bbox=dict(boxstyle='round,pad=0.15', fc='white',
                        alpha=0.8, lw=0))

for st in SUBTYPES_ORDER:
    cx,cy = centroids[st]
    n = n_cells[st]
    ax2a.text(cx, cy, f'{st}\n(n={n})', fontsize=8.5, fontweight='bold',
              ha='center', va='center', color='white',
              bbox=dict(boxstyle='round,pad=0.3', fc=SUBTYPE_PALETTE[st],
                        alpha=0.9, lw=0.5, edgecolor='white'))

ax2a.set_title('Mono/Mac Differentiation Trajectory\n'
               '(PAGA + DPT, arrows = direction)',
               fontsize=12, fontweight='bold')
ax2a.set_xlabel('UMAP1', fontsize=10)
ax2a.set_ylabel('UMAP2', fontsize=10)
ax2a.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax2a.spines.values(): sp.set_visible(False)

# ── Panel B：伪时间梯度 ───────────────────────────────────────────────────────
ax2b = axes2[1]

sc2b = ax2b.scatter(umap_coords[:,0], umap_coords[:,1],
                    c=pt_arr, cmap='viridis', s=20, alpha=0.85,
                    linewidths=0, rasterized=True, vmin=0, vmax=1)

for src, tgt, w in directed_edges:
    x0,y0 = centroids[src]; x1,y1 = centroids[tgt]
    dx,dy = x1-x0, y1-y0
    shrink = 0.22
    xs,ys = x0+dx*shrink, y0+dy*shrink
    xe,ye = x1-dx*shrink, y1-dy*shrink
    ax2b.annotate('', xy=(xe,ye), xytext=(xs,ys),
                  arrowprops=dict(arrowstyle='->', color='white',
                                  lw=3.5, mutation_scale=20))
    ax2b.annotate('', xy=(xe,ye), xytext=(xs,ys),
                  arrowprops=dict(arrowstyle='->', color='#ff4444',
                                  lw=2.0, mutation_scale=18))

for st in SUBTYPES_ORDER:
    cx,cy = centroids[st]
    pt_val = pt_means[st]
    ax2b.text(cx, cy, f'{st}\nPT={pt_val:.2f}', fontsize=8, fontweight='bold',
              ha='center', va='center', color='white',
              bbox=dict(boxstyle='round,pad=0.25', fc=SUBTYPE_PALETTE[st],
                        alpha=0.9, lw=0.5, edgecolor='white'))

cbar2 = fig2.colorbar(sc2b, ax=ax2b, shrink=0.7, pad=0.02)
cbar2.set_label('DPT Pseudotime', fontsize=9)
cbar2.ax.tick_params(labelsize=8)
ax2b.set_title('DPT Pseudotime Projection\n'
               '(red arrows = differentiation direction)',
               fontsize=12, fontweight='bold')
ax2b.set_xlabel('UMAP1', fontsize=10)
ax2b.set_ylabel('UMAP2', fontsize=10)
ax2b.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax2b.spines.values(): sp.set_visible(False)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig2.savefig(os.path.join(TRAJ_DIR, f'monomac_trajectory_directed.{fmt}'),
                 dpi=DPI, bbox_inches='tight', format=fmt)
plt.close(fig2)
print("  ✓ monomac_trajectory_directed saved")


# ══════════════════════════════════════════════════════════════════════════════
# Fig3 前置参数：Gene Dynamics Heatmap Along Pseudotime
# ══════════════════════════════════════════════════════════════════════════════

P03 = {

    # ── 基因分组 ──────────────────────────────────────────────────────────────
    "gene_groups": {
        'Mono':    ['FCN1', 'S100A9', 'S100A8'],
        'IFN_TAM': ['MT1H', 'CCL8',   'CCL2'],
        'LA_TAM':  ['APOE', 'CCL18',  'CTSD'],
        'DC_2':    ['CST3', 'CLEC9A', 'DNASE1L3'],
        'DC_1':    ['CD1C', 'CLEC10A','LTB'],
    },

    # ── 分组颜色 ──────────────────────────────────────────────────────────────
    "group_colors": {
        'Mono':    '#3498DB',
        'IFN_TAM': '#E67E22',
        'LA_TAM':  '#2ECC71',
        'DC_2':    '#9B59B6',
        'DC_1':    '#E74C3C',
    },

    # ── 伪时间分箱 ────────────────────────────────────────────────────────────
    "n_bins"           : 80,             # 分箱数量

    # ── 平滑参数 ──────────────────────────────────────────────────────────────
    "sigma"            : 2.5,            # gaussian_filter1d 标准差（越大越平滑）

    # ── 热图颜色映射 ──────────────────────────────────────────────────────────
    "heatmap_cmap"     : 'RdBu_r',
    "heatmap_vmin"     : -2,             # Z-score 色阶下限
    "heatmap_vmax"     :  2,             # Z-score 色阶上限
    "heatmap_interp"   : 'bilinear',

    # ── 亚型色条 ──────────────────────────────────────────────────────────────
    "subtype_default_color": '#CCCCCC',  # 未在 SUBTYPE_PALETTE 中的亚型颜色

    # ── 伪时间渐变条颜色 ──────────────────────────────────────────────────────
    "pt_bar_cmap"      : 'viridis',

    # ── 图例（亚型色条上方）──────────────────────────────────────────────────
    "legend_fontsize"       : 12,
    "legend_bbox"           : (0.23, 1.05),   # 紧贴色条上方
    "legend_framealpha"     : 0.9,
    "legend_handlelength"   : 1.0,
    "legend_handletextpad"  : 0.4,
    "legend_columnspacing"  : 0.8,

    # ── 右侧分组标注 ──────────────────────────────────────────────────────────
    "group_label_fontsize"  : 15,
    "group_label_x_offset"  : 1.5,       # 标注文字距热图右边缘的偏移（data 坐标）
    "group_sep_color"       : 'white',   # 分组分隔线颜色
    "group_sep_lw"          : 1.5,       # 分组分隔线宽度

    # ── 轴标签与刻度 ──────────────────────────────────────────────────────────
    "yticklabel_fontsize"   : 15,
    "subtype_ylabel"        : 'Subtype',
    "pt_ylabel"             : 'Pseudotime',
    "xtick_fontsize"        : 10,
    "xlabel_fontsize"       : 15,
    "cbar_label_fontsize"   : 11,
    "cbar_tick_fontsize"    : 10,

    # ── 伪时间 x 轴刻度位置与标签 ────────────────────────────────────────────
    "xtick_ratios"          : [0, 0.25, 0.5, 0.75, 1.0],
    "xtick_labels"          : ['0', '0.25', '0.5', '0.75', '1.0'],

    # ── 总标题 ────────────────────────────────────────────────────────────────
    "suptitle"         : '',
    "suptitle_fontsize": 13,
    "suptitle_y"       : 0.98,

    # ── 图像尺寸与布局 ────────────────────────────────────────────────────────
    "figsize"          : (14, 8),
    "fig_facecolor"    : 'white',
    "height_ratios"    : [0.10, 0.10, 1],
    "width_ratios"     : [1, 0.1, 0.01],  # 主图 | 分组标注留白 | colorbar
    "hspace"           : 0.04,
    "wspace"           : 0.02,

    # ── 整体边距 ──────────────────────────────────────────────────────────────
    "margin_left"      : 0.10,
    "margin_right"     : 0.92,
    "margin_top"       : 0.91,
    "margin_bottom"    : 0.08,

    # ── 输出 ──────────────────────────────────────────────────────────────────
    "output_name"      : 'monomac_gene_dynamics',
    "output_dir"       : TRAJ_DIR,
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig3 绘图
# ══════════════════════════════════════════════════════════════════════════════

print("[Fig 3] Gene dynamics heatmap...")

# ── 基因过滤 ──────────────────────────────────────────────────────────────────
all_genes = [g for genes in P03["gene_groups"].values()
             for g in genes if g in adata_mac.var_names]
print(f"  Genes available: {len(all_genes)}")

# ── 表达矩阵 ──────────────────────────────────────────────────────────────────
X_mac = adata_mac.X
if sps.issparse(X_mac):
    X_mac = X_mac.toarray()

gene_idx = {g: list(adata_mac.var_names).index(g) for g in all_genes}
pt_mac   = adata_mac.obs['dpt_pseudotime'].values
st_mac   = adata_mac.obs['mac_subtype'].values

# ── 分箱统计 ──────────────────────────────────────────────────────────────────
n_bins        = P03["n_bins"]
bin_edges     = np.linspace(0, 1, n_bins + 1)
bin_centers   = (bin_edges[:-1] + bin_edges[1:]) / 2
expr_binned   = np.zeros((len(all_genes), n_bins))
bin_subtype   = []

for b in range(n_bins):
    mask = (pt_mac >= bin_edges[b]) & (pt_mac < bin_edges[b + 1])
    if mask.sum() > 0:
        for gi, gene in enumerate(all_genes):
            expr_binned[gi, b] = X_mac[mask, gene_idx[gene]].mean()
        bin_subtype.append(Counter(st_mac[mask]).most_common(1)[0][0])
    else:
        if b > 0:
            expr_binned[:, b] = expr_binned[:, b - 1]
        bin_subtype.append(bin_subtype[-1] if bin_subtype else 'Mono')

# ── 平滑与 Z-score ────────────────────────────────────────────────────────────
expr_smooth = np.array([gaussian_filter1d(expr_binned[gi], P03["sigma"])
                        for gi in range(len(all_genes))])
expr_z      = ((expr_smooth - expr_smooth.mean(axis=1, keepdims=True)) /
               (expr_smooth.std(axis=1,  keepdims=True) + 1e-9))

# ── 画布与 GridSpec ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=P03["figsize"])
fig.patch.set_facecolor(P03["fig_facecolor"])

gs3 = fig.add_gridspec(
    3, 3,
    height_ratios = P03["height_ratios"],
    width_ratios  = P03["width_ratios"],
    hspace        = P03["hspace"],
    wspace        = P03["wspace"],
)

ax_ct3  = fig.add_subplot(gs3[0, 0])
ax_pt3  = fig.add_subplot(gs3[1, 0])
ax_hm3  = fig.add_subplot(gs3[2, 0])
ax_cbar = fig.add_subplot(gs3[2, 2])

fig.suptitle(P03["suptitle"],
             fontsize=P03["suptitle_fontsize"],
             fontweight='bold',
             y=P03["suptitle_y"])

# ── 亚型色条 ──────────────────────────────────────────────────────────────────
ct_colors3 = np.array([[mcolors.to_rgb(
                            SUBTYPE_PALETTE.get(st, P03["subtype_default_color"]))
                         for st in bin_subtype]])
ax_ct3.imshow(ct_colors3, aspect='auto', interpolation='nearest')
ax_ct3.set_yticks([0])
ax_ct3.set_yticklabels([P03["subtype_ylabel"]], fontsize=P03["yticklabel_fontsize"])
ax_ct3.set_xticks([])
ax_ct3.set_xlim(-0.5, n_bins - 0.5)

# ── 图例 ──────────────────────────────────────────────────────────────────────
legend_handles3 = [mpatches.Patch(color=SUBTYPE_PALETTE[st], label=st)
                   for st in SUBTYPES_ORDER]
title_handle    = mpatches.Patch(color='none', label='Subtype:')
ax_ct3.legend(
    handles       = [title_handle] + legend_handles3,
    fontsize      = P03["legend_fontsize"],
    loc           = 'lower left',
    bbox_to_anchor= P03["legend_bbox"],
    framealpha    = P03["legend_framealpha"],
    edgecolor     = 'none',
    ncol          = len(SUBTYPES_ORDER) + 1,
    handlelength  = P03["legend_handlelength"],
    handletextpad = P03["legend_handletextpad"],
    columnspacing = P03["legend_columnspacing"],
)

# ── 伪时间渐变条 ──────────────────────────────────────────────────────────────
pt_bar3 = np.linspace(0, 1, n_bins).reshape(1, -1)
ax_pt3.imshow(pt_bar3, aspect='auto',
              cmap=P03["pt_bar_cmap"], interpolation='nearest')
ax_pt3.set_yticks([0])
ax_pt3.set_yticklabels([P03["pt_ylabel"]], fontsize=P03["yticklabel_fontsize"])
ax_pt3.set_xticks([])
ax_pt3.set_xlim(-0.5, n_bins - 0.5)

# ── 热图 ──────────────────────────────────────────────────────────────────────
im3 = ax_hm3.imshow(expr_z, aspect='auto',
                    cmap          = P03["heatmap_cmap"],
                    vmin          = P03["heatmap_vmin"],
                    vmax          = P03["heatmap_vmax"],
                    interpolation = P03["heatmap_interp"])
ax_hm3.set_xlim(-0.5, n_bins - 0.5)
ax_hm3.set_yticks(range(len(all_genes)))
ax_hm3.set_yticklabels(all_genes, fontsize=P03["yticklabel_fontsize"])

gene_to_group = {g: grp for grp, genes in P03["gene_groups"].items() for g in genes}
for tick, gene in zip(ax_hm3.get_yticklabels(), all_genes):
    grp = gene_to_group.get(gene, 'Mono')
    tick.set_color(P03["group_colors"][grp])
    tick.set_fontweight('bold')

# ── 右侧分组标注 ──────────────────────────────────────────────────────────────
group_sizes3 = [len(v) for v in P03["gene_groups"].values()]
pos3 = 0
for k, (grp, size) in enumerate(zip(P03["gene_groups"].keys(), group_sizes3)):
    mid3 = pos3 + size / 2 - 0.5
    ax_hm3.annotate(
        grp,
        xy             = (n_bins - 0.5, mid3),
        xytext         = (n_bins + P03["group_label_x_offset"], mid3),
        xycoords       = 'data',
        textcoords     = 'data',
        fontsize       = P03["group_label_fontsize"],
        color          = P03["group_colors"][grp],
        fontweight     = 'bold',
        ha             = 'left',
        va             = 'center',
        annotation_clip= False,
        arrowprops     = None,
    )
    if k < len(group_sizes3) - 1:
        ax_hm3.axhline(pos3 + size - 0.5,
                       color=P03["group_sep_color"],
                       lw   =P03["group_sep_lw"],
                       zorder=5)
    pos3 += size

# ── x 轴刻度 ──────────────────────────────────────────────────────────────────
tick_pos3 = [min(int(n_bins * t), n_bins - 1) for t in P03["xtick_ratios"]]
ax_hm3.set_xticks(tick_pos3)
ax_hm3.set_xticklabels(P03["xtick_labels"], fontsize=P03["xtick_fontsize"])
ax_hm3.set_xlabel('Pseudotime', fontsize=P03["xlabel_fontsize"])

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar3 = fig.colorbar(im3, cax=ax_cbar)
cbar3.set_label('Z-score', fontsize=P03["cbar_label_fontsize"])
cbar3.ax.tick_params(labelsize=P03["cbar_tick_fontsize"])

# ── 整体边距 ──────────────────────────────────────────────────────────────────
fig.subplots_adjust(
    left   = P03["margin_left"],
    right  = P03["margin_right"],
    top    = P03["margin_top"],
    bottom = P03["margin_bottom"],
)

# ── 输出 ──────────────────────────────────────────────────────────────────────
for fmt in ['png', 'svg']:
    fig.savefig(
        os.path.join(P03["output_dir"], f'{P03["output_name"]}.{fmt}'),
        dpi=DPI, bbox_inches='tight', format=fmt,
    )
plt.close(fig)
print("  ✓ monomac_gene_dynamics saved")


# =============================================================================
#  图 4/5：SCENIC TF 动态图 + 热图
# =============================================================================
auc_path  = os.path.join(SCENIC_DIR, 'data_AUC_matrix.csv')
meta_path = os.path.join(SCENIC_DIR, 'data_cell_metadata.csv')

if os.path.exists(auc_path) and os.path.exists(meta_path):
    print("\n[Fig 4/5] SCENIC TF plots...")

    auc_mtx      = pd.read_csv(auc_path, index_col=0)
    cell_meta_raw = pd.read_csv(meta_path, index_col='cell_id')
    cell_meta     = cell_meta_raw.rename(columns={'subtype': 'mac_subtype'})[
        ['mac_subtype', 'dpt_pseudotime']].copy()

    common_cells = auc_mtx.index.intersection(cell_meta.index)
    print(f"  Common cells: {len(common_cells)}")
    auc_sub  = auc_mtx.loc[common_cells]
    meta_sub = cell_meta.loc[common_cells]
    df_scenic = auc_sub.join(meta_sub).sort_values('dpt_pseudotime')

    # 差异 TF：IFN_TAM vs others
    results_sc = []
    for tf in auc_mtx.columns:
        i_vals = auc_sub.loc[meta_sub['mac_subtype'] == 'IFN_TAM', tf].values
        o_vals = auc_sub.loc[meta_sub['mac_subtype'] != 'IFN_TAM', tf].values
        if len(i_vals) == 0 or len(o_vals) == 0: continue
        _, p = mannwhitneyu(i_vals, o_vals, alternative='greater')
        results_sc.append({'TF': tf, 'mean_IFN_TAM': i_vals.mean(),
                           'mean_other': o_vals.mean(), 'pval': p})

    df_diff_sc = pd.DataFrame(results_sc)
    _, padj_sc, _, _ = multipletests(df_diff_sc['pval'], method='fdr_bh')
    df_diff_sc['padj'] = padj_sc
    sig_tfs_sc  = df_diff_sc[df_diff_sc['padj'] < 0.05].sort_values('padj')
    selected_tfs = sig_tfs_sc.head(10)['TF'].tolist()
    print(f"  Significant TFs: {len(sig_tfs_sc)}  |  Selected top-10: "
          f"{[t.replace('(+)','') for t in selected_tfs]}")

    # 平滑 TF 活性
    n_bins_sc  = 99; sigma_sc = 3
    pt_min_sc  = df_scenic['dpt_pseudotime'].min()
    pt_max_sc  = df_scenic['dpt_pseudotime'].max()
    pt_bins_sc = np.linspace(pt_min_sc, pt_max_sc, n_bins_sc)
    pt_ctrs_sc = (pt_bins_sc[:-1] + pt_bins_sc[1:]) / 2

    bin_subtypes_sc = []
    for i in range(len(pt_bins_sc) - 1):
        mask_b = ((df_scenic['dpt_pseudotime'] >= pt_bins_sc[i]) &
                  (df_scenic['dpt_pseudotime'] <  pt_bins_sc[i + 1]))
        bin_subtypes_sc.append(
            df_scenic.loc[mask_b, 'mac_subtype'].value_counts().idxmax()
            if mask_b.sum() > 0 else 'Mono')

    smoothed_sc = {}
    for tf in selected_tfs:
        vals = []
        for i in range(len(pt_bins_sc) - 1):
            mask_b = ((df_scenic['dpt_pseudotime'] >= pt_bins_sc[i]) &
                      (df_scenic['dpt_pseudotime'] <  pt_bins_sc[i + 1]))
            vals.append(df_scenic.loc[mask_b, tf].mean()
                        if mask_b.sum() > 0 else np.nan)
        vals = pd.Series(vals).interpolate().values
        smoothed_sc[tf] = gaussian_filter1d(vals, sigma=sigma_sc)

    # ── 图 4：TF 动态折线图 ────────────────────────────────────────────────────
    if need_redraw('scenic_tf_dynamics'):
        fig4 = plt.figure(figsize=(12, 9))
        fig4.patch.set_facecolor('white')
        gs4 = gridspec.GridSpec(2, 1, height_ratios=[0.06, 1], hspace=0.08)
        ax_band4 = fig4.add_subplot(gs4[0])
        ax_main4 = fig4.add_subplot(gs4[1])

        bin_w4 = pt_ctrs_sc[1] - pt_ctrs_sc[0]
        for i, st in enumerate(bin_subtypes_sc):
            ax_band4.axvspan(pt_ctrs_sc[i] - bin_w4/2,
                             pt_ctrs_sc[i] + bin_w4/2,
                             color=SUBTYPE_COLORS_SC.get(st, '#CCCCCC'),
                             alpha=0.85, linewidth=0)
        for st in ['Mono', 'IFN_TAM', 'LA_TAM']:
            st_pts = df_scenic[df_scenic['mac_subtype'] == st]['dpt_pseudotime']
            if len(st_pts) == 0: continue
            ax_band4.text(st_pts.median(), 0.5, st, ha='center', va='center',
                          fontsize=8.5, fontweight='bold', color='white',
                          transform=ax_band4.get_xaxis_transform())
        ax_band4.set_xlim(pt_ctrs_sc[0], pt_ctrs_sc[-1])
        ax_band4.set_xticks([]); ax_band4.set_yticks([])
        ax_band4.set_ylabel('Subtype', fontsize=9, rotation=0,
                            ha='right', va='center', labelpad=40)
        for sp4 in ax_band4.spines.values(): sp4.set_visible(False)

        for idx, tf in enumerate(selected_tfs):
            vals      = smoothed_sc[tf]
            vals_norm = (vals - vals.min()) / (vals.max() - vals.min() + 1e-9)
            ax_main4.plot(pt_ctrs_sc, vals_norm,
                          color=TF_PALETTE[idx % len(TF_PALETTE)],
                          linewidth=2.2, label=tf.replace('(+)', ''), alpha=0.9)

        for st in ['Mono', 'IFN_TAM', 'LA_TAM']:
            st_pts = df_scenic[df_scenic['mac_subtype'] == st]['dpt_pseudotime'].values
            if len(st_pts) == 0: continue
            ax_main4.axvspan(st_pts.min(), st_pts.max(),
                             alpha=0.06, color=SUBTYPE_COLORS_SC[st], zorder=0)
            ax_main4.axvline(np.median(st_pts), color=SUBTYPE_COLORS_SC[st],
                             ls='--', lw=1.0, alpha=0.6, zorder=1)

        ax_main4.set_xlabel('DPT Pseudotime', fontsize=12)
        ax_main4.set_ylabel('Normalized TF Activity (AUC)', fontsize=12)
        ax_main4.set_xlim(pt_ctrs_sc[0], pt_ctrs_sc[-1])
        ax_main4.set_ylim(-0.05, 1.15)
        ax_main4.legend(title='Transcription Factor', loc='upper left',
                        fontsize=9, title_fontsize=9.5, frameon=True,
                        framealpha=0.9, ncol=2, bbox_to_anchor=(0.01, 0.99))
        sns.despine(ax=ax_main4)
        fig4.suptitle(
            'TF Regulon Activity Along Mono → IFN_TAM → LA_TAM Trajectory\n'
            '(pySCENIC AUCell; top 10 IFN_TAM-enriched TFs, padj<0.05)',
            fontsize=13, fontweight='bold', y=1.01)
        plt.tight_layout()
        save_fig(fig4, 'scenic_tf_dynamics')

    # ── 图 5：TF 热图 ──────────────────────────────────────────────────────────
    if need_redraw('scenic_tf_heatmap'):
        fig5, ax5 = plt.subplots(figsize=(max(10, len(sig_tfs_sc) * 0.32), 3.8))
        fig5.patch.set_facecolor('white')

        mean_auc5 = (auc_sub.join(meta_sub)
                     .groupby('mac_subtype')[sig_tfs_sc['TF'].tolist()]
                     .mean()
                     .loc[['Mono', 'IFN_TAM', 'LA_TAM']])
        mean_auc5.columns = [c.replace('(+)', '') for c in mean_auc5.columns]
        mean_auc5_z = (mean_auc5.apply(zscore, axis=0)
                       .T.sort_values('IFN_TAM', ascending=False).T)

        sns.heatmap(mean_auc5_z, ax=ax5, cmap='RdBu_r', center=0,
                    vmin=-2, vmax=2, linewidths=0.4, linecolor='white',
                    cbar_kws={'label': 'Z-score (AUC)', 'shrink': 0.6})

        for tick, st in zip(ax5.get_yticklabels(), ['Mono', 'IFN_TAM', 'LA_TAM']):
            tick.set_color(SUBTYPE_COLORS_SC[st])
            tick.set_fontweight('bold'); tick.set_fontsize(10)
        ax5.set_xticklabels(ax5.get_xticklabels(),
                            rotation=45, ha='right', fontsize=8.5)
        ax5.set_xlabel(''); ax5.set_ylabel('')
        ax5.set_title(
            'IFN_TAM-enriched TF Regulon Activity  (pySCENIC AUCell, padj<0.05)',
            fontsize=11, fontweight='bold', pad=10)
        plt.tight_layout()
        save_fig(fig5, 'scenic_tf_heatmap')

else:
    print(f"  SKIP Fig 4/5: SCENIC files not found at {SCENIC_DIR}")

# =============================================================================
#  图 6/7：Receptor 基因伪时间热图 + 折线图
# =============================================================================
volcano_path = os.path.join(CELLCOMM_DIR,
    'volcano_data_malignant_vs_normal_monomac.csv')

if os.path.exists(volcano_path):
    print("\n[Fig 6/7] Receptor pseudotime plots...")

    lr_df  = pd.read_csv(volcano_path)
    mal_lr = lr_df[lr_df['category'] == 'Malignant-enriched'].copy()

    def extract_genes(s):
        return [p.strip() for p in str(s).split('_')
                if p.strip() and p.strip() != 'nan']

    receptor_genes_all = []
    for rec in mal_lr['receptor_complex']:
        receptor_genes_all.extend(extract_genes(rec))
    receptor_genes = list(dict.fromkeys(receptor_genes_all))
    available_rec  = [g for g in receptor_genes if g in adata_path.var_names]
    print(f"  Receptor genes available: {len(available_rec)}")

    _X6 = adata_path[:, available_rec].layers['log1p_norm']
    if sps.issparse(_X6): _X6 = _X6.toarray()
    expr_path6 = pd.DataFrame(_X6, index=adata_path.obs_names,
                               columns=available_rec)

    corr_results6 = []
    for gene in available_rec:
        r, p = spearmanr(pt_sorted, expr_path6[gene].values)
        corr_results6.append({'gene': gene, 'spearman_r': r, 'pval': p})
    corr_res6 = pd.DataFrame(corr_results6)
    _, padj6, _, _ = multipletests(corr_res6['pval'], method='fdr_bh')
    corr_res6['padj'] = padj6
    corr_dict6 = dict(zip(corr_res6['gene'], corr_res6['spearman_r']))
    padj_dict6 = dict(zip(corr_res6['gene'], corr_res6['padj']))

    pos_genes6 = corr_res6[(corr_res6['spearman_r'] > 0) &
                            (corr_res6['padj'] <= 0.05)].sort_values(
                                'spearman_r', ascending=False)
    print(f"  Positively correlated genes: {len(pos_genes6)}")

    # ── 图 6：热图 ────────────────────────────────────────────────────────────
    if need_redraw('receptor_pseudotime_heatmap_path'):
        n_bins6 = 100; sigma6 = 3
        bin_edges6 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins6 + 1)
        bin_ctrs6  = (bin_edges6[:-1] + bin_edges6[1:]) / 2

        expr_binned6  = np.zeros((len(available_rec), n_bins6))
        subtype_bins6 = []
        for b in range(n_bins6):
            mask_b = ((pt_sorted >= bin_edges6[b]) &
                      (pt_sorted <  bin_edges6[b + 1]))
            if mask_b.sum() > 0:
                expr_binned6[:, b] = expr_path6[available_rec].values[mask_b].mean(axis=0)
                subtype_bins6.append(
                    Counter(subtype_sorted[mask_b]).most_common(1)[0][0])
            else:
                expr_binned6[:, b] = 0
                subtype_bins6.append(subtype_bins6[-1] if subtype_bins6 else 'Mono')

        expr_smooth6 = np.array([gaussian_filter1d(expr_binned6[i], sigma6)
                                  for i in range(len(available_rec))])
        expr_z6 = np.apply_along_axis(
            lambda x: zscore(x) if x.std() > 0 else x, 1, expr_smooth6)

        peak_idx6   = np.argmax(expr_z6, axis=1)
        gene_order6 = np.argsort(peak_idx6)
        genes_sorted6  = [available_rec[i] for i in gene_order6]
        expr_z_sorted6 = expr_z6[gene_order6]

        def sig_marker6(gene):
            r = corr_dict6.get(gene, 0)
            q = padj_dict6.get(gene, 1)
            if q <= 0.05:
                return ('▲', '#D4603A') if r > 0 else ('▼', '#2E7DC9')
            return ('', '#888888')

        sub_num6  = np.array([[{'Mono': 0, 'IFN_TAM': 1, 'LA_TAM': 2}
                                .get(s, 3) for s in subtype_bins6]])
        cmap_sub6 = ListedColormap(['#E9C46A', '#E63946', '#F4A261', '#CCCCCC'])

        fig6 = plt.figure(figsize=(14, 16))
        fig6.patch.set_facecolor('white')
        gs6 = gridspec.GridSpec(3, 2, height_ratios=[0.04, 1, 0.02],
                                 width_ratios=[1, 0.03], hspace=0.02, wspace=0.02)
        ax_bar6  = fig6.add_subplot(gs6[0, 0])
        ax_heat6 = fig6.add_subplot(gs6[1, 0])
        ax_cbar6 = fig6.add_subplot(gs6[1, 1])

        ax_bar6.imshow(sub_num6, aspect='auto', cmap=cmap_sub6,
                       vmin=0, vmax=3, interpolation='nearest')
        ax_bar6.set_xticks([]); ax_bar6.set_yticks([0])
        ax_bar6.set_yticklabels(['Subtype'], fontsize=8)
        ax_bar6.tick_params(left=False)
        ax_bar6.legend(
            handles=[Patch(color='#E9C46A', label='Mono'),
                     Patch(color='#E63946', label='IFN_TAM'),
                     Patch(color='#F4A261', label='LA_TAM')],
            loc='upper right', bbox_to_anchor=(1.12, 2.5),
            fontsize=8, frameon=True, ncol=3)

        cmap_heat6 = LinearSegmentedColormap.from_list(
            'expr', ['#2166AC', '#F7F7F7', '#D6604D'])
        im6 = ax_heat6.imshow(expr_z_sorted6, aspect='auto', cmap=cmap_heat6,
                               vmin=-2, vmax=2, interpolation='nearest')

        ytick_labels6, ytick_colors6 = [], []
        for g in genes_sorted6:
            marker, color = sig_marker6(g)
            ytick_labels6.append(f'{g} {marker}' if marker else g)
            ytick_colors6.append(color)
        ax_heat6.set_yticks(range(len(genes_sorted6)))
        ax_heat6.set_yticklabels(ytick_labels6, fontsize=7)
        for tick, color in zip(ax_heat6.get_yticklabels(), ytick_colors6):
            tick.set_color(color)
            if color != '#888888': tick.set_fontweight('bold')

        xtick_pos6  = np.linspace(0, n_bins6 - 1, 6)
        xtick_vals6 = np.linspace(pt_sorted.min(), pt_sorted.max(), 6)
        ax_heat6.set_xticks(xtick_pos6)
        ax_heat6.set_xticklabels([f'{v:.2f}' for v in xtick_vals6], fontsize=8.5)
        ax_heat6.set_xlabel('Pseudotime  (Mono → IFN_TAM → LA_TAM)', fontsize=10)
        ax_heat6.set_ylabel('Receptor genes (sorted by peak expression)', fontsize=9)

        prev6 = subtype_bins6[0]
        for b, st in enumerate(subtype_bins6):
            if st != prev6:
                ax_heat6.axvline(b - 0.5, color='white', lw=1.5,
                                 ls='--', alpha=0.8)
                prev6 = st

        plt.colorbar(im6, cax=ax_cbar6, label='Z-score')
        ax_cbar6.tick_params(labelsize=7)
        ax_heat6.text(1.08, -0.02,
                      '▲ pos. corr.  ▼ neg. corr.  (padj≤0.05)',
                      transform=ax_heat6.transAxes, fontsize=7.5,
                      ha='right', va='top', color='#555555')
        fig6.suptitle(
            f'Receptor Gene Expression Along Pseudotime\n'
            f'Mono → IFN_TAM → LA_TAM  ({len(available_rec)} receptor genes)',
            fontsize=12, fontweight='bold', y=1.005)
        plt.tight_layout()
        save_fig(fig6, 'receptor_pseudotime_heatmap_path')

    # ── 图 7：折线图（top 4 正相关基因）────────────────────────────────────────
    if need_redraw('receptor_pseudotime_lineplots_path'):
        plot_genes7 = pos_genes6['gene'].tolist()[:4]
        print(f"  Top positive genes: {plot_genes7}")

        mono_end7   = adata_path.obs[
            adata_path.obs['mac_subtype'] == 'Mono']['dpt_pseudotime'].max()
        ifntam_end7 = adata_path.obs[
            adata_path.obs['mac_subtype'] == 'IFN_TAM']['dpt_pseudotime'].max()

        n_bins7    = 80; sigma7 = 3
        bin_edges7 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins7 + 1)
        bin_ctrs7  = (bin_edges7[:-1] + bin_edges7[1:]) / 2

        fig7, axes7 = plt.subplots(2, 2, figsize=(13, 9))
        fig7.patch.set_facecolor('white')

        for ax7, gene in zip(axes7.flatten(), plot_genes7):
            expr7 = expr_path6[gene].values
            bm7 = np.zeros(n_bins7); bs7 = np.zeros(n_bins7)
            for b in range(n_bins7):
                mask_b = ((pt_sorted >= bin_edges7[b]) &
                          (pt_sorted <  bin_edges7[b + 1]))
                if mask_b.sum() > 0:
                    vals7  = expr7[mask_b]
                    bm7[b] = vals7.mean()
                    bs7[b] = vals7.std() / np.sqrt(len(vals7)) if len(vals7) > 1 else 0
            smooth7 = gaussian_filter1d(bm7, sigma=sigma7)
            sem7    = gaussian_filter1d(bs7, sigma=sigma7)

            ax7.axvspan(pt_sorted.min(), mono_end7,
                        alpha=0.12, color=PATH_COLORS['Mono'],    zorder=0)
            ax7.axvspan(mono_end7, ifntam_end7,
                        alpha=0.12, color=PATH_COLORS['IFN_TAM'], zorder=0)
            ax7.axvspan(ifntam_end7, pt_sorted.max(),
                        alpha=0.12, color=PATH_COLORS['LA_TAM'],  zorder=0)
            ax7.axvline(mono_end7,   color='#888888', lw=0.8, ls='--', zorder=1)
            ax7.axvline(ifntam_end7, color='#888888', lw=0.8, ls='--', zorder=1)

            for st7, sc_c7 in PATH_COLORS.items():
                mask_st7 = subtype_sorted == st7
                ax7.scatter(pt_sorted[mask_st7], expr7[mask_st7],
                            c=sc_c7, s=8, alpha=0.35, zorder=2, linewidths=0)

            ax7.plot(bin_ctrs7, smooth7, color='#333333', lw=2.2, zorder=4)
            ax7.fill_between(bin_ctrs7, smooth7 - sem7, smooth7 + sem7,
                             color='#333333', alpha=0.15, zorder=3)

            r_val7 = corr_dict6[gene]
            p_val7 = padj_dict6[gene]
            p_str7 = f'{p_val7:.2e}' if p_val7 >= 1e-10 else '<1e-10'
            ax7.text(0.97, 0.95,
                     f'r = {r_val7:.3f}\npadj = {p_str7}',
                     transform=ax7.transAxes, ha='right', va='top',
                     fontsize=9, color='#D4603A', fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.3', fc='white',
                               ec='#D4603A', alpha=0.85))
            ax7.set_title(gene, fontsize=13, fontweight='bold', color='#D4603A')
            ax7.set_xlabel('Pseudotime', fontsize=9)
            ax7.set_ylabel('log-normalized expression', fontsize=9)
            ax7.spines[['top', 'right']].set_visible(False)
            ax7.set_xlim(pt_sorted.min() - 0.005, pt_sorted.max() + 0.005)

        for ax7, gene in zip(axes7.flatten(), plot_genes7):
            ymax7 = ax7.get_ylim()[1]
            for st7, x7 in [
                ('Mono',    (pt_sorted.min() + mono_end7)   / 2),
                ('IFN_TAM', (mono_end7 + ifntam_end7)       / 2),
                ('LA_TAM',  (ifntam_end7 + pt_sorted.max()) / 2),
            ]:
                ax7.text(x7, ymax7 * 0.97, st7, ha='center', va='top',
                         fontsize=8, color=PATH_COLORS[st7], fontweight='bold')

        legend7 = (
            [Patch(color=c, alpha=0.5, label=s) for s, c in PATH_COLORS.items()] +
            [Line2D([0], [0], color='#333333', lw=2, label='Smoothed mean')]
        )
        fig7.legend(handles=legend7, loc='lower center', ncol=4,
                    fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))
        fig7.suptitle(
            'Receptor Genes Positively Correlated with Pseudotime\n'
            'Mono → IFN_TAM → LA_TAM  (Spearman, BH-corrected, padj≤0.05)',
            fontsize=12, fontweight='bold', y=1.01)
        plt.tight_layout()
        save_fig(fig7, 'receptor_pseudotime_lineplots_path')

else:
    print(f"  SKIP Fig 6/7: {volcano_path} not found")

# =============================================================================
#  图 8/9：IFN_TAM vs LA_TAM DEG 火山图 + 富集分析条形图
# =============================================================================
print("\n[Fig 8/9] DEG volcano + enrichment...")

adata_sub8 = adata_mac[adata_mac.obs['mac_subtype'].isin(
    ['IFN_TAM', 'LA_TAM'])].copy()
expr8 = (adata_sub8.layers['log1p_norm'].toarray()
         if sps.issparse(adata_sub8.layers['log1p_norm'])
         else adata_sub8.layers['log1p_norm'])

ifntam_mask8 = (adata_sub8.obs['mac_subtype'] == 'IFN_TAM').values
latam_mask8  = (adata_sub8.obs['mac_subtype'] == 'LA_TAM').values
expr_i8 = expr8[ifntam_mask8]
expr_l8 = expr8[latam_mask8]
genes8  = adata_sub8.var_names.tolist()
print(f"  Testing {len(genes8)} genes")

pvals8, log2fcs8 = [], []
for i in range(len(genes8)):
    g_i = expr_i8[:, i]; g_l = expr_l8[:, i]
    if g_i.max() == 0 and g_l.max() == 0:
        pvals8.append(1.0); log2fcs8.append(0.0)
    else:
        _, p = mannwhitneyu(g_i, g_l, alternative='two-sided')
        log2fcs8.append(np.log2((g_i.mean() + 1e-9) / (g_l.mean() + 1e-9)))
        pvals8.append(p)

_, padj8, _, _ = multipletests(pvals8, method='fdr_bh')
deg_df8 = pd.DataFrame({'gene': genes8, 'log2FC': log2fcs8,
                         'pval': pvals8, 'padj': padj8})

sig_df8  = deg_df8[(deg_df8['padj'] <= 0.05) &
                   (deg_df8['log2FC'].abs() >= 0.5)].copy()
up_df8   = sig_df8[sig_df8['log2FC'] > 0].sort_values('log2FC', ascending=False)
down_df8 = sig_df8[sig_df8['log2FC'] < 0].sort_values('log2FC')
print(f"  Up in IFN_TAM: {len(up_df8)}, Up in LA_TAM: {len(down_df8)}")

os.makedirs(DEG_DIR, exist_ok=True)
sig_df8.sort_values('log2FC', ascending=False).to_csv(
    os.path.join(DEG_DIR, 'DEG_IFN_TAM_vs_LA_TAM_full.csv'), index=False)

# ── 图 8：火山图 ──────────────────────────────────────────────────────────────
if need_redraw('volcano_IFN_TAM_vs_LA_TAM', DEG_DIR):
    try:
        from adjustText import adjust_text
        HAS_ADJUSTTEXT = True
    except ImportError:
        HAS_ADJUSTTEXT = False
        print("  WARNING: adjustText not installed")

    plot_df8 = deg_df8.copy()
    plot_df8['-log10padj']     = -np.log10(plot_df8['padj'].clip(lower=1e-30))
    plot_df8['log2FC_clipped'] = plot_df8['log2FC'].clip(-8, 8)

    def classify8(row):
        if row['padj'] <= 0.05 and row['log2FC'] >= 0.5:  return 'Up'
        if row['padj'] <= 0.05 and row['log2FC'] <= -0.5: return 'Down'
        return 'NS'

    plot_df8['class'] = plot_df8.apply(classify8, axis=1)
    colors8 = {'Up': '#D4603A', 'Down': '#2E7DC9', 'NS': '#CCCCCC'}
    sizes8  = {'Up': 18, 'Down': 18, 'NS': 6}
    alphas8 = {'Up': 0.75, 'Down': 0.75, 'NS': 0.25}

    top_up8   = plot_df8[plot_df8['class'] == 'Up'].nlargest(20, 'log2FC')
    top_down8 = plot_df8[plot_df8['class'] == 'Down'].nsmallest(15, 'log2FC')
    label_genes8 = pd.concat([top_up8, top_down8])

    fig8, ax8 = plt.subplots(figsize=(9, 7))
    fig8.patch.set_facecolor('white'); ax8.set_facecolor('white')

    for cls in ['NS', 'Up', 'Down']:
        sub = plot_df8[plot_df8['class'] == cls]
        ax8.scatter(sub['log2FC_clipped'], sub['-log10padj'],
                    c=colors8[cls], s=sizes8[cls], alpha=alphas8[cls],
                    linewidths=0, zorder=2 if cls != 'NS' else 1)

    ax8.axhline(-np.log10(0.05), color='#888888', lw=0.8, ls='--', zorder=0)
    ax8.axvline( 0.5, color='#888888', lw=0.8, ls='--', zorder=0)
    ax8.axvline(-0.5, color='#888888', lw=0.8, ls='--', zorder=0)

    texts8 = []
    for _, row in label_genes8.iterrows():
        t = ax8.text(row['log2FC_clipped'], row['-log10padj'], row['gene'],
                     fontsize=7.5, color=colors8[row['class']],
                     fontweight='bold', zorder=5)
        texts8.append(t)

    if HAS_ADJUSTTEXT:
        adjust_text(texts8, ax=ax8,
                    arrowprops=dict(arrowstyle='-', color='#999999', lw=0.5),
                    expand=(1.2, 1.4), force_text=(0.3, 0.5))

    n_up8   = (plot_df8['class'] == 'Up').sum()
    n_down8 = (plot_df8['class'] == 'Down').sum()
    ax8.text(0.97, 0.97, f'Up in IFN_TAM: {n_up8}',
             transform=ax8.transAxes, ha='right', va='top',
             fontsize=10, color=colors8['Up'], fontweight='bold')
    ax8.text(0.03, 0.97, f'Up in LA_TAM: {n_down8}',
             transform=ax8.transAxes, ha='left', va='top',
             fontsize=10, color=colors8['Down'], fontweight='bold')

    ax8.set_xlabel('log₂ Fold Change  (IFN_TAM vs LA_TAM)', fontsize=11)
    ax8.set_ylabel('−log₁₀(adjusted p-value)', fontsize=11)
    ax8.set_title('IFN_TAM vs LA_TAM — Differential Expression\n'
                  '(Wilcoxon, BH-corrected; padj≤0.05, |log₂FC|≥0.5)',
                  fontsize=11, fontweight='bold')
    legend8 = [
        mpatches.Patch(color=colors8['Up'],   label=f'Up in IFN_TAM ({n_up8})'),
        mpatches.Patch(color=colors8['Down'], label=f'Up in LA_TAM ({n_down8})'),
        mpatches.Patch(color=colors8['NS'],   label='Not significant'),
    ]
    ax8.legend(handles=legend8, fontsize=9, frameon=True,
               loc='lower right', framealpha=0.9)
    ax8.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    save_fig(fig8, 'volcano_IFN_TAM_vs_LA_TAM', DEG_DIR)

# ── 图 9：富集分析条形图 ──────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
# Fig9 前置参数：Pathway Enrichment Barplot — IFN_TAM vs LA_TAM
# ══════════════════════════════════════════════════════════════════════════════

P09 = {

    # ── 数据输入 ──────────────────────────────────────────────────────────────
    "deg_dir"              : DEG_DIR,             # 富集结果 CSV 的读写目录

    # ── Enrichr 基因集 ────────────────────────────────────────────────────────
    "gs_kegg"              : 'KEGG_2021_Human',
    "gs_gobp"              : 'GO_Biological_Process_2021',
    "organism"             : 'human',
    "padj_cutoff"          : 0.05,                # 显著性阈值
    "padj_clip"            : 1e-20,               # -log10 计算时的 p 值下限（防止 inf）

    # ── 每个面板展示的最多条目数 ──────────────────────────────────────────────
    "top_n"                : 10,

    # ── Term 名称截断长度 ─────────────────────────────────────────────────────
    "term_max_len"         : 32,                  # 超过此长度截断并加 '...'

    # ── 颜色映射（上调 / 下调各一套）─────────────────────────────────────────
    "cmap_up"              : plt.cm.YlGnBu,       # IFN_TAM 上调（蓝绿色系）
    "cmap_dn"              : plt.cm.YlOrRd,       # LA_TAM  上调（橙红色系）

    # ── 条形图样式 ────────────────────────────────────────────────────────────
    "bar_height"           : 0.75,                # 条形高度
    "bar_edgecolor"        : 'white',
    "bar_edgelw"           : 0.5,

    # ── 条形内 n= 标注 ────────────────────────────────────────────────────────
    "n_label_x_ratio"      : 0.97,                # 标注 x 位置（条宽 × ratio）
    "n_label_fontsize"     : 16,
    "n_label_color"        : 'white',
    "n_label_fontweight"   : 'bold',

    # ── 轴标签与刻度 ──────────────────────────────────────────────────────────
    "xlabel_fontsize"      : 14,
    "yticklabel_fontsize"  : 13.5,
    "title_fontsize"       : 18,
    "title_pad"            : 8,
    "no_sig_fontsize"      : 16,                  # "No significant terms" 提示字号
    "no_sig_color"         : '#888888',

    # ── 色条（colorbar）────────────────────────────────────────────────────────
    "cbar_shrink"          : 0.4,
    "cbar_pad"             : 0.02,
    "cbar_aspect"          : 15,
    "cbar_label_fontsize"  : 12,
    "cbar_tick_fontsize"   : 10,

    # ── 图像尺寸 ──────────────────────────────────────────────────────────────
    "figsize"              : (15, 9),
    "fig_facecolor"        : 'white',

    # ── 四个面板的标题 ────────────────────────────────────────────────────────
    "panel_titles"         : [
        'KEGG Pathways\n(IFN_TAM upregulated)',
        'GO Biological Process\n(IFN_TAM upregulated)',
        'KEGG Pathways\n(LA_TAM upregulated)',
        'GO Biological Process\n(LA_TAM upregulated)',
    ],

    # ── 总标题 ────────────────────────────────────────────────────────────────
    "suptitle_fontsize"    : 11,
    "suptitle_y"           : 1.01,

    # ── 输出文件名 ────────────────────────────────────────────────────────────
    "output_name"          : 'enrichment_barplot_IFN_TAM_vs_LA_TAM',
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig9 辅助函数
# ══════════════════════════════════════════════════════════════════════════════

def load_or_run_enrichr9(gene_list, gs_name, short_name, label, out_dir,
                          padj_cutoff, organism):
    """读取缓存 CSV，若不存在则调用 Enrichr 并保存结果。"""
    csv_path = os.path.join(out_dir, f'enrichment_{label}_{short_name}.csv')
    if os.path.exists(csv_path):
        df  = pd.read_csv(csv_path)
        sig = df[df['Adjusted P-value'] <= padj_cutoff].sort_values('Adjusted P-value')
        print(f"    Loaded {os.path.basename(csv_path)}: {len(sig)} sig terms")
        return sig
    try:
        import gseapy as gp
        print(f"    Running Enrichr ({label}/{gs_name})...")
        enr = gp.enrichr(gene_list=gene_list, gene_sets=gs_name,
                         organism=organism, outdir=None, verbose=False)
        sig = enr.results[enr.results['Adjusted P-value'] <= padj_cutoff].sort_values(
            'Adjusted P-value')
        sig.to_csv(csv_path, index=False)
        print(f"    Done: {len(sig)} sig terms")
        return sig
    except Exception as e:
        print(f"    Enrichr failed: {e}")
        return pd.DataFrame()


def clean_term9(term, max_len):
    term = re.sub(r'\s*\(GO:\d+\)', '', term)
    return term[:max_len - 3] + '...' if len(term) > max_len else term


def parse_overlap9(s):
    try:    return int(str(s).split('/')[0])
    except: return 0


def prep9(df, n, max_len, padj_clip):
    if df is None or len(df) == 0:
        return pd.DataFrame()
    df = df.head(n).copy()
    df['Term_clean']   = df['Term'].apply(lambda t: clean_term9(t, max_len))
    df['-log10padj']   = -np.log10(df['Adjusted P-value'].clip(lower=padj_clip))
    df['n_genes']      = df['Overlap'].apply(parse_overlap9)
    return df.sort_values('-log10padj').reset_index(drop=True)

# ══════════════════════════════════════════════════════════════════════════════
# Fig9 绘图
# ══════════════════════════════════════════════════════════════════════════════

print("\n[Fig 9] Drawing pathway enrichment barplot (IFN_TAM vs LA_TAM)...")

# ── 运行 / 加载富集分析 ───────────────────────────────────────────────────────
_kw = dict(out_dir      = P09["deg_dir"],
           padj_cutoff  = P09["padj_cutoff"],
           organism     = P09["organism"])

up_kegg9 = load_or_run_enrichr9(up_df8['gene'].tolist(),
                                 P09["gs_kegg"], 'KEGG', 'up', **_kw)
up_gobp9 = load_or_run_enrichr9(up_df8['gene'].tolist(),
                                 P09["gs_gobp"], 'GOBP', 'up', **_kw)
dn_kegg9 = load_or_run_enrichr9(down_df8['gene'].tolist(),
                                 P09["gs_kegg"], 'KEGG', 'dn', **_kw)
dn_gobp9 = load_or_run_enrichr9(down_df8['gene'].tolist(),
                                 P09["gs_gobp"], 'GOBP', 'dn', **_kw)

# ── 数据预处理 ────────────────────────────────────────────────────────────────
_prep_kw = dict(n=P09["top_n"], max_len=P09["term_max_len"],
                padj_clip=P09["padj_clip"])

top_up_kegg9 = prep9(up_kegg9, **_prep_kw)
top_up_gobp9 = prep9(up_gobp9, **_prep_kw)
top_dn_kegg9 = prep9(dn_kegg9, **_prep_kw)
top_dn_gobp9 = prep9(dn_gobp9, **_prep_kw)

# ── 画布 ──────────────────────────────────────────────────────────────────────
fig9, axes9 = plt.subplots(2, 2, figsize=P09["figsize"])
fig9.patch.set_facecolor(P09["fig_facecolor"])

panel_configs9 = [
    (axes9[0, 0], top_up_kegg9, P09["panel_titles"][0], P09["cmap_up"]),
    (axes9[0, 1], top_up_gobp9, P09["panel_titles"][1], P09["cmap_up"]),
    (axes9[1, 0], top_dn_kegg9, P09["panel_titles"][2], P09["cmap_dn"]),
    (axes9[1, 1], top_dn_gobp9, P09["panel_titles"][3], P09["cmap_dn"]),
]

for ax9, df9, title9, cmap9 in panel_configs9:

    # 无显著条目时显示占位文字
    if len(df9) == 0:
        ax9.text(0.5, 0.5, 'No significant terms',
                 ha='center', va='center',
                 fontsize=P09["no_sig_fontsize"],
                 transform=ax9.transAxes,
                 color=P09["no_sig_color"])
        ax9.set_title(title9,
                      fontsize=P09["title_fontsize"], fontweight='bold')
        ax9.axis('off')
        continue

    y_pos9   = np.arange(len(df9))
    norm9    = Normalize(vmin=df9['-log10padj'].min(),
                         vmax=df9['-log10padj'].max())
    bar_clrs = [cmap9(norm9(v)) for v in df9['-log10padj']]

    bars9 = ax9.barh(y_pos9, df9['-log10padj'],
                     color=bar_clrs,
                     edgecolor=P09["bar_edgecolor"],
                     linewidth=P09["bar_edgelw"],
                     height=P09["bar_height"])

    # 条形内 n= 标注
    for bar9, cnt9 in zip(bars9, df9['n_genes']):
        w9 = bar9.get_width()
        ax9.text(w9 * P09["n_label_x_ratio"],
                 bar9.get_y() + bar9.get_height() / 2,
                 f'n={cnt9}',
                 ha='right', va='center',
                 fontsize=P09["n_label_fontsize"],
                 color=P09["n_label_color"],
                 fontweight=P09["n_label_fontweight"])

    ax9.set_yticks(y_pos9)
    ax9.set_yticklabels(df9['Term_clean'], fontsize=P09["yticklabel_fontsize"])
    ax9.set_xlabel('−log10(adjusted p-value)', fontsize=P09["xlabel_fontsize"])
    ax9.set_title(title9,
                  fontsize=P09["title_fontsize"], fontweight='bold',
                  pad=P09["title_pad"])
    ax9.spines[['top', 'right']].set_visible(False)

    # 色条
    sm9 = plt.cm.ScalarMappable(cmap=cmap9, norm=norm9)
    sm9.set_array([])
    cbar9 = plt.colorbar(sm9, ax=ax9,
                          shrink=P09["cbar_shrink"],
                          pad=P09["cbar_pad"],
                          aspect=P09["cbar_aspect"])
    cbar9.set_label('−log10(padj)', fontsize=P09["cbar_label_fontsize"])
    cbar9.ax.tick_params(labelsize=P09["cbar_tick_fontsize"])

# ── 总标题 ────────────────────────────────────────────────────────────────────
fig9.suptitle(
    f'',
    fontsize=P09["suptitle_fontsize"], fontweight='bold',
    y=P09["suptitle_y"],
)

plt.tight_layout()
save_fig(fig9, P09["output_name"], P09["deg_dir"])
plt.close(fig9)
print("  ✓ enrichment_barplot_IFN_TAM_vs_LA_TAM saved")

# =============================================================================
#  图 10/11：CellRank 基因接力热图 + fate trend 折线图
# =============================================================================
print("\n[Fig 10/11] CellRank plots...")

try:
    import cellrank as cr
    CELLRANK_OK = True
    print(f"  CellRank version: {cr.__version__}")
except ImportError:
    CELLRANK_OK = False
    print("  WARNING: cellrank not installed → pip install cellrank")

if CELLRANK_OK and os.path.exists(volcano_path):

    # ── CellRank 核心计算（健壮版）────────────────────────────────────────────────
    # 替换原来的 sc.pp.neighbors 调用
    sc.pp.neighbors(
        adata_mac,
        n_neighbors=15,
        n_pcs=20,
        method='umap',          # 明确指定，不让 scanpy 自动选
        random_state=0,         # 锁定种子
        metric='euclidean'
    )
    params = adata_mac.uns['neighbors']['params']
    print("neighbors params:", params)
    # scanpy 1.9.8 + annoy 已安装时，method 应显示 'gauss' 或 params 里含 annoy 相关信息


    adata_path.obs['mac_subtype'] = (
        adata_path.obs['mac_subtype'].astype('category'))
    
    pk = cr.kernels.PseudotimeKernel(adata_path, time_key='dpt_pseudotime')
    pk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)

    g2 = cr.estimators.GPCCA(pk)

    # ── Step 1：先尝试 n_states=3，若找不到 IFN_TAM 则升至 4/5 ──────────────────
    for n_st in [3, 4, 5]:
        g2.fit(n_states=n_st, cluster_key='mac_subtype')
        actual_names = list(g2.macrostates_memberships.names)
        print(f"  n_states={n_st} → macrostates: {actual_names}")

    # 检查三个目标亚型是否都能找到匹配
        def find_state(target, names):
            """大小写不敏感的模糊匹配"""
            t = target.lower()
            for n in names:
                if t in n.lower():
                    return n
            return None

        matched = {st: find_state(st, actual_names)
                   for st in ['Mono', 'IFN_TAM', 'LA_TAM']}
        print(f"  Matched: {matched}")

        if all(v is not None for v in matched.values()):
            break
        print(f"  Not all targets found, trying n_states={n_st + 1}...")

    # ── Step 2：用实际名称设置 terminal states ───────────────────────────────────
    terminal_actual = [v for v in matched.values() if v is not None]
    print(f"  Setting terminal states: {terminal_actual}")
    g2.set_terminal_states(terminal_actual, cluster_key='mac_subtype')
    g2.compute_fate_probabilities()

    # ── Step 3：构建 fate probability DataFrame，列名统一映射回原始亚型名 ─────────
    fp_arr  = np.array(g2.fate_probabilities.X)
    fp_cols = list(g2.fate_probabilities.names)          # 实际列名，如 Mono_1
    fp_df   = pd.DataFrame(fp_arr,
                            index=adata_path.obs_names,
                            columns=fp_cols)

    # 反向映射：实际名 → 原始亚型名
    reverse_map = {v: k for k, v in matched.items() if v is not None}
    fp_df = fp_df.rename(columns=reverse_map)            # 列名改为 Mono/IFN_TAM/LA_TAM

    # 若某亚型未匹配到，用 0 填充
    for st in ['Mono', 'IFN_TAM', 'LA_TAM']:
        if st not in fp_df.columns:
            fp_df[st] = 0.0
            print(f"  WARNING: {st} not found in fate probabilities, filled with 0")

    adata_path.obs['fate_LA_TAM']  = fp_df['LA_TAM'].values
    adata_path.obs['fate_IFN_TAM'] = fp_df['IFN_TAM'].values
    adata_path.obs['fate_Mono']    = fp_df['Mono'].values
    adata_path.obsm['lineages_fwd'] = g2.fate_probabilities

    print("  Fate probabilities assigned:")
    print(adata_path.obs[['fate_Mono', 'fate_IFN_TAM', 'fate_LA_TAM']].describe().round(3))

    # ── Step 4：Lineage drivers（用实际列名）────────────────────────────────────
    la_actual  = matched.get('LA_TAM')
    ifn_actual = matched.get('IFN_TAM')

    drivers_la = drivers_ifn = None

    if la_actual:
        try:
            drivers_la = g2.compute_lineage_drivers(
                lineages=la_actual, use_raw=False,
                layer='log1p_norm', return_drivers=True)
            # 统一列名
            drivers_la = drivers_la.rename(
                columns={f'{la_actual}_corr': 'LA_TAM_corr',
                         f'{la_actual}_pval': 'LA_TAM_pval'})
            print(f"  LA_TAM drivers computed: {drivers_la.shape}")
        except Exception as e:
            print(f"  WARNING: LA_TAM drivers failed: {e}")

    if ifn_actual:
        try:
            drivers_ifn = g2.compute_lineage_drivers(
                lineages=ifn_actual, use_raw=False,
                layer='log1p_norm', return_drivers=True)
            drivers_ifn = drivers_ifn.rename(
                columns={f'{ifn_actual}_corr': 'IFN_TAM_corr',
                         f'{ifn_actual}_pval': 'IFN_TAM_pval'})
            print(f"  IFN_TAM drivers computed: {drivers_ifn.shape}")
        except Exception as e:
            print(f"  WARNING: IFN_TAM drivers failed: {e}")

    # ── Step 5：筛选 receptor driver 基因（兼容 drivers 为 None 的情况）──────────
    receptor_drivers_la  = ([g for g in available_rec if g in drivers_la.index]
                             if drivers_la is not None else [])
    receptor_drivers_ifn = ([g for g in available_rec if g in drivers_ifn.index]
                             if drivers_ifn is not None else [])

    top_la = (drivers_la.loc[receptor_drivers_la, 'LA_TAM_corr']
              .sort_values(ascending=False).head(10)
              if drivers_la is not None and receptor_drivers_la
              else pd.Series(dtype=float))

    top_ifn = (drivers_ifn.loc[receptor_drivers_ifn, 'IFN_TAM_corr']
               .sort_values(ascending=False).head(10)
               if drivers_ifn is not None and receptor_drivers_ifn
               else pd.Series(dtype=float))

    # 若 receptor driver 不足，用伪时间相关性 top 基因补充
    if len(top_la) + len(top_ifn) < 4:
        print("  Falling back to pseudotime-correlated receptor genes for relay plot")
        fallback = pos_genes6['gene'].tolist()[:10]
        top_la  = pd.Series(
            {g: corr_dict6[g] for g in fallback[:5] if g in corr_dict6},
            dtype=float)
        top_ifn = pd.Series(
            {g: corr_dict6[g] for g in fallback[5:10] if g in corr_dict6},
            dtype=float)

    relay_genes = list(dict.fromkeys(
        top_la.index.tolist() + top_ifn.index.tolist()))[:20]
    print(f"  Final relay genes ({len(relay_genes)}): {relay_genes}")


    # ── 图 10：基因接力热图 ───────────────────────────────────────────────────
    if need_redraw('cellrank_gene_relay_heatmap'):
        n_bins10  = 100; sigma10 = 3
        bin_edges10 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins10 + 1)
        bin_ctrs10  = (bin_edges10[:-1] + bin_edges10[1:]) / 2

        relay_avail = [g for g in relay_genes if g in adata_path.var_names]
        _X10 = adata_path[:, relay_avail].layers['log1p_norm']
        if sps.issparse(_X10): _X10 = _X10.toarray()
        expr10 = pd.DataFrame(_X10, index=adata_path.obs_names,
                               columns=relay_avail)

        expr_binned10 = np.zeros((len(relay_avail), n_bins10))
        subtype_bins10 = []
        for b in range(n_bins10):
            mask_b = ((pt_sorted >= bin_edges10[b]) &
                      (pt_sorted <  bin_edges10[b + 1]))
            if mask_b.sum() > 0:
                expr_binned10[:, b] = expr10[relay_avail].values[
                    mask_b].mean(axis=0)
                subtype_bins10.append(
                    Counter(subtype_sorted[mask_b]).most_common(1)[0][0])
            else:
                expr_binned10[:, b] = (expr_binned10[:, b-1]
                                       if b > 0 else 0)
                subtype_bins10.append(
                    subtype_bins10[-1] if subtype_bins10 else 'Mono')

        expr_smooth10 = np.array([
            gaussian_filter1d(expr_binned10[i], sigma10)
            for i in range(len(relay_avail))])
        expr_z10 = np.apply_along_axis(
            lambda x: zscore(x) if x.std() > 0 else x, 1, expr_smooth10)

        peak_idx10   = np.argmax(expr_z10, axis=1)
        order10      = np.argsort(peak_idx10)
        genes_ord10  = [relay_avail[i] for i in order10]
        expr_z_ord10 = expr_z10[order10]

        # 标注 lineage 归属
        gene_lineage10 = {}
        for g in genes_ord10:
            in_la  = g in top_la.index
            in_ifn = g in top_ifn.index
            if in_la and in_ifn:
                gene_lineage10[g] = 'Both'
            elif in_la:
                gene_lineage10[g] = 'LA_TAM'
            elif in_ifn:
                gene_lineage10[g] = 'IFN_TAM'
            else:
                gene_lineage10[g] = 'Other'

        lineage_colors10 = {
            'LA_TAM':  '#2ECC71', 'IFN_TAM': '#E67E22',
            'Both':    '#9B59B6', 'Other':   '#AAAAAA',
        }

        sub_num10  = np.array([[{'Mono': 0, 'IFN_TAM': 1, 'LA_TAM': 2}
                                 .get(s, 3) for s in subtype_bins10]])
        cmap_sub10 = ListedColormap(['#E9C46A', '#E63946', '#F4A261', '#CCCCCC'])

        fig10 = plt.figure(figsize=(14, max(8, len(genes_ord10) * 0.55 + 2)))
        fig10.patch.set_facecolor('white')
        gs10 = gridspec.GridSpec(2, 2,
                                  height_ratios=[0.04, 1],
                                  width_ratios=[1, 0.03],
                                  hspace=0.02, wspace=0.02)
        ax_bar10  = fig10.add_subplot(gs10[0, 0])
        ax_heat10 = fig10.add_subplot(gs10[1, 0])
        ax_cbar10 = fig10.add_subplot(gs10[1, 1])

        ax_bar10.imshow(sub_num10, aspect='auto', cmap=cmap_sub10,
                        vmin=0, vmax=3, interpolation='nearest')
        ax_bar10.set_xticks([]); ax_bar10.set_yticks([0])
        ax_bar10.set_yticklabels(['Subtype'], fontsize=8)
        ax_bar10.tick_params(left=False)
        ax_bar10.legend(
            handles=[Patch(color='#E9C46A', label='Mono'),
                     Patch(color='#E63946', label='IFN_TAM'),
                     Patch(color='#F4A261', label='LA_TAM')],
            loc='upper right', bbox_to_anchor=(1.12, 2.8),
            fontsize=8, frameon=True, ncol=3)

        cmap_heat10 = LinearSegmentedColormap.from_list(
            'relay', ['#2166AC', '#F7F7F7', '#D6604D'])
        im10 = ax_heat10.imshow(expr_z_ord10, aspect='auto',
                                 cmap=cmap_heat10, vmin=-2, vmax=2,
                                 interpolation='nearest')

        ax_heat10.set_yticks(range(len(genes_ord10)))
        ax_heat10.set_yticklabels(genes_ord10, fontsize=8)
        for tick, gene in zip(ax_heat10.get_yticklabels(), genes_ord10):
            lin = gene_lineage10.get(gene, 'Other')
            tick.set_color(lineage_colors10[lin])
            if lin != 'Other': tick.set_fontweight('bold')

        # 右侧 lineage 标注条
        for idx, gene in enumerate(genes_ord10):
            lin = gene_lineage10.get(gene, 'Other')
            ax_heat10.annotate(
                '', xy=(n_bins10 + 1.5, idx),
                xytext=(n_bins10 + 0.5, idx),
                xycoords='data',
                arrowprops=dict(arrowstyle='-',
                                color=lineage_colors10[lin], lw=3))

        xtick_pos10  = np.linspace(0, n_bins10 - 1, 6)
        xtick_vals10 = np.linspace(pt_sorted.min(), pt_sorted.max(), 6)
        ax_heat10.set_xticks(xtick_pos10)
        ax_heat10.set_xticklabels(
            [f'{v:.2f}' for v in xtick_vals10], fontsize=8.5)
        ax_heat10.set_xlabel(
            'Pseudotime  (Mono → IFN_TAM / LA_TAM)', fontsize=10)
        ax_heat10.set_ylabel(
            'Receptor genes (CellRank lineage drivers)', fontsize=9)

        prev10 = subtype_bins10[0]
        for b, st in enumerate(subtype_bins10):
            if st != prev10:
                ax_heat10.axvline(b - 0.5, color='white',
                                   lw=1.5, ls='--', alpha=0.8)
                prev10 = st

        plt.colorbar(im10, cax=ax_cbar10, label='Z-score')
        ax_cbar10.tick_params(labelsize=7)

        # 图例
        fig10.legend(
            handles=[Patch(color=c, label=l)
                     for l, c in lineage_colors10.items()
                     if l != 'Other'],
            title='Lineage driver', title_fontsize=8.5,
            fontsize=8, loc='lower right',
            bbox_to_anchor=(0.88, 0.02), framealpha=0.9)

        fig10.suptitle(
            f'CellRank Lineage Driver Receptor Genes — Gene Relay Heatmap\n'
            f'Mono → IFN_TAM / LA_TAM  ({len(genes_ord10)} genes)',
            fontsize=12, fontweight='bold', y=1.005)
        plt.tight_layout()
        save_fig(fig10, 'cellrank_gene_relay_heatmap')

    # ── 图 11：Fate trend 折线图（top 4 receptor drivers）────────────────────
    if need_redraw('cellrank_fate_trend_plots'):
        # 选 top 4：LA_TAM + IFN_TAM 各取 top 2，去重
        trend_genes = list(dict.fromkeys(
            top_la.index[:2].tolist() +
            top_ifn.index[:2].tolist()))[:4]
        # 如不足 4 个则补充
        if len(trend_genes) < 4:
            extra = [g for g in relay_genes
                     if g not in trend_genes]
            trend_genes += extra[:4 - len(trend_genes)]
        print(f"  Fate trend genes: {trend_genes}")

        n_bins11  = 80; sigma11 = 3
        bin_edges11 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins11 + 1)
        bin_ctrs11  = (bin_edges11[:-1] + bin_edges11[1:]) / 2

        trend_avail = [g for g in trend_genes if g in adata_path.var_names]
        _X11 = adata_path[:, trend_avail].layers['log1p_norm']
        if sps.issparse(_X11): _X11 = _X11.toarray()
        expr11 = pd.DataFrame(_X11, index=adata_path.obs_names,
                               columns=trend_avail)

        fate_la_arr  = adata_path.obs['fate_LA_TAM'].values
        fate_ifn_arr = adata_path.obs['fate_IFN_TAM'].values

        fig11, axes11 = plt.subplots(2, 2, figsize=(14, 10))
        fig11.patch.set_facecolor('white')

        for ax11, gene in zip(axes11.flatten(), trend_avail):
            expr_g = expr11[gene].values

            # 平滑表达
            bm11 = np.zeros(n_bins11); bs11 = np.zeros(n_bins11)
            for b in range(n_bins11):
                mask_b = ((pt_sorted >= bin_edges11[b]) &
                          (pt_sorted <  bin_edges11[b + 1]))
                if mask_b.sum() > 0:
                    v = expr_g[mask_b]
                    bm11[b] = v.mean()
                    bs11[b] = (v.std() / np.sqrt(len(v))
                               if len(v) > 1 else 0)
            smooth11 = gaussian_filter1d(bm11, sigma=sigma11)
            sem11    = gaussian_filter1d(bs11, sigma=sigma11)

            # 平滑 fate probability
            fp_la_bin  = np.zeros(n_bins11)
            fp_ifn_bin = np.zeros(n_bins11)
            for b in range(n_bins11):
                mask_b = ((pt_sorted >= bin_edges11[b]) &
                          (pt_sorted <  bin_edges11[b + 1]))
                if mask_b.sum() > 0:
                    fp_la_bin[b]  = fate_la_arr[mask_b].mean()
                    fp_ifn_bin[b] = fate_ifn_arr[mask_b].mean()
            fp_la_sm  = gaussian_filter1d(fp_la_bin,  sigma=sigma11)
            fp_ifn_sm = gaussian_filter1d(fp_ifn_bin, sigma=sigma11)

            # 背景色带
            mono_end11   = adata_path.obs[
                adata_path.obs['mac_subtype'] == 'Mono'][
                'dpt_pseudotime'].max()
            ifntam_end11 = adata_path.obs[
                adata_path.obs['mac_subtype'] == 'IFN_TAM'][
                'dpt_pseudotime'].max()

            ax11.axvspan(pt_sorted.min(), mono_end11,
                         alpha=0.10, color=PATH_COLORS['Mono'],    zorder=0)
            ax11.axvspan(mono_end11, ifntam_end11,
                         alpha=0.10, color=PATH_COLORS['IFN_TAM'], zorder=0)
            ax11.axvspan(ifntam_end11, pt_sorted.max(),
                         alpha=0.10, color=PATH_COLORS['LA_TAM'],  zorder=0)
            ax11.axvline(mono_end11,   color='#AAAAAA',
                         lw=0.8, ls='--', zorder=1)
            ax11.axvline(ifntam_end11, color='#AAAAAA',
                         lw=0.8, ls='--', zorder=1)

            # 散点
            for st11, sc_c11 in PATH_COLORS.items():
                mask_st11 = subtype_sorted == st11
                ax11.scatter(pt_sorted[mask_st11], expr_g[mask_st11],
                             c=sc_c11, s=7, alpha=0.30,
                             zorder=2, linewidths=0)

            # 表达曲线（主轴）
            ax11.plot(bin_ctrs11, smooth11,
                      color='#333333', lw=2.2, zorder=4,
                      label='Expression')
            ax11.fill_between(bin_ctrs11,
                              smooth11 - sem11, smooth11 + sem11,
                              color='#333333', alpha=0.12, zorder=3)

            # fate probability（副轴）
            ax11b = ax11.twinx()
            ax11b.plot(bin_ctrs11, fp_la_sm,
                       color=SUBTYPE_COLORS_SC['LA_TAM'],
                       lw=1.8, ls='--', alpha=0.85,
                       label='Fate: LA_TAM')
            ax11b.plot(bin_ctrs11, fp_ifn_sm,
                       color=SUBTYPE_COLORS_SC['IFN_TAM'],
                       lw=1.8, ls='--', alpha=0.85,
                       label='Fate: IFN_TAM')
            ax11b.set_ylabel('Fate Probability', fontsize=8,
                             color='#666666')
            ax11b.tick_params(axis='y', labelsize=7, colors='#666666')
            ax11b.set_ylim(0, 1.05)
            ax11b.spines[['top']].set_visible(False)

            # 相关性标注
            r_la,  _ = spearmanr(pt_sorted, fate_la_arr)
            r_ifn, _ = spearmanr(pt_sorted, fate_ifn_arr)
            r_expr,_ = spearmanr(pt_sorted, expr_g)
            ax11.text(0.97, 0.95,
                      f'expr r={r_expr:.2f}'
                      f'LA_TAM fate r={r_la:.2f}'
                      f'IFN_TAM fate r={r_ifn:.2f}',
                      transform=ax11.transAxes, ha='right', va='top',
                      fontsize=8, color='#333333',
                      bbox=dict(boxstyle='round,pad=0.3', fc='white',
                                ec='#BBBBBB', alpha=0.88))

            ax11.set_title(gene, fontsize=13, fontweight='bold',
                           color='#333333')
            ax11.set_xlabel('Pseudotime', fontsize=9)
            ax11.set_ylabel('log-normalized expression', fontsize=9)
            ax11.spines[['top', 'right']].set_visible(False)
            ax11.set_xlim(pt_sorted.min() - 0.005,
                          pt_sorted.max() + 0.005)

            # 亚型标注
            ymax11 = ax11.get_ylim()[1]
            for st11, x11 in [
                ('Mono',    (pt_sorted.min() + mono_end11)    / 2),
                ('IFN_TAM', (mono_end11 + ifntam_end11)       / 2),
                ('LA_TAM',  (ifntam_end11 + pt_sorted.max())  / 2),
            ]:
                ax11.text(x11, ymax11 * 0.97, st11,
                          ha='center', va='top', fontsize=8,
                          color=PATH_COLORS[st11], fontweight='bold')

            # 合并图例
            handles_expr = [
                Line2D([0],[0], color='#333333', lw=2,
                       label='Expression (smoothed)'),
            ]
            handles_fate = [
                Line2D([0],[0], color=SUBTYPE_COLORS_SC['LA_TAM'],
                       lw=1.8, ls='--',
                       label='(solid = TF AUC, dashed = receptor expr)'),
                Line2D([0],[0], color=SUBTYPE_COLORS_SC['IFN_TAM'],
                       lw=1.8, ls='--', label='Fate: IFN_TAM'),
            ]
            ax11.legend(handles=handles_expr + handles_fate,
                        fontsize=7.5, loc='upper left',
                        framealpha=0.85, edgecolor='none')

        fig11.suptitle(
            'CellRank Fate Trend — Receptor Lineage Driver Genes'
            'Mono → IFN_TAM / LA_TAM  '
            '(solid = expression, dashed = fate probability)',
            fontsize=12, fontweight='bold', y=1.01)

        # 全局图例
        global_handles = [
            Patch(color=c, alpha=0.5, label=s)
            for s, c in PATH_COLORS.items()
        ] + [
            Line2D([0],[0], color='#333333', lw=2,
                   label='Expression (smoothed)'),
            Line2D([0],[0], color=SUBTYPE_COLORS_SC['LA_TAM'],
                   lw=1.8, ls='--', label='Fate: LA_TAM'),
            Line2D([0],[0], color=SUBTYPE_COLORS_SC['IFN_TAM'],
                   lw=1.8, ls='--', label='Fate: IFN_TAM'),
        ]
        fig11.legend(handles=global_handles,
                     loc='lower center', ncol=6,
                     fontsize=8.5, frameon=True,
                     bbox_to_anchor=(0.5, -0.03),
                     framealpha=0.92, edgecolor='#CCCCCC')

        plt.tight_layout()
        save_fig(fig11, 'cellrank_fate_trend_plots')

else:
    if not CELLRANK_OK:
        print("  SKIP Fig 10/11: cellrank not available")
    else:
        print(f"  SKIP Fig 10/11: {volcano_path} not found")

# =============================================================================
#  完成
# =============================================================================
print()
print("=" * 60)
print("ALL DONE — output files:")
for fname in [
    'monomac_paga_trajectory',
    'monomac_trajectory_directed',
    'monomac_gene_dynamics',
    'scenic_tf_dynamics',
    'scenic_tf_heatmap',
    'receptor_pseudotime_heatmap_path',
    'receptor_pseudotime_lineplots_path',
]:
    for fmt in ['png', 'svg']:
        p = os.path.join(TRAJ_DIR, f'{fname}.{fmt}')
        status = '✓' if os.path.exists(p) else '✗ MISSING'
        print(f"  [{status}] {p}")

for fname in [
    'volcano_IFN_TAM_vs_LA_TAM',
    'enrichment_barplot_IFN_TAM_vs_LA_TAM',
]:
    for fmt in ['png', 'svg']:
        p = os.path.join(DEG_DIR, f'{fname}.{fmt}')
        status = '✓' if os.path.exists(p) else '✗ MISSING'
        print(f"  [{status}] {p}")

for fname in [
    'cellrank_gene_relay_heatmap',
    'cellrank_fate_trend_plots',
]:
    for fmt in ['png', 'svg']:
        p = os.path.join(TRAJ_DIR, f'{fname}.{fmt}')
        status = '✓' if os.path.exists(p) else '✗ MISSING'
        print(f"  [{status}] {p}")

print("=" * 60)

# =============================================================================
#  图 10/11：CellRank 基因接力热图 + fate trend 折线图
#  本地 Windows 版本 — 路径、预处理、iroot 全部对齐 paga2.txt
# =============================================================================

import os, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore, spearmanr
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import scanpy as sc
import cellrank as cr
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── 路径配置（本地 Windows）──────────────────────────────────────────────────
BASE_DIR  = translate(r'D:/bulk-download')
DATA_DIR  = os.path.join(BASE_DIR, 'GSE182434')
TRAJ_DIR  = os.path.join(DATA_DIR, 'trajectory')
DEG_DIR   = os.path.join(DATA_DIR, 'IFN_TAM_vs_LA_TAM_DEG')
H5AD_PATH = translate(os.path.join(DATA_DIR, 'adata_processed.h5ad'))
CELLCOMM_DIR = os.path.join(DATA_DIR, 'cellcomm')
DPI = 300
os.makedirs(TRAJ_DIR, exist_ok=True)

# ── 颜色常量（与 paga2.txt 完全一致）────────────────────────────────────────
SUBTYPE_PALETTE = {
    'Mono': '#3498DB', 'DC_1': '#E74C3C', 'LA_TAM': '#2ECC71',
    'IFN_TAM': '#E67E22', 'DC_2': '#9B59B6',
}
SUBTYPES_ORDER    = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
PATH_COLORS       = {'Mono': '#E9C46A', 'IFN_TAM': '#E63946', 'LA_TAM': '#F4A261'}
SUBTYPE_COLORS_SC = {'Mono': '#4878CF', 'IFN_TAM': '#E8601C', 'LA_TAM': '#7BAE7F'}

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
#  STEP 1：数据加载与预处理（完全对齐 paga2.txt 原始逻辑）
# =============================================================================
print("=" * 60)
print("[CellRank] Loading h5ad and computing pseudotime...")
print("=" * 60)

adata = sc.read_h5ad(H5AD_PATH)
adata_mac = adata[adata.obs['cell_type'] == 'Monocytes/Macrophages'].copy()
print(f"Mono/Mac cells: {adata_mac.shape[0]}")

# ── 亚型评分（与 paga2.txt 完全一致）────────────────────────────────────────
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

# ── 关键：显式赋值 .X，清除 log1p 元数据（与 paga2.txt 完全一致）────────────
_layer = adata_mac.layers['log1p_norm']
adata_mac.X = (_layer.tocsr().copy() if sps.issparse(_layer)
               else sps.csr_matrix(_layer.copy()))
print(f"  .X assigned: shape={adata_mac.X.shape}, max={adata_mac.X.max():.2f}")

if 'log1p' in adata_mac.uns:
    del adata_mac.uns['log1p']

# ── HVG + scale + PCA（与 paga2.txt 完全一致）────────────────────────────────
sc.pp.highly_variable_genes(adata_mac, n_top_genes=2000, flavor='seurat')
sc.pp.scale(adata_mac, max_value=10)
sc.tl.pca(adata_mac, n_comps=30, use_highly_variable=True)
sc.pp.neighbors(adata_mac, n_neighbors=15, n_pcs=20)
sc.tl.diffmap(adata_mac, n_comps=15)

# ── categories 顺序强制对齐 SUBTYPES_ORDER（与 paga2.txt 完全一致）──────────
adata_mac.obs['mac_subtype'] = pd.Categorical(
    adata_mac.obs['mac_subtype'],
    categories=SUBTYPES_ORDER, ordered=False)

# ── DPT：iroot 用 np.where 整数位置（与 paga2.txt 完全一致）─────────────────
mono_mask     = adata_mac.obs['mac_subtype'] == 'Mono'
mono_idx      = np.where(mono_mask)[0]
dc1_scores    = adata_mac.obs['score_Mono'].values
root_cell_idx = mono_idx[np.argmax(dc1_scores[mono_idx])]
adata_mac.uns['iroot'] = root_cell_idx
sc.tl.dpt(adata_mac, n_dcs=10)

print(f"DPT pseudotime range: "
      f"{adata_mac.obs['dpt_pseudotime'].min():.3f} – "
      f"{adata_mac.obs['dpt_pseudotime'].max():.3f}")
print("\nMean pseudotime per subtype:")
print(adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()
      .sort_values().round(3))

# ── 路径子集（Mono + IFN_TAM + LA_TAM）──────────────────────────────────────
path_mask  = adata_mac.obs['mac_subtype'].isin(['Mono', 'IFN_TAM', 'LA_TAM'])
adata_path = adata_mac[path_mask].copy()
adata_path = adata_path[adata_path.obs['dpt_pseudotime'].argsort()].copy()
adata_path.obs['mac_subtype'] = adata_path.obs['mac_subtype'].astype('category')

pt_sorted      = adata_path.obs['dpt_pseudotime'].values
subtype_sorted = adata_path.obs['mac_subtype'].values
print(f"Path cells (Mono+IFN_TAM+LA_TAM): {adata_path.n_obs}")

# =============================================================================
#  STEP 2：CellRank 核心计算
# =============================================================================
print("\n[CellRank] Building kernel and computing fate probabilities...")

sc.pp.neighbors(adata_path, n_neighbors=15, n_pcs=20)

pk = cr.kernels.PseudotimeKernel(adata_path, time_key='dpt_pseudotime')
pk.compute_transition_matrix(threshold_scheme='soft', nu=0.5)
print(f"Transition matrix computed: {pk.transition_matrix.shape}")

g2 = cr.estimators.GPCCA(pk)

# ── 自动尝试 n_states=3/4/5，直到三个目标亚型都能匹配──────────────────────
def find_state(target, names):
    t = target.lower()
    for n in names:
        if t in n.lower():
            return n
    return None

matched = {}
for n_st in [3, 4, 5]:
    g2.fit(n_states=n_st, cluster_key='mac_subtype')
    actual_names = list(g2.macrostates_memberships.names)
    print(f"  n_states={n_st} → macrostates: {actual_names}")
    matched = {st: find_state(st, actual_names)
               for st in ['Mono', 'IFN_TAM', 'LA_TAM']}
    print(f"  Matched: {matched}")
    if all(v is not None for v in matched.values()):
        break
    print(f"  Not all targets found, trying n_states={n_st + 1}...")

# ── 设置 terminal states 并计算 fate probabilities ──────────────────────────
terminal_actual = [v for v in matched.values() if v is not None]
print(f"  Setting terminal states: {terminal_actual}")
g2.set_terminal_states(terminal_actual, cluster_key='mac_subtype')
g2.compute_fate_probabilities()

# ── 构建 fate probability DataFrame，列名映射回原始亚型名 ────────────────────
fp_arr  = np.array(g2.fate_probabilities.X)
fp_cols = list(g2.fate_probabilities.names)
fp_df   = pd.DataFrame(fp_arr,
                        index=adata_path.obs_names,
                        columns=fp_cols)

reverse_map = {v: k for k, v in matched.items() if v is not None}
fp_df = fp_df.rename(columns=reverse_map)

for st in ['Mono', 'IFN_TAM', 'LA_TAM']:
    if st not in fp_df.columns:
        fp_df[st] = 0.0
        print(f"  WARNING: {st} not found in fate probabilities, filled with 0")

adata_path.obs['fate_LA_TAM']  = fp_df['LA_TAM'].values
adata_path.obs['fate_IFN_TAM'] = fp_df['IFN_TAM'].values
adata_path.obs['fate_Mono']    = fp_df['Mono'].values
adata_path.obsm['lineages_fwd'] = g2.fate_probabilities

print("  Fate probabilities assigned:")
print(adata_path.obs[['fate_Mono', 'fate_IFN_TAM', 'fate_LA_TAM']]
      .describe().round(3))

# ── Lineage drivers ──────────────────────────────────────────────────────────
la_actual  = matched.get('LA_TAM')
ifn_actual = matched.get('IFN_TAM')

drivers_la = drivers_ifn = None

if la_actual:
    try:
        drivers_la = g2.compute_lineage_drivers(
            lineages=la_actual, use_raw=False,
            layer='log1p_norm', return_drivers=True)
        drivers_la = drivers_la.rename(columns={
            f'{la_actual}_corr': 'LA_TAM_corr',
            f'{la_actual}_pval': 'LA_TAM_pval'})
        print(f"  LA_TAM drivers computed: {drivers_la.shape}")
    except Exception as e:
        print(f"  WARNING: LA_TAM drivers failed: {e}")

if ifn_actual:
    try:
        drivers_ifn = g2.compute_lineage_drivers(
            lineages=ifn_actual, use_raw=False,
            layer='log1p_norm', return_drivers=True)
        drivers_ifn = drivers_ifn.rename(columns={
            f'{ifn_actual}_corr': 'IFN_TAM_corr',
            f'{ifn_actual}_pval': 'IFN_TAM_pval'})
        print(f"  IFN_TAM drivers computed: {drivers_ifn.shape}")
    except Exception as e:
        print(f"  WARNING: IFN_TAM drivers failed: {e}")

# ── 读取 receptor 基因列表 ────────────────────────────────────────────────────
volcano_path = os.path.join(CELLCOMM_DIR,
                             'volcano_data_malignant_vs_normal_monomac.csv')
lr_df   = pd.read_csv(volcano_path)
mal_lr  = lr_df[lr_df['category'] == 'Malignant-enriched'].copy()

def extract_genes(complex_str):
    return [p.strip() for p in str(complex_str).split('_')
            if p.strip() and p.strip() != 'nan']

receptor_genes_all = []
for rec in mal_lr['receptor_complex']:
    receptor_genes_all.extend(extract_genes(rec))
receptor_genes = list(dict.fromkeys(receptor_genes_all))
available_rec  = [g for g in receptor_genes if g in adata_path.var_names]
print(f"  Receptor genes available: {len(available_rec)}")

# ── 构建 receptor_drivers 表 ─────────────────────────────────────────────────
_X_rec = adata_path[:, available_rec].layers['log1p_norm']
if sps.issparse(_X_rec): _X_rec = _X_rec.toarray()
expr_rec = pd.DataFrame(_X_rec, index=adata_path.obs_names,
                         columns=available_rec)

if drivers_la is not None:
    receptor_drivers = drivers_la.loc[
        drivers_la.index.isin(available_rec)].copy()
    receptor_drivers.index.name = 'gene'
    receptor_drivers = receptor_drivers.reset_index()
    receptor_drivers = receptor_drivers.sort_values(
        'LA_TAM_corr', ascending=False)

    # 补充 qval 列（若不存在则用 pval 替代）
    if 'LA_TAM_qval' not in receptor_drivers.columns:
        if 'LA_TAM_pval' in receptor_drivers.columns:
            from statsmodels.stats.multitest import multipletests
            _, qvals, _, _ = multipletests(
                receptor_drivers['LA_TAM_pval'].fillna(1), method='fdr_bh')
            receptor_drivers['LA_TAM_qval'] = qvals
        else:
            receptor_drivers['LA_TAM_qval'] = 1.0

    for st in ['Mono', 'IFN_TAM', 'LA_TAM']:
        mask = adata_path.obs['mac_subtype'] == st
        receptor_drivers[f'mean_expr_{st}'] = receptor_drivers['gene'].map(
            expr_rec.loc[mask].mean())

    receptor_drivers['direction'] = receptor_drivers['LA_TAM_corr'].apply(
        lambda r: 'positive' if r > 0 else 'negative')
    receptor_drivers['significant'] = receptor_drivers['LA_TAM_qval'] <= 0.05

    receptor_drivers.to_csv(
        os.path.join(TRAJ_DIR, 'cellrank_LA_TAM_fate_drivers_receptors.csv'),
        index=False)
    print(f"  Saved receptor drivers: {len(receptor_drivers)} genes")
else:
    # fallback：用伪时间 Spearman 相关性替代
    print("  WARNING: drivers_la is None, falling back to pseudotime correlation")
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests
    rows = []
    for gene in available_rec:
        r, p = spearmanr(pt_sorted, expr_rec[gene].values)
        rows.append({'gene': gene, 'LA_TAM_corr': r, 'LA_TAM_pval': p})
    receptor_drivers = pd.DataFrame(rows)
    _, qvals, _, _ = multipletests(
        receptor_drivers['LA_TAM_pval'].fillna(1), method='fdr_bh')
    receptor_drivers['LA_TAM_qval'] = qvals
    receptor_drivers['significant'] = receptor_drivers['LA_TAM_qval'] <= 0.05
    receptor_drivers = receptor_drivers.sort_values(
        'LA_TAM_corr', ascending=False)

# ── top_la / top_ifn（用于图 10 基因接力）────────────────────────────────────
top_la = (drivers_la.loc[
    [g for g in available_rec if g in drivers_la.index], 'LA_TAM_corr']
    .sort_values(ascending=False).head(10)
    if drivers_la is not None else pd.Series(dtype=float))

top_ifn = (drivers_ifn.loc[
    [g for g in available_rec if g in drivers_ifn.index], 'IFN_TAM_corr']
    .sort_values(ascending=False).head(10)
    if drivers_ifn is not None else pd.Series(dtype=float))

# fallback：不足时用伪时间相关性补充
if len(top_la) + len(top_ifn) < 4:
    print("  Falling back to pseudotime-correlated genes for relay plot")
    from scipy.stats import spearmanr
    fallback_rows = []
    for gene in available_rec:
        r, _ = spearmanr(pt_sorted, expr_rec[gene].values)
        fallback_rows.append((gene, r))
    fallback_sorted = sorted(fallback_rows, key=lambda x: -x[1])
    top_la  = pd.Series({g: r for g, r in fallback_sorted[:5]}, dtype=float)
    top_ifn = pd.Series({g: r for g, r in fallback_sorted[5:10]}, dtype=float)

relay_genes = list(dict.fromkeys(
    top_la.index.tolist() + top_ifn.index.tolist()))[:20]
print(f"  Relay genes ({len(relay_genes)}): {relay_genes}")

# =============================================================================
#  图 10：CellRank 基因接力热图
# =============================================================================
if need_redraw('cellrank_gene_relay_heatmap'):
    print("\n[Fig 10] CellRank gene relay heatmap...")

    n_bins10    = 100; sigma10 = 3
    bin_edges10 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins10 + 1)

    relay_avail = [g for g in relay_genes if g in adata_path.var_names]
    _X10 = adata_path[:, relay_avail].layers['log1p_norm']
    if sps.issparse(_X10): _X10 = _X10.toarray()
    expr10 = pd.DataFrame(_X10, index=adata_path.obs_names,
                           columns=relay_avail)

    expr_binned10  = np.zeros((len(relay_avail), n_bins10))
    subtype_bins10 = []
    for b in range(n_bins10):
        mask_b = ((pt_sorted >= bin_edges10[b]) &
                  (pt_sorted <  bin_edges10[b + 1]))
        if mask_b.sum() > 0:
            expr_binned10[:, b] = expr10[relay_avail].values[
                mask_b].mean(axis=0)
            subtype_bins10.append(
                Counter(subtype_sorted[mask_b]).most_common(1)[0][0])
        else:
            expr_binned10[:, b] = (expr_binned10[:, b-1] if b > 0 else 0)
            subtype_bins10.append(
                subtype_bins10[-1] if subtype_bins10 else 'Mono')

    expr_smooth10 = np.array([
        gaussian_filter1d(expr_binned10[i], sigma10)
        for i in range(len(relay_avail))])
    expr_z10 = np.apply_along_axis(
        lambda x: zscore(x) if x.std() > 0 else x, 1, expr_smooth10)

    peak_idx10  = np.argmax(expr_z10, axis=1)
    order10     = np.argsort(peak_idx10)
    genes_ord10 = [relay_avail[i] for i in order10]
    expr_z_ord10 = expr_z10[order10]

    lineage_colors10 = {
        'LA_TAM': '#2ECC71', 'IFN_TAM': '#E67E22',
        'Both':   '#9B59B6', 'Other':   '#AAAAAA',
    }
    gene_lineage10 = {}
    for g in genes_ord10:
        in_la  = g in top_la.index
        in_ifn = g in top_ifn.index
        if in_la and in_ifn:   gene_lineage10[g] = 'Both'
        elif in_la:            gene_lineage10[g] = 'LA_TAM'
        elif in_ifn:           gene_lineage10[g] = 'IFN_TAM'
        else:                  gene_lineage10[g] = 'Other'

    sub_num10  = np.array([[{'Mono': 0, 'IFN_TAM': 1, 'LA_TAM': 2}
                             .get(s, 3) for s in subtype_bins10]])
    cmap_sub10 = ListedColormap(['#E9C46A', '#E63946', '#F4A261', '#CCCCCC'])

    fig10 = plt.figure(figsize=(14, max(8, len(genes_ord10) * 0.55 + 2)))
    fig10.patch.set_facecolor('white')
    gs10 = gridspec.GridSpec(2, 2,
                              height_ratios=[0.04, 1],
                              width_ratios=[1, 0.03],
                              hspace=0.02, wspace=0.02)
    ax_bar10  = fig10.add_subplot(gs10[0, 0])
    ax_heat10 = fig10.add_subplot(gs10[1, 0])
    ax_cbar10 = fig10.add_subplot(gs10[1, 1])

    ax_bar10.imshow(sub_num10, aspect='auto', cmap=cmap_sub10,
                    vmin=0, vmax=3, interpolation='nearest')
    ax_bar10.set_xticks([]); ax_bar10.set_yticks([0])
    ax_bar10.set_yticklabels(['Subtype'], fontsize=8)
    ax_bar10.tick_params(left=False)
    ax_bar10.legend(
        handles=[Patch(color='#E9C46A', label='Mono'),
                 Patch(color='#E63946', label='IFN_TAM'),
                 Patch(color='#F4A261', label='LA_TAM')],
        loc='upper right', bbox_to_anchor=(1.12, 2.8),
        fontsize=8, frameon=True, ncol=3)

    cmap_heat10 = LinearSegmentedColormap.from_list(
        'relay', ['#2166AC', '#F7F7F7', '#D6604D'])
    im10 = ax_heat10.imshow(expr_z_ord10, aspect='auto',
                             cmap=cmap_heat10, vmin=-2, vmax=2,
                             interpolation='nearest')

    ax_heat10.set_yticks(range(len(genes_ord10)))
    ax_heat10.set_yticklabels(genes_ord10, fontsize=8)
    for tick, gene in zip(ax_heat10.get_yticklabels(), genes_ord10):
        lin = gene_lineage10.get(gene, 'Other')
        tick.set_color(lineage_colors10[lin])
        if lin != 'Other': tick.set_fontweight('bold')

    for idx, gene in enumerate(genes_ord10):
        lin = gene_lineage10.get(gene, 'Other')
        ax_heat10.annotate(
            '', xy=(n_bins10 + 1.5, idx),
            xytext=(n_bins10 + 0.5, idx),
            xycoords='data',
            arrowprops=dict(arrowstyle='-',
                            color=lineage_colors10[lin], lw=3))

    xtick_pos10  = np.linspace(0, n_bins10 - 1, 6)
    xtick_vals10 = np.linspace(pt_sorted.min(), pt_sorted.max(), 6)
    ax_heat10.set_xticks(xtick_pos10)
    ax_heat10.set_xticklabels(
        [f'{v:.2f}' for v in xtick_vals10], fontsize=8.5)
    ax_heat10.set_xlabel(
        'Pseudotime  (Mono → IFN_TAM / LA_TAM)', fontsize=10)
    ax_heat10.set_ylabel(
        'Receptor genes (CellRank lineage drivers)', fontsize=9)

    prev10 = subtype_bins10[0]
    for b, st in enumerate(subtype_bins10):
        if st != prev10:
            ax_heat10.axvline(b - 0.5, color='white',
                               lw=1.5, ls='--', alpha=0.8)
            prev10 = st

    plt.colorbar(im10, cax=ax_cbar10, label='Z-score')
    ax_cbar10.tick_params(labelsize=7)

    fig10.legend(
        handles=[Patch(color=c, label=l)
                 for l, c in lineage_colors10.items() if l != 'Other'],
        title='Lineage driver', title_fontsize=8.5,
        fontsize=8, loc='lower right',
        bbox_to_anchor=(0.88, 0.02), framealpha=0.9)

    fig10.suptitle(
        f'CellRank Lineage Driver Receptor Genes — Gene Relay Heatmap\n'
        f'Mono → IFN_TAM / LA_TAM  ({len(genes_ord10)} genes)',
        fontsize=12, fontweight='bold', y=1.005)
    plt.tight_layout()
    save_fig(fig10, 'cellrank_gene_relay_heatmap')

# =============================================================================
#  图 11：Fate trend 折线图（top 4 receptor drivers）
# =============================================================================
if need_redraw('cellrank_fate_trend_plots'):
    print("\n[Fig 11] CellRank fate trend plots...")

    trend_genes = list(dict.fromkeys(
        top_la.index[:2].tolist() +
        top_ifn.index[:2].tolist()))[:4]
    if len(trend_genes) < 4:
        extra = [g for g in relay_genes if g not in trend_genes]
        trend_genes += extra[:4 - len(trend_genes)]
    print(f"  Fate trend genes: {trend_genes}")

    n_bins11    = 80; sigma11 = 3
    bin_edges11 = np.linspace(pt_sorted.min(), pt_sorted.max(), n_bins11 + 1)
    bin_ctrs11  = (bin_edges11[:-1] + bin_edges11[1:]) / 2

    trend_avail = [g for g in trend_genes if g in adata_path.var_names]
    _X11 = adata_path[:, trend_avail].layers['log1p_norm']
    if sps.issparse(_X11): _X11 = _X11.toarray()
    expr11 = pd.DataFrame(_X11, index=adata_path.obs_names,
                           columns=trend_avail)

    fate_la_arr  = adata_path.obs['fate_LA_TAM'].values
    fate_ifn_arr = adata_path.obs['fate_IFN_TAM'].values

    mono_end11   = adata_path.obs[
        adata_path.obs['mac_subtype'] == 'Mono']['dpt_pseudotime'].max()
    ifntam_end11 = adata_path.obs[
        adata_path.obs['mac_subtype'] == 'IFN_TAM']['dpt_pseudotime'].max()

    fig11, axes11 = plt.subplots(2, 2, figsize=(14, 10))
    fig11.patch.set_facecolor('white')

    for ax11, gene in zip(axes11.flatten(), trend_avail):
        expr_g = expr11[gene].values

        bm11 = np.zeros(n_bins11); bs11 = np.zeros(n_bins11)
        for b in range(n_bins11):
            mask_b = ((pt_sorted >= bin_edges11[b]) &
                      (pt_sorted <  bin_edges11[b + 1]))
            if mask_b.sum() > 0:
                v = expr_g[mask_b]
                bm11[b] = v.mean()
                bs11[b] = v.std() / np.sqrt(len(v)) if len(v) > 1 else 0
        smooth11 = gaussian_filter1d(bm11, sigma=sigma11)
        sem11    = gaussian_filter1d(bs11, sigma=sigma11)

        fp_la_bin  = np.zeros(n_bins11)
        fp_ifn_bin = np.zeros(n_bins11)
        for b in range(n_bins11):
            mask_b = ((pt_sorted >= bin_edges11[b]) &
                      (pt_sorted <  bin_edges11[b + 1]))
            if mask_b.sum() > 0:
                fp_la_bin[b]  = fate_la_arr[mask_b].mean()
                fp_ifn_bin[b] = fate_ifn_arr[mask_b].mean()
        fp_la_sm  = gaussian_filter1d(fp_la_bin,  sigma=sigma11)
        fp_ifn_sm = gaussian_filter1d(fp_ifn_bin, sigma=sigma11)

        ax11.axvspan(pt_sorted.min(), mono_end11,
                     alpha=0.10, color=PATH_COLORS['Mono'],    zorder=0)
        ax11.axvspan(mono_end11, ifntam_end11,
                     alpha=0.10, color=PATH_COLORS['IFN_TAM'], zorder=0)
        ax11.axvspan(ifntam_end11, pt_sorted.max(),
                     alpha=0.10, color=PATH_COLORS['LA_TAM'],  zorder=0)
        ax11.axvline(mono_end11,   color='#AAAAAA', lw=0.8, ls='--', zorder=1)
        ax11.axvline(ifntam_end11, color='#AAAAAA', lw=0.8, ls='--', zorder=1)

        for st11, sc_c11 in PATH_COLORS.items():
            mask_st11 = subtype_sorted == st11
            ax11.scatter(pt_sorted[mask_st11], expr_g[mask_st11],
                         c=sc_c11, s=7, alpha=0.30,
                         zorder=2, linewidths=0)

        ax11.plot(bin_ctrs11, smooth11,
                  color='#333333', lw=2.2, zorder=4, label='Expression')
        ax11.fill_between(bin_ctrs11,
                          smooth11 - sem11, smooth11 + sem11,
                          color='#333333', alpha=0.12, zorder=3)

        ax11b = ax11.twinx()
        ax11b.plot(bin_ctrs11, fp_la_sm,
                   color=SUBTYPE_COLORS_SC['LA_TAM'],
                   lw=1.8, ls='--', alpha=0.85, label='Fate: LA_TAM')
        ax11b.plot(bin_ctrs11, fp_ifn_sm,
                   color=SUBTYPE_COLORS_SC['IFN_TAM'],
                   lw=1.8, ls='--', alpha=0.85, label='Fate: IFN_TAM')
        ax11b.set_ylabel('Fate Probability', fontsize=8, color='#666666')
        ax11b.tick_params(axis='y', labelsize=7, colors='#666666')
        ax11b.set_ylim(0, 1.05)
        ax11b.spines[['top']].set_visible(False)

        r_expr, _ = spearmanr(pt_sorted, expr_g)
        r_la,   _ = spearmanr(pt_sorted, fate_la_arr)
        r_ifn,  _ = spearmanr(pt_sorted, fate_ifn_arr)
        ax11.text(0.97, 0.95,
                  f'expr r={r_expr:.2f}\nLA_TAM fate r={r_la:.2f}\nIFN_TAM fate r={r_ifn:.2f}',
                  transform=ax11.transAxes, ha='right', va='top',
                  fontsize=8, color='#333333',
                  bbox=dict(boxstyle='round,pad=0.3', fc='white',
                            ec='#BBBBBB', alpha=0.88))

        ax11.set_title(gene, fontsize=13, fontweight='bold', color='#333333')
        ax11.set_xlabel('Pseudotime', fontsize=9)
        ax11.set_ylabel('log-normalized expression', fontsize=9)
        ax11.spines[['top', 'right']].set_visible(False)
        ax11.set_xlim(pt_sorted.min() - 0.005, pt_sorted.max() + 0.005)

        ymax11 = ax11.get_ylim()[1]
        for st11, x11 in [
            ('Mono',    (pt_sorted.min() + mono_end11)    / 2),
            ('IFN_TAM', (mono_end11 + ifntam_end11)       / 2),
            ('LA_TAM',  (ifntam_end11 + pt_sorted.max())  / 2),
        ]:
            ax11.text(x11, ymax11 * 0.97, st11,
                      ha='center', va='top', fontsize=8,
                      color=PATH_COLORS[st11], fontweight='bold')

    global_handles = [
        Patch(color=c, alpha=0.5, label=s)
        for s, c in PATH_COLORS.items()
    ] + [
        Line2D([0], [0], color='#333333', lw=2,
               label='Expression (smoothed)'),
        Line2D([0], [0], color=SUBTYPE_COLORS_SC['LA_TAM'],
               lw=1.8, ls='--', label='Fate: LA_TAM'),
        Line2D([0], [0], color=SUBTYPE_COLORS_SC['IFN_TAM'],
               lw=1.8, ls='--', label='Fate: IFN_TAM'),
    ]
    fig11.legend(handles=global_handles,
                 loc='lower center', ncol=6,
                 fontsize=8.5, frameon=True,
                 bbox_to_anchor=(0.5, -0.03),
                 framealpha=0.92, edgecolor='#CCCCCC')

    fig11.suptitle(
        'CellRank Fate Trend — Receptor Lineage Driver Genes\n'
        'Mono → IFN_TAM / LA_TAM\n'
        '(solid = expression, dashed = fate probability)',
        fontsize=12, fontweight='bold', y=1.01)

    plt.tight_layout()
    save_fig(fig11, 'cellrank_fate_trend_plots')

# =============================================================================
#  完成
# =============================================================================
print()
print("=" * 60)
print("CellRank DONE — output files:")
for fname in ['cellrank_gene_relay_heatmap', 'cellrank_fate_trend_plots']:
    for fmt in ['png', 'svg']:
        p = os.path.join(TRAJ_DIR, f'{fname}.{fmt}')
        status = '✓' if os.path.exists(p) else '✗ MISSING'
        print(f"  [{status}] {p}")
print("=" * 60)


