# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

# ── Load adata ─────────────────────────────────────────────────────────────────

adata = sc.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_processed.h5ad'))

adata_mac = adata[adata.obs['cell_type'] == 'Monocytes/Macrophages'].copy()
print(f"Mono/Mac cells: {adata_mac.shape[0]}")

# ── Marker genes (from notebook Cell 733) ─────────────────────────────────────

subtype_markers = {
    'Mono':    ['FCN1', 'S100A9', 'S100A8', 'S100A4', 'APOBEC3A'],
    'DC_1':    ['LTB', 'CLEC10A', 'CD1C', 'JAML', 'CD1E'],
    'LA_TAM':  ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD'],
    'IFN_TAM': ['MT1H', 'MT1G', 'CCL8', 'CCL2', 'MT1X'],
    'DC_2':    ['DNASE1L3', 'RGCC', 'CST3', 'CLEC9A', 'SNX3'],
}

# Check availability
print("\nMarker gene availability:")
for ct, genes in subtype_markers.items():
    present = [g for g in genes if g in adata_mac.var_names]
    missing = [g for g in genes if g not in adata_mac.var_names]
    print(f"  {ct:10s}: {len(present)}/{len(genes)} present | missing: {missing if missing else 'none'}")

# ── Score each subtype ─────────────────────────────────────────────────────────

for ct, genes in subtype_markers.items():
    present_genes = [g for g in genes if g in adata_mac.var_names]
    sc.tl.score_genes(adata_mac, gene_list=present_genes,
                      score_name=f'score_{ct}', use_raw=False)

# ── Assign subtype by highest score ───────────────────────────────────────────

subtypes   = list(subtype_markers.keys())
score_cols = [f'score_{ct}' for ct in subtypes]
score_mat  = adata_mac.obs[score_cols].values
best_idx   = np.argmax(score_mat, axis=1)

adata_mac.obs['mac_subtype'] = pd.Categorical(
    [subtypes[i] for i in best_idx], categories=subtypes
)

print("\nSubtype distribution:")
print(adata_mac.obs['mac_subtype'].value_counts())

print("\nMean scores per assigned subtype:")
print(adata_mac.obs[score_cols + ['mac_subtype']].groupby('mac_subtype')[score_cols].mean().round(3))
import h5py

def clean_uns(d):
    """递归清理 uns 中无法序列化的 h5py.Empty，替换为 None"""
    cleaned = {}
    for k, v in d.items():
        if isinstance(v, h5py._hl.base.Empty):
            print(f'  清理 Empty: {k}')
            cleaned[k] = None          # 或直接 continue 跳过不保存
        elif isinstance(v, dict):
            cleaned[k] = clean_uns(v)
        else:
            cleaned[k] = v
    return cleaned

adata_mac.uns = clean_uns(adata_mac.uns)
# ── 保存注释好的巨噬细胞亚群 adata ──────────────────────────────────────────

adata_mac.write_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac_subtyped.h5ad'))

print(f"✓ 已保存: adata_mac_subtyped.h5ad")
print(f"  细胞数: {adata_mac.shape[0]}, 基因数: {adata_mac.shape[1]}")
print(f"  亚群分布:\n{adata_mac.obs['mac_subtype'].value_counts()}")


import scanpy as sc
import numpy as np
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── Re-process from log-norm .X ───────────────────────────────────────────────
# .X is already log-normalized (max ~7.5)
# Store scores before processing
score_backup = adata_mac.obs[['mac_subtype'] + [f'score_{ct}' for ct in subtypes]].copy()

# HVG selection
sc.pp.highly_variable_genes(adata_mac, n_top_genes=2000, flavor='seurat')
print(f"HVGs: {adata_mac.var['highly_variable'].sum()}")

# Scale (for PCA)
sc.pp.scale(adata_mac, max_value=10)

# PCA
sc.tl.pca(adata_mac, n_comps=30, use_highly_variable=True)

# Restore scores (scaling overwrites obs? No, obs is preserved)
# Verify
print(f"mac_subtype still in obs: {'mac_subtype' in adata_mac.obs.columns}")

# Neighbors graph
sc.pp.neighbors(adata_mac, n_neighbors=15, n_pcs=20)

# Diffusion map (required for DPT pseudotime)
sc.tl.diffmap(adata_mac, n_comps=15)
print("Diffusion map computed")

# UMAP (standard, for comparison)
sc.tl.umap(adata_mac, min_dist=0.3)
print("UMAP computed")

# Leiden clustering for PAGA
sc.tl.leiden(adata_mac, resolution=0.5, key_added='leiden_mac')
print(f"Leiden clusters: {adata_mac.obs['leiden_mac'].nunique()}")


import scanpy as sc
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── PAGA on mac_subtype ────────────────────────────────────────────────────────
sc.tl.paga(adata_mac, groups='mac_subtype')
print("PAGA connectivity matrix:")
conn = adata_mac.uns['paga']['connectivities'].toarray()
subtypes_order = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
conn_df = pd.DataFrame(conn, index=subtypes_order, columns=subtypes_order)
print(conn_df.round(3))

# ── DPT pseudotime: root = Mono (most progenitor-like) ────────────────────────
# Find the Mono cell closest to the diffusion map centroid as root
mono_mask = adata_mac.obs['mac_subtype'] == 'Mono'
mono_idx = np.where(mono_mask)[0]

# Use the cell with highest DC1 score (most "monocyte-like") as root
dc1_scores = adata_mac.obs['score_Mono'].values
root_cell_idx = mono_idx[np.argmax(dc1_scores[mono_idx])]
print(f"\nRoot cell index: {root_cell_idx} (Mono, highest score_Mono)")

adata_mac.uns['iroot'] = root_cell_idx
sc.tl.dpt(adata_mac, n_dcs=10)
print(f"DPT pseudotime range: {adata_mac.obs['dpt_pseudotime'].min():.3f} – {adata_mac.obs['dpt_pseudotime'].max():.3f}")

# Mean pseudotime per subtype
print("\nMean pseudotime per subtype:")
pt_means = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean().sort_values()
print(pt_means.round(3))


import scanpy as sc
import numpy as np
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ── PAGA-initialized UMAP ─────────────────────────────────────────────────────
# This uses PAGA positions as initialization → layout reflects differentiation topology
sc.pl.paga(adata_mac, show=False)  # compute PAGA positions
sc.tl.umap(adata_mac, init_pos='paga', min_dist=0.5, random_state=35)
print("PAGA-initialized UMAP computed")
print(f"UMAP shape: {adata_mac.obsm['X_umap'].shape}")

# Store both UMAPs
adata_mac.obsm['X_umap_paga'] = adata_mac.obsm['X_umap'].copy()
print("Stored as X_umap_paga")


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# ── Color palette ──────────────────────────────────────────────────────────────
subtype_palette = {
    'Mono':    '#3498DB',
    'DC_1':    '#E74C3C',
    'LA_TAM':  '#2ECC71',
    'IFN_TAM': '#E67E22',
    'DC_2':    '#9B59B6',
}
subtypes_order = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']

umap_coords = adata_mac.obsm['X_umap_paga']
subtypes_arr = adata_mac.obs['mac_subtype'].values
pt_arr       = adata_mac.obs['dpt_pseudotime'].values

# ── Figure layout: 1×3 ────────────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 6.5))
gs  = fig.add_gridspec(1, 3, wspace=0.35)

# ─────────────────────────────────────────────────────────────────────────────
# Panel A: PAGA topology graph (manual drawing for full control)
# ─────────────────────────────────────────────────────────────────────────────
ax_paga = fig.add_subplot(gs[0])

# Get PAGA connectivity
conn = adata_mac.uns['paga']['connectivities'].toarray()
# Get PAGA positions (from sc.pl.paga call above)
paga_pos = adata_mac.uns['paga']['pos']  # shape (n_groups, 2)

# Node sizes proportional to cell count
n_cells = {st: (adata_mac.obs['mac_subtype'] == st).sum() for st in subtypes_order}
node_sizes = np.array([n_cells[st] for st in subtypes_order])
node_sizes_scaled = (node_sizes / node_sizes.max()) * 1200 + 300

# Draw edges (thickness = connectivity weight)
threshold = 0.15  # only draw edges above this
for i, st_i in enumerate(subtypes_order):
    for j, st_j in enumerate(subtypes_order):
        if j <= i:
            continue
        w = conn[i, j]
        if w < threshold:
            continue
        x0, y0 = paga_pos[i]
        x1, y1 = paga_pos[j]
        lw = w * 8
        alpha = min(0.9, 0.3 + w * 0.7)
        ax_paga.plot([x0, x1], [y0, y1], '-',
                     color='#888888', lw=lw, alpha=alpha, zorder=1)
        # Edge weight label
        mid_x, mid_y = (x0 + x1) / 2, (y0 + y1) / 2
        ax_paga.text(mid_x, mid_y, f'{w:.2f}', fontsize=7,
                     ha='center', va='center', color='#555',
                     bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.7, lw=0))

# Draw nodes
for i, st in enumerate(subtypes_order):
    x, y = paga_pos[i]
    color = subtype_palette[st]
    ax_paga.scatter(x, y, s=node_sizes_scaled[i], c=color,
                    zorder=3, edgecolors='white', linewidths=2)
    ax_paga.text(x, y + 0.08, st, ha='center', va='bottom',
                 fontsize=10, fontweight='bold', color=color,
                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8, lw=0))
    ax_paga.text(x, y - 0.08, f'n={n_cells[st]}', ha='center', va='top',
                 fontsize=8, color='#555')

ax_paga.set_title('PAGA Topology Graph\n(edge weight = connectivity)', 
                  fontsize=12, fontweight='bold')
ax_paga.set_xlim(paga_pos[:, 0].min() - 0.3, paga_pos[:, 0].max() + 0.3)
ax_paga.set_ylim(paga_pos[:, 1].min() - 0.3, paga_pos[:, 1].max() + 0.3)
ax_paga.axis('off')

# ─────────────────────────────────────────────────────────────────────────────
# Panel B: PAGA-initialized UMAP colored by subtype
# ─────────────────────────────────────────────────────────────────────────────
ax_sub = fig.add_subplot(gs[1])

for st in subtypes_order:
    mask = subtypes_arr == st
    ax_sub.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                   c=subtype_palette[st], s=18, alpha=0.75,
                   linewidths=0, label=st, rasterized=True)

# Centroid labels
for st in subtypes_order:
    mask = subtypes_arr == st
    cx, cy = umap_coords[mask, 0].mean(), umap_coords[mask, 1].mean()
    ax_sub.text(cx, cy, st, fontsize=8.5, fontweight='bold',
                ha='center', va='center', color='white',
                bbox=dict(boxstyle='round,pad=0.25', fc=subtype_palette[st],
                          alpha=0.85, lw=0))

ax_sub.set_title('PAGA-initialized UMAP\n(Cell Subtypes)', fontsize=12, fontweight='bold')
ax_sub.set_xlabel('UMAP1', fontsize=10)
ax_sub.set_ylabel('UMAP2', fontsize=10)
ax_sub.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_sub.spines.values():
    sp.set_visible(False)

legend_handles = [mpatches.Patch(color=subtype_palette[st], label=st) for st in subtypes_order]
ax_sub.legend(handles=legend_handles, fontsize=8, loc='lower right',
              framealpha=0.85, edgecolor='none', title='Subtype', title_fontsize=8.5)

# ─────────────────────────────────────────────────────────────────────────────
# Panel C: PAGA-initialized UMAP colored by DPT pseudotime
# ─────────────────────────────────────────────────────────────────────────────
ax_pt = fig.add_subplot(gs[2])

sc_pt = ax_pt.scatter(umap_coords[:, 0], umap_coords[:, 1],
                      c=pt_arr, cmap='viridis', s=18, alpha=0.85,
                      linewidths=0, rasterized=True)

# Subtype centroid arrows showing pseudotime order
pt_means = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()
centroids = {st: umap_coords[subtypes_arr == st].mean(axis=0) for st in subtypes_order}

# Draw trajectory arrows: Mono → IFN_TAM → LA_TAM, Mono → DC_2 → DC_1
trajectory_edges = [
    ('Mono', 'IFN_TAM'),
    ('IFN_TAM', 'LA_TAM'),
    ('Mono', 'DC_2'),
    ('DC_2', 'DC_1'),
]
for src, tgt in trajectory_edges:
    x0, y0 = centroids[src]
    x1, y1 = centroids[tgt]
    ax_pt.annotate('', xy=(x1, y1), xytext=(x0, y0),
                   arrowprops=dict(arrowstyle='->', color='white',
                                   lw=2.5, mutation_scale=18))
    ax_pt.annotate('', xy=(x1, y1), xytext=(x0, y0),
                   arrowprops=dict(arrowstyle='->', color='#333333',
                                   lw=1.5, mutation_scale=16))

# Subtype centroid labels
for st in subtypes_order:
    cx, cy = centroids[st]
    ax_pt.text(cx, cy, st, fontsize=8, fontweight='bold',
               ha='center', va='center', color='white',
               bbox=dict(boxstyle='round,pad=0.2', fc=subtype_palette[st],
                         alpha=0.85, lw=0))

cbar = fig.colorbar(sc_pt, ax=ax_pt, shrink=0.7, pad=0.02)
cbar.set_label('Pseudotime\n(DPT)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax_pt.set_title('PAGA-initialized UMAP\n(DPT Pseudotime)', fontsize=12, fontweight='bold')
ax_pt.set_xlabel('UMAP1', fontsize=10)
ax_pt.set_ylabel('UMAP2', fontsize=10)
ax_pt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_pt.spines.values():
    sp.set_visible(False)

# ── Save ──────────────────────────────────────────────────────────────────────
import os
os.makedirs('trajectory', exist_ok=True)

plt.savefig('trajectory/monomac_paga_trajectory.png',
            dpi=150, bbox_inches='tight')
plt.savefig('trajectory/monomac_paga_trajectory.svg',
            bbox_inches='tight')
plt.close()
print("Figure saved")


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import scipy.sparse as sp
import warnings
warnings.filterwarnings('ignore')

subtype_palette = {
    'Mono':    '#3498DB',
    'DC_1':    '#E74C3C',
    'LA_TAM':  '#2ECC71',
    'IFN_TAM': '#E67E22',
    'DC_2':    '#9B59B6',
}
subtypes_order = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']

umap_coords = adata_mac.obsm['X_umap_paga']
subtypes_arr = adata_mac.obs['mac_subtype'].values
pt_arr       = adata_mac.obs['dpt_pseudotime'].values

# ── Compute centroids ──────────────────────────────────────────────────────────
centroids = {st: umap_coords[subtypes_arr == st].mean(axis=0) for st in subtypes_order}
pt_means  = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()

# ── Trajectory edges based on PAGA connectivity + pseudotime direction ─────────
conn = adata_mac.uns['paga']['connectivities'].toarray()
st_idx = {st: i for i, st in enumerate(subtypes_order)}

# Build directed edges: draw arrow from lower PT → higher PT if connectivity > threshold
threshold = 0.25
directed_edges = []
for i, st_i in enumerate(subtypes_order):
    for j, st_j in enumerate(subtypes_order):
        if j <= i:
            continue
        w = conn[i, j]
        if w < threshold:
            continue
        # Direction: lower pseudotime → higher pseudotime
        if pt_means[st_i] < pt_means[st_j]:
            directed_edges.append((st_i, st_j, w))
        else:
            directed_edges.append((st_j, st_i, w))

print("Directed edges (src → tgt, weight):")
for src, tgt, w in sorted(directed_edges, key=lambda x: -x[2]):
    print(f"  {src} → {tgt}  (w={w:.3f}, ΔPT={pt_means[tgt]-pt_means[src]:.3f})")

# ── Figure: 2-panel (PAGA directed + pseudotime UMAP) ─────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

# ─── Panel A: PAGA-UMAP with directed arrows ──────────────────────────────────
ax = axes[0]

# Background scatter
for st in subtypes_order:
    mask = subtypes_arr == st
    ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
               c=subtype_palette[st], s=20, alpha=0.65,
               linewidths=0, rasterized=True)

# Draw directed arrows between centroids
for src, tgt, w in directed_edges:
    x0, y0 = centroids[src]
    x1, y1 = centroids[tgt]
    # Shorten arrow to not overlap node labels
    dx, dy = x1 - x0, y1 - y0
    dist = np.sqrt(dx**2 + dy**2)
    shrink = 0.22  # fraction to shrink from each end
    xs = x0 + dx * shrink
    ys = y0 + dy * shrink
    xe = x1 - dx * shrink
    ye = y1 - dy * shrink
    lw = 1.5 + w * 4
    # White outline
    ax.annotate('', xy=(xe, ye), xytext=(xs, ys),
                arrowprops=dict(arrowstyle='->', color='white',
                                lw=lw + 1.5, mutation_scale=20))
    # Colored arrow
    ax.annotate('', xy=(xe, ye), xytext=(xs, ys),
                arrowprops=dict(arrowstyle='->', color='#333333',
                                lw=lw, mutation_scale=18))
    # Weight label at midpoint
    mx, my = (xs + xe) / 2, (ys + ye) / 2
    ax.text(mx, my, f'{w:.2f}', fontsize=7.5, ha='center', va='center',
            color='#333', bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.8, lw=0))

# Centroid labels
for st in subtypes_order:
    cx, cy = centroids[st]
    n = (subtypes_arr == st).sum()
    ax.text(cx, cy, f'{st}\n(n={n})', fontsize=8.5, fontweight='bold',
            ha='center', va='center', color='white',
            bbox=dict(boxstyle='round,pad=0.3', fc=subtype_palette[st],
                      alpha=0.9, lw=0.5, edgecolor='white'))

ax.set_title('Mono/Mac Differentiation Trajectory\n(PAGA + DPT, arrows = direction)',
             fontsize=12, fontweight='bold')
ax.set_xlabel('UMAP1', fontsize=10)
ax.set_ylabel('UMAP2', fontsize=10)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax.spines.values():
    sp.set_visible(False)

# Pseudotime colorbar legend
sm = cm.ScalarMappable(cmap='viridis', norm=mcolors.Normalize(0, 1))
sm.set_array([])

# ─── Panel B: Pseudotime gradient ─────────────────────────────────────────────
ax2 = axes[1]

sc_pt = ax2.scatter(umap_coords[:, 0], umap_coords[:, 1],
                    c=pt_arr, cmap='viridis', s=20, alpha=0.85,
                    linewidths=0, rasterized=True,
                    vmin=0, vmax=1)

# Directed arrows on pseudotime panel too
for src, tgt, w in directed_edges:
    x0, y0 = centroids[src]
    x1, y1 = centroids[tgt]
    dx, dy = x1 - x0, y1 - y0
    shrink = 0.22
    xs = x0 + dx * shrink
    ys = y0 + dy * shrink
    xe = x1 - dx * shrink
    ye = y1 - dy * shrink
    ax2.annotate('', xy=(xe, ye), xytext=(xs, ys),
                 arrowprops=dict(arrowstyle='->', color='white',
                                 lw=3.5, mutation_scale=20))
    ax2.annotate('', xy=(xe, ye), xytext=(xs, ys),
                 arrowprops=dict(arrowstyle='->', color='#ff4444',
                                 lw=2.0, mutation_scale=18))

# Centroid labels with PT value
for st in subtypes_order:
    cx, cy = centroids[st]
    pt_val = pt_means[st]
    ax2.text(cx, cy, f'{st}\nPT={pt_val:.2f}', fontsize=8, fontweight='bold',
             ha='center', va='center', color='white',
             bbox=dict(boxstyle='round,pad=0.25', fc=subtype_palette[st],
                       alpha=0.9, lw=0.5, edgecolor='white'))

cbar = fig.colorbar(sc_pt, ax=ax2, shrink=0.7, pad=0.02)
cbar.set_label('DPT Pseudotime', fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax2.set_title('DPT Pseudotime Projection\n(red arrows = differentiation direction)',
              fontsize=12, fontweight='bold')
ax2.set_xlabel('UMAP1', fontsize=10)
ax2.set_ylabel('UMAP2', fontsize=10)
ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax2.spines.values():
    sp.set_visible(False)

plt.tight_layout()
plt.savefig('trajectory/monomac_trajectory_directed.png',
            dpi=150, bbox_inches='tight')
plt.savefig('trajectory/monomac_trajectory_directed.svg',
            bbox_inches='tight')
plt.close()
print("Directed trajectory figure saved")


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# ── Load data ──────────────────────────────────────────────────────────────────

import scanpy as sc
import h5py
import shutil

src = translate(r'D:\bulk-download\GSE182434\trajectory\adata_mac_annotated.h5ad')
dst = translate(r'D:\bulk-download\GSE182434\trajectory\adata_mac_annotated_patched.h5ad')

shutil.copy2(src, dst)

with h5py.File(dst, 'a') as f:
    if 'uns/log1p' in f and 'base' in f['uns/log1p']:
        del f['uns/log1p/base']
        print("已修复 uns/log1p/base")

print("完成，使用 patched 文件读取")


adata_mac = sc.read_h5ad(translate(r'D:\bulk-download\GSE182434\trajectory\adata_mac_annotated_patched.h5ad'))
print(adata_mac)

subtype_palette = {
    'Mono':    '#3498DB',
    'DC_1':    '#E74C3C',
    'LA_TAM':  '#2ECC71',
    'IFN_TAM': '#E67E22',
    'DC_2':    '#9B59B6',
}
subtypes_order = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']

umap_coords  = adata_mac.obsm['X_umap_paga']
subtypes_arr = adata_mac.obs['mac_subtype'].values
pt_arr       = adata_mac.obs['dpt_pseudotime'].values
conn         = adata_mac.uns['paga']['connectivities'].toarray()
paga_pos     = adata_mac.uns['paga']['pos']

n_cells   = {st: (subtypes_arr == st).sum() for st in subtypes_order}
centroids = {st: umap_coords[subtypes_arr == st].mean(axis=0) for st in subtypes_order}
pt_means  = adata_mac.obs.groupby('mac_subtype')['dpt_pseudotime'].mean()

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: 3-panel PAGA trajectory
# ══════════════════════════════════════════════════════════════════════════════

fig = plt.figure(figsize=(20, 6.5))
gs  = fig.add_gridspec(1, 3, wspace=0.35)

# ─── Panel A: PAGA topology ────────────────────────────────────────────────────

ax_paga = fig.add_subplot(gs[0])

node_sizes        = np.array([n_cells[st] for st in subtypes_order])
node_sizes_scaled = (node_sizes / node_sizes.max()) * 1200 + 300

threshold = 0.15
for i, st_i in enumerate(subtypes_order):
    for j, st_j in enumerate(subtypes_order):
        if j <= i: continue
        w = conn[i, j]
        if w < threshold: continue
        x0, y0 = paga_pos[i]; x1, y1 = paga_pos[j]
        ax_paga.plot([x0, x1], [y0, y1], '-', color='#888888',
                     lw=w*8, alpha=min(0.9, 0.3+w*0.7), zorder=1)
        ax_paga.text((x0+x1)/2, (y0+y1)/2, f'{w:.2f}', fontsize=7,
                     ha='center', va='center', color='#555',
                     bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.7, lw=0))

for i, st in enumerate(subtypes_order):
    x, y = paga_pos[i]
    ax_paga.scatter(x, y, s=node_sizes_scaled[i], c=subtype_palette[st],
                    zorder=3, edgecolors='white', linewidths=2)
    ax_paga.text(x, y+0.08, st, ha='center', va='bottom', fontsize=10,
                 fontweight='bold', color=subtype_palette[st],
                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8, lw=0))
    ax_paga.text(x, y-0.08, f'n={n_cells[st]}', ha='center', va='top',
                 fontsize=8, color='#555')

ax_paga.set_title('PAGA Topology Graph\n(edge weight = connectivity)',
                  fontsize=12, fontweight='bold')
ax_paga.set_xlim(paga_pos[:,0].min()-0.3, paga_pos[:,0].max()+0.3)
ax_paga.set_ylim(paga_pos[:,1].min()-0.3, paga_pos[:,1].max()+0.3)
ax_paga.axis('off')

# ─── Panel B: subtype UMAP ────────────────────────────────────────────────────

ax_sub = fig.add_subplot(gs[1])

for st in subtypes_order:
    mask = subtypes_arr == st
    ax_sub.scatter(umap_coords[mask,0], umap_coords[mask,1],
                   c=subtype_palette[st], s=18, alpha=0.75,
                   linewidths=0, label=st, rasterized=True)

for st in subtypes_order:
    cx, cy = centroids[st]
    ax_sub.text(cx, cy, st, fontsize=14, fontweight='bold',
                ha='center', va='center', color='white',
                bbox=dict(boxstyle='round,pad=0.25', fc=subtype_palette[st], alpha=0.85, lw=0))

ax_sub.set_title('PAGA-initialized UMAP\n(Cell Subtypes)', fontsize=12, fontweight='bold')
ax_sub.set_xlabel('UMAP1', fontsize=12); ax_sub.set_ylabel('UMAP2', fontsize=12)
ax_sub.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_sub.spines.values(): sp.set_visible(False)

legend_handles = [mpatches.Patch(color=subtype_palette[st], label=st) for st in subtypes_order]
ax_sub.legend(handles=legend_handles, fontsize=10, loc='lower right',
              framealpha=0.85, edgecolor='none', title='Subtype', title_fontsize=12)

# ─── Panel C: DPT pseudotime ──────────────────────────────────────────────────

ax_pt = fig.add_subplot(gs[2])

sc_pt = ax_pt.scatter(umap_coords[:,0], umap_coords[:,1],
                      c=pt_arr, cmap='viridis', s=18, alpha=0.85,
                      linewidths=0, rasterized=True)

trajectory_edges = [
    ('Mono',    'IFN_TAM'),
    ('IFN_TAM', 'LA_TAM'),
    ('Mono',    'DC_2'),
    ('DC_2',    'DC_1'),
]

for src, tgt in trajectory_edges:
    x0, y0 = centroids[src]
    x1, y1 = centroids[tgt]
    ax_pt.annotate('', xy=(x1, y1), xytext=(x0, y0),
                   arrowprops=dict(arrowstyle='->', color='white',
                                   lw=2.5, mutation_scale=18))
    ax_pt.annotate('', xy=(x1, y1), xytext=(x0, y0),
                   arrowprops=dict(arrowstyle='->', color='#333333',
                                   lw=1.5, mutation_scale=16))

for st in subtypes_order:
    cx, cy = centroids[st]
    ax_pt.text(cx, cy, st, fontsize=8, fontweight='bold',
               ha='center', va='center', color='white',
               bbox=dict(boxstyle='round,pad=0.2', fc=subtype_palette[st], alpha=0.85, lw=0))

cbar = fig.colorbar(sc_pt, ax=ax_pt, shrink=0.7, pad=0.02)
cbar.set_label('Pseudotime\n(DPT)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax_pt.set_title('PAGA-initialized UMAP\n(DPT Pseudotime)', fontsize=12, fontweight='bold')
ax_pt.set_xlabel('UMAP1', fontsize=10); ax_pt.set_ylabel('UMAP2', fontsize=10)
ax_pt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax_pt.spines.values(): sp.set_visible(False)

OUT_DIR = translate(r'D:\bulk-download\GSE182434\trajectory')

plt.savefig(fr'{OUT_DIR}\monomac_paga_trajectory_rawpt.png',
            dpi=150, bbox_inches='tight')
plt.savefig(fr'{OUT_DIR}\monomac_paga_trajectory_rawpt.svg',
            bbox_inches='tight')
plt.close()
print("3-panel PAGA figure saved")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: 2-panel directed trajectory
# ══════════════════════════════════════════════════════════════════════════════

threshold = 0.25
directed_edges = []
for i, st_i in enumerate(subtypes_order):
    for j, st_j in enumerate(subtypes_order):
        if j <= i: continue
        w = conn[i, j]
        if w < threshold: continue
        if pt_means[st_i] < pt_means[st_j]:
            directed_edges.append((st_i, st_j, w))
        else:
            directed_edges.append((st_j, st_i, w))

print("Directed edges (src → tgt, weight):")
for src, tgt, w in sorted(directed_edges, key=lambda x: -x[2]):
    print(f"  {src} → {tgt}  (w={w:.3f}, ΔPT={pt_means[tgt]-pt_means[src]:.3f})")

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

# ─── Panel A: subtype colors + directed arrows ────────────────────────────────

ax = axes[0]

for st in subtypes_order:
    mask = subtypes_arr == st
    ax.scatter(umap_coords[mask,0], umap_coords[mask,1],
               c=subtype_palette[st], s=20, alpha=0.65,
               linewidths=0, rasterized=True)

for src, tgt, w in directed_edges:
    x0, y0 = centroids[src]; x1, y1 = centroids[tgt]
    dx, dy = x1-x0, y1-y0
    shrink = 0.22
    xs, ys = x0+dx*shrink, y0+dy*shrink
    xe, ye = x1-dx*shrink, y1-dy*shrink
    lw = 1.5 + w*4
    ax.annotate('', xy=(xe,ye), xytext=(xs,ys),
                arrowprops=dict(arrowstyle='->', color='white', lw=lw+1.5, mutation_scale=20))
    ax.annotate('', xy=(xe,ye), xytext=(xs,ys),
                arrowprops=dict(arrowstyle='->', color='#333333', lw=lw, mutation_scale=18))
    mx, my = (xs+xe)/2, (ys+ye)/2
    ax.text(mx, my, f'{w:.2f}', fontsize=12, ha='center', va='center',
            color='#333', bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.8, lw=0))

for st in subtypes_order:
    cx, cy = centroids[st]
    n = (subtypes_arr == st).sum()
    ax.text(cx, cy, f'{st}\n(n={n})', fontsize=14, fontweight='bold',
            ha='center', va='center', color='white',
            bbox=dict(boxstyle='round,pad=0.3', fc=subtype_palette[st],
                      alpha=0.9, lw=0.5, edgecolor='white'))

ax.set_title('Mono/Mac Differentiation Trajectory\n(PAGA + DPT, arrows = direction)',
             fontsize=12, fontweight='bold')
ax.set_xlabel('UMAP1', fontsize=10); ax.set_ylabel('UMAP2', fontsize=10)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax.spines.values(): sp.set_visible(False)

# ─── Panel B: DPT gradient + 红色箭头 ────────────────────────────────────────

ax2 = axes[1]

sc_pt = ax2.scatter(umap_coords[:,0], umap_coords[:,1],
                    c=pt_arr, cmap='viridis', s=20, alpha=0.85,
                    linewidths=0, rasterized=True, vmin=0, vmax=1)

for src, tgt, w in directed_edges:
    x0, y0 = centroids[src]; x1, y1 = centroids[tgt]
    dx, dy = x1-x0, y1-y0
    shrink = 0.22
    xs, ys = x0+dx*shrink, y0+dy*shrink
    xe, ye = x1-dx*shrink, y1-dy*shrink
    ax2.annotate('', xy=(xe,ye), xytext=(xs,ys),
                 arrowprops=dict(arrowstyle='->', color='white', lw=3.5, mutation_scale=20))
    ax2.annotate('', xy=(xe,ye), xytext=(xs,ys),
                 arrowprops=dict(arrowstyle='->', color='#ff4444', lw=2.0, mutation_scale=18))

for st in subtypes_order:
    cx, cy = centroids[st]
    pt_val = pt_means[st]
    ax2.text(cx, cy, f'{st}\nPT={pt_val:.2f}', fontsize=14, fontweight='bold',
             ha='center', va='center', color='white',
             bbox=dict(boxstyle='round,pad=0.25', fc=subtype_palette[st],
                       alpha=0.9, lw=0.5, edgecolor='white'))

cbar = fig.colorbar(sc_pt, ax=ax2, shrink=0.7, pad=0.02)
cbar.set_label('DPT Pseudotime', fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax2.set_title('DPT Pseudotime Projection\n(red arrows = differentiation direction)',
              fontsize=12, fontweight='bold')
ax2.set_xlabel('UMAP1', fontsize=10); ax2.set_ylabel('UMAP2', fontsize=10)
ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
for sp in ax2.spines.values(): sp.set_visible(False)

plt.tight_layout()
plt.savefig(fr'{OUT_DIR}\monomac_trajectory_directed_rawpt.png',
            dpi=300, bbox_inches='tight')
plt.savefig(fr'{OUT_DIR}\monomac_trajectory_directed_rawpt.svg',
            bbox_inches='tight')
plt.close()
print("2-panel directed trajectory figure saved")
