# -*- coding: utf-8 -*-
"""Fig1b marker dotplot - 复刻原版样式(Reds/基因在Y/细胞类型在X)，仅改尺寸与竖直图例"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import translate

OUT = os.path.join('results', 'uav_panels', 'fig1b_v2')
os.makedirs(OUT, exist_ok=True)

CT_ORDER = ['CD8 T cells', 'CD4 T cells', 'Tregs', 'TFH cells', 'NK cells',
            'Naive B cells', 'GC B cells', 'Proliferative GC B cells',
            'Memory B cells', 'Age-associated B cells', 'B cells (other)',
            'Monocytes/Macrophages', 'Plasma cells', 'pDC/Other']

adata = sc.read_h5ad(translate('D:/bulk-download/GSE182434/adata_processed.h5ad'))
markers = pd.read_csv('D:/bulk-download/GSE182434/comparison/cluster_markers.csv')

top = markers.groupby('group').head(2).groupby('group')['names'].apply(list).to_dict()
cts = [c for c in CT_ORDER if c in top]
genes, groups = [], []
for c in cts:
    gs = [g for g in top[c] if g not in genes]
    genes += gs
    groups.append((c, len(gs)))
print('types:', len(cts), 'genes:', len(genes))

X = adata[:, genes]
M = X.X.toarray() if hasattr(X.X, 'toarray') else np.asarray(X.X)
cats = adata.obs['cell_type'].astype('category')
means = np.zeros((len(cts), len(genes)))
fracs = np.zeros_like(means)
for i, c in enumerate(cts):
    m = (cats == c).values
    if m.sum() == 0:
        continue
    sub = M[m]
    means[i] = sub.mean(axis=0)
    fracs[i] = (sub > 0).mean(axis=0)
gmin, gmax = means.min(axis=0, keepdims=True), means.max(axis=0, keepdims=True)
scaled = (means - gmin) / np.maximum(gmax - gmin, 1e-9)

FS = 6.0
fig = plt.figure(figsize=(3.6, 3.1))
ax = fig.add_axes([0.13, 0.165, 0.665, 0.80])  # ★ 加宽7% 加高4%
cmap = plt.get_cmap('Reds')
for i in range(len(cts)):          # x = 细胞类型
    for j in range(len(genes)):    # y = 基因
        f = fracs[i, j]; v = scaled[i, j]
        if f <= 0:
            continue
        ax.scatter(i, -j, s=10 + 40 * f, c=[cmap(0.08 + 0.92 * v)],
                   edgecolors='#666666', linewidths=0.25, zorder=3)
ax.set_xlim(-0.7, len(cts) - 0.3)
ax.set_ylim(-len(genes) + 0.6, 0.7)
ax.set_yticks([-j for j in range(len(genes))])
ax.set_yticklabels(genes, fontsize=FS)
ax.tick_params(axis='y', length=0, pad=1)
ax.set_xticks(range(len(cts)))
ax.set_xticklabels([c.replace(' cells', '').replace('Monocytes/Macrophages', 'Mono/Mac')
                    .replace('Proliferative GC B', 'Prolif. GC B')
                    .replace('Age-associated B', 'Age-assoc. B') for c in cts],
                   rotation=45, ha='right', fontsize=FS)  # ★ 45° 省高度
ax.tick_params(axis='x', length=0, pad=1)
for sp in ax.spines.values():
    sp.set_visible(False)

lax = fig.add_axes([0.785, 0.60, 0.20, 0.175])  # ★ 右移 ~3mm
lax.set_xlim(0, 1); lax.set_ylim(0, 1)
lax.axis('off')
lax.text(0.10, 1.05, 'Fraction (%)', fontsize=FS * 0.85, ha='left', va='bottom')  # ★ 左对齐
for k, f in enumerate([0.25, 0.6, 1.0]):
    lax.scatter(0.18, 0.72 - k * 0.28, s=10 + 40 * f, c='0.92',
                edgecolors='#666666', linewidths=0.25)
    lax.text(0.38, 0.72 - k * 0.28, str(int(f * 100)) + '%', fontsize=FS * 0.85, va='center')

cax = fig.add_axes([0.807, 0.30, 0.028, 0.20])  # ★ 右移 0.3cm 定位
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1)),
                  cax=cax, orientation='vertical')
cb.set_ticks([0, 0.5, 1.0])
cb.set_label('Mean\nexpression', fontsize=FS * 0.85)  # ★ 两短行居中于色条
cb.ax.tick_params(labelsize=FS * 0.8)
cb.outline.set_linewidth(0.4)

fig.savefig(os.path.join(OUT, 'fig1b.pdf'), format='pdf')
print('saved', os.path.join(OUT, 'fig1b.pdf'))
