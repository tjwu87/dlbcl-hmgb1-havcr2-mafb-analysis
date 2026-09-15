# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  GSE232853 —— GeoMx Q3 归一化 + Harmony 批次校正 + PCA/UMAP
#                                                     （section5 上游 · 第 2 步）
#
#  来源：Biomni 平台分析记录 `14_空间转录组分析.md`
#        · L13815-13929  Q3 归一化 + HVG + PCA + Harmony + UMAP
#        · L14423-14477  Harmony 输出形状的最终修正（Z_corr 已是 (样本, PC)，**不要转置**）
#        · L14011-14015  导出 UMAP 坐标
#        · L12630-12710  导出 PCA 坐标 / HVG（下游 06_01 读这些 CSV）
#
#  GeoMx Q3 归一化（官方推荐）：
#        X_norm = X / Q3_i * median(Q3)   ，Q3_i = 每个样本非零基因的 75 分位数
#        再做 log1p。
#
#  ⚠ Harmony 踩坑（原文记录）：`ho.Z_corr` 在容器里返回的已经是 (572, 50)，
#     直接 `np.ascontiguousarray` 赋值即可；早期版本误加了 `.T` 导致形状错乱。
#     另外要把 PCA 矩阵转成 C 连续（`np.ascontiguousarray`）才能喂给 harmonypy。
#
#  产出（写入 <数据根>/GSE232853_v2/）：
#        expression_q3_lognorm.csv      572 样本 × 16 560 基因（Q3 + log1p）
#        Q3_normalization_stats.csv     每样本 Q3 值与缩放系数
#        task2_PCA_coordinates.csv      前 20 个 PC + 元数据
#        task2_UMAP_coordinates.csv     UMAP_1/2 + 元数据
#        task2_highly_variable_genes.csv
#        GSE232853_adata_v2.h5ad        供后续步骤复用（落到数据根的 03_中间/ 下）
#
#  环境：pandas / numpy / scanpy / anndata / harmonypy / seaborn / matplotlib
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import harmonypy as hm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

OUT = translate(r'D:\bulk-download\GSE232853_v2')     # -> <数据根>/GSE232853_v2
os.makedirs(OUT, exist_ok=True)

# ── ① 读入原始（已 QC）矩阵 ───────────────────────────────────────────────────
df_raw = pd.read_csv(os.path.join(OUT, 'expression_raw_filtered.csv'), index_col=0)
meta = pd.read_csv(os.path.join(OUT, 'meta_filtered.csv'), index_col='Sample_ID')
expr = df_raw.T                                  # 转成 样本 × 基因
X_raw = expr.values.astype(np.float64)
print(f'Expression: {X_raw.shape[0]} samples x {X_raw.shape[1]} genes')

# ── ② GeoMx Q3 归一化 ─────────────────────────────────────────────────────────
q3_per_sample = np.array([
    np.percentile(row[row > 0], 75) if (row > 0).any() else 1.0
    for row in X_raw
])
median_q3 = np.median(q3_per_sample)
print(f'Q3 range: {q3_per_sample.min():.1f} – {q3_per_sample.max():.1f}')
print(f'Median Q3 (scale factor): {median_q3:.1f}')

X_q3 = X_raw / q3_per_sample[:, np.newaxis] * median_q3
X_log = np.log1p(X_q3).astype(np.float32)
print(f'Q3-normalized log1p range: {X_log.min():.3f} – {X_log.max():.3f}')

pd.DataFrame({'Sample_ID': expr.index,
              'Q3_value': q3_per_sample,
              'scale_factor': median_q3 / q3_per_sample}
             ).to_csv(os.path.join(OUT, 'Q3_normalization_stats.csv'), index=False)

# 全量 Q3-lognorm 矩阵（下游直接读这一份）
pd.DataFrame(X_log, index=expr.index, columns=expr.columns
             ).to_csv(os.path.join(OUT, 'expression_q3_lognorm.csv'))

# ── ③ 构建 AnnData ────────────────────────────────────────────────────────────
adata = ad.AnnData(X=X_log.copy())
adata.obs_names = expr.index.tolist()
adata.var_names = expr.columns.tolist()
adata.layers['q3_lognorm'] = X_log.copy()
adata.layers['raw_counts'] = X_raw.astype(np.float32)
for col in ['Slide_ID', 'ROI_Num', 'Mask', 'Unique_ROI', 'Tissue_Type']:
    adata.obs[col] = meta.loc[adata.obs_names, col].values
print(f'AnnData: {adata.shape}')

# ── ④ HVG + scale + PCA ───────────────────────────────────────────────────────
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
print(f"HVGs: {int(adata.var['highly_variable'].sum())}")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

# ── ⑤ Harmony 批次校正（按 Slide_ID）─────────────────────────────────────────
pca_mat = np.ascontiguousarray(adata.obsm['X_pca'])
ho = hm.run_harmony(pca_mat, adata.obs, 'Slide_ID',
                    max_iter_harmony=20, random_state=42)
# Z_corr 已是 (n_samples, n_pcs)，直接用，不要转置
adata.obsm['X_pca_harmony'] = np.ascontiguousarray(ho.Z_corr)
print(f"Harmony done: {adata.obsm['X_pca_harmony'].shape}")

# ── ⑥ 邻居 + UMAP ─────────────────────────────────────────────────────────────
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep='X_pca_harmony')
sc.tl.umap(adata, min_dist=0.35, spread=1.2, random_state=42)

# ── ⑦ 导出下游所需 CSV ───────────────────────────────────────────────────────
umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names,
                       columns=['UMAP_1', 'UMAP_2'])
umap_df = umap_df.join(adata.obs[['Slide_ID', 'Mask', 'Tissue_Type', 'Unique_ROI']])
umap_df.index.name = 'Sample_ID'
umap_df.to_csv(os.path.join(OUT, 'task2_UMAP_coordinates.csv'))

pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :20], index=adata.obs_names,
                      columns=[f'PC{i+1}' for i in range(20)])
pca_df = pca_df.join(adata.obs[['Slide_ID', 'Mask', 'Tissue_Type', 'Unique_ROI']])
pca_df.index.name = 'Sample_ID'
pca_df.to_csv(os.path.join(OUT, 'task2_PCA_coordinates.csv'))

hvg_df = adata.var[['highly_variable', 'means', 'dispersions',
                    'dispersions_norm']].copy()
hvg_df.index.name = 'gene'
hvg_df.to_csv(os.path.join(OUT, 'task2_highly_variable_genes.csv'))

# ── ⑧ 等高线 UMAP 图（论文 Fig5-a 的底稿，供核对）────────────────────────────
umap = adata.obsm['X_umap']
MASK_COLORS = {'CD20': '#2166AC', 'CD68': '#D6604D'}
TISSUE_COLORS = {'DLBCL': '#B2182B', 'Normal': '#4393C3'}

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.patch.set_facecolor('#F8F8F8')
for ax, color_by, cmap, title in zip(
        axes, ['Mask', 'Tissue_Type'], [MASK_COLORS, TISSUE_COLORS],
        ['Cell Compartment (CD20 / CD68)', 'Tissue Type (DLBCL / Normal)']):
    ax.set_facecolor('#F0F0F0')
    ax.scatter(umap[:, 0], umap[:, 1], c='#DDDDDD', s=6, alpha=0.3,
               linewidths=0, zorder=1, rasterized=True)
    for grp in sorted(adata.obs[color_by].unique()):
        color = cmap.get(grp, '#888888')
        idx = adata.obs[color_by].values == grp
        x_g, y_g = umap[idx, 0], umap[idx, 1]
        ax.scatter(x_g, y_g, c=color, s=14, alpha=0.55, linewidths=0,
                   zorder=2, rasterized=True, label=f'{grp} (n={int(idx.sum())})')
        try:
            sns.kdeplot(x=x_g, y=y_g, ax=ax, levels=5, color=color,
                        linewidths=[0.8, 1.2, 1.6, 2.0, 2.4], alpha=0.85,
                        zorder=4, fill=False)
            sns.kdeplot(x=x_g, y=y_g, ax=ax, levels=3, color=color,
                        alpha=0.08, zorder=3, fill=True)
        except Exception as e:
            print(f'  KDE failed for {grp}: {e}')
    ax.set_xlabel('UMAP 1', fontsize=11)
    ax.set_ylabel('UMAP 2', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
    ax.legend(fontsize=9.5, framealpha=0.9, markerscale=1.8,
              loc='upper right', edgecolor='#CCCCCC')
    ax.tick_params(labelsize=8)
    for sp in ax.spines.values():
        sp.set_edgecolor('#CCCCCC')
fig.text(0.5, -0.02,
         f"Harmony batch correction: {adata.obs['Slide_ID'].nunique()} slides",
         ha='center', fontsize=9, color='#555555', style='italic')
fig.suptitle('GeoMx DSP — Topographic Contour UMAP\n'
             'GSE232853 · DLBCL · Q3-Normalized · 572 ROIs · 16,560 genes',
             fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'task2_contour_UMAP.png'), dpi=150,
            bbox_inches='tight', facecolor='#F8F8F8')
plt.close()

# ── ⑨ 保存 AnnData 供后续步骤复用 ────────────────────────────────────────────
adata.write_h5ad(os.path.join(OUT, 'GSE232853_adata_v2.h5ad'))
print('[OK] Q3 stats / 坐标 / HVG / AnnData 已写出')
