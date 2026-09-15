# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# =============================================================================
# GSE182434 单细胞 RNA-seq 分析流程
# 数据集：DLBCL（弥漫性大B细胞淋巴瘤）vs 扁桃体（Tonsil）
# 流程：数据下载 → QC过滤 → 归一化 → HVG → PCA → Harmony批次校正
#        → 聚类 → UMAP → CellTypist注释 → 组成分析 → 差异基因
#
# 运行环境：Windows 10，Python 3.9+
# 依赖安装（在cmd/PowerShell中执行，运行代码前完成）：
#   pip install GEOparse scanpy anndata scipy pandas numpy matplotlib
#   pip install harmonypy celltypist scrublet leidenalg
# =============================================================================

import os
import gc
import warnings
import urllib.request

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import anndata as ad
import harmonypy as hm
import celltypist
from celltypist import models
import scrublet as scr
import GEOparse

warnings.filterwarnings('ignore')
matplotlib.use('Agg')  # 非交互式后端，适合保存图片而不弹窗

# =============================================================================
# 0. 路径配置（Windows 本地路径，按需修改）
# =============================================================================

DATA_DIR   = translate(r"D:\scRNA\GSE182434\data")
RESULT_DIR = translate(r"D:\scRNA\GSE182434\results")

QC_DIR    = os.path.join(RESULT_DIR, "qc")
UMAP_DIR  = os.path.join(RESULT_DIR, "umap")
ANNOT_DIR = os.path.join(RESULT_DIR, "annotation")
COMP_DIR  = os.path.join(RESULT_DIR, "comparison")
TMP_DIR   = os.path.join(RESULT_DIR, "tmp")

for d in [DATA_DIR, QC_DIR, UMAP_DIR, ANNOT_DIR, COMP_DIR, TMP_DIR]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# 全局常量：细胞类型顺序与配色（T细胞 → B细胞 → 固有免疫）
# =============================================================================

CT_ORDER = [
    'CD8 T cells', 'CD4 T cells', 'Tregs', 'TFH cells', 'NK cells',
    'Naive B cells', 'GC B cells', 'Proliferative GC B cells',
    'Memory B cells', 'Age-associated B cells', 'B cells (other)',
    'Monocytes/Macrophages', 'Plasma cells', 'pDC/Other'
]

CT_PALETTE = {
    'CD8 T cells':               '#1f77b4',
    'CD4 T cells':               '#aec7e8',
    'Tregs':                     '#9467bd',
    'TFH cells':                 '#c5b0d5',
    'NK cells':                  '#17becf',
    'Naive B cells':             '#2ca02c',
    'GC B cells':                '#98df8a',
    'Proliferative GC B cells':  '#d62728',
    'Memory B cells':            '#ff9896',
    'Age-associated B cells':    '#8c564b',
    'B cells (other)':           '#c49c94',
    'Monocytes/Macrophages':     '#e377c2',
    'Plasma cells':              '#ff7f0e',
    'pDC/Other':                 '#7f7f7f',
}

TISSUE_COLORS  = {'DLBCL': '#0279EE', 'Tonsil': '#FF9400'}
PATIENT_COLORS = {
    'DLBCL002': '#1f77b4', 'DLBCL007': '#ff7f0e',
    'DLBCL008': '#2ca02c', 'DLBCL111': '#d62728', 'T2': '#9467bd'
}

# =============================================================================
# 1. 下载 GEO 元数据与补充文件
# =============================================================================

print("=" * 60)
print("Step 1: 下载 GEO 元数据与数据文件")
print("=" * 60)

gse = pd.read_csv( translate("d:/BadiduNetdiskDownload/R/PTL_M1_M2/GSE182434_raw_count_matrix.txt.gz"), sep="\t", index_col=0, compression="gzip")
# 兼容性修正：此处读入的是表达矩阵 DataFrame，不是 GEOparse 对象，
# 原代码访问 .metadata / .gsms / .gpls 会 AttributeError。
print(f"表达矩阵: {getattr(gse, 'shape', 'N/A')}")

SUPPL_URLS = [
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182434/suppl/GSE182434_raw_count_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182434/suppl/GSE182434_cell_annotation.txt.gz"
]

for url in SUPPL_URLS:
    fname = url.split("/")[-1]
    dest  = os.path.join(DATA_DIR, fname)
    if not os.path.exists(dest):
        print(f"正在下载 {fname} ...")
        # Windows 下使用 Python 内置 urllib 替代 wget
        urllib.request.urlretrieve(url, dest)
        print(f"  ✓ 下载完成: {os.path.getsize(dest)/1e6:.1f} MB")
    else:
        print(f"  ✓ 文件已存在: {fname} ({os.path.getsize(dest)/1e6:.1f} MB)")

# =============================================================================
# 2. 加载数据，构建 AnnData 对象
# =============================================================================

print("\n" + "=" * 60)
print("Step 2: 加载数据，构建 AnnData")
print("=" * 60)

# 加载细胞注释（列包含 Tissue / Patient / CellType / Sample）
annot = pd.read_csv(
    os.path.join(DATA_DIR, "GSE182434_cell_annotation.txt.gz"),
    sep="\t", index_col=0
)
print(f"注释表维度: {annot.shape}  列: {annot.columns.tolist()}")

# 仅保留 DLBCL 和 Tonsil 两组
annot_filtered = annot[annot['Tissue'].isin(['DLBCL', 'Tonsil'])].copy()
print(f"过滤后细胞数: {len(annot_filtered)}")
print(annot_filtered['Tissue'].value_counts().to_string())

# 加载原始计数矩阵（行=基因，列=细胞）
print("\n正在加载计数矩阵（可能需要几分钟）...")
counts = pd.read_csv(
    os.path.join(DATA_DIR, "GSE182434_raw_count_matrix.txt.gz"),
    sep="\t", index_col=0
)
print(f"计数矩阵维度 (基因 × 细胞): {counts.shape}")

# 取目标细胞交集
target_cells = annot_filtered.index.tolist()
overlap      = set(target_cells) & set(counts.columns)
print(f"目标细胞: {len(target_cells)}  矩阵中细胞: {len(counts.columns)}  交集: {len(overlap)}")

# 子集矩阵，转置为细胞×基因，构建稀疏矩阵
counts_sub = counts[target_cells]
X = sp.csr_matrix(counts_sub.values.T)

adata = ad.AnnData(
    X=X,
    obs=annot_filtered.loc[target_cells].copy(),
    var=pd.DataFrame(index=counts_sub.index)
)
adata.var_names_make_unique()

# 用 Patient 作为批次键（Harmony 批次校正所需）
adata.obs['batch'] = adata.obs['Patient'].astype(str)
# 保存原始计数，后续步骤不会覆盖
adata.layers['counts'] = adata.X.copy()

del counts, counts_sub
gc.collect()

print(f"\nAnnData 构建完成: {adata.shape}")
print(f"批次分布:\n{adata.obs['batch'].value_counts().to_string()}")

# =============================================================================
# 3. QC 质控：计算指标 → MAD 异常值检测 → Scrublet 双细胞检测 → 过滤
# =============================================================================

print("\n" + "=" * 60)
print("Step 3: QC 质控")
print("=" * 60)

# 标记线粒体基因（人类线粒体基因以 MT- 开头）
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True
)
print("QC 指标统计:")
print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe().round(2))

def batch_mad_outlier(adata, batch_key, metrics, nmads=5):
    """
    批次感知 MAD 异常值检测：
    对每个批次独立计算中位数和 MAD，标记超出 nmads 倍 MAD 的细胞。
    线粒体比例仅检测上界（高线粒体 = 细胞受损）。
    """
    outlier = pd.Series(False, index=adata.obs_names)
    for batch in adata.obs[batch_key].unique():
        mask = adata.obs[batch_key] == batch
        for metric in metrics:
            vals = adata.obs.loc[mask, metric]
            med  = vals.median()
            mad  = (vals - med).abs().median()
            hi   = med + nmads * mad
            lo   = med - nmads * mad
            if 'mt' in metric:
                outlier[mask] |= (vals > hi)
            else:
                outlier[mask] |= (vals < lo) | (vals > hi)
    adata.obs['outlier'] = outlier
    return adata

adata = batch_mad_outlier(
    adata, batch_key='batch',
    metrics=['log1p_total_counts', 'log1p_n_genes_by_counts', 'pct_counts_mt'],
    nmads=5
)
print(f"MAD 异常值: {adata.obs['outlier'].sum()} / {adata.n_obs} "
      f"({100 * adata.obs['outlier'].mean():.1f}%)")

# Scrublet 双细胞检测（按批次独立运行）
predicted_doublet = pd.Series(False, index=adata.obs_names)
doublet_score     = pd.Series(0.0,   index=adata.obs_names)

for batch in adata.obs['batch'].unique():
    mask  = adata.obs['batch'] == batch
    X_raw = adata.layers['counts'][mask]
    scrub  = scr.Scrublet(X_raw, expected_doublet_rate=0.06)
    scores, calls = scrub.scrub_doublets(
        min_counts=2, min_cells=3, n_prin_comps=30, verbose=False
    )
    predicted_doublet[mask] = calls
    doublet_score[mask]     = scores

adata.obs['predicted_doublet'] = predicted_doublet.values
adata.obs['doublet_score']     = doublet_score.values
print(f"双细胞: {adata.obs['predicted_doublet'].sum()} "
      f"({100 * adata.obs['predicted_doublet'].mean():.1f}%)")

# 过滤前绘制 QC 小提琴图
sc.settings.figdir = QC_DIR
sc.pl.violin(
    adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='batch', rotation=45,
    save="_before_filter.png", show=False
)

# 过滤：去除 MAD 异常值 + 双细胞（score > 0.25）
n_before = adata.n_obs
keep = (~adata.obs['outlier'] &
        ~adata.obs['predicted_doublet'] &
        (adata.obs['doublet_score'] <= 0.25))
adata = adata[keep].copy()
print(f"\n过滤结果: {n_before} → {adata.n_obs} 细胞 "
      f"(保留 {100 * adata.n_obs / n_before:.1f}%)")
print(adata.obs['Tissue'].value_counts().to_string())

# =============================================================================
# 4. 归一化 → 高变基因 → 缩放 → PCA
# =============================================================================

print("\n" + "=" * 60)
print("Step 4: 归一化、HVG 筛选、PCA")
print("=" * 60)

# 每细胞总计数标准化到 10000，再取 log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# 保存 log 归一化数据，供 CellTypist / 差异基因分析使用
adata.layers['log1p_norm'] = adata.X.copy()
print("✓ 归一化完成，log1p_norm layer 已保存")

# 批次感知高变基因筛选（seurat_v3：每批次独立计算方差后取并集）
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor='seurat_v3',
    batch_key='batch', layer='log1p_norm'
)
print(f"✓ 高变基因: {adata.var['highly_variable'].sum()} 个")

# 回归去除技术变异（总计数 + 线粒体比例），再缩放（截断最大值为 10）
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# PCA（50 个主成分，仅对高变基因）
sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, save="_variance.png", show=False)
print("✓ PCA 完成")

# =============================================================================
# 5. Harmony 批次校正
# =============================================================================

print("\n" + "=" * 60)
print("Step 5: Harmony 批次校正")
print("=" * 60)

# 取前 30 个 PC 输入 Harmony
pca_30 = adata.obsm['X_pca'][:, :30]
harmony_out = hm.run_harmony(
    pca_30, adata.obs, 'batch',
    theta=2.0, max_iter_harmony=10, random_state=42, verbose=True
)

# Z_corr 形状为 (n_pcs, n_cells)，转置后存入 obsm
adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
print(f"✓ Harmony 完成，X_pca_harmony 形状: {adata.obsm['X_pca_harmony'].shape}")

adata.uns['harmony_params'] = {
    'batch_key': 'batch', 'theta': 2.0,
    'n_pcs': 30, 'n_batches': adata.obs['batch'].nunique()
}

# =============================================================================
# 6. KNN 图 → Leiden 聚类 → UMAP
# =============================================================================

print("\n" + "=" * 60)
print("Step 6: 构建邻域图、Leiden 聚类、UMAP 降维")
print("=" * 60)

# 基于 Harmony 校正后的表示构建 KNN 图
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca_harmony', random_state=42)
print("✓ KNN 图构建完成")

# 多分辨率 Leiden 聚类（同时运行 4 个分辨率，便于后续选择）
for res in [0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', random_state=42)
    print(f"  leiden_{res}: {adata.obs[f'leiden_{res}'].nunique()} 个簇")

sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=42)
print(f"✓ UMAP 完成，X_umap 形状: {adata.obsm['X_umap'].shape}")

# =============================================================================
# 7. UMAP 可视化（概览：组织 / 批次 / 聚类 / 细胞类型）
# =============================================================================

print("\n" + "=" * 60)
print("Step 7: UMAP 概览图")
print("=" * 60)

umap = adata.obsm['X_umap']

def scatter_by_group(ax, labels, palette, title):
    """通用分组散点图辅助函数"""
    for label, color in palette.items():
        mask = labels == label
        if mask.sum() > 0:
            ax.scatter(umap[mask, 0], umap[mask, 1],
                       c=color, s=2, alpha=0.4, rasterized=True, label=label)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
    ax.legend(markerscale=4, frameon=False, fontsize=9)
    ax.set_aspect('equal')

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('GSE182434: DLBCL vs Tonsil scRNA-seq\n(Harmony 批次校正)',
             fontsize=14, fontweight='bold')

scatter_by_group(axes[0, 0], adata.obs['Tissue'],  TISSUE_COLORS,  'Tissue')
scatter_by_group(axes[0, 1], adata.obs['Patient'], PATIENT_COLORS, 'Patient (Batch)')

# Leiden 0.6 聚类
ax = axes[1, 0]
clusters  = adata.obs['leiden_0.6'].astype(str)
unique_cl = sorted(clusters.unique(), key=int)
cmap_tab  = plt.cm.get_cmap('tab20', len(unique_cl))
for i, cl in enumerate(unique_cl):
    mask = clusters == cl
    ax.scatter(umap[mask, 0], umap[mask, 1],
               c=[cmap_tab(i)], s=2, alpha=0.5, rasterized=True, label=cl)
ax.set_title('Leiden 聚类 (res=0.6)', fontsize=12)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(markerscale=4, frameon=False, fontsize=7, ncol=2, title='Cluster')
ax.set_aspect('equal')

# GEO 原始细胞类型注释
geo_ct_colors = {
    'B cells': '#1f77b4', 'T cells CD8': '#d62728', 'T cells CD4': '#ff7f0e',
    'Tregs': '#9467bd', 'TFH': '#8c564b', 'Monocytes and Macrophages': '#e377c2',
    'NK cells': '#7f7f7f', 'Plasma cells': '#bcbd22', 'Others': '#17becf'
}
scatter_by_group(axes[1, 1], adata.obs['CellType'], geo_ct_colors, 'Cell Type (GEO 注释)')

plt.tight_layout()
for fmt in ['png', 'svg']:
    plt.savefig(os.path.join(UMAP_DIR, f'umap_overview.{fmt}'),
                dpi=150, bbox_inches='tight', format=fmt)
plt.close()
print("✓ UMAP 概览图已保存")

# =============================================================================
# 8. CellTypist 自动注释 + 精细化细胞类型
# =============================================================================

print("\n" + "=" * 60)
print("Step 8: CellTypist 细胞类型注释")
print("=" * 60)

# 下载免疫细胞模型（首次运行需联网，后续自动使用缓存）
models.download_models(model='Immune_All_Low.pkl', force_update=False)

# CellTypist 要求 log1p 归一化（target_sum=10000）的数据
adata_ct = adata.copy()
adata_ct.X = adata_ct.layers['log1p_norm'].copy()

predictions = celltypist.annotate(
    adata_ct, model='Immune_All_Low.pkl', majority_voting=True
)
del adata_ct

adata.obs['celltypist_label']    = predictions.predicted_labels['predicted_labels'].values
adata.obs['celltypist_majority'] = predictions.predicted_labels['majority_voting'].values
adata.obs['celltypist_conf']     = predictions.probability_matrix.max(axis=1).values
print("Top 15 预测细胞类型:")
print(adata.obs['celltypist_majority'].value_counts().head(15).to_string())

# B 细胞亚型映射（CellTypist 标签 → 统一命名）
B_SUBTYPE_MAP = {
    'Naive B cells':                         'Naive B cells',
    'Germinal center B cells':               'GC B cells',
    'Proliferative germinal center B cells': 'Proliferative GC B cells',
    'Memory B cells':                        'Memory B cells',
    'Age-associated B cells':                'Age-associated B cells',
    'Plasma cells':                          'Plasma cells',
}

# 其他细胞类型映射（GEO 标签 → 统一命名）
GEO_TO_REFINED = {
    'TFH':                       'TFH cells',
    'Tregs':                     'Tregs',
    'T cells CD8':               'CD8 T cells',
    'T cells CD4':               'CD4 T cells',
    'NK cells':                  'NK cells',
    'Monocytes and Macrophages': 'Monocytes/Macrophages',
    'Plasma cells':              'Plasma cells',
    'Others':                    'pDC/Other',
}

def refine_cell_type(row):
    """
    精细化注释策略：
    - B cells：使用 CellTypist 亚型（生发中心、记忆等）
    - 其他类型：使用 GEO 标注映射到统一命名
    """
    geo = row['CellType']
    ct  = row['celltypist_majority']
    if geo == 'B cells':
        return B_SUBTYPE_MAP.get(ct, 'B cells (other)')
    return GEO_TO_REFINED.get(geo, geo)

adata.obs['cell_type'] = adata.obs.apply(refine_cell_type, axis=1)
print(f"\n精细化细胞类型分布 ({adata.obs['cell_type'].nunique()} 种):")
print(adata.obs['cell_type'].value_counts().to_string())

# =============================================================================
# 9. 细胞类型 UMAP 图 + 标志基因表达图
# =============================================================================

print("\n" + "=" * 60)
print("Step 9: 细胞类型 UMAP 与标志基因可视化")
print("=" * 60)

# --- 细胞类型 UMAP ---
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
fig.suptitle('GSE182434: 细胞类型注释\n(GEO 标签 + CellTypist B 细胞亚型细化)',
             fontsize=13, fontweight='bold')

ax = axes[0]
for ct, color in CT_PALETTE.items():
    mask = adata.obs['cell_type'] == ct
    if mask.sum() > 0:
        ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=2, alpha=0.5,
                   rasterized=True, label=f'{ct} ({mask.sum()})')
ax.set_title('精细化细胞类型', fontsize=12)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(markerscale=5, frameon=False, fontsize=7.5,
          loc='upper left', bbox_to_anchor=(1.01, 1), title='Cell Type (n)')
ax.set_aspect('equal')

ax = axes[1]
for tissue, color in TISSUE_COLORS.items():
    mask = adata.obs['Tissue'] == tissue
    ax.scatter(umap[mask, 0], umap[mask, 1], c=color, s=2, alpha=0.35,
               rasterized=True, label=f'{tissue} (n={mask.sum()})')
ax.set_title('组织（DLBCL vs Tonsil）', fontsize=12)
ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2')
ax.legend(markerscale=5, frameon=False, fontsize=10, title='Tissue')
ax.set_aspect('equal')

plt.tight_layout()
for fmt in ['png', 'svg']:
    plt.savefig(os.path.join(ANNOT_DIR, f'umap_celltypes.{fmt}'),
                dpi=150, bbox_inches='tight', format=fmt)
plt.close()
print("✓ 细胞类型 UMAP 已保存")

# --- 标志基因表达图 ---
MARKERS = {
    'CD3D':  'T cells',       'CD8A':  'CD8 T cells',
    'CD4':   'CD4 T cells',   'FOXP3': 'Tregs',
    'CXCR5': 'TFH',           'MS4A1': 'B cells (CD20)',
    'CD38':  'Plasma/GC B',   'MKI67': 'Proliferating',
    'LYZ':   'Monocytes/Mac', 'GNLY':  'NK cells',
    'IGHG1': 'Plasma cells',  'BCL6':  'GC B cells',
}

X_log      = adata.layers['log1p_norm']
gene_names = adata.var_names.tolist()

fig, axes = plt.subplots(3, 4, figsize=(18, 13))
fig.suptitle('标志基因表达（log 归一化）', fontsize=13, fontweight='bold')

for ax, (gene, label) in zip(axes.flatten(), MARKERS.items()):
    if gene in gene_names:
        idx  = gene_names.index(gene)
        expr = (np.array(X_log[:, idx].todense()).flatten()
                if sp.issparse(X_log) else X_log[:, idx])
        order    = np.argsort(expr)  # 高表达细胞绘制在上层
        sc_plot  = ax.scatter(umap[order, 0], umap[order, 1],
                              c=expr[order], s=1.5, alpha=0.7,
                              cmap='YlOrRd', rasterized=True, vmin=0)
        plt.colorbar(sc_plot, ax=ax, shrink=0.7, pad=0.02)
        ax.set_title(f'{gene}\n({label})', fontsize=9)
    else:
        ax.set_title(f'{gene} (未找到)', fontsize=9)
        ax.text(0.5, 0.5, 'Not in dataset', ha='center', va='center',
                transform=ax.transAxes)
    ax.set_xlabel('UMAP1', fontsize=7); ax.set_ylabel('UMAP2', fontsize=7)
    ax.tick_params(labelsize=6); ax.set_aspect('equal')

plt.tight_layout()
for fmt in ['png', 'svg']:
    plt.savefig(os.path.join(ANNOT_DIR, f'marker_genes.{fmt}'),
                dpi=150, bbox_inches='tight', format=fmt)
plt.close()
print("✓ 标志基因图已保存")

# =============================================================================
# 10. 细胞组成分析（DLBCL vs Tonsil）
# =============================================================================

print("\n" + "=" * 60)
print("Step 10: 细胞类型组成分析")
print("=" * 60)

# 计算每个样本各细胞类型的比例
obs = adata.obs[['Sample', 'Tissue', 'Patient', 'cell_type']].copy()
counts_df = obs.groupby(['Sample', 'cell_type']).size().reset_index(name='n')
totals    = counts_df.groupby('Sample')['n'].sum().rename('total')
counts_df = counts_df.join(totals, on='Sample')
counts_df['proportion'] = counts_df['n'] / counts_df['total']

# 验证每个样本比例之和为 1
check = counts_df.groupby('Sample')['proportion'].sum()
assert (check.round(3) == 1.0).all(), "比例计算错误！"

# 转为宽格式，添加元数据
wide = counts_df.pivot_table(
    index='Sample', columns='cell_type',
    values='proportion', aggfunc='first', fill_value=0
)
meta = obs[['Sample', 'Tissue', 'Patient']].drop_duplicates().set_index('Sample')
wide = wide.join(meta).reset_index()
wide = wide.sort_values(['Tissue', 'Patient']).reset_index(drop=True)

ct_cols   = [c for c in CT_ORDER if c in wide.columns]
n         = len(wide)
n_dlbcl   = (wide['Tissue'] == 'DLBCL').sum()

def sample_label(row):
    """生成简洁的样本标签，区分 B（肿瘤）和 NB（非肿瘤）活检"""
    suffix = 'NB' if row['Sample'].endswith('NB') else 'B'
    return f"{row['Patient']}\n({suffix})"

sample_labels = [sample_label(row) for _, row in wide.iterrows()]

# --- 图1：每样本堆叠柱状图 ---
fig1, ax = plt.subplots(figsize=(12, 5.5))
fig1.suptitle('Per-Sample Cell Type Composition: DLBCL vs Tonsil',
              fontsize=13, fontweight='bold')

x      = np.arange(n)
bottom = np.zeros(n)
for ct in ct_cols:
    vals = wide[ct].values
    ax.bar(x, vals, bottom=bottom, color=CT_PALETTE[ct],
           label=ct, width=0.72, linewidth=0)
    bottom += vals

# 组织背景色
ax.axvspan(-0.5, n_dlbcl - 0.5, alpha=0.07, color='#0279EE', zorder=0)
ax.axvspan(n_dlbcl - 0.5, n - 0.5, alpha=0.07, color='#FF9400', zorder=0)
ax.axvline(x=n_dlbcl - 0.5, color='black', linewidth=1.5, linestyle='--', alpha=0.7)
ax.text((n_dlbcl - 1) / 2, 1.04, 'DLBCL',
        ha='center', fontsize=11, fontweight='bold', color='#0279EE')
ax.text(n_dlbcl + (n - n_dlbcl - 1) / 2, 1.04, 'Tonsil',
        ha='center', fontsize=11, fontweight='bold', color='#FF9400')

ax.set_xticks(x)
ax.set_xticklabels(sample_labels, fontsize=9.5, ha='center')
ax.set_ylabel('Cell Type Proportion', fontsize=11)
ax.set_ylim(0, 1.10); ax.set_xlim(-0.5, n - 0.5)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

handles = [mpatches.Patch(color=CT_PALETTE[ct], label=ct) for ct in ct_cols]
ax.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left',
          fontsize=8.5, frameon=False, title='Cell Type', title_fontsize=9.5)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig1.savefig(os.path.join(COMP_DIR, f'celltype_composition_per_sample.{fmt}'),
                 dpi=150, bbox_inches='tight', format=fmt)
plt.close(fig1)
print("✓ 每样本堆叠柱状图已保存")

# --- 图2：DLBCL vs Tonsil 平均比例对比图 ---
mean_by_tissue = wide.groupby('Tissue')[ct_cols].mean()
mean_pivot     = mean_by_tissue.T.reindex([c for c in CT_ORDER if c in mean_by_tissue.columns])

fig2, ax2 = plt.subplots(figsize=(8, 7))
fig2.suptitle('Mean Cell Type Proportion: DLBCL vs Tonsil',
              fontsize=13, fontweight='bold')

y_pos       = np.arange(len(mean_pivot))
dlbcl_vals  = mean_pivot.get('DLBCL',  pd.Series(0, index=mean_pivot.index)).values
tonsil_vals = mean_pivot.get('Tonsil', pd.Series(0, index=mean_pivot.index)).values

bars_d = ax2.barh(y_pos - 0.2, dlbcl_vals,  height=0.35,
                  color='#0279EE', alpha=0.85, label='DLBCL')
bars_t = ax2.barh(y_pos + 0.2, tonsil_vals, height=0.35,
                  color='#FF9400', alpha=0.85, label='Tonsil')

# 柱端数值标注
for bar, val in zip(bars_d, dlbcl_vals):
    if val > 0.01:
        ax2.text(val + 0.004, bar.get_y() + bar.get_height() / 2,
                 f'{val:.2f}', va='center', fontsize=7.5, color='#0055bb')
for bar, val in zip(bars_t, tonsil_vals):
    if val > 0.01:
        ax2.text(val + 0.004, bar.get_y() + bar.get_height() / 2,
                 f'{val:.2f}', va='center', fontsize=7.5, color='#cc6600')

ax2.set_yticks(y_pos)
ax2.set_yticklabels(mean_pivot.index, fontsize=10)
ax2.set_xlabel('Mean Proportion', fontsize=11)
ax2.legend(fontsize=10, frameon=False, loc='lower right')
ax2.invert_yaxis()
ax2.set_xlim(0, max(dlbcl_vals.max(), tonsil_vals.max()) * 1.3)
for spine in ['top', 'right']:
    ax2.spines[spine].set_visible(False)

plt.tight_layout()
for fmt in ['png', 'svg']:
    fig2.savefig(os.path.join(COMP_DIR, f'celltype_mean_proportion.{fmt}'),
                 dpi=150, bbox_inches='tight', format=fmt)
plt.close(fig2)
print("✓ 平均比例对比图已保存")

# =============================================================================
# 11. 差异基因分析（Wilcoxon，按细胞类型）
# =============================================================================

print("" + "=" * 60)
print("Step 11: 差异基因分析")
print("=" * 60)

# 差异基因分析使用 log1p 归一化数据
adata.X = adata.layers['log1p_norm'].copy()

sc.tl.rank_genes_groups(
    adata, groupby='cell_type', method='wilcoxon',
    pts=True, min_in_group_fraction=0.1,
    min_fold_change=0.25, key_added='rank_genes_celltypes'
)
print("✓ 差异基因分析完成")

markers_df = sc.get.rank_genes_groups_df(
    adata, group=None, key='rank_genes_celltypes',
    pval_cutoff=0.05, log2fc_min=0.25
)
markers_df.to_csv(os.path.join(COMP_DIR, 'cluster_markers.csv'), index=False)
print(f"✓ 差异基因已保存: {len(markers_df)} 个显著基因")

print("每种细胞类型 Top 3 标志基因:")
top3 = markers_df.groupby('group').head(3)[['group', 'names', 'logfoldchanges', 'pvals_adj']]
print(top3.to_string(index=False))

# --- Dot plot：各细胞类型 Top4 标志基因 ---
top_markers = (
    markers_df.groupby('group').head(4)
    .groupby('group')['names'].apply(list).to_dict()
)
ct_order_present = [c for c in CT_ORDER if c in top_markers]
genes_ordered    = []
for ct in ct_order_present:
    genes_ordered.extend(top_markers[ct])

# 去重保留顺序
seen       = set()
genes_plot = [g for g in genes_ordered if not (g in seen or seen.add(g))]
print(f"Dot plot 基因数: {len(genes_plot)}")

sc.settings.figdir = COMP_DIR
with plt.rc_context({'figure.figsize': (14, 7)}):
    dp = sc.pl.dotplot(
        adata, var_names=genes_plot, groupby='cell_type',
        categories_order=ct_order_present,
        standard_scale='var',
        colorbar_title='Mean expr\n(scaled)',
        size_title='Fraction\nexpressing',
        show=False, return_fig=True
    )
    dp.savefig(os.path.join(COMP_DIR, 'marker_dotplot.png'),
               dpi=150, bbox_inches='tight')
    dp.savefig(os.path.join(COMP_DIR, 'marker_dotplot.svg'),
               bbox_inches='tight')
    plt.close('all')
print("✓ Dot plot 已保存")

# =============================================================================
# 12. 保存处理后的 AnnData
# =============================================================================

print(" " + "=" * 60)
print("Step 12: 保存 AnnData")
print("=" * 60)

# 避免 h5ad 保存时斜杠导致路径问题
adata.obs['cell_type_safe'] = adata.obs['cell_type'].str.replace('/', '_', regex=False)

# 差异基因结果已导出为 CSV，uns 中的大型结果可删除以减小文件体积
if 'rank_genes_celltypes' in adata.uns:
    del adata.uns['rank_genes_celltypes']

out_path = os.path.join(RESULT_DIR, 'adata_processed.h5ad')
adata.write_h5ad(out_path, compression='gzip')
print(f"✓ AnnData 已保存: {adata.shape} → {os.path.getsize(out_path)/1e6:.1f} MB")
print(f"obs 列: {adata.obs.columns.tolist()}")
print(f"obsm:   {list(adata.obsm.keys())}")
print(f"layers: {list(adata.layers.keys())}")
