# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  恶性 vs 正常 B 细胞判定（inferCNVpy）      （section1 上游 · 恶性判定）
#
#  来源：Biomni 平台分析记录 `01_对话记录/01_B细胞轨迹与恶性鉴定.md`
#        · L2625-2650  取 B 细胞子集 + 设 cnv_group（Tonsil 为参考）
#        · L2822-2887  用 scanpy biomart 取基因染色体坐标（chr1–chr22 + X）
#        · L3394-3411  染色体命名修正（inferCNVpy 要求 'chr1' 而非 '1'）
#        · L3472-3482  cnv.tl.infercnv 实际参数
#        · L3714-3746  cnv_score + 阈值（Tonsil 99 分位）→ malignancy 标签
#        · L6215-6216  落盘
#
#  产出（写入 <数据根>/GSE182434/cnv/）：
#        malignancy_classification.csv
#            cell_type / Tissue / Patient / cnv_group / cnv_score / malignancy
#        chromosome_heatmap.png          染色体 CNV 热图（供核对）
#
#  ⚠ 原记录的两个坑（照做可避）：
#     1. 基因染色体名必须是 'chr1'/'chr2'…，inferCNVpy 内部用
#        `x.startswith("chr")` 过滤，'1' 会被整条丢掉 → 必须显式加前缀。
#     2. 每条染色体基因数要 ≥ window_size(200)，否则 `_running_mean_by_chromosome`
#        报错；chrY/chrM 等小染色体建议直接排除。
#
#  环境：scanpy（含 sc.queries.biomart_annotations，需联网）/ infercnvpy / anndata
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import warnings
warnings.filterwarnings('ignore')

GSE = translate(r'D:\bulk-download\GSE182434')          # -> <数据根>/GSE182434
CNV_DIR = os.path.join(GSE, 'cnv')
os.makedirs(CNV_DIR, exist_ok=True)

# ── ① 读入已整合/注释的单细胞数据 ─────────────────────────────────────────────
adata = sc.read_h5ad(os.path.join(GSE, 'adata_processed.h5ad'))
print(f'adata: {adata.shape}')

B_TYPES = ['Naive B cells', 'GC B cells', 'Proliferative GC B cells',
           'Memory B cells', 'Age-associated B cells', 'B cells (other)',
           'Plasma cells']
adata_b = adata[adata.obs['cell_type'].isin(B_TYPES)].copy()
print(f'Total B cells: {adata_b.shape[0]}')

# inferCNVpy 要求原始 counts
adata_b.X = adata_b.layers['counts'].copy()

# 参考组：扁桃体（Tonsil）为正常基线
adata_b.obs['cnv_group'] = 'DLBCL_B'
adata_b.obs.loc[adata_b.obs['Tissue'] == 'Tonsil', 'cnv_group'] = 'Tonsil_ref'
print(adata_b.obs['cnv_group'].value_counts().to_string())

# ── ② 基因染色体坐标（Ensembl biomart，带磁盘缓存）───────────────────────────
annot = sc.queries.biomart_annotations(
    'hsapiens',
    ['ensembl_gene_id', 'external_gene_name', 'chromosome_name',
     'start_position', 'end_position'],
    use_cache=True,
)
valid_chroms = [str(i) for i in range(1, 23)] + ['X']
annot_filt = annot[annot['chromosome_name'].isin(valid_chroms)].drop_duplicates(
    'external_gene_name')
gene_pos = annot_filt.set_index('external_gene_name')[
    ['chromosome_name', 'start_position', 'end_position']]

adata_b.var['chromosome'] = gene_pos['chromosome_name'].reindex(adata_b.var_names)
adata_b.var['start'] = gene_pos['start_position'].reindex(adata_b.var_names)
adata_b.var['end'] = gene_pos['end_position'].reindex(adata_b.var_names)
print(f"Annotated genes: {adata_b.var['chromosome'].notna().sum()} / {adata_b.n_vars}")

# 只保留有坐标的基因
adata_cnv = adata_b[:, adata_b.var['chromosome'].notna()].copy()

# 关键修正：加 'chr' 前缀（否则 inferCNVpy 内部会整条过滤掉）
adata_cnv.var['chromosome'] = adata_cnv.var['chromosome'].apply(lambda c: f'chr{c}')
print(f"Chromosome format fixed: {adata_cnv.var['chromosome'].unique()[:5]}")

# ── ③ 推断 CNV ────────────────────────────────────────────────────────────────
print('Running inferCNVpy ...')
cnv.tl.infercnv(
    adata_cnv,
    reference_key='cnv_group',
    reference_cat='Tonsil_ref',
    window_size=200,
    step=10,
    lfc_clip=3.0,
    dynamic_threshold=1.5,
    chunksize=5000,
    n_jobs=1,
)

# ── ④ 每细胞 CNV 得分 ─────────────────────────────────────────────────────────
cnv.tl.cnv_score(adata_cnv)
print(adata_cnv.obs['cnv_score'].describe().to_string())

# ── ⑤ 阈值判定：Tonsil 参考的 99 分位 ─────────────────────────────────────────
dlbcl_scores = adata_cnv.obs.loc[adata_cnv.obs['cnv_group'] == 'DLBCL_B', 'cnv_score']
tonsil_scores = adata_cnv.obs.loc[adata_cnv.obs['cnv_group'] == 'Tonsil_ref', 'cnv_score']
threshold = tonsil_scores.quantile(0.99)
print(f'Threshold (Tonsil 99th pct): {threshold:.4f}')

adata_cnv.obs['malignancy'] = 'Normal B cell'          # Tonsil 全部为正常参考
dlbcl_mask = adata_cnv.obs['cnv_group'] == 'DLBCL_B'
malignant_mask = dlbcl_mask & (adata_cnv.obs['cnv_score'] > threshold)
normal_dlbcl_mask = dlbcl_mask & (adata_cnv.obs['cnv_score'] <= threshold)
adata_cnv.obs.loc[malignant_mask, 'malignancy'] = 'Malignant B cell'
adata_cnv.obs.loc[normal_dlbcl_mask, 'malignancy'] = 'Normal B cell (DLBCL)'
print(adata_cnv.obs['malignancy'].value_counts().to_string())

# ── ⑥ 落盘 ────────────────────────────────────────────────────────────────────
cols = ['cell_type', 'Tissue', 'Patient', 'cnv_group', 'cnv_score', 'malignancy']
cnv_meta = adata_cnv.obs[[c for c in cols if c in adata_cnv.obs.columns]].copy()
cnv_meta.to_csv(os.path.join(CNV_DIR, 'malignancy_classification.csv'))
print(f'[OK] malignancy_classification.csv  ({len(cnv_meta)} cells)')

# 染色体热图（仅作核对，不参与统计）
try:
    import matplotlib
    matplotlib.use('Agg')
    axes_dict = cnv.pl.chromosome_heatmap(
        adata_cnv, groupby='malignancy', show=False)
    fig = list(axes_dict.values())[0].figure
    fig.savefig(os.path.join(CNV_DIR, 'chromosome_heatmap.png'),
                dpi=150, bbox_inches='tight')
    print('[OK] chromosome_heatmap.png')
except Exception as e:
    print(f'[warn] chromosome_heatmap 绘制失败（不影响 CSV）: {e}')
