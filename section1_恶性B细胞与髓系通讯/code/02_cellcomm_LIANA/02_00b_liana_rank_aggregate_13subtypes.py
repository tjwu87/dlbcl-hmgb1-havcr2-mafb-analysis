# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR
# <<< 自动注入结束 >>>

# =============================================================================
#  LIANA 细胞通讯本体（rank_aggregate）—— 论文 section1 上游
#
#  来源：Biomni 平台分析记录（02_细胞通讯LIANA全细胞 / 05_细胞通讯v2_13亚型）
#  说明：本脚本重算 LIANA，产出 02_00_liana_run.py / 02_01 / 02_02 所读的 CSV。
#        运行约需数分钟（14 细胞类型）到十余分钟（13 亚型，含 100 次置换）。
#
#  环境：scanpy + anndata + liana(>=1.1) + pandas + numpy
#        adata 的 .X 必须是 log1p 归一化后的表达（liana 会按 log-normalized 处理）
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import liana
import scipy.sparse as sp
import warnings
warnings.filterwarnings('ignore')

GSE = translate(r'D:\bulk-download\GSE182434')   # -> <数据根>/GSE182434
CELLCOMM = os.path.join(GSE, 'cellcomm')
os.makedirs(CELLCOMM, exist_ok=True)

print(f'LIANA version: {liana.__version__}')
adata = sc.read_h5ad(os.path.join(GSE, 'adata_processed.h5ad'))
print(f'adata: {adata.shape}')

# =============================================================================
#  ② 13 个细胞亚型（Malignant B / Normal B ＋ 髓系细分）—— 论文主图所用
#
#     分类规则（照 Biomni 记录 05_细胞通讯v2_13亚型 原文）：
#       · T/NK/pDC 由 cell_type 直接映射
#       · B 细胞按 malignancy_classification.csv 拆成 Malignant B / Normal B
#       · Monocytes/Macrophages 按 5 个 marker 打分取最大 → Mono/DC_1/DC_2/ IFN_TAM/LA_TAM
# =============================================================================
MAC = os.path.join(GSE, 'cnv', 'malignancy_classification.csv')
mal_labels = pd.read_csv(MAC, index_col=0)
print(f'Malignancy labels: {len(mal_labels)}')

# ── 髓系 5 亚型打分（marker 取均值，argmax 定类）──────────────────────────────
adata_mac = adata[adata.obs['cell_type'] == 'Monocytes/Macrophages'].copy()
subtype_markers = {
    'Mono':    ['FCN1', 'S100A9', 'S100A8', 'S100A4', 'APOBEC3A'],
    'DC_1':    ['LTB', 'CLEC10A', 'CD1C', 'JAML', 'CD1E'],
    'LA_TAM':  ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD'],
    'IFN_TAM': ['MT1H', 'MT1G', 'CCL8', 'CCL2', 'MT1X'],
    'DC_2':    ['DNASE1L3', 'RGCC', 'CST3', 'CLEC9A', 'SNX3'],
}
subtypes = list(subtype_markers.keys())
for subtype, genes in subtype_markers.items():
    present = [g for g in genes if g in adata_mac.var_names]
    idx = [adata_mac.var_names.get_loc(g) for g in present]
    expr = adata_mac.X[:, idx]
    if sp.issparse(expr):
        expr = expr.toarray()
    adata_mac.obs[f'score_{subtype}'] = expr.mean(axis=1) if idx else 0.0
best = np.argmax(adata_mac.obs[[f'score_{s}' for s in subtypes]].values, axis=1)
adata_mac.obs['mac_subtype'] = pd.Categorical(
    [subtypes[i] for i in best], categories=subtypes)

# ── 拼出统一的 cell_subtype ───────────────────────────────────────────────────
subtype_labels = pd.Series('Unknown', index=adata.obs_names)

cell_type_map = {
    'CD8 T cells':  'CD8 T',   'CD4 T cells': 'CD4 T',  'Tregs':      'Tregs',
    'TFH cells':    'TFH',     'NK cells':    'NK',     'Plasma cells': 'Plasma cells',
    'pDC/Other':    'pDC',
}
for ct, label in cell_type_map.items():
    subtype_labels[adata.obs['cell_type'] == ct] = label

b_cell_types = ['Naive B cells', 'GC B cells', 'Proliferative GC B cells',
                'Memory B cells', 'Age-associated B cells', 'B cells (other)',
                'Plasma cells']
b_idx = adata.obs[adata.obs['cell_type'].isin(b_cell_types)].index
for cell in b_idx.intersection(mal_labels.index):
    subtype_labels[cell] = ('Malignant B'
                            if mal_labels.loc[cell, 'malignancy'] == 'Malignant B cell'
                            else 'Normal B')
plasma_no_label = adata.obs[
    (adata.obs['cell_type'] == 'Plasma cells') & (~adata.obs_names.isin(mal_labels.index))
].index
subtype_labels[plasma_no_label] = 'Plasma cells'

for cell in adata_mac.obs_names:
    subtype_labels[cell] = adata_mac.obs.loc[cell, 'mac_subtype']

adata.obs['cell_subtype'] = pd.Categorical(subtype_labels)
print('Subtypes:', sorted(adata.obs['cell_subtype'].unique().tolist()))

# ── 运行 13 亚型 LIANA ────────────────────────────────────────────────────────
liana.mt.rank_aggregate(
    adata,
    groupby='cell_subtype',
    expr_prop=0.1,
    min_cells=5,
    use_raw=False,
    verbose=False,
    n_perms=100,
    seed=42,
    key_added='liana_res_subtype',
)

liana_res_13 = adata.uns['liana_res_subtype'].copy()
print(f'13-subtype LIANA: {liana_res_13.shape}')
print('Sources:', sorted(liana_res_13['source'].unique().tolist()))

liana_res_13.to_csv(os.path.join(CELLCOMM, 'liana_results_13subtypes_full.csv'),
                    index=False)
liana_sig = liana_res_13[liana_res_13['magnitude_rank'] < 0.05].copy()
liana_sig.to_csv(os.path.join(CELLCOMM, 'liana_results_13subtypes_significant.csv'),
                 index=False)
print(f'[OK] 13subtypes_full={len(liana_res_13):,}  significant={len(liana_sig):,}')
