# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  构建 pySCENIC 输入 loom（论文 section3 上游 · 第 1 步）
#
#  来源：Biomni 平台分析记录（19_PAGA与注释重做）
#
#  规则（照原文）：
#    · 细胞 = Monocytes/Macrophages 经 mac_subtype 注释后的 327 个
#      （Mono 90 / IFN_TAM 161 / LA_TAM 76）
#    · 表达取 adata 原始 .X（log-normalized，max ≈ 7.5）
#    · 基因过滤：至少在 5% 细胞中表达 → (327, 9118)
#    · loom 要求 genes × cells
#
#  环境：scanpy + anndata + loompy + numpy
#  产出：<SCENIC_DB>/mono_ifntam_latam.loom（≈3.3 MB）
#        以及 <数据根>/GSE182434/trajectory/ 下的持久副本
# =============================================================================

import os
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import loompy

GSE = translate(r'D:\bulk-download\GSE182434')          # -> <数据根>/GSE182434
SCENIC_DB = os.environ.get('SCENIC_DB', os.path.join(GSE, 'scenic_db'))
TRAJ = os.path.join(GSE, 'trajectory')
os.makedirs(SCENIC_DB, exist_ok=True)
os.makedirs(TRAJ, exist_ok=True)

# ── 加载已注释的 adata（mac_subtype 由 04_00b 之前的注释步骤写入）──────────────
adata = sc.read_h5ad(os.path.join(GSE, 'adata_mac_annotated_patched.h5ad'))
adata_3 = adata[adata.obs['mac_subtype'].isin(['Mono', 'IFN_TAM', 'LA_TAM'])].copy()
print('Subtypes:', adata_3.obs['mac_subtype'].value_counts().to_dict())

X = adata_3.X.toarray() if sp.issparse(adata_3.X) else np.asarray(adata_3.X)

# ── 基因过滤：≥5% 细胞表达 ────────────────────────────────────────────────────
mask = (X > 0).mean(axis=0) >= 0.05
X_filtered = X[:, mask]
gene_names = np.asarray(adata_3.var_names)[mask]
cell_names = np.asarray(adata_3.obs_names)
print(f'Genes expressed in >=5% cells: {int(mask.sum())} / {len(mask)}')
print(f'Final matrix: {X_filtered.shape}  (cells x genes)')

# ── 写 loom（genes x cells）────────────────────────────────────────────────────
loom_path = os.path.join(SCENIC_DB, 'mono_ifntam_latam.loom')
if os.path.exists(loom_path):
    os.remove(loom_path)

loompy.create(
    loom_path,
    X_filtered.T,
    {'Gene': gene_names},
    {'CellID': cell_names,
     'mac_subtype': adata_3.obs['mac_subtype'].values.astype(str)},
)
print(f'Loom created: {loom_path}  ({os.path.getsize(loom_path) / 1e6:.1f} MB)')

with loompy.connect(loom_path, 'r') as ds:
    print(f'Verify — genes={ds.shape[0]}, cells={ds.shape[1]}')

shutil.copy2(loom_path, os.path.join(TRAJ, 'mono_ifntam_latam.loom'))
print(f'[OK] 持久副本 -> {TRAJ}')
