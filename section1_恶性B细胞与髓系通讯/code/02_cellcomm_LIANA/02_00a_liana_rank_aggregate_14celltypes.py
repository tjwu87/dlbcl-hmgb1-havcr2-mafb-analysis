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
#  ① 全细胞（14 个 cell_type）—— 产出 liana_results_full / significant / 权重矩阵
# =============================================================================
print(f"\nCells: {adata.n_obs}, Cell types: {adata.obs['cell_type'].nunique()}")

liana.mt.rank_aggregate(
    adata,
    groupby='cell_type',
    expr_prop=0.1,       # 基因至少在 10% 细胞中表达
    min_cells=5,         # 每组至少 5 个细胞
    use_raw=False,       # .X 已是 log1p 归一化
    verbose=True,
    n_perms=100,         # 置换次数（用于显著性）
    seed=42,
)

liana_res = adata.uns['liana_res'].copy()
print(f"LIANA results shape: {liana_res.shape}")
print(f"Unique source cell types: {liana_res['source'].nunique()}")

# ── 显著交互：magnitude_rank < 0.05 ───────────────────────────────────────────
sig = liana_res[liana_res['magnitude_rank'] < 0.05].copy()
print(f"Significant interactions (magnitude_rank < 0.05): {len(sig)}")

# ── source × target 通讯权重矩阵：Σ lr_means ──────────────────────────────────
weight_matrix = (sig.groupby(['source', 'target'])['lr_means']
                 .sum()
                 .unstack(fill_value=0))
all_ct = sorted(liana_res['source'].unique())
weight_matrix = weight_matrix.reindex(index=all_ct, columns=all_ct, fill_value=0)
print(f"Weight matrix shape: {weight_matrix.shape}")

weight_matrix.to_csv(os.path.join(CELLCOMM, 'interaction_weight_matrix.csv'))
liana_res.to_csv(os.path.join(CELLCOMM, 'liana_results_full.csv'), index=False)
sig.to_csv(os.path.join(CELLCOMM, 'liana_results_significant.csv'), index=False)
print('\n[OK] 已写出 interaction_weight_matrix.csv / liana_results_full.csv / '
      'liana_results_significant.csv')
