# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import pandas as pd

# 侦察1：16基因评分文件结构
latam = pd.read_csv(translate(r'D:\bulk-download\GSE232853_v2\fig4_paired_ROI_16gene_LATAM.csv'))
print(latam.columns.tolist())
print(latam.head(3))
print('行数:', len(latam))

# 侦察2：meta注释结构（确认compartment/tissue列名）
meta = pd.read_csv(translate(r'D:\bulk-download\GSE232853_v2\meta_df.csv'), index_col=0)
print(meta.columns.tolist())

# 侦察3：Q3矩阵结构
expr = pd.read_csv(translate(r'D:\bulk-download\GSE232853_v2\expression_q3_lognorm.csv'), index_col=0, nrows=5)
print('Q3矩阵形状(前5行×全部列):', expr.shape)

import pandas as pd
import numpy as np
from scipy.stats import spearmanr

# ---------- 1) 读入 + 转置 ----------
latam = pd.read_csv(translate(r'D:\bulk-download\GSE232853_v2\fig4_paired_ROI_16gene_LATAM.csv'))
expr  = pd.read_csv(translate(r'D:\bulk-download\GSE232853_v2\expression_q3_lognorm.csv'), index_col=0)

if 'ZNF83' in expr.index:
    expr = expr.T

# ---------- 2) 解析AOI索引（Index→Series再split）----------
idx = expr.index.to_series().astype(str)          # 先转Series
parts = idx.str.split(r'\s*\|\s*', expand=True)   # 此时返回DataFrame

expr_chk = expr.copy()
expr_chk['Slide_ID'] = parts[0].str.strip()
expr_chk['ROI_Num']  = parts[1].str.strip().str.lstrip('0')
expr_chk['ROI_Num']  = expr_chk['ROI_Num'].replace('', np.nan)
expr_chk['Mask']     = parts[2].str.strip()

print('Mask分布:'); print(expr_chk['Mask'].value_counts())


# ---------- 3) 取CD68侧 + 构造对齐键 ----------
cd68 = expr_chk[expr_chk['Mask'] == 'CD68'].copy()
cd68['key'] = cd68['Slide_ID'] + '_' + cd68['ROI_Num'].astype(str)

latam = latam.copy()
latam['key'] = latam['Unique_ROI'].str.replace(
    r'_(\d+)$', lambda m: '_' + str(int(m.group(1))), regex=True)

common = latam['key'].isin(cd68['key'])
print(f'\nlatam可匹配CD68 AOI: {common.sum()}/{len(latam)}')

if common.sum() == len(latam):
    matched = latam.set_index('key')
    cd68m = cd68.set_index('key').loc[matched.index]
else:
    # 打印不匹配的键，看差异在哪
    miss = latam.loc[~common, ['Unique_ROI', 'key']]
    print('未匹配示例:')
    print(miss.head(10).to_string())
    print('CD68键示例:', cd68['key'].head(5).tolist())
    matched = latam.loc[common].set_index('key')
    cd68m = cd68.set_index('key').loc[matched.index]

# ---------- 4) 基因面板 + 自校验 ----------
genes16 = ['PTGDS','CCL18','APOE','CHI3L1','CTSD','GPNMB','APOC1','PLA2G2D',
           'CAPG','MMP9','NUPR1','CTSL','RARRES1','IL32','FUCA1','LGMN']
genes5  = ['APOE','APOC1','ACP5','CCL18','CTSD']

missing = [g for g in genes16 + genes5 if g not in expr_chk.columns]
print('缺失基因:', missing or '无')

recomp16 = cd68m[genes16].mean(axis=1)
r_chk, _ = spearmanr(recomp16, matched['LA_TAM_Score'])
print(f'\n[自校验] 重算16基因 vs LA_TAM_Score: rho={r_chk:.4f}, '
      f'最大差={np.abs(recomp16 - matched["LA_TAM_Score"]).max():.4f}')

# ---------- 5) 主结果 ----------
score5 = cd68m[genes5].mean(axis=1)
r, p = spearmanr(score5, matched['LA_TAM_Score'])
print(f'\n16 vs 5基因 LA-TAM评分, 配对CD68+ ROIs (n={len(matched)}): rho={r:.3f}, p={p:.2g}')

for tis in ['DLBCL', 'Normal']:
    sel = (matched['Tissue_Type'] == tis).values
    if sel.sum() > 3:
        r_t, p_t = spearmanr(score5[sel], matched.loc[sel, 'LA_TAM_Score'])
        print(f'  {tis} (n={sel.sum()}): rho={r_t:.3f}, p={p_t:.2g}')
