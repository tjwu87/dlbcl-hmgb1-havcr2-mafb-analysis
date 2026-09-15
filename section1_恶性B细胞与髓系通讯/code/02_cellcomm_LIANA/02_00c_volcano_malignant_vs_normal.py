# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  恶性 vs 正常 B 细胞 → 巨噬细胞的差异通讯（LIANA + 火山图数据）
#                                       （section1 上游 · 图 S1/正文火山图）
#
#  来源：Biomni 平台分析记录 `01_对话记录/02_细胞通讯LIANA全细胞.md`
#        · L2420-2456  用 CNV 判定结果改写 cell_type_mal，重跑 LIANA
#        · L6041-6096  取 Malignant/Normal B → Mono/Mac 两组，按 LR 对合并
#        · L6158-6214  Delta_Magnitude / p 值 / 四类标签
#        · L6383-6387  火山图数据落盘
#
#  依赖上一步：`cnv/malignancy_classification.csv`（由 01_00 产出）
#
#  产出（写入 <数据根>/GSE182434/cellcomm/）：
#        volcano_data_malignant_vs_normal_monomac.csv
#            LR_pair / ligand_complex / receptor_complex / strength_mal /
#            strength_nor / Delta_Magnitude / pval / neg_log10_pval / category
#        volcano_malignant_vs_normal_monomac.png/svg   火山图
#
#  定义（照原文）：
#        strength      = 1 − magnitude_rank（越大越强）
#        Delta         = strength_mal − strength_nor（>0 表示恶性侧更强）
#        pval          = 取 Malignant 组的 cellphone_pvals（clip 到 1e-10）
#        category      = 显著(p<0.05) 且 Δ>+0.05 → Malignant-enriched
#                        显著 且 Δ<−0.05        → Normal-enriched
#                        显著但 |Δ|≤0.05         → Significant (small Δ)
#                        否则                    → Not significant
#        只在一组出现的 LR 对：另一组 magnitude_rank 补 1.0（= 最弱）
#
#  环境：scanpy / anndata / liana / pandas / numpy / matplotlib / seaborn / adjustText
# =============================================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import liana
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

GSE = translate(r'D:\bulk-download\GSE182434')          # -> <数据根>/GSE182434
CELLCOMM = os.path.join(GSE, 'cellcomm')
os.makedirs(CELLCOMM, exist_ok=True)

# ── ① 用 CNV 判定改写细胞类型标签 ─────────────────────────────────────────────
mal_df = pd.read_csv(os.path.join(GSE, 'cnv', 'malignancy_classification.csv'),
                     index_col=0)
print('Malignancy labels:', mal_df['malignancy'].value_counts().to_dict())

adata_mal = sc.read_h5ad(os.path.join(GSE, 'adata_processed.h5ad'))
adata_mal.obs['cell_type_mal'] = adata_mal.obs['cell_type'].astype(str)

for idx in adata_mal.obs_names:
    if idx in mal_df.index:
        status = str(mal_df.loc[idx, 'malignancy'])
        if status == 'Malignant B cell':
            adata_mal.obs.loc[idx, 'cell_type_mal'] = 'Malignant B cell'
        elif 'Normal B cell' in status:
            adata_mal.obs.loc[idx, 'cell_type_mal'] = 'Normal B cell'

print(adata_mal.obs['cell_type_mal'].value_counts().to_string())

# ── ② 重跑 LIANA（13 亚型那套参数）────────────────────────────────────────────
liana.mt.rank_aggregate(
    adata_mal,
    groupby='cell_type_mal',
    expr_prop=0.1,
    min_cells=5,
    use_raw=False,
    verbose=False,
    n_perms=100,
    seed=42,
)
liana_res_mal = adata_mal.uns['liana_res'].copy()
print(f'LIANA done: {len(liana_res_mal):,} interactions')
print('Sources:', sorted(liana_res_mal['source'].unique().tolist()))

# ── ③ 取两组子表并按 LR 对合并 ────────────────────────────────────────────────
mal_sub = liana_res_mal[(liana_res_mal['source'] == 'Malignant B cell') &
                        (liana_res_mal['target'] == 'Monocytes/Macrophages')].copy()
nor_sub = liana_res_mal[(liana_res_mal['source'] == 'Normal B cell') &
                        (liana_res_mal['target'] == 'Monocytes/Macrophages')].copy()
print(f'Malignant B -> Mono/Mac: {len(mal_sub)} LR pairs')
print(f'Normal    B -> Mono/Mac: {len(nor_sub)} LR pairs')

cols = ['magnitude_rank', 'lrscore', 'scaled_weight', 'cellphone_pvals']
mal_sub['lr_key'] = mal_sub['ligand_complex'] + '__' + mal_sub['receptor_complex']
nor_sub['lr_key'] = nor_sub['ligand_complex'] + '__' + nor_sub['receptor_complex']

merged = pd.merge(
    mal_sub[['lr_key', 'ligand_complex', 'receptor_complex'] + cols],
    nor_sub[['lr_key'] + cols],
    on='lr_key', suffixes=('_mal', '_nor'), how='outer',
)
# 只在一组出现的 LR 对：另一组补"最弱"值
merged['magnitude_rank_mal'] = merged['magnitude_rank_mal'].fillna(1.0)
merged['magnitude_rank_nor'] = merged['magnitude_rank_nor'].fillna(1.0)
merged['lrscore_mal'] = merged['lrscore_mal'].fillna(0.0)
merged['lrscore_nor'] = merged['lrscore_nor'].fillna(0.0)
print(f'Merged table: {len(merged)} LR pairs')

# ── ④ 强度差 / p 值 / 分类 ────────────────────────────────────────────────────
merged['strength_mal'] = 1 - merged['magnitude_rank_mal']
merged['strength_nor'] = 1 - merged['magnitude_rank_nor']
merged['Delta_Magnitude'] = merged['strength_mal'] - merged['strength_nor']

# p 值取 Malignant 组（该组覆盖全部 LR 对）
merged['pval'] = merged['cellphone_pvals_mal'].fillna(1.0).clip(lower=1e-10)
merged['neg_log10_pval'] = -np.log10(merged['pval'])
merged['LR_pair'] = merged['ligand_complex'] + ' \u2192 ' + merged['receptor_complex']

PVAL_THRESH, DELTA_THRESH = 0.05, 0.05


def classify(row):
    sig = row['pval'] < PVAL_THRESH
    up = row['Delta_Magnitude'] > DELTA_THRESH
    dn = row['Delta_Magnitude'] < -DELTA_THRESH
    if sig and up:
        return 'Malignant-enriched'
    if sig and dn:
        return 'Normal-enriched'
    if sig:
        return 'Significant (small \u0394)'
    return 'Not significant'


merged['category'] = merged.apply(classify, axis=1)
print(merged['category'].value_counts().to_string())

# ── ⑤ 落盘（火山图数据）──────────────────────────────────────────────────────
out_cols = ['LR_pair', 'ligand_complex', 'receptor_complex', 'strength_mal',
            'strength_nor', 'Delta_Magnitude', 'pval', 'neg_log10_pval', 'category']
table = merged[out_cols].sort_values('Delta_Magnitude', ascending=False)
table.to_csv(os.path.join(CELLCOMM,
                          'volcano_data_malignant_vs_normal_monomac.csv'),
             index=False)
print(f'[OK] volcano_data_malignant_vs_normal_monomac.csv  ({len(table)} rows)')

# ── ⑥ 火山图 ─────────────────────────────────────────────────────────────────
palette = {'Malignant-enriched': '#d62728', 'Normal-enriched': '#2ca02c',
           'Significant (small \u0394)': '#ff7f0e', 'Not significant': '#bbbbbb'}
size_map = {'Malignant-enriched': 60, 'Normal-enriched': 60,
            'Significant (small \u0394)': 40, 'Not significant': 20}
alpha_map = {'Malignant-enriched': 0.85, 'Normal-enriched': 0.85,
             'Significant (small \u0394)': 0.75, 'Not significant': 0.45}

sns.set_theme(style='ticks', font_scale=1.0)
fig, ax = plt.subplots(figsize=(11, 8))
for cat in palette:
    sub = table[table['category'] == cat]
    if sub.empty:
        continue
    ax.scatter(sub['Delta_Magnitude'], sub['neg_log10_pval'], s=size_map[cat],
               c=palette[cat], alpha=alpha_map[cat],
               edgecolors='white' if cat == 'Not significant' else '#333333',
               linewidths=0.3 if cat == 'Not significant' else 0.5,
               label=f'{cat} (n={len(sub)})', zorder=3)

pval_line = -np.log10(0.05)
for y in (pval_line,):
    ax.axhline(y, color='#555555', linewidth=1.0, linestyle='--', alpha=0.7, zorder=2)
for x in (0.05, -0.05):
    ax.axvline(x, color='#555555', linewidth=1.0, linestyle='--', alpha=0.7, zorder=2)
ax.axvline(0, color='#999999', linewidth=0.8, alpha=0.5, zorder=1)

top_mal = table[table['category'] == 'Malignant-enriched'].nlargest(10, 'neg_log10_pval')
top_nor = table[table['category'] == 'Normal-enriched'].nlargest(5, 'neg_log10_pval')
texts = []
for _, row in pd.concat([top_mal, top_nor]).iterrows():
    color = '#d62728' if row['category'] == 'Malignant-enriched' else '#2ca02c'
    texts.append(ax.text(row['Delta_Magnitude'], row['neg_log10_pval'], row['LR_pair'],
                         fontsize=7.5, color=color, fontweight='bold', zorder=5,
                         path_effects=[pe.withStroke(linewidth=2.0, foreground='white')]))
try:
    from adjustText import adjust_text
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6))
except ImportError:
    pass

ax.set_xlabel('\u0394Magnitude  (Malignant B \u2212 Normal B interaction strength)',
              fontsize=10.5, labelpad=8)
ax.set_ylabel('\u2212log\u2081\u2080(p-value)  [CellPhoneDB permutation test]',
              fontsize=10.5, labelpad=8)
ax.set_title('Differential Cell-Cell Communication: Malignant B cell vs Normal B cell\n'
             'Target: Monocytes/Macrophages  |  LIANA rank_aggregate',
             fontsize=12, fontweight='bold', pad=14)
ax.legend(loc='upper left', fontsize=8.5, frameon=True, framealpha=0.9,
          title='Category', title_fontsize=9)
sns.despine(ax=ax)
plt.tight_layout()
for fmt in ('png', 'svg'):
    fig.savefig(os.path.join(CELLCOMM, f'volcano_malignant_vs_normal_monomac.{fmt}'),
                dpi=150, bbox_inches='tight', format=fmt)
plt.close(fig)
print('[OK] volcano_malignant_vs_normal_monomac.png/.svg')
