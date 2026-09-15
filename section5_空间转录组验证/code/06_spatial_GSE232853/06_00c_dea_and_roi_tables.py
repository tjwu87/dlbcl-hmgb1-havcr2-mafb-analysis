# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  GSE232853 —— 差异表达 + 配对 ROI 表 + 分组表
#                                                     （section5 上游 · 第 3 步）
#
#  来源：Biomni 平台分析记录 `14_空间转录组分析.md`
#        · L15079-15129  逐基因 Mann-Whitney（CD68+ 内 DLBCL vs Normal）
#        · L17304-17353  scanpy Wilcoxon rank_genes_groups（同上，scanpy 版）
#        · L20460-20577  fig2 Wilcoxon 统计 + 配对 ROI 表 + 相关性表 + 分组表
#
#  产出（写入 <数据根>/GSE232853_v2/，均为下游 06_01 直接读取）：
#        fig2_wilcoxon_HMGB1_HAVCR2_MAFB.csv     CD20 vs CD68 三基因 Wilcoxon
#        fig3_DEA_CD68_scanpy.csv                scanpy Wilcoxon DEA（火山图数据）
#        task4_DEA_CD68_DLBCL_vs_Normal.csv      逐基因 MWU + BH 校正
#        fig4_paired_ROI_16gene_LATAM.csv        配对 ROI 表（含 16 基因 LA_TAM_Score）
#        fig4_correlation_stats.csv              Fig4 A/B 的 Pearson/Spearman
#        fig5_DLBCL_LATAM_group_assignment.csv   DLBCL 按 LA_TAM_Score 中位数分 Rich/Poor
#
#  关键定义（照原文）：
#        LA_TAM_GENES（16 个，FTL 不在 panel 故剔除）
#        LA_TAM_Score = CD68+ ROI 这 16 基因的均值
#        Rich/Poor  = 在 DLBCL 样本内按 LA_TAM_Score 中位数二分（各 57 个）
#
#  环境：pandas / numpy / scanpy / anndata / scipy / statsmodels
# =============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

OUT = translate(r'D:\bulk-download\GSE232853_v2')     # -> <数据根>/GSE232853_v2
adata = ad.read_h5ad(os.path.join(OUT, 'GSE232853_adata_v2.h5ad'))
expr = pd.DataFrame(adata.layers['q3_lognorm'],
                    index=adata.obs_names, columns=adata.var_names)
meta = adata.obs.copy()

LATAM_GENES = ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD', 'GPNMB', 'APOC1',
               'PLA2G2D', 'CAPG', 'MMP9', 'NUPR1', 'CTSL', 'RARRES1', 'IL32',
               'FUCA1', 'LGMN']

# ══════════════════════════════════════════════════════════════════════════════
# ① fig2 —— HMGB1 / HAVCR2 / MAFB 在 CD20+ vs CD68+ 的 Wilcoxon 统计
# ══════════════════════════════════════════════════════════════════════════════
cd20_idx = meta['Mask'] == 'CD20'
cd68_idx = meta['Mask'] == 'CD68'

rows = []
for g in ['HMGB1', 'HAVCR2', 'MAFB']:
    a = expr.loc[cd20_idx, g].values
    b = expr.loc[cd68_idx, g].values
    stat, pval = mannwhitneyu(a, b, alternative='two-sided')
    rows.append({'Gene': g, 'n_CD20': int(cd20_idx.sum()), 'n_CD68': int(cd68_idx.sum()),
                 'median_CD20': np.median(a), 'median_CD68': np.median(b),
                 'mean_CD20': np.mean(a), 'mean_CD68': np.mean(b),
                 'MWU_stat': stat, 'p_value': pval})
df_w = pd.DataFrame(rows)
df_w['padj_BH'] = multipletests(df_w['p_value'], method='fdr_bh')[1]
df_w.to_csv(os.path.join(OUT, 'fig2_wilcoxon_HMGB1_HAVCR2_MAFB.csv'), index=False)
print('  [OK] fig2_wilcoxon_HMGB1_HAVCR2_MAFB.csv')

# ══════════════════════════════════════════════════════════════════════════════
# ② fig3 —— scanpy Wilcoxon DEA（CD68+ 内 DLBCL vs Normal），供火山图
# ══════════════════════════════════════════════════════════════════════════════
adata_cd68 = adata[adata.obs['Mask'] == 'CD68'].copy()
adata_cd68.X = adata_cd68.layers['q3_lognorm'].copy()

sc.tl.rank_genes_groups(adata_cd68, groupby='Tissue_Type', groups=['DLBCL'],
                        reference='Normal', method='wilcoxon',
                        key_added='wilcoxon_DLBCL_vs_Normal', pts=True)
result = sc.get.rank_genes_groups_df(adata_cd68, group='DLBCL',
                                     key='wilcoxon_DLBCL_vs_Normal')
result.columns = ['gene', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj',
                  'pts', 'pts_rest']
result = result.dropna(subset=['pvals_adj'])
result['-log10_padj'] = -np.log10(result['pvals_adj'].clip(lower=1e-300))
result.to_csv(os.path.join(OUT, 'fig3_DEA_CD68_scanpy.csv'), index=False)
print(f'  [OK] fig3_DEA_CD68_scanpy.csv  ({len(result)} genes, '
      f"sig={int(((result['pvals_adj'] <= 0.05) & (result['logfoldchanges'].abs() >= 0.5)).sum())})")

# ══════════════════════════════════════════════════════════════════════════════
# ③ task4 —— 逐基因 Mann-Whitney + BH（另一套 DEA 口径，论文补图用）
# ══════════════════════════════════════════════════════════════════════════════
expr_cd68 = expr.loc[cd68_idx]
meta_cd68 = meta.loc[cd68_idx]
dl = meta_cd68['Tissue_Type'] == 'DLBCL'
nr = meta_cd68['Tissue_Type'] == 'Normal'
print(f'  CD68+ ROIs — DLBCL: {int(dl.sum())}, Normal: {int(nr.sum())}')

res = []
for gene in expr_cd68.columns:
    a = expr_cd68.loc[dl, gene].values
    b = expr_cd68.loc[nr, gene].values
    if np.std(a) == 0 and np.std(b) == 0:
        continue
    stat, pval = mannwhitneyu(a, b, alternative='two-sided')
    res.append({'Gene': gene, 'LogFC_DLBCL_vs_Normal': np.mean(a) - np.mean(b),
                'Mean_DLBCL': np.mean(a), 'Mean_Normal': np.mean(b),
                'MWU_stat': stat, 'p_value': pval})
dea = pd.DataFrame(res)
dea['padj_BH'] = multipletests(dea['p_value'], method='fdr_bh')[1]
dea = dea.sort_values('padj_BH')
dea.to_csv(os.path.join(OUT, 'task4_DEA_CD68_DLBCL_vs_Normal.csv'), index=False)
print(f"  [OK] task4_DEA_CD68_DLBCL_vs_Normal.csv  ({len(dea)} genes, "
      f"sig={int(((dea['padj_BH'] < 0.05) & (dea['LogFC_DLBCL_vs_Normal'].abs() > 0.5)).sum())})")

# ══════════════════════════════════════════════════════════════════════════════
# ④ 配对 ROI 表：CD20+ 与 CD68+ 按 Unique_ROI 内连接
# ══════════════════════════════════════════════════════════════════════════════
cd20_meta = meta[cd20_idx].copy()
cd68_meta = meta[cd68_idx].copy()
cd20_expr = expr.loc[cd20_meta.index].copy()
cd68_expr = expr.loc[cd68_meta.index].copy()
cd20_expr.index = cd20_meta['Unique_ROI'].values
cd68_expr.index = cd68_meta['Unique_ROI'].values
cd20_meta = cd20_meta.set_index('Unique_ROI')
cd68_meta = cd68_meta.set_index('Unique_ROI')

common_rois = cd20_expr.index.intersection(cd68_expr.index)
cd20_paired = cd20_expr.loc[common_rois]
cd68_paired = cd68_expr.loc[common_rois]
tissue_type = cd20_meta.loc[common_rois, 'Tissue_Type']
slide_id = cd20_meta.loc[common_rois, 'Slide_ID']

la_tam_score = cd68_paired[LATAM_GENES].mean(axis=1)

paired_df = pd.DataFrame({
    'Unique_ROI': common_rois,
    'Slide_ID': slide_id.values,
    'Tissue_Type': tissue_type.values,
    'LA_TAM_Score': la_tam_score.values,
    'B_HMGB1': cd20_paired['HMGB1'].values,
    'B_HAVCR2': cd20_paired['HAVCR2'].values,
    'Mac_HAVCR2': cd68_paired['HAVCR2'].values,
    'Mac_MAFB': cd68_paired['MAFB'].values,
    'Mac_APOE': cd68_paired['APOE'].values,
})
paired_df.to_csv(os.path.join(OUT, 'fig4_paired_ROI_16gene_LATAM.csv'), index=False)
print(f'  [OK] fig4_paired_ROI_16gene_LATAM.csv  ({len(paired_df)} paired ROIs)')

# ══════════════════════════════════════════════════════════════════════════════
# ⑤ 相关性统计表（Fig4 A: B_HMGB1 vs Mac_HAVCR2；B: B_HMGB1 vs LA_TAM_Score）
# ══════════════════════════════════════════════════════════════════════════════
b_hmgb1 = cd20_paired['HMGB1'].values
m_havcr2 = cd68_paired['HAVCR2'].values
la_tam = la_tam_score.values
tissue = tissue_type.values

corr_rows = []
for subset, mask in [('All', np.ones(len(tissue), dtype=bool)),
                     ('DLBCL', tissue == 'DLBCL'),
                     ('Normal', tissue == 'Normal')]:
    for fig, x, y, ylabel in [('Fig_A', b_hmgb1[mask], m_havcr2[mask], 'Mac_HAVCR2'),
                              ('Fig_B', b_hmgb1[mask], la_tam[mask], 'LA_TAM_Score_16gene')]:
        r, pr = pearsonr(x, y)
        rho, ps = spearmanr(x, y)
        corr_rows.append({'Figure': fig, 'Subset': subset, 'n': int(mask.sum()),
                          'X': 'B_HMGB1', 'Y': ylabel,
                          'Pearson_r': round(r, 4), 'Pearson_p': pr,
                          'Spearman_rho': round(rho, 4), 'Spearman_p': ps})
pd.DataFrame(corr_rows).to_csv(os.path.join(OUT, 'fig4_correlation_stats.csv'),
                               index=False)
print('  [OK] fig4_correlation_stats.csv')

# ══════════════════════════════════════════════════════════════════════════════
# ⑥ DLBCL 按 LA_TAM_Score 中位数分 LA_TAM_Rich / Poor
# ══════════════════════════════════════════════════════════════════════════════
dlbcl_mask = tissue == 'DLBCL'
la_tam_dlbcl = la_tam_score.loc[tissue_type == 'DLBCL']
median_score = la_tam_dlbcl.median()
print(f'  LA_TAM_Score median (DLBCL): {median_score:.4f}')

group_df = pd.DataFrame({
    'Unique_ROI': common_rois[dlbcl_mask],
    'Slide_ID': slide_id.values[dlbcl_mask],
    'LA_TAM_Score': la_tam[dlbcl_mask],
    'LA_TAM_Group': np.where(la_tam[dlbcl_mask] >= median_score,
                             'LA_TAM_Rich', 'LA_TAM_Poor'),
})
group_df.to_csv(os.path.join(OUT, 'fig5_DLBCL_LATAM_group_assignment.csv'), index=False)
print('  [OK] fig5_DLBCL_LATAM_group_assignment.csv  '
      f"{group_df['LA_TAM_Group'].value_counts().to_dict()}")
print('[OK] DEA 与 ROI 相关表全部写出')
