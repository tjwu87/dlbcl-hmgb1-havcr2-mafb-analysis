# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  GSE232853 —— GSEA prerank + ssGSEA（Hallmark）
#                                                     （section5 上游 · 第 4 步）
#
#  来源：Biomni 平台分析记录 `14_空间转录组分析.md`
#        · L19603-19664  重建配对数据 + 16 基因 LA_TAM_Score + 中位数分组
#        · L19746-19818  GSEA prerank（CD20+ 内 Rich vs Poor）
#        · L20685-20879  ssGSEA（每个 CD20+ ROI 一套 Hallmark 打分）
#
#  两种富集口径（论文 Fig5-f 与 Fig6 各用其一）：
#    · GSEA prerank：先把 Rich/Poor 的均值差作为基因排序，再跑 prerank
#    · ssGSEA      ：对每个样本单独打通路分，再比较两组分布
#
#  产出（写入 <数据根>/GSE232853_v2/）：
#        fig5_GSEA_Hallmark_LATAM_Rich_vs_Poor_CD20.csv   prerank 结果（50 条通路）
#        fig5_GSEA_ranking_CD20_DLBCL.csv                 排序表（Rich - Poor 均值差）
#        fig6_ssGSEA_Hallmark_CD20_DLBCL.csv              样本 × 通路 NES 矩阵
#
#  环境：pandas / numpy / anndata / gseapy（首次跑会自动拉 MSigDB Hallmark）
# =============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

OUT = translate(r'D:\bulk-download\GSE232853_v2')     # -> <数据根>/GSE232853_v2
adata = ad.read_h5ad(os.path.join(OUT, 'GSE232853_adata_v2.h5ad'))
expr = pd.DataFrame(adata.layers['q3_lognorm'],
                    index=adata.obs_names, columns=adata.var_names)
meta = adata.obs.copy()

LATAM_GENES = ['PTGDS', 'CCL18', 'APOE', 'CHI3L1', 'CTSD', 'GPNMB', 'APOC1',
               'PLA2G2D', 'CAPG', 'MMP9', 'NUPR1', 'CTSL', 'RARRES1', 'IL32',
               'FUCA1', 'LGMN']

# ── ① 重建配对数据 + LA_TAM_Score + Rich/Poor 分组 ───────────────────────────
cd20_meta = meta[meta['Mask'] == 'CD20'].copy()
cd68_meta = meta[meta['Mask'] == 'CD68'].copy()
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

la_tam_score = cd68_paired[LATAM_GENES].mean(axis=1)
dlbcl_mask = tissue_type == 'DLBCL'
cd20_dlbcl = cd20_paired.loc[dlbcl_mask]
la_tam_dlbcl = la_tam_score.loc[dlbcl_mask]
median_score = la_tam_dlbcl.median()
rich_rois = la_tam_dlbcl[la_tam_dlbcl >= median_score].index
poor_rois = la_tam_dlbcl[la_tam_dlbcl < median_score].index
print(f"DLBCL paired ROIs: {int(dlbcl_mask.sum())}  "
      f"(Rich {len(rich_rois)} / Poor {len(poor_rois)})")

# ── ② 排序表：CD20+ 内 Rich 均值 − Poor 均值（log 空间，≈ log2FC）────────────
ranking = (cd20_dlbcl.loc[rich_rois].mean(axis=0)
           - cd20_dlbcl.loc[poor_rois].mean(axis=0)).sort_values(ascending=False)
rnk_df = ranking.reset_index()
rnk_df.columns = ['gene', 'MeanDiff_Rich_minus_Poor']
rnk_df.to_csv(os.path.join(OUT, 'fig5_GSEA_ranking_CD20_DLBCL.csv'), index=False)
print(f'  Ranking: {len(rnk_df)} genes, range {ranking.min():.3f} – {ranking.max():.3f}')

# ── ③ GSEA prerank（Hallmark）─────────────────────────────────────────────────
print('Running GSEA prerank (Hallmark)...')
res = gp.prerank(rnk=rnk_df, gene_sets='MSigDB_Hallmark_2020',
                 min_size=10, max_size=500, permutation_num=1000,
                 seed=42, verbose=False, outdir=None)
gsea_df = res.res2d.copy().sort_values('NES', ascending=False)
gsea_df.to_csv(os.path.join(OUT, 'fig5_GSEA_Hallmark_LATAM_Rich_vs_Poor_CD20.csv'),
               index=False)
print(f"  [OK] prerank: {len(gsea_df)} pathways, "
      f"FDR<0.25 = {int((gsea_df['FDR q-val'] < 0.25).sum())}")
print(gsea_df[gsea_df['NES'] > 0][['Term', 'NES', 'FDR q-val']].head(6).to_string(index=False))

# ── ④ ssGSEA（Hallmark）───────────────────────────────────────────────────────
# 用已经分好组的表，只取有分组的 CD20+ DLBCL ROI
group_df = pd.read_csv(os.path.join(OUT, 'fig5_DLBCL_LATAM_group_assignment.csv'))
roi_to_group = group_df.set_index('Unique_ROI')['LA_TAM_Group'].to_dict()

cd20_only = expr.loc[cd20_meta.index].copy()
cd20_only.index = cd20_meta['Unique_ROI'].values
cd20_dlbcl_expr = cd20_only.loc[cd20_only.index.isin(roi_to_group.keys())]
groups = pd.Series([roi_to_group[r] for r in cd20_dlbcl_expr.index],
                   index=cd20_dlbcl_expr.index)
print(f"\nssGSEA input: {cd20_dlbcl_expr.shape[0]} ROIs "
      f"({groups.value_counts().to_dict()})")

# ssGSEA 要求 基因 × 样本
ss = gp.ssgsea(data=cd20_dlbcl_expr.T, gene_sets='MSigDB_Hallmark_2020',
               sample_norm_method='rank', min_size=10, max_size=500,
               scale=True, outdir=None, verbose=False)

# 结果展平为 样本 × 通路（取 NES）
scores_wide = ss.res2d.pivot(index='Name', columns='Term', values='NES')
scores_wide['LA_TAM_Group'] = scores_wide.index.map(roi_to_group)
scores_wide.to_csv(os.path.join(OUT, 'fig6_ssGSEA_Hallmark_CD20_DLBCL.csv'))
print(f'  [OK] ssGSEA: {scores_wide.shape[0]} ROIs x '
      f'{scores_wide.shape[1] - 1} pathways')
print('[OK] GSEA / ssGSEA 全部写出')
