# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  GSE232853 —— Target_Gene_Score + 空间共变主表
#                                                     （section5 上游 · 第 5 步）
#
#  来源：Biomni 平台分析记录 `14_空间转录组分析.md`
#        · L21568-21608  在 CD68+ ROI 上用 scanpy.score_genes 打 Target_Gene_Score
#        · L21645-21678  与配对 ROI 表合并成 fig7 主表
#
#  Target_Gene_Score 用的是 **MAFB 扰动下游的 23 个靶基因**（与 section4 的
#  CellOracle / section6 的预后模型同一套基因），把「MAFB 活性」投射到空间上。
#
#  产出（写入 <数据根>/GSE232853_v2/）：
#        fig7_master_spatial_coevolution.csv
#            ROI_ID / Tissue_Type / CD20_Ligand(=B_HMGB1) / CD68_Receptor(=Mac_HAVCR2)
#            / CD68_TF(=Mac_MAFB) / CD68_Target_Gene_Score
#
#  环境：pandas / numpy / scanpy / anndata
# =============================================================================

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

OUT = translate(r'D:\bulk-download\GSE232853_v2')     # -> <数据根>/GSE232853_v2

# MAFB 扰动下游 23 靶基因（与 section4 / section6 同一套）
TARGET_23 = ['CEP192', 'CHCHD10', 'NAGK', 'PPRC1', 'BASP1', 'MTR', 'TIMM17A',
             'AAK1', 'MGAT1', 'USP12', 'FAM20A', 'APOL3', 'GSDMD', 'ELMO1',
             'PICALM', 'RASSF4', 'IFIT3', 'NAE1', 'NAP1L4', 'PPP2R5A',
             'MRGBP', 'PPP1R18', 'CTSL']

adata = ad.read_h5ad(os.path.join(OUT, 'GSE232853_adata_v2.h5ad'))

# ── ① CD68+ 子集上打分 ────────────────────────────────────────────────────────
adata_cd68 = adata[adata.obs['Mask'] == 'CD68'].copy()
adata_cd68.X = adata_cd68.layers['q3_lognorm'].copy()
print(f'CD68+ subset: {adata_cd68.shape}')

missing = [g for g in TARGET_23 if g not in adata_cd68.var_names]
if missing:
    print(f'  WARNING: {len(missing)} target genes absent from panel: {missing}')

sc.tl.score_genes(adata_cd68, gene_list=TARGET_23,
                  score_name='Target_Gene_Score')
score = adata_cd68.obs['Target_Gene_Score']
print(f'Target_Gene_Score range: {score.min():.4f} – {score.max():.4f} '
      f'(mean {score.mean():.4f})')

score_series = pd.Series(score.values,
                         index=adata_cd68.obs['Unique_ROI'].values,
                         name='CD68_Target_Gene_Score')

# ── ② 与配对 ROI 表合并成 fig7 主表 ──────────────────────────────────────────
paired_df = pd.read_csv(os.path.join(OUT, 'fig4_paired_ROI_16gene_LATAM.csv'))
score_df = score_series.reset_index()
score_df.columns = ['Unique_ROI', 'CD68_Target_Gene_Score']

master = paired_df.merge(score_df, on='Unique_ROI', how='inner')
print(f'Master DataFrame: {master.shape}')

master = master.rename(columns={'Unique_ROI': 'ROI_ID',
                                'B_HMGB1': 'CD20_Ligand',
                                'Mac_HAVCR2': 'CD68_Receptor',
                                'Mac_MAFB': 'CD68_TF'})
master_clean = master[['ROI_ID', 'Tissue_Type', 'CD20_Ligand',
                       'CD68_Receptor', 'CD68_TF',
                       'CD68_Target_Gene_Score']].copy()
print(master_clean.describe().to_string())

master_clean.to_csv(os.path.join(OUT, 'fig7_master_spatial_coevolution.csv'),
                    index=False)
print(f'[OK] fig7_master_spatial_coevolution.csv  ({len(master_clean)} ROIs)')
