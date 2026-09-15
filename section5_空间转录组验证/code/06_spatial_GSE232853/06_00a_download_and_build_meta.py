# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  GSE232853 GeoMx DSP —— 下载原始数据 + 构建样本元数据 + 严格 QC
#                                                     （section5 上游 · 第 1 步）
#
#  来源：Biomni 平台分析记录 `14_空间转录组分析.md`
#        · L10733-10795  GEO 文件清单与下载
#        · L11025-11096  构建 meta_df（解析列名）
#        · L13621-13681  严格 QC（剔除 Full ROI）+ 抽出 Tissue_Type
#
#  数据来源（GEO 公开，无需注册）：
#        GSE232853_Processed_data_CD20_CD68_final.csv.gz
#        ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232853/suppl/
#
#  实验设计：GeoMx DSP，同一张切片上相邻 ROI 分别用 CD20 / CD68 抗体掩膜
#        （Mask），因此每个 ROI 有 CD20+ 与 CD68+ 两个配对样本。
#        原始矩阵 16 560 基因 × 578 样本，其中 6 个是 "Full ROI"（无掩膜）需剔除。
#
#  产出（写入 <数据根>/GSE232853_v2/）：
#        meta_df.csv               578 样本完整元数据（含 Full ROI）
#        meta_filtered.csv         572 样本（CD20 253 + CD68 319）
#        expression_raw_filtered.csv   16 560 基因 × 572 样本的原始计数
#
#  环境：pandas / numpy（下载用 urllib，无需 GEOparse）
# =============================================================================

import os
import urllib.request
import numpy as np
import pandas as pd

# ── 路径 ──────────────────────────────────────────────────────────────────────
GEO_URL = ('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232853/suppl/'
           'GSE232853_Processed_data_CD20_CD68_final.csv.gz')
OUT = translate(r'D:\bulk-download\GSE232853_v2')          # -> <数据根>/GSE232853_v2
RAW_DIR = translate(r'D:\bulk-download\external\GSE232853')  # 原始下载件
os.makedirs(OUT, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)
RAW = os.path.join(RAW_DIR, 'GSE232853_Processed_data_CD20_CD68_final.csv.gz')

# ── ① 下载（已存在则跳过）─────────────────────────────────────────────────────
if not os.path.exists(RAW):
    print(f'Downloading {GEO_URL} ...', flush=True)
    urllib.request.urlretrieve(GEO_URL, RAW)
    print('Done.', flush=True)
else:
    print(f'Already downloaded: {RAW}', flush=True)

df_raw = pd.read_csv(RAW, index_col=0)
print(f'Raw matrix: {df_raw.shape[0]} genes x {df_raw.shape[1]} samples')

# 列名不应重复；如有则去重（照原文处理）
dup = df_raw.columns[df_raw.columns.duplicated()].tolist()
if dup:
    print(f'Duplicate columns: {len(dup)} -> dropping')
    df_raw = df_raw.loc[:, ~df_raw.columns.duplicated()]

# ── ② 解析列名构建 meta_df ────────────────────────────────────────────────────
# 列名格式："Slide_ID | ROI_Num | Mask"
rows = []
for col in df_raw.columns:
    parts = col.split(' | ')
    if len(parts) == 3:
        slide_id, roi_num, mask = (p.strip() for p in parts)
    elif len(parts) == 2:
        slide_id, roi_num = (p.strip() for p in parts)
        mask = 'Unknown'
    else:
        slide_id, roi_num, mask = col, 'NA', 'Unknown'
    rows.append({'Sample_ID': col, 'Slide_ID': slide_id, 'ROI_Num': roi_num,
                 'Mask': mask, 'Unique_ROI': f'{slide_id}_{roi_num}'})

meta_df = pd.DataFrame(rows).set_index('Sample_ID')
print(f'meta_df: {meta_df.shape}')
print(meta_df['Mask'].value_counts().to_string())

# 目标基因必须在 panel 里
for g in ['HMGB1', 'HAVCR2', 'MAFB']:
    print(f'  {g} in panel: {g in df_raw.index}')
print(f'Missing values: {int(df_raw.isnull().sum().sum())}')
meta_df.to_csv(os.path.join(OUT, 'meta_df.csv'))

# ── ③ 严格 QC：只保留 CD20 / CD68 掩膜 ───────────────────────────────────────
keep = meta_df['Mask'].isin(['CD20', 'CD68'])
meta_filtered = meta_df[keep].copy()
print(f"\nRemoved {int((~keep).sum())} samples: "
      f"{meta_df.loc[~keep, 'Mask'].value_counts().to_dict()}")
print(f'Retained: {len(meta_filtered)} samples')

# ── ④ 从 Slide_ID 推断组织类型 ────────────────────────────────────────────────
def infer_tissue(slide_id: str) -> str:
    """含 tonsil 且不含 dlbcl → Normal；含 dlbcl → DLBCL；否则 Unknown。"""
    s = slide_id.lower()
    if 'tonsil' in s and 'dlbcl' not in s:
        return 'Normal'
    if 'dlbcl' in s:
        return 'DLBCL'
    return 'Unknown'


meta_filtered['Tissue_Type'] = meta_filtered['Slide_ID'].apply(infer_tissue)
unknown = meta_filtered[meta_filtered['Tissue_Type'] == 'Unknown']
if len(unknown):
    print(f"WARNING: {len(unknown)} samples with Unknown Tissue_Type: "
          f"{unknown['Slide_ID'].unique()}")
else:
    print('All samples have valid Tissue_Type.')

print('\nFinal distribution (Mask x Tissue_Type):')
print(meta_filtered.groupby(['Mask', 'Tissue_Type']).size().to_string())

# ── ⑤ 同步裁剪表达矩阵并保存 ──────────────────────────────────────────────────
df_filtered = df_raw[meta_filtered.index.tolist()]
print(f'\nFiltered expression: {df_filtered.shape[0]} genes x '
      f'{df_filtered.shape[1]} samples')
meta_filtered.to_csv(os.path.join(OUT, 'meta_filtered.csv'))
df_filtered.to_csv(os.path.join(OUT, 'expression_raw_filtered.csv'))
print('[OK] meta_df.csv / meta_filtered.csv / expression_raw_filtered.csv')
