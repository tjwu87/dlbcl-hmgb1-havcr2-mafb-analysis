# >>> 自动注入：统一路径配置 <<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT
# <<< 自动注入结束 >>>

# =============================================================================
#  pySCENIC 三步：GRNBoost2 → RcisTarget → AUCell，并导出论文所用 CSV
#                                                                   （section3 上游 · 第 2 步）
#
#  来源：Biomni 平台分析记录
#        · 19_PAGA与注释重做      —— GRNBoost2 / RcisTarget / AUCell 原始参数
#        · 09_CellRank命运与Regulon精修 —— 173 regulon 重算与下游 CSV 导出
#
#  参数（照原文，未改）：
#        GRNBoost2   seed=42，TF=出现在数据中的 834 个
#        RcisTarget  rank_threshold=1500, auc_threshold=0.05, nes_threshold=3.0,
#                    motif_similarity_fdr=0.001, weighted_recovery=False,
#                    filter_for_annotation=True, custom_multiprocessing x4
#        AUCell      auc_threshold=0.05, noweights=False
#
#  环境：pyscenic + arboreto + ctxcore + loompy + pandas + numpy
#  产出（写入 <数据根>/GSE182434/scenic/，即 04_01 下游脚本所读目录）：
#        data_AUC_matrix.csv                        (cells x regulons)
#        data_RSS_scores.csv                        (regulons x 3 亚型)
#        data_GRNBoost2_adjacencies.csv
#        data_RcisTarget_regulon_targets.csv
#        data_core_TF_AUC_and_receptor_expression.csv
#        data_TF_receptor_intersection.csv
#        data_pseudotime_spearman_correlations.csv
#        data_top30_RSS_per_subtype.csv
#        data_cell_metadata.csv
#
#  注：CLI 形式（pyscenic grn/ctx/aucell）见 04_00a_setup_scenic_db.sh 注释；
#      本脚本用 Python API，因为 CLI 的 dask 后端在容器内不稳定。
# =============================================================================

import os
import pickle
import numpy as np
import pandas as pd
import loompy
from scipy.stats import spearmanr
from scipy.spatial.distance import jensenshannon
import warnings
warnings.filterwarnings('ignore')

from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

np.object = object          # pyscenic 旧版本兼容补丁

GSE = translate(r'D:\bulk-download\GSE182434')
SCENIC_DB = os.environ.get('SCENIC_DB', os.path.join(GSE, 'scenic_db'))
OUT = os.path.join(GSE, 'scenic')
os.makedirs(OUT, exist_ok=True)

PATH_SUBTYPES = ['Mono', 'IFN_TAM', 'LA_TAM']

# ── ① 读入表达矩阵与 TF 列表 ──────────────────────────────────────────────────
with loompy.connect(os.path.join(SCENIC_DB, 'mono_ifntam_latam.loom'), 'r') as ds:
    cell_ids = list(ds.ca['CellID'])
    subtypes = list(ds.ca['mac_subtype'])
    genes = list(ds.ra['Gene'])
    raw = ds[:, :]                                   # genes x cells

expr = pd.DataFrame(raw.T, index=cell_ids, columns=genes)   # cells x genes
tf_all = [l.strip() for l in
          open(os.path.join(SCENIC_DB, 'allTFs_hg38.txt')).readlines() if l.strip()]
tfs_final = [t for t in tf_all if t in expr.columns]
print(f'Expression: {expr.shape}   TFs present: {len(tfs_final)} / {len(tf_all)}')

# ── ② GRNBoost2：TF→target 共表达网络 ────────────────────────────────────────
print('Running GRNBoost2 ...')
adj = grnboost2(expression_data=expr, tf_names=tfs_final, verbose=True, seed=42)
adj.to_csv(os.path.join(SCENIC_DB, 'adjacencies_new.tsv'), sep='\t', index=False)
print(f'Adjacencies: {len(adj):,} pairs')

# ── ③ 共表达模块 → RcisTarget motif 富集 → regulon ───────────────────────────
print('Building modules ...')
modules = list(modules_from_adjacencies(
    adjacencies=adj, ex_mtx=expr,
    rho_mask_dropouts=False, min_genes=5, keep_only_activating=True))
print(f'Modules: {len(modules)}')

db = RankingDatabase(
    fname=os.path.join(
        SCENIC_DB,
        'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'),
    name='hg38_10kbp')

print('Running RcisTarget (custom_multiprocessing) ...')
df_motifs = prune2df(
    rnkdbs=[db],
    modules=modules,
    motif_annotations_fname=os.path.join(
        SCENIC_DB, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'),
    rank_threshold=1500,
    auc_threshold=0.05,
    nes_threshold=3.0,
    motif_similarity_fdr=0.001,
    orthologuous_identity_threshold=0.0,
    weighted_recovery=False,
    client_or_address='custom_multiprocessing',
    num_workers=4,
    filter_for_annotation=True,
)
df_motifs.to_csv(os.path.join(SCENIC_DB, 'motifs_enriched_new.csv'))
print(f'Enriched motif-module pairs: {len(df_motifs)}')

regulons = df2regulons(df_motifs)
with open(os.path.join(SCENIC_DB, 'regulons_new.pkl'), 'wb') as f:
    pickle.dump(regulons, f)
regulon_targets = {r.name: sorted(r.genes) for r in regulons}
print(f'Motif-pruned regulons: {len(regulons)}')

# ── ④ AUCell：每细胞每 regulon 活性 ──────────────────────────────────────────
print('Running AUCell ...')
auc_mtx = aucell(exp_mtx=expr, signatures=regulons,
                 auc_threshold=0.05, num_workers=4, noweights=False)
auc_mtx.to_csv(os.path.join(SCENIC_DB, 'auc_matrix_new.csv'))
print(f'AUC matrix: {auc_mtx.shape}')

# ── ⑤ 导出下游 CSV ───────────────────────────────────────────────────────────
dom = [r.name.split('(')[0] for r in regulons]
auc_mtx.columns = [f'{t}(+)' for t in dom]
sub = pd.Series(subtypes, index=cell_ids)

def compute_rss(auc_matrix, cell_labels, cell_types):
    """RSS(regulon, celltype) = 1 - JS_divergence(regulon 分布, 细胞类型指示分布)"""
    out = {ct: {} for ct in cell_types}
    for reg in auc_matrix.columns:
        vals = auc_matrix[reg].values.astype(float)
        p = vals / (vals.sum() + 1e-10)
        for ct in cell_types:
            ind = (cell_labels == ct).astype(float)
            q = ind / (ind.sum() + 1e-10)
            out[ct][reg] = 1 - jensenshannon(p, q)
    return pd.DataFrame(out)

rss = compute_rss(auc_mtx, sub.values, PATH_SUBTYPES).dropna()
print(f'RSS: {rss.shape}')

# 细胞 × 核心 TF AUC + 受体表达 + 伪时间
CORE_TFS = ['ETV5', 'MITF', 'MAFB', 'XBP1', 'FOSB', 'KLF2', 'KLF4', 'IRF1', 'STAT1']
RECEPTORS = ['AXL', 'HAVCR2', 'CD74', 'HLA-DPA1']
core = pd.DataFrame(index=auc_mtx.index)
for tf in CORE_TFS:
    col = [c for c in auc_mtx.columns if c.split('(')[0] == tf]
    if col:
        core[f'AUC_{tf}'] = auc_mtx[col[0]]
for g in RECEPTORS:
    if g in expr.columns:
        core[g] = expr.loc[core.index, g]
core['mac_subtype'] = sub.loc[core.index]

# TF→受体 交集（regulon 内含受体基因者）
inter = []
for tf, targets in regulon_targets.items():
    for rec in RECEPTORS:
        if rec in targets:
            inter.append({'TF': tf, 'receptor': rec})
inter = pd.DataFrame(inter, columns=['TF', 'receptor'])

# 伪时间 Spearman（若 metadata 提供 dpt_pseudotime）
meta_p = os.path.join(GSE, 'trajectory', 'cell_metadata.csv')
if os.path.exists(meta_p):
    meta = pd.read_csv(meta_p, index_col=0)
    if 'dpt_pseudotime' in meta.columns:
        common = core.index.intersection(meta.index)
        rows = []
        pt = meta.loc[common, 'dpt_pseudotime'].values
        for c in core.columns:
            if c in ('mac_subtype',):
                continue
            rho, pv = spearmanr(pt, core.loc[common, c].values)
            rows.append({'feature': c, 'rho': rho, 'pval': pv})
        pd.DataFrame(rows).to_csv(
            os.path.join(OUT, 'data_pseudotime_spearman_correlations.csv'), index=False)

auc_mtx.to_csv(os.path.join(OUT, 'data_AUC_matrix.csv'))
rss.to_csv(os.path.join(OUT, 'data_RSS_scores.csv'))
adj.to_csv(os.path.join(OUT, 'data_GRNBoost2_adjacencies.csv'), index=False)
pd.DataFrame([{'TF': k, 'target': t} for k, v in regulon_targets.items()
              for t in v]).to_csv(
    os.path.join(OUT, 'data_RcisTarget_regulon_targets.csv'), index=False)
core.to_csv(os.path.join(OUT, 'data_core_TF_AUC_and_receptor_expression.csv'))
inter.to_csv(os.path.join(OUT, 'data_TF_receptor_intersection.csv'), index=False)
rss.apply(lambda s: s.sort_values(ascending=False).head(30).index.tolist()).to_csv(
    os.path.join(OUT, 'data_top30_RSS_per_subtype.csv'))
sub.to_frame('mac_subtype').to_csv(
    os.path.join(OUT, 'data_cell_metadata.csv'))

print(f'[OK] {len([f for f in os.listdir(OUT) if f.startswith("data_")])} '
      f'个 data_*.csv 已写入 {OUT}')
