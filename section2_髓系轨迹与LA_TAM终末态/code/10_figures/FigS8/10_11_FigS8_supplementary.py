# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import anndata as ad
A = ad.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac_subtyped.h5ad'))
print('形状:', A.shape)                    # (细胞数, 基因数) —— 需要全转录组
print('\nobs列名:', A.obs.columns.tolist())
print('\n前3行obs:')
print(A.obs.head(3).to_string())
print('\nX类型:', type(A.X), '| 是否稀疏:', hasattr(A.X, 'toarray'))



# ══════════════════════════════════════════════════════════════════
# Fig S8 — External LAM signature benchmarking of LA_TAM
# Requires: adata_mac_subtyped.h5ad, scanpy>=1.9, matplotlib, scipy
# ══════════════════════════════════════════════════════════════════
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

mpl.rcParams['svg.fonttype'] = 'none'
OUT = translate(r'D:\bulk-download\GSE182434')

# ── Step 0: load & sanity-check normalization ────────────────────
A = ad.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac_subtyped.h5ad'))
subtypes = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
A.obs['mac_subtype'] = pd.Categorical(A.obs['mac_subtype'], categories=subtypes, ordered=True)

# 粗判 X 是否 raw counts（若最大值远超 log-normalized 范围则先归一化）
xmax = A.X.max() if not hasattr(A.X, 'toarray') else A.X.max()
if xmax > 30:
    print(f'Max X = {xmax:.1f} → 判定为 raw counts，执行 CP10K + log1p')
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
else:
    print(f'Max X = {xmax:.2f} → 已是 log-normalized')

# ── Step 1: external signatures ──────────────────────────────────
sig_jaitin = ['TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5',
              'LGALS1','LGALS3','CD9','CD36']                      # Jaitin 2019, Cell
sig_klooster = ['GPNMB','FABP5','HMOX1','SPP1','ARG1']             # Kloosterman 2024, Cell
signatures = {'LAM_Jaitin2019': sig_jaitin, 'LLM_Kloosterman2024': sig_klooster}

# ── Step 2: detection filter + overlap bookkeeping ───────────────
det = pd.Series(np.asarray((A.X > 0).mean(axis=0)).ravel(), index=A.var_names)
la_tam_defining = ['PTGDS','CCL18','APOE','CHI3L1','CTSD','GPNMB','APOC1',
                   'PLA2G2D','CAPG','MMP9','NUPR1','CTSL','RARRES1',
                   'IL32','FUCA1','LGMN','FTL']                     # Dai 2024 LA_TAM markers

records = []
for name, genes in signatures.items():
    present = [g for g in genes if g in A.var_names]
    kept    = [g for g in present if det[g] > 0]
    overlap = [g for g in kept if g in la_tam_defining]
    records.append({'signature': name, 'n_input': len(genes),
                    'n_in_var': len(present), 'n_detectable': len(kept),
                    'missing': [g for g in genes if g not in A.var_names],
                    'zero_detect': [g for g in present if g not in kept],
                    'overlap_with_LAm_defining': overlap})
    sc.tl.score_genes(A, gene_list=kept, score_name=f'ext_{name}',
                      ctrl_size=min(50, len(kept)), n_bins=25,
                      use_raw=False, random_state=0)
overlap_df = pd.DataFrame(records)
overlap_df.to_csv(f'{OUT}/S8_signature_filtering.csv', index=False)
print(overlap_df.to_string())

# ── Step 3: statistics — LA_TAM vs all others (one-sided MWU) ────
stat_rows = []
for name in signatures:
    col = f'ext_{name}'
    la = A.obs.loc[A.obs.mac_subtype == 'LA_TAM', col].values
    ot = A.obs.loc[A.obs.mac_subtype != 'LA_TAM', col].values
    u, p = stats.mannwhitneyu(la, ot, alternative='greater')
    # AUROC = U / (n1*n2)
    auc = u / (len(la) * len(ot))
    d = (la.mean() - ot.mean()) / np.sqrt((la.std(ddof=1)**2 + ot.std(ddof=1)**2) / 2)
    stat_rows.append({'signature': name, 'median_LA_TAM': np.median(la),
                      'median_others': np.median(ot), 'U': u, 'p_value': p,
                      'AUROC': auc, 'cohens_d': d})
    # Kruskal–Wallis across all 5 states (context)
    kw, kwp = stats.kruskal(*[A.obs.loc[A.obs.mac_subtype == s, col] for s in subtypes])
    stat_rows[-1].update({'KruskalWallas_H': kw, 'KW_p': kwp})
stat_df = pd.DataFrame(stat_rows)
stat_df.to_csv(f'{OUT}/S8_statistics.csv', index=False)
print(stat_df.round(4).to_string())

# ── Step 3b: sensitivity — drop overlap genes, re-score ──────────
sens_rows = []
for name, genes in signatures.items():
    ov = [g for g in genes if g in la_tam_defining and det.get(g, 0) > 0]
    kept = [g for g in genes if g in A.var_names and det[g] > 0 and g not in la_tam_defining]
    sc.tl.score_genes(A, gene_list=kept, score_name=f'ext_{name}_sens',
                      ctrl_size=min(50, len(kept)), n_bins=25, random_state=0)
    col = f'ext_{name}_sens'
    la = A.obs.loc[A.obs.mac_subtype == 'LA_TAM', col].values
    ot = A.obs.loc[A.obs.mac_subtype != 'LA_TAM', col].values
    u, p = stats.mannwhitneyu(la, ot, alternative='greater')
    sens_rows.append({'signature': name, 'dropped': ov, 'p_sensitivity': p,
                      'AUROC_sens': u / (len(la) * len(ot))})
sens_df = pd.DataFrame(sens_rows)
sens_df.to_csv(f'{OUT}/S8_sensitivity.csv', index=False)
print(sens_df.to_string())

# ── Step 4: QC check — MAFB & CD163 per subtype ──────────────────
qc_genes = [g for g in ['MAFB', 'CD163', 'HAVCR2', 'HMGB1'] if g in A.var_names]
qc = A.obs.groupby('mac_subtype', observed=True).apply(
    lambda df: pd.Series({g: np.asarray(A[df.index, g].X.todense()).mean() for g in qc_genes}
))
qc.to_csv(f'{OUT}/S8_MAFB_CD163_check.csv')
print(qc.round(3).to_string())

# ── Step 5: Figure S8 ────────────────────────────────────────────
sub_cols = {'Mono': '#2980B9', 'DC_1': '#9B59B6', 'LA_TAM': '#C0392B',
            'IFN_TAM': '#F39C12', 'DC_2': '#E74C3C'}
sig_titles = {'LAM_Jaitin2019': 'Jaitin LAM signature (2019, Cell)',
              'LLM_Kloosterman2024': 'Lipid-laden macrophage signature (Kloosterman 2024, Cell)'}

fig = plt.figure(figsize=(12, 8))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.3)

# Panels a–b: violin per signature
for i, name in enumerate(signatures):
    ax = fig.add_subplot(gs[0, i])
    data = [A.obs.loc[A.obs.mac_subtype == s, f'ext_{name}'].values for s in subtypes]
    parts = ax.violinplot(data, positions=range(5), showextrema=False, widths=0.85)
    for pc, s in zip(parts['bodies'], subtypes):
        pc.set_facecolor(sub_cols[s]); pc.set_alpha(0.75 if s == 'LA_TAM' else 0.4)
    for j, d in enumerate(data):
        ax.hlines(np.median(d), j - 0.28, j + 0.28, color='k', lw=1.4)
        ax.scatter(np.full(len(d), j) + np.random.uniform(-0.15, 0.15, len(d)), d,
                   s=2, color='k', alpha=0.15, zorder=0)
    pv = stat_df.loc[stat_df.signature == name, 'p_value'].iloc[0]
    auc = stat_df.loc[stat_df.signature == name, 'AUROC'].iloc[0]
    ax.set_title(f'{sig_titles[name]}\nLA_TAM vs others: p = {pv:.1e}, AUROC = {auc:.2f}',
                 fontsize=9)
    ax.set_xticks(range(5)); ax.set_xticklabels(subtypes, rotation=30, ha='right')
    ax.set_ylabel('Signature score' if i == 0 else '')
    ax.spines[['top', 'right']].set_visible(False)

# Panel c: dot plot of shared/Jaitin LAM genes across subtypes
genes_dot = [g for g in sig_jaitin if g in A.var_names and det[g] > 0]
axc = fig.add_subplot(gs[1, :])
rows = []
for s in subtypes:
    mask = (A.obs.mac_subtype == s).values
    sub = A[mask, genes_dot].X
    sub = sub.toarray() if hasattr(sub, 'toarray') else sub
    for j, g in enumerate(genes_dot):
        rows.append({'subtype': s, 'gene': g, 'mean_exp': sub[:, j].mean(),
                     'pct': (sub[:, j] > 0).mean() * 100})
dd = pd.DataFrame(rows)
for g in genes_dot:
    v = dd.loc[dd.gene == g, 'mean_exp']
    dd.loc[dd.gene == g, 'mean_norm'] = (v - v.min()) / (v.max() - v.min() + 1e-9)
sc_map = axc.scatter(dd.gene, dd.subtype, s=dd.pct / 100 * 380 + 8,
                     c=dd.mean_norm, cmap='YlOrRd', edgecolors='#444', lw=0.4)
plt.colorbar(sc_map, ax=axc, label='Mean expression (scaled)', shrink=0.7)
axc.set_xticklabels(genes_dot, rotation=45, ha='right', fontsize=9)
axc.invert_yaxis()
axc.set_title('LAM signature gene expression across myeloid states', fontsize=10)
axc.spines[['top', 'right']].set_visible(False)

for fmt in ['png', 'svg']:
    fig.savefig(f'{OUT}/FigS8.{fmt}', dpi=300, bbox_inches='tight', format=fmt)
plt.close(fig)
print(f'✓ FigS8 saved to {OUT}')
