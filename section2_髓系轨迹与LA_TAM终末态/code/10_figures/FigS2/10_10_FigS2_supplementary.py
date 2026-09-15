# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

# ══════════════════════════════════════════════════════════════════
# Fig S8 — External LAM signature benchmarking of LA_TAM (FINAL)
# 输出文件：
#   S8_signature_filtering.csv   过滤与重叠记账（含 per-gene 检出率）
#   S8_statistics.csv            two-sided MWU + BH + AUROC + d + KW
#   S8_sensitivity.csv           剔除重叠基因敏感性
#   S8_perpatient.csv            per-patient 明细（内部QC/回复信用，不入图）
#   S8_leave_one_patient_out.csv LOPO（图 panel d 的数据）
#   S8_rho.csv                   外部签名 vs 5标记 ρ（两个口径）
#   S8_MAFB_CD163_check.csv      内部QC（不入补充材料）
#   FigS8.png / FigS8.svg
# ══════════════════════════════════════════════════════════════════
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

mpl.rcParams['svg.fonttype'] = 'none'
OUT      = translate(r'D:\bulk-download\GSE182434')
DET_MIN  = 0.01
MARKERS5 = ['APOE','APOC1','ACP5','CCL18','CTSD']
MIN_LA   = 5                      # per-patient 可评估的最低 LA_TAM 细胞数
PAT_COL  = 'Patient'              # 已确认：529 细胞全部来自 4 位 DLBCL 患者
RNG      = np.random.default_rng(0)

def bh_adjust(pvals):
    p = np.asarray(pvals, float); n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1)
    return out

# ── Step 0: load & sanity-check ───────────────────────────────────
A = ad.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac_subtyped.h5ad'))
subtypes = ['Mono', 'DC_1', 'LA_TAM', 'IFN_TAM', 'DC_2']
A.obs['mac_subtype'] = pd.Categorical(A.obs['mac_subtype'],
                                      categories=subtypes, ordered=True)
X = A.X.toarray() if hasattr(A.X, 'toarray') else np.asarray(A.X)
assert X.min() >= 0, \
    f'X.min()={X.min():.3f} < 0 → 疑似 scaled 数据，请换回 log-normalized 层！'
xmax = X.max()
if xmax > 30:
    print(f'Max X = {xmax:.1f} → raw counts，执行 CP10K + log1p')
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
else:
    print(f'Max X = {xmax:.2f} → 已是 log-normalized')
assert A.obs[PAT_COL].nunique() == 4, '预期 4 位 DLBCL 患者，请检查数据'

# ── Step 1: external signatures ───────────────────────────────────
sig_jaitin   = ['TREM2','LIPA','LPL','CTSB','CTSL','FABP4','FABP5',
                'LGALS1','LGALS3','CD9','CD36']          # Jaitin 2019 — primary
sig_klooster = ['GPNMB','FABP5','HMOX1','SPP1','ARG1']   # Kloosterman 2024 — secondary
signatures = {'LAM_Jaitin2019': sig_jaitin,
              'LLM_Kloosterman2024': sig_klooster}

# ── Step 2: detection filter (>=1%) + overlap bookkeeping ────────
det = pd.Series(np.asarray((X > 0).mean(axis=0)).ravel(), index=A.var_names)
la_tam_defining = ['PTGDS','CCL18','APOE','CHI3L1','CTSD','GPNMB','APOC1',
                   'PLA2G2D','CAPG','MMP9','NUPR1','CTSL','RARRES1',
                   'IL32','FUCA1','LGMN','FTL']          # = Table S21 LA_TAM 行
kept_map, records = {}, []
for name, genes in signatures.items():
    present = [g for g in genes if g in A.var_names]
    kept    = [g for g in present if det[g] >= DET_MIN]
    kept_map[name] = kept
    records.append({'signature': name, 'n_input': len(genes),
                    'n_in_var': len(present), 'n_detectable': len(kept),
                    'missing': ';'.join(g for g in genes if g not in A.var_names),
                    'low_detect(<1%)': ';'.join(g for g in present if g not in kept),
                    'overlap_defining': ';'.join(g for g in kept
                                                 if g in la_tam_defining),
                    'per_gene_det': {g: round(float(det[g]), 4) for g in present}})
pd.DataFrame(records).to_csv(f'{OUT}/S8_signature_filtering.csv', index=False)
print(pd.DataFrame([{k: v for k, v in r.items() if k != 'per_gene_det'}
                    for r in records]).to_string())

for name, kept in kept_map.items():
    sc.tl.score_genes(A, gene_list=kept, score_name=f'ext_{name}',
                      ctrl_size=min(50, len(kept)), n_bins=25,
                      use_raw=False, random_state=0)

# ── Step 3: LA_TAM vs others — TWO-SIDED MWU + BH ─────────────────
stat_rows = []
for name in signatures:
    col = f'ext_{name}'
    la = A.obs.loc[A.obs.mac_subtype == 'LA_TAM', col].values
    ot = A.obs.loc[A.obs.mac_subtype != 'LA_TAM', col].values
    u, p = stats.mannwhitneyu(la, ot, alternative='two-sided')
    auc = u / (len(la) * len(ot))
    d   = (la.mean() - ot.mean()) / np.sqrt(
              (la.std(ddof=1)**2 + ot.std(ddof=1)**2) / 2)
    kw, kwp = stats.kruskal(*[A.obs.loc[A.obs.mac_subtype == s, col]
                              for s in subtypes])
    stat_rows.append({'signature': name, 'median_LA_TAM': np.median(la),
                      'median_others': np.median(ot), 'U': u,
                      'p_two_sided': p, 'p_BH': np.nan, 'AUROC': auc,
                      'cohens_d': d, 'KW_H': kw, 'KW_p': kwp,
                      'n_kept_genes': len(kept_map[name])})
stat_df = pd.DataFrame(stat_rows)
stat_df['p_BH'] = bh_adjust(stat_df['p_two_sided'].values)
stat_df.to_csv(f'{OUT}/S8_statistics.csv', index=False)
print(stat_df.round(6).to_string())

# ── Step 3b: sensitivity — drop overlap genes, re-score ──────────
sens_rows = []
for name, genes in signatures.items():
    dropped = [g for g in genes if g in la_tam_defining
               and g in A.var_names and det[g] >= DET_MIN]
    kept_s  = [g for g in genes if g in A.var_names
               and det[g] >= DET_MIN and g not in la_tam_defining]
    sc.tl.score_genes(A, gene_list=kept_s, score_name=f'ext_{name}_sens',
                      ctrl_size=min(50, len(kept_s)), n_bins=25, random_state=0)
    col = f'ext_{name}_sens'
    la = A.obs.loc[A.obs.mac_subtype == 'LA_TAM', col].values
    ot = A.obs.loc[A.obs.mac_subtype != 'LA_TAM', col].values
    u, p = stats.mannwhitneyu(la, ot, alternative='two-sided')
    sens_rows.append({'signature': name, 'dropped': ';'.join(dropped),
                      'n_kept_sens': len(kept_s), 'kept_sens': ';'.join(kept_s),
                      'p_sens_two_sided': p, 'AUROC_sens': u/(len(la)*len(ot))})
sens_df = pd.DataFrame(sens_rows)
sens_df.to_csv(f'{OUT}/S8_sensitivity.csv', index=False)
print(sens_df.to_string())

# ── Step 3c: per-patient 明细（内部QC/回复信用，不入图）────────────
col = 'ext_LAM_Jaitin2019'
pp_rows = []
for pat, sub in A.obs.groupby(PAT_COL, observed=True):
    la = sub.loc[sub.mac_subtype == 'LA_TAM', col].values
    ot = sub.loc[sub.mac_subtype != 'LA_TAM', col].values
    p  = np.nan
    if len(la) >= 3 and len(ot) >= 10:
        _, p = stats.mannwhitneyu(la, ot, alternative='two-sided')
    pp_rows.append({'patient': pat, 'n_total': len(sub),
                    'n_LA_TAM': len(la), 'n_others': len(ot),
                    'median_LA_TAM': np.median(la) if len(la) else np.nan,
                    'median_others': np.median(ot) if len(ot) else np.nan,
                    'direction': ('LA_TAM_higher' if len(la) and
                                  np.median(la) > np.median(ot)
                                  else 'not_evaluable_or_not_higher'),
                    'p': p})
pp_df = pd.DataFrame(pp_rows)
pp_df.to_csv(f'{OUT}/S8_perpatient.csv', index=False)
print(pp_df.round(4).to_string())
evaluable = pp_df[pp_df['n_LA_TAM'] >= MIN_LA]
print(f'\n可评估患者(≥{MIN_LA} LA_TAM): {list(evaluable.patient)}；'
      f'方向一致 {(evaluable.direction=="LA_TAM_higher").sum()}/{len(evaluable)}')
print('不可评估患者:', pp_df[pp_df['n_LA_TAM'] < MIN_LA].patient.tolist())

# ── Step 3d: external score vs 5-marker ρ（两个口径）──────────────
mk = [g for g in MARKERS5 if g in A.var_names]
Xm = A[:, mk].X.toarray() if hasattr(A[:, mk].X, 'toarray') else np.asarray(A[:, mk].X)
score5 = Xm.mean(axis=1)
rho_rows = []
for name in signatures:
    s = A.obs[f'ext_{name}'].values
    r1, p1 = stats.spearmanr(s, score5)
    m = (A.obs.mac_subtype != 'LA_TAM').values
    r2, p2 = stats.spearmanr(s[m], score5[m])
    rho_rows.append({'signature': name, 'markers': ';'.join(mk),
                     'rho_all': r1, 'p_all': p1,
                     'rho_nonLATAM': r2, 'p_nonLATAM': p2})
rho_df = pd.DataFrame(rho_rows)
rho_df.to_csv(f'{OUT}/S8_rho.csv', index=False)
print(rho_df.round(4).to_string())
# 引用口径：图注/回复信优先用 rho_nonLATAM（防同义反复）

# ── Step 3e: leave-one-patient-out（panel d 的数据，版本B核心证据）──
lopo_rows = []
for pat in A.obs[PAT_COL].unique():
    m = (A.obs[PAT_COL] != pat).values
    la = A.obs.loc[m & (A.obs.mac_subtype == 'LA_TAM'), col].values
    ot = A.obs.loc[m & (A.obs.mac_subtype != 'LA_TAM'), col].values
    u, p = stats.mannwhitneyu(la, ot, alternative='two-sided')
    lopo_rows.append({'excluded_patient': pat, 'n_LA_TAM_left': len(la),
                      'AUROC': u/(len(la)*len(ot)), 'p': p})
lopo_df = pd.DataFrame(lopo_rows)
lopo_df.to_csv(f'{OUT}/S8_leave_one_patient_out.csv', index=False)
print('\n=== LOPO（关键判定：剔除任一患者后是否仍显著）===')
print(lopo_df.round(6).to_string())

# ── Step 4: QC — MAFB/CD163（仅内部核对，不入补充材料）────────────
qc_genes = [g for g in ['MAFB','CD163','HAVCR2','HMGB1'] if g in A.var_names]
qc = pd.DataFrame({g: A.obs.groupby('mac_subtype', observed=True).apply(
        lambda df, g=g: X[A.obs.index.isin(df.index),
                          A.var_names.get_loc(g)].mean())
        for g in qc_genes})
qc.to_csv(f'{OUT}/S8_MAFB_CD163_check.csv')
print(qc.round(3).to_string())

# ── Step 5: Figure S8（a/b violin + c dot plot + d LOPO + e ρ）────
sub_cols = {'Mono': '#2980B9', 'DC_1': '#9B59B6', 'LA_TAM': '#C0392B',
            'IFN_TAM': '#F39C12', 'DC_2': '#E74C3C'}
sig_titles = {'LAM_Jaitin2019':
              'a  Jaitin LAM signature (2019, Cell) — primary',
              'LLM_Kloosterman2024':
              'b  Lipid-laden macrophage signature (2024, Cell) — complementary'}

fig = plt.figure(figsize=(13, 14))
gs  = fig.add_gridspec(3, 2, hspace=0.5, wspace=0.28)

# (a,b) violins
for i, name in enumerate(signatures):
    ax = fig.add_subplot(gs[0, i])
    data = [A.obs.loc[A.obs.mac_subtype == s, f'ext_{name}'].values
            for s in subtypes]
    parts = ax.violinplot(data, positions=range(5), showextrema=False,
                          widths=0.85)
    for pc, s in zip(parts['bodies'], subtypes):
        pc.set_facecolor(sub_cols[s]); pc.set_alpha(0.8 if s=='LA_TAM' else 0.4)
    for j, d_ in enumerate(data):
        ax.hlines(np.median(d_), j-0.28, j+0.28, color='k', lw=1.4)
        ax.scatter(np.full(len(d_), j) + RNG.uniform(-0.15, 0.15, len(d_)),
                   d_, s=2, color='k', alpha=0.12, zorder=0)
    r = stat_df[stat_df.signature == name].iloc[0]
    ax.set_title(f"{sig_titles[name]}\np = {r['p_two_sided']:.1e} "
                 f"(BH {r['p_BH']:.1e}), AUROC = {r['AUROC']:.2f}", fontsize=9)
    ax.set_xticks(range(5)); ax.set_xticklabels(subtypes, rotation=30, ha='right')
    ax.set_ylabel('Signature score' if i == 0 else '')
    ax.spines[['top','right']].set_visible(False)

# (c) dot plot — 两个签名基因 + 分块标注
genes_dot_j = [g for g in sig_jaitin if g in A.var_names and det[g] >= DET_MIN]
genes_dot_k = [g for g in sig_klooster if g in A.var_names and det[g] >= DET_MIN]
genes_all = genes_dot_j + genes_dot_k
axc = fig.add_subplot(gs[1, :])
rows = []
for s in subtypes:
    m = (A.obs.mac_subtype == s).values
    sub = A[m, genes_all].X
    sub = sub.toarray() if hasattr(sub, 'toarray') else sub
    for j, g in enumerate(genes_all):
        rows.append({'subtype': s, 'gene': g, 'mean_exp': sub[:, j].mean(),
                     'pct': (sub[:, j] > 0).mean()*100})
dd = pd.DataFrame(rows)
for g in genes_all:
    v = dd.loc[dd.gene == g, 'mean_exp']
    dd.loc[dd.gene == g, 'mean_norm'] = (v-v.min())/(v.max()-v.min()+1e-9)
sc_map = axc.scatter(dd.gene, dd.subtype, s=dd.pct/100*380+8, c=dd.mean_norm,
                     cmap='YlOrRd', edgecolors='#444', lw=0.4)
plt.colorbar(sc_map, ax=axc, label='Mean expression (scaled)', shrink=0.7)
axc.axvline(len(genes_dot_j)-0.5, color='#888', ls='--', lw=0.8)
axc.text(len(genes_dot_j)/2-0.5, -0.85, 'Jaitin LAM genes',
         ha='center', fontsize=8, style='italic')
axc.text(len(genes_dot_j)+len(genes_dot_k)/2-0.5, -0.85,
         'Kloosterman LLM genes', ha='center', fontsize=8, style='italic')
axc.set_xticks(range(len(genes_all)))
axc.set_xticklabels(genes_all, rotation=45, ha='right', fontsize=9)
axc.invert_yaxis(); axc.set_ylim(4.9, -1.4)
axc.set_title('c  External signature gene expression across myeloid states '
              '(dot size = % detected)', fontsize=10)
axc.spines[['top','right']].set_visible(False)

# (d) LOPO — bar of AUROC per excluded patient
axd = fig.add_subplot(gs[2, 0])
lopo_sorted = lopo_df.set_index('excluded_patient').loc[
    sorted(lopo_df.excluded_patient)]
bars = axd.bar(range(len(lopo_sorted)), lopo_sorted['AUROC'],
               color='#5DADE2', edgecolor='#2E86C1')
for j, (idx, r) in enumerate(lopo_sorted.iterrows()):
    star = '***' if r['p'] < 0.001 else ('**' if r['p'] < 0.01
           else ('*' if r['p'] < 0.05 else 'ns'))
    axd.text(j, r['AUROC'] + 0.01, f"{star}\np={r['p']:.1e}",
             ha='center', fontsize=7)
axd.axhline(0.5, color='#C0392B', ls=':', lw=1)
axd.text(len(lopo_sorted)-0.5, 0.505, 'chance', color='#C0392B',
         fontsize=7, ha='right')
axd.set_xticks(range(len(lopo_sorted)))
axd.set_xticklabels(lopo_sorted.index, rotation=20, ha='right', fontsize=8)
axd.set_ylim(0.4, 1.0)
axd.set_ylabel('AUROC (LA-TAM vs others)')
axd.set_title('d  Leave-one-patient-out sensitivity\n'
              '(Jaitin LAM score; each bar excludes one patient)', fontsize=9)
axd.spines[['top','right']].set_visible(False)

# (e) ρ — Jaitin score vs 5-marker mean
axe = fig.add_subplot(gs[2, 1])
axe.scatter(A.obs[col].values, score5, s=4, alpha=0.3, color='#566573')
rr = rho_df[rho_df.signature == 'LAM_Jaitin2019'].iloc[0]
axe.set_title(f"e  Jaitin LAM score vs 5-marker LA-TAM definition\n"
              f"ρ(all cells) = {rr['rho_all']:.2f}; "
              f"ρ(non-LA-TAM) = {rr['rho_nonLATAM']:.2f}", fontsize=9)
axe.set_xlabel('Jaitin LAM signature score')
axe.set_ylabel('Mean of APOE/APOC1/ACP5/\nCCL18/CTSD (log-normalized)')
axe.spines[['top','right']].set_visible(False)

for fmt in ['png', 'svg']:
    fig.savefig(f'{OUT}/FigS8.{fmt}', dpi=300, bbox_inches='tight', format=fmt)
plt.close(fig)
print('✓ FigS8 final saved')
