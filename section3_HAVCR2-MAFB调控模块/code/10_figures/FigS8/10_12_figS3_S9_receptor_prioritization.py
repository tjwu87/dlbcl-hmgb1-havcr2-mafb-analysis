# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr

ad   = sc.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac.h5ad'))
core = pd.read_csv(translate(r'D:\bulk-download\GSE182434\scenic\data_core_TF_AUC_and_receptor_expression.csv'), index_col=0)
auc  = pd.read_csv(translate(r'D:\bulk-download\GSE182434\scenic\data_AUC_matrix.csv'), index_col=0)

# ---------- 直接按原名对齐，不做任何映射 ----------
common = core.index.intersection(ad.obs.index)
print(f'直接匹配: {len(common)}/{len(core)}')          # 应为 272/272
print('ad.obs的列:', ad.obs.columns.tolist())            # 看有没有subtype列

# ---------- 提取四个受体 ----------
genes = ['HAVCR2','TLR2','TLR4','AGER']
X = ad[common, genes].to_df()

# ---------- ★质控：HAVCR2必须复现 ----------
hav_ref = core['expr_HAVCR2']
r, p = spearmanr(X['HAVCR2'], hav_ref.loc[common])
print(f'质控 spearman={r:.4f}, 最大绝对差={np.abs(X["HAVCR2"]-hav_ref.loc[common]).max():.4f}')
assert r > 0.999, '口径不一致，先看下面的提示'

# ---------- 组装统计用表 ----------
df = X.copy()
df['subtype'] = core.loc[common, 'subtype'] if 'subtype' in core.columns else ad.obs.loc[common, 'subtype']
df['pt']     = core.loc[common, 'dpt_pseudotime']
df['MAFB']   = auc['MAFB(+)'].reindex(common).values

# ---------- S19三个维度的统计 ----------
print('\n== 1) 亚型均值 ==')
print(df.groupby('subtype')[genes].mean().round(2))

print('\n== 2) vs pseudotime ==')
for g in genes:
    r, p = spearmanr(df[g], df['pt'])
    print(f'{g:8s} rho={r:+.3f}, p={p:.2g}')

print('\n== 3) vs MAFB(+) AUC ==')
for g in genes:
    r, p = spearmanr(df[g], df['MAFB'])
    print(f'{g:8s} rho={r:+.3f}, p={p:.2g}')

meta = pd.read_csv(translate(r'D:\bulk-download\GSE182434\cell_metadata.csv'), index_col=0)
print(len(meta))
# 以下为早期手工 QC：仅当该版 cell_metadata.csv 同时含 expr_HAVCR2 与 dpt_pseudotime 时才校验。
# 论文最终版数据里 cell_metadata.csv 无 dpt_pseudotime 列（伪时序存于 trajectory/ 下），
# 故此处改为条件执行，避免历史断言误伤复现。
if {'expr_HAVCR2', 'dpt_pseudotime'} <= set(meta.columns):
    r, p = spearmanr(meta['expr_HAVCR2'], meta['dpt_pseudotime'])
    print(f'[QC] cell_metadata expr_HAVCR2 vs dpt_pseudotime: rho={r:+.4f}, p={p:.2g}')
else:
    print('[QC] 跳过：cell_metadata.csv 不含 expr_HAVCR2/dpt_pseudotime（该列在 trajectory/ 伪时序表中）')

ad0 = sc.read_h5ad(translate(r'D:\bulk-download\GSE182434\adata_mac.h5ad'))
common0 = core.index.intersection(ad0.obs.index)
print(f'adata_mac.h5ad 匹配: {len(common0)}/{len(core)}')
X0 = ad0[common0, ['HAVCR2']].to_df()
r0, _ = spearmanr(X0['HAVCR2'], core['expr_HAVCR2'].loc[common0])
print(f'adata_mac.h5ad 质控: rho={r0:.4f}')


# ============================================================
# Fig S9: Prioritization of candidate HMGB1 receptors
# Panel A: quadrant scatter — pseudotime rho vs MAFB(+)-AUC rho
#          (point size = LA_TAM/Mono expression ratio)
# Panel B: differential engagement magnitude (Table S2)
# 数据来源: adata_mac.h5ad 272-cell trajectory core + Table S2
# ============================================================
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 9,
    'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# ---------- 数据（与 Table S19 一一对应） ----------
recs = ['HAVCR2', 'TLR2', 'TLR4', 'RAGE/AGER']
rho_pt    = [0.235, -0.122, -0.007, -0.144]   # expression vs pseudotime
rho_mafb  = [0.277,  0.020,  0.037,  0.015]   # expression vs MAFB(+) AUC
pvals     = [9.1e-05, 0.044, 0.91, 0.017]
expr_LA   = [0.34, 0.44, 0.22, 0.01]          # log-normalized mean, LA_TAM
expr_Mono = [0.16, 0.60, 0.31, 0.06]          # log-normalized mean, Mono
delta     = [0.064, 0.077, 0.086, np.nan]     # Table S2 (AGER not detected)

# 点面积：LA_TAM/Mono 表达倍数（给个底面积，避免 0.01/0.06 太小）
ratio = np.array(expr_LA) / np.array(expr_Mono)
size  = 150 * (ratio / ratio.max()) + 40      # 40–190 pt^2

# 颜色：HAVCR2 高亮，其余灰阶
colors   = ['#D62728', '#7F7F7F', '#7F7F7F', '#7F7F7F']
edgecols = ['#D62728', '#4D4D4D', '#4D4D4D', '#4D4D4D']

# ---------- 画布 ----------
fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(7.2, 3.2),
    gridspec_kw={'width_ratios': [1.5, 1]})

# ---------- Panel A: 四象限散点 ----------
axA.axhline(0, color='#BBBBBB', lw=0.8, ls='--', zorder=1)
axA.axvline(0, color='#BBBBBB', lw=0.8, ls='--', zorder=1)

for i, r in enumerate(recs):
    axA.scatter(rho_pt[i], rho_mafb[i], s=size[i],
                c=colors[i], edgecolors=edgecols[i], linewidth=1.2,
                alpha=0.9 if i == 0 else 0.55, zorder=3 if i == 0 else 2)
    # 标签位置微调，避免重叠
    dx, dy = {
        'HAVCR2':   (0.015,  0.020),
        'TLR2':     (0.015, -0.030),
        'TLR4':     (0.015,  0.025),
        'RAGE/AGER':(-0.020, -0.035),
    }[r]
    ha = 'left' if dx > 0 else 'right'
    axA.text(rho_pt[i] + dx, rho_mafb[i] + dy, r,
             fontsize=9, ha=ha, va='center',
             fontweight='bold' if i == 0 else 'normal',
             color=colors[i])

# 右上象限标注
axA.text(0.145, 0.26, 'LA_TAM-enriched\n& MAFB-coupled',
         fontsize=7.5, style='italic', color='#D62728',
         ha='center', va='center', alpha=0.8)

axA.set_xlabel('Expression vs pseudotime\n(Spearman $\\rho$)', fontsize=9)
axA.set_ylabel('Expression vs MAFB(+) regulon AUC\n(Spearman $\\rho$)', fontsize=9)
axA.set_xlim(-0.22, 0.33)
axA.set_ylim(-0.08, 0.33)
axA.text(-0.03, 1.02, 'a', transform=axA.transAxes,
         fontsize=12, fontweight='bold', va='bottom')

# 点大小说明（图内小图例）
for s_lab, r_lab, pos in [(size[0], ratio[0], (0.30, -0.065)),
                          (size[2], ratio[2], (0.30, -0.045))]:
    pass  # 如需图例可解开下面注释
# 图例：LA_TAM/Mono expression ratio
h1 = axA.scatter([], [], s=size[0], c='grey', alpha=0.5, edgecolors='k', lw=0.8,
                 label=f'LA_TAM/Mono\nexpression ratio\n({ratio[0]:.1f}$\\times$, HAVCR2)')
h2 = axA.scatter([], [], s=size[3], c='grey', alpha=0.5, edgecolors='k', lw=0.8,
                 label=f'({ratio[3]:.1f}$\\times$, RAGE/AGER)')
axA.legend(handles=[h1, h2], loc='lower left', fontsize=6.5,
           frameon=False, handletextpad=0.8, borderpad=0.2)

# ---------- Panel B: Δ engagement magnitude ----------
d_ok  = [d for d in delta if not np.isnan(d)]
r_ok  = [r for r, d in zip(recs, delta) if not np.isnan(d)]
c_ok  = [c for c, d in zip(colors, delta) if not np.isnan(d)]
ypos = np.arange(len(d_ok))[::-1]

axB.barh(ypos, d_ok, color=c_ok, alpha=0.85, edgecolor='white', height=0.55)
for y, d in zip(ypos, d_ok):
    axB.text(d + 0.002, y, f'{d:.3f}', va='center', fontsize=8)

axB.set_yticks(ypos)
axB.set_yticklabels(r_ok, fontsize=9)
axB.set_xlabel('Differential engagement $\\Delta$\n(malignant vs normal B, Table S2)', fontsize=9)
axB.set_xlim(0, 0.105)
axB.text(-0.18, 1.02, 'b', transform=axB.transAxes,
         fontsize=12, fontweight='bold', va='bottom')

#AGER 注记
axB.text(0.052, ypos[0] if len(ypos) > 0 else 0, '', fontsize=1)
axB.text(0.05, ypos[-1] - 0.35, 'RAGE/AGER:\nnot detected\namong prioritized pairs',
         fontsize=7, style='italic', color='#7F7F7F', va='center')

# 去顶右边框
for ax in (axA, axB):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=0.8, length=3)

fig.suptitle('Candidate HMGB1 receptor prioritization across '
             'expression, trajectory, and regulon-coupling dimensions',
             fontsize=9.5, y=1.02, style='italic')

plt.tight_layout()
plt.savefig('Fig_S9_receptor_prioritization.png', dpi=300, bbox_inches='tight')
plt.savefig('Fig_S9_receptor_prioritization.pdf',       bbox_inches='tight')
plt.savefig('Fig_S9_receptor_prioritization.svg',       bbox_inches='tight')
plt.show()
print('输出: Fig_S9_receptor_prioritization.png / .pdf / .svg')

# ============================================================
# Fig S3: Prioritization of candidate HMGB1 receptors
# 配色对齐原有图系统: REC_COLORS['HAVCR2'] = #1F77B4
# Panel A: quadrant scatter; Panel B: Δ engagement bars
# ============================================================
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 9,
    'axes.linewidth': 0.8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# ---------- 数据（与 Table S19 一致） ----------
recs      = ['HAVCR2', 'TLR2', 'TLR4', 'RAGE/AGER']
rho_pt    = [0.235, -0.122, -0.007, -0.144]
rho_mafb  = [0.277,  0.020,  0.037,  0.015]
expr_LA   = [0.34, 0.44, 0.22, 0.01]
expr_Mono = [0.16, 0.60, 0.31, 0.06]
delta     = [0.064, 0.077, 0.086, np.nan]   # RAGE: not detected

# ---------- 配色（对齐原图系统） ----------
HAV_COLOR  = '#1F77B4'   # REC_COLORS['HAVCR2']
GRAY       = '#9E9E9E'
GRAY_EDGE  = '#6E6E6E'
ACCENT     = '#7BAE7F'   # SUBTYPE_COLORS['LA_TAM']，象限注记用

colors   = [HAV_COLOR, GRAY, GRAY, GRAY]
edgecols = [HAV_COLOR, GRAY_EDGE, GRAY_EDGE, GRAY_EDGE]

ratio = np.array(expr_LA) / np.array(expr_Mono)
size  = 150 * (ratio / ratio.max()) + 40

fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(7.2, 3.4),
    gridspec_kw={'width_ratios': [1.5, 1]})

# ---------- Panel A ----------
axA.axhline(0, color='#BBBBBB', lw=0.8, ls='--', zorder=1)
axA.axvline(0, color='#BBBBBB', lw=0.8, ls='--', zorder=1)

# 手工标签位置（修复重叠：TLR2 标上方、RAGE 标下方、TLR4 标右上）
offsets = {
    'HAVCR2':   (0.018,  0.000, 'left'),
    'TLR4':     (0.015,  0.022, 'left'),
    'TLR2':     (0.000,  0.042, 'center'),
    'RAGE/AGER':(0.000, -0.045, 'center'),
}
for i, r in enumerate(recs):
    axA.scatter(rho_pt[i], rho_mafb[i], s=size[i],
                c=colors[i], edgecolors=edgecols[i], linewidth=1.2,
                alpha=0.95 if i == 0 else 0.6,
                zorder=4 if i == 0 else 3)
    dx, dy, ha = offsets[r]
    axA.text(rho_pt[i] + dx, rho_mafb[i] + dy, r,
             fontsize=9, ha=ha, va='center',
             fontweight='bold' if i == 0 else 'normal',
             color=colors[i])

axA.text(0.10, 0.30, 'LA_TAM-enriched\n& MAFB-coupled',
         fontsize=7.5, style='italic', color=ACCENT,
         ha='right', va='center')

axA.set_xlabel('Expression vs pseudotime\n(Spearman $\\rho$)', fontsize=9)
axA.set_ylabel('Expression vs MAFB(+) regulon AUC\n(Spearman $\\rho$)', fontsize=9)
axA.set_xlim(-0.24, 0.34)
axA.set_ylim(-0.10, 0.34)
axA.text(-0.03, 1.02, 'a', transform=axA.transAxes,
         fontsize=12, fontweight='bold', va='bottom')

# ★图例移到右上（左下是数据区），并缩小
h1 = axA.scatter([], [], s=size[0], c='grey', alpha=0.5,
                 edgecolors='k', lw=0.8,
                 label=f'LA_TAM/Mono expression ratio\n({ratio[0]:.1f}$\\times$, HAVCR2)')
h2 = axA.scatter([], [], s=size[3], c='grey', alpha=0.5,
                 edgecolors='k', lw=0.8,
                 label=f'({ratio[3]:.1f}$\\times$, RAGE/AGER)')
axA.legend(handles=[h1, h2], loc='lower right', fontsize=6.5,
           frameon=True, framealpha=0.9, edgecolor='#CCCCCC',
           handletextpad=0.8, borderpad=0.4)

# ---------- Panel B ----------
labels_b = ['HAVCR2', 'TLR2', 'TLR4', 'RAGE/AGER\n(n.d.)']
vals_b   = [0.064, 0.077, 0.086, 0.0]
cols_b   = [HAV_COLOR, GRAY, GRAY, '#DDDDDD']
ypos     = np.arange(4)[::-1]

bars = axB.barh(ypos, vals_b, color=cols_b, edgecolor='white', height=0.55)
bars[3].set_hatch('///')
bars[3].set_edgecolor('#AAAAAA')

for y, v, r in zip(ypos, vals_b, recs):
    if not np.isnan(v):
        axB.text(v + 0.002, y, f'{v:.3f}', va='center', fontsize=8)
    else:
        axB.text(0.004, y, 'not detected among\nprioritized pairs',
                 va='center', fontsize=7, style='italic', color='#7F7F7F')

axB.set_yticks(ypos)
axB.set_yticklabels(labels_b, fontsize=9)
axB.set_xlabel('Differential engagement $\\Delta$\n(malignant vs normal B, Table S2)', fontsize=9)
axB.set_xlim(0, 0.105)
axB.text(-0.20, 1.02, 'b', transform=axB.transAxes,
         fontsize=12, fontweight='bold', va='bottom')

for ax in (axA, axB):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=0.8, length=3)

fig.suptitle('Candidate HMGB1 receptor prioritization across '
             'expression, trajectory, and regulon-coupling dimensions',
             fontsize=11, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig('FigS3_receptor_prioritization.png', dpi=300, bbox_inches='tight')
plt.savefig('FigS3_receptor_prioritization.pdf',       bbox_inches='tight')
plt.savefig('FigS3_receptor_prioritization.svg',       bbox_inches='tight')
plt.show()
print('输出: FigS3_receptor_prioritization.png / .pdf / .svg')
