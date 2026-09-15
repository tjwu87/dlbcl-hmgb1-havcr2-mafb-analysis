# =============================================================================
#  pySCENIC 分析 — 论文版（7 张图，尺寸 ×0.6，字号放大，统一配色）
# =============================================================================

import os, re, warnings
import numpy as np
import pandas as pd
from scipy.stats import zscore, spearmanr
from scipy.ndimage import gaussian_filter1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
import seaborn as sns
warnings.filterwarnings('ignore')

# =============================================================================
#  路径配置
# =============================================================================
BASE_DIR   = r'D:\bulk-download'
SCENIC_DIR = os.path.join(BASE_DIR, 'GSE182434', 'scenic')
REPO_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # 仓库根
OUT_DIR    = os.environ.get(
    'FIG3D_OUT',
    os.path.join(REPO_DIR, 'results', 'uav_panels', 'suoxiao_rss'),
)
DPI        = 300
os.makedirs(OUT_DIR, exist_ok=True)

CSV = {
    'auc':       os.path.join(SCENIC_DIR, 'data_AUC_matrix.csv'),
    'meta':      os.path.join(SCENIC_DIR, 'data_cell_metadata.csv'),
    'rss':       os.path.join(SCENIC_DIR, 'data_RSS_scores.csv'),
    'adj':       os.path.join(SCENIC_DIR, 'data_GRNBoost2_adjacencies.csv'),
    'reg':       os.path.join(SCENIC_DIR, 'data_RcisTarget_regulon_targets.csv'),
    'core':      os.path.join(SCENIC_DIR, 'data_core_TF_AUC_and_receptor_expression.csv'),
    'int':       os.path.join(SCENIC_DIR, 'data_TF_receptor_intersection.csv'),
    'spear':     os.path.join(SCENIC_DIR, 'data_pseudotime_spearman_correlations.csv'),
    'top_rss':   os.path.join(SCENIC_DIR, 'data_top30_RSS_per_subtype.csv'),
    'enr_gobp':  os.path.join(SCENIC_DIR, 'MAFB_gseapy_GO_BP_enrichment.csv'),
    'enr_kegg':  os.path.join(SCENIC_DIR, 'MAFB_gseapy_KEGG_enrichment.csv'),
    'enr_react': os.path.join(SCENIC_DIR, 'MAFB_gseapy_Reactome_enrichment.csv'),
}

def check_csv(key):
    path = CSV[key]
    if not os.path.exists(path):
        print(f'  [MISSING] {os.path.basename(path)}')
        return False
    return True

def need_redraw(fname):
    p = os.path.join(OUT_DIR, f'{fname}.png')
    if os.path.exists(p):
        print(f'  SKIP {fname} (already exists)')
        return False
    return True

def save_fig(fig, fname):
    for fmt in ['png', 'svg']:
        fig.savefig(os.path.join(OUT_DIR, f'{fname}.{fmt}'),
                    dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  [OK] {fname} saved')

# =============================================================================
#  ★ 统一配色（替换为 SUBTYPE_PALETTE 风格）
# =============================================================================
SUBTYPE_COLORS = {
    'Mono':    '#3498DB',   # ← 统一风格
    'IFN_TAM': '#E67E22',
    'LA_TAM':  '#2ECC71',
}
AXIS_COLORS = {
    'LA_TAM':  '#2ECC71',   # ← 统一风格
    'IFN_TAM': '#E67E22',
}
TF_COLORS = {
    'ETV5': '#D6604D',  'MITF': '#27AE60', 'MAFB': '#1E8449',
    'XBP1': '#A06080', 'FOSB': '#E67E22', 'KLF2': '#D35400', 'KLF4': '#E59866',
    'IRF1': '#CA6F1E', 'STAT1': '#3498DB',
}
REC_COLORS = {
    'AXL': '#2CA02C', 'HAVCR2': '#1F77B4',
    'CD74': '#2166AC', 'HLA-DPA1': '#4393C3',
}
path_subtypes = ['Mono', 'IFN_TAM', 'LA_TAM']

TOP5 = {
    'ETV5':  ['MYNN', 'UBE2O', 'GTF2H4', 'ZFAND2A', 'C12orf4'],
    'MITF':  ['THAP5', 'STMN1', 'DAB2', 'CCDC93', 'MCFD2'],
    'MAFB':  ['CEP192', 'CHCHD10', 'NAGK', 'PPRC1', 'BASP1'],
    'XBP1':  ['BLOC1S6', 'PDSS2', 'RBM12B', 'SSR3', 'CDC123'],
    'STAT1': ['UBD', 'GBP1', 'HLA-DQA2', 'CD86', 'IL32'],
    'IRF1':  ['CXCL10', 'GRAMD1A', 'PTPN2', 'IKBKE', 'SLC31A2'],
    'KLF4':  ['FOSB*', 'DUSP1', 'SRSF5', 'IER2', 'PHF19'],
    'KLF2':  ['YTHDC1', 'ZNF644', 'CTNNB1', 'CAMKK2', 'TMEM175'],
    'FOSB':  ['VCAN', 'CCNL1', 'JUND', 'JUN', 'UBXN2A'],
}
LA_TFS  = ['ETV5', 'MITF', 'MAFB', 'XBP1']
IFN_TFS = ['IRF1', 'KLF4', 'KLF2', 'FOSB']
SHARED  = ['STAT1']

print('Loading RSS CSV...')
rss_df = pd.read_csv(CSV['rss'], index_col=0)
print('rss_df:', rss_df.shape)

if True:
    print('\n[Fig 2] RSS scatter...')
    if rss_df is None:
        print('  SKIP: missing RSS CSV')
    else:
        try:
            from adjustText import adjust_text
            HAS_AT = True
        except ImportError:
            HAS_AT = False

        rss_la      = rss_df['LA_TAM'].values
        rss_ifn     = rss_df['IFN_TAM'].values
        reg_names   = rss_df.index.tolist()
        clean_names = [r.replace('(+)', '').replace('(-)', '') for r in reg_names]

        la_core_set  = {'ETV5',  'MITF', 'MAFB', 'XBP1'}
        ifn_core_set = {'FOSB', 'KLF2', 'KLF4', 'IRF1', 'STAT1'}
        threshold    = 0.33

        la_specific  = (rss_la > threshold) & (rss_la > rss_ifn + 0.02)
        ifn_specific = (rss_ifn > threshold) & (rss_ifn > rss_la + 0.02)
        shared       = ((rss_la > threshold) & (rss_ifn > threshold) &
                        ~la_specific & ~ifn_specific)
        other        = ~(la_specific | ifn_specific | shared)

        # ★ 尺寸 ×0.6；点颜色换统一风格
        fig, ax = plt.subplots(figsize=(5.4, 4.8))
        fig.patch.set_facecolor('white'); ax.set_facecolor('white')

        ax.scatter(rss_ifn[other], rss_la[other],
                   c='#CCCCCC', s=22, alpha=0.5, linewidths=0, zorder=2)
        for mask, color, label in [
            (la_specific,  AXIS_COLORS['LA_TAM'],  'LA_TAM specific'),   # ★
            (ifn_specific, AXIS_COLORS['IFN_TAM'], 'IFN_TAM specific'),  # ★
            (shared,       '#8B5CF6', 'Shared'),
        ]:
            ax.scatter(rss_ifn[mask], rss_la[mask], c=color, s=55, alpha=0.85,
                       linewidths=0.5, edgecolors='white', zorder=3, label=label)

        for i, (name, la, ifn) in enumerate(zip(clean_names, rss_la, rss_ifn)):
            if name in la_core_set:
                ax.scatter(ifn, la, c=AXIS_COLORS['LA_TAM'], s=160,
                           linewidths=1.5, edgecolors='#333', zorder=5, marker='*')
            elif name in ifn_core_set:
                ax.scatter(ifn, la, c=AXIS_COLORS['IFN_TAM'], s=160,
                           linewidths=1.5, edgecolors='#333', zorder=5, marker='*')

        lim_max = float(max(rss_la.max(), rss_ifn.max())) + 0.04
        ax.plot([0, lim_max], [0, lim_max], color='#AAAAAA', lw=0.8, ls='--', zorder=1)
        ax.axhline(threshold, color=AXIS_COLORS['LA_TAM'],  lw=0.8, ls=':', alpha=0.7)
        ax.axvline(threshold, color=AXIS_COLORS['IFN_TAM'], lw=0.8, ls=':', alpha=0.7)

        label_mask = la_specific | ifn_specific | shared
        texts = []
        for i, (name, la, ifn, lbl) in enumerate(
                zip(clean_names, rss_la, rss_ifn, label_mask)):
            if lbl or name in la_core_set or name in ifn_core_set:
                is_core = name in la_core_set or name in ifn_core_set
                color = (AXIS_COLORS['LA_TAM']  if la_specific[i] or name in la_core_set
                         else AXIS_COLORS['IFN_TAM'] if ifn_specific[i] or name in ifn_core_set
                         else '#8B5CF6')
                t = ax.text(float(ifn), float(la), name,
                            fontsize=9, color=color,   # ★ 7.5→9
                            fontweight='bold' if is_core else 'normal', zorder=6)
                texts.append(t)

        if HAS_AT:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.5),
                        expand=(1.2, 1.5), force_text=(0.3, 0.5))

        n_la  = int(la_specific.sum())
        n_ifn = int(ifn_specific.sum())
        n_sh  = int(shared.sum())
        ax.text(0.97, 0.03,
                f'LA_TAM specific: {n_la}\nIFN_TAM specific: {n_ifn}\nShared: {n_sh}',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=10,   # ★ 9→10
                bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#CCCCCC', alpha=0.9))

        legend_elements = [
            mpatches.Patch(color=AXIS_COLORS['LA_TAM'],  label='LA_TAM specific'),
            mpatches.Patch(color=AXIS_COLORS['IFN_TAM'], label='IFN_TAM specific'),
            mpatches.Patch(color='#8B5CF6', label='Shared'),
            mpatches.Patch(color='#CCCCCC', label='Other'),
            Line2D([0], [0], marker='*', color='w', markerfacecolor='#333',
                   markersize=10, label='Core axis TF'),
        ]
        ax.legend(handles=legend_elements, fontsize=8.5, frameon=True,
                  loc='lower left', framealpha=0.9, edgecolor='#CCCCCC',
                  handlelength=1.2, handletextpad=0.4, labelspacing=0.3,
                  borderpad=0.4)  # ★ 紧凑小图例
        ax.set_xlabel('RSS — IFN_TAM', fontsize=12)   # ★ 11→12
        ax.set_ylabel('RSS — LA_TAM',  fontsize=12)
        ax.tick_params(labelsize=10)                   # ★ 新增
        ax.set_title(
            f'Regulon Specificity Score: LA_TAM vs IFN_TAM',
            fontsize=12, fontweight='bold')   # ★ 标题缩为一行
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_xlim(-0.01, lim_max); ax.set_ylim(-0.01, lim_max)
        plt.tight_layout()
        save_fig(fig, 'Fig2_RSS_scatter_paper')

# =============================================================================
#  Fig 3：TF-Target Network
#  原尺寸 (13, 11) → (7.8, 6.6)
print('DONE')
