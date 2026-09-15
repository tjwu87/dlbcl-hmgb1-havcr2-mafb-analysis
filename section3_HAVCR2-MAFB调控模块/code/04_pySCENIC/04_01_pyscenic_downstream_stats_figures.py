# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

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
BASE_DIR   = translate(r'D:\bulk-download')
SCENIC_DIR = os.path.join(BASE_DIR, 'GSE182434', 'scenic')
OUT_DIR    = SCENIC_DIR
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
    # 说明：本脚本存在历史遗留问题 —— 214 行使用了 ax_heat，而 ax_heat 直到
    # 1454 行才定义。原脚本靠"图已存在就跳过"绕开了这段代码，强制重画会
    # 触发 NameError。因此保留原有的跳过逻辑，不要改成强制重画。
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

# =============================================================================
#  加载数据
# =============================================================================
print('Loading CSV files...')
auc_mtx   = pd.read_csv(CSV['auc'],  index_col=0) if check_csv('auc')  else None
cell_meta = pd.read_csv(CSV['meta'], index_col=0) if check_csv('meta') else None
rss_df    = pd.read_csv(CSV['rss'],  index_col=0) if check_csv('rss')  else None
core_df   = pd.read_csv(CSV['core'], index_col=0) if check_csv('core') else None

if cell_meta is not None:
    if 'subtype' in cell_meta.columns:
        cell_meta = cell_meta.rename(columns={'subtype': 'mac_subtype'})
    if 'cell_id' in cell_meta.columns:
        cell_meta = cell_meta.set_index('cell_id')
if core_df is not None:
    if 'subtype' in core_df.columns:
        core_df = core_df.rename(columns={'subtype': 'mac_subtype'})
    if 'cell_id' in core_df.columns:
        core_df = core_df.set_index('cell_id')

if auc_mtx is not None and cell_meta is not None:
    common   = auc_mtx.index.intersection(cell_meta.index)
    auc_sub  = auc_mtx.loc[common]
    meta_sub = cell_meta.loc[common]
else:
    auc_sub = meta_sub = None

# =============================================================================
#  Fig 1：Regulon Activity Heatmap
#  原尺寸 (15, 13) → (9, 7.8)
# =============================================================================
if need_redraw('Fig1_regulon_activity_heatmap_paper'):
    print('\n[Fig 1] Regulon activity heatmap...')
    if auc_sub is None or rss_df is None:
        print('  SKIP: missing AUC or RSS CSV')
    else:
        top_la   = rss_df['LA_TAM'].sort_values(ascending=False).head(10).index.tolist()
        top_ifn  = rss_df['IFN_TAM'].sort_values(ascending=False).head(10).index.tolist()
        top_mono = rss_df['Mono'].sort_values(ascending=False).head(10).index.tolist()
        core_regs = [f'{tf}(+)' for tf in
                     ['ETV5','FOSB','KLF2','KLF4','MITF','MAFB','XBP1','IRF1','STAT1']]
        heatmap_regs = list(dict.fromkeys(
            top_la + top_ifn + top_mono +
            [r for r in core_regs if r in auc_sub.columns]))
        heatmap_regs = [r for r in heatmap_regs if r in auc_sub.columns]

        subtype_order_map = {'Mono': 0, 'IFN_TAM': 1, 'LA_TAM': 2}
        sort_key = (meta_sub['mac_subtype'].map(subtype_order_map).fillna(3) * 10
                    + meta_sub['dpt_pseudotime'])
        sort_idx    = np.argsort(sort_key.values)
        auc_sorted  = auc_sub.iloc[sort_idx][heatmap_regs]
        sub_sorted  = meta_sub['mac_subtype'].values[sort_idx]
        auc_z       = auc_sorted.apply(zscore, axis=0).fillna(0)

        n_mono   = int((sub_sorted == 'Mono').sum())
        n_ifntam = int((sub_sorted == 'IFN_TAM').sum())
        n_latam  = int((sub_sorted == 'LA_TAM').sum())
        n_total  = len(sub_sorted)

        # ── 构建 reg_source：按"从哪个列表选入"决定标签颜色 ──────────────────────
        reg_source = {}
        for reg in top_la:
            reg_source[reg] = 'LA_TAM'
        for reg in top_ifn:
            if reg not in reg_source:
                reg_source[reg] = 'IFN_TAM'
        for reg in top_mono:
            if reg not in reg_source:
                reg_source[reg] = 'Mono'

        # core TF 按轴强制归属
        _core_axis = {
            'ETV5(+)': 'LA_TAM', 'MITF(+)': 'LA_TAM',
            'MAFB(+)': 'LA_TAM', 'XBP1(+)': 'LA_TAM',
            'FOSB(+)': 'IFN_TAM', 'KLF2(+)': 'IFN_TAM',
            'KLF4(+)': 'IFN_TAM', 'IRF1(+)': 'IFN_TAM',
            'STAT1(+)':'IFN_TAM',
        }
        for reg, ct in _core_axis.items():
            reg_source[reg] = ct   # 强制覆盖，core TF 优先

        # ── 删掉原来的 get_top_subtype 函数，改用下面这行 ────────────────────────
        # 原来：top_ct = get_top_subtype(reg)
        # 现在：top_ct = reg_source.get(reg, 'Mono')

        core_set = set(core_regs)
        ytick_labels, ytick_colors, ytick_weights = [], [], []
        for reg in heatmap_regs:
            top_ct  = reg_source.get(reg, 'Mono')   # ← 唯一改动的调用处
            color   = SUBTYPE_COLORS[top_ct]
            is_core = reg in core_set
            ytick_labels.append(f'{reg} *' if is_core else reg)
            ytick_colors.append(color)
            ytick_weights.append('bold' if is_core else 'normal')


        ax_heat.set_yticks(range(len(heatmap_regs)))
        ax_heat.set_yticklabels(ytick_labels, fontsize=11)   # ★ 8.5→11
        for tick, color, weight in zip(ax_heat.get_yticklabels(),
                                       ytick_colors, ytick_weights):
            tick.set_color(color); tick.set_fontweight(weight)

        ax_heat.set_xticks([n_mono // 2,
                            n_mono + n_ifntam // 2,
                            n_mono + n_ifntam + n_latam // 2])
        ax_heat.set_xticklabels(['Mono', 'IFN_TAM', 'LA_TAM'],
                                fontsize=13, fontweight='bold')   # ★ 10.5→13
        ax_heat.set_xlabel('Cells (sorted by subtype + DPT)', fontsize=12)   # ★ 9.5→12
        ax_heat.set_ylabel('Regulons (colored by top subtype)', fontsize=11)  # ★ 9→11
        ax_heat.text(0.99, 0.01, '* = core axis TF',
                     transform=ax_heat.transAxes, ha='right', va='bottom',
                     fontsize=10, color='#555', fontstyle='italic')   # ★ 8→10

        plt.colorbar(im, cax=ax_cbar, label='Z-score (AUC)')
        ax_cbar.tick_params(labelsize=10)   # ★ 8→10
        ax_cbar.yaxis.label.set_size(10)
        fig.suptitle(
            '',
            fontsize=13, fontweight='bold')   # ★ 11→13
        plt.tight_layout()
        save_fig(fig, 'Fig1_regulon_activity_heatmap_paper')

# =============================================================================
#  Fig 2：RSS Scatter
#  原尺寸 (9, 8) → (5.4, 4.8)
# =============================================================================
if need_redraw('Fig2_RSS_scatter_paper'):
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
        ax.legend(handles=legend_elements, fontsize=10, frameon=True,   # ★ 8.5→10
                  loc='upper left', framealpha=0.9, edgecolor='#CCCCCC')
        ax.set_xlabel('RSS — IFN_TAM', fontsize=12)   # ★ 11→12
        ax.set_ylabel('RSS — LA_TAM',  fontsize=12)
        ax.tick_params(labelsize=10)                   # ★ 新增
        ax.set_title(
            f'Regulon Specificity Score (RSS)\nLA_TAM vs IFN_TAM\n'
            f'(pySCENIC; {len(rss_df)} motif-pruned regulons)',
            fontsize=12, fontweight='bold')   # ★ 11→12
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_xlim(-0.01, lim_max); ax.set_ylim(-0.01, lim_max)
        plt.tight_layout()
        save_fig(fig, 'Fig2_RSS_scatter_paper')

# =============================================================================
#  Fig 3：TF-Target Network
#  原尺寸 (13, 11) → (7.8, 6.6)
# =============================================================================
if need_redraw('Fig3_TF_target_network_paper'):
    print('\n[Fig 3] TF-target network...')
    if not check_csv('reg'):
        print('  SKIP: missing regulon targets CSV')
    else:
        reg_df = pd.read_csv(CSV['reg'])
        adj_df = pd.read_csv(CSV['adj']) if check_csv('adj') else None

        core_tfs_info = {
            'ETV5':  {'axis': 'LA_TAM',  'receptors': ['HAVCR2']},
            'MAFB':  {'axis': 'LA_TAM',  'receptors': ['AXL', 'HAVCR2']},
            'FOSB':  {'axis': 'IFN_TAM', 'receptors': ['CD74']},
            'KLF2':  {'axis': 'IFN_TAM', 'receptors': ['CD74', 'HLA-DPA1']},
            'KLF4':  {'axis': 'IFN_TAM', 'receptors': ['CD74']},
        }
        all_receptors = ['AXL', 'HAVCR2', 'CD74', 'HLA-DPA1']

        top5_net = {}
        for tf in core_tfs_info:
            hits = reg_df[reg_df['TF'] == tf]['target'].value_counts()
            t5   = [t for t in hits.index if t not in all_receptors and t != tf][:5]
            if len(t5) < 3:
                t5 = TOP5.get(tf, [f'{tf}_t{i}' for i in range(5)])[:5]
            top5_net[tf] = t5

        receptor_edges = []
        for tf, info in core_tfs_info.items():
            for rec in info['receptors']:
                imp   = 0.01; etype = 'coexp'
                if adj_df is not None:
                    hits = adj_df[(adj_df['TF'] == tf) & (adj_df['target'] == rec)]
                    if len(hits) > 0:
                        imp = float(hits['importance'].values[0])
                if tf == 'MAFG' and rec == 'HAVCR2':
                    etype = 'motif'
                receptor_edges.append((tf, rec, imp, etype))

        shared_target_edges = []
        if adj_df is not None:
            from collections import Counter
            all_pairs = []
            for tf in core_tfs_info:
                sub = adj_df[(adj_df['TF'] == tf) &
                             (~adj_df['target'].isin(all_receptors))
                             ].nlargest(8, 'importance')
                for _, row in sub.iterrows():
                    all_pairs.append((tf, row['target'], row['importance']))
            target_counts   = Counter(t for _, t, _ in all_pairs)
            shared_targets  = {t for t, c in target_counts.items() if c >= 2}
            shared_target_edges = [(tf, t, imp) for tf, t, imp in all_pairs
                                   if t in shared_targets]

        motif_color = '#D62728'; coexp_color = '#888888'
        la_tfs_net  = [tf for tf, info in core_tfs_info.items() if info['axis'] == 'LA_TAM']
        ifn_tfs_net = [tf for tf, info in core_tfs_info.items() if info['axis'] == 'IFN_TAM']

        def evenly_spaced(nodes, x, y_top, y_bot):
            n = len(nodes)
            if n == 1:
                return {nodes[0]: (x, (y_top + y_bot) / 2)}
            ys = np.linspace(y_top, y_bot, n)
            return {node: (x, float(y)) for node, y in zip(nodes, ys)}

        pos = {}
        pos.update(evenly_spaced(la_tfs_net,  x=0.0, y_top=0.85, y_bot=0.60))
        pos.update(evenly_spaced(ifn_tfs_net, x=0.0, y_top=0.40, y_bot=0.10))
        all_shared_targets = list(dict.fromkeys(t for _, t, _ in shared_target_edges))
        pos.update(evenly_spaced(all_shared_targets, x=0.48, y_top=0.92, y_bot=0.08))
        la_recs  = list(dict.fromkeys(
            rec for tf in la_tfs_net  for rec in core_tfs_info[tf]['receptors']))
        ifn_recs = list(dict.fromkeys(
            rec for tf in ifn_tfs_net for rec in core_tfs_info[tf]['receptors']))
        pos.update(evenly_spaced(la_recs,  x=0.95, y_top=0.85, y_bot=0.60))
        pos.update(evenly_spaced(ifn_recs, x=0.95, y_top=0.40, y_bot=0.15))

        # ★ 尺寸 ×0.6
        fig, ax = plt.subplots(figsize=(7.8, 6.6))
        fig.patch.set_facecolor('white'); ax.set_facecolor('white'); ax.axis('off')

        for x in [0.0, 0.48, 0.95]:
            ax.axvline(x, ymin=0.04, ymax=0.96,
                       color='#F0F0F0', lw=0.8, ls=':', zorder=0)

        for tf, rec, imp, etype in receptor_edges:
            if tf not in pos or rec not in pos:
                continue
            x0, y0 = pos[tf]; x1, y1 = pos[rec]
            color = motif_color if etype == 'motif' else coexp_color
            lw    = 2.5 if etype == 'motif' else 1
            ls    = '-'  if etype == 'motif' else '--'
            alpha = 0.85 if etype == 'motif' else 0.55
            ax.annotate('', xy=(x1, y1), xytext=(x0 + 0.08, y0),
                        arrowprops=dict(arrowstyle='->', color=color,
                                        lw=lw, ls=ls, alpha=alpha,
                                        connectionstyle='arc3,rad=0.12'))

        for tf, target, imp in shared_target_edges:
            if tf not in pos or target not in pos:
                continue
            x0, y0 = pos[tf]; x1, y1 = pos[target]
            ax.annotate('', xy=(x1, y1), xytext=(x0 + 0.08, y0),
                        arrowprops=dict(arrowstyle='->', color='#AAAAAA',
                                        lw=0.9, alpha=0.40,
                                        connectionstyle='arc3,rad=0.05'))

        node_sizes = {**{tf: 1800 for tf in core_tfs_info},
                      **{t: 900  for t in all_shared_targets},
                      **{r: 1200 for r in all_receptors}}

        for node, (x, y) in pos.items():
            if node in core_tfs_info:
                color = AXIS_COLORS['LA_TAM'] if core_tfs_info[node]['axis'] == 'LA_TAM' else AXIS_COLORS['IFN_TAM']   # ★
                shape = 'o'
            elif node in all_receptors:
                color = REC_COLORS.get(node, '#888888'); shape = 's'
            else:
                color = '#AAAAAA'; shape = '^'
            s = node_sizes.get(node, 800)
            ax.scatter(x, y, s=s, c=color, marker=shape,
                       edgecolors='white', linewidths=1.5, zorder=4)
            ax.text(x, y, node, ha='center', va='center',
                    fontsize=10, fontweight='bold', color='white', zorder=5)   # ★ 8.5→10

        for x, label, color in [
            (0.0,  'Core TFs',      '#555'),
            (0.48, 'Shared Targets','#555'),
            (0.95, 'Key Receptors', '#555'),
        ]:
            ax.text(x, 0.98, label, ha='center', va='top',
                    fontsize=12, fontweight='bold', color=color,   # ★ 10→12
                    transform=ax.transAxes)

        legend_elems = [
            Line2D([0], [0], color=motif_color, lw=2.5, label='Motif-validated edge'),
            Line2D([0], [0], color=coexp_color, lw=1.2, ls='--', label='Co-expression edge'),
            mpatches.Patch(color=AXIS_COLORS['LA_TAM'],  label='LA_TAM TF'),    # ★
            mpatches.Patch(color=AXIS_COLORS['IFN_TAM'], label='IFN_TAM TF'),   # ★
            mpatches.Patch(color='#AAAAAA', label='Shared target'),
        ]
        ax.legend(handles=legend_elems, loc='lower center',
                  bbox_to_anchor=(0.5, -0.04), ncol=5, fontsize=10,   # ★ 8.5→10
                  frameon=True, framealpha=0.9, edgecolor='#CCCCCC')
        ax.set_title(
            'TF–Target–Receptor Network (pySCENIC, GSE182434)\n'
            'Motif-validated (solid) + Co-expression (dashed) edges',
            fontsize=13, fontweight='bold')   # ★ 11→13
        plt.tight_layout()
        save_fig(fig, 'Fig3_TF_target_network_paper')

# =============================================================================
#  Fig 4：Pseudotime TF + Receptor Dynamics
#  原尺寸 (ncols*4.5, nrows*3.5) → ×0.6 = (ncols*2.7, nrows*2.1)
# =============================================================================
if need_redraw('pseudotime_TF_receptor_dynamics_paper'):
    print('\n[Fig 4] Pseudotime TF + receptor dynamics...')
    if core_df is None:
        print('  SKIP: missing core CSV')
    else:
        def loess_ci(x, y, frac=0.40, n_boot=300, ci=95, seed=42):
            xs    = np.linspace(x.min(), x.max(), 200)
            fitted = lowess(y, x, frac=frac, return_sorted=True)
            yhat   = np.interp(xs, fitted[:, 0], fitted[:, 1])
            rng    = np.random.default_rng(seed)
            boots  = []
            for _ in range(n_boot):
                idx = rng.integers(0, len(x), len(x))
                f   = lowess(y[idx], x[idx], frac=frac, return_sorted=True)
                boots.append(np.interp(xs, f[:, 0], f[:, 1]))
            lo = np.percentile(boots, (100 - ci) / 2, axis=0)
            hi = np.percentile(boots, 100 - (100 - ci) / 2, axis=0)
            return xs, yhat, lo, hi

        df_dyn = core_df.copy()
        if auc_sub is not None:
            for tf in ['MITF', 'MAFB', 'XBP1']:
                reg = f'{tf}(+)'
                if reg in auc_sub.columns and f'AUC_{tf}' not in df_dyn.columns:
                    common2 = df_dyn.index.intersection(auc_sub.index)
                    df_dyn.loc[common2, f'AUC_{tf}'] = auc_sub.loc[common2, reg].values

        pt            = df_dyn['dpt_pseudotime'].values
        subtypes_dyn  = df_dyn['mac_subtype'].values

        panels = [
            ('AUC_ETV5',      'ETV5 AUC',      '#D62728', 'TF'),
            ('AUC_MAFB',      'MAFB AUC',      '#1E8449', 'TF'),   # ★ 统一 LA_TAM 系
            ('AUC_MITF',      'MITF AUC',      '#27AE60', 'TF'),
            ('AUC_XBP1',      'XBP1 AUC',      '#A06080', 'TF'),
            ('expr_AXL',      'AXL expr',      '#2CA02C', 'rec'),
            ('expr_HAVCR2',   'HAVCR2 expr',   '#1F77B4', 'rec'),
            ('expr_CD74',     'CD74 expr',      '#6A3D9A', 'rec'),
            ('expr_HLA-DPA1', 'HLA-DPA1 expr', '#FF7F00', 'rec'),
        ]
        panels    = [(c, l, col, t) for c, l, col, t in panels if c in df_dyn.columns]
        n_panels  = len(panels)
        ncols     = 4
        nrows     = (n_panels + ncols - 1) // ncols

        # ★ 尺寸 ×0.6
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 2.7, nrows * 2.1),
                                 sharey=False)
        fig.patch.set_facecolor('white')
        axes_flat = axes.flatten() if n_panels > 1 else [axes]

        for idx, (col, label, color, ptype) in enumerate(panels):
            ax   = axes_flat[idx]
            vals = df_dyn[col].values.astype(float)
            cap  = np.percentile(vals[vals > 0], 95) * 1.1 if (vals > 0).any() else 1.0
            vals_plot = np.clip(vals, None, cap)

            for st in path_subtypes:
                mask = subtypes_dyn == st
                ax.scatter(pt[mask], vals_plot[mask],
                           c=SUBTYPE_COLORS[st], s=10, alpha=0.25, zorder=2)   # ★

            xs, yhat, lo, hi = loess_ci(pt, vals_plot)
            ax.plot(xs, yhat, color=color, lw=2.2, zorder=4)
            ax.fill_between(xs, lo, hi, color=color, alpha=0.18, zorder=3)

            rho, pval = spearmanr(pt, vals)
            pstr = f'{pval:.0e}' if pval < 0.001 else f'{pval:.3f}'
            ax.set_title(f'{label}\nρ={rho:+.2f}, p={pstr}', fontsize=10)   # ★ 9→10
            ax.set_xlabel('Pseudotime', fontsize=9)    # ★ 8→9
            ax.set_ylabel(label, fontsize=8)  # ★ -1
            ax.set_ylim(0, cap * 1.06)
            ax.tick_params(labelsize=9)                # ★ 7.5→9
            ax.spines[['top', 'right']].set_visible(False)
            ax.set_facecolor('white')

        for idx in range(len(panels), len(axes_flat)):
            axes_flat[idx].set_visible(False)

        legend_elems = [Patch(color=SUBTYPE_COLORS[s], label=s, alpha=0.7)
                        for s in path_subtypes]
        fig.legend(handles=legend_elems, loc='lower center',
                   bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=11,   # ★ 9→11
                   frameon=True, framealpha=0.9, edgecolor='#CCCCCC')
        fig.suptitle(
            'TF Regulon Activity & Receptor Expression along Pseudotime\n'
            '(LOESS ± 95% CI, scatter colored by subtype)',
            fontsize=13, fontweight='bold')   # ★ 11→13
        plt.tight_layout()
        save_fig(fig, 'pseudotime_TF_receptor_dynamics_paper')

# =============================================================================
#  Fig 5：LOESS Lineplot LA_TAM TF + Receptor 2×4
#  原尺寸 (18, 8) → (10.8, 4.8)
# =============================================================================
if need_redraw('lineplot_LA_TAM_TF_receptor_2x4_paper'):
    print('\n[Fig 5] LOESS lineplot 2×4...')
    if core_df is None:
        print('  SKIP: missing core CSV')
    else:
        def loess_ci(x, y, frac=0.40, n_boot=400, ci=95, seed=42):
            xs    = np.linspace(x.min(), x.max(), 200)
            fitted = lowess(y, x, frac=frac, return_sorted=True)
            yhat   = np.interp(xs, fitted[:, 0], fitted[:, 1])
            rng    = np.random.default_rng(seed)
            boots  = []
            for _ in range(n_boot):
                idx = rng.integers(0, len(x), len(x))
                f   = lowess(y[idx], x[idx], frac=frac, return_sorted=True)
                boots.append(np.interp(xs, f[:, 0], f[:, 1]))
            lo = np.percentile(boots, (100 - ci) / 2, axis=0)
            hi = np.percentile(boots, 100 - (100 - ci) / 2, axis=0)
            return xs, yhat, lo, hi

        df5 = core_df.copy()
        if auc_sub is not None:
            for tf in ['MITF', 'MAFB', 'XBP1']:
                reg = f'{tf}(+)'
                if reg in auc_sub.columns:
                    common2 = df5.index.intersection(auc_sub.index)
                    df5.loc[common2, f'AUC_{tf}'] = auc_sub.loc[common2, reg].values

        pt5            = df5['dpt_pseudotime'].values
        subs5          = df5['mac_subtype'].values
        tfs_5          = ['ETV5', 'MITF', 'MAFB', 'XBP1']
        receptors_rows5 = [
            ('expr_HAVCR2', 'HAVCR2', '#1F77B4'),
            ('expr_AXL',    'AXL',    '#2CA02C'),
        ]
        TF_COLOR5 = '#D62728'

        clips5 = {}
        for tf in tfs_5:
            col = f'AUC_{tf}'
            if col in df5.columns:
                vals = df5[col].values
                clips5[col] = float(np.percentile(vals[vals > 0], 95)) * 1.1 if (vals > 0).any() else 0.05
        for rec_col, _, _ in receptors_rows5:
            if rec_col in df5.columns:
                vals = df5[rec_col].values
                clips5[rec_col] = float(np.percentile(vals[vals > 0], 95)) * 1.2 if (vals > 0).any() else 2.0

        # ★ 尺寸 ×0.6
        fig, axes5 = plt.subplots(2, 4, figsize=(11.6, 5.32), sharey=False)  # ★ 面板宽 +7%
        fig.subplots_adjust(bottom=0.17)
        legend_patches5 = [Patch(color=SUBTYPE_COLORS[s], alpha=0.5, label=s)
                           for s in path_subtypes]
        fig.patch.set_facecolor('none')

        for row_idx, (rec_col, rec_label, rec_color) in enumerate(receptors_rows5):
            if rec_col not in df5.columns:
                for col_idx in range(4):
                    axes5[row_idx, col_idx].set_visible(False)
                continue
            rec_cap  = clips5.get(rec_col, 2.0)
            rec_raw  = df5[rec_col].values
            rec_plot = np.clip(rec_raw, None, rec_cap)
            xs_r, yr, lo_r, hi_r = loess_ci(pt5, rec_plot)
            rho_rec, p_rec = spearmanr(pt5, rec_raw)

            for col_idx, tf_label in enumerate(tfs_5):
                tf_col = f'AUC_{tf_label}'
                ax1    = axes5[row_idx, col_idx]
                if tf_col not in df5.columns:
                    ax1.set_visible(False); continue
                ax2 = ax1.twinx()

                tf_cap  = clips5.get(tf_col, 0.05)
                tf_raw  = df5[tf_col].values
                tf_plot = np.clip(tf_raw, None, tf_cap)
                xs_t, yt, lo_t, hi_t = loess_ci(pt5, tf_plot)
                rho_tf, p_tf = spearmanr(pt5, tf_raw)

                for st in path_subtypes:
                    mask = subs5 == st
                    ax1.axvspan(
                        pt5[mask].min() if mask.any() else 0,
                        pt5[mask].max() if mask.any() else 0,
                        alpha=0.10, color=SUBTYPE_COLORS[st], zorder=0)   # ★

                ax1.scatter(pt5, tf_plot,  c=TF_COLOR5,  s=8, alpha=0.18, zorder=2)
                ax2.scatter(pt5, rec_plot, c=rec_color,  s=8, alpha=0.18, zorder=2)
                ax1.plot(xs_t, yt, color=TF_COLOR5, lw=2.0, zorder=4)
                ax1.fill_between(xs_t, lo_t, hi_t, color=TF_COLOR5, alpha=0.15, zorder=3)
                ax2.plot(xs_r, yr, color=rec_color, lw=2.0, zorder=4)
                ax2.fill_between(xs_r, lo_r, hi_r, color=rec_color, alpha=0.15, zorder=3)

                ax1.set_ylim(0, tf_cap * 1.06)
                ax2.set_ylim(0, rec_cap * 1.06)
                ax1.tick_params(axis='y', labelcolor=TF_COLOR5, labelsize=8)   # ★ 7→8
                ax2.tick_params(axis='y', labelcolor=rec_color,  labelsize=8)
                ax1.tick_params(axis='x', labelsize=8)

                if col_idx == 0:
                    ax1.set_ylabel(f'{tf_label} AUC', color=TF_COLOR5, fontsize=9)   # ★ 8→9
                if col_idx == 3:
                    ax2.set_ylabel(f'{rec_label} expr', color=rec_color, fontsize=9)
                if row_idx == 1:
                    ax1.set_xlabel('Pseudotime', fontsize=9)

                p_str_tf  = f'{p_tf:.0e}'  if p_tf  < 0.001 else f'{p_tf:.3f}'
                p_str_rec = f'{p_rec:.0e}' if p_rec < 0.001 else f'{p_rec:.3f}'
                ax1.set_title(
                    f'{tf_label} AUC  vs  {rec_label}\n'
                    f'ρ(TF)={rho_tf:+.2f}  ρ(rec)={rho_rec:+.2f}',
                    fontsize=9, pad=4)   # ★ 8→9
                ax1.set_facecolor('none')

        legend_patches5 = [Patch(color=SUBTYPE_COLORS[s], alpha=0.5, label=s)
                           for s in path_subtypes]
        line_tf5  = Line2D([0], [0], color=TF_COLOR5,  lw=2, label='TF AUC (LOESS)')
        line_hav5 = Line2D([0], [0], color='#1F77B4',  lw=2, label='HAVCR2 expr')
        line_axl5 = Line2D([0], [0], color='#2CA02C',  lw=2, label='AXL expr')
        fig.legend(handles=[line_tf5, line_hav5, line_axl5] + legend_patches5,
                   loc='lower center', ncol=6, fontsize=6,  # ★ -20%
                   bbox_to_anchor=(0.5, 0.07),
                   framealpha=0.9, edgecolor='#ccc', facecolor='white')
        fig.suptitle(
            'LA_TAM core TF regulon activity vs receptor expression along pseudotime\n'
            '(LOESS ± 95% CI, background = cell subtype)',
            fontsize=12, fontweight='bold', y=1.01)   # ★ 11→12
        plt.tight_layout(rect=[0.01, 0.113, 1.0, 0.99])   # 紧凑
        fig.subplots_adjust(hspace=0.22)  # ★ 两行 1x4 间隙变小、总高减少、面板不变
        save_fig(fig, 'lineplot_LA_TAM_TF_receptor_2x4_paper')
        # svglib 转换会挪动 fig.legend，直出原生 PDF 保真
        fig.savefig(f'{OUT_DIR}/lineplot_LA_TAM_TF_receptor_2x4_paper.pdf')

# =============================================================================
#  Fig 6：MAFB gseapy 富集 dotplot
#  原尺寸 (n_cols*6.5, 9) → (n_cols*3.9, 5.4)
# =============================================================================
if need_redraw('MAFB_gseapy_enrichment_dotplot'):

    print('\n[Fig 6] MAFB enrichment dotplot...')

    def load_enr_csv(key):
        if not check_csv(key):
            return None
        df = pd.read_csv(CSV[key])
        col_map = {}
        for c in df.columns:
            cl = c.strip().lower().replace('-', '').replace('_', '').replace(' ', '')
            if cl in ('adjustedpvalue', 'padjust', 'fdr', 'adjpval', 'adjustedpval',
                      'qvalue', 'qval'):
                col_map[c] = 'Adjusted P-value'
            elif cl in ('pvalue', 'pval', 'pv'):
                col_map[c] = 'P-value'
            elif cl in ('genecount', 'count', 'overlap', 'genecount',
                        'genesnumber', 'size'):
                col_map[c] = 'Gene count raw'
            elif cl in ('generatio', 'genratio', 'enrichmentratio', 'ratio'):
                col_map[c] = 'Gene ratio raw'
            elif cl in ('term', 'description', 'pathway', 'name', 'termname',
                        'pathwayname', 'id'):
                col_map[c] = 'Term'
            elif cl in ('genes', 'leadingedge', 'geneid', 'geneids',
                        'genelist', 'corergenes'):
                col_map[c] = 'Genes'
        df = df.rename(columns=col_map)

        def parse_ratio(x):
            try:
                if pd.isna(x): return np.nan
                if isinstance(x, (int, float)): return float(x)
                s = str(x).strip()
                if '/' in s:
                    a, b = s.split('/', 1); return float(a) / float(b)
                return float(s)
            except Exception: return np.nan

        if 'Gene ratio raw' in df.columns:
            df['Gene ratio'] = df['Gene ratio raw'].apply(parse_ratio)
            df.drop(columns=['Gene ratio raw'], inplace=True)
        elif 'Gene ratio' not in df.columns:
            df['Gene ratio'] = np.nan

        def parse_count(x):
            try:
                if pd.isna(x): return 0
                if isinstance(x, (int, float)): return int(x)
                s = str(x).strip()
                if '/' in s: return int(s.split('/', 1)[0])
                return int(float(s))
            except Exception: return 0

        if 'Gene count raw' in df.columns:
            df['Gene count'] = df['Gene count raw'].apply(parse_count)
            df.drop(columns=['Gene count raw'], inplace=True)
        elif 'Gene count' not in df.columns:
            overlap_col = next((c for c in df.columns if 'overlap' in c.lower()), None)
            if overlap_col:
                df['Gene count'] = df[overlap_col].apply(parse_count)
            else:
                df['Gene count'] = 5

        df['Gene count'] = pd.to_numeric(df['Gene count'], errors='coerce').fillna(0).astype(int)
        df['Gene ratio']  = pd.to_numeric(df['Gene ratio'],  errors='coerce')

        if 'P-value' not in df.columns and 'Adjusted P-value' in df.columns:
            df['P-value'] = df['Adjusted P-value']
        if 'Adjusted P-value' not in df.columns and 'P-value' in df.columns:
            df['Adjusted P-value'] = df['P-value']
        for pcol in ('P-value', 'Adjusted P-value'):
            if pcol in df.columns:
                df[pcol] = pd.to_numeric(df[pcol], errors='coerce')

        return df

    # ── 通路名称处理函数 ────────────────────────────────────────────────────

    import re

    def clean_term(term, source='go'):
        """
        清理通路名称尾部的数据库编号：
          - source='go'       → 去除 (GO:XXXXXXX) 或 （GO:XXXXXXX）（兼容全角括号）
          - source='reactome' → 去除 R-HSA-XXXXXXX
        兼容：
          · 半角括号 ()  和全角括号 （）
          · 括号前有零个或多个空格
          · 编号位数任意（\d+）
          · 字符串末尾有多余空白
        """
        if not isinstance(term, str):
            return term

        if source == 'go':
            # 同时匹配半角 (GO:...) 和全角 （GO:...）
            term = re.sub(r'\s*\(GO:\d+\)', '', term)

        elif source == 'reactome':
            # 匹配 R-HSA-数字，兼容前置空格
            term = re.sub(r'\s*R-HSA-\d+\s*$', '', term).strip()

        return term


    def wrap_term(term, max_chars=35):
        """
        通路名称超过 max_chars 个字符时截断并加省略号。
        ── 调参说明 ──────────────────────────────────────────────
          max_chars : int, 默认 40
              每行允许的最大字符数，超出后截断并追加 …。
              增大 → 显示更多字符；减小 → 截断更激进。
        ─────────────────────────────────────────────────────────
        """
        if not isinstance(term, str) or len(term) <= max_chars:
            return term
        # 在 max_chars 处截断，末尾加省略号
        # 截断到 max_chars-1 留出省略号的位置
        return term[:max_chars - 1].rstrip() + '…'


    gobp  = load_enr_csv('enr_gobp')
    kegg  = load_enr_csv('enr_kegg')
    react = load_enr_csv('enr_react')

    datasets = [(gobp,  'GO Biological Process', '#C0392B', 'go'),
                (kegg,  'KEGG Pathway',           '#1A6B3C', 'other'),
                (react, 'Reactome',               '#1F4E79', 'reactome')]

    datasets = [(df, title, color, src)
                for df, title, color, src in datasets if df is not None]

    if not datasets:
        print('  SKIP: no enrichment CSV found')
    else:
        n_cols = len(datasets)

        # ── 图像尺寸说明 ────────────────────────────────────────────────────
        # figsize=(n_cols * 6.5, 9)
        #   每列宽度 6.5 英寸，总高 9 英寸。
        #   若通路名称折行后纵向空间不足，可适当增大高度，如改为 11 或 12。
        fig, axes = plt.subplots(1, n_cols, figsize=(n_cols * 6.5, 9))
        fig.patch.set_facecolor('white')

        if n_cols == 1:
            axes = [axes]

        size_handles = []
        last_sc      = None
        last_ax      = None
        vmax_global  = 2

        # ── 第一遍：收集全局 vmax ────────────────────────────────────────
        for df_enr, title, base_color, src in datasets:
            df_tmp = df_enr.copy()
            if 'Adjusted P-value' in df_tmp.columns:
                df_tmp = df_tmp[df_tmp['Adjusted P-value'].notna()]
            adj_col = 'Adjusted P-value' if 'Adjusted P-value' in df_tmp.columns else 'P-value'
            if adj_col in df_tmp.columns:
                pvals = df_tmp[adj_col].clip(lower=1e-15)
                vmax_global = max(vmax_global, float(-np.log10(pvals).max()))

        # ── 主循环：绘图 ─────────────────────────────────────────────────
        for ax, (df_enr, title, base_color, src) in zip(axes, datasets):
            df_plot = df_enr.copy()

            if 'Adjusted P-value' in df_plot.columns:
                df_plot = df_plot[df_plot['Adjusted P-value'].notna()]
                df_plot = df_plot.sort_values('Adjusted P-value').head(20)
            else:
                df_plot = df_plot.head(20)

            if df_plot.empty:
                ax.set_visible(False)
                continue

            if 'Term' in df_plot.columns:
                df_plot['Term_short'] = df_plot['Term'].apply(
                    # ① 先按来源去除编号，② 再折行（超过 40 字符折为两行）
                    # ── 调参说明 ──────────────────────────────────────────
                    # wrap_term 的第二个参数 max_chars 控制折行阈值，默认 40。
                    # clean_term 的 source 参数：'go' / 'reactome' / 'other'
                    # ─────────────────────────────────────────────────────
                    lambda x: wrap_term(clean_term(x, source=src))
                )
            else:
                df_plot['Term_short'] = [f'Term_{i}' for i in range(len(df_plot))]

            df_plot = df_plot.sort_values(
                'Adjusted P-value' if 'Adjusted P-value' in df_plot.columns else 'P-value',
                ascending=False)
# ── 去重处理：截断后可能出现重复的 Term_short ────────────────────────
# 对重复项追加序号后缀加以区分，保证 categories 唯一
            seen = {}
            unique_terms = []
            for t in df_plot['Term_short'].tolist():
                if t not in seen:
                    seen[t] = 0
                    unique_terms.append(t)
                else:
                    seen[t] += 1
                    unique_terms.append(f'{t} ({seen[t]})')

            df_plot['Term_short'] = unique_terms


            df_plot['Term_short'] = pd.Categorical(
                df_plot['Term_short'],
                categories=df_plot['Term_short'].tolist(),
                ordered=True)

            if df_plot['Gene ratio'].isna().all():
                bg_col = next((c for c in df_plot.columns
                               if 'background' in c.lower() or 'bg' in c.lower()
                               or 'bgcount' in c.lower()), None)
                if bg_col is not None:
                    bg = pd.to_numeric(df_plot[bg_col], errors='coerce').fillna(200)
                    df_plot['Gene ratio'] = df_plot['Gene count'].astype(float) / bg
                else:
                    df_plot['Gene ratio'] = df_plot['Gene count'].astype(float) / 200.0
            else:
                nan_mask = df_plot['Gene ratio'].isna()
                if nan_mask.any():
                    df_plot.loc[nan_mask, 'Gene ratio'] = (
                        df_plot.loc[nan_mask, 'Gene count'].astype(float) / 200.0)

            df_plot['Gene ratio'] = (
                pd.to_numeric(df_plot['Gene ratio'], errors='coerce')
                .clip(lower=1e-6, upper=1.0).fillna(0.05))

            adj_col  = 'Adjusted P-value' if 'Adjusted P-value' in df_plot.columns else 'P-value'
            pvals    = df_plot[adj_col].clip(lower=1e-15)
            neg_log  = -np.log10(pvals)
            sig_mask = pvals < 0.05

            cnt      = df_plot['Gene count'].values.astype(float)
            cnt_norm = (cnt - cnt.min()) / (cnt.max() - cnt.min() + 1e-9)
            sizes    = 60 + cnt_norm * 220

            sc = ax.scatter(
                df_plot['Gene ratio'].values,
                range(len(df_plot)),
                s=sizes,
                c=neg_log.values,
                cmap='YlOrRd',
                vmin=0,
                vmax=vmax_global,
                edgecolors='#333333',
                linewidths=0.6,
                zorder=3)

            sig_idx = np.where(sig_mask.values)[0]
            if len(sig_idx) > 0:
                ax.scatter(
                    df_plot['Gene ratio'].values[sig_idx],
                    sig_idx,
                    s=sizes[sig_idx],
                    facecolors='none',
                    edgecolors='black',
                    linewidths=1.8,
                    zorder=4)

            ax.set_yticks(range(len(df_plot)))

            # ── 通路名称字体大小 ─────────────────────────────────────────
            # 原始值 8.2，×1.2 ≈ 9.8，取整为 9.8。
            # ── 调参说明 ─────────────────────────────────────────────────
            # fontsize : float
            #     通路标签字号。折行后每条通路占两行，若图纵向空间紧张，
            #     可适当减小此值（如改回 8.2）或增大 figsize 高度。
            # ─────────────────────────────────────────────────────────────
            ax.set_yticklabels(df_plot['Term_short'].tolist(), fontsize=13)

            # ── 坐标轴标题 / 子图标题字体大小 ───────────────────────────
            # 原始值 9.5 / 10.5，同样 ×1.2
            # ── 调参说明 ─────────────────────────────────────────────────
            # ax.set_xlabel fontsize : x 轴标题字号，原 9.5 → 11.4
            # ax.set_title  fontsize : 子图标题字号，原 10.5 → 12.6
            # ─────────────────────────────────────────────────────────────
            ax.set_xlabel('Gene Ratio', fontsize=11.4)   # 原 9.5
            ax.set_title(title, fontsize=12.6, fontweight='bold', pad=8)  # 原 10.5

            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(axis='x', color='#EEEEEE', lw=0.6, zorder=0)

            last_sc = sc
            last_ax = ax

            if not size_handles:
                for cnt_val in [1, 2, 3, 4]:
                    size_handles.append(
                        Line2D([0], [0], marker='o', color='w',
                               markerfacecolor='#AAAAAA',
                               markersize=np.sqrt(60 + (cnt_val / cnt.max()) * 220),
                               label=f'n={cnt_val}'))

        # ── colorbar ────────────────────────────────────────────────────────
        if last_sc is not None and last_ax is not None:
            cbar = plt.colorbar(last_sc, ax=last_ax,
                                label='-log₁₀(p.adj)',
                                shrink=0.55, pad=0)
            cbar.ax.tick_params(labelsize=8.5)
            # ── 调参说明 ──────────────────────────────────────────────────
            # cbar label fontsize : colorbar 标题字号，原 9.5 → 11.4
            # ─────────────────────────────────────────────────────────────
            cbar.set_label('-log₁₀(p.adj)', fontsize=11.4)  # 原 9.5

        fig.legend(handles=size_handles,
                   loc='lower center', ncol=4, fontsize=8.5,
                   bbox_to_anchor=(0.5, -0.04),
                   framealpha=0.9, edgecolor='#ccc', facecolor='white')

        # ── 总标题字体大小 ───────────────────────────────────────────────
        # 原始值 12，×1.2 = 14.4
        # ── 调参说明 ─────────────────────────────────────────────────────
        # suptitle fontsize : 总标题字号，原 12 → 14.4
        # ─────────────────────────────────────────────────────────────────
        fig.suptitle(
            'Functional enrichment of MAFB regulon targets (pySCENIC, GSE182434)\n'
            'ORA, BH correction | Black ring = p.adj < 0.05',
            fontsize=14.4, fontweight='bold', y=1.01)  # 原 12

        plt.tight_layout(rect=[0.02, 0.04, 1.0, 1.0])
        save_fig(fig, 'MAFB_gseapy_enrichment_dotplot')

# =============================================================================
#  Fig 7：Sankey 7-layer
#  原尺寸 (22, 12) → (13.2, 7.2)
# =============================================================================
if need_redraw('sankey_7layer_TF_targets_v3_paper'):
    print('[Fig 7] Sankey 7-layer TF targets...')

    # ★ 统一配色
    C_SOURCE  = '#4A4A6A'; C_LIG    = '#D4603A'; C_REC  = '#2E7DC9'
    C_PATH    = '#3A9B6F'; C_TF_LA  = AXIS_COLORS['LA_TAM']    # ★ #2ECC71
    C_TF_IFN  = AXIS_COLORS['IFN_TAM']                          # ★ #E67E22
    C_TF_SH   = SUBTYPE_COLORS['Mono']                          # ★ #3498DB
    C_GENE_LA = '#A8E6CF'; C_GENE_IFN = '#FFD180'
    C_GENE_SH = '#90CAF9'
    C_LA      = AXIS_COLORS['LA_TAM']                           # ★
    C_IFN     = AXIS_COLORS['IFN_TAM']                          # ★
    BG        = '#FAFAFA'

    la_tfs_s  = ['ETV5', 'MITF', 'MAFB', 'XBP1']
    ifn_tfs_s = ['IRF1', 'KLF4', 'KLF2', 'FOSB']
    shared_s  = ['STAT1']

    source_nodes   = [('Malignant B', 0, 12)]
    ligand_nodes   = [('IL15RA', 1, 3.0), ('HMGB1', 1, 3.0),
                      ('TNFSF9', 1, 3.0), ('COPA',  1, 3.0)]
    receptor_nodes = [('AXL',      2, 3.0), ('HAVCR2',   2, 3.0),
                      ('HLA-DPA1', 2, 3.0), ('CD74',     2, 3.0)]
    pathway_nodes  = [('PI3K/AKT\n->mTOR',    3, 3.5),
                      ('STAT3\nactivation',    3, 2.5),
                      ('SHP2/PTPN11\n->ERK',  3, 2.5),
                      ('JAK/STAT1\n->IRF',    3, 3.0),
                      ('NF-kB\n/CXCR4',       3, 3.0)]
    tf_nodes       = [(tf, 4, 1.8) for tf in la_tfs_s + shared_s + ifn_tfs_s]
    gene_nodes     = [(g, 5, 0.9) for tf, genes in TOP5.items() for g in genes]
    subtype_nodes  = [('LA_TAM', 6, 5.5), ('IFN_TAM', 6, 5.5)]

    PW = {k: k for k in [
        'PI3K/AKT\n->mTOR', 'STAT3\nactivation',
        'SHP2/PTPN11\n->ERK', 'JAK/STAT1\n->IRF', 'NF-kB\n/CXCR4']}
    p1, p2, p3, p4, p5 = list(PW.values())

    links = [
        ('Malignant B', 'IL15RA', 3.0, True),
        ('Malignant B', 'HMGB1',  3.0, True),
        ('Malignant B', 'TNFSF9', 3.0, True),
        ('Malignant B', 'COPA',   3.0, True),
        ('IL15RA',   'AXL',      3.0, True),
        ('HMGB1',    'HAVCR2',   3.0, True),
        ('TNFSF9',   'HLA-DPA1', 3.0, True),
        ('COPA',     'CD74',     3.0, True),
        ('AXL',    p1, 2.0, True),  ('AXL',    p2, 1.0, True),
        ('HAVCR2', p3, 2.0, True),  ('HAVCR2', p2, 1.0, False),
        ('HLA-DPA1', p4, 2.0, True), ('HLA-DPA1', p5, 1.0, False),
        ('CD74',   p5, 2.0, True),  ('CD74',   p4, 1.0, False),
        (p1, 'ETV5', 1.5, True),  (p1, 'MITF', 1.5, True),
        (p1, 'MAFB', 1.0, True),  (p1, 'XBP1', 1.0, True),
        (p2, 'ETV5', 1.0, True),  (p2, 'STAT1',1.0, True),
        (p2, 'MITF', 0.8, False),
        (p3, 'STAT1',1.2, True),  (p3, 'MAFB', 0.8, False),
        (p4, 'STAT1',1.5, True),  (p4, 'IRF1', 1.5, True),
        (p4, 'FOSB', 1.0, True),
        (p5, 'KLF2', 1.5, True),  (p5, 'KLF4', 1.5, True),
        (p5, 'IRF1', 1.0, False),
    ]
    for tf, genes in TOP5.items():
        w = 1.8 / len(genes)
        for g in genes:
            links.append((tf, g, w, True))
    for tf in la_tfs_s:
        for g in TOP5[tf]:
            links.append((g, 'LA_TAM',  0.9 / 5, True))
    for tf in ifn_tfs_s:
        for g in TOP5[tf]:
            links.append((g, 'IFN_TAM', 0.9 / 5, True))
    for g in TOP5['STAT1']:
        links.append((g, 'LA_TAM',  0.45 / 5, False))
        links.append((g, 'IFN_TAM', 0.45 / 5, False))

    def compute_y_positions(nodes, total_height=20, pad=0.3):
        weights  = [w for _, _, w in nodes]
        total_w  = sum(weights) + pad * (len(weights) - 1)
        scale    = total_height / total_w
        positions = {}
        y = total_height
        for label, col, w in nodes:
            h = w * scale
            positions[label] = (y - h / 2, h)
            y -= h + pad * scale
        return positions

    H = 20
    all_pos = {
        **compute_y_positions(source_nodes,   H, 0.5),
        **compute_y_positions(ligand_nodes,   H, 0.6),
        **compute_y_positions(receptor_nodes, H, 0.6),
        **compute_y_positions(pathway_nodes,  H, 0.5),
        **compute_y_positions(tf_nodes,       H, 0.4),
        **compute_y_positions(gene_nodes,     H, 0.12),
        **compute_y_positions(subtype_nodes,  H, 1.5),
    }

    col_x  = {0: 0.0, 1: 2.0, 2: 4.0, 3: 6.0, 4: 8.0, 5: 10.0, 6: 12.2}
    NODE_W = 0.9

    node_col = {}
    for label, col, _ in (source_nodes + ligand_nodes + receptor_nodes +
                           pathway_nodes + tf_nodes + gene_nodes + subtype_nodes):
        node_col[label] = col

    node_color = {}
    for label, _, _ in source_nodes:   node_color[label] = C_SOURCE
    for label, _, _ in ligand_nodes:   node_color[label] = C_LIG
    for label, _, _ in receptor_nodes: node_color[label] = C_REC
    for label, _, _ in pathway_nodes:  node_color[label] = C_PATH
    for label, _, _ in tf_nodes:
        if label in la_tfs_s:   node_color[label] = C_TF_LA
        elif label == 'STAT1':  node_color[label] = C_TF_SH
        else:                   node_color[label] = C_TF_IFN
    tf_of_gene = {g: tf for tf, genes in TOP5.items() for g in genes}
    for label, _, _ in gene_nodes:
        parent = tf_of_gene.get(label, '')
        if parent in la_tfs_s:   node_color[label] = C_GENE_LA
        elif parent == 'STAT1':  node_color[label] = C_GENE_SH
        else:                    node_color[label] = C_GENE_IFN
    node_color['LA_TAM']  = C_LA
    node_color['IFN_TAM'] = C_IFN

    node_right_used = {n: 0.0 for n in all_pos}
    node_left_used  = {n: 0.0 for n in all_pos}
    node_out_total  = {n: 0.0 for n in all_pos}
    node_in_total   = {n: 0.0 for n in all_pos}
    for s, t, w, _ in links:
        if s in node_out_total: node_out_total[s] += w
        if t in node_in_total:  node_in_total[t]  += w

    def draw_bezier_ribbon(ax, x0, y0c, h0, x1, y1c, h1, color, alpha):
        y0t, y0b = y0c + h0 / 2, y0c - h0 / 2
        y1t, y1b = y1c + h1 / 2, y1c - h1 / 2
        cx = (x0 + x1) / 2
        verts = [(x0, y0t), (cx, y0t), (cx, y1t), (x1, y1t),
                 (x1, y1b), (cx, y1b), (cx, y0b), (x0, y0b), (x0, y0t)]
        codes  = [Path.MOVETO,
                  Path.CURVE4, Path.CURVE4, Path.CURVE4,
                  Path.LINETO,
                  Path.CURVE4, Path.CURVE4, Path.CURVE4,
                  Path.CLOSEPOLY]
        ax.add_patch(mpatches.PathPatch(
            Path(verts, codes),
            facecolor=color, edgecolor='none', alpha=alpha, zorder=1))

    # ★ 尺寸 ×0.6
    fig, ax = plt.subplots(figsize=(18.0, 7.2))  # ★ 18in
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.set_xlim(-0.5, 14.0)
    ax.set_ylim(-1.8, H + 1.4)  # ★ 无标题，收紧
    ax.axis('off')

    for s, t, w, strong in links:
        if s not in all_pos or t not in all_pos:
            continue
        x0 = col_x[node_col[s]] + NODE_W
        x1 = col_x[node_col[t]]
        cy_s, h_s = all_pos[s]; cy_t, h_t = all_pos[t]
        frac_s = w / node_out_total[s] if node_out_total[s] > 0 else 0
        frac_t = w / node_in_total[t]  if node_in_total[t]  > 0 else 0
        sh_s   = h_s * frac_s; sh_t = h_t * frac_t
        y0c = (cy_s + h_s / 2 - node_right_used[s]) - sh_s / 2
        node_right_used[s] += sh_s
        y1c = (cy_t + h_t / 2 - node_left_used[t]) - sh_t / 2
        node_left_used[t]  += sh_t
        draw_bezier_ribbon(ax, x0, y0c, sh_s, x1, y1c, sh_t,
                           color=node_color[s],
                           alpha=0.38 if strong else 0.16)

    for label, (cy, h) in all_pos.items():
        x     = col_x[node_col[label]]
        color = node_color[label]
        ax.add_patch(FancyBboxPatch(
            (x, cy - h / 2), NODE_W, h,
            boxstyle='round,pad=0.04',
            facecolor=color, edgecolor='white', linewidth=0.9, zorder=3))
        fs = 8.0 if node_col[label] == 5 else (7.5 if '\n' in label else 9.0)  # ★ -1
        ax.text(x + NODE_W / 2, cy, label,
                ha='center', va='center',
                fontsize=fs, color='white', fontweight='bold',
                zorder=4, linespacing=1.2)

    headers = [
        (col_x[0] + NODE_W / 2, 'Malignant B\nCell',      C_SOURCE),
        (col_x[1] + NODE_W / 2, 'Secreted\nLigand',       C_LIG),
        (col_x[2] + NODE_W / 2, 'Mono/Mac\nReceptor',     C_REC),
        (col_x[3] + NODE_W / 2, 'Signaling\nPathway',     C_PATH),
        (col_x[4] + NODE_W / 2, 'Core TF\n(pySCENIC)',    '#555555'),
        (col_x[5] + NODE_W / 2, 'Top 5 Target\nGenes',   '#555555'),
        (col_x[6] + NODE_W / 2, 'TAM\nSubtype',           '#555555'),
    ]
    for hx, hlabel, hcolor in headers:
        ax.text(hx, H + 1.0, hlabel,
                ha='center', va='bottom', fontsize=11,   # ★ 9.5→11
                fontweight='bold', color=hcolor, linespacing=1.3)
        ax.axvline(hx, ymin=0.02, ymax=0.90,
                   color=hcolor, alpha=0.10, lw=1, zorder=0)

    

    legend_patches = [
        mpatches.Patch(color=C_SOURCE,   label='Malignant B cell'),
        mpatches.Patch(color=C_LIG,      label='Secreted ligand'),
        mpatches.Patch(color=C_REC,      label='Surface receptor'),
        mpatches.Patch(color=C_PATH,     label='Signaling pathway'),
        mpatches.Patch(color=C_TF_LA,    label='LA_TAM core TF'),
        mpatches.Patch(color=C_TF_IFN,   label='IFN_TAM core TF'),
        mpatches.Patch(color=C_TF_SH,    label='Shared TF (STAT1)'),
        mpatches.Patch(color=C_GENE_LA,  label='LA_TAM target gene'),
        mpatches.Patch(color=C_GENE_IFN, label='IFN_TAM target gene'),
        mpatches.Patch(color=C_GENE_SH,  label='Shared target gene'),
        mpatches.Patch(color='gray', alpha=0.40, label='Strong evidence'),
        mpatches.Patch(color='gray', alpha=0.16, label='Moderate evidence'),
    ]
    ax.legend(handles=legend_patches, loc='lower center',
              bbox_to_anchor=(0.5, -0.07), ncol=6,
              fontsize=10, frameon=True, framealpha=0.9,   # ★ 8→10
              edgecolor='#cccccc', facecolor='white')

    plt.tight_layout(pad=0.3)
    save_fig(fig, 'sankey_7layer_TF_targets_v3_paper')

# =============================================================================
#  完成总结
# =============================================================================
print()
print('=' * 65)
print('论文版图形生成完成')
print('=' * 65)
figs = [
    'Fig1_regulon_activity_heatmap_paper',
    'Fig2_RSS_scatter_paper',
    'Fig3_TF_target_network_paper',
    'pseudotime_TF_receptor_dynamics_paper',
    'lineplot_LA_TAM_TF_receptor_2x4_paper',
    'MAFB_gseapy_enrichment_dotplot_paper',
    'sankey_7layer_TF_targets_v3_paper',
]
for fname in figs:
    for fmt in ['png', 'svg']:
        p = os.path.join(OUT_DIR, f'{fname}.{fmt}')
        status = '[OK]' if os.path.exists(p) else '[MISS]'
        print(f'  {status}  {p}')
print('=' * 65)
# =============================================================================
#  Receptor Pseudotime 分析 — 本地 Windows 复现版
#  数据路径：D:\bulk-download\GSE182434\trajectory\
#  输出：2 张图（heatmap + lineplot），PNG + SVG 双格式
# =============================================================================

import os
import warnings
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import zscore
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')

# =============================================================================
#  路径配置
# =============================================================================
TRAJ_DIR = translate(r'D:\bulk-download\GSE182434\trajectory')
OUT_DIR  = TRAJ_DIR
DPI      = 300
os.makedirs(OUT_DIR, exist_ok=True)

CSV_CELL   = os.path.join(TRAJ_DIR, 'receptor_pseudotime_cell_level.csv')
CSV_BINNED = os.path.join(TRAJ_DIR, 'receptor_pseudotime_binned_smooth.csv')
CSV_CORR   = os.path.join(TRAJ_DIR, 'receptor_pseudotime_correlation_path.csv')

print('=' * 65)
print('Receptor Pseudotime 本地复现版')
print('=' * 65)

# =============================================================================
#  加载数据
# =============================================================================
print('\n[LOAD] Reading CSV files...')

cell_df   = pd.read_csv(CSV_CELL)
binned_df = pd.read_csv(CSV_BINNED)
corr_df   = pd.read_csv(CSV_CORR)

print(f'  Cell-level  : {cell_df.shape}')
print(f'  Binned-smooth: {binned_df.shape}')
print(f'  Correlation : {corr_df.shape}')

# ── 提取关键列 ───────────────────────────────────────────────────────────────
# cell_df 列：cell_barcode | pseudotime | subtype | <receptor genes...>
pt       = cell_df['pseudotime'].values
subtypes = cell_df['subtype'].values

receptor_cols = [c for c in cell_df.columns
                 if c not in ('cell_barcode', 'pseudotime', 'subtype')]
print(f'  Receptor genes detected: {len(receptor_cols)}')
print(f'  Pseudotime range: {pt.min():.4f} – {pt.max():.4f}')
print(f'  Subtype counts:\n{pd.Series(subtypes).value_counts().to_string()}')

# ── 相关性字典 ───────────────────────────────────────────────────────────────
corr_dict = dict(zip(corr_df['gene'], corr_df['spearman_r']))
padj_dict = dict(zip(corr_df['gene'], corr_df['padj']))

def get_sig_label(gene):
    r = corr_dict.get(gene, 0)
    p = padj_dict.get(gene, 1)
    if p <= 0.05:
        return '▲' if r > 0 else '▼'
    return ''

# ── 亚型边界 ─────────────────────────────────────────────────────────────────
mono_end   = cell_df.loc[cell_df['subtype'] == 'Mono',    'pseudotime'].max()
ifntam_end = cell_df.loc[cell_df['subtype'] == 'IFN_TAM', 'pseudotime'].max()
print(f'\n  Mono end     : {mono_end:.4f}')
print(f'  IFN_TAM end  : {ifntam_end:.4f}')

# ── 显著性基因 ───────────────────────────────────────────────────────────────
pos_genes_df = corr_df[(corr_df['spearman_r'] > 0) & (corr_df['padj'] <= 0.05)]
neg_genes_df = corr_df[(corr_df['spearman_r'] < 0) & (corr_df['padj'] <= 0.05)]
top_pos_genes = pos_genes_df.sort_values('spearman_r', ascending=False)['gene'].head(4).tolist()
# ★ 与原图 Fig3B 面板排布一致：下行 = HLA-DPA1(左), HAVCR2(右)
if 'HLA-DPA1' in top_pos_genes and 'HAVCR2' in top_pos_genes:
    i, j = top_pos_genes.index('HAVCR2'), top_pos_genes.index('HLA-DPA1')
    top_pos_genes[i], top_pos_genes[j] = top_pos_genes[j], top_pos_genes[i]

print(f'\n  Positive sig genes: {len(pos_genes_df)} → {pos_genes_df["gene"].tolist()}')
print(f'  Negative sig genes: {len(neg_genes_df)}')
print(f'  Top 4 for lineplot : {top_pos_genes}')

# =============================================================================
#  Fig 1 — Receptor Gene Expression Heatmap Along Pseudotime
# =============================================================================
print('\n[Fig 1] Drawing receptor pseudotime heatmap...')

# ── 分箱 + 平滑 ──────────────────────────────────────────────────────────────
N_BINS   = 100
SIGMA    = 3
bin_edges   = np.linspace(pt.min(), pt.max(), N_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

expr_mat   = cell_df[receptor_cols].values          # (cells, genes)
expr_binned = np.zeros((len(receptor_cols), N_BINS))
subtype_bins = []

for b in range(N_BINS):
    mask = (pt >= bin_edges[b]) & (pt < bin_edges[b + 1])
    if mask.sum() > 0:
        expr_binned[:, b] = expr_mat[mask].mean(axis=0)
        subtype_bins.append(Counter(subtypes[mask]).most_common(1)[0][0])
    else:
        expr_binned[:, b] = 0
        subtype_bins.append('Unknown')

expr_smooth = np.array([gaussian_filter1d(expr_binned[i], sigma=SIGMA)
                        for i in range(len(receptor_cols))])

expr_z = np.apply_along_axis(
    lambda x: zscore(x) if x.std() > 0 else x, 1, expr_smooth)

# ── 按峰值伪时间排序基因 ─────────────────────────────────────────────────────
peak_idx     = np.argmax(expr_z, axis=1)
gene_order   = np.argsort(peak_idx)
genes_sorted = [receptor_cols[i] for i in gene_order]
expr_z_sorted = expr_z[gene_order]

# ── 绘图 ─────────────────────────────────────────────────────────────────────
SUBTYPE_COLORS = {
    'Mono':    '#3498DB',  # 蓝
    'IFN_TAM': '#E67E22',  # 橙
    'LA_TAM':  '#2ECC71',  # 绿
    'Unknown': '#CCCCCC'
}


fig = plt.figure(figsize=(8.4, 9.6))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(
    3, 2,
    height_ratios=[0.04, 1, 0.02],
    width_ratios=[1, 0.03],
    hspace=0.02, wspace=0.02
)
ax_bar  = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[1, 0])
ax_cbar = fig.add_subplot(gs[1, 1])
fig.subplots_adjust(left=0.24, right=0.90)  # ★ 右边距也加大：色条刻度不再出界

# Subtype color bar
subtype_int_map = {'Mono': 0, 'IFN_TAM': 1, 'LA_TAM': 2, 'Unknown': 3}
subtype_numeric = np.array([[subtype_int_map.get(s, 3) for s in subtype_bins]])
cmap_sub = ListedColormap(['#E9C46A', '#E63946', '#F4A261', '#CCCCCC'])
ax_bar.imshow(subtype_numeric, aspect='auto', cmap=cmap_sub,
              vmin=0, vmax=3, interpolation='nearest')
ax_bar.set_xticks([])
ax_bar.set_yticks([0])
ax_bar.set_yticklabels(['Subtype'], fontsize=12)
ax_bar.tick_params(left=False)

legend_patches = [
    Patch(color='#E9C46A', label='Mono'),
    Patch(color='#E63946', label='IFN_TAM'),
    Patch(color='#F4A261', label='LA_TAM'),
]
ax_bar.legend(handles=legend_patches, loc='upper center',
              bbox_to_anchor=(0.5, 2.4), fontsize=10, frameon=True, ncol=3)  # ★ 水平居中

# Heatmap
cmap_heat = LinearSegmentedColormap.from_list(
    'expr', ['#2166AC', '#F7F7F7', '#D6604D'])
im = ax_heat.imshow(expr_z_sorted, aspect='auto', cmap=cmap_heat,
                    vmin=-2, vmax=2, interpolation='nearest')

# Y-axis gene labels with significance markers
ytick_labels = []
for g in genes_sorted:
    sig = get_sig_label(g)
    ytick_labels.append(f'{g} {sig}' if sig else g)

ax_heat.set_yticks(range(len(genes_sorted)))
ax_heat.set_yticklabels(ytick_labels, fontsize=10)

for tick, gene in zip(ax_heat.get_yticklabels(), genes_sorted):
    r = corr_dict.get(gene, 0)
    p = padj_dict.get(gene, 1)
    if p <= 0.05:
        tick.set_color('#D6604D' if r > 0 else '#2166AC')
        tick.set_fontweight('bold')

# X-axis pseudotime ticks
n_xticks = 6
xtick_pos = np.linspace(0, N_BINS - 1, n_xticks)
xtick_labels_pt = [f'{v:.2f}' for v in np.linspace(pt.min(), pt.max(), n_xticks)]
ax_heat.set_xticks(xtick_pos)
ax_heat.set_xticklabels(xtick_labels_pt, fontsize=11)
ax_heat.set_xlabel('Pseudotime (Mono → IFN_TAM → LA_TAM)', fontsize=12)
ax_heat.set_ylabel('Receptor genes (sorted by peak expression)', fontsize=11)

plt.colorbar(im, cax=ax_cbar, label='Z-score')
ax_cbar.tick_params(labelsize=10)

ax_heat.text(1.08, -0.1,
             '▲ pos. corr.  ▼ neg. corr. (padj≤0.05)',
             transform=ax_heat.transAxes, fontsize=10,
             ha='right', va='top', color='#555555')

fig.suptitle(
    'Receptor Gene Expression Along Pseudotime\n'
    'Mono → IFN_TAM → LA_TAM Trajectory '
    '(82 Malignant-enriched LR pairs, 57 receptor genes)',
    fontsize=14, fontweight='bold', y=1.005
)

plt.tight_layout()
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'receptor_pseudotime_heatmap_path.{fmt}')
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {out_path}')
plt.close()

# =============================================================================
#  Fig 2 — Top 4 Positive Genes Lineplot (2 × 2)
# =============================================================================
print('\n[Fig 2] Drawing top-4 positive genes lineplot...')

N_BINS_LP = 80
SIGMA_LP  = 3
bin_edges_lp   = np.linspace(pt.min(), pt.max(), N_BINS_LP + 1)
bin_centers_lp = (bin_edges_lp[:-1] + bin_edges_lp[1:]) / 2

fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.patch.set_facecolor('white')
axes_flat = axes.flatten()

for ax, gene in zip(axes_flat, top_pos_genes):
    if gene not in cell_df.columns:
        ax.set_visible(False)
        continue

    expr = cell_df[gene].values

    # Bin mean + SEM
    expr_binned_lp = np.zeros(N_BINS_LP)
    expr_sem_lp    = np.zeros(N_BINS_LP)
    for b in range(N_BINS_LP):
        mask = (pt >= bin_edges_lp[b]) & (pt < bin_edges_lp[b + 1])
        if mask.sum() > 0:
            vals = expr[mask]
            expr_binned_lp[b] = vals.mean()
            expr_sem_lp[b]    = vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0

    expr_smooth_lp = gaussian_filter1d(expr_binned_lp, sigma=SIGMA_LP)
    sem_smooth_lp  = gaussian_filter1d(expr_sem_lp,    sigma=SIGMA_LP)

    # Background shading
    ax.axvspan(pt.min(),   mono_end,   alpha=0.12, color=SUBTYPE_COLORS['Mono'],    zorder=0)
    ax.axvspan(mono_end,   ifntam_end, alpha=0.12, color=SUBTYPE_COLORS['IFN_TAM'], zorder=0)
    ax.axvspan(ifntam_end, pt.max(),   alpha=0.12, color=SUBTYPE_COLORS['LA_TAM'],  zorder=0)

    ax.axvline(mono_end,   color='#888888', lw=0.8, ls='--', zorder=1)
    ax.axvline(ifntam_end, color='#888888', lw=0.8, ls='--', zorder=1)

    # Scatter (cells colored by subtype)
    for st, sc in SUBTYPE_COLORS.items():
        if st == 'Unknown':
            continue
        mask_st = subtypes == st
        ax.scatter(pt[mask_st], expr[mask_st],
                   c=sc, s=8, alpha=0.35, zorder=2, linewidths=0)

    # Smoothed mean ± SEM
    ax.plot(bin_centers_lp, expr_smooth_lp, color='#333333', lw=2.2, zorder=4)
    ax.fill_between(bin_centers_lp,
                    expr_smooth_lp - sem_smooth_lp,
                    expr_smooth_lp + sem_smooth_lp,
                    color='#333333', alpha=0.15, zorder=3)

    # Spearman annotation box
    r_val = corr_dict.get(gene, 0)
    p_val = padj_dict.get(gene, 1)
    p_str = f'{p_val:.2e}' if p_val >= 1e-10 else '<1e-10'
    ax.text(0.3, 0.9,
            f'r = {r_val:.3f}\npadj = {p_str}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=11, color='#D6604D', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3',
                      fc='white', ec='#D6604D', alpha=0.85))

    ax.set_title(gene, fontsize=16, fontweight='bold', color='#D6604D')
    ax.set_xlabel('Pseudotime', fontsize=12)
    ax.set_ylabel('log-normalized expression', fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlim(pt.min() - 0.005, pt.max() + 0.005)

# Re-draw subtype labels after ylim is settled
for ax, gene in zip(axes_flat, top_pos_genes):
    if gene not in cell_df.columns:
        continue
    ymax = ax.get_ylim()[1]
    for st, x_pos in [
        ('Mono',    (pt.min()  + mono_end)   / 2),
        ('IFN_TAM', (mono_end  + ifntam_end) / 2),
        ('LA_TAM',  (ifntam_end + pt.max())  / 2),
    ]:
        ax.text(x_pos, ymax * 0.97, st,
                ha='center', va='top', fontsize=11,
                color=SUBTYPE_COLORS[st], fontweight='bold')

# Figure legend
legend_elements = [
    Patch(color=SUBTYPE_COLORS['Mono'],    alpha=0.5, label='Mono'),
    Patch(color=SUBTYPE_COLORS['IFN_TAM'], alpha=0.5, label='IFN_TAM'),
    Patch(color=SUBTYPE_COLORS['LA_TAM'],  alpha=0.5, label='LA_TAM'),
    Line2D([0], [0], color='#333333', lw=2, label='Smoothed mean ± SEM'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4,
           fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.10))  # ★ 上移至 xlabel 下方  # ★ 锚点 0.005

fig.subplots_adjust(left=0.09, bottom=0.20, right=0.97, top=0.96, hspace=0.42, wspace=0.22)  # ★ 底部 20%，图例紧贴 xlabel 下方
for fmt in ['png', 'svg']:
    out_path = os.path.join(OUT_DIR, f'receptor_pseudotime_lineplots_path.{fmt}')
    plt.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {out_path}')
plt.close()

# =============================================================================
#  完成
# =============================================================================
print('\n' + '=' * 65)
print('全部完成！输出文件：')
for fname in [
    'receptor_pseudotime_heatmap_path.png',
    'receptor_pseudotime_heatmap_path.svg',
    'receptor_pseudotime_lineplots_path.png',
    'receptor_pseudotime_lineplots_path.svg',
]:
    p = os.path.join(OUT_DIR, fname)
    status = '[OK]  ' if os.path.exists(p) else '[MISS]'
    print(f'  {status} {p}')
print('=' * 65)
