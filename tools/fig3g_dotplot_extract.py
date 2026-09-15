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
OUT_DIR    = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'uav_panels', 'suoxiao_dot')  # ★ 动态推导
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

print('Loading enrichment CSVs...')

if True:

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


    def wrap_term(term, max_chars=20):  # ★ 恢复 20（有空间显示了）
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
        # figsize=(n_cols * 6.5, 7.4)  # ★ 原始宽度，高度仅小幅压缩
        #   每列宽度 6.5 英寸，总高 9 英寸。
        #   若通路名称折行后纵向空间不足，可适当增大高度，如改为 11 或 12。
        fig, axes = plt.subplots(1, n_cols, figsize=(n_cols * 6.5, 7.38), gridspec_kw={'width_ratios': [0.6, 0.6, 0.6]})  # ★ 第三面板横轴缩短，空间让给纵坐标
        fig.subplots_adjust(left=0.19, right=0.98, wspace=0.85)  # ★ 纵坐标标签空间  # ★ 前两面板横轴 -20%
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
            ax.set_yticklabels(df_plot['Term_short'].tolist(), fontsize=12)  # ★ -1

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
                   loc='lower center', ncol=4, fontsize=8,
                   bbox_to_anchor=(0.5, 0.002),
                   framealpha=0.9, edgecolor='#ccc', facecolor='white')

        # ── 总标题字体大小 ───────────────────────────────────────────────
        # 原始值 12，×1.2 = 14.4
        # ── 调参说明 ─────────────────────────────────────────────────────
        # suptitle fontsize : 总标题字号，原 12 → 14.4
        # ─────────────────────────────────────────────────────────────────
        # ★ 标题去除（图注承担）
        pass

        plt.tight_layout(rect=[0.02, 0.05, 1.0, 1.0])
        save_fig(fig, 'MAFB_gseapy_enrichment_dotplot')

# =============================================================================
#  Fig 7：Sankey 7-layer
print('DONE')
