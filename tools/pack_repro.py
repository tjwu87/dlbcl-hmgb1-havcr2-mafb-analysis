# -*- coding: utf-8 -*-
"""把 fig1/fig3 复现所需的全部文件收集到 fig1_fig3_repro_pack/ 目录。

用法（仓库根目录下）：
  python tools/pack_repro.py             # 代码 + 布局 + 静态资产（不含数据）
  python tools/pack_repro.py --with-data # 额外复制 D:/bulk-download 中被引用的数据文件

之后把整个 fig1_fig3_repro_pack/ 拷到新电脑，按 docs/复现指南_fig1_fig3.md 操作。
"""
import os
import shutil
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = 'D:/bulk-download'
PACK = os.path.join(REPO, 'fig1_fig3_repro_pack')

CODE = [
    'tools/rebuild_fig1.py', 'tools/rebuild_fig3.py', 'tools/vector_crop_util.py',
    'tools/fig_compose_pdf.py', 'tools/04_01_gate_runner.py', 'tools/run_with_retrofit.py',
    'tools/10_00_truncated_fig1a_fig1c.py', 'tools/fig1b_dotplot.py',
    'tools/fig1f_volcano.py', 'tools/fig3d_rss_extract.py',
    'tools/fig3g_dotplot_extract.py',
    'config/paths.py', 'config/retrofit.py',
    'results/assembled/fig1.json', 'results/assembled/fig3.json',
    'docs/复现指南_fig1_fig3.md',
]
ASSETS = [
    'results/assembled/fig1/fig1d.pdf',      # 1d 静态矢量资产
    'results/assembled/fig1/fig1e.pdf',      # 1e 静态矢量资产
    'results/uav_panels/pngpdf/Fig1_regulon_activity_heatmap_paper.pdf',   # 3c 静态
    'results/uav_panels/suoxiao_rss/Fig2_RSS_scatter_paper.pdf',           # 3d 静态
]
# 3f 的透明背景版输入 PDF（fig1f_volcano.py 的提取源不在本仓库，见指南"静态资产"节）
ASSETS += ['results/uav_panels/fig1f_v2/volcano_malignant_vs_normal_monomac.pdf']

DATA_FILES = [
    'GSE182434/adata_processed.h5ad',
    'GSE182434/comparison/cluster_markers.csv',
    'GSE182434/annotation',
    'GSE182434/umap',
    'GSE182434/cnv/aucell_metabolic_scores.csv',
    'GSE182434/cnv/malignancy_classification.csv',
    'GSE182434/cellcomm/volcano_data_malignant_vs_normal_monomac.csv',
    'GSE182434/trajectory/receptor_pseudotime_cell_level.csv',
    'GSE182434/scenic/data_AUC_matrix.csv',
    'GSE182434/scenic/data_cell_metadata.csv',
    'GSE182434/scenic/data_RSS_scores.csv',
    'GSE182434/scenic/data_GRNBoost2_adjacencies.csv',
    'GSE182434/scenic/data_RcisTarget_regulon_targets.csv',
    'GSE182434/scenic/data_core_TF_AUC_and_receptor_expression.csv',
    'GSE182434/scenic/data_TF_receptor_intersection.csv',
    'GSE182434/scenic/data_pseudotime_spearman_correlations.csv',
    'GSE182434/scenic/data_top30_RSS_per_subtype.csv',
    'GSE182434/scenic/MAFB_gseapy_GO_BP_enrichment.csv',
    'GSE182434/scenic/MAFB_gseapy_KEGG_enrichment.csv',
    'GSE182434/scenic/MAFB_gseapy_Reactome_enrichment.csv',
    'Macro_mono_cellmarker.xlsx',
]


def copy(src_rel, root_from=REPO, root_to=None):
    src = os.path.join(root_from, src_rel)
    if not os.path.exists(src):
        print('[MISS]', src_rel)
        return False
    dst = os.path.join(root_to or PACK, src_rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.isdir(src):
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)
    return True


def main():
    with_data = '--with-data' in sys.argv
    missing = []
    for f in CODE + ASSETS:
        if not copy(f):
            missing.append(f)
    if with_data:
        for f in DATA_FILES:
            if not copy(f, root_from=DATA):
                missing.append('DATA: ' + f)
    print()
    print('打包完成 ->', os.path.abspath(PACK))
    if missing:
        print('缺失清单（需人工补）:')
        for m in missing:
            print('  ', m)


if __name__ == '__main__':
    main()
