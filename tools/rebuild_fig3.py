# -*- coding: utf-8 -*-
"""Fig3 一键重建：面板重渲 -> 矢量裁切 -> 合成（布局由 results/assembled/fig3.json 固定）。

用法（在仓库根目录 DLBCL_HMGB1_HAVCR2_MAFB 下）：
  python tools/rebuild_fig3.py
环境要求：DLBCL_DATA_ROOT=D:/bulk-download（脚本内已设默认）

面板来源（经 tools/04_01_gate_runner.py 强制单一 gate 重渲，替代早期两份 1650 行副本）：
  a = gate:sankey -> receptor_pseudotime_heatmap_path.pdf（left=0.24/right=0.90 边距）
  b = gate:lineplot -> receptor_pseudotime_lineplots_path.pdf（subplots_adjust 绝对定位）
  c = pngpdf/Fig1_regulon_activity_heatmap_paper.pdf（静态资产）
  d = gate:Fig2_RSS -> Fig2_RSS_scatter_paper.pdf（备用静态：suoxiao_rss 目录）
  e = gate:sankey -> sankey_7layer_TF_targets_v3_paper.pdf（18in 宽画布）
  f = gate:lineplot -> lineplot_LA_TAM_TF_receptor_2x4_paper.pdf（11.6x5.32in、hspace=0.13、
      图例 6pt@bbox(0.5,0.07)、透明背景）
  g = gate:MAFB_dotplot -> MAFB_gseapy_enrichment_dotplot.pdf（ratios [0.6,0.6,0.6]、wspace 0.85，
      输出重定向 suoxiao_dot 目录，与 fig3.json 一致）
g 在 fig3.json 中 x_mm 已含 -3mm 左移；f-g 行距 2.5mm。
"""
import os
import sys
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY310 = r'C:/Users/admin/AppData/Local/Programs/Python/Python310/python.exe'
if not os.path.exists(PY310):
    PY310 = sys.executable  # ★ 新机器：回退到运行本脚本的解释器（需已装 scanpy 等）
RETRO = os.path.join(REPO, 'tools', 'run_with_retrofit.py')
os.chdir(REPO)
os.environ.setdefault('DLBCL_DATA_ROOT', 'D:/bulk-download')

sys.path.insert(0, os.path.join(REPO, 'tools'))
from vector_crop_util import vector_crop

V3 = 'results/uav_panels/04_01_v3'


def run(cmd, env_extra=None):
    print('>>', ' '.join(str(c) for c in cmd), flush=True)
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    r = subprocess.run(cmd, cwd=REPO, env=env)
    if r.returncode != 0:
        raise SystemExit(f'[FAILED] exit={r.returncode}: {cmd}')


# ── 1) 面板重渲（gate runner 只强制重渲指定节）──
GATES = [
    ('sankey_7layer_TF_targets_v3_paper', None),      # 3a + 3e
    ('lineplot_LA_TAM_TF_receptor_2x4_paper', None),  # 3b + 3f
    ('MAFB_gseapy_enrichment_dotplot', 'results/uav_panels/suoxiao_dot'),  # 3g
]
# 3d 用已确认版式的单节提取件（主文件 gate 版式与论文定稿不同）
run([PY310, RETRO, 'tools/fig3d_rss_extract.py', '--plain',
     '--out', 'results/uav_panels/suoxiao_rss'])
for gate, out_override in GATES:
    out = out_override or V3
    run([PY310, RETRO, 'tools/04_01_gate_runner.py', '6.3', '8.0',
         '--design', '--strip-titles', '--out', out],
        env_extra={'PYSCENIC_FORCE_GATE': gate})

# ── 2) 矢量裁切（定稿参数）──
vector_crop(f'{V3}/receptor_pseudotime_heatmap_path.pdf',
            f'{V3}/heatmap_pad.pdf', pad_pt=12)                    # 3a 呼吸空间
vector_crop(f'{V3}/lineplot_LA_TAM_TF_receptor_2x4_paper.pdf',
            f'{V3}/lineplot_2x4_trim.pdf', pad_pt=4)               # 3f 顶部残带

# ── 3) 合成 ──
run([sys.executable, 'tools/fig_compose_pdf.py', 'results/assembled/fig3.json'])
print('[DONE] Fig3 -> results/assembled/fig3/Fig3_composed.pdf')
