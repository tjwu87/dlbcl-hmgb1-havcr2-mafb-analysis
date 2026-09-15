# -*- coding: utf-8 -*-
"""Fig1 一键重建：面板重渲 -> 矢量裁切 -> 合成（布局由 results/assembled/fig1.json 固定）。

用法（在仓库根目录 DLBCL_HMGB1_HAVCR2_MAFB 下）：
  python tools/rebuild_fig1.py          # 全量：重渲 1a/1b/1c/1f + 裁切 + 合成
  python tools/rebuild_fig1.py --fast   # 复用已有 1a/1c 原始输出，只重跑 1b/1f + 裁切 + 合成

环境要求：
  - DLBCL_DATA_ROOT=D:/bulk-download（脚本内已设默认）
  - 1a/1b/1c 用 Python310 系统解释器（scanpy 环境）；其余与主脚本同解释器
面板来源：
  a = 10_00_fig1a_only.py -> umap_celltypes2.pdf（矢量）      [慢，全量 SC 流程]
  b = fig1b_dotplot.py -> fig1b.pdf（手写点图，纯矢量）
  c = 10_00 截断脚本 -> chromosome_heatmap.pdf（矢量）
  d = assembled/fig1/fig1d.pdf（svglib 矢量资产，静态）
  e = assembled/fig1/fig1e.pdf（静态资产）
  f = fig1f_volcano.py -> volcano_malignant_vs_normal_monomac.pdf（02_02 Fig6 提取重渲）
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

fast = '--fast' in sys.argv
V2 = 'results/uav_panels/fig1a_v2'


def run(cmd):
    print('>>', ' '.join(str(c) for c in cmd), flush=True)
    r = subprocess.run(cmd, cwd=REPO)
    if r.returncode != 0:
        raise SystemExit(f'[FAILED] exit={r.returncode}: {cmd}')


# ── 1) 面板重渲 ──
if not fast:
    # 1a + 1c 原始输出（10_00 是截断脚本：umap 保存后继续画 chromosome_heatmap 再退出）
    run([PY310, RETRO, 'tools/10_00_truncated_fig1a_fig1c.py', '6.3', '8.0',
         '--design', '--out', V2])

run([PY310, 'tools/fig1b_dotplot.py'])                       # 1b
run([PY310, RETRO, 'tools/fig1f_volcano.py', '--plain',      # 1f
     '--out', 'results/uav_panels/fig1f_v2'])

# ── 2) 矢量裁切（参数是定稿值，勿改：1a 顶 pad 大防截断）──
vector_crop(f'{V2}/umap_celltypes2.pdf', f'{V2}/umap_celltypes2_crop.pdf',
            pad_pt=10, th=252)
vector_crop(f'{V2}/chromosome_heatmap.pdf',
            'results/assembled/fig1/fig1c_tight.pdf', pad_pt=4)
vector_crop('results/assembled/fig1/fig1d.pdf',
            'results/assembled/fig1/fig1d_tight.pdf', pad_pt=4)

# ── 3) 合成（布局固定在 fig1.json）──
run([sys.executable, 'tools/fig_compose_pdf.py', 'results/assembled/fig1.json'])
print('[DONE] Fig1 -> results/assembled/fig1/Fig1_composed.pdf')
