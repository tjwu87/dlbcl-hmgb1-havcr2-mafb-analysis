# -*- coding: utf-8 -*-
"""以"强制单一 gate"方式运行 04_pySCENIC/04_01_pyscenic_downstream_stats_figures.py 的指定小节。

背景：主脚本 04_01 含 7 个 need_redraw gate（按产物缓存跳过）。本运行器把指定 gate
强制为 True 后整体执行主脚本，从而只重渲该节——替代早期维护的两份 1650 行全文副本
（04_01_sankey_only.py / 04_01_lineplot_only.py，已删）。

用法（经 retrofit 包装，与运行主脚本完全一致）：
    # bash:  export PYSCENIC_FORCE_GATE=sankey_7layer_TF_targets_v3_paper
    # cmd :  set PYSCENIC_FORCE_GATE=sankey_7layer_TF_targets_v3_paper
    python tools/run_with_retrofit.py tools/04_01_gate_runner.py 6.3 8.0 \
        --design --strip-titles --out results/uav_panels/04_01_v3

可用 gate 名（= 04_01 的 need_redraw 名，即产物名）：
    Fig1_regulon_activity_heatmap_paper | Fig2_RSS_scatter_paper |
    Fig3_TF_target_network_paper | pseudotime_TF_receptor_dynamics_paper |
    lineplot_LA_TAM_TF_receptor_2x4_paper | MAFB_gseapy_enrichment_dotplot |
    sankey_7layer_TF_targets_v3_paper
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# 主脚本位置解析（整理后已移入 section3，保留旧结构回退 + 全局兜底搜索）─────
_CANDS = [
    os.path.join(REPO, '04_pySCENIC', '04_01_pyscenic_downstream_stats_figures.py'),
    os.path.join(REPO, 'section3_HAVCR2-MAFB调控模块', 'code', '04_pySCENIC',
                 '04_01_pyscenic_downstream_stats_figures.py'),
]
MAIN = next((p for p in _CANDS if os.path.exists(p)), None)
if MAIN is None:
    import glob
    _hits = [
        p for p in glob.glob(
            os.path.join(REPO, '**', '04_pySCENIC',
                         '04_01_pyscenic_downstream_stats_figures.py'),
            recursive=True)
        if '历史版本' not in p and 'zenodo' not in p and '_backup' not in p
    ]
    MAIN = _hits[0] if _hits else None
if MAIN is None:
    sys.exit('[04_01_gate_runner] 找不到主脚本 04_01_pyscenic_downstream_stats_figures.py')

gate = os.environ.get('PYSCENIC_FORCE_GATE', '')
if not gate:
    sys.exit('[04_01_gate_runner] 未设置 PYSCENIC_FORCE_GATE 环境变量')

src = open(MAIN, encoding='utf-8').read()
target = f"if need_redraw('{gate}'):"
if target not in src:
    sys.exit(f'[04_01_gate_runner] gate 不存在: {gate}')
src = src.replace(target, 'if True:', 1)

# 主脚本内部用 __file__ 定位自身目录 -> 指回主脚本真实路径
g = {'__name__': '__main__', '__file__': MAIN}
print(f'[04_01_gate_runner] 强制 gate: {gate}')
exec(compile(src, MAIN, 'exec'), g)
