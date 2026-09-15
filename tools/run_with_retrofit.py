"""
以出图改造层运行既有画图脚本（路线 B 的执行入口）
================================================================================
在不改动原脚本内部逻辑的前提下，把画布缩到版面尺寸并同步放大字号，
使输出的图在 LaTeX 中 1:1 显示时，文字就是设定的字号。

用法：
    python tools/run_with_retrofit.py <script.py> [target_width_in] [font_pt]

示例：
    # Fig4：按两栏布局（每个 panel 3.2 in 宽），最终字号 7 pt
    python tools/run_with_retrofit.py 05_CellOracle/05_90_reproduce_figures_local.py 3.2 7

    # 单栏大图（panel 6.3 in 宽），字号 8 pt
    python tools/run_with_retrofit.py 04_pySCENIC/04_01_pyscenic_regulons.py 6.3 8
================================================================================
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    target = Path(sys.argv[1])
    if not target.exists():
        print(f"脚本不存在：{target}")
        sys.exit(1)

    # 解析参数：位置参数 [宽度] [字号] [输出目录]，另支持 --out <目录>
    pos, out_flag, design, plain, strip = [], None, False, False, False
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--out" and i + 1 < len(sys.argv):
            out_flag = sys.argv[i + 1]
            i += 2
            continue
        if sys.argv[i] == "--design":
            design = True
            i += 1
            continue
        if sys.argv[i] == "--plain":
            plain = True
            i += 1
            continue
        if sys.argv[i] == "--strip-titles":
            strip = True
            i += 1
            continue
        pos.append(sys.argv[i])
        i += 1

    # 第 1 个位置参数 = 这张图在论文中的最终显示宽度（inch），不是画布宽度。
    # 画布保持脚本原尺寸，字号与图形元素按此反推同步放大。
    width = float(pos[1]) if len(pos) > 1 else 6.3
    font = float(pos[2]) if len(pos) > 2 else 7.0

    # 输出目录：默认 results/figures/<脚本名>，可覆盖（--out 优先）。
    # 目的是绝不覆盖用户已有的原始图片。
    out = (out_flag
           if out_flag
           else (pos[3] if len(pos) > 3
                 else str(REPO / "results" / "figures" / target.stem)))

    import config.retrofit as retrofit
    retrofit.install(target_display_in=width, target_font_pt=font,
                     output_dir=out, design=design, plain=plain,
                     strip_titles=strip)

    print(f"[retrofit] 最终显示宽 {width} in，排版后字号 {font} pt")
    print("[retrofit] 模式：" + ("原样输出（保持脚本原始尺寸与字号）" if plain else
                                ("定尺设计（画布=最终显示尺寸）" if design else "缩画布保字号")))
    print(f"[retrofit] 输出重定向到 {out}")
    print(f"[retrofit] 运行 {target}")
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()
