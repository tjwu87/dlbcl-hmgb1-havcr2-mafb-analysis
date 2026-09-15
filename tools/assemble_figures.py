"""
主图组装（路线 B 的最后一环）
================================================================================
把重出的 panel 级矢量图按网格组装成主图的 LaTeX 代码。

为什么用 LaTeX 而不是拼成一张大图：
  - panel 保持矢量（SVG→PDF），放大不失真
  - 每个 panel 的物理尺寸固定为 3.2 in，排版时 1:1，文字就是设定的字号
  - 调整顺序/增删 panel 只需改几行，不必重新跑脚本

用法：
    python tools/assemble_figures.py --scan          # 扫描已重出的 panel，生成配置模板
    python tools/assemble_figures.py                 # 按 config/figure_layout.yaml 生成 LaTeX
================================================================================
"""

from __future__ import annotations

import argparse
from pathlib import Path

import re
import yaml

REPO = Path(__file__).resolve().parents[1]
FIG_ROOT = REPO / "results" / "figures"
LAYOUT = REPO / "config" / "figure_layout.yaml"
OUT_DIR = REPO / "docs" / "figure_latex"

# 脚本输出目录 → 主图（用于 --scan 时给出建议）
SCRIPT_TO_FIGURE = {
    "02_02_liana_visualization": "Fig1",
    "10_00_all_panels_GSE182434": "Fig1",
    "03_01_paga_trajectory": "Fig2",
    "04_01_pyscenic_regulons": "Fig3",
    "celloracle": "Fig4",
    "05_90_reproduce_figures_local": "Fig4",
    "06_01_spatial_analysis": "Fig5",
    "07_01_prognosis_main": "Fig6",
    "07_02_lasso_cv_curves": "Fig6",
    "08_01_ic50_and_boltz2": "Fig7",
    "01_01_qc_integration": "FigS1",
    "10_11_FigS8_supplementary": "FigS8",
    "10_10_FigS2_supplementary": "FigS2",
}


def classify(stem: str, script: str) -> str:
    """判断一个 panel 属于哪张图。

    优先级：文件名里的补充图编号（如 S5A_xxx → FigS5）
            > 脚本与图的固定映射
    """
    m = re.match(r"^S(\d)[A-Za-z]?[_-]", stem)
    if m:
        return f"FigS{m.group(1)}"
    m = re.match(r"^Fig[S]?(\d)", stem)
    if m and stem.startswith("FigS"):
        return f"FigS{m.group(1)}"
    return SCRIPT_TO_FIGURE.get(script, "未指定")


def scan_panels() -> dict:
    """扫描已重出的 panel（svg 与 png 都算，按 stem 去重），按图归类。"""
    out: dict = {}
    for d in sorted(FIG_ROOT.iterdir()):
        if not d.is_dir():
            continue
        stems = sorted({p.stem for p in list(d.glob("*.svg")) + list(d.glob("*.png"))})
        for stem in stems:
            fig = classify(stem, d.name)
            out.setdefault(fig, {}).setdefault(d.name, []).append(stem)
    return out


def default_config() -> dict:
    """按扫描结果生成默认配置：panel 按文件名字母序，2 列。

    注意：真实顺序需对照图注（docs/figure_panel_map.md）人工调整。
    """
    found = scan_panels()
    cfg = {}
    for fig, scripts in found.items():
        if fig == "未指定":
            continue
        panels = []
        for script in sorted(scripts):
            for stem in scripts[script]:
                panels.append(f"results/figures/{script}/{stem}")
        cfg[fig] = {
            "cols": 2,
            "width_in": 3.2,
            "note": "panel 顺序为文件名字母序，需对照图注调整",
            "panels": panels,
        }
    return cfg


def build_latex(fig: str, panels: list, cols: int, width_in: float = 3.2,
                label: str | None = None) -> str:
    """生成单张主图的 LaTeX（minipage 网格）。"""
    if not panels:
        return ""
    label = label or f"fig:{fig.lower()}"
    col_w = f"{width_in}in"
    lines = [
        f"% ===== {fig} （{len(panels)} 个 panel，{cols} 列）=====",
        "\\begin{figure}[p]",
        "\\centering",
    ]
    for i, p in enumerate(panels):
        # panel 统一引用 pdf/ 子目录下的 PDF（由 tools/vectorize_panels.py 生成）
        pp = Path(p)
        # LaTeX 路径一律用正斜杠
        src = str(pp.parent / "pdf" / f"{pp.stem}.pdf").replace("\\", "/")
        lines.append(f"\\begin{{minipage}}[b]{{{col_w}}}")
        lines.append("\\centering")
        lines.append(f"  \\includegraphics[width=\\linewidth]{{{src}}}")
        lines.append(f"\\end{{minipage}}\\hfill" if (i + 1) % cols else "\\end{minipage}")
        if (i + 1) % cols == 0 and i + 1 < len(panels):
            lines.append("\\vspace{0.35em}")
            lines.append("")
    lines += [
        "\\vspace{0.4em}",
        "%% 图注：请按 panel 实际顺序填写（对照 docs/figure_panel_map.md）",
        "\\figurecaptiontext{...}",
        f"\\label{{{label}}}",
        "\\end{figure}",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true", help="只扫描并生成配置模板")
    args = ap.parse_args()

    if args.scan:
        found = scan_panels()
        print("已重出的 panel（按主图分组）：\n")
        for fig in sorted(found):
            print(f"## {fig}")
            for script, files in found[fig].items():
                print(f"   {script}: {len(files)} 个")
                for f in files[:6]:
                    print(f"      - {f}")
                if len(files) > 6:
                    print(f"      - ... 共 {len(files)} 个")
            print()
        return

    if not LAYOUT.exists():
        cfg = default_config()
        LAYOUT.parent.mkdir(parents=True, exist_ok=True)
        LAYOUT.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=True),
                          encoding="utf-8")
        print(f"已生成默认配置 {LAYOUT}（panel 按字母序，需人工调整顺序）")
    else:
        cfg = yaml.safe_load(LAYOUT.read_text(encoding="utf-8"))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for fig, spec in cfg.items():
        tex = build_latex(fig, spec.get("panels", []),
                          spec.get("cols", 2),
                          spec.get("width_in", 3.2),
                          spec.get("label"))
        if not tex:
            continue
        (OUT_DIR / f"{fig}.tex").write_text(tex, encoding="utf-8")
        print(f"[生成] {OUT_DIR / fig}.tex  ({len(spec.get('panels', []))} panel)")


if __name__ == "__main__":
    main()
