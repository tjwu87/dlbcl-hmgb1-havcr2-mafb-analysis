"""
统一出图规范（解决审稿人"图模糊/字太小"的核心文件）
================================================================================
背景
--------------------------------------------------------------------------------
期刊模板 sn-jnl 的版面尺寸为 text={160mm, 216mm}，即
    \\textwidth  = 160 mm = 6.30 inch
    \\textheight = 216 mm = 8.50 inch

旧脚本的做法是：figsize 开到 13×13 ~ 20×6.5 inch，字号用 matplotlib 默认值
(10 pt)，再以 300 dpi 导出。这样的图在 LaTeX 里被 width=\\textwidth 缩放到
6.3 inch，缩放系数约 0.116，图上 10 pt 的字在 PDF 中只剩约 1.2 pt —— 这才是
"模糊"的真正原因，与像素多少无关（旧图 6520 px 但元数据只有 120 dpi）。

正确做法
--------------------------------------------------------------------------------
figsize 直接按最终版面尺寸给，字号按最终字号给，1:1 不缩放。

    from config.plot_style import apply_paper_style, panel_figsize, save_figure

    apply_paper_style()                       # 字号 8 pt 起
    fig, axes = plt.subplots(2, 3, figsize=panel_figsize(3, 2))
    ...
    save_figure(fig, FIG_DIR / "Fig1")        # 同时写出 PDF(矢量) + PNG(600 dpi)

导出的 PDF 是矢量格式，放大任意倍数都不失真，优先用于投稿。
================================================================================
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, Tuple

import matplotlib as mpl

# ── 版面尺寸（来自 sn-jnl.cls: text={160mm, 216mm}）──────────────────────────
TEXTWIDTH_IN: float = 160 / 25.4      # 6.30 inch
TEXTHEIGHT_IN: float = 216 / 25.4     # 8.50 inch

# 图注与正文留白：单张图建议不超过这个高度，否则需要整页浮动
SAFE_HEIGHT_IN: float = 7.60

# ── 字号（pt，即 PDF 中的真实字号，缩放后不变）──────────────────────────────
BASE_FONT_PT: float = 8.0
MIN_FONT_PT: float = 6.0              # 低于此值在印刷品上不可读，不要再调小

FONT_SIZES = {
    "panel_label": BASE_FONT_PT + 3,  # a / b / c 面板标签，加粗
    "title": BASE_FONT_PT + 1,
    "axis_label": BASE_FONT_PT + 0.5,
    "tick": BASE_FONT_PT - 0.5,
    "legend": BASE_FONT_PT - 0.5,
    "legend_title": BASE_FONT_PT,
    "annotation": BASE_FONT_PT - 0.5,  # 统计标注 P 值等，仍须 ≥ MIN_FONT_PT
}

# ── 导出 ─────────────────────────────────────────────────────────────────────
DEFAULT_DPI: int = 600                # 位图分辨率；矢量 PDF 不受此影响
DEFAULT_FORMATS: Tuple[str, ...] = ("pdf", "png")


def apply_paper_style(base_pt: float = BASE_FONT_PT, font_family: str = "sans-serif") -> None:
    """把 matplotlib 全局样式设为论文出图样式。"""
    scale = base_pt / BASE_FONT_PT
    sizes = {k: max(v * scale, MIN_FONT_PT) for k, v in FONT_SIZES.items()}

    mpl.rcParams.update({
        "font.family": font_family,
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": base_pt,
        "axes.labelsize": sizes["axis_label"],
        "axes.titlesize": sizes["title"],
        "axes.linewidth": 0.6,
        "axes.labelpad": 3,
        "axes.titlepad": 6,
        "xtick.labelsize": sizes["tick"],
        "ytick.labelsize": sizes["tick"],
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": sizes["legend"],
        "legend.title_fontsize": sizes["legend_title"],
        "legend.frameon": False,
        "legend.handlelength": 1.2,
        "legend.borderpad": 0.3,
        "figure.dpi": 100,             # 屏幕显示用，不影响 savefig
        "savefig.dpi": DEFAULT_DPI,
        "savefig.bbox": "tight",       # 仅裁白边，不改变字号
        "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42,            # 嵌入 TrueType，避免字体替换
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "patch.linewidth": 0.5,
        "grid.linewidth": 0.4,
        "image.interpolation": "none",
    })
    return None


def panel_figsize(
    cols: int,
    rows: int,
    panel_aspect: float = 1.0,
    width_frac: float = 1.0,
    wspace: float = 0.55,
    hspace: float = 0.65,
    top_space_in: float = 0.25,
) -> Tuple[float, float]:
    """按 panel 网格反推 figsize，保证每个 panel 有足够物理尺寸。

    参数
    ----
    cols, rows : panel 网格列数/行数
    panel_aspect : 单个 panel 的 宽/高 比（1.0 方形；1.4 常见横向；0.8 竖向）
    width_frac : 占 \\textwidth 的比例（默认铺满）
    wspace/hspace : 以字号为单位的留白系数，越大越宽松
    """
    total_w = TEXTWIDTH_IN * width_frac
    # 留白按字号折算成 inch：1 pt = 1/72 inch
    gap_x = (BASE_FONT_PT * wspace) / 72.0 * (cols - 1)
    gap_y = (BASE_FONT_PT * hspace) / 72.0 * (rows - 1)
    panel_w = (total_w - gap_x) / cols
    panel_h = panel_w / panel_aspect
    total_h = panel_h * rows + gap_y + top_space_in
    return (round(total_w, 3), round(total_h, 3))


def check_figsize(size: Sequence[float], name: str = "") -> None:
    """出图前自检：宽度不超版面，高度不超安全高度。返回提示字符串。"""
    w, h = size
    msgs = []
    if w > TEXTWIDTH_IN * 1.02:
        msgs.append(f"宽度 {w:.2f}in 超过 \\textwidth({TEXTWIDTH_IN:.2f}in)，会被缩小")
    if h > SAFE_HEIGHT_IN:
        msgs.append(f"高度 {h:.2f}in 超过安全高度 {SAFE_HEIGHT_IN:.2f}in，建议拆页或改用整页浮动")
    if msgs:
        print(f"[plot_style] {name} " + "; ".join(msgs))
    return "; ".join(msgs)


def save_figure(
    fig,
    out_stem,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = DEFAULT_DPI,
    verbose: bool = True,
) -> list:
    """按统一规格导出。out_stem 不带扩展名，例如 FIG_DIR/'Fig1'。

    PDF 为矢量，用于投稿；PNG 为 600 dpi 位图，用于单独上传与预览。
    """
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = out_stem.with_suffix("." + fmt)
        fig.savefig(path, dpi=dpi, transparent=False)
        written.append(path)
        if verbose:
            size_in = "×".join(f"{v:.2f}in" for v in fig.get_size_inches())
            print(f"[saved] {path.name}  ({size_in}, {dpi} dpi)")
    return written


def add_panel_labels(axes, labels=None, x=-0.12, y=1.08, fontweight="bold"):
    """给每个 panel 加 a/b/c 粗体标签。axes 可以是 Axes 数组或列表。"""
    import numpy as np

    axs = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
    if labels is None:
        labels = [chr(ord("a") + i) for i in range(len(axs))]
    for ax, lab in zip(axs, labels):
        ax.text(x, y, lab, transform=ax.transAxes,
                fontsize=FONT_SIZES["panel_label"], fontweight=fontweight,
                va="top", ha="right")
    return axs
