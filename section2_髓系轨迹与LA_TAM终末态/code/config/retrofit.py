"""
出图改造层（路线 B·第三版）
================================================================================
本版解决三件事
--------------------------------------------------------------------------------
1. **矢量 + 文字可提取**
   `svg.fonttype='none'`：SVG 里的文字保留为 <text> 而非路径，
   这样 (a) 缩放不失真，(b) 程序可以读出图上写了什么，用于建立
   「代码 ↔ 子图」的对应关系。同时 pdf.fonttype=42 保证 PDF 字体嵌入。

2. **字号锚定到目标值**
   矢量图缩放时字号会等比变化，所以「字号固定」意味着「图尺寸被锁定」。
   这里反推：若画布宽 W、当前字号 F，希望排版后字号达到 target_pt，
   则画布宽应为  W_need = F × 显示宽 / target_pt。
   画布宽超过 W_need 就按比例缩小画布（**字号不动**）——
   缩小画布会让文字相对图变大，从而使缩进版面后仍然清晰。

3. **只输出矢量**
   默认只写 SVG/PDF；PNG 仅在显式需要时生成（用于快速预览）。
   栅格图缩进版面后字形边缘必然发糊，不能作为投稿图。

用法
--------------------------------------------------------------------------------
    import config.retrofit as retrofit
    retrofit.install(target_display_in=6.3, target_font_pt=7.0)

环境变量
--------------------------------------------------------------------------------
    DLBCL_FIG_DISPLAY_W  版面可用宽度（inch），默认 6.3（= \\textwidth）
    DLBCL_FIG_FONT_PT    排版后期望字号（pt），默认 7.0
    DLBCL_FIG_VECTOR_ONLY  设为 1 时不写 PNG（默认 1）
================================================================================
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

# ★ 保护脚本显式的 subplots_adjust 不被 design 模式的 tight_layout 覆盖
import matplotlib.figure as _mpl_fig
_orig_subplots_adjust = _mpl_fig.Figure.subplots_adjust
def _subplots_adjust_guarded(self, *a, **k):
    self._user_adjusted = True
    return _orig_subplots_adjust(self, *a, **k)
matplotlib.figure.Figure.subplots_adjust = _subplots_adjust_guarded

from matplotlib import pyplot as plt

_DISPLAY_W = float(os.environ.get("DLBCL_FIG_DISPLAY_W", 6.3))
_FONT_PT = float(os.environ.get("DLBCL_FIG_FONT_PT", 7.0))
_VECTOR_ONLY = os.environ.get("DLBCL_FIG_VECTOR_ONLY", "1") == "1"

_out_dir = ""
_installed = False
_DESIGN = False            # True = 画布直接设为最终显示尺寸（1:1 排入版面）
_PLAIN = False             # True = 按脚本原始尺寸/字号原样出图，不做任何缩放
_STRIP = False             # True = 去掉描述性长标题（拼图时由图注承担）
_BASE: dict = {}           # rcParams 基线快照
_BASE_FONT = 10.0
_resized: list = []          # 记录被缩放的画布，供出图后核对


def install(target_display_in: float | None = None,
            target_font_pt: float | None = None,
            output_dir: str | None = None,
            vector_only: bool | None = None,
            design: bool = False,
            plain: bool = False,
            strip_titles: bool = False) -> None:
    global _DISPLAY_W, _FONT_PT, _out_dir, _installed, _VECTOR_ONLY, _DESIGN, _PLAIN, _STRIP
    _DESIGN = design
    _PLAIN = plain
    _STRIP = strip_titles
    if target_display_in:
        _DISPLAY_W = target_display_in
    if target_font_pt:
        _FONT_PT = target_font_pt
    if output_dir:
        _out_dir = output_dir
    if vector_only is not None:
        _VECTOR_ONLY = vector_only

    if _installed:
        return
    _installed = True

    rc = matplotlib.rcParams
    rc["svg.fonttype"] = "none"      # ★ 文字保留为文本，可提取、可缩放
    rc["pdf.fonttype"] = 42
    rc["ps.fonttype"] = 42

    # 快照基线：字号类与线宽类都要存，否则 'large'/'medium' 这类字符串
    # 会因缺少基线而被跳过，字号不受控（实测标题仍为 11 px）。
    for k in _FONT_KEYS + _SIZE_KEYS:
        _BASE.setdefault(k, rc.get(k))

    _orig_figure = plt.figure
    _orig_subplots = plt.subplots

    def patched_figure(*args, **kwargs):
        nf = _resize_figsize(kwargs.get("figsize"))
        if nf:
            kwargs["figsize"] = nf
        return _orig_figure(*args, **kwargs)

    def patched_subplots(*args, **kwargs):
        # 注意：不要在这里缩放！plt.subplots 内部会调用 plt.figure，
        # 而 plt.figure 也被打了补丁 → 同一个 figsize 会被缩放两次
        # （实测 0.65^2，图被多缩一轮）。缩放统一放在 patched_figure。
        return _orig_subplots(*args, **kwargs)

    plt.figure = patched_figure
    plt.subplots = patched_subplots
    matplotlib.pyplot.figure = patched_figure
    matplotlib.pyplot.subplots = patched_subplots

    _patch_savefig()
    _patch_h5ad_reader()


def required_canvas_width() -> float:
    """在当前字号下，画布宽最多应为多少，才能让排版后字号达到目标值。

    W_need = F × 显示宽 / 目标字号
    """
    f = float(matplotlib.rcParams.get("font.size", 10.0))
    return f * _DISPLAY_W / _FONT_PT


def _apply_design_scale(s: float) -> None:
    """设计模式：几何量按 s 缩放，字号单独提升到目标值。

    s = 最终显示宽 / 原始画布宽。字号若按 s 缩会掉到 3-4 pt（审稿人抱怨的原因），
    所以字号只按「目标字号 / 原始字号」缩放，而线宽、点径、刻度按几何比例缩放，
    这样线条粗细观感与原图一致，只有文字被放大到可读。
    """
    rc = matplotlib.rcParams
    fk = _FONT_PT / _BASE_FONT
    from matplotlib import font_manager as _fm
    for k in _FONT_KEYS:
        v = _BASE.get(k, rc.get(k))
        if v is None:
            continue
        # matplotlib 的默认字号多为字符串（'large' / 'medium' / 'xx-small'）。
        # 直接 float() 会抛异常被跳过，导致标题、轴标签、刻度字号**不受控**
        # （实测：画布缩到 3.12 in 后标题仍是 12 pt，把坐标区挤没）。
        # 这里按 font_scalings 折算成绝对 pt，再乘以字号系数。
        if isinstance(v, str):
            sc = _fm.font_scalings.get(v)
            if sc is None:
                continue
            rc[k] = sc * _BASE_FONT * fk
        else:
            try:
                rc[k] = float(v) * fk
            except (TypeError, ValueError):
                pass
    rc["font.size"] = _FONT_PT
    for k in _SIZE_KEYS:
        if _BASE.get(k) is not None:
            try:
                rc[k] = float(_BASE[k]) * s
            except (TypeError, ValueError):
                pass


_FONT_KEYS = ("font.size", "axes.labelsize", "axes.titlesize",
              "xtick.labelsize", "ytick.labelsize",
              "legend.fontsize", "legend.title_fontsize",
              "figure.titlesize")
_SIZE_KEYS = ("lines.linewidth", "axes.linewidth", "patch.linewidth",
              "grid.linewidth", "xtick.major.width", "ytick.major.width",
              "xtick.minor.width", "ytick.minor.width",
              "lines.markersize", "xtick.major.size", "ytick.major.size",
              "xtick.minor.size", "ytick.minor.size")


def _resize_figsize(figsize):
    """设计模式：画布直接设为最终显示尺寸；否则沿用「超宽才缩」。"""
    if not figsize:
        return None
    try:
        w, h = float(figsize[0]), float(figsize[1])
    except (TypeError, ValueError, IndexError):
        return None
    if _PLAIN:
        return None          # 原样：不改 figsize
    if _DESIGN:
        # 作者要求：**整张图跟文字按同一比例等比缩放**。
        # 比例只由字号决定 s = 目标字号 / 原字号，
        # 这样布局、线条、标记的相对关系与原设计完全一致，只是物理尺寸变小，
        # 文字、图例、坐标区之间绝不可能重叠。
        # （之前用 `_DISPLAY_W / w` 把画布宽度锁到版面宽，而字号用另一个比例，
        #  两者不一致 → 文字相对图形变大 → 图例/标题盖住坐标区。）
        s = _FONT_PT / _BASE_FONT
        _resized.append((round(w, 2), round(h, 2),
                         round(w * s, 2), round(h * s, 2)))
        return (w * s, h * s)
    need = required_canvas_width()
    if w <= need:
        return None
    k = need / w
    _resized.append((round(w, 2), round(h, 2), round(w * k, 2), round(h * k, 2)))
    return (w * k, h * k)


def _scale_figure_objects(fig, fk: float, s: float) -> None:
    """把图内所有文字/线条/标记按比例缩放，并重排布局。

    为什么必须这么做：许多脚本把 fontsize= / linewidth= **硬编码**在调用里
    （实测 06_01 有 94 处 fontsize=11/10/9…），rcParams 完全管不到。
    只有在 savefig 前直接改图形对象的属性，才能统一控制。
    fk = 目标字号 / 原设计字号；s = 几何缩放（画布缩放比）。
    """
    if getattr(fig, "_retrofit_scaled", False):
        return
    fig._retrofit_scaled = True

    def _fs(obj):
        try:
            v = obj.get_fontsize()
            if isinstance(v, (int, float)) and v > 0:
                obj.set_fontsize(v * fk)
        except Exception:                       # noqa: BLE001
            pass

    for ax in list(fig.axes):
        objs = [ax.title, ax.xaxis.label, ax.yaxis.label]
        objs += [getattr(ax, "_left_title", None), getattr(ax, "_right_title", None)]
        objs += list(ax.texts)
        objs += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        for o in objs:
            if o is not None:
                _fs(o)
        for leg in ([ax.get_legend()] if ax.get_legend() else []):
            for o in list(leg.get_texts()) + ([leg.get_title()] if leg.get_title() else []):
                _fs(o)
        # 线条 / 标记 / 边框 —— 按几何比例，保持与原设计一致的粗细观感
        for ln in ax.get_lines():
            try:
                ln.set_linewidth(ln.get_linewidth() * s)
                ln.set_markersize(ln.get_markersize() * s)
            except Exception:                   # noqa: BLE001
                pass
        for sp in ax.spines.values():
            try:
                sp.set_linewidth(sp.get_linewidth() * s)
            except Exception:                   # noqa: BLE001
                pass
        # 散点（PathCollection）的 s 参数是**面积**（pt^2）；
        # 要让点的直径按 s 缩放，面积需按 s^2 缩放。
        for coll in ax.collections:
            try:
                sz = coll.get_sizes()
                if len(sz):
                    coll.set_sizes(sz * (s ** 2))
            except Exception:                   # noqa: BLE001
                pass
        try:
            ax.tick_params(width=1.0 * s, length=3.5 * s)
        except Exception:                       # noqa: BLE001
            pass

    for o in list(fig.texts):
        _fs(o)
    if getattr(fig, "_suptitle", None) is not None:
        _fs(fig._suptitle)

    # 去掉描述性长标题：拼进主图后这些信息由图注（caption）承担，
    # 留着只会占用版面并造成 panel 之间拥挤。
    # 规则：图级 suptitle 一律去掉；axes 标题若含换行或长度 > 40 字符也去掉，
    # 短标题（如 1a 的 "Refined Cell Types" / "Tissue (DLBCL vs Tonsil)"）保留，
    # 因为它们是区分两个子坐标区的必要信息。
    if _STRIP:
        if getattr(fig, "_suptitle", None) is not None:
            fig._suptitle.set_text("")
        for ax in fig.axes:
            ti = ax.get_title()
            if ti and (("\n" in ti) or len(ti) > 40):
                ax.set_title("")

    # 字号变小后重新排布。**不要用 constrained**：脚本里有 fig.add_axes 手动放
    # colorbar 的图，constrained 会把主坐标区压扁（实测 joyplot 的 18 行全挤成一团）。
    # tight_layout 与脚本原本的用法一致，更安全。
    try:
        if not getattr(fig, '_user_adjusted', False):  # ★ 脚本用过 subplots_adjust 则不覆盖
            fig.tight_layout()
    except Exception:                           # noqa: BLE001
        pass



def _qa_bbox_check(fig, fname) -> None:
    """figure-style §9.1：保存前程序化检测文字/边框 bbox 碰撞（Agg 渲染器）。

    只报告「可见且非空」的文本框两两重叠、以及非刻度文本压到坐标区边框的情况。
    刻度标签贴着自己的轴不算问题。结果打到 stdout（进 runner 日志），并累计到
    _QA_VIOLATIONS 供批量结束后汇总。
    """
    try:
        import matplotlib.text as _mt
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        texts = [(t, t.get_window_extent(r)) for t in fig.findobj(_mt.Text)
                 if t.get_text().strip() and t.get_visible()]
        spines = [(s, s.get_window_extent(r)) for ax in fig.axes
                  for s in ax.spines.values() if s.get_visible()]
        tick = {id(ax): set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
                for ax in fig.axes}
        ov = [(a.get_text()[:30], b.get_text()[:30])
              for i, (a, ba) in enumerate(texts) for b, bb in texts[i+1:] if ba.overlaps(bb)]
        ov += [(t.get_text()[:30], "<spine>") for t, bt in texts for s, bs in spines
               if bt.overlaps(bs) and t not in tick.get(id(s.axes), set())]
        # 图例框压坐标区（skill 弱点补丁：之前只查文字互撞）
        for ax in fig.axes:
            leg = ax.get_legend()
            if leg is not None and leg.get_visible():
                lb = leg.get_window_extent(r)
                ab = ax.get_window_extent(r)
                if lb.overlaps(ab):
                    ov.append(("<legend>", "<axes>"))
        for leg in getattr(fig, "legends", []):
            lb = leg.get_window_extent(r)
            for ax in fig.axes:
                if lb.overlaps(ax.get_window_extent(r)):
                    ov.append(("<fig-legend>", "<axes>"))
        W, H = fig.canvas.get_width_height()
        out = [t.get_text()[:30] for t, bt in texts
               if bt.x0 < -2 or bt.y0 < -2 or bt.x1 > W + 2 or bt.y1 > H + 2]
        name = Path(str(fname)).name
        if ov or out:
            _QA_VIOLATIONS.append((name, ov, out))
            print(f"[QA] {name}: {len(ov)} 处重叠, {len(out)} 处出画布")
            for a, b in ov[:5]:
                print(f"     overlap: {a!r} <-> {b!r}")
            for o in out[:5]:
                print(f"     outside: {o!r}")
        else:
            print(f"[QA] {name}: OK")
    except Exception as e:                       # noqa: BLE001
        print(f"[QA] {name if 'name' in dir() else '?'} 检查失败: {str(e)[:80]}")

_QA_VIOLATIONS: list = []

def _patch_savefig() -> None:
    orig_fig = matplotlib.figure.Figure.savefig

    def patched_fig(self, fname, *args, **kwargs):
        fname = _redirect(fname)
        kwargs.pop("bbox_inches", None)
        if _DESIGN:
            _s = _FONT_PT / _BASE_FONT
            _scale_figure_objects(self, _s, _s)
        kwargs.setdefault("pad_inches", 0.0)
        if _VECTOR_ONLY and str(fname).lower().endswith(".png"):
            _qa_bbox_check(self, fname)
            return None                      # 矢量优先，跳过位图
        _qa_bbox_check(self, fname)
        r = orig_fig(self, fname, *args, **kwargs)
        # ★ svg→pdf 走 svglib 会丢 fill-opacity（透明底色变实心）和文字旋转，
        #   因此凡脚本写 svg 的同时，用 matplotlib 原生直出同名 pdf。
        if str(fname).lower().endswith(".svg"):
            pdf_name = str(fname)[:-4] + ".pdf"
            try:
                pdf_kwargs = dict(kwargs)
                pdf_kwargs["format"] = "pdf"   # ★ 脚本可能显式传 format='svg'，侧写必须覆盖
                orig_fig(self, pdf_name, *args, **pdf_kwargs)
            except Exception as e:
                print(f"[retrofit][warn] native PDF failed for {pdf_name}: {e}")
        return r

    matplotlib.figure.Figure.savefig = patched_fig

    orig_plt = plt.savefig

    def patched_plt(fname, *args, **kwargs):
        fname = _redirect(fname)
        kwargs.pop("bbox_inches", None)
        _f = plt.gcf()
        if _DESIGN:
            _s = _FONT_PT / _BASE_FONT
            _scale_figure_objects(_f, _s, _s)
        kwargs.setdefault("pad_inches", 0.0)
        if _VECTOR_ONLY and str(fname).lower().endswith(".png"):
            _qa_bbox_check(_f, fname)
            return None
        _qa_bbox_check(_f, fname)
        r = orig_plt(fname, *args, **kwargs)
        if str(fname).lower().endswith(".svg"):
            pdf_name = str(fname)[:-4] + ".pdf"
            try:
                pdf_kwargs = dict(kwargs)
                pdf_kwargs["format"] = "pdf"   # ★ 同 patched_fig
                orig_plt(pdf_name, *args, **pdf_kwargs)
            except Exception as e:
                print(f"[retrofit][warn] native PDF failed for {pdf_name}: {e}")
        return r

    plt.savefig = patched_plt
    matplotlib.pyplot.savefig = patched_plt


def _redirect(fname):
    if not _out_dir:
        return fname
    try:
        dest = Path(_out_dir) / Path(str(fname)).name
        dest.parent.mkdir(parents=True, exist_ok=True)
        return str(dest)
    except (TypeError, ValueError, OSError):
        return fname


def _patch_h5ad_reader() -> None:
    """旧版 h5ad 的 uns/log1p 带 null 编码，新版 anndata 读不了（见 tools/fix_h5ad.py）。"""
    fixed_dir = Path(__file__).resolve().parents[1] / "results" / "data_fixed"
    if not fixed_dir.is_dir():
        return
    try:
        import anndata as ad
    except ImportError:
        return

    orig = ad.read_h5ad

    def patched(filename, *args, **kwargs):
        try:
            cand = fixed_dir / Path(str(filename)).name
            if cand.exists():
                filename = str(cand)
        except (TypeError, ValueError, OSError):
            pass
        return orig(filename, *args, **kwargs)

    ad.read_h5ad = patched
    try:
        import scanpy as sc
        sc.read_h5ad = patched
    except ImportError:
        pass


def report_resized() -> str:
    """输出被缩放过的画布清单，便于核对。"""
    if not _resized:
        return "没有画布被缩放（都已在目标字号允许的尺寸内）"
    lines = [f"{w}x{h} → {nw}x{nh}" for w, h, nw, nh in _resized]
    return f"共缩放 {len(lines)} 个画布：\n  " + "\n  ".join(lines)
