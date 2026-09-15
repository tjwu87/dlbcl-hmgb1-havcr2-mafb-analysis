# -*- coding: utf-8 -*-
"""矢量 PDF 合成器（figure-composer 规范，作者标准）：
  - 各面板 PDF 用 pypdf 矢量合并（缩放/平移不失真，文字仍可选可选中）
  - 两种排版模式：
    (1) rows 模式：每行面板等行高、按纵横比铺满行宽（零裁切、零竖向留白）
    (2) panels 模式（显式矩形）：config["panels"] = [{"letter","pdf","x_mm","y_mm","w_mm"}]
        —— y_mm 自页顶计，高由纵横比推得；用于复刻原图不等宽双列版式
  - 面板字母 a/b/c... 用 reportlab 矢量叠加（16pt 加粗、透明背景、面板左侧专用槽）
  - 合成后自动逐面板边缘复检：对比"合成前 vs 合成后"四边贴边率，变差即报警
用法: python tools/fig_compose_pdf.py <config.json>
"""
import sys, json, io
from pathlib import Path
import matplotlib, os

MM = 25.4
W_MM, GUT_MM = 174.0, 4.0
LET_MM, LGAP_MM = 5.0, 1.0   # 每个面板左侧的字母专用槽宽 + 槽后间隙
mm2pt = lambda v: v * 72 / MM

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas as rl_canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import pypdfium2 as pdfium

_font_cands = [
    os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf"),
    r"C:\Windows\Fonts\arialbd.ttf",
]
FT_PATH = next((f for f in _font_cands if os.path.exists(f)), None)
if FT_PATH:
    pdfmetrics.registerFont(TTFont("PanelBold", FT_PATH))
BOLD = "PanelBold" if FT_PATH else "Helvetica-Bold"

def edge_touch(img):
    """四边最外 3px 非白像素占比 (top,right,bottom,left)"""
    import numpy as np
    a = np.array(img.convert("L"))
    f = lambda s: float((s <= 245).mean())
    return (f(a[:3, :]), f(a[:, :3]), f(a[-3:, :]), f(a[:, -3:]))

def render_page(pdf_path, width_px=2000):
    pg = pdfium.PdfDocument(str(pdf_path))[0]
    return pg.render(scale=width_px / pg.get_size()[0]).to_pil().convert("RGB")

def _aspect(pdf):
    box = PdfReader(str(pdf)).pages[0].mediabox
    return float(box.width) / float(box.height)

def main(cfg_path):
    cfg = json.loads(Path(cfg_path).read_text(encoding="utf-8"))
    out = Path(cfg["out"]); out.parent.mkdir(parents=True, exist_ok=True)
    seq_letters = "abcdefghijklmnopqrstuvwxyz"
    # 页宽：优先用配置里的 page_mm[0]（布局脚本按面板宽度反算出来的）；
    #       配置没写才退回默认双栏整页宽 174mm。
    W = mm2pt(float(cfg["page_mm"][0]) if cfg.get("page_mm") else W_MM)
    LET = mm2pt(LET_MM); LGAP = mm2pt(LGAP_MM); GUT = mm2pt(GUT_MM)
    # 上下页边距（mm，可选）。整页左右本来有 LET/LGAP 边距，上下却是 0
    # （面板紧贴页顶/页底 → 顶部标题、底部 x 轴标题看起来像被裁掉）。
    # 配置里写 vpad_mm = N 即上下各留 N mm；不写 / 写 0 = 维持旧行为
    # （逐像素与历史产物完全一致）。
    VPAD = mm2pt(float(cfg.get("vpad_mm") or 0.0))

    # ---- 归一化为 placed = [(letter, pdf, x_pt, y_top_pt, w_pt, h_pt)]，y 自页顶 ----
    placed = []
    if "panels" in cfg:                       # 显式矩形模式
        for p in cfg["panels"]:
            asp = _aspect(p["pdf"])
            w = mm2pt(p["w_mm"]); h = mm2pt(p["h_mm"]) if p.get("h_mm") else w / asp
            x = mm2pt(p["x_mm"]); ytop = mm2pt(p["y_mm"]) + VPAD
            L = p.get("letter") or seq_letters[len(placed)]
            placed.append((L, p["pdf"], x, ytop, w, h))
        H_pt = max(y + h for _, _, _, y, _, h in placed) + VPAD
    else:                                     # rows 模式
        rows = cfg["rows"]
        row_h_pt = []
        for row in rows:
            avail = W - LET * len(row) - LGAP * (len(row) - 1)
            row_h_pt.append(avail / sum(_aspect(p[0]) if isinstance(p, (list, tuple)) else _aspect(p) for p in row))
        H_pt = sum(row_h_pt) + GUT * (len(rows) - 1) + 2 * VPAD
        ytop = VPAD
        idx = 0
        for ri, row in enumerate(rows):
            h = row_h_pt[ri]; x = 0.0
            for pdf in row:
                shift_mm = 0.0
                if isinstance(pdf, (list, tuple)):
                    pdf, shift_mm = pdf[0], float(pdf[1])
                L = seq_letters[idx]; idx += 1
                w = _aspect(pdf) * h
                placed.append((L, pdf, x + LET + mm2pt(shift_mm), ytop, w, h))
                x += LET + w + LGAP
            ytop += h + GUT

    writer = PdfWriter()
    page = writer.add_blank_page(width=W, height=H_pt)

    # 字母叠加层（reportlab 矢量）
    buf = io.BytesIO()
    c = rl_canvas.Canvas(buf, pagesize=(W, H_pt))
    for L, pdf, x, ytop, w, h in placed:
        src = PdfReader(str(pdf)).pages[0]
        sx = w / float(src.mediabox.width)
        sy = h / float(src.mediabox.height)   # ★ h_mm 时 sy!=sx（非等比）
        s = sx
        y = H_pt - ytop - h
        print(f"    merge {L}: src={Path(pdf).name} s={s:.3f} x={x:.0f} y={y:.0f} "
              f"src_wh=({float(src.mediabox.width):.0f},{float(src.mediabox.height):.0f})")
        page.merge_transformed_page(src, Transformation(ctm=(sx, 0, 0, sy, x, y)))
        # 字母：面板左侧专用槽内、顶对齐；透明背景，绝不压图
        c.setFillColorRGB(0, 0, 0)
        c.setFont(BOLD, 16)
        c.drawString(max(x - LET + mm2pt(0.6), mm2pt(0.6)), y + h - mm2pt(5.0), L)
    c.save()
    buf.seek(0)
    overlay = PdfReader(buf).pages[0]
    page.merge_page(overlay)          # 字母叠加（同尺寸页面，恒等变换）

    try:
        with open(out, "wb") as f:
            writer.write(f)
    except PermissionError:
        out = out.with_name(out.stem + "_new.pdf")
        with open(out, "wb") as f:
            writer.write(f)
        print(f"[warn] 目标被占用，改写 {out}")
    print(f"合成(矢量): {out}  {W/72*25.4:.0f}x{H_pt/72*25.4:.0f} mm")

    # ---- 合成后逐面板边缘复检：合成前 vs 合成后 ----
    import numpy as np
    print("边缘复检（合成前 -> 合成后，任一边恶化>0.10 即 ⚠️）:")
    problems = []
    doc = pdfium.PdfDocument(str(out))[0]
    full = doc.render(scale=1400 / W).to_pil().convert("RGB")
    for L, pdf, x, ytop, w, h in placed:
        y = H_pt - ytop - h
        before = edge_touch(render_page(pdf, max(300, int(w * 1400 / W))))
        px_x0 = int(x / W * full.size[0]); px_x1 = int((x + w) / W * full.size[0])
        px_y0 = int(ytop / H_pt * full.size[1]); px_y1 = int((ytop + h) / H_pt * full.size[1])
        after = edge_touch(full.crop((px_x0, px_y0, px_x1, px_y1)))
        worse = max(a2 - a1 for a1, a2 in zip(before, after))
        mark = " ⚠️" if worse > 0.10 else ""
        print(f"  {L}: 前{tuple(round(v,2) for v in before)} -> 后{tuple(round(v,2) for v in after)}{mark}")
        if worse > 0.10:
            problems.append((L, before, after))
    if problems:
        print(f"⚠️ {len(problems)} 个面板合成后边缘变差: {[p[0] for p in problems]}")
    else:
        print("✅ 全部面板合成前后边缘一致，无新增裁切")
    return problems

if __name__ == "__main__":
    main(sys.argv[1])
