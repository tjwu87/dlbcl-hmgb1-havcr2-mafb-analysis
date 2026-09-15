# -*- coding: utf-8 -*-
"""SVG -> 矢量 PDF 高保真转换（Chrome headless print-to-pdf）。

背景（2026-09-13 复现确认）：
    svglib 转 matplotlib SVG 会丢 fill-opacity / stroke-opacity（小提琴半透明填充、
    低 alpha 连线整体消失），还会丢粗体、错位背景 patch —— 完全不可用。
    本模块改用 Chrome 的 SVG 引擎渲染再打印为 PDF，保真度与浏览器所见一致，
    产出为矢量（Form XObject + 文字路径），源 SVG 里 rasterized 的散点保留为位图。

用法:
    python tools/svg2pdf_chrome.py <out_dir> <svg1> [svg2 ...]
    python tools/svg2pdf_chrome.py <out_dir> --list list.txt

要点:
    - SVG 用 @page size = 原始 pt 尺寸 + margin 0 内联进包装 HTML，保证 1:1。
    - Chrome 输出尺寸存在约 0.1~0.2% 的内部舍入，转完用 pypdf 做一次
      ctm 微缩放把页面精确校正回 SVG 原始尺寸（避免非等比形变）。
"""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation

CHROME_CANDIDATES = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]


def find_chrome():
    for p in CHROME_CANDIDATES:
        if os.path.exists(p):
            return p
    raise RuntimeError("未找到 Chrome / Edge 可执行文件")


def svg_size(svg_text):
    """取 SVG 根元素的 width/height（pt）。"""
    m = re.search(r"<svg[^>]*>", svg_text, re.S)
    if not m:
        raise ValueError("找不到 <svg> 根元素")
    tag = m.group(0)
    w = re.search(r'\bwidth="([\d.]+)(pt)?"', tag)
    h = re.search(r'\bheight="([\d.]+)(pt)?"', tag)
    if w and h:
        return float(w.group(1)), float(h.group(1))
    vb = re.search(r'viewBox="([\d.\s-]+)"', tag)
    if vb:
        v = vb.group(1).split()
        return float(v[2]), float(v[3])
    raise ValueError("SVG 无 width/height/viewBox")


def convert(svg_path, out_pdf, chrome=None):
    """单个 SVG -> 矢量 PDF，返回 (w_pt, h_pt)。"""
    svg_path = Path(svg_path)
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    chrome = chrome or find_chrome()

    raw = svg_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"(<svg.*?</svg>)", raw, re.S)
    if not m:
        raise ValueError(f"{svg_path.name}: 找不到完整 <svg>...</svg>")
    body = m.group(1)
    W, H = svg_size(raw)

    html = (
        '<!DOCTYPE html><html><head><meta charset="utf-8"><style>\n'
        f"@page {{ size: {W}pt {H}pt; margin: 0; }}\n"
        "html,body{margin:0;padding:0;background:#fff;}\n"
        "svg{display:block;}\n"
        "</style></head><body>" + body + "</body></html>"
    )

    with tempfile.TemporaryDirectory() as td:
        html_path = Path(td) / (svg_path.stem + ".html")
        html_path.write_text(html, encoding="utf-8")
        tmp_pdf = Path(td) / (svg_path.stem + ".pdf")
        cmd = [
            chrome, "--headless=new", "--disable-gpu", "--no-sandbox",
            "--no-pdf-header-footer", "--run-all-compositor-stages-before-draw",
            "--virtual-time-budget=5000",
            f"--print-to-pdf={tmp_pdf}",
            html_path.as_uri(),
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if not tmp_pdf.exists():
            raise RuntimeError(f"Chrome 未产出 PDF: {r.stderr[-500:]}")
        tmp_copy = out_pdf.parent / ("_tmp_" + out_pdf.name)
        tmp_copy.write_bytes(tmp_pdf.read_bytes())

    # Chrome 页面尺寸校正 -> SVG 原始 pt（消除 0.1~0.2% 内部舍入）
    reader = PdfReader(str(tmp_copy))
    pg = reader.pages[0]
    cw, ch = float(pg.mediabox.width), float(pg.mediabox.height)
    if abs(cw - W) > 0.01 or abs(ch - H) > 0.01:
        sx, sy = W / cw, H / ch
        pg.add_transformation(Transformation(ctm=(sx, 0, 0, sy, 0, 0)))
        pg.mediabox.lower_left = (0, 0)
        pg.mediabox.upper_right = (W, H)
    w = PdfWriter()
    w.add_page(pg)
    with open(out_pdf, "wb") as f:
        w.write(f)
    tmp_copy.unlink(missing_ok=True)
    return W, H


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    out_dir = Path(sys.argv[1])
    args = sys.argv[2:]
    svgs = []
    if args and args[0] == "--list":
        svgs = [ln.strip() for ln in Path(args[1]).read_text(encoding="utf-8").splitlines() if ln.strip()]
    else:
        svgs = args
    chrome = find_chrome()
    print(f"[chrome] {chrome}")
    for s in svgs:
        s = Path(s)
        dst = out_dir / (s.stem + ".pdf")
        W, H = convert(s, dst, chrome)
        r = PdfReader(str(dst))
        b = r.pages[0].mediabox
        print(f"[ok] {s.name:52s} {W:7.2f}x{H:7.2f} -> {float(b.width):7.2f}x{float(b.height):7.2f} pt  {dst.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
