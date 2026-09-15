# -*- coding: utf-8 -*-
"""矢量裁切公共函数：测内容 bbox -> 平移内容到原点 + 重设 mediabox（保持矢量）。

被 rebuild_fig1.py / rebuild_fig3.py 调用。
注意：
  - pad_pt 会四边各留白 pad_pt 磅；顶部内容贴边的图可加大 pad（如 1a 用 10）。
  - th 是"非白"灰度阈值；浅灰文字（灰阶 235-250）需用 250-252 才能检测到。
  - 必须同时平移内容并重设 mediabox（只改 mediabox 会导致内容错位/溢出——已踩坑）。
"""
import numpy as np
import pypdfium2 as pdfium
from pypdf import PdfReader, PdfWriter, Transformation


def vector_crop(src, dst, pad_pt=4, th=250):
    im = pdfium.PdfDocument(src)[0].render(scale=3.0).to_pil().convert('L')
    a = np.array(im)
    mask = a < th
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    reader = PdfReader(src)
    pg = reader.pages[0]
    W, H = float(pg.mediabox.width), float(pg.mediabox.height)
    ph, pw = im.height, im.width
    x0 = max(0, cols[0] / pw * W - pad_pt)
    x1 = min(W, cols[-1] / pw * W + pad_pt)
    y_top = rows[0] / ph * H
    y_bot = rows[-1] / ph * H
    y0 = max(0, H - y_bot - pad_pt)
    y1 = min(H, H - y_top + pad_pt)
    pg.add_transformation(Transformation().translate(tx=-x0, ty=-y0))
    pg.mediabox.lower_left = (0, 0)
    pg.mediabox.upper_right = (x1 - x0, y1 - y0)
    w = PdfWriter()
    w.add_page(pg)
    with open(dst, 'wb') as f:
        w.write(f)
    bb = PdfReader(dst).pages[0].mediabox
    print('[crop]', src, '->', dst,
          round(float(bb.width) / 72 * 25.4, 1), 'x',
          round(float(bb.height) / 72 * 25.4, 1), 'mm')
    return float(bb.width) / float(bb.height)
