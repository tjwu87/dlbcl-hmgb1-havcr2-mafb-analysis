# -*- coding: utf-8 -*-
"""PNG 面板 -> 单页 PDF（300dpi 定尺寸，可选裁切/自动去白边）
用法:
  python tools/png2pdf.py <in.png> <out.pdf> [--autocrop] [--crop x0,y0,x1,y1]
  --crop 用像素坐标（基于原图），先裁切再 autocrop
"""
import sys
import numpy as np
from PIL import Image
from pypdf import PdfWriter, PdfReader
import io

def autocrop(im, thresh=245):
    a = np.array(im.convert("L"))
    mask = a < thresh
    if not mask.any():
        return im
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    pad = 4
    y0 = max(0, rows[0] - pad); y1 = min(a.shape[0], rows[-1] + 1 + pad)
    x0 = max(0, cols[0] - pad); x1 = min(a.shape[1], cols[-1] + 1 + pad)
    return im.crop((x0, y0, x1, y1))

def main():
    src = sys.argv[1]; dst = sys.argv[2]
    im = Image.open(src).convert("RGB")
    if "--crop" in sys.argv:
        box = [int(v) for v in sys.argv[sys.argv.index("--crop") + 1].split(",")]
        im = im.crop(tuple(box))
    if "--autocrop" in sys.argv:
        im = autocrop(im)
    w_in, h_in = im.size[0] / 300.0, im.size[1] / 300.0
    buf = io.BytesIO()
    im.save(buf, format="PDF", resolution=300.0)
    buf.seek(0)
    w = PdfWriter(); w.append(buf)
    with open(dst, "wb") as f:
        w.write(f)
    print(f"{dst}: {im.size}px -> {w_in*25.4:.0f}x{h_in*25.4:.0f} mm")

if __name__ == "__main__":
    main()
