#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
01_05_fig1b_marker_dotplot.py —— 重跑 Fig1b（marker 点图）

为什么单独写：原 `01_01_qc_integration.py` 要读 GSE182434_raw_count_matrix.txt.gz
（58 MB gz，解压后极大），在本机反复段错误（exit=139）。
而 `adata_processed.h5ad`（17484 细胞 × 49632 基因，obs 含 cell_type）本地就有，
可直接做 Wilcoxon + sc.pl.dotplot，得到与原图一致的「每类 top4 标记基因」点图。

输出：marker_dotplot.svg / .pdf（配 retrofit 后画布=目标显示宽，字号≈6.5 pt）
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.paths import translate          # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

OUT_DIR = os.environ.get("DLBCL_RETRO_OUT") or str(
    Path(__file__).resolve().parents[1] / "results" / "fig1_panels" / "01b_marker_dotplot")
os.makedirs(OUT_DIR, exist_ok=True)

# 14 种细胞类型的展示顺序（与原图 1a/1b 图例一致）
CT_ORDER = [
    "CD8 T cells", "CD4 T cells", "Tregs", "TFH cells", "NK cells",
    "Naive B cells", "GC B cells", "Proliferative GC B cells",
    "Memory B cells", "Age-associated B cells", "B cells (other)",
    "Monocytes/Macrophages", "Plasma cells", "pDC/Other",
]

print("读取 adata_processed.h5ad …")
adata = sc.read_h5ad(translate(r"D:\bulk-download\GSE182434\adata_processed.h5ad"))
print(f"  {adata.shape}")

# 表达值优先用 raw（计数经 log1p 归一化后更适合做标记基因排序）
X = adata.raw[:, adata.var_names].X if adata.raw is not None else adata.X
import anndata as ad
adata = ad.AnnData(X=X, obs=adata.obs[["cell_type"]].copy(), var=adata.var[[]])
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].cat.set_categories(
    [c for c in CT_ORDER if c in adata.obs["cell_type"].cat.categories])
adata = adata[adata.obs["cell_type"].notna()].copy()

print("Wilcoxon rank_genes_groups（每类 vs 其余）…")
sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon",
                        pts=True, key_added="rg_celltypes")

markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rg_celltypes")
if "pvals_adj" in markers_df.columns:
    markers_df = markers_df[markers_df["pvals_adj"] < 0.01]
print(f"  显著基因 {len(markers_df)} 个")

# 每类取 top4，按 CT_ORDER 排列并去重 → 与原脚本逻辑一致
present = [c for c in CT_ORDER if c in markers_df["group"].unique()]
top_markers = (markers_df.groupby("group", sort=False).head(2)
               .groupby("group", sort=False)["names"].apply(list).to_dict())
genes_ordered = []
for ct in present:
    genes_ordered.extend(top_markers.get(ct, []))
genes_plot = [g for g in genes_ordered if not (g in genes_ordered[:genes_ordered.index(g)]
                                               or False)]
seen, genes_plot = set(), []
for g in genes_ordered:
    if g not in seen:
        seen.add(g)
        genes_plot.append(g)
print(f"Dot plot 基因数: {len(genes_plot)}")

# 原图 1B 是「宽扁」条状（2792x839 px，比例 3.3:1）：细胞类型在 Y、基因在 X。
# 必须显式给定 figsize，否则 scanpy 会按 52 个基因自算出一个 18 in 高的画布。
dp = sc.pl.dotplot(
    adata, var_names=genes_plot, groupby="cell_type",
    standard_scale="var", swap_axes=True, return_fig=True,
    figsize=(9.8, 7.6))
dp.style(cmap="Reds", dot_edge_color="#555555", dot_edge_lw=0.4)
# 竖排 90° 的细胞名过长会被画布底边裁掉（QA 实测 11 处出界）。
# 改 45° + 右对齐，让标签斜向伸入下方留白，同时缩个字号。
# 注意：DotPlot.savefig 每次都会调用 make_figure() **重建图形**（无幂等守卫），
# 任何 savefig 之前的刻度修改都会被抹掉。所以这里：先显式 make_figure()，
# 再改刻度，最后直接用 fig.savefig 保存（绕过 DotPlot.savefig）。
dp.make_figure()
for _ax in dp.get_axes().values():
    _xl = _ax.get_xticklabels()
    if len(_xl) >= 5 and any(t.get_text().strip() for t in _xl):
        _ax.tick_params(axis="x", labelsize=9)
        for tick in _xl:
            tick.set_ha("right")
            tick.set_rotation(45)
_fig = plt.gcf()
_fig.subplots_adjust(bottom=0.30)   # 给 45° 细胞名标签预留底部空间（否则被画布裁掉）
_fig.savefig(os.path.join(OUT_DIR, "marker_dotplot.svg"))
_fig.savefig(os.path.join(OUT_DIR, "marker_dotplot.pdf"))
plt.close("all")
print(f"✓ 已保存: {OUT_DIR}/marker_dotplot.(svg|pdf)")
