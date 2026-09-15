# Fig1 / Fig3 从零复现指南（新电脑）

本包 `fig1_fig3_repro_pack/` 包含在另一台电脑上从零重建 Fig1 与 Fig3 所需的
**全部代码、布局与静态资产**（`--with-data` 打包时还含输入数据）。

---

## 1. 环境要求

- **Python 3.10**（无需 R——Fig1/Fig3 全流程均不使用 R）
- Python 包（pip 安装）：
  ```
  pip install scanpy anndata h5py infercnvpy matplotlib numpy pandas scipy
  pip install seaborn statsmodels adjustText pypdf pypdfium2 reportlab Pillow openpyxl
  ```
- 数据盘：默认放在 **`D:/bulk-download/`**（`config/paths.py` 的 translate 按此根解析；
  如放别的盘，改 `config/paths.py` 里的根路径即可，全部脚本都经它取数）

## 2. 目录摆放

```
D:/bulk-download/                 ← 数据（--with-data 打包时在包内 bulk_download/ 下，拷到 D: 即可）
任意位置/
  └─ 仓库/                        ← 包内除 bulk_download 外的所有内容（保持相对结构）
      ├─ tools/  config/  results/assembled/fig1.json fig3.json ...
      ├─ docs/复现指南_fig1_fig3.md（本文件）
```

## 3. 一键运行

在**仓库根目录**下：

```bash
python tools/rebuild_fig1.py          # Fig1 全量（1a/1b/1c/1f 重渲 + 裁切 + 合成）
python tools/rebuild_fig1.py --fast   # 快速：复用已有 1a/1c 原始输出，只跑 1b/1f + 合成
python tools/rebuild_fig3.py          # Fig3 全量（3a+3e / 3b+3f / 3g 三节重渲 + 裁切 + 合成）
```

产出：
- `results/assembled/fig1/Fig1_composed.pdf`（174×214 mm）
- `results/assembled/fig3/Fig3_composed.pdf`（174×337 mm）

## 4. 面板 → 生成器对照（改内容时用）

| 面板 | 生成器 | 输入数据 |
|---|---|---|
| Fig1a | tools/10_00_fig1a_only.py（截断脚本） | GSE182434/adata_processed.h5ad、annotation/、cnv/、Macro_mono_cellmarker.xlsx |
| Fig1b | tools/fig1b_dotplot.py（手写点图） | adata_processed.h5ad、comparison/cluster_markers.csv |
| Fig1c | tools/10_00_fig1a_only.py（同一次运行） | 同 Fig1a（inferCNVpy） |
| Fig1d/1e | **静态资产**（svglib 转换件，不重渲） | — |
| Fig1f | tools/fig1f_volcano.py（02_02 火山图提取重渲） | GSE182434/cellcomm/volcano_data_….csv |
| Fig3a+3e | tools/04_01_sankey_only.py（强制副本） | GSE182434/scenic/*.csv、trajectory/*.csv |
| Fig3b+3f | tools/04_01_lineplot_only.py（强制副本） | 同上 |
| Fig3c/3d | **静态资产** | — |
| Fig3g | tools/suoxiao_dot_only.py | GSE182434/scenic/MAFB_gseapy_*.csv |

## 5. 布局修改

面板坐标、行距、字母位置全部固定在 `results/assembled/fig1.json` / `fig3.json`
（panels 显式坐标模式，x_mm/y_mm 自页顶/页左计）。**改布局=改 json，改内容=改面板脚本**。
合成器：`python tools/fig_compose_pdf.py results/assembled/fig1.json`。
注意 fig3 的 g 面板 x_mm 已含 -3mm 左移；f-g 行距 2.5mm。

## 6. 已知注意事项

- 重建脚本里的 `PY310` 变量若在本机不存在，会自动回退到运行脚本的解释器；
  请确保该解释器装有上列全部包。
- `vector_crop` 的 pad/th 参数是定稿值（1a pad=10/th=252 防顶部截断；浅灰文字需 th≥250）。
- 合成器 `fig_compose_pdf.py` 的 design 模式会跳过已显式调用 `subplots_adjust`
  的脚本的 `tight_layout`（保护补丁已内置于 config/retrofit.py）。
- 目标 PDF 被查看器占用时，合成器自动改写 `*_new.pdf`。
