# DLBCL HMGB1–HAVCR2–MAFB 多组学分析流程

弥漫大 B 细胞淋巴瘤（DLBCL）肿瘤微环境中 HMGB1、HAVCR2、MAFB 与 T 细胞功能关联的
**多组学分析代码**，按论文的 Results 小节组织：

| 部分 | 内容 | 图件 |
|---|---|---|
| `section1_恶性B细胞与髓系通讯` | 单细胞 QC/整合/注释、恶性判定、LIANA 细胞通讯 | Figure 1 (a–h)、S1 |
| `section2_髓系轨迹与LA_TAM终末态` | 髓系再分型、PAGA/CellRank 轨迹、LA-TAM 打分 | Figure 2 (a–f)、S2、S8 |
| `section3_HAVCR2-MAFB调控模块` | pySCENIC 调控子活性、RSS 受体优先级、伪时序相关 | Figure 3 (a–g)、S3、S9 |
| `section4_MAFB计算机扰动` | CellOracle：GRN → 扰动模拟 → 命运概率 → 敏感性 | Figure 4 (a–h) |
| `section5_空间转录组验证` | GeoMx：Q3 归一化 → 配对 ROI → DEA → GSEA → 空间共变 | Figure 5 (a–f)、S5 |
| `section6_MAFB预后评分模型` | bulk 多队列 LASSO-Cox、KM/校准、CIBERSORT 免疫浸润 | Figure 6 (a–g)、S6 + 补充表 |

> ⚠ **本版本只包含 Results 1–6 的代码。** 论文已删除
> 「HAVCR2 结构与成药性」（原 Figure 7 / S7：GDSC IC50、Boltz-2、分子对接、MD），
> 因此该部分代码不在本仓库内。
>
> ⚠ 本研究**全部为计算分析**，不含实验验证。文中"轴 / 驱动 / 靶向"等表述均为
> **假说性**结论，请按原文的 "candidate / associated with" 理解。

---

## 目录结构

```
.
├── run_all.py              # 一键复现入口（--stage 1..6 / --list / --check / --dry-run）
├── config/                 # 路径解析（paths.py/paths.R）、绘图规范、渲染包装、参数
├── section1..section6/     # 各 Results 小节：README.md + code/<模块>/
│   └── code/config/        # 各 section 的 config 副本（便于单独运行该 section；
│                           #   内容与根 config/ 一致）
├── tools/                  # 图面复现流水线（Fig1/Fig3 一键重建、通用合成、面板提取件）
├── docs/                   # panel ↔ 代码映射、版本审计、复现指南与状态报告
├── 00_download/            # 公开数据下载脚本
├── environment.yml         # conda 主环境（+ environment_celloracle/pyscenic.yml）
├── install_R_packages.R    # R 依赖
└── data/                   # 数据目录占位（不随仓库分发，见 data/README.md）
```

---

## 快速开始

```bash
# 1) 环境
conda env create -f environment.yml && conda activate dlbcl
Rscript install_R_packages.R
#   CellOracle / pySCENIC 依赖较苛刻，建议单独建环境：
#   conda env create -f environment_celloracle.yml
#   conda env create -f environment_pyscenic.yml

# 2) 数据（约 8 GB，不含在本仓库内）—— 见 data/README.md
export DLBCL_DATA_ROOT=/path/to/data      # Windows: set DLBCL_DATA_ROOT=D:\bulk-download
python 00_download/download_all.py --list
python 00_download/download_all.py --all

# 3) 一键复现
python run_all.py                 # 全部（section1 → section6）
python run_all.py --stage 4       # 只跑 CellOracle 扰动
python run_all.py --list          # 列出各节脚本与就绪状态
python run_all.py --check         # 只做输入数据预检
python run_all.py --dry-run       # 只打印将要执行的命令
```

`DLBCL_DATA_ROOT` 是唯一需要配置的路径。所有脚本通过 `config/paths.py`
（R 侧 `config/paths.R`）解析路径；代码中保留的历史绝对路径
（`D:\bulk-download\…`、`/mnt/results/…`）由 `config.paths.translate()`
在运行时映射到当前数据根，**无需手工改脚本**。

### ⚠ 已知环境坑：旧版 h5ad 读不了

随项目分发的 h5ad 由旧版 anndata 写出，`uns/log1p` 带 `encoding_type='null'`；
anndata 0.11+ 取消了对应读方法，会直接报
`IORegistryError: No read method registered for IOSpec(encoding_type='null', ...)`。
执行 `python tools/fix_h5ad.py` 即可修复：它把 h5ad **复制**到 `results/data_fixed/`
并删掉该键（**原始数据不动**），`config/paths.py` 会自动优先使用修补版。

---

## 从代码到图

主图由 panel 级脚本产出后拼装。一键重建与合成：

```bash
python tools/rebuild_fig1.py --fast                 # Figure 1
python tools/rebuild_fig3.py                        # Figure 3
python tools/fig_compose_pdf.py results/assembled/figN.json    # N = 2,4,5,6
```

| 图 | 主要绘图脚本 | 数据来源 |
|---|---|---|
| Fig 1 | `section1/code/10_figures/10_00_umap_cnv_myeloid_panels.py` + `tools/rebuild_fig1.py` | `section1/code/01_*`, `02_*` |
| Fig 2 | `section2/code/03_trajectory_PAGA/03_01_paga_trajectory.py` | `section2/code/01_*` |
| Fig 3 | `section3/code/04_pySCENIC/04_01_pyscenic_downstream_stats_figures.py` | pySCENIC 产物 CSV |
| Fig 4 | `section4/code/05_CellOracle/05_90_reproduce_figures_local.py` | `section4/code/step*` |
| Fig 5 | `section5/code/06_spatial_GSE232853/06_01_spatial_analysis.py` | `data/GSE232853_v2/` |
| Fig 6 | `section6/code/07_bulk_prognosis/07_01_prognosis_main.R`、`07_02_lasso_cv_curves.R` | `data/DLBCL_prognosis/` |

panel 与最终主图的对应关系见 `docs/figure_panel_map.md`、
`docs/figure_code_map_final.md`；各 section 的复现状态见
`docs/复现状态报告_20260913.md`。

---

## 出图规范

期刊版面 `\textwidth = 160 mm`、`\textheight = 216 mm`。
旧脚本把画布开到 13×13 in 甚至 20×6.5 in，字号却用默认 10 pt，
排版缩放后图上文字只剩 1–2 pt —— 这才是"图模糊"的真正原因，与像素数无关。

现在统一在 `config/plot_style.py`（R 侧 `config/plot_style.R`）中约定：

```python
from config.plot_style import apply_paper_style, panel_figsize, save_figure

apply_paper_style()                                    # 字号 8 pt 起，最小 6 pt
fig, axes = plt.subplots(2, 3, figsize=panel_figsize(3, 2))
save_figure(fig, FIG_DIR / "Fig1")                     # 输出 PDF(矢量) + PNG(600 dpi)
```

**矢量 PDF 优先**，放大任意倍数不失真。

对既有脚本的批量改造（不改原脚本，用改造层运行）：

```bash
python tools/run_with_retrofit.py <脚本.py> <显示宽in> <字号pt>
Rscript tools/run_with_retrofit.R <脚本.R> <显示宽in> <字号pt>
```

原理见 `config/retrofit.py` / `config/retrofit.R`：拦截 `plt.subplots` / `ggsave`
的画布尺寸并同步放大字号，输出统一重定向到 `results/`，**不覆盖原始图片**。

---

## 分析参数在哪里

分析代码与超参数**内嵌在管线脚本中**（非集中式配置），直接打开即可看到：
QC/双联去噪（scrublet 阈值）、归一化与 HVG（`n_top_genes`）、整合（harmonypy）、
聚类（leiden resolution）、PAGA/DPT 参数、pySCENIC 下游统计、CellOracle 逐步参数
（step00–step09：GRN 过滤、扰动模拟、敏感性网格）、空间分析（Q3 校正）、
LASSO/Cox（λ 网格、交叉验证折数）。

随机性：`seed = 42`（见 `config/params.yaml`）。审稿人关注的阈值
（LIANA `magnitude_rank < 0.05`、GRN 边数、HVG 数量等）集中在 `config/params.yaml`，
并附敏感性分析脚本（`step08*`）。

---

## 已知缺口（诚实说明）

> 本仓库**只含代码**。"下载仓库 → 直接出图"是**不成立**的：
> 脚本读的大部分输入是**分析中间产物**（约 8 GB），不在任何公开下载源里。
> 下面把缺口按性质列清（审计于 2026-09-14，逐条可追溯到文件路径）。

### A. 需要你自备的数据（不在仓库）

| 类别 | 内容 | 获取方式 |
|---|---|---|
| 原始公开数据 | GSE182434 原始 count 矩阵、GSE10846/GSE87371/GSE11318/GSE181063 series matrix、GeoMx GSE232853 原始 RLT 矩阵、GDSC2 训练矩阵 | `python 00_download/download_all.py --list`（该脚本给出的队列清单**与脚本实际使用的不完全一致**，以各 section 代码中出现的 `GSE*` 为准） |
| 分析中间产物 | `adata_processed.h5ad`、`adata_mac_*`、pySCENIC 产物、`celloracle0331/`、`GSE232853_v2/` 的 `expression_raw_filtered.csv`/`Q3_normalization_stats.csv`/`task2_*`、`DLBCL_prognosis/tables/*.csv`、`GSE10846_sur_model.Rdata` 等 | **需作者侧提供**（无法从 GEO 下载） |
| 需注册获取 | CIBERSORT 的 `LM22.txt`、CellOracle 的 `base_GRN_human_promoter.csv`、`Macro_mono_cellmarker.xlsx` | 见 `data/README.md` |

### B. 上游脚本未收录（原 4 处，**已全部补齐**）

1. ~~**pySCENIC 本体**~~ **已补齐**：原 `section3` 只有下游统计脚本
   `04_01_pyscenic_downstream_stats_figures.py`。现已从 Biomni 分析记录中找回本体，
   收录为 `section3/code/04_pySCENIC/` 下的三个文件：
   `04_00a_setup_scenic_db.sh`（参考数据库下载 + 官方 CLI 三步命令）、
   `04_00b_build_loom.py`（构建 327 细胞 × 9118 基因的 loom）、
   `04_00c_pyscenic_run_and_export.py`（GRNBoost2 → RcisTarget NES≥3.0 → AUCell，
   并导出下游所需 `data_*.csv`）。参数、数据库 URL、随机种子均为**原文照录**。
   详见 `section3/code/04_pySCENIC/README_recovered.md`。
2. ~~**LIANA 计算本体**~~ **已补齐**：`section1/code/02_cellcomm_LIANA/02_00_liana_run.py`
   名为 "run" 但实际是**绘图脚本**（读三张结果表）。现补入两个本体脚本：
   `02_00a_liana_rank_aggregate_14celltypes.py`（全细胞 14 类型）与
   `02_00b_liana_rank_aggregate_13subtypes.py`（13 亚型，含 Malignant/Normal B 拆分规则）。
   参数 `expr_prop=0.1 / min_cells=5 / n_perms=100 / seed=42` 为原文照录。
   详见 `section1/code/02_cellcomm_LIANA/README_recovered.md`。
3. ~~**空间转录组的上游处理**~~ **已补齐**：`section5` 的
   `06_01_spatial_analysis.py` 读的 11 张 CSV 现已有产出者 ——
   `section5/code/06_spatial_GSE232853/` 下补入 5 个脚本：
   `06_00a`（下载 GEO 原始矩阵 + 样本元数据 + 严格 QC）、
   `06_00b`（GeoMx Q3 归一化 + Harmony 批次校正 + PCA/UMAP）、
   `06_00c`（两套 DEA 口径 + 配对 ROI 表 + LA_TAM 分组表）、
   `06_00d`（GSEA prerank + ssGSEA）、
   `06_00e`（MAFB 靶基因打分 + 空间共变主表）。
   原始数据可从 GEO 公开下载（`GSE232853_Processed_data_CD20_CD68_final.csv.gz`）。
   详见 `section5/code/06_spatial_GSE232853/README_recovered.md`。
4. ~~**预后模型的若干中间表**~~ **已补齐**：`07_01_prognosis_main.R` 读的
   `Cox_univariable_multivariable_results.csv`、`Calibration_data_1_3_5yr.csv`、
   `Immune_infiltration_scores.csv`、`All_cohorts_risk_scores.csv`、
   `GSE87371/GSE11318/GSE181063/TCGA_*_risk_scores.csv` 等，
   其产出者是作者原始母脚本，现已收录到
   `section6/code/legacy/yuhou.R`（25 基因版，**产出数据根里现有的那批表**）
   与 `section6/code/legacy/yuhou2.R`（23 基因版平行支线，含论文未用的 NCICCR 队列）。
   逐表 ↔ 行号对照、以及"可跳过的探索段"清单见 `section6/code/legacy/README.md`。
   `07_02_lasso_cv_curves.R` 另产出 `LASSO_prognostic_formula.csv` 与 `GSE10846_risk_scores.csv`。

### C. 代码层面的已知问题（已在归档副本中修好或记录）

1. ✅ **`07_05_cibersort.R` 原缺路径迁移**（源项目里没有 `source(config/paths.R)`，
   3 处裸绝对路径）→ 归档副本已补注入并把路径包进 `translate_path()`。
2. ✅ **`config/paths.py` / `paths.R` 的 `LEGACY_MAP` 缺 2 个前缀**
   （`D:/BadiduNetdiskDownload/R/Tcell`、`…/R/GSE10846survive`）→ 已补映射到
   `data/external/Tcell` 与 `data/external/GSE10846survive`。
3. ⚠️ **CIBERSORT 重算需要 `exp.txt`**：`07_05` 第 117 行
   `cibersort(lm22f, "exp.txt", perm = 1000, QN = T)` 里的 `exp.txt` **脚本不生成、
   仓库与数据里也没有**；它有 `if (!file.exists("GSE10846_ciber.Rdata"))` 短路——
   即只有缓存了 `GSE10846_ciber.Rdata` 时才跑得通。要真正重算，需自备 `exp.txt`
   （基因×样本的表达矩阵）放到 `data/external/Tcell/`。
   ⚠️ 逆向地看：**CIBERSORT 的 `LM22.txt` 是从 R 包 `CIBERSORT` 的 `extdata/` 读的**
   （`system.file("extdata","LM22.txt", package="CIBERSORT")`），不是从工作目录读，
   所以需要先安装该 R 包（`install_R_packages.R` 无法自动安装，需手工获取）。
4. ⚠️ **`00_download/download_all.py` 的队列清单与实际不符**：脚本里列的是
   GSE32918 / GSE4475 / GSE23501，而 `section6` 代码实际用的是
   GSE10846 / GSE87371 / GSE11318 / GSE181063 / TCGA-DLBC。**以代码为准**。
   该脚本还错标 GEO 提供 `adata_processed.h5ad`（实为作者处理产物）。
5. ⚠️ **`/tmp/…` 交互路径**：`section4` 的 `step*` 用 `/tmp/oracle_after_02B.pkl` 等
   做步骤间暂存，在非 Unix 环境（Windows）会解析成 `	mp`。Linux/macOS 下正常。
6. ⚠️ **`config/paths.py` 在源码项目里依赖 Windows junction**（`section*/code/config`
   是指向根 `config/` 的目录联接）。归档已实体化并改为"向上自定位仓库根"，
   便于跨平台与 git 管理。

### D. 其他

- `docs/` 里的审计与映射文档记录了整理过程中的版本差异（例如某些图的面板源图
  来自素材库的早期版本），仅供追溯，不影响流程运行。
- Fig4 的面板源图由 matplotlib 直接输出矢量 PDF；若改用其它平台渲染同一份 SVG，
  需注意半透明填充（fill-opacity）可能丢失。

---

## 数据可用性

原始数据全部来自公开数据库（GEO / GDC / GDSC / MSigDB），
清单见 `python 00_download/download_all.py --list`。
本仓库**不包含**任何原始数据或大体积中间产物。

---

## 引用与许可

- 代码许可：**MIT**（见 `LICENSE`）
- 引用信息：见 `CITATION.cff`
- 若使用本代码，请引用原论文（信息待补充）。
