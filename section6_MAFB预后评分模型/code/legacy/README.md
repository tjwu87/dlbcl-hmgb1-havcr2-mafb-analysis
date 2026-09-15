# code/legacy/ —— 归档件（不参与 run_all.py 流程）

这两个文件是作者当年的**原始母脚本**，从 `Tcell/` 按原样收录（仅追加了一段
"归档收录说明"头 + 把服务器绝对路径包进 `translate_path()`，正文逻辑未改）。

## 为什么收录

`section6` 的 `07_01_prognosis_main.R` 会**读**下面这些表，但 `07_01`~`07_05`
只读不产 —— 仓库原先没有任何脚本能产出它们。这两个文件就是产出者：

| 中间表 | `yuhou.R`（25 基因版） | `yuhou2.R`（23 基因版） |
|---|---|---|
| `LASSO_prognostic_formula.csv` | L289 | L80 |
| `GSE10846_risk_scores.csv` | L304 | L96 |
| `Cox_univariable_multivariable_results.csv` | **L957** | L411 |
| `Calibration_data_1_3_5yr.csv` | **L1571** | L528 |
| `Immune_infiltration_scores.csv` | **L1731** | — |
| `Immune_RiskScore_correlations.csv` | **L1734** | — |
| `All_cohorts_risk_scores.csv`（5 队列） | **L3356** | — |
| `MultiCohort_KM_statistics.csv` | **L3404**（5 队列） | L1995 起（6 队列） |
| `TCGA_DLBC_risk_scores.csv` | **L3414** | — |
| `GSE87371/11318/181063_risk_scores.csv` | **L3421/3426/3431** | — |

> **判定依据**：数据根 `DLBCL_prognosis/tables/` 里的
> `Cox_univariable_multivariable_results.csv` 列名为
> `Variable,label,type,HR,CI_low,CI_high,pval,pval_str,hr_str`，
> 与 `yuhou.R` L957 的 `select(Variable, label, type=section, HR, CI_low, CI_high,
> pval, pval_str, hr_str)` **逐字一致**；`All_cohorts_risk_scores.csv` 的 5 个队列
> （412/221/199/882/45）在 `yuhou.R` L3350-3351 明确列出。
> 文件顶部的 `LASSO_prognostic_formula.csv` 含 15 个入选基因，全部落在
> `07_02_lasso_cv_curves.R` 的 23 基因候选集内。

## 两者的关系

| | `yuhou.R` | `yuhou2.R` |
|---|---|---|
| 基因集 | 25 个（HMGB1+HAVCR2+23） | **23 个**（去掉 HMGB1/HAVCR2）＝ `07_02` 的候选集 |
| 输出目录 | `/mnt/results/DLBCL_prognosis/` → `data/DLBCL_prognosis/` | `/mnt/results/DLBCL_prognosis_v2/` → 同上（LEGACY_MAP 归并） |
| 队列 | **5 个**（GSE10846/87371/11318/181063/TCGA）＝ 论文所用 | 6 个（多一个 **NCICCR-DLBCL**） |
| 与论文关系 | **产出数据根里现有的那批表** | 平行支线；LASSO 部分已由 `07_02` 覆盖 |

## 可以跳过 / 非必需的段落（行号＝原文行号）

`yuhou.R`
- L225-268 —— 探索：改用全 420 样本、试 elastic net(α=0.5)、试纯 LASSO
- L374-442 / L455-541 / L542-720 —— CV 曲线与系数轨迹的**三版重绘尝试**，后者取代前者
  （对应图现在由 `07_01` 的 `A0*` 产出）
- L2064-2205 —— 变量/对象的一次性 reload 调试段
- L2897-2955 —— GSE87371 方向性排查（**但 L2992 的 `cens_os` 方向修正是必需的，别跳**）

`yuhou2.R`
- L685-1150 —— **NCICCR-DLBCL** 队列（论文未使用；含 GDC 下载与两套缩放方案对比）
- L1538-1690 —— TCGA 走 TCGAbiolinks 的实现（`yuhou.R` 走 GDC API，二选一即可）
- L2020-2394 —— 6 队列 KM 的多轮配色/排布重绘

## 运行提示

- ⚠️ **这两个文件是 R + 内嵌 Python 的混合工作文件**：正文里夹着
  `import requests` / `pd.read_csv(...)` 等 Python 片段（`yuhou.R` 约 L2371 起、
  `yuhou2.R` 约 L681 起是 TCGA/GDC 交互）。因此**整文件无法直接 `source()` 运行**
  ——用 `Rscript -e 'parse("<file>")'` 会在 Python 段报错（**原文即如此，不是收录时改坏的**）。
  当年应是分块手动执行。收录目的是**保留产出逻辑与参数**，不是提供一键入口。
- 它们**依赖上游产物**：`GSE10846_series_matrix.txt.gz`（GEO 可下，见 `data/README.md`）、
  `LASSO_prognostic_formula.csv` 等；不是"下载仓库就能直接跑"的入口。
- 原文用 `/workspace` 作暂存目录，已映射到 `<数据根>/external/workspace`，
  首次运行前建议先建好该目录。
- 依赖 R 包：`GEOquery`、`Biobase`、`hgu133plus2.db`、`illuminaHumanv4.db`、
  `glmnet`、`survival`、`survminer`、`rms`、`ggplot2`、`patchwork`、`cowplot`、
  `dplyr`、`stringr`、`httr`/`jsonlite`（TCGA GDC API）；`yuhou2.R` 另需
  `reticulate` 或手工分块跑 Python 段。

## 收录时的改动（可审计）

相对 `Tcell/` 下的原文，归档件**只**多了两处改动，已用行级 diff 逐行核对：

1. 文件头追加了一段「归档收录说明」注释（`yuhou.R` 29 行 / `yuhou2.R` 25 行）；
2. 把服务器绝对路径字面量包进 `translate_path()`（`yuhou.R` 80 处 / `yuhou2.R` 33 处）。

除此之外**没有改动任何一行正文逻辑**。
