# code/04_pySCENIC/ —— 上游本体脚本（找回件）

`04_01_pyscenic_downstream_stats_figures.py` 只做**统计与出图**，它读的是
`<数据根>/GSE182434/scenic/data_*.csv`。这些 CSV 由 pySCENIC 本体产生，
原先仓库里没有产出它们的代码。下面三个脚本从 Biomni 平台的分析记录中找回，
参数与数据库地址均为**原文照录**。

| 脚本 | 作用 |
| --- | --- |
| `04_00a_setup_scenic_db.sh` | 下载 pySCENIC 三个参考数据库（TF 列表 / hg38 motif 排名 / motif-TF 注释），并在注释中给出官方 CLI 三步命令 |
| `04_00b_build_loom.py` | 构建 loom 输入（Mono+IFN_TAM+LA_TAM = 327 细胞 × 9118 基因，≥5% 表达过滤） |
| `04_00c_pyscenic_run_and_export.py` | GRNBoost2 → 共表达模块 → RcisTarget（NES≥3.0）→ AUCell，并导出下游所需 `data_*.csv` |

## 环境

`pyscenic`、`arboreto`、`ctxcore`、`loompy`，另需 `scanpy`/`anndata` 读 h5ad。

## 为什么用 Python API 而不是 CLI

Biomni 记录显示 CLI 的 dask 后端在容器里不稳定（`pyscenic ctx` 长时间输出 0 字节，
且 `aggregate_func` 与新版 dask 的 `from_delayed` 不兼容）。最终采用
`client_or_address='custom_multiprocessing'` 的 Python API 绕过 dask，
阈值参数与 CLI 完全一致。CLI 等价命令写在 `04_00a` 的文件头注释里。

## 前置条件

`04_00b` 需要已注释髓系亚型的 h5ad（含 `mac_subtype` 列，
Mono/IFN_TAM/LA_TAM 之外的细胞会被排除）。该注释属单细胞上游步骤，
不在本仓库范围内。
