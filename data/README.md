# 数据目录（占位，不随仓库分发）

本仓库**只包含代码**。原始数据全部来自公开数据库（GEO / GDC / GDSC / MSigDB），
体积约 8 GB，请自行下载后放到这里，或用环境变量指向已有位置：

```bash
# 方式一：放到本目录（推荐）
#   把各数据集解压/下载到 data/ 下，目录名与 run_all.py --check 的提示一致：
#   data/GSE182434/       单细胞发现队列
#   data/GSE232853_v2/    GeoMx 空间转录组（重处理版）
#   data/celloracle0331/  CellOracle 扰动产物
#   data/DLBCL_prognosis/ bulk 预后建模
#   data/GDSC2/           药物敏感性
#   data/spatial_analysis/

# 方式二：指向任意已有目录
export DLBCL_DATA_ROOT=/path/to/data      # Linux / macOS
set    DLBCL_DATA_ROOT=D:\bulk-download   # Windows

# 查看需要哪些数据
python 00_download/download_all.py --list
```

`config/paths.py` 的解析顺序：`DLBCL_DATA_ROOT` → `<仓库>/data/` → 仓库内 `data/`。

需注册后手动获取的参考文件（不便随包分发）：

| 文件 | 放到 | 说明 |
|---|---|---|
| `GSE10846_series_matrix.txt.gz` | `data/external/GSE10846survive/` | GEO 公开可下；`section6` 的 R 脚本按此路径查找 |
| `GSE10846_ciber.Rdata` | `data/external/Tcell/` | CIBERSORT 结果的**缓存**（有它就跑得通；删掉则需重算） |
| `immu_check_point.txt` | `data/external/Tcell/` | 免疫检查点基因列表（脚本按工作目录读） |
| `DLBCL_MRGs_Cluster_PD1.csv` | `data/external/Tcell/` | TIDE 输入表 |
| `exp.txt` | `data/external/Tcell/` | **仅重算 CIBERSORT 时必需**（基因×样本表达矩阵；脚本不生成） |
| `CIBERSORT` R 包 | R 库 | 提供 `cibersort()` 与内置 `extdata/LM22.txt`，需从 <https://ciberx.stanford.edu/> 注册获取（仅限学术使用） |
| `base_GRN_human_promoter.csv` | `data/GSE182434/celloracle_rerun_20260324_035651/` | CellOracle 先验 GRN，见 <https://github.com/morris-lab/CellOracle> |
| `Macro_mono_cellmarker.xlsx` | `data/` | 髓系 marker 表（Fig1a/1c 与髓系面板用） |

> `data/external/` 这两个子目录名对应 `config/paths.py` / `paths.R` 中 LEGACY_MAP
> 的 `external/Tcell`、`external/GSE10846survive` 两项——脚本里的历史绝对路径
> （`D:/BadiduNetdiskDownload/R/…`）会在运行时被翻译到这里，无需改代码。
