"""
统一路径配置
================================================================================
本仓库中所有脚本都必须通过本模块获取路径，禁止再硬编码
    D:\\bulk-download\\...        (本地 Windows 画图环境)
    /mnt/results/...             (服务器 Linux 分析环境)
    /data/...                    (CellOracle 容器内路径)

用法
--------------------------------------------------------------------------------
    from config.paths import GSE182434, FIG_DIR, out_figure
    adata_path = GSE182434 / "adata_processed.h5ad"

切换数据位置（不改动任何代码）：
    Windows : set DLBCL_DATA_ROOT=D:\\bulk-download
    Linux   : export DLBCL_DATA_ROOT=/mnt/results
    Python  : os.environ["DLBCL_DATA_ROOT"] = "..."

本文件只做路径解析，不做任何 IO 副作用（不自动建目录），
以保证 import 时不会在用户磁盘上产生意外写入。
================================================================================
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Union

# ── 仓库根与数据根 ────────────────────────────────────────────────────────────
def _find_repo_root(start: Path) -> Path:
    """向上查找仓库根（同时含 `run_all.py` 与 `config/` 的目录）。

    兼容两种布局，避免依赖符号链接 / Windows junction：

        <repo>/config/paths.py                 ← 根配置（权威源）
        <repo>/sectionN/code/config/paths.py   ← 各 section 的副本，
                                                  便于单独运行该 section

    两种情形都会解析到同一个 <repo>，因此 DATA_ROOT / FIG_DIR 始终一致。
    """
    start = start.resolve()
    for cand in (start, *start.parents):
        if (cand / "run_all.py").is_file() and (cand / "config").is_dir():
            return cand
    # 兜底：保持与旧实现一致的行为
    return start.parents[1] if len(start.parents) > 1 else start


REPO_ROOT: Path = _find_repo_root(Path(__file__).parent)

# 数据根解析顺序：
#   1) 环境变量 DLBCL_DATA_ROOT（最高优先级，用于指向任意位置的大体积数据）
#   2) 归档内公用数据目录 00_共用/数据/（精简复制后的自包含数据根）
#   3) 仓库内 data/（仅放小体积参考数据与 .gitkeep）
DATA_ROOT: Path = Path(
    os.environ.get(
        "DLBCL_DATA_ROOT",
        (REPO_ROOT / "00_共用" / "数据")
        if (REPO_ROOT / "00_共用" / "数据").is_dir()
        else (REPO_ROOT / "data"),
    )
).expanduser()


def sub(name: str) -> Path:
    """返回数据根下的子目录（不创建）。"""
    return DATA_ROOT / name


# ── 各数据集目录（与论文结果部分一一对应）────────────────────────────────────
GSE182434: Path = sub("GSE182434")              # 单细胞发现队列 (scRNA-seq)
CELLORACLE: Path = sub("celloracle0331")        # CellOracle 扰动
GSE232853: Path = sub("GSE232853")              # GeoMx 空间转录组 (原始)
# 注：原始版未被论文使用，归档中不保留；若不存在则回退到重处理版 _v2，
# 避免任何残留引用静默读到空目录。
if not GSE232853.is_dir() and (DATA_ROOT / "GSE232853_v2").is_dir():
    GSE232853 = DATA_ROOT / "GSE232853_v2"
GSE232853_V2: Path = sub("GSE232853_v2")        # GeoMx 空间转录组 (重处理)
PROGNOSIS: Path = sub("DLBCL_prognosis")        # bulk 预后建模
GDSC2: Path = sub("GDSC2")                      # 药物敏感性
SPATIAL: Path = sub("spatial_analysis")         # 空间分析附加结果
CNV: Path = sub("CNV")                          # inferCNV

# ── 输出目录 ──────────────────────────────────────────────────────────────────
RESULTS: Path = Path(os.environ.get("DLBCL_OUT_ROOT", REPO_ROOT / "results"))
FIG_DIR: Path = RESULTS / "figures"
TAB_DIR: Path = RESULTS / "tables"

# ── 历史路径 → 新路径映射 ─────────────────────────────────────────────────────
# 用于把旧脚本里的绝对路径批量翻译成新结构，见 translate()
LEGACY_MAP = {
    # 本地 Windows 画图环境
    "d:/bulk-download/gse182434": "GSE182434",
    "d:/scrna/gse182434": "GSE182434",
    "d:/badidunetdiskdownload/r/ptl_m1_m2": "GSE182434",
    # 需注册 / 手工获取的外部文件（放数据根的 external/ 下，见 data/README.md）
    "d:/badidunetdiskdownload/r/gse10846survive": "external/GSE10846survive",
    "d:/badidunetdiskdownload/r/tcell": "external/Tcell",
    # 历史服务器暂存目录（原始母脚本用；见 section6/code/legacy/README.md）
    "/workspace": "external/workspace",
    "d:/bulk-download/celloracle0331": "CELLORACLE",
    "d:/bulk-download/gse232853_v2": "GSE232853_V2",
    "d:/bulk-download/gse232853": "GSE232853",
    "d:/bulk-download/dlbcl_prognosis": "PROGNOSIS",
    "d:/bulk-download/gse10846_expression_gene_level.csv": "GSE10846_expression_gene_level.csv",
    "d:/bulk-download": ".",
    # 服务器 Linux 分析环境
    "/mnt/results/dlbcl_prognosis_v2": "PROGNOSIS",
    "/mnt/results/dlbcl_prognosis": "PROGNOSIS",
    "/mnt/results/gse182434": "GSE182434",
    "/mnt/results": ".",
    # CellOracle 容器
    "/data": "CELLORACLE",
}

_NAMED = {
    "GSE182434": GSE182434,
    "CELLORACLE": CELLORACLE,
    "GSE232853_V2": GSE232853_V2,
    "GSE232853": GSE232853,
    "PROGNOSIS": PROGNOSIS,
    "CNV": CNV,
}


def translate(path: Union[str, Path]) -> Path:
    """把历史硬编码绝对路径翻译成当前数据根下的路径。

    未命中的路径原样返回（转为 Path），便于逐步迁移。

    特例：`adata_processed.h5ad` 用旧版 anndata 写出，其中
    `uns/log1p/base` 带 `encoding_type='null'`，新版 anndata（0.11+）没有对应
    读方法，读取即报 IORegistryError。修补版（删除该键）位于
    `results/data_fixed/`，若存在则自动改用，避免修改原始数据。
    """
    s = str(path).replace("\\", "/")
    low = s.lower().rstrip("/")
    for prefix, target in LEGACY_MAP.items():
        if low == prefix:
            base = _NAMED.get(target, DATA_ROOT / target)
            return _prefer_fixed(base)
        if low.startswith(prefix + "/"):
            rest = s[len(prefix) + 1:]
            base = _NAMED.get(target, DATA_ROOT / target)
            return _prefer_fixed(base / rest)
    return _prefer_fixed(Path(s))


_FIXED_DIR = REPO_ROOT / "results" / "data_fixed"
_FIXED_H5AD = _FIXED_DIR / "adata_processed.h5ad"


def _prefer_fixed(p: Path) -> Path:
    """若 results/data_fixed/ 下存在同名修补文件，优先使用。

    背景：这些 h5ad 由旧版 anndata 写出，`uns/log1p` 带 `encoding_type='null'`，
    新版 anndata（0.11+）读取即报 IORegistryError。修补版删除了该键，
    避免修改用户的原始数据。见 tools/fix_h5ad.py。
    """
    try:
        cand = _FIXED_DIR / p.name
        if cand.exists():
            return cand
    except (AttributeError, OSError, TypeError):
        pass
    return p


def ensure(*dirs: Union[str, Path]) -> None:
    """按需创建目录（仅在脚本显式调用时才建）。"""
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)


def legacy_roots() -> dict:
    """返回历史根目录清单，供 README 说明迁移方式。"""
    return {
        "windows_local": "D:/bulk-download",
        "linux_server": "/mnt/results",
        "celloracle_container": "/data",
    }
