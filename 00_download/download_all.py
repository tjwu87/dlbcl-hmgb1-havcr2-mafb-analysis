"""
数据获取
================================================================================
本仓库不包含任何原始数据或中间产物。运行本脚本只下载**公开原始数据**；
分析中间产物（adata_*.h5ad、pySCENIC/ CellOracle 产物、DLBCL_prognosis 的
tables/ 等，约 8 GB）**无法从 GEO 下载**，需由作者侧提供 —— 见仓库 README
「已知缺口」与 data/README.md。

    python 00_download/download_all.py --list          # 查看全部数据集
    python 00_download/download_all.py --dataset GSE182434
    python 00_download/download_all.py --all           # 全部下载（耗时且体积大）

下载位置由环境变量 DLBCL_DATA_ROOT 决定（默认 ./data）。
================================================================================
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from config.paths import DATA_ROOT  # noqa: E402

# ── 数据集清单 ────────────────────────────────────────────────────────────────
DATASETS = {
    # ── 单细胞（GEO）─────────────────────────────────────────────────────────
    "GSE182434": {
        "desc": "单细胞 RNA-seq 发现队列（4 例 DLBCL + 1 例扁桃体）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182434",
        "target": "GSE182434/",
        "size": "~2.3 GB",
        "used_in": "section1-2",
        "files": [
            "GSE182434_raw_count_matrix.txt.gz",   # ← GEO 提供（01_01 的输入）
            "GSE182434_cell_annotation.txt.gz",    # ← GEO 提供
            # 注：adata_processed.h5ad / adata_mac_*.h5ad 是**分析产物**，
            #     由 section1/code/01_01_qc_integration.py 等产出，GEO 没有。
        ],
    },
    # ── 空间转录组（GEO）─────────────────────────────────────────────────────
    "GSE232853": {
        "desc": "GeoMx DSP 空间转录组（DLBCL 组织芯片）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232853",
        "target": "GSE232853_v2/",
        "size": "~200 MB",
        "used_in": "section5",
        "files": [
            # 注意：论文用的是 Q3 重处理版（_v2）。原始 RLT 矩阵下载后需要
            # 先做过滤与 Q3 归一化 —— 脚本 06_01_spatial_analysis.py 内含
            # Q3 校正，但其上游 expression_raw_filtered.csv / meta_filtered.csv
            # 需另行准备（见 README 已知缺口 B.3）。
            "meta_filtered.csv",
        ],
    },
    # ── bulk 预后队列（GEO）—— 与 section6 代码实际使用的一致 ────────────────
    #   训练队列 + 4 个验证队列；05_90 之前代码里出现的 GSE32918/GSE4475/
    #   GSE23501 已核对为**非论文所用**，已从清单移除。
    "GSE10846": {
        "desc": "bulk 表达谱 —— 预后模型**训练队列**（GPL570；n=412 有生存信息）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846",
        "target": "DLBCL_prognosis/",
        "size": "~48 MB",
        "used_in": "section6",
        "files": [
            "GSE10846_series_matrix.txt.gz",       # ← 07_02/07_03 的输入
            "GSE10846_expression_gene_level.csv",  # 由 series matrix 整理得到
        ],
    },
    "GSE87371": {
        "desc": "bulk 表达谱 —— 验证队列（GPL570；n=221）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87371",
        "target": "DLBCL_prognosis/",
        "size": "~63 MB",
        "used_in": "section6",
        "files": ["GSE87371_series_matrix.txt.gz"],
    },
    "GSE11318": {
        "desc": "bulk 表达谱 —— 验证队列（GPL570；n=199）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11318",
        "target": "DLBCL_prognosis/",
        "size": "~20 MB",
        "used_in": "section6",
        "files": ["GSE11318_series_matrix.txt.gz"],
    },
    "GSE181063": {
        "desc": "bulk 表达谱 —— 验证队列（Illumina HT-12 v4 / GPL14951；n=882）",
        "source": "GEO",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181063",
        "target": "DLBCL_prognosis/",
        "size": "~99 MB",
        "used_in": "section6",
        "files": ["GSE181063_series_matrix.txt.gz"],
    },
    "TCGA-DLBC": {
        "desc": "TCGA 弥漫大 B 细胞淋巴瘤（验证队列，n=45）",
        "source": "GDC",
        "url": "https://portal.gdc.cancer.gov/projects/TCGA-DLBC",
        "target": "DLBCL_prognosis/TCGA/",
        "size": "~100 MB",
        "used_in": "section6",
        "files": [],
    },
    # ── 基因集（gseapy 自动下载）─────────────────────────────────────────────
    "MSigDB": {
        "desc": "Hallmark / GO / Reactome 基因集（gseapy 自动下载）",
        "source": "MSigDB",
        "url": "https://www.gsea-msigdb.org/gsea/msigdb/",
        "target": "reference/",
        "size": "~50 MB",
        "used_in": "section3,5,6",
        "files": [],
    },
    # ── 以下属已被论文删除的第七部分（HAVCR2 结构与成药性），本归档不含其代码 ──
    "GDSC2": {
        "desc": "药物敏感性参考数据（IC50 与细胞系表达谱）—— **第七部分用，本归档不含**",
        "source": "GDSC",
        "url": "https://www.cancerrxgene.org/downloads",
        "target": "GDSC2/",
        "size": "~260 MB",
        "used_in": "（不含）",
        "files": ["GDSC2_RNAseq_log2TPM_training_cells.csv"],
    },
}

# ── 参考文件（需手动获取，许可限制）────────────────────────────────────────────
MANUAL = {
    "CIBERSORT": {
        "desc": "免疫浸润反卷积 —— R 包（内含 extdata/LM22.txt）",
        "url": "https://ciberx.stanford.edu/",
        "note": "需注册；仅限学术研究使用。注意 section6 的 07_05_cibersort.R 是从"
                "**包内** extdata/LM22.txt 读（system.file(...)），不是从工作目录读；"
                "重算还需自备 exp.txt，见 data/README.md。",
    },
    "base_GRN_human_promoter.csv": {
        "desc": "CellOracle 启动子先验基因调控网络",
        "url": "https://github.com/morris-lab/CellOracle",
        "note": "从 CellOracle 仓库获取；放 data/GSE182434/celloracle_rerun_20260324_035651/",
    },
    "Macro_mono_cellmarker.xlsx": {
        "desc": "髓系 marker 表（Fig1a/1c 与髓系面板用）",
        "url": "—",
        "note": "作者侧参考表；放 data/ 下",
    },
    "GSE10846_series_matrix.txt.gz": {
        "desc": "GSE10846 series matrix（GEO 公开可下，非注册资源）",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846",
        "note": "脚本按 data/external/GSE10846survive/ 查找（LEGACY_MAP 映射）",
    },
}

# ── 无法下载、必须由作者提供的分析中间产物（约 8 GB）──────────────────────────
NOT_DOWNLOADABLE = {
    "GSE182434/adata_processed.h5ad 等": "section1/2 的起点（由 01_01 从原始矩阵重算）",
    "GSE182434/cellcomm/liana_results_*.csv": "LIANA 计算结果（计算本体不在仓库内）",
    "GSE182434/scenic/data_*.csv": "pySCENIC 产物（上游命令行不在仓库内）",
    "celloracle0331/* 与 GSE182434/celloracle_*": "CellOracle 产物（可由 section4 的 step* 重算）",
    "GSE232853_v2/expression_raw_filtered.csv 等": "空间数据上游预处理产物",
    "DLBCL_prognosis/tables/*.csv": "预后模型中间表（可由 section6/code/legacy/yuhou.R 产出）",
}


def list_datasets() -> None:
    print(f"{'名称':<14}{'体积':<10}{'用于':<12}说明")
    print("-" * 78)
    for name, meta in DATASETS.items():
        print(f"{name:<14}{meta['size']:<10}{meta['used_in']:<12}{meta['desc']}")
    print("\n需手动获取：")
    for name, meta in MANUAL.items():
        print(f"  {name}: {meta['url']}  ({meta['note']})")
    print("\n⚠ 本脚本只覆盖公开原始数据；下列分析中间产物**无法下载**，")
    print("  必须由作者侧提供（详见仓库 README「已知缺口」）：")
    for k, v in NOT_DOWNLOADABLE.items():
        print(f"  - {k}  —— {v}")


def download_geo(acc: str, meta: dict) -> bool:
    """用 GEOparse 下载 GEO 数据集。"""
    target = DATA_ROOT / meta["target"]
    target.mkdir(parents=True, exist_ok=True)
    print(f"[GEO] {acc} -> {target}")
    print(f"      {meta['url']}")
    try:
        import GEOparse
        GEOparse.get_GEO(acc, destdir=str(target))
        return True
    except ImportError:
        print("      未安装 GEOparse，改用浏览器/Aspera 手动下载：")
        print(f"      {meta['url']}")
        return False


def download_gdc(meta: dict) -> bool:
    target = DATA_ROOT / meta["target"]
    target.mkdir(parents=True, exist_ok=True)
    print(f"[GDC] TCGA-DLBC -> {target}")
    print("      推荐使用 TCGAbiolinks (R)：")
    print('      query <- GDCquery(project = "TCGA-DLBC",')
    print('                        data.category = "Transcriptome Profiling",')
    print('                        data.type = "Gene Expression Quantification")')
    print("      GDCdownload(query); data <- GDCprepare(query)")
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dataset")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    print(f"数据根目录: {DATA_ROOT}")

    if args.list or not (args.dataset or args.all):
        list_datasets()
        print(f"\n下载示例：python 00_download/download_all.py --dataset GSE182434")
        return

    names = list(DATASETS) if args.all else [args.dataset]
    for name in names:
        meta = DATASETS.get(name)
        if not meta:
            print(f"[warn] 未知数据集 {name}")
            continue
        if meta["source"] == "GEO":
            download_geo(name, meta)
        elif meta["source"] == "GDC":
            download_gdc(meta)
        else:
            print(f"[{meta['source']}] {name}: 请从 {meta['url']} 手动下载到 {meta['target']}")

    print("\n下载完成后确认目录结构与 config/paths.py 中的约定一致。")


if __name__ == "__main__":
    main()
