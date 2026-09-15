"""
修补旧版 h5ad，使其可被新版 anndata 读取
================================================================================
问题
--------------------------------------------------------------------------------
这些 h5ad 由旧版 anndata 写出，其中 `uns/log1p` 带有 `encoding_type='null'`。
新版 anndata（0.11+）取消了对应的读方法，读取时直接报：

    IORegistryError: No read method registered for IOSpec(
        encoding_type='null', encoding_version='0.1.0')

后果：PAGA 轨迹（Fig2）与单细胞主面板（Fig1）脚本全部无法运行。

做法
--------------------------------------------------------------------------------
把 h5ad **复制**到 `results/data_fixed/`（不改动原始数据），删除 `uns/log1p`——
该键只记录 log1p 变换参数，不参与任何后续绘图。复制件通过
`config/paths.py` 的 `_prefer_fixed()` 自动生效，脚本无需改动。

用法
--------------------------------------------------------------------------------
    python tools/fix_h5ad.py                 # 修补默认列表
    python tools/fix_h5ad.py <file.h5ad>     # 修补指定文件
================================================================================
"""

from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

import h5py

REPO = Path(__file__).resolve().parents[1]
FIXED_DIR = REPO / "results" / "data_fixed"

DEFAULT_TARGETS = [
    "GSE182434/adata_processed.h5ad",
    "GSE182434/adata_mac.h5ad",
    "GSE182434/adata_mac_subtyped.h5ad",
]


def needs_fix(path: Path) -> bool:
    """检查 uns/log1p 是否存在且带 null 编码。"""
    try:
        with h5py.File(path, "r") as h:
            if "uns" not in h or "log1p" not in h["uns"]:
                return False
            g = h["uns"]["log1p"]
            et = g.attrs.get("encoding-type")
            if isinstance(et, bytes):
                et = et.decode()
            # 空 group 或显式 null 都需要处理
            return et == "null" or len(g) == 0
    except Exception:
        return False


def fix(src: Path, dst: Path) -> bool:
    t0 = time.time()
    dst.parent.mkdir(parents=True, exist_ok=True)
    print(f"复制 {src.name} ({src.stat().st_size / 1048576:.0f} MB) …")
    shutil.copy2(src, dst)
    with h5py.File(dst, "a") as h:
        if "uns" in h and "log1p" in h["uns"]:
            del h["uns"]["log1p"]
    print(f"  已删除 uns/log1p，用时 {time.time() - t0:.0f}s -> {dst}")
    return True


def main() -> None:
    data_root = Path(__file__).resolve().parents[1].parent / "bulk-download"
    args = sys.argv[1:]
    targets = [Path(a) for a in args] if args else [data_root / t for t in DEFAULT_TARGETS]

    for src in targets:
        if not src.exists():
            print(f"[跳过] 不存在 {src}")
            continue
        dst = FIXED_DIR / src.name
        if dst.exists():
            print(f"[已存在] {dst.name}")
            continue
        if not needs_fix(src):
            print(f"[无需修补] {src.name}")
            continue
        fix(src, dst)

    print("\n修补后校验：")
    import anndata as ad

    for f in sorted(FIXED_DIR.glob("*.h5ad")):
        try:
            a = ad.read_h5ad(f)
            print(f"  OK  {f.name}: {a.shape[0]} cells × {a.shape[1]} genes")
        except Exception as e:
            print(f"  失败 {f.name}: {str(e)[:80]}")


if __name__ == "__main__":
    main()
