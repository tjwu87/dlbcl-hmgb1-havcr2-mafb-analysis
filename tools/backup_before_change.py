"""
改动前备份
================================================================================
**规矩：任何会改变图形尺寸、重出图、覆盖产物的操作之前，先跑这个脚本。**

被保护的对象（默认）
--------------------------------------------------------------------------------
  1. 论文结果/主图/          原始主图（Fig1–Fig7 的位图原件）
  2. 论文结果/附图/          原始补充图
  3. 论文结果/主图_重排版/    第一轮重排版产物
  4. results/figures/        重出的矢量 panel（每跑一次覆盖一次）
  5. results/data_fixed/     修补版 h5ad（体积大，默认跳过）

用法
--------------------------------------------------------------------------------
    python tools/backup_before_change.py --tag 重出图前
    python tools/backup_before_change.py --tag 调字号前 --include-data
    python tools/backup_before_change.py --list        # 只列出现有备份

备份位置：Tcell/_backup/<日期_时间>_<tag>/
================================================================================
"""

from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]      # DLBCL_HMGB1_HAVCR2_MAFB
TCELL = REPO.parent
BACKUP_ROOT = TCELL / "_backup"

TARGETS = [
    ("论文结果/主图", "主图_原始"),
    ("论文结果/附图", "附图_原始"),
    ("论文结果/主图_重排版", "主图_重排版"),
    ("DLBCL_HMGB1_HAVCR2_MAFB/results/figures", "矢量产出"),
]
DATA_TARGET = ("DLBCL_HMGB1_HAVCR2_MAFB/results/data_fixed", "修补数据_h5ad")


def size_mb(p: Path) -> float:
    if not p.exists():
        return 0.0
    if p.is_file():
        return p.stat().st_size / 1048576
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) / 1048576


def do_backup(tag: str, include_data: bool) -> None:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    dest_root = BACKUP_ROOT / f"{stamp}_{tag}"
    dest_root.mkdir(parents=True, exist_ok=True)

    targets = list(TARGETS)
    if include_data:
        targets.append(DATA_TARGET)

    print(f"备份到 {dest_root}\n")
    total = 0.0
    for rel, name in targets:
        src = TCELL / rel
        if not src.exists():
            print(f"  [跳过] {rel}（不存在）")
            continue
        mb = size_mb(src)
        dst = dest_root / name
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)
        total += mb
        print(f"  [完成] {name:16s} {mb:8.1f} MB  ← {rel}")

    manifest = dest_root / "备份说明.txt"
    manifest.write_text(
        f"备份时间：{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"标签：{tag}\n"
        f"内容：\n"
        + "".join(f"  - {n}  ← {r}\n" for r, n in targets)
        + f"\n合计约 {total:.1f} MB\n"
        "\n还原方式：把对应目录复制回原位即可。\n",
        encoding="utf-8")
    print(f"\n合计 {total:.1f} MB，清单写入 {manifest.name}")


def list_backups() -> None:
    if not BACKUP_ROOT.exists():
        print("暂无备份")
        return
    for d in sorted(BACKUP_ROOT.iterdir()):
        if d.is_dir():
            print(f"{d.name}   {sum(size_mb(x) for x in d.iterdir()):.1f} MB")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="手动备份", help="备份标签，写进目录名")
    ap.add_argument("--include-data", action="store_true",
                    help="同时备份修补版 h5ad（体积大）")
    ap.add_argument("--list", action="store_true", help="只列出现有备份")
    args = ap.parse_args()

    if args.list:
        list_backups()
        return
    do_backup(args.tag, args.include_data)


if __name__ == "__main__":
    main()
