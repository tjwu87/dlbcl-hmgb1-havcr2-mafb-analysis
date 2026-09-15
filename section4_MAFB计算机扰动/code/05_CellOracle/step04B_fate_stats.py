# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 04B: fate score 统计比较
- 以 cell_id (obs_names) 为键对齐四个条件，确认顺序完全一致
- Wilcoxon signed-rank test (paired) + BH 校正
- 重点：Mono KO vs WT, KO vs Null; IFN_TAM 同上; LA_TAM KO vs WT (supplementary)
- delta 方向明确定义为 median_KO - median_WT (或 median_KO - median_Null)
输出：
  04_fate/04_fate_stats.csv
  04_fate/04_fate_log.txt
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_04 = os.path.join(OUTDIR, "04_fate")
os.makedirs(OUTDIR_04, exist_ok=True)

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 04B: fate score 统计比较")
log(f"Start: {t0.isoformat()}")
log(f"OUTDIR: {OUTDIR}")
log("="*60)

# ── 加载四个条件的 fate score CSV ─────────────────────────────
log(f"\n--- 加载 fate score CSV ---")
FILE_MAP = {
    "WT":      "cell_fate_scores_WT.csv",
    "MAFB_KO": "cell_fate_scores_KO.csv",
    "MAFB_OE": "cell_fate_scores_OE.csv",
    "Null":    "cell_fate_scores_Null.csv",
}
dfs = {}
for cond, fname in FILE_MAP.items():
    path = os.path.join(OUTDIR_04, fname)
    assert os.path.exists(path), f"FAIL: {path} not found"
    df = pd.read_csv(path)
    assert "cell_id" in df.columns and "subtype" in df.columns and "fate_score" in df.columns, \
        f"FAIL: {fname} missing required columns"
    dfs[cond] = df.set_index("cell_id")
    log(f"  [PASS] {cond}: n={len(df)}, cols={list(df.columns)}")

# ── 对齐检查：以 cell_id 为键，确认顺序完全一致 ──────────────
log(f"\n--- 对齐检查: cell_id 顺序一致性 ---")
ref_ids = list(dfs["WT"].index)
for cond, df in dfs.items():
    ids = list(df.index)
    if ids != ref_ids:
        log(f"  [FAIL] {cond} cell_id 顺序与 WT 不一致！")
        if set(ids) == set(ref_ids):
            log(f"    集合相同但顺序不同 — 尝试重新对齐")
            dfs[cond] = df.loc[ref_ids]
            log(f"    [PASS] {cond} 已按 WT 顺序重新对齐")
        else:
            log(f"    集合也不同 — 严重错误，停止")
            sys.exit(1)
    else:
        log(f"  [PASS] {cond}: cell_id 顺序与 WT 完全一致 (n={len(ids)})")

# 提取 numpy arrays（已对齐）
scores = {cond: dfs[cond]["fate_score"].values for cond in FILE_MAP}
subtypes = dfs["WT"]["subtype"].values
n_cells = len(ref_ids)
log(f"  [PASS] 对齐完成。n_cells={n_cells}")

# ── Cohen's d 计算函数 ────────────────────────────────────────
def cohens_d(a, b):
    """Paired Cohen's d = mean(diff) / std(diff)"""
    diff = a - b
    if diff.std() == 0:
        return np.nan
    return diff.mean() / diff.std()

# ── 统计比较定义 ──────────────────────────────────────────────
# 格式: (subtype, cond_A_label, cond_B_label, cond_A_key, cond_B_key, role)
# delta = median_A - median_B（A 通常是 KO）
COMPARISONS = [
    # Mono — 主要 readout
    ("Mono",    "KO", "WT",   "MAFB_KO", "WT",      "primary"),
    ("Mono",    "KO", "Null", "MAFB_KO", "Null",    "primary"),
    # IFN_TAM — 次级
    ("IFN_TAM", "KO", "WT",   "MAFB_KO", "WT",      "secondary"),
    ("IFN_TAM", "KO", "Null", "MAFB_KO", "Null",    "secondary"),
    # LA_TAM — supplementary only
    ("LA_TAM",  "KO", "WT",   "MAFB_KO", "WT",      "supplementary"),
]

log(f"\n--- 统计比较 ---")
log(f"  delta 方向定义：delta = median_A - median_B（A=KO, B=WT 或 Null）")
log(f"  Wilcoxon signed-rank test (paired) + BH 校正")
log(f"  LA_TAM 仅作 supplementary，在结果中明确标注")

rows = []
for subtype, label_A, label_B, cond_A, cond_B, role in COMPARISONS:
    mask = (subtypes == subtype)
    n = mask.sum()
    a = scores[cond_A][mask]
    b = scores[cond_B][mask]

    # Wilcoxon signed-rank (paired)
    stat, pval = stats.wilcoxon(a, b, alternative='two-sided')

    med_A = np.median(a)
    med_B = np.median(b)
    delta = med_A - med_B  # delta = median_KO - median_WT (or Null)
    cd    = cohens_d(a, b)

    rows.append({
        "subtype":      subtype,
        "comparison":   f"{label_A}_vs_{label_B}",
        "cond_A":       cond_A,
        "cond_B":       cond_B,
        "median_A":     round(med_A, 6),
        "median_B":     round(med_B, 6),
        "delta":        round(delta, 6),
        "delta_direction": f"median_{label_A} - median_{label_B}",
        "p_value":      pval,
        "padj":         np.nan,   # filled after BH
        "cohens_d":     round(cd, 6) if not np.isnan(cd) else np.nan,
        "n_cells":      n,
        "role":         role,
    })
    log(f"  {subtype} {label_A} vs {label_B}: "
        f"median_A={med_A:.4f}, median_B={med_B:.4f}, "
        f"delta={delta:.4f}, p={pval:.4e}, d={cd:.3f}, n={n} [{role}]")

# ── BH 校正（全部比较一起校正）────────────────────────────────
log(f"\n--- BH 多重检验校正 ---")
pvals = [r["p_value"] for r in rows]
_, padjs, _, _ = multipletests(pvals, method="fdr_bh")
for i, r in enumerate(rows):
    r["padj"] = round(padjs[i], 6)
    log(f"  {r['subtype']} {r['comparison']}: p={r['p_value']:.4e} → padj={r['padj']:.4e}")

# ── 保存 04_fate_stats.csv ────────────────────────────────────
stats_df = pd.DataFrame(rows, columns=[
    "subtype", "comparison", "cond_A", "cond_B",
    "median_A", "median_B", "delta", "delta_direction",
    "p_value", "padj", "cohens_d", "n_cells", "role"
])
out_stats = os.path.join(OUTDIR_04, "04_fate_stats.csv")
stats_df.to_csv(out_stats, index=False)
log(f"\n  [PASS] 04_fate_stats.csv saved: {out_stats}")

# ── 打印完整统计表 ────────────────────────────────────────────
log(f"\n{'='*60}")
log(f"完整统计表：")
log(f"{'亚群':<10} {'比较':<15} {'median_A':>9} {'median_B':>9} {'delta':>8} {'p_value':>10} {'padj':>10} {'Cohen_d':>8} {'n':>4} {'role'}")
log(f"{'-'*100}")
for _, r in stats_df.iterrows():
    sig = "**" if r['padj'] < 0.01 else ("*" if r['padj'] < 0.05 else "ns")
    log(f"{r['subtype']:<10} {r['comparison']:<15} {r['median_A']:>9.4f} {r['median_B']:>9.4f} "
        f"{r['delta']:>8.4f} {r['p_value']:>10.4e} {r['padj']:>10.4e} "
        f"{r['cohens_d']:>8.3f} {r['n_cells']:>4}  {r['role']} {sig}")

# ── 保存 04_fate_log.txt ──────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 04B COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if total/60>45 else 'NOT triggered'}")
log(f"stall: NO")
log("="*60)
log(f"\n停止等待确认。")

log_path = os.path.join(OUTDIR_04, "04_fate_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 04_fate_log.txt saved")
