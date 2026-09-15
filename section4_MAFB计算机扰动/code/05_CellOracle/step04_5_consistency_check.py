# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 04.5: 一致性验证（强制门控）
V1: SNR 排序检查
V2: Mono 主 readout 检查（必须 STOP 条件）
V3: IFN_TAM 检查（WARN 条件）
V4: LA_TAM 自身 fate 检查（supplementary）
V5: SNR 幅度变化检查（Phase B k=40 参考值）

Phase B (k=40) 参考值：Mono≈2.05, IFN_TAM≈1.84, LA_TAM≈3.10
输出：04_fate/consistency_check.txt
"""
import os, sys
import numpy as np, pandas as pd
from datetime import datetime

t0 = datetime.now()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_04 = os.path.join(OUTDIR, "04_fate")
log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 04.5: 一致性验证（强制门控）")
log(f"Start: {t0.isoformat()}")
log("="*60)

# ── 加载 fate score CSV ───────────────────────────────────────
FILE_MAP = {
    "WT":      "cell_fate_scores_WT.csv",
    "MAFB_KO": "cell_fate_scores_KO.csv",
}
dfs = {}
for cond, fname in FILE_MAP.items():
    path = os.path.join(OUTDIR_04, fname)
    df = pd.read_csv(path).set_index("cell_id")
    dfs[cond] = df
    log(f"  [PASS] {cond}: n={len(df)}")

# ── 加载 04_fate_stats.csv ────────────────────────────────────
stats_path = os.path.join(OUTDIR_04, "04_fate_stats.csv")
stats_df = pd.read_csv(stats_path)
log(f"  [PASS] 04_fate_stats.csv loaded: {len(stats_df)} rows")

# ── 计算 SNR（per subtype）────────────────────────────────────
# SNR = mean_KO_fate / std_KO_fate (per subtype)
# 计划书未明确 SNR 公式，使用标准定义：mean / std
log(f"\n--- SNR 计算 ---")
subtypes_order = ["Mono", "IFN_TAM", "LA_TAM"]
snr = {}
for st in subtypes_order:
    mask = dfs["MAFB_KO"]["subtype"] == st
    vals = dfs["MAFB_KO"].loc[mask, "fate_score"].values
    mean_v = vals.mean()
    std_v  = vals.std()
    snr_v  = mean_v / std_v if std_v > 0 else np.nan
    snr[st] = snr_v
    log(f"  {st}: mean={mean_v:.4f}, std={std_v:.4f}, SNR={snr_v:.4f}")

# Phase B (k=40) 参考值
PHASE_B_REF = {"Mono": 2.05, "IFN_TAM": 1.84, "LA_TAM": 3.10}

# ── 从 stats_df 提取关键数值 ──────────────────────────────────
def get_stat(subtype, comparison):
    row = stats_df[(stats_df["subtype"]==subtype) & (stats_df["comparison"]==comparison)]
    if len(row) == 0:
        return None
    return row.iloc[0]

mono_ko_wt   = get_stat("Mono",    "KO_vs_WT")
mono_ko_null = get_stat("Mono",    "KO_vs_Null")
ifn_ko_wt    = get_stat("IFN_TAM", "KO_vs_WT")
ifn_ko_null  = get_stat("IFN_TAM", "KO_vs_Null")
la_ko_wt     = get_stat("LA_TAM",  "KO_vs_WT")

# ── 状态追踪 ─────────────────────────────────────────────────
status_flags = []   # list of (check, status, detail)
STOP_triggered = False
WARN_triggered = False

log(f"\n{'='*60}")
log(f"V1: SNR 排序检查")
log(f"{'='*60}")
# 要求：LA_TAM SNR 必须是三个亚群中最高值
la_snr_max = snr["LA_TAM"] == max(snr.values())
log(f"  SNR: Mono={snr['Mono']:.4f}, IFN_TAM={snr['IFN_TAM']:.4f}, LA_TAM={snr['LA_TAM']:.4f}")
if not la_snr_max:
    log(f"  [STOP] LA_TAM SNR 不是最高值！触发 STOP。")
    status_flags.append(("V1_LA_TAM_SNR_max", "STOP", f"LA_TAM SNR={snr['LA_TAM']:.4f} 不是最高"))
    STOP_triggered = True
else:
    log(f"  [PASS] LA_TAM SNR 是最高值 ✓")
    status_flags.append(("V1_LA_TAM_SNR_max", "PASS", f"LA_TAM SNR={snr['LA_TAM']:.4f} 最高"))

# Mono vs IFN_TAM 相对顺序
if snr["IFN_TAM"] > snr["Mono"]:
    log(f"  [WARN] IFN_TAM SNR ({snr['IFN_TAM']:.4f}) > Mono SNR ({snr['Mono']:.4f}) — 顺序翻转")
    log(f"    Phase B 参考值两者差距极小，轻微数值波动可导致翻转，生物学意义有限")
    status_flags.append(("V1_Mono_IFN_order", "WARN", f"IFN_TAM SNR > Mono SNR (差距小，生物学意义有限)"))
    WARN_triggered = True
else:
    log(f"  [PASS] Mono SNR ({snr['Mono']:.4f}) > IFN_TAM SNR ({snr['IFN_TAM']:.4f}) ✓")
    status_flags.append(("V1_Mono_IFN_order", "PASS", f"Mono SNR > IFN_TAM SNR"))

log(f"\n{'='*60}")
log(f"V2: Mono 主 readout 检查（必须 STOP 条件）")
log(f"{'='*60}")
# 要求：Mono KO vs WT delta 必须为下降方向，且 padj < 0.05
mono_delta = mono_ko_wt["delta"]
mono_padj  = mono_ko_wt["padj"]
log(f"  Mono KO vs WT: delta={mono_delta:.4f}, padj={mono_padj:.4e}")
if mono_delta >= 0:
    log(f"  [STOP] Mono KO vs WT delta 方向错误（应为下降，实为 {mono_delta:.4f}）")
    status_flags.append(("V2_Mono_direction", "STOP", f"delta={mono_delta:.4f} 方向错误"))
    STOP_triggered = True
elif mono_padj >= 0.05:
    log(f"  [STOP] Mono KO vs WT padj={mono_padj:.4e} >= 0.05，不显著")
    status_flags.append(("V2_Mono_significance", "STOP", f"padj={mono_padj:.4e} >= 0.05"))
    STOP_triggered = True
else:
    log(f"  [PASS] Mono KO vs WT: delta={mono_delta:.4f} (下降) ✓, padj={mono_padj:.4e} < 0.05 ✓")
    status_flags.append(("V2_Mono_readout", "PASS",
        f"delta={mono_delta:.4f} 下降, padj={mono_padj:.4e}"))

log(f"\n{'='*60}")
log(f"V3: IFN_TAM 检查")
log(f"{'='*60}")
# 若 IFN_TAM 出现强显著结果（padj < 0.05 且效应明显），标记 WARN
ifn_delta = ifn_ko_wt["delta"]
ifn_padj  = ifn_ko_wt["padj"]
ifn_d     = ifn_ko_wt["cohens_d"]
log(f"  IFN_TAM KO vs WT: delta={ifn_delta:.4f}, padj={ifn_padj:.4e}, Cohen's d={ifn_d:.3f}")
if ifn_padj < 0.05 and abs(ifn_d) > 0.3:
    log(f"  [WARN] IFN_TAM 出现强显著结果 (padj={ifn_padj:.4e}, |d|={abs(ifn_d):.3f})")
    log(f"    需在最终报告中解释其与历史不稳定性的关系")
    status_flags.append(("V3_IFN_TAM", "WARN",
        f"padj={ifn_padj:.4e} < 0.05 且 |d|={abs(ifn_d):.3f} > 0.3"))
    WARN_triggered = True
else:
    log(f"  [PASS] IFN_TAM 无强显著结果 (padj={ifn_padj:.4e}, |d|={abs(ifn_d):.3f})")
    status_flags.append(("V3_IFN_TAM", "PASS",
        f"padj={ifn_padj:.4e}, |d|={abs(ifn_d):.3f}"))

log(f"\n{'='*60}")
log(f"V4: LA_TAM 自身 fate 检查（supplementary）")
log(f"{'='*60}")
la_delta = la_ko_wt["delta"]
la_padj  = la_ko_wt["padj"]
log(f"  LA_TAM KO vs WT: delta={la_delta:.4f}, padj={la_padj:.4e}")
# 若方向与 Mono 不同，标记 WARN
if (la_delta > 0) == (mono_delta > 0):
    log(f"  [WARN] LA_TAM delta 方向与 Mono 相同（均为 {'上升' if la_delta>0 else '下降'}）")
    log(f"    LA_TAM 自身 fate 上升（KO 后 LA_TAM 细胞更倾向 LA_TAM fate）是 state-dependent observation")
    log(f"    不得自动升级为主结论，仅在 supplementary 中降级解释")
    status_flags.append(("V4_LA_TAM_direction", "WARN",
        f"LA_TAM delta={la_delta:.4f} 与 Mono delta={mono_delta:.4f} 方向相同"))
    WARN_triggered = True
else:
    log(f"  [PASS] LA_TAM delta 方向与 Mono 不同（LA_TAM={la_delta:.4f}, Mono={mono_delta:.4f}）")
    status_flags.append(("V4_LA_TAM_direction", "PASS",
        f"LA_TAM delta={la_delta:.4f} 与 Mono delta={mono_delta:.4f} 方向不同"))

log(f"\n{'='*60}")
log(f"V5: SNR 幅度变化检查（Phase B k=40 参考值）")
log(f"{'='*60}")
log(f"  Phase B (k=40) 参考值: Mono={PHASE_B_REF['Mono']}, IFN_TAM={PHASE_B_REF['IFN_TAM']}, LA_TAM={PHASE_B_REF['LA_TAM']}")
log(f"  注意：主分析使用 knn_imputation_k=20，与 Phase B k=40 不同，")
log(f"        V5 仅用于量级参考，不作为严格数值一致性判定依据。")
log(f"        优先比较：(1) SNR 排序稳定性；(2) Mono fate 主结论稳定性；(3) 主要图形结构。")
for st in subtypes_order:
    ref = PHASE_B_REF[st]
    cur = snr[st]
    rel_change = abs(cur - ref) / ref * 100
    if rel_change > 50:
        flag = "HIGH-WARN"
        WARN_triggered = True
    elif rel_change > 20:
        flag = "WARN"
        WARN_triggered = True
    else:
        flag = "OK"
    log(f"  {st}: current={cur:.4f}, ref={ref:.2f}, rel_change={rel_change:.1f}% → [{flag}]")
    status_flags.append((f"V5_SNR_{st}", flag,
        f"current={cur:.4f}, ref={ref:.2f}, rel_change={rel_change:.1f}%"))
    if flag in ("WARN", "HIGH-WARN"):
        log(f"    → 必须在最终报告中单列解释原因（k=20 vs k=40 差异）")

# ── 总状态判定 ────────────────────────────────────────────────
log(f"\n{'='*60}")
log(f"总状态判定")
log(f"{'='*60}")
if STOP_triggered:
    overall = "STOP"
elif WARN_triggered:
    overall = "WARN"
else:
    overall = "PASS"

log(f"\n  总状态: {overall}")
log(f"\n  各检查项汇总：")
for check, status, detail in status_flags:
    log(f"    {check:<30} [{status:<9}] {detail}")

if overall == "STOP":
    log(f"\n  [STOP] 触发 STOP 的检查项：")
    for check, status, detail in status_flags:
        if status == "STOP":
            log(f"    {check}: {detail}")
    log(f"\n  可能原因（至少3个）：")
    log(f"    1. knn_imputation_k 参数设置影响了 transition probability 计算")
    log(f"    2. GRN 拟合质量不足，导致 MAFB 扰动效应不显著")
    log(f"    3. 细胞数量过少（n=272），统计功效不足")
    log(f"\n  备选方案：")
    log(f"    方案 X：保留严格主分析但降级结论（仅报告探索性结果）")
    log(f"    方案 Y：回到 STEP 02 检查主参数是否被错误实现")
    log(f"\n  等待人工确认。")
elif overall == "WARN":
    log(f"\n  [WARN] 存在警告项，但允许继续。")
    log(f"  必须在最终报告中说明所有 WARN 项。")
    log(f"  继续执行后续步骤。")
else:
    log(f"\n  [PASS] 所有检查通过，继续执行后续步骤。")

# ── 保存 consistency_check.txt ────────────────────────────────
out_path = os.path.join(OUTDIR_04, "consistency_check.txt")
with open(out_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
    f.write(f"\nOVERALL_STATUS: {overall}\n")
log(f"\n  [PASS] consistency_check.txt saved: {out_path}")

total = (datetime.now()-t0).total_seconds()
log(f"\nSTEP 04.5 COMPLETE. Total: {total:.1f}s")
log(f"OVERALL_STATUS: {overall}")

# 若 STOP，以非零退出码退出
if overall == "STOP":
    sys.exit(2)
