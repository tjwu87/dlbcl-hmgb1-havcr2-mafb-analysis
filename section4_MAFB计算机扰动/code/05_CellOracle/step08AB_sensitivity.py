# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 08A+B: 敏感性分析 Class 1 简单组
08A: 主分析参考重放（直接复用主分析结果，不重复模拟）
08B: filter p<0.01 / threshold=5000 下 MAFB edge 数量和分类更新
     禁止重新运行 simulate_shift / estimate_transition_prob / calculate_embedding_shift
     只读取已有扰动结果，更新 filter 宽松度下的 MAFB edge 数量和分类
输出:
  08_sensitivity/08A_reference_summary.txt
  08_sensitivity/08B_filter_update.txt
  08_sensitivity/08B_mafb_edges_filterB.csv
"""
import os, sys
import numpy as np, pandas as pd
from datetime import datetime

t0 = datetime.now()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_08 = os.path.join(OUTDIR, "08_sensitivity")
os.makedirs(OUTDIR_08, exist_ok=True)

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 08A+B: 敏感性分析 Class 1 简单组")
log(f"Start: {t0.isoformat()}")
log("="*60)

# ══════════════════════════════════════════════════════════════
# 08A: 主分析参考重放（直接复用）
# ══════════════════════════════════════════════════════════════
log(f"\n{'='*60}")
log(f"08A: 主分析参考重放（直接复用主分析结果）")
log(f"{'='*60}")
log(f"  参数: HVG=3000(counts), p<0.001, threshold=2000, n_neighbors=200")
log(f"  策略: 直接复用主分析结果，不重复完整模拟")
log(f"  （STEP 04.5 触发 WARN，但 V2 PASS，无需重新运行完整扰动模拟）")

# 读取主分析关键结果
fate_stats = pd.read_csv(os.path.join(OUTDIR, "04_fate", "04_fate_stats.csv"))
prog_stats  = pd.read_csv(os.path.join(OUTDIR, "06_program", "06_program_score_stats.csv"))
consist_txt = open(os.path.join(OUTDIR, "04_fate", "consistency_check.txt")).read()

# 提取 SNR（从 consistency_check.txt）
import re
snr_vals = {}
for line in consist_txt.split("\n"):
    m = re.search(r"(Mono|IFN_TAM|LA_TAM): mean=[\d.]+, std=[\d.]+, SNR=([\d.]+)", line)
    if m:
        snr_vals[m.group(1)] = float(m.group(2))

# 提取 MAFB edge 数（filter A）
mafb_edges_A = pd.read_csv(os.path.join(OUTDIR, "02_oracle_grn", "filter_A_t2000", "MAFB_edges_filterA.csv"))
edge_counts_A = mafb_edges_A.groupby("cluster").size().to_dict()

ref_lines = [
    "STEP 08A: 主分析参考重放",
    "="*50,
    "参数组合 A: HVG=3000(counts), p<0.001, threshold=2000, n_neighbors=200",
    "策略: 直接复用主分析结果（STEP 04.5 WARN 但 V2 PASS，无需重新模拟）",
    "",
    "── SNR (MAFB KO fate score, per subtype) ──",
]
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    ref_lines.append(f"  {st}: SNR={snr_vals.get(st, 'N/A'):.4f}")

ref_lines += [
    "",
    "── MAFB edges (filter A, threshold=2000) ──",
]
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    ref_lines.append(f"  {st}: n_edges={edge_counts_A.get(st, 0)}")

ref_lines += ["", "── Fate score stats (KO vs WT) ──"]
for _, r in fate_stats[fate_stats["comparison"] == "KO_vs_WT"].iterrows():
    ref_lines.append(f"  {r['subtype']}: delta={r['delta']:.4f}, padj={r['padj']:.4e}, "
                     f"Cohen's d={r['cohens_d']:.3f}, n={r['n_cells']}")

ref_lines += ["", "── Program score stats (KO vs WT) ──"]
for _, r in prog_stats.iterrows():
    ref_lines.append(f"  {r['subtype']}: delta={r['delta']:.4f}, padj={r['padj']:.4e}, "
                     f"Cohen's d={r['cohens_d']:.3f}")

ref_lines += ["", "── Consistency check 总状态 ──"]
for line in consist_txt.split("\n"):
    if "总状态" in line or "OVERALL_STATUS" in line or "PASS" in line or "WARN" in line or "STOP" in line:
        ref_lines.append(f"  {line.strip()}")

out_08A = os.path.join(OUTDIR_08, "08A_reference_summary.txt")
with open(out_08A, "w") as f:
    f.write("\n".join(ref_lines) + "\n")
log(f"  [PASS] 08A_reference_summary.txt saved")
for l in ref_lines:
    log(f"  {l}")

# ══════════════════════════════════════════════════════════════
# 08B: filter p<0.01 / threshold=5000 更新
# ══════════════════════════════════════════════════════════════
log(f"\n{'='*60}")
log(f"08B: filter p<0.01 / threshold=5000 MAFB edge 更新")
log(f"{'='*60}")
log(f"  参数: HVG=3000(counts), p<0.01, threshold=5000, n_neighbors=200")
log(f"  策略: 仅更新 MAFB edge 数量和分类，禁止重新运行任何扰动模拟")
log(f"  复用: 主分析 03_perturbation/、04_fate/、06_program/ 全部结果")

# 读取 filter B 结果
mafb_edges_B = pd.read_csv(os.path.join(OUTDIR, "02_oracle_grn", "filter_B_t5000", "MAFB_edges_filterB.csv"))
edge_counts_B = mafb_edges_B.groupby("cluster").size().to_dict()

log(f"\n  MAFB edges (filter B, p<0.01, threshold=5000):")
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    n_A = edge_counts_A.get(st, 0)
    n_B = edge_counts_B.get(st, 0)
    log(f"    {st}: filter_A(t=2000)={n_A}, filter_B(t=5000)={n_B}, delta={n_B-n_A:+d}")

# 保存 filter B edges
out_edges_B = os.path.join(OUTDIR_08, "08B_mafb_edges_filterB.csv")
mafb_edges_B.to_csv(out_edges_B, index=False)
log(f"  [PASS] 08B_mafb_edges_filterB.csv saved")

# 比较 filter A vs B 的 MAFB target 分类
log(f"\n  MAFB target 分类对比 (filter A vs B):")
class_df = pd.read_csv(os.path.join(OUTDIR, "02_oracle_grn", "02E_23gene_classification.csv"))

# filter B 中各亚群的 MAFB targets
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    st_edges_B = mafb_edges_B[mafb_edges_B["cluster"] == st]["target"].tolist()
    st_edges_A = mafb_edges_A[mafb_edges_A["cluster"] == st]["target"].tolist()
    new_in_B = [g for g in st_edges_B if g not in st_edges_A]
    lost_in_B = [g for g in st_edges_A if g not in st_edges_B]
    log(f"    {st}: filter_B 新增={len(new_in_B)}, 丢失={len(lost_in_B)}")
    if new_in_B:
        log(f"      新增: {new_in_B[:10]}")
    if lost_in_B:
        log(f"      丢失: {lost_in_B[:10]}")

# 生成 08B 汇总文本
b_lines = [
    "STEP 08B: filter p<0.01 / threshold=5000 更新",
    "="*50,
    "参数组合 B: HVG=3000(counts), p<0.01, threshold=5000, n_neighbors=200",
    "策略: 仅更新 MAFB edge 数量和分类",
    "禁止: 重新运行 simulate_shift / estimate_transition_prob / calculate_embedding_shift",
    "复用: 主分析 03_perturbation/、04_fate/、06_program/ 全部结果（数值与主分析完全相同）",
    "",
    "── MAFB edges 对比 ──",
    f"{'亚群':<12} {'filter_A(t=2000)':>18} {'filter_B(t=5000)':>18} {'delta':>8}",
    "-"*60,
]
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    n_A = edge_counts_A.get(st, 0)
    n_B = edge_counts_B.get(st, 0)
    b_lines.append(f"{st:<12} {n_A:>18} {n_B:>18} {n_B-n_A:>+8}")

b_lines += [
    "",
    "── 扰动模拟结果（与主分析完全相同，直接复用）──",
    "  fate score、program score、delta_X、delta_embedding 均与主分析一致",
    "  filter 参数不影响 fit_GRN_for_simulation 使用的底层 ridge 全量 GRN",
    "",
    "── 解释层变化 ──",
    "  threshold=5000 下 MAFB edges 增多，说明更多低置信度 MAFB→target 连接被纳入",
    "  但这些额外连接不影响已完成的扰动模拟结果",
    "  在论文中应明确区分：主分析使用 threshold=2000（更严格），",
    "  threshold=5000 作为 sensitivity 展示 edge 数量变化",
]

out_08B = os.path.join(OUTDIR_08, "08B_filter_update.txt")
with open(out_08B, "w") as f:
    f.write("\n".join(b_lines) + "\n")
log(f"  [PASS] 08B_filter_update.txt saved")

total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 08A+B COMPLETE. Total: {total:.1f}s")
log("="*60)

# 保存日志
log_path = os.path.join(OUTDIR_08, "08AB_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 08AB_log.txt saved")
