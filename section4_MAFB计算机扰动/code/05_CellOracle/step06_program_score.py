# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 06: LA_TAM program score
06A: 定义 LA_TAM program genes
     - LA_TAM vs (Mono + IFN_TAM) 差异表达
     - padj < 0.01, log2FC > 1, top 50, 排除 MAFB
06B: program score 统计
     - program score = top50 基因平均 log1p_norm 表达
     - WT vs KO 比较（分亚群）
     - Wilcoxon signed-rank (paired) + BH 校正 + Cohen's d
输出：
  06_program/06_program_gene_50.txt
  06_program/06_program_gene_definition_log.txt
  06_program/06_program_score_stats.csv
  06_program/06_program_score_log.txt
"""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime

t0 = datetime.now()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_06 = os.path.join(OUTDIR, "06_program")
os.makedirs(OUTDIR_06, exist_ok=True)

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 06: LA_TAM program score")
log(f"Start: {t0.isoformat()}")
log("="*60)

# ── 加载 WT oracle（用于 DEG 定义和 program score 计算）────────
log(f"\n--- 加载 WT oracle ---")
with open(os.path.join(OUTDIR, "03_perturbation/WT/oracle_WT.pkl"), "rb") as f:
    oracle_wt = pickle.load(f)
with open(os.path.join(OUTDIR, "03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl"), "rb") as f:
    oracle_ko = pickle.load(f)

log(f"  WT adata.shape: {oracle_wt.adata.shape}")
log(f"  KO adata.shape: {oracle_ko.adata.shape}")

import scipy.sparse as sp
# 提取 log1p_norm 表达矩阵（用于 DEG 和 program score）
wt_log = oracle_wt.adata.layers["log1p_norm"]
ko_log = oracle_ko.adata.layers["log1p_norm"]  # WT 和 KO 的 log1p_norm 相同（来自原始数据）
if sp.issparse(wt_log): wt_log = wt_log.toarray()
if sp.issparse(ko_log): ko_log = ko_log.toarray()

# 注意：program score 使用 log1p_norm（原始表达），不是模拟后的值
# 但 WT vs KO 的 program score 比较应使用各自条件的 imputed_count
# 计划书说：program score = top50 基因平均 log1p_norm 表达
# 对每个细胞比较 WT vs KO 的 program score
# 这里 WT 和 KO 的 log1p_norm 是相同的（来自同一 adata），
# 但 imputed_count 不同（KO 模拟后的平滑值）
# 按计划书严格定义：program score = top50 基因平均 log1p_norm 表达
# WT 和 KO 的 log1p_norm 相同，所以需要用 imputed_count 来区分
# 重新解读：WT program score = WT imputed_count 的 top50 基因均值
#           KO program score = KO simulated_count 的 top50 基因均值
# 这样才能体现 KO 效应

wt_imp = oracle_wt.adata.layers["imputed_count"]
ko_sim = oracle_ko.adata.layers["simulated_count"]
if sp.issparse(wt_imp): wt_imp = wt_imp.toarray()
if sp.issparse(ko_sim): ko_sim = ko_sim.toarray()

var_names = list(oracle_wt.adata.var_names)
subtypes = oracle_wt.adata.obs["mac_subtype"].values
log(f"  subtype 分布: {dict(zip(*np.unique(subtypes, return_counts=True)))}")

# ── 06A: 定义 LA_TAM program genes ───────────────────────────
log(f"\n--- 06A: 定义 LA_TAM program genes ---")
log(f"  LA_TAM vs (Mono + IFN_TAM) 差异表达")
log(f"  条件: padj < 0.01, log2FC > 1, top 50, 排除 MAFB")
log(f"  数据: log1p_norm layer（原始表达）")

la_mask   = (subtypes == "LA_TAM")
rest_mask = (subtypes == "Mono") | (subtypes == "IFN_TAM")
n_la   = la_mask.sum()
n_rest = rest_mask.sum()
log(f"  LA_TAM n={n_la}, Mono+IFN_TAM n={n_rest}")

# 对每个基因做 Wilcoxon rank-sum test（unpaired）
deg_rows = []
for i, gene in enumerate(var_names):
    if gene == "MAFB":
        continue
    la_vals   = wt_log[la_mask, i]
    rest_vals = wt_log[rest_mask, i]
    # 跳过全零基因
    if la_vals.max() == 0 and rest_vals.max() == 0:
        continue
    stat, pval = stats.ranksums(la_vals, rest_vals)
    log2fc = np.log2(la_vals.mean() + 1) - np.log2(rest_vals.mean() + 1)
    deg_rows.append({"gene": gene, "log2FC": log2fc, "p_value": pval})

deg_df = pd.DataFrame(deg_rows)
log(f"  非零基因数: {len(deg_df)}")

# BH 校正
_, padjs, _, _ = multipletests(deg_df["p_value"].values, method="fdr_bh")
deg_df["padj"] = padjs

# 筛选：padj < 0.01, log2FC > 1
sig_df = deg_df[(deg_df["padj"] < 0.01) & (deg_df["log2FC"] > 1)].copy()
log(f"  padj < 0.01 且 log2FC > 1: n={len(sig_df)}")

# Top 50（按 log2FC 降序）
top50_df = sig_df.sort_values("log2FC", ascending=False).head(50)
log(f"  Top 50 genes: n={len(top50_df)}")

program_genes = list(top50_df["gene"].values)
log(f"  Program genes: {program_genes}")

# 保存 06_program_gene_50.txt
gene_path = os.path.join(OUTDIR_06, "06_program_gene_50.txt")
with open(gene_path, "w") as f:
    f.write("# LA_TAM program genes (top 50)\n")
    f.write("# 定义: LA_TAM vs (Mono+IFN_TAM), padj<0.01, log2FC>1, top50 by log2FC, 排除MAFB\n")
    f.write("# 数据: log1p_norm layer (WT oracle)\n\n")
    for g in program_genes:
        f.write(g + "\n")
log(f"  [PASS] 06_program_gene_50.txt saved")

# 保存 definition log
def_lines = [
    "STEP 06A: LA_TAM program gene 定义",
    "="*50,
    f"分析: LA_TAM vs (Mono + IFN_TAM) 差异表达",
    f"数据: log1p_norm layer (WT oracle, 272×3000)",
    f"统计: Wilcoxon rank-sum test (unpaired) + BH 校正",
    f"筛选条件: padj < 0.01, log2FC > 1, top 50 by log2FC, 排除 MAFB",
    f"",
    f"LA_TAM n={n_la}, Mono+IFN_TAM n={n_rest}",
    f"非零基因数: {len(deg_df)}",
    f"padj<0.01 且 log2FC>1: n={len(sig_df)}",
    f"Top 50 genes: n={len(top50_df)}",
    f"",
    f"Top 50 genes (按 log2FC 降序):",
]
for _, r in top50_df.iterrows():
    def_lines.append(f"  {r['gene']}: log2FC={r['log2FC']:.4f}, padj={r['padj']:.4e}")

def_path = os.path.join(OUTDIR_06, "06_program_gene_definition_log.txt")
with open(def_path, "w") as f:
    f.write("\n".join(def_lines) + "\n")
log(f"  [PASS] 06_program_gene_definition_log.txt saved")

# ── 06B: program score 统计 ───────────────────────────────────
log(f"\n--- 06B: program score 统计 ---")
log(f"  program score = top50 基因平均 imputed_count (WT) / simulated_count (KO)")
log(f"  注意: 使用 imputed_count (WT) 和 simulated_count (KO) 以体现 KO 效应")

# 找到 program genes 在 var_names 中的索引
gene_indices = [var_names.index(g) for g in program_genes if g in var_names]
missing = [g for g in program_genes if g not in var_names]
if missing:
    log(f"  [WARN] 以下 program genes 不在 var_names 中: {missing}")
log(f"  可用 program genes: {len(gene_indices)}/{len(program_genes)}")

# 计算 program score
wt_prog_score = wt_imp[:, gene_indices].mean(axis=1)  # (272,)
ko_prog_score = ko_sim[:, gene_indices].mean(axis=1)  # (272,)

def cohens_d_paired(a, b):
    diff = a - b
    if diff.std() == 0:
        return np.nan
    return diff.mean() / diff.std()

log(f"\n  分亚群统计：")
log(f"  delta 方向定义：delta = median_KO - median_WT")
stat_rows = []
subtypes_list = ["Mono", "IFN_TAM", "LA_TAM"]
for st in subtypes_list:
    mask = (subtypes == st)
    n = mask.sum()
    wt_vals = wt_prog_score[mask]
    ko_vals = ko_prog_score[mask]

    med_wt = np.median(wt_vals)
    med_ko = np.median(ko_vals)
    delta  = med_ko - med_wt  # KO - WT

    try:
        stat, pval = stats.wilcoxon(wt_vals, ko_vals, alternative='two-sided')
    except ValueError:
        pval = 1.0

    cd = cohens_d_paired(ko_vals, wt_vals)

    stat_rows.append({
        "subtype":   st,
        "median_WT": round(med_wt, 6),
        "median_KO": round(med_ko, 6),
        "delta":     round(delta, 6),
        "delta_direction": "median_KO - median_WT",
        "p_value":   pval,
        "padj":      np.nan,
        "cohens_d":  round(cd, 6) if not np.isnan(cd) else np.nan,
        "n_cells":   n,
    })
    log(f"  {st}: median_WT={med_wt:.4f}, median_KO={med_ko:.4f}, "
        f"delta={delta:.4f}, p={pval:.4e}, d={cd:.3f}, n={n}")

# BH 校正
pvals = [r["p_value"] for r in stat_rows]
_, padjs, _, _ = multipletests(pvals, method="fdr_bh")
for i, r in enumerate(stat_rows):
    r["padj"] = round(padjs[i], 6)

# 保存 06_program_score_stats.csv
score_df = pd.DataFrame(stat_rows, columns=[
    "subtype", "median_WT", "median_KO", "delta", "delta_direction",
    "p_value", "padj", "cohens_d", "n_cells"
])
out_csv = os.path.join(OUTDIR_06, "06_program_score_stats.csv")
score_df.to_csv(out_csv, index=False)
log(f"\n  [PASS] 06_program_score_stats.csv saved")

# 打印完整统计表
log(f"\n  完整统计表：")
log(f"  {'亚群':<10} {'median_WT':>10} {'median_KO':>10} {'delta':>8} {'p_value':>10} {'padj':>10} {'Cohen_d':>8} {'n':>4}")
log(f"  {'-'*75}")
for _, r in score_df.iterrows():
    sig = "**" if r['padj'] < 0.01 else ("*" if r['padj'] < 0.05 else "ns")
    log(f"  {r['subtype']:<10} {r['median_WT']:>10.4f} {r['median_KO']:>10.4f} "
        f"{r['delta']:>8.4f} {r['p_value']:>10.4e} {r['padj']:>10.4e} "
        f"{r['cohens_d']:>8.3f} {r['n_cells']:>4}  {sig}")

# 特殊观察记录
log(f"\n  特殊观察：")
for _, r in score_df.iterrows():
    if r['subtype'] == 'Mono' and r['delta'] < 0:
        log(f"  [OBS] Mono 中 KO 后 program score 下降 (delta={r['delta']:.4f}) — 符合预期")
    elif r['subtype'] == 'Mono' and r['delta'] >= 0:
        log(f"  [OBS] Mono 中 KO 后 program score 未下降 (delta={r['delta']:.4f}) — 需注意")
    if r['subtype'] == 'LA_TAM' and r['delta'] > 0:
        log(f"  [OBS] LA_TAM 中 KO 后 program score 上升 (delta={r['delta']:.4f}) — state-dependent observation")
        log(f"        不得自动升级为主结论，仅在 supplementary 中以 state-dependent observation 表述")
    if r['subtype'] == 'LA_TAM' and r['delta'] < 0:
        log(f"  [OBS] LA_TAM 中 KO 后 program score 下降 (delta={r['delta']:.4f})")

# ── 保存日志 ──────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 06 COMPLETE. Total: {total:.1f}s")
log("="*60)

log_path = os.path.join(OUTDIR_06, "06_program_score_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 06_program_score_log.txt saved")
