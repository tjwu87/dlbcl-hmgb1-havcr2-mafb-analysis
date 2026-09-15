# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 05: 23 个 pySCENIC 靶基因 KO 效应分析
05A: 确定可评估与不可评估基因
05B: WT imputed_count vs KO simulated_count 基因层面统计
     - Wilcoxon signed-rank test (paired)
     - BH 校正 padj
     - log2FC
     - Cohen's d（若方差极小则标记 numeric_instability）
输出：
  05_targets/05_gene_evaluation_scope.txt
  05_targets/05_scenic_target_stats.csv
  05_targets/05_scenic_target_log.txt
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

OUTDIR_05 = os.path.join(OUTDIR, "05_targets")
os.makedirs(OUTDIR_05, exist_ok=True)

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 05: 23 个 pySCENIC 靶基因 KO 效应分析")
log(f"Start: {t0.isoformat()}")
log("="*60)

# ── 加载 23 基因分类表 ────────────────────────────────────────
class_path = os.path.join(OUTDIR, "02_oracle_grn", "02E_23gene_classification.csv")
class_df = pd.read_csv(class_path)
log(f"\n--- 加载 23 基因分类表 ---")
log(f"  n_genes: {len(class_df)}")
log(f"  columns: {list(class_df.columns)}")

# ── 05A: 确定可评估与不可评估基因 ────────────────────────────
log(f"\n--- 05A: 确定可评估与不可评估基因 ---")
# 可评估 = in_HVG3000=True（class C 或 indirect，且在 oracle 中）
# 不可评估 = in_HVG3000=False（class E）
evaluable = class_df[class_df["in_HVG3000"] == True].copy()
not_evaluable = class_df[class_df["in_HVG3000"] == False].copy()

log(f"  可评估基因 (in_HVG3000=True): n={len(evaluable)}")
for _, r in evaluable.iterrows():
    log(f"    {r['gene']}: class={r['class']}, in_oracle={r['in_oracle']}")
log(f"  不可评估基因 (in_HVG3000=False): n={len(not_evaluable)}")
for _, r in not_evaluable.iterrows():
    log(f"    {r['gene']}: class={r['class']}")

# 保存 05_gene_evaluation_scope.txt
scope_lines = [
    "STEP 05A: 基因可评估性分类",
    "="*50,
    f"总基因数: {len(class_df)}",
    f"可评估基因 (in_HVG3000=True): n={len(evaluable)}",
    "  定义: 在 HVG3000 中，可在 oracle adata 中找到对应表达值",
    "  class C (in_HVG, not TF) 或 class E (not in HVG, not TF)",
    "  注意: 所有 23 基因均为 class C 或 E（非 TF），无 strict/supported/indirect 分类",
    "",
    "可评估基因列表:",
]
for _, r in evaluable.iterrows():
    scope_lines.append(f"  {r['gene']}: class={r['class']}, "
                       f"MAFB_target_IFN_TAM={r['MAFB_target_IFN_TAM_filterA']}, "
                       f"MAFB_target_LA_TAM={r['MAFB_target_LA_TAM_filterA']}, "
                       f"MAFB_target_Mono={r['MAFB_target_Mono_filterA']}")
scope_lines += [
    "",
    "不可评估基因列表:",
]
for _, r in not_evaluable.iterrows():
    scope_lines.append(f"  {r['gene']}: class={r['class']}")

scope_path = os.path.join(OUTDIR_05, "05_gene_evaluation_scope.txt")
with open(scope_path, "w") as f:
    f.write("\n".join(scope_lines) + "\n")
log(f"  [PASS] 05_gene_evaluation_scope.txt saved")

# ── 加载 WT 和 KO oracle ─────────────────────────────────────
log(f"\n--- 加载 WT 和 KO oracle ---")
with open(os.path.join(OUTDIR, "03_perturbation/WT/oracle_WT.pkl"), "rb") as f:
    oracle_wt = pickle.load(f)
with open(os.path.join(OUTDIR, "03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl"), "rb") as f:
    oracle_ko = pickle.load(f)

log(f"  WT adata.shape: {oracle_wt.adata.shape}")
log(f"  KO adata.shape: {oracle_ko.adata.shape}")

# 对齐检查
wt_ids = list(oracle_wt.adata.obs_names)
ko_ids = list(oracle_ko.adata.obs_names)
if wt_ids != ko_ids:
    log(f"  [FAIL] WT 和 KO obs_names 不一致！")
    sys.exit(1)
log(f"  [PASS] WT 和 KO obs_names 完全一致 (n={len(wt_ids)})")

# 提取 subtype 信息
subtypes = oracle_wt.adata.obs["mac_subtype"].values
log(f"  subtype 分布: {dict(zip(*np.unique(subtypes, return_counts=True)))}")

# ── 05B: 基因层面统计 ─────────────────────────────────────────
log(f"\n--- 05B: WT imputed_count vs KO simulated_count 基因层面统计 ---")
log(f"  WT 数据: imputed_count layer")
log(f"  KO 数据: simulated_count layer")
log(f"  统计: Wilcoxon signed-rank (paired) + BH 校正 + log2FC + Cohen's d")

# 提取表达矩阵
wt_mat = oracle_wt.adata.layers["imputed_count"]   # (272, 3000)
ko_mat = oracle_ko.adata.layers["simulated_count"] # (272, 3000)
var_names = list(oracle_wt.adata.var_names)

# 转为 dense array（若为稀疏矩阵）
import scipy.sparse as sp
if sp.issparse(wt_mat): wt_mat = wt_mat.toarray()
if sp.issparse(ko_mat): ko_mat = ko_mat.toarray()

log(f"  wt_mat shape: {wt_mat.shape}, ko_mat shape: {ko_mat.shape}")

def cohens_d_paired(a, b):
    diff = a - b
    if diff.std() == 0:
        return np.nan
    return diff.mean() / diff.std()

NUMERIC_INSTABILITY_THRESHOLD = 10.0  # |Cohen's d| > 10 且 std 极小

rows = []
subtypes_list = ["Mono", "IFN_TAM", "LA_TAM"]

for gene in evaluable["gene"].values:
    if gene not in var_names:
        log(f"  [WARN] {gene} not in var_names — 跳过")
        continue
    gene_idx = var_names.index(gene)
    gene_row = class_df[class_df["gene"] == gene].iloc[0]

    for st in subtypes_list:
        mask = (subtypes == st)
        n = mask.sum()
        wt_vals = wt_mat[mask, gene_idx]
        ko_vals = ko_mat[mask, gene_idx]

        # log2FC = log2(mean_KO + 1) - log2(mean_WT + 1)
        log2fc = np.log2(ko_vals.mean() + 1) - np.log2(wt_vals.mean() + 1)

        # Wilcoxon signed-rank
        try:
            stat, pval = stats.wilcoxon(wt_vals, ko_vals, alternative='two-sided')
        except ValueError as e:
            # 若所有差值为0
            pval = 1.0
            stat = 0.0

        # Cohen's d
        cd = cohens_d_paired(ko_vals, wt_vals)
        numeric_instability = False
        if not np.isnan(cd) and abs(cd) > NUMERIC_INSTABILITY_THRESHOLD:
            diff = ko_vals - wt_vals
            if diff.std() < 1e-6:
                numeric_instability = True

        rows.append({
            "gene":                    gene,
            "subtype":                 st,
            "class_t2000":             gene_row["class"],
            "class_t5000":             gene_row["class"],  # same classification
            "evaluated_or_not":        "evaluated",
            "MAFB_target_filterA":     int(gene_row[f"MAFB_target_{st}_filterA"]),
            "log2FC":                  round(log2fc, 6),
            "p_value":                 pval,
            "padj":                    np.nan,
            "cohens_d":                round(cd, 6) if not np.isnan(cd) else np.nan,
            "n_cells":                 n,
            "numeric_instability_flag": numeric_instability,
        })

# 不可评估基因也加入结果表（标注 not_evaluated）
for gene in not_evaluable["gene"].values:
    gene_row = class_df[class_df["gene"] == gene].iloc[0]
    for st in subtypes_list:
        rows.append({
            "gene":                    gene,
            "subtype":                 st,
            "class_t2000":             gene_row["class"],
            "class_t5000":             gene_row["class"],
            "evaluated_or_not":        "not_evaluated",
            "MAFB_target_filterA":     0,
            "log2FC":                  np.nan,
            "p_value":                 np.nan,
            "padj":                    np.nan,
            "cohens_d":                np.nan,
            "n_cells":                 np.nan,
            "numeric_instability_flag": False,
        })

# BH 校正（仅对 evaluated 基因）
eval_mask = [r["evaluated_or_not"] == "evaluated" for r in rows]
eval_pvals = [rows[i]["p_value"] for i, m in enumerate(eval_mask) if m]
if len(eval_pvals) > 0:
    _, padjs, _, _ = multipletests(eval_pvals, method="fdr_bh")
    j = 0
    for i, m in enumerate(eval_mask):
        if m:
            rows[i]["padj"] = round(padjs[j], 6)
            j += 1

# 保存结果
result_df = pd.DataFrame(rows, columns=[
    "gene", "subtype", "class_t2000", "class_t5000", "evaluated_or_not",
    "MAFB_target_filterA", "log2FC", "p_value", "padj", "cohens_d",
    "n_cells", "numeric_instability_flag"
])
out_csv = os.path.join(OUTDIR_05, "05_scenic_target_stats.csv")
result_df.to_csv(out_csv, index=False)
log(f"  [PASS] 05_scenic_target_stats.csv saved: {out_csv}")

# ── 打印摘要 ──────────────────────────────────────────────────
log(f"\n--- 结果摘要 ---")
log(f"  总行数: {len(result_df)}")
eval_df = result_df[result_df["evaluated_or_not"] == "evaluated"]
log(f"  可评估基因×亚群: {len(eval_df)}")
sig_df = eval_df[eval_df["padj"] < 0.05]
log(f"  padj < 0.05: {len(sig_df)}")
instab_df = eval_df[eval_df["numeric_instability_flag"] == True]
log(f"  numeric_instability_flag=True: {len(instab_df)}")

log(f"\n  可评估基因统计（按亚群）：")
log(f"  {'基因':<12} {'亚群':<10} {'log2FC':>8} {'p_value':>10} {'padj':>10} {'Cohen_d':>8} {'MAFB_target':>12} {'instab'}")
log(f"  {'-'*80}")
for _, r in eval_df.sort_values(["subtype","gene"]).iterrows():
    sig = "**" if r['padj'] < 0.01 else ("*" if r['padj'] < 0.05 else "ns")
    log(f"  {r['gene']:<12} {r['subtype']:<10} {r['log2FC']:>8.4f} "
        f"{r['p_value']:>10.4e} {r['padj']:>10.4e} "
        f"{r['cohens_d']:>8.3f} {r['MAFB_target_filterA']:>12} "
        f"{'YES' if r['numeric_instability_flag'] else 'no':>6} {sig}")

# ── 保存日志 ──────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 05 COMPLETE. Total: {total:.1f}s")
log("="*60)

log_path = os.path.join(OUTDIR_05, "05_scenic_target_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 05_scenic_target_log.txt saved")
