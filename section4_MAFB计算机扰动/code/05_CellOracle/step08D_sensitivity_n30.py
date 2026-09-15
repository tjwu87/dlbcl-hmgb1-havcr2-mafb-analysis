# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 08D: 敏感性分析 — n_neighbors=30
Class 1: 复用主 Oracle 和主 fitted GRN
仅重跑 estimate_transition_prob(n_neighbors=30) + calculate_embedding_shift + fate score
禁止重新初始化 Oracle 或重新拟合 GRN
输出:
  08_sensitivity/08D_n30_fate_scores.csv
  08_sensitivity/08D_n30_fate_stats.csv
  08_sensitivity/08D_n30_log.txt
"""
import os, sys, pickle, copy
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

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
log("STEP 08D: 敏感性分析 — n_neighbors=30")
log(f"Start: {t0.isoformat()}")
log("="*60)
log(f"  策略: 复用主 Oracle + 主 fitted GRN")
log(f"  仅重跑: estimate_transition_prob(n_neighbors=30) + calculate_embedding_shift + fate score")
log(f"  禁止: 重新初始化 Oracle 或重新拟合 GRN")

# ── 加载主 Oracle (02F checkpoint，已 fit_GRN) ────────────────
log(f"\n--- 加载主 Oracle (02F checkpoint) ---")
main_pkl = os.path.join(OUTDIR, "02_oracle_grn", "oracle_checkpoint_02F.pkl")
assert os.path.exists(main_pkl), f"FAIL: {main_pkl} not found"

with open(main_pkl, "rb") as f:
    oracle = pickle.load(f)
log(f"  [PASS] oracle loaded: adata.shape={oracle.adata.shape}")
log(f"  [PASS] 主 Oracle 已包含 fit_GRN 结果，不重新拟合")

# 验证 GRN 已拟合
assert hasattr(oracle, 'coef_matrix_per_cluster'), "FAIL: coef_matrix_per_cluster not found"
log(f"  [PASS] coef_matrix_per_cluster 存在，GRN 已拟合")

subtypes = oracle.adata.obs["mac_subtype"].values
log(f"  subtype 分布: {dict(zip(*np.unique(subtypes, return_counts=True)))}")

# ── MAFB KO 扰动模拟 (复用 GRN，仅重跑 transition prob) ──────
log(f"\n--- MAFB KO 扰动模拟 (n_neighbors=30) ---")
log(f"  simulate_shift: MAFB=0")
oracle.simulate_shift(
    perturb_condition={"MAFB": 0},
    n_propagation=3,
)
log(f"  [PASS] simulate_shift done. [{elapsed()}]")

log(f"  estimate_transition_prob: n_neighbors=30, knn_random=True, sampled_fraction=1")
oracle.estimate_transition_prob(
    n_neighbors=30,
    knn_random=True,
    sampled_fraction=1,
    random_seed=42,
)
log(f"  [PASS] estimate_transition_prob done. [{elapsed()}]")

log(f"  calculate_embedding_shift: sigma_corr=0.05")
oracle.calculate_embedding_shift(sigma_corr=0.05)
log(f"  [PASS] calculate_embedding_shift done. [{elapsed()}]")

# ── 计算 LA_TAM fate score ────────────────────────────────────
log(f"\n--- 计算 LA_TAM fate score (n_neighbors=30) ---")
tp = oracle.transition_prob  # (272, 272)
la_mask = (subtypes == "LA_TAM")
la_indices = np.where(la_mask)[0]
fate_scores = tp[:, la_indices].sum(axis=1)
log(f"  transition_prob shape: {tp.shape}")
log(f"  LA_TAM indices: {len(la_indices)}")
log(f"  fate_score range: [{fate_scores.min():.4f}, {fate_scores.max():.4f}]")

# 保存 fate scores
fate_df = pd.DataFrame({
    "cell_id":    list(oracle.adata.obs_names),
    "subtype":    subtypes,
    "fate_score": fate_scores,
    "condition":  "MAFB_KO_n30",
})
out_fate = os.path.join(OUTDIR_08, "08D_n30_fate_scores.csv")
fate_df.to_csv(out_fate, index=False)
log(f"  [PASS] 08D_n30_fate_scores.csv saved")

# ── 与主分析 WT fate score 比较 ───────────────────────────────
log(f"\n--- 与主分析 WT fate score 比较 ---")
wt_fate = pd.read_csv(os.path.join(OUTDIR, "04_fate", "cell_fate_scores_WT.csv")).set_index("cell_id")

# 对齐
ref_ids = list(wt_fate.index)
ko_ids  = list(fate_df.set_index("cell_id").index)
assert ref_ids == ko_ids, "FAIL: cell_id 顺序不一致"
log(f"  [PASS] cell_id 顺序一致 (n={len(ref_ids)})")

wt_scores = wt_fate["fate_score"].values
ko_scores = fate_scores

def cohens_d_paired(a, b):
    diff = a - b
    return diff.mean() / diff.std() if diff.std() > 0 else np.nan

stat_rows = []
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    mask = (subtypes == st)
    n = mask.sum()
    wt_v = wt_scores[mask]
    ko_v = ko_scores[mask]
    med_wt = np.median(wt_v)
    med_ko = np.median(ko_v)
    delta  = med_ko - med_wt
    snr    = ko_v.mean() / ko_v.std() if ko_v.std() > 0 else np.nan
    try:
        _, pval = stats.wilcoxon(wt_v, ko_v, alternative='two-sided')
    except:
        pval = 1.0
    cd = cohens_d_paired(ko_v, wt_v)
    stat_rows.append({
        "subtype": st, "n_cells": n,
        "median_WT": round(med_wt, 6), "median_KO_n30": round(med_ko, 6),
        "delta": round(delta, 6), "delta_direction": "median_KO_n30 - median_WT",
        "SNR_KO_n30": round(snr, 4) if not np.isnan(snr) else np.nan,
        "p_value": pval, "padj": np.nan,
        "cohens_d": round(cd, 6) if not np.isnan(cd) else np.nan,
    })
    log(f"  {st}: median_WT={med_wt:.4f}, median_KO_n30={med_ko:.4f}, "
        f"delta={delta:.4f}, SNR={snr:.4f}, p={pval:.4e}, d={cd:.3f}, n={n}")

# BH 校正
pvals = [r["p_value"] for r in stat_rows]
_, padjs, _, _ = multipletests(pvals, method="fdr_bh")
for i, r in enumerate(stat_rows):
    r["padj"] = round(padjs[i], 6)

stats_df = pd.DataFrame(stat_rows)
out_stats = os.path.join(OUTDIR_08, "08D_n30_fate_stats.csv")
stats_df.to_csv(out_stats, index=False)
log(f"  [PASS] 08D_n30_fate_stats.csv saved")

# ── 与主分析 (n_neighbors=200) 对比 ──────────────────────────
log(f"\n--- 与主分析 (n_neighbors=200) 对比 ---")
main_stats = pd.read_csv(os.path.join(OUTDIR, "04_fate", "04_fate_stats.csv"))
main_ko_wt = main_stats[main_stats["comparison"] == "KO_vs_WT"].set_index("subtype")

log(f"  {'亚群':<10} {'主分析 delta':>14} {'n30 delta':>12} {'方向一致':>10} {'主分析 padj':>14} {'n30 padj':>12}")
log(f"  {'-'*75}")
for _, r in stats_df.iterrows():
    st = r["subtype"]
    main_delta = main_ko_wt.loc[st, "delta"] if st in main_ko_wt.index else np.nan
    main_padj  = main_ko_wt.loc[st, "padj"]  if st in main_ko_wt.index else np.nan
    direction_ok = "YES" if (r["delta"] * main_delta > 0) else "NO"
    log(f"  {st:<10} {main_delta:>14.4f} {r['delta']:>12.4f} {direction_ok:>10} "
        f"{main_padj:>14.4e} {r['padj']:>12.4e}")

total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 08D COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

log_path = os.path.join(OUTDIR_08, "08D_n30_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 08D_n30_log.txt saved")
