# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 08E: k sensitivity (knn_imputation_k = 15 / 30 / 40)
主分析已用 k=20，此处对 k=15/30/40 计算关键 readout：
  - SNR (Mono / IFN_TAM / LA_TAM)
  - Mono fate: KO vs WT delta, padj, Cohen's d
  - IFN_TAM fate: KO vs WT delta, padj
  - LA_TAM fate: KO vs WT delta, padj (supplementary)
  - Mono program score: KO vs WT delta, padj
策略: 从 02C checkpoint (get_links 后，fit_GRN 前) 加载，
      重新 knn_imputation(k=X) → fit_GRN → simulate_shift → transition_prob → fate score
      每个 k 独立运行，避免超时
输出:
  08_sensitivity/08E_k{k}_fate_scores.csv
  08_sensitivity/08E_k_sensitivity_results.csv
  08_sensitivity/08E_k_sensitivity_log.txt
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

# 从命令行获取 k 值（默认运行所有）
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--k", type=int, default=None, help="knn_imputation k value (15/30/40)")
args = parser.parse_args()
K_VALUES = [args.k] if args.k else [15, 30, 40]

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log(f"STEP 08E: k sensitivity — k={K_VALUES}")
log(f"Start: {t0.isoformat()}")
log("="*60)
log(f"  主分析 k=20（已完成），此处计算 k={K_VALUES}")
log(f"  策略: 从 02C checkpoint 加载 → knn_imputation(k=X) → fit_GRN → simulate_shift → fate score")

# 加载 program genes（用于 program score 计算）
prog_genes_path = os.path.join(OUTDIR, "06_program", "06_program_gene_50.txt")
program_genes = []
with open(prog_genes_path) as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith("#"):
            program_genes.append(line)
log(f"  Program genes: {len(program_genes)} — {program_genes}")

# 加载主分析 WT fate scores（用于比较）
wt_fate_df = pd.read_csv(os.path.join(OUTDIR, "04_fate", "cell_fate_scores_WT.csv")).set_index("cell_id")

def cohens_d_paired(a, b):
    diff = a - b
    return diff.mean() / diff.std() if diff.std() > 0 else np.nan

all_results = []

for k in K_VALUES:
    log(f"\n{'='*60}")
    log(f"  k = {k}")
    log(f"{'='*60}")

    # 从 02C checkpoint 加载（get_links 后，fit_GRN 前）
    ckpt_02C = os.path.join(OUTDIR, "02_oracle_grn", "oracle_checkpoint_02C.pkl")
    assert os.path.exists(ckpt_02C), f"FAIL: {ckpt_02C} not found"
    with open(ckpt_02C, "rb") as f:
        oracle = pickle.load(f)
    log(f"  [PASS] oracle_02C loaded: adata.shape={oracle.adata.shape} [{elapsed()}]")

    subtypes = oracle.adata.obs["mac_subtype"].values

    # knn_imputation
    log(f"  knn_imputation(k={k})...")
    n_cells = oracle.adata.shape[0]
    b_sight = min(k*8, n_cells - 1)
    b_maxl  = min(k*4, n_cells - 1)
    oracle.knn_imputation(
        n_pca_dims=50,
        k=k,
        balanced=True,
        b_sight=b_sight,
        b_maxl=b_maxl,
        n_jobs=4,
    )
    log(f"  knn_imputation params: k={k}, b_sight={b_sight}, b_maxl={b_maxl}")
    log(f"  [PASS] knn_imputation done. [{elapsed()}]")

    # fit_GRN
    log(f"  fit_GRN_for_simulation...")
    oracle.fit_GRN_for_simulation(
        alpha=10,
        use_cluster_specific_TFdict=False,
        GRN_unit="cluster",
    )
    log(f"  [PASS] fit_GRN done. [{elapsed()}]")

    # ── WT baseline: estimate_transition_prob ────────────────
    log(f"  WT: simulate_shift (empty dict = WT baseline)...")
    oracle.simulate_shift(
        perturb_condition={},
        n_propagation=3,
    )
    oracle.estimate_transition_prob(
        n_neighbors=200,
        knn_random=True,
        sampled_fraction=1,
        random_seed=42,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    tp_wt = oracle.transition_prob.copy()
    log(f"  [PASS] WT transition_prob done. [{elapsed()}]")

    # WT fate score
    la_mask = (subtypes == "LA_TAM")
    la_indices = np.where(la_mask)[0]
    wt_fate_k = tp_wt[:, la_indices].sum(axis=1)

    # WT imputed_count for program score
    import scipy.sparse as sp
    wt_imp = oracle.adata.layers["imputed_count"]
    if sp.issparse(wt_imp): wt_imp = wt_imp.toarray()
    var_names = list(oracle.adata.var_names)
    prog_indices = [var_names.index(g) for g in program_genes if g in var_names]
    wt_prog = wt_imp[:, prog_indices].mean(axis=1)

    # ── KO: simulate_shift MAFB=0 ────────────────────────────
    log(f"  KO: simulate_shift (MAFB=0)...")
    oracle.simulate_shift(
        perturb_condition={"MAFB": 0},
        n_propagation=3,
    )
    oracle.estimate_transition_prob(
        n_neighbors=200,
        knn_random=True,
        sampled_fraction=1,
        random_seed=42,
    )
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    tp_ko = oracle.transition_prob.copy()
    log(f"  [PASS] KO transition_prob done. [{elapsed()}]")

    # KO fate score
    ko_fate_k = tp_ko[:, la_indices].sum(axis=1)

    # KO simulated_count for program score
    ko_sim = oracle.adata.layers["simulated_count"]
    if sp.issparse(ko_sim): ko_sim = ko_sim.toarray()
    ko_prog = ko_sim[:, prog_indices].mean(axis=1)

    # 保存 fate scores
    fate_out = pd.DataFrame({
        "cell_id":    list(oracle.adata.obs_names),
        "subtype":    subtypes,
        "wt_fate":    wt_fate_k,
        "ko_fate":    ko_fate_k,
        "k":          k,
    })
    fate_out.to_csv(os.path.join(OUTDIR_08, f"08E_k{k}_fate_scores.csv"), index=False)
    log(f"  [PASS] 08E_k{k}_fate_scores.csv saved")

    # ── 计算关键 readout ──────────────────────────────────────
    log(f"\n  关键 readout (k={k}):")
    for st in ["Mono", "IFN_TAM", "LA_TAM"]:
        mask = (subtypes == st)
        n = mask.sum()
        wt_v = wt_fate_k[mask]
        ko_v = ko_fate_k[mask]
        snr  = ko_v.mean() / ko_v.std() if ko_v.std() > 0 else np.nan
        delta = np.median(ko_v) - np.median(wt_v)
        try:
            _, pval = stats.wilcoxon(wt_v, ko_v, alternative='two-sided')
        except:
            pval = 1.0
        cd = cohens_d_paired(ko_v, wt_v)

        # program score (Mono only for primary readout)
        prog_delta = np.nan
        prog_pval  = np.nan
        if st == "Mono":
            wt_p = wt_prog[mask]
            ko_p = ko_prog[mask]
            prog_delta = np.median(ko_p) - np.median(wt_p)
            try:
                _, prog_pval = stats.wilcoxon(wt_p, ko_p, alternative='two-sided')
            except:
                prog_pval = 1.0

        all_results.append({
            "k": k,
            "subtype": st,
            "SNR_KO": round(snr, 4) if not np.isnan(snr) else np.nan,
            "fate_delta": round(delta, 6),
            "fate_pval": pval,
            "fate_cohens_d": round(cd, 6) if not np.isnan(cd) else np.nan,
            "prog_delta": round(prog_delta, 6) if not np.isnan(prog_delta) else np.nan,
            "prog_pval": prog_pval if not np.isnan(prog_pval) else np.nan,
            "n_cells": n,
            "wt_baseline_source": "per-subtype WT transition_prob (k-specific)",
        })
        log(f"    {st}: SNR={snr:.4f}, fate_delta={delta:.4f}, p={pval:.4e}, d={cd:.3f}"
            + (f", prog_delta={prog_delta:.4f}" if st == "Mono" else ""))

    log(f"  k={k} done. [{elapsed()}]")

# ── 汇总结果 ──────────────────────────────────────────────────
log(f"\n{'='*60}")
log(f"汇总 k sensitivity 结果")
log(f"{'='*60}")

# 加入主分析 k=20 结果
main_fate = pd.read_csv(os.path.join(OUTDIR, "04_fate", "04_fate_stats.csv"))
main_prog = pd.read_csv(os.path.join(OUTDIR, "06_program", "06_program_score_stats.csv"))
import re
consist_txt = open(os.path.join(OUTDIR, "04_fate", "consistency_check.txt")).read()
snr_main = {}
for line in consist_txt.split("\n"):
    m = re.search(r"(Mono|IFN_TAM|LA_TAM): mean=[\d.]+, std=[\d.]+, SNR=([\d.]+)", line)
    if m:
        snr_main[m.group(1)] = float(m.group(2))

for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    row_fate = main_fate[(main_fate["subtype"]==st) & (main_fate["comparison"]=="KO_vs_WT")]
    row_prog = main_prog[main_prog["subtype"]==st]
    all_results.append({
        "k": 20,
        "subtype": st,
        "SNR_KO": snr_main.get(st, np.nan),
        "fate_delta": row_fate["delta"].values[0] if len(row_fate) > 0 else np.nan,
        "fate_pval":  row_fate["p_value"].values[0] if len(row_fate) > 0 else np.nan,
        "fate_cohens_d": row_fate["cohens_d"].values[0] if len(row_fate) > 0 else np.nan,
        "prog_delta": row_prog["delta"].values[0] if (st=="Mono" and len(row_prog)>0) else np.nan,
        "prog_pval":  row_prog["p_value"].values[0] if (st=="Mono" and len(row_prog)>0) else np.nan,
        "n_cells": row_fate["n_cells"].values[0] if len(row_fate) > 0 else np.nan,
        "wt_baseline_source": "per-subtype WT transition_prob (main analysis k=20)",
    })

results_df = pd.DataFrame(all_results).sort_values(["subtype","k"])
out_results = os.path.join(OUTDIR_08, "08E_k_sensitivity_results.csv")
results_df.to_csv(out_results, index=False)
log(f"  [PASS] 08E_k_sensitivity_results.csv saved")

# 打印汇总表
log(f"\n  k sensitivity 汇总（Mono fate delta 和 SNR）：")
log(f"  {'k':>4} {'亚群':<10} {'SNR_KO':>8} {'fate_delta':>12} {'fate_pval':>12} {'Cohen_d':>9} {'prog_delta':>12} {'方向稳健'}")
log(f"  {'-'*85}")
for _, r in results_df.iterrows():
    direction_ok = "YES" if r["fate_delta"] < 0 and r["subtype"] == "Mono" else \
                   ("YES" if r["fate_delta"] > 0 and r["subtype"] == "LA_TAM" else "-")
    log(f"  {int(r['k']):>4} {r['subtype']:<10} {r['SNR_KO']:>8.4f} {r['fate_delta']:>12.4f} "
        f"{r['fate_pval']:>12.4e} {r['fate_cohens_d']:>9.3f} "
        f"{str(r['prog_delta'])[:8]:>12} {direction_ok}")

# 方向稳健性判断
log(f"\n  方向稳健性判断：")
mono_results = results_df[results_df["subtype"] == "Mono"]
all_mono_down = all(mono_results["fate_delta"] < 0)
log(f"  Mono fate delta 全部为负（KO 后下降）: {'YES ✓' if all_mono_down else 'NO ✗'}")
snr_order_ok = all(
    results_df[results_df["subtype"]=="LA_TAM"]["SNR_KO"].values >
    results_df[results_df["subtype"]=="Mono"]["SNR_KO"].values
)
log(f"  LA_TAM SNR 始终高于 Mono SNR: {'YES ✓' if snr_order_ok else 'NO ✗'}")

total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 08E COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

log_path = os.path.join(OUTDIR_08, "08E_k_sensitivity_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 08E_k_sensitivity_log.txt saved")
