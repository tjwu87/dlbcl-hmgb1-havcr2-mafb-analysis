# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
Sub-step 04C / STEP 07C: Fig23 前景基因集准备 + ORA 富集分析
- 读取 WT imputed_count 和 KO simulated_count
- 分亚群 (Mono / IFN_TAM / LA_TAM) 做 MAFB KO vs WT Wilcoxon signed-rank (paired) + BH
- log2FC = log2((mean_KO + 1e-9) / (mean_WT + 1e-9))
- 前景集: padj < 0.05 且 |log2FC| > 0.3，分 up / down
- gseapy enrichr ORA: GO_Biological_Process_2023, KEGG_2021_Human, MSigDB_Hallmark_2020
- 背景集: HVG3000 全部基因
输出:
  07_enrichment/07_ko_response_deg_per_subtype.csv
  07_enrichment/07_ko_response_enrichment_results.csv  (覆盖旧版，含 direction 列)
  07_enrichment/07_ko_response_enrichment_log.txt
脚本: scripts/step07C_ko_response_enrichment.py
"""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp
from datetime import datetime

t0 = datetime.now()
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_07 = os.path.join(OUTDIR, "07_enrichment")
os.makedirs(OUTDIR_07, exist_ok=True)

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg, flush=True)

log("="*60)
log("STEP 07C (Sub-step 04C): Fig23 前景基因集 + ORA 富集分析")
log(f"Start: {t0.isoformat()}")
log("="*60)

GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]
SUBTYPES  = ["Mono", "IFN_TAM", "LA_TAM"]

# ── 加载 WT 和 KO oracle ─────────────────────────────────────
log(f"\n--- 加载 oracle ---")

# 优先从 03_perturbation 加载
wt_pkl  = os.path.join(OUTDIR, "03_perturbation/WT/oracle_WT.pkl")
ko_pkl  = os.path.join(OUTDIR, "03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl")

with open(wt_pkl, "rb") as f:
    oracle_wt = pickle.load(f)
with open(ko_pkl, "rb") as f:
    oracle_ko = pickle.load(f)

log(f"  WT adata.shape: {oracle_wt.adata.shape}")
log(f"  KO adata.shape: {oracle_ko.adata.shape}")

# 对齐检查
wt_ids = list(oracle_wt.adata.obs_names)
ko_ids = list(oracle_ko.adata.obs_names)
assert wt_ids == ko_ids, "FAIL: WT 和 KO obs_names 不一致"
log(f"  [PASS] obs_names 完全一致 (n={len(wt_ids)})")

import scipy.sparse as sp

# WT: imputed_count layer
wt_mat = oracle_wt.adata.layers["imputed_count"]
if sp.issparse(wt_mat): wt_mat = wt_mat.toarray()

# KO: 使用 simulated_count（KO 模拟后的表达值）
# 注意: simulation_input 与 WT imputed_count 完全相同（差值为0），不能用于比较
# simulated_count 才是 MAFB KO 扰动后的模拟结果
ko_mat = oracle_ko.adata.layers["simulated_count"]
ko_layer_used = "simulated_count"
if sp.issparse(ko_mat): ko_mat = ko_mat.toarray()

log(f"  WT layer: imputed_count")
log(f"  KO layer: {ko_layer_used}")
log(f"  wt_mat shape: {wt_mat.shape}, ko_mat shape: {ko_mat.shape}")

var_names = list(oracle_wt.adata.var_names)
subtypes  = oracle_wt.adata.obs["mac_subtype"].values
log(f"  subtype 分布: {dict(zip(*np.unique(subtypes, return_counts=True)))}")

# ── 分亚群 Wilcoxon + BH + log2FC ────────────────────────────
log(f"\n--- 分亚群 Wilcoxon signed-rank (paired) + BH + log2FC ---")
log(f"  log2FC = log2((mean_KO + 1e-9) / (mean_WT + 1e-9))")
log(f"  前景集: padj < 0.05 且 |log2FC| > 0.3")

all_deg_rows = []

for st in SUBTYPES:
    log(f"\n  --- {st} ---")
    mask = (subtypes == st)
    n = mask.sum()
    wt_sub = wt_mat[mask, :]  # (n, 3000)
    ko_sub = ko_mat[mask, :]  # (n, 3000)

    # log2FC
    mean_wt = wt_sub.mean(axis=0)
    mean_ko = ko_sub.mean(axis=0)
    log2fc  = np.log2((mean_ko + 1e-9) / (mean_wt + 1e-9))

    # Wilcoxon signed-rank (paired) per gene
    pvals = []
    for i in range(len(var_names)):
        diff = ko_sub[:, i] - wt_sub[:, i]
        if diff.std() == 0:
            pvals.append(1.0)
        else:
            try:
                _, p = stats.wilcoxon(wt_sub[:, i], ko_sub[:, i], alternative='two-sided')
                pvals.append(p)
            except:
                pvals.append(1.0)

    pvals = np.array(pvals)
    _, padjs, _, _ = multipletests(pvals, method="fdr_bh")

    # 前景集
    fg_mask = (padjs < 0.05) & (np.abs(log2fc) > 0.3)
    n_fg = fg_mask.sum()
    log(f"  n_cells={n}, 前景集 (padj<0.05, |log2FC|>0.3): n={n_fg}")

    for i, gene in enumerate(var_names):
        direction = "ns"
        if padjs[i] < 0.05 and log2fc[i] > 0.3:
            direction = "up"
        elif padjs[i] < 0.05 and log2fc[i] < -0.3:
            direction = "down"
        all_deg_rows.append({
            "gene":    gene,
            "subtype": st,
            "log2FC":  round(float(log2fc[i]), 6),
            "p_value": float(pvals[i]),
            "padj":    round(float(padjs[i]), 6),
            "direction": direction,
        })

# 保存 07_ko_response_deg_per_subtype.csv
deg_df = pd.DataFrame(all_deg_rows, columns=["gene","subtype","log2FC","p_value","padj","direction"])
out_deg = os.path.join(OUTDIR_07, "07_ko_response_deg_per_subtype.csv")
deg_df.to_csv(out_deg, index=False)
log(f"\n  [PASS] 07_ko_response_deg_per_subtype.csv saved: {len(deg_df)} rows")

# 摘要
for st in SUBTYPES:
    sub = deg_df[deg_df["subtype"] == st]
    n_up   = (sub["direction"] == "up").sum()
    n_down = (sub["direction"] == "down").sum()
    log(f"  {st}: up={n_up}, down={n_down}, total_fg={n_up+n_down}")

# ── ORA 富集分析 ──────────────────────────────────────────────
log(f"\n--- ORA 富集分析 (gseapy enrichr) ---")
log(f"  数据库: {GENE_SETS}")
log(f"  背景集: HVG3000 全部基因 (n={len(var_names)})")

enr_results = []
enr_log = []

for st in SUBTYPES:
    sub = deg_df[deg_df["subtype"] == st]
    for direction in ["up", "down"]:
        fg_genes = list(sub[sub["direction"] == direction]["gene"].values)
        n_fg = len(fg_genes)
        log(f"\n  {st} {direction}: n_fg={n_fg}")

        if n_fg < 10:
            msg = f"  {st} {direction}: 前景集 n={n_fg} < 10，如实记录，不强行解释"
            log(msg)
            enr_log.append(msg)
            if n_fg == 0:
                continue

        if n_fg < 15:
            note = f"  {st} {direction}: n={n_fg} < 15，ORA 结果标注 exploratory，不作为主文证据"
            log(note)
            enr_log.append(note)

        for gs in GENE_SETS:
            try:
                enr = gp.enrichr(
                    gene_list=fg_genes,
                    gene_sets=gs,
                    background=var_names,
                    outdir=None,
                    verbose=False,
                )
                res = enr.results
                if res is not None and len(res) > 0:
                    res = res.copy()
                    res["subtype"]    = st
                    res["direction"]  = direction
                    res["gene_set"]   = gs
                    res["n_fg_genes"] = n_fg
                    res["n_bg_genes"] = len(var_names)
                    res["exploratory"] = (n_fg < 15)
                    enr_results.append(res)
                    sig = res[res["Adjusted P-value"] < 0.05]
                    log(f"    {gs}: {len(res)} terms, {len(sig)} padj<0.05")
                else:
                    log(f"    {gs}: no results")
            except Exception as e:
                log(f"    {gs}: ERROR - {e}")
                enr_log.append(f"{st} {direction} {gs}: ERROR - {e}")

# 保存 07_ko_response_enrichment_results.csv（覆盖旧版，含 direction 列）
if enr_results:
    enr_df = pd.concat(enr_results, ignore_index=True)
    out_enr = os.path.join(OUTDIR_07, "07_ko_response_enrichment_results.csv")
    enr_df.to_csv(out_enr, index=False)
    log(f"\n  [PASS] 07_ko_response_enrichment_results.csv saved: {len(enr_df)} rows")
    # 打印显著结果摘要
    sig_df = enr_df[enr_df["Adjusted P-value"] < 0.05]
    log(f"  显著富集 (padj<0.05): {len(sig_df)} 条目")
    for st in SUBTYPES:
        for direction in ["up", "down"]:
            sub_sig = sig_df[(sig_df["subtype"]==st) & (sig_df["direction"]==direction)]
            if len(sub_sig) > 0:
                log(f"    {st} {direction}: {len(sub_sig)} 显著条目")
                for _, r in sub_sig.head(3).iterrows():
                    log(f"      {r['gene_set']}: {r['Term'][:60]} (padj={r['Adjusted P-value']:.3e})")
else:
    log(f"\n  [WARN] 无 ORA 富集结果")
    pd.DataFrame().to_csv(os.path.join(OUTDIR_07, "07_ko_response_enrichment_results.csv"), index=False)

# 保存日志
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 07C COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

enr_log_path = os.path.join(OUTDIR_07, "07_ko_response_enrichment_log.txt")
with open(enr_log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 07_ko_response_enrichment_log.txt saved")
