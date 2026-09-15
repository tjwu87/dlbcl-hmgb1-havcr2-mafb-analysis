# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 07: 富集分析（严格背景定义）
07A-B: Regulon-based enrichment
  背景: threshold=2000 下 MAFB direct target 全部基因（各亚群）
  前景A: STEP 05 中 KO 后 padj < 0.05 的已评估基因
  前景B: threshold=2000 direct targets 中 |log2FC| > 0.3 的基因
07C: 全局 KO 响应基因集富集分析
  前景: HVG3000 中 MAFB KO vs WT padj < 0.05 且 |median_log2FC| > 0.3
  背景: HVG3000 全部基因
  方法: ORA (gseapy enrichr) + GSEA preranked
  数据库: GO_Biological_Process_2023, KEGG_2021_Human, MSigDB_Hallmark_2020
输出：
  07_enrichment/07_regulon_based_enrichment_results.csv
  07_enrichment/07_regulon_based_enrichment_log.txt
  07_enrichment/07_ko_response_enrichment_results.csv
  07_enrichment/07_ko_response_enrichment_log.txt
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
log("STEP 07: 富集分析")
log(f"Start: {t0.isoformat()}")
log("="*60)

GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "MSigDB_Hallmark_2020"]
SUBTYPES  = ["Mono", "IFN_TAM", "LA_TAM"]

# ── 加载必要文件 ──────────────────────────────────────────────
log(f"\n--- 加载必要文件 ---")

# 05 scenic target stats
scenic_stats = pd.read_csv(os.path.join(OUTDIR, "05_targets", "05_scenic_target_stats.csv"))
log(f"  scenic_target_stats: {len(scenic_stats)} rows")

# filter_A MAFB edges
mafb_edges_path = os.path.join(OUTDIR, "02_oracle_grn", "filter_A_t2000", "MAFB_edges_filterA.csv")
mafb_edges = pd.read_csv(mafb_edges_path)
log(f"  MAFB_edges_filterA: {len(mafb_edges)} rows, cols={list(mafb_edges.columns)}")

# HVG3000 var_names
with open(os.path.join(OUTDIR, "03_perturbation/WT/oracle_WT.pkl"), "rb") as f:
    oracle_wt = pickle.load(f)
with open(os.path.join(OUTDIR, "03_perturbation/MAFB_KO/oracle_MAFB_KO.pkl"), "rb") as f:
    oracle_ko = pickle.load(f)

var_names = list(oracle_wt.adata.var_names)
subtypes  = oracle_wt.adata.obs["mac_subtype"].values
log(f"  HVG3000 var_names: {len(var_names)}")

import scipy.sparse as sp
wt_imp = oracle_wt.adata.layers["imputed_count"]
ko_sim = oracle_ko.adata.layers["simulated_count"]
if sp.issparse(wt_imp): wt_imp = wt_imp.toarray()
if sp.issparse(ko_sim): ko_sim = ko_sim.toarray()

# ── 07A-B: Regulon-based enrichment ──────────────────────────
log(f"\n{'='*60}")
log(f"07A-B: Regulon-based enrichment")
log(f"{'='*60}")

# 背景基因集：各亚群 threshold=2000 下 MAFB direct target 全部基因
# 前景A: STEP 05 padj < 0.05 的已评估基因
# 前景B: threshold=2000 direct targets 中 |log2FC| > 0.3

# 检查 MAFB_edges_filterA.csv 结构
log(f"  MAFB_edges_filterA columns: {list(mafb_edges.columns)}")
log(f"  MAFB_edges_filterA head:\n{mafb_edges.head()}")

regulon_results = []
regulon_log = []

for st in SUBTYPES:
    log(f"\n  --- {st} ---")
    # 背景：该亚群的 MAFB direct targets (filter A)
    if "cluster" in mafb_edges.columns:
        bg_genes = list(mafb_edges[mafb_edges["cluster"] == st]["target"].unique())
    elif "subtype" in mafb_edges.columns:
        bg_genes = list(mafb_edges[mafb_edges["subtype"] == st]["target"].unique())
    else:
        # 尝试从 filtered_A_{st}.csv 读取
        fa_path = os.path.join(OUTDIR, "02_oracle_grn", "filter_A_t2000", f"filtered_A_{st}.csv")
        if os.path.exists(fa_path):
            fa_df = pd.read_csv(fa_path)
            log(f"    filtered_A_{st}.csv cols: {list(fa_df.columns)}")
            # 找 target 列
            target_col = [c for c in fa_df.columns if "target" in c.lower() or "gene" in c.lower()]
            if target_col:
                bg_genes = list(fa_df[target_col[0]].unique())
            else:
                bg_genes = []
        else:
            bg_genes = []

    log(f"  背景基因 (MAFB direct targets, filter A): n={len(bg_genes)}")
    if len(bg_genes) == 0:
        log(f"  [WARN] 背景基因集为空，跳过 {st}")
        regulon_log.append(f"{st}: 背景基因集为空，跳过")
        continue

    # 前景A: STEP 05 padj < 0.05 的已评估基因（该亚群）
    st_scenic = scenic_stats[(scenic_stats["subtype"] == st) &
                              (scenic_stats["evaluated_or_not"] == "evaluated") &
                              (scenic_stats["padj"] < 0.05)]
    fg_A = list(st_scenic["gene"].values)
    log(f"  前景A (STEP05 padj<0.05): n={len(fg_A)}, genes={fg_A}")

    # 前景B: threshold=2000 direct targets 中 |log2FC| > 0.3
    st_scenic_all = scenic_stats[(scenic_stats["subtype"] == st) &
                                  (scenic_stats["evaluated_or_not"] == "evaluated")]
    # 只看 MAFB_target_filterA=1 的基因
    mafb_target_genes = st_scenic_all[st_scenic_all["MAFB_target_filterA"] == 1]
    fg_B = list(mafb_target_genes[mafb_target_genes["log2FC"].abs() > 0.3]["gene"].values)
    log(f"  前景B (MAFB target |log2FC|>0.3): n={len(fg_B)}, genes={fg_B}")

    for fg_label, fg_genes in [("foreground_A", fg_A), ("foreground_B", fg_B)]:
        if len(fg_genes) < 3:
            msg = f"  {st} {fg_label}: 前景集 n={len(fg_genes)} < 3，跳过 ORA"
            log(msg)
            regulon_log.append(msg)
            if len(fg_genes) > 0:
                note = "前景集过小，ORA 结果仅供参考，不得作为论文主文证据（supplementary exploratory）"
                log(f"    注意: {note}")
                regulon_log.append(f"    {note}")
            continue

        if len(fg_genes) < 15:
            log(f"  [WARN] {st} {fg_label}: 前景集 n={len(fg_genes)} < 15，ORA 结果仅供参考，supplementary exploratory")
            regulon_log.append(f"{st} {fg_label}: n={len(fg_genes)} < 15，exploratory only")

        for gs in GENE_SETS:
            try:
                enr = gp.enrichr(
                    gene_list=fg_genes,
                    gene_sets=gs,
                    background=bg_genes,
                    outdir=None,
                    verbose=False,
                )
                res = enr.results
                if res is not None and len(res) > 0:
                    res = res.copy()
                    res["subtype"]    = st
                    res["fg_type"]    = fg_label
                    res["gene_set"]   = gs
                    res["n_fg_genes"] = len(fg_genes)
                    res["n_bg_genes"] = len(bg_genes)
                    regulon_results.append(res)
                    sig = res[res["Adjusted P-value"] < 0.05]
                    log(f"    {gs}: {len(res)} terms, {len(sig)} padj<0.05")
                else:
                    log(f"    {gs}: no results")
            except Exception as e:
                log(f"    {gs}: ERROR - {e}")
                regulon_log.append(f"{st} {fg_label} {gs}: ERROR - {e}")

# 保存 regulon-based 结果
if regulon_results:
    reg_df = pd.concat(regulon_results, ignore_index=True)
    out_reg = os.path.join(OUTDIR_07, "07_regulon_based_enrichment_results.csv")
    reg_df.to_csv(out_reg, index=False)
    log(f"\n  [PASS] 07_regulon_based_enrichment_results.csv saved: {len(reg_df)} rows")
else:
    log(f"\n  [WARN] 无 regulon-based 富集结果（前景集过小）")
    pd.DataFrame().to_csv(os.path.join(OUTDIR_07, "07_regulon_based_enrichment_results.csv"), index=False)

reg_log_path = os.path.join(OUTDIR_07, "07_regulon_based_enrichment_log.txt")
with open(reg_log_path, "w") as f:
    f.write("STEP 07A-B: Regulon-based enrichment log\n")
    f.write("="*50 + "\n")
    f.write("\n".join(regulon_log) + "\n")
log(f"  [PASS] 07_regulon_based_enrichment_log.txt saved")

# ── 07C: 全局 KO 响应基因集富集分析 ──────────────────────────
log(f"\n{'='*60}")
log(f"07C: 全局 KO 响应基因集富集分析")
log(f"{'='*60}")
log(f"  前景: HVG3000 中 MAFB KO vs WT padj < 0.05 且 |median_log2FC| > 0.3")
log(f"  背景: HVG3000 全部基因 (n={len(var_names)})")

ko_response_results = []
ko_response_log = []

for st in SUBTYPES:
    log(f"\n  --- {st} ---")
    mask = (subtypes == st)
    n = mask.sum()

    # 对每个 HVG3000 基因计算 WT imputed_count vs KO simulated_count
    log(f"  计算 {n} 个细胞的 {len(var_names)} 个基因统计...")
    wt_vals = wt_imp[mask, :]  # (n, 3000)
    ko_vals = ko_sim[mask, :]  # (n, 3000)

    # 批量计算 log2FC 和 Wilcoxon
    log2fc_arr = np.log2(ko_vals.mean(axis=0) + 1) - np.log2(wt_vals.mean(axis=0) + 1)

    pvals = []
    for i in range(len(var_names)):
        wt_g = wt_vals[:, i]
        ko_g = ko_vals[:, i]
        diff = ko_g - wt_g
        if diff.std() == 0:
            pvals.append(1.0)
        else:
            try:
                _, p = stats.wilcoxon(wt_g, ko_g, alternative='two-sided')
                pvals.append(p)
            except:
                pvals.append(1.0)

    pvals = np.array(pvals)
    _, padjs, _, _ = multipletests(pvals, method="fdr_bh")

    # 前景集
    fg_mask = (padjs < 0.05) & (np.abs(log2fc_arr) > 0.3)
    fg_genes = [var_names[i] for i in range(len(var_names)) if fg_mask[i]]
    log(f"  前景集 (padj<0.05, |log2FC|>0.3): n={len(fg_genes)}")

    if len(fg_genes) < 10:
        msg = f"  {st}: 前景集 n={len(fg_genes)} < 10，如实记录，不强行解释"
        log(msg)
        ko_response_log.append(msg)
        if len(fg_genes) == 0:
            continue

    if len(fg_genes) < 15:
        log(f"  [WARN] {st}: 前景集 n={len(fg_genes)} < 15，ORA 结果仅供参考，supplementary exploratory")
        ko_response_log.append(f"{st}: n={len(fg_genes)} < 15，exploratory only")

    # ORA
    log(f"  运行 ORA...")
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
                res["analysis"]   = "ORA"
                res["gene_set"]   = gs
                res["n_fg_genes"] = len(fg_genes)
                res["n_bg_genes"] = len(var_names)
                ko_response_results.append(res)
                sig = res[res["Adjusted P-value"] < 0.05]
                log(f"    ORA {gs}: {len(res)} terms, {len(sig)} padj<0.05")
            else:
                log(f"    ORA {gs}: no results")
        except Exception as e:
            log(f"    ORA {gs}: ERROR - {e}")
            ko_response_log.append(f"{st} ORA {gs}: ERROR - {e}")

    # GSEA preranked（按 signed log2FC 排序）
    # 有效排序基因数 = 非零 log2FC 的基因数
    nonzero_mask = log2fc_arr != 0
    n_valid = nonzero_mask.sum()
    log(f"  有效排序基因数 (log2FC != 0): {n_valid}")

    if n_valid < 15:
        msg = f"  {st}: 有效排序基因数 {n_valid} < 15，跳过 GSEA preranked，仅运行 ORA"
        log(msg)
        ko_response_log.append(msg)
    else:
        # 构建 ranked gene list
        ranked = pd.Series(log2fc_arr, index=var_names).sort_values(ascending=False)
        ranked = ranked[ranked != 0]  # 去除零值
        log(f"  运行 GSEA preranked (n={len(ranked)})...")
        for gs in GENE_SETS:
            try:
                pre_res = gp.prerank(
                    rnk=ranked,
                    gene_sets=gs,
                    outdir=None,
                    permutation_num=100,
                    seed=42,
                    verbose=False,
                )
                res = pre_res.res2d
                if res is not None and len(res) > 0:
                    res = res.copy()
                    res["subtype"]    = st
                    res["analysis"]   = "GSEA_preranked"
                    res["gene_set"]   = gs
                    res["n_fg_genes"] = len(fg_genes)
                    res["n_bg_genes"] = len(var_names)
                    ko_response_results.append(res)
                    sig = res[res["FDR q-val"] < 0.05] if "FDR q-val" in res.columns else pd.DataFrame()
                    log(f"    GSEA {gs}: {len(res)} terms, {len(sig)} FDR<0.05")
                else:
                    log(f"    GSEA {gs}: no results")
            except Exception as e:
                log(f"    GSEA {gs}: ERROR - {e}")
                ko_response_log.append(f"{st} GSEA {gs}: ERROR - {e}")

# 保存 07C 结果
if ko_response_results:
    ko_df = pd.concat(ko_response_results, ignore_index=True)
    out_ko = os.path.join(OUTDIR_07, "07_ko_response_enrichment_results.csv")
    ko_df.to_csv(out_ko, index=False)
    log(f"\n  [PASS] 07_ko_response_enrichment_results.csv saved: {len(ko_df)} rows")
else:
    log(f"\n  [WARN] 无 KO 响应富集结果")
    pd.DataFrame().to_csv(os.path.join(OUTDIR_07, "07_ko_response_enrichment_results.csv"), index=False)

ko_log_path = os.path.join(OUTDIR_07, "07_ko_response_enrichment_log.txt")
with open(ko_log_path, "w") as f:
    f.write("STEP 07C: KO response enrichment log\n")
    f.write("="*50 + "\n")
    f.write("\n".join(ko_response_log) + "\n")
log(f"  [PASS] 07_ko_response_enrichment_log.txt saved")

# ── 保存总日志 ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 07 COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

log_path = os.path.join(OUTDIR_07, "07_enrichment_main_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 07_enrichment_main_log.txt saved")
