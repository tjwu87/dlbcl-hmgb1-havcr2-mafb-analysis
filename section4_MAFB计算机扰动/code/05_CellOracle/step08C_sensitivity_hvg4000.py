# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 08C: 敏感性分析 — HVG=4000 (Class 2 rebuild)
因 HVG 改为 4000，允许重建 sensitivity-specific Oracle
执行最小流程: HVG4000 → Oracle → knn_imputation(k=20) → get_links → fit_GRN → simulate_shift → fate score
所有对象明确标注 sensitivity-specific，不覆盖主 Oracle
输出:
  08_sensitivity/08C_hvg4000_fate_scores.csv
  08_sensitivity/08C_hvg4000_fate_stats.csv
  08_sensitivity/08C_hvg4000_log.txt
"""
import os, sys, pickle
import numpy as np, pandas as pd
import scanpy as sc
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
log("STEP 08C: 敏感性分析 — HVG=4000 (Class 2 rebuild)")
log(f"Start: {t0.isoformat()}")
log("="*60)
log(f"  策略: 重建 sensitivity-specific Oracle (HVG4000)")
log(f"  标注: 所有对象明确为 sensitivity-specific，不覆盖主 Oracle")

import celloracle as co

# ── 加载原始 adata ────────────────────────────────────────────
log(f"\n--- 加载原始 adata ---")
adata_path = "/tmp/adata_mac_annotated_fixed.h5ad"
assert os.path.exists(adata_path), f"FAIL: {adata_path} not found"
adata_full = sc.read_h5ad(adata_path)
log(f"  adata_full.shape: {adata_full.shape}")

# 筛选 mac 细胞（与主分析相同）
mac_subtypes = ["Mono", "IFN_TAM", "LA_TAM"]
adata_mac = adata_full[adata_full.obs["mac_subtype"].isin(mac_subtypes)].copy()
log(f"  adata_mac.shape: {adata_mac.shape}")
log(f"  subtype 分布: {dict(adata_mac.obs['mac_subtype'].value_counts())}")

# ── HVG 4000 (counts layer, seurat_v3) ───────────────────────
log(f"\n--- HVG 4000 (counts layer, seurat_v3) ---")
adata_hvg = adata_mac.copy()
adata_hvg.X = adata_hvg.layers["counts"].copy()
sc.pp.highly_variable_genes(
    adata_hvg,
    n_top_genes=4000,
    flavor="seurat_v3",
    layer="counts",
)
hvg_genes = adata_hvg.var_names[adata_hvg.var["highly_variable"]].tolist()
log(f"  HVG4000: n={len(hvg_genes)}")

# 子集到 HVG4000
adata_hvg4000 = adata_mac[:, hvg_genes].copy()
# 设置 X = raw_count（celloracle 要求）
adata_hvg4000.X = adata_hvg4000.layers["counts"].copy()
log(f"  adata_hvg4000.shape: {adata_hvg4000.shape}")

# ── 构建 sensitivity-specific Oracle ─────────────────────────
log(f"\n--- 构建 sensitivity-specific Oracle ---")
base_grn_path = translate("/mnt/results/GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv")
assert os.path.exists(base_grn_path), f"FAIL: {base_grn_path} not found"
grn = pd.read_csv(base_grn_path, index_col=0).reset_index()
log(f"  base GRN loaded: {grn.shape}")

oracle_c = co.Oracle()
oracle_c.import_anndata_as_raw_count(
    adata=adata_hvg4000,
    cluster_column_name="mac_subtype",
    embedding_name="X_umap",
)
oracle_c.import_TF_data(TF_info_matrix=grn)
log(f"  [PASS] Oracle built. adata.shape={oracle_c.adata.shape} [{elapsed()}]")

# PCA + knn_imputation
oracle_c.perform_PCA()
log(f"  [PASS] PCA done. [{elapsed()}]")

n_cells = oracle_c.adata.shape[0]
k = 20
b_sight = min(k*8, n_cells - 1)
b_maxl  = min(k*4, n_cells - 1)
oracle_c.knn_imputation(
    n_pca_dims=50,
    k=k,
    balanced=True,
    b_sight=b_sight,
    b_maxl=b_maxl,
    n_jobs=4,
)
log(f"  [PASS] knn_imputation(k={k}) done. [{elapsed()}]")

# get_links
log(f"  get_links...")
links_c = oracle_c.get_links(
    cluster_name_for_GRN_unit="mac_subtype",
    alpha=10,
    verbose_level=10,
    test_mode=False,
)
log(f"  [PASS] get_links done. [{elapsed()}]")

# filter_links (p<0.001, threshold=2000)
links_c.filter_links(
    p=0.001,
    weight="coef_abs",
    threshold_number=2000,
)
log(f"  [PASS] filter_links done. [{elapsed()}]")

# fit_GRN
oracle_c.get_cluster_specific_TFdict_from_Links(links_object=links_c)
oracle_c.fit_GRN_for_simulation(
    alpha=10,
    use_cluster_specific_TFdict=False,
    GRN_unit="cluster",
)
log(f"  [PASS] fit_GRN done. [{elapsed()}]")

# 保存 sensitivity-specific Oracle checkpoint
ckpt_c = os.path.join(OUTDIR_08, "08C_oracle_hvg4000_sensitivity.pkl")
with open(ckpt_c, "wb") as f:
    pickle.dump(oracle_c, f)
log(f"  [PASS] 08C_oracle_hvg4000_sensitivity.pkl saved [{elapsed()}]")

subtypes = oracle_c.adata.obs["mac_subtype"].values

# ── WT baseline ───────────────────────────────────────────────
log(f"\n--- WT baseline (empty dict) ---")
oracle_c.simulate_shift(perturb_condition={}, n_propagation=3)
oracle_c.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1, random_seed=42)
oracle_c.calculate_embedding_shift(sigma_corr=0.05)
tp_wt = oracle_c.transition_prob.copy()
la_mask = (subtypes == "LA_TAM")
la_indices = np.where(la_mask)[0]
wt_fate = tp_wt[:, la_indices].sum(axis=1)
log(f"  [PASS] WT done. [{elapsed()}]")

# ── KO: MAFB=0 ───────────────────────────────────────────────
log(f"\n--- KO: MAFB=0 ---")
oracle_c.simulate_shift(perturb_condition={"MAFB": 0}, n_propagation=3)
oracle_c.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1, random_seed=42)
oracle_c.calculate_embedding_shift(sigma_corr=0.05)
tp_ko = oracle_c.transition_prob.copy()
ko_fate = tp_ko[:, la_indices].sum(axis=1)
log(f"  [PASS] KO done. [{elapsed()}]")

# 保存 fate scores
fate_df = pd.DataFrame({
    "cell_id":    list(oracle_c.adata.obs_names),
    "subtype":    subtypes,
    "wt_fate":    wt_fate,
    "ko_fate":    ko_fate,
    "sensitivity": "HVG4000",
})
fate_df.to_csv(os.path.join(OUTDIR_08, "08C_hvg4000_fate_scores.csv"), index=False)
log(f"  [PASS] 08C_hvg4000_fate_scores.csv saved")

# ── 统计比较 ──────────────────────────────────────────────────
log(f"\n--- 统计比较 ---")
def cohens_d_paired(a, b):
    diff = a - b
    return diff.mean() / diff.std() if diff.std() > 0 else np.nan

stat_rows = []
for st in ["Mono", "IFN_TAM", "LA_TAM"]:
    mask = (subtypes == st)
    n = mask.sum()
    wt_v = wt_fate[mask]
    ko_v = ko_fate[mask]
    snr  = ko_v.mean() / ko_v.std() if ko_v.std() > 0 else np.nan
    delta = np.median(ko_v) - np.median(wt_v)
    try:
        _, pval = stats.wilcoxon(wt_v, ko_v, alternative='two-sided')
    except:
        pval = 1.0
    cd = cohens_d_paired(ko_v, wt_v)
    stat_rows.append({
        "subtype": st, "n_cells": n,
        "median_WT": round(np.median(wt_v), 6),
        "median_KO": round(np.median(ko_v), 6),
        "delta": round(delta, 6),
        "SNR_KO": round(snr, 4) if not np.isnan(snr) else np.nan,
        "p_value": pval, "padj": np.nan,
        "cohens_d": round(cd, 6) if not np.isnan(cd) else np.nan,
        "sensitivity": "HVG4000",
    })
    log(f"  {st}: SNR={snr:.4f}, delta={delta:.4f}, p={pval:.4e}, d={cd:.3f}, n={n}")

pvals = [r["p_value"] for r in stat_rows]
_, padjs, _, _ = multipletests(pvals, method="fdr_bh")
for i, r in enumerate(stat_rows):
    r["padj"] = round(padjs[i], 6)

stats_df = pd.DataFrame(stat_rows)
stats_df.to_csv(os.path.join(OUTDIR_08, "08C_hvg4000_fate_stats.csv"), index=False)
log(f"  [PASS] 08C_hvg4000_fate_stats.csv saved")

# 与主分析对比
log(f"\n--- 与主分析 (HVG3000) 对比 ---")
main_stats = pd.read_csv(os.path.join(OUTDIR, "04_fate", "04_fate_stats.csv"))
main_ko_wt = main_stats[main_stats["comparison"] == "KO_vs_WT"].set_index("subtype")
log(f"  {'亚群':<10} {'主分析 delta':>14} {'HVG4000 delta':>15} {'方向一致':>10}")
log(f"  {'-'*55}")
for _, r in stats_df.iterrows():
    st = r["subtype"]
    main_delta = main_ko_wt.loc[st, "delta"] if st in main_ko_wt.index else np.nan
    direction_ok = "YES" if (r["delta"] * main_delta > 0) else "NO"
    log(f"  {st:<10} {main_delta:>14.4f} {r['delta']:>15.4f} {direction_ok:>10}")

total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 08C COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

log_path = os.path.join(OUTDIR_08, "08C_hvg4000_log.txt")
with open(log_path, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"  [PASS] 08C_hvg4000_log.txt saved")
