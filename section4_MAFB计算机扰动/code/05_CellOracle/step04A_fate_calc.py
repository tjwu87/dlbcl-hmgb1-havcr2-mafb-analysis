# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 04A: LA_TAM-directed fate score 计算
唯一正式定义：
  fate_score[i] = sum of oracle.transition_prob[i, j] for all j where mac_subtype[j] == "LA_TAM"
来源：oracle.transition_prob (272×272), oracle.adata.obs['mac_subtype']
四个条件：WT / MAFB_KO / MAFB_OE / Null
前置检查：四个条件的 obs_names 顺序必须完全一致，否则立即报错停止
输出：
  04_fate/cell_fate_scores_WT.csv   (列: cell_id, subtype, fate_score)
  04_fate/cell_fate_scores_KO.csv
  04_fate/cell_fate_scores_OE.csv
  04_fate/cell_fate_scores_Null.csv
  04_fate/04_fate_definition.txt
"""
import os, sys, pickle
import numpy as np, pandas as pd
import scipy.sparse as sp
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

OUTDIR_03  = os.path.join(OUTDIR, "03_perturbation")
OUTDIR_04  = os.path.join(OUTDIR, "04_fate")
os.makedirs(OUTDIR_04, exist_ok=True)

def log(msg):
    print(msg, flush=True)

log("="*60)
log("STEP 04A: LA_TAM-directed fate score 计算")
log(f"Start: {t0.isoformat()}")
log(f"OUTDIR: {OUTDIR}")
log("="*60)

# ── Oracle pkl 路径 ───────────────────────────────────────────
PKL_MAP = {
    "WT":      os.path.join(OUTDIR_03, "WT",      "oracle_WT.pkl"),
    "MAFB_KO": os.path.join(OUTDIR_03, "MAFB_KO", "oracle_MAFB_KO.pkl"),
    "MAFB_OE": os.path.join(OUTDIR_03, "MAFB_OE", "oracle_MAFB_OE.pkl"),
    "Null":    os.path.join(OUTDIR_03, "Null",    "oracle_Null.pkl"),
}

for cond, path in PKL_MAP.items():
    assert os.path.exists(path), f"FAIL: {path} not found"
    log(f"  [OK] {cond}: {path} ({os.path.getsize(path)//1024//1024}MB)")

# ── 加载所有 oracle ───────────────────────────────────────────
log(f"\n--- Loading oracle pkl files ---")
oracles = {}
for cond, path in PKL_MAP.items():
    with open(path, "rb") as f:
        oracles[cond] = pickle.load(f)
    log(f"  [PASS] {cond} loaded. shape={oracles[cond].adata.shape} [{elapsed()}]")

# ── 前置检查：obs_names 顺序完全一致 ─────────────────────────
log(f"\n--- 前置检查: obs_names 顺序一致性 ---")
ref_cond = "WT"
ref_obs  = list(oracles[ref_cond].adata.obs_names)
for cond, oracle in oracles.items():
    obs = list(oracle.adata.obs_names)
    if obs != ref_obs:
        log(f"  [FAIL] {cond} obs_names 与 WT 不一致！")
        log(f"    WT[:5]:   {ref_obs[:5]}")
        log(f"    {cond}[:5]: {obs[:5]}")
        # 检查是否只是顺序不同（集合相同）
        if set(obs) == set(ref_obs):
            log(f"    [INFO] 集合相同但顺序不同 — 这是严重错误，停止")
        else:
            log(f"    [INFO] 集合也不同 — 更严重错误，停止")
        sys.exit(1)
    log(f"  [PASS] {cond}: obs_names 与 WT 完全一致 (n={len(obs)})")

cell_ids = ref_obs
n_cells  = len(cell_ids)
log(f"  [PASS] 所有条件 obs_names 一致。n_cells={n_cells}")

# ── 获取 LA_TAM 列索引（从 WT oracle 读取，所有条件一致）────
log(f"\n--- 确定 LA_TAM 列索引 ---")
subtypes = oracles["WT"].adata.obs["mac_subtype"].values
la_tam_mask = (subtypes == "LA_TAM")
la_tam_idx  = np.where(la_tam_mask)[0]
log(f"  LA_TAM 细胞数: {la_tam_mask.sum()} (列索引 {la_tam_idx[:5]}...)")
log(f"  亚群分布: {dict(pd.Series(subtypes).value_counts())}")

# ── 计算 fate score ───────────────────────────────────────────
log(f"\n--- 计算 fate score ---")
fate_results = {}

for cond, oracle in oracles.items():
    tp = oracle.transition_prob
    if sp.issparse(tp):
        tp = tp.toarray()
    assert tp.shape == (n_cells, n_cells), \
        f"FAIL: {cond} transition_prob shape={tp.shape}, expected ({n_cells},{n_cells})"

    # fate_score[i] = sum of tp[i, j] for j in LA_TAM
    fate_score = tp[:, la_tam_idx].sum(axis=1)
    assert fate_score.shape == (n_cells,), f"FAIL: fate_score shape={fate_score.shape}"

    fate_results[cond] = fate_score
    log(f"  {cond}: fate_score shape={fate_score.shape}, "
        f"min={fate_score.min():.4f}, max={fate_score.max():.4f}, "
        f"mean={fate_score.mean():.4f} [{elapsed()}]")

# ── 保存 CSV ──────────────────────────────────────────────────
log(f"\n--- 保存 fate score CSV ---")
COND_FILE_MAP = {
    "WT":      "cell_fate_scores_WT.csv",
    "MAFB_KO": "cell_fate_scores_KO.csv",
    "MAFB_OE": "cell_fate_scores_OE.csv",
    "Null":    "cell_fate_scores_Null.csv",
}

for cond, fname in COND_FILE_MAP.items():
    df = pd.DataFrame({
        "cell_id":    cell_ids,
        "subtype":    subtypes,
        "fate_score": fate_results[cond],
    })
    out = os.path.join(OUTDIR_04, fname)
    df.to_csv(out, index=False)
    log(f"  [PASS] {fname} saved (n={len(df)})")

# ── 保存 fate_definition.txt ──────────────────────────────────
log(f"\n--- 保存 04_fate_definition.txt ---")
defn_lines = [
    "LA_TAM-directed Fate Score — 唯一正式定义",
    "="*60,
    "",
    "定义：",
    "  fate_score[i] = sum_{j: mac_subtype[j]=='LA_TAM'} transition_prob[i, j]",
    "",
    "实现步骤：",
    "  1. 从各条件 oracle pkl 读取 oracle.transition_prob (272×272 numpy array)",
    "  2. 从 oracle.adata.obs['mac_subtype'] 读取细胞类型标签",
    "  3. 找出所有 mac_subtype == 'LA_TAM' 的细胞列索引",
    "  4. 对每个细胞 i，对 LA_TAM 列求和得到 fate_score[i]",
    "  5. 前置检查：四个条件的 obs_names 顺序完全一致（已验证）",
    "",
    "参数来源：",
    "  - transition_prob 由 estimate_transition_prob(n_neighbors=200, knn_random=True,",
    "    sampled_fraction=1, random_seed=42) 计算",
    "  - 细胞类型标签来自 oracle.adata.obs['mac_subtype']",
    "",
    f"数据规格：",
    f"  - n_cells = {n_cells}",
    f"  - n_LA_TAM = {la_tam_mask.sum()}",
    f"  - 亚群分布: {dict(pd.Series(subtypes).value_counts())}",
    "",
    "条件说明：",
    "  WT:      真实GRN，无扰动 (perturb_condition={})",
    "  MAFB_KO: 真实GRN，MAFB=0 (perturb_condition={'MAFB': 0.0})",
    "  MAFB_OE: 真实GRN，MAFB=2 (perturb_condition={'MAFB': 2.0})",
    "  Null:    随机化GRN，MAFB=0 (perturb_condition={'MAFB': 0}, use_randomized_GRN=True)",
    "           用途：SNR 分母，测量随机网络下 MAFB KO 的背景噪声水平",
    "",
    f"生成时间: {t0.isoformat()}",
    f"脚本: scripts/step04A_fate_calc.py",
]
defn_path = os.path.join(OUTDIR_04, "04_fate_definition.txt")
with open(defn_path, "w") as f:
    f.write("\n".join(defn_lines) + "\n")
log(f"  [PASS] 04_fate_definition.txt saved")

# ── 打印每个亚群在每个条件下的 fate score 均值 ───────────────
log(f"\n{'='*60}")
log(f"每个亚群在每个条件下的 fate score 均值：")
log(f"{'亚群':<12} {'WT':>10} {'MAFB_KO':>10} {'MAFB_OE':>10} {'Null':>10}")
log(f"{'-'*52}")
for cl in ["Mono", "IFN_TAM", "LA_TAM"]:
    mask = (subtypes == cl)
    vals = {cond: fate_results[cond][mask].mean() for cond in ["WT","MAFB_KO","MAFB_OE","Null"]}
    log(f"{cl:<12} {vals['WT']:>10.4f} {vals['MAFB_KO']:>10.4f} {vals['MAFB_OE']:>10.4f} {vals['Null']:>10.4f}")

total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 04A COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log(f"45min check: {'YES' if total/60>45 else 'NOT triggered'}")
log(f"stall: NO")
log("="*60)
log(f"\n停止等待确认。")
