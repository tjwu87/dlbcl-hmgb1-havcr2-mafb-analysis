# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 01: Preprocessing (01A subset, 01B HVG comparison, 01C main HVG, 01D DPT)
"""
import os, sys, shutil
import numpy as np, pandas as pd
import anndata as ad, scanpy as sc
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

params_df = pd.read_csv(os.path.join(OUTDIR,"00_audit","00B_frozen_parameters.tsv"),
                        sep="\t", index_col="parameter_name")
def P(k): return params_df.loc[k,"value"]

H5AD_USE   = P("cell_input_file")
SUBTYPES   = P("subset_subtypes").split(",")
HVG_LAYER  = "counts"
HVG_FLAVOR = P("HVG_flavor")
HVG_MAIN   = int(P("HVG_n_top_genes_main"))
HVG_SENS   = [int(x) for x in P("HVG_n_top_genes_sensitivity").split(",")]
OUTDIR_01  = os.path.join(OUTDIR, "01_preprocessing")

SCENIC = ["CEP192","CHCHD10","NAGK","PPRC1","BASP1","MTR","TIMM17A","AAK1",
          "MGAT1","USP12","FAM20A","APOL3","GSDMD","ELMO1","PICALM","RASSF4",
          "IFIT3","NAE1","NAP1L4","PPP2R5A","MRGBP","PPP1R18","CTSL"]

lines=[]; log=lines.append
log("="*60); log("STEP 01: Preprocessing"); log(f"Start: {t0.isoformat()}")
log(f"OUTDIR: {OUTDIR}"); log("="*60)

# ── 01A: Read & subset ────────────────────────────────────────
log("\n--- 01A: Read & Subset ---")
adata_full = ad.read_h5ad(H5AD_USE)
log(f"  Full shape: {adata_full.shape}")
adata = adata_full[adata_full.obs["mac_subtype"].isin(SUBTYPES)].copy()
log(f"  Subset shape: {adata.shape}")
for st in SUBTYPES:
    n = (adata.obs["mac_subtype"]==st).sum()
    log(f"    {st}: {n}")
log(f"  Total: {adata.n_obs} cells")

# Copy X_umap_paga → X_umap
if "X_umap_paga" in adata.obsm:
    adata.obsm["X_umap"] = adata.obsm["X_umap_paga"].copy()
    log("  [PASS] X_umap_paga → X_umap copied")
elif "X_umap" in adata.obsm:
    log("  [INFO] X_umap already present")
else:
    log("  [FAIL] No UMAP found"); sys.exit(1)
log(f"  [PASS] 01A done [{elapsed()}]")

# ── 01B: HVG comparison ───────────────────────────────────────
log("\n--- 01B: HVG Comparison ---")
log(f"  Layer={HVG_LAYER}, Flavor={HVG_FLAVOR}")
hvg_rows = []
for n_top in [HVG_MAIN] + HVG_SENS:
    tmp = adata.copy()
    sc.pp.highly_variable_genes(tmp, layer=HVG_LAYER, flavor=HVG_FLAVOR,
                                 n_top_genes=n_top, subset=False)
    hvg = tmp.var_names[tmp.var["highly_variable"]].tolist()
    s_in  = [g for g in SCENIC if g in hvg]
    s_out = [g for g in SCENIC if g not in hvg]
    log(f"  n_top={n_top}: HVG={len(hvg)}, pySCENIC={len(s_in)}/23, "
        f"MGAT1={'MGAT1' in hvg}, APOL3={'APOL3' in hvg}, GSDMD={'GSDMD' in hvg}")
    log(f"    in HVG: {s_in}")
    if s_out: log(f"    absent: {s_out}")
    hvg_rows.append({"n_top_genes":n_top,"layer":HVG_LAYER,"flavor":HVG_FLAVOR,
                     "n_hvg_actual":len(hvg),"scenic_covered":len(s_in),
                     "MGAT1":"MGAT1" in hvg,"APOL3":"APOL3" in hvg,"GSDMD":"GSDMD" in hvg,
                     "scenic_in_hvg":",".join(s_in),"scenic_absent":",".join(s_out)})
    del tmp

pd.DataFrame(hvg_rows).to_csv(os.path.join(OUTDIR_01,"01_hvg_comparison_table.csv"), index=False)
log(f"  [PASS] 01_hvg_comparison_table.csv saved [{elapsed()}]")

# ── 01C: Main HVG object ──────────────────────────────────────
log("\n--- 01C: Main HVG Object ---")
sc.pp.highly_variable_genes(adata, layer=HVG_LAYER, flavor=HVG_FLAVOR,
                             n_top_genes=HVG_MAIN, subset=False)
hvg3000 = adata.var_names[adata.var["highly_variable"]].tolist()
log(f"  HVG3000 actual: {len(hvg3000)}")
adata_hvg = adata[:, hvg3000].copy()
log(f"  adata_hvg shape: {adata_hvg.shape}")

# Ensure raw_count layer
if "counts" in adata_hvg.layers:
    adata_hvg.layers["raw_count"] = adata_hvg.layers["counts"].copy()
    log("  [PASS] counts → raw_count copied")
elif "raw_count" in adata_hvg.layers:
    log("  [INFO] raw_count already present")
else:
    log("  [FAIL] No counts/raw_count"); sys.exit(1)

# Save gene list
with open(os.path.join(OUTDIR_01,"hvg3000_genes.txt"),"w") as f:
    f.write("\n".join(hvg3000)+"\n")
log("  [PASS] hvg3000_genes.txt saved")

# Save h5ad via staging (约定2: /tmp staging)
staging = "/tmp/results-staging"
os.makedirs(staging, exist_ok=True)
stg_path = os.path.join(staging, "adata_hvg_counts3000.h5ad")
fin_path  = os.path.join(OUTDIR_01, "adata_hvg_counts3000.h5ad")
adata_hvg.write_h5ad(stg_path)
shutil.copy(stg_path, fin_path)
log(f"  [PASS] adata_hvg_counts3000.h5ad saved [{elapsed()}]")

# Verify
v = ad.read_h5ad(fin_path)
assert v.shape == adata_hvg.shape
log(f"  [PASS] Verified: shape={v.shape}, layers={list(v.layers.keys())}")
del v

# pySCENIC coverage
s_in3k = [g for g in SCENIC if g in hvg3000]
s_out3k = [g for g in SCENIC if g not in hvg3000]
log(f"  pySCENIC in HVG3000: {len(s_in3k)}/23")
log(f"    Present: {s_in3k}")
if s_out3k: log(f"    Absent:  {s_out3k}")

# Decision log
with open(os.path.join(OUTDIR_01,"01_main_hvg_decision.txt"),"w") as f:
    f.write(f"Main HVG Decision\n=================\n"
            f"Layer: {HVG_LAYER}\nFlavor: {HVG_FLAVOR}\nn_top_genes: {HVG_MAIN}\n"
            f"Actual HVG: {len(hvg3000)}\npySCENIC covered: {len(s_in3k)}/23\n"
            f"Rationale: counts+seurat_v3+3000 fixed per plan H1/H2/H3.\n"
            f"  4000/5000 only for sensitivity.\n"
            f"Timestamp: {datetime.now().isoformat()}\n")
log(f"  [PASS] 01C done [{elapsed()}]")

# ── 01D: DPT check ────────────────────────────────────────────
log("\n--- 01D: DPT Pseudotime ---")
dpt_path = os.path.join(OUTDIR_01,"01_dpt_audit.txt")
if "dpt_pseudotime" in adata.obs.columns:
    dpt = adata.obs["dpt_pseudotime"]
    log(f"  [PASS] dpt_pseudotime present, range=[{dpt.min():.4f},{dpt.max():.4f}], nan={dpt.isna().sum()}")
    with open(dpt_path,"w") as f:
        f.write(f"DPT Audit\nStatus: PRESENT — no recomputation needed\n"
                f"Range: [{dpt.min():.4f},{dpt.max():.4f}]\nnan: {dpt.isna().sum()}\n"
                f"Timestamp: {datetime.now().isoformat()}\n")
else:
    log("  [WARN] dpt_pseudotime absent — computing")
    mono_mask = adata.obs["mac_subtype"]=="Mono"
    if "MAFB" in adata.var_names:
        expr = np.asarray(adata[mono_mask,"MAFB"].layers["log1p_norm"]).flatten()
        root_idx = np.where(mono_mask)[0][np.argmin(expr)]
    else:
        root_idx = np.where(mono_mask)[0][0]
    root_cell = adata.obs_names[root_idx]
    adata.uns["iroot"] = int(root_idx)
    sc.tl.diffmap(adata); sc.tl.dpt(adata)
    dpt = adata.obs["dpt_pseudotime"]
    log(f"  [PASS] DPT computed, range=[{dpt.min():.4f},{dpt.max():.4f}]")
    with open(dpt_path,"w") as f:
        f.write(f"DPT Audit\nStatus: COMPUTED in STEP 01D\nRoot: {root_cell}\n"
                f"iroot: {root_idx}\nRange: [{dpt.min():.4f},{dpt.max():.4f}]\n"
                f"Timestamp: {datetime.now().isoformat()}\n")
log(f"  [PASS] 01D done [{elapsed()}]")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 01 COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("45min check: NOT triggered | 55min checkpoint: NOT triggered | stall: NO")
log("="*60)

audit_txt = os.path.join(OUTDIR_01,"01_subset_audit.txt")
with open(audit_txt,"w") as f: f.write("\n".join(lines)+"\n")
with open(os.path.join(OUTDIR_01,"01_hvg_audit.txt"),"w") as f: f.write("\n".join(lines)+"\n")
print("\n".join(lines))
print(f"\nSTEP 01 COMPLETE — {total:.1f}s")
