# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

#!/usr/bin/env python3
"""
STEP 02D: Filter links (threshold A=2000, B=5000) + MAFB edge audit
STEP 02E: 23-gene 5-class classification table
"""
import os, sys, pickle, time
import numpy as np, pandas as pd
from datetime import datetime

t0 = datetime.now()
def elapsed(): return f"{(datetime.now()-t0).total_seconds():.1f}s"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(os.path.dirname(SCRIPT_DIR), "tmp_manifest", "outdir.txt")) as f:
    OUTDIR = f.read().strip()

params_df = pd.read_csv(os.path.join(OUTDIR,"00_audit","00B_frozen_parameters.tsv"),
                        sep="\t", index_col="parameter_name")
def P(k): return params_df.loc[k,"value"]

OUTDIR_02 = os.path.join(OUTDIR, "02_oracle_grn")
T_MAIN    = int(P("filter_threshold_main"))      # 2000
T_SUPP    = int(P("filter_threshold_supported")) # 5000
P_MAIN    = float(P("filter_p_main"))            # 0.001

PYSCENIC_23 = [
    "BASP1","CEP192","CHCHD10","NAGK","PPRC1","TIMM17A","MTR","MGAT1",
    "USP12","FAM20A","AAK1","ELMO1","PICALM","NAP1L4","PPP2R5A","MRGBP",
    "APOL3","GSDMD","RASSF4","IFIT3","NAE1","PPP1R18","CTSL"
]

lines=[]; log=lines.append
log("="*60); log("STEP 02D+02E: Filter Audit + 23-gene Classification")
log(f"Start: {t0.isoformat()}"); log(f"OUTDIR: {OUTDIR}")
log(f"Params: T_MAIN={T_MAIN}, T_SUPP={T_SUPP}, P_MAIN={P_MAIN}")
log("="*60)

# ── Load links ────────────────────────────────────────────────
log(f"\n--- Loading links_unfiltered ---")
with open("/tmp/links_unfiltered.pkl","rb") as f:
    links = pickle.load(f)
clusters = list(links.links_dict.keys())
log(f"  [PASS] Clusters: {clusters}")

# ── 02D: Filter A (threshold_number=2000) ────────────────────
log(f"\n--- 02D-A: Filter A (threshold_number={T_MAIN}, p<{P_MAIN}) ---")
import celloracle as co

# Filter A
links.filter_links(p=P_MAIN, weight="coef_abs", threshold_number=T_MAIN)
log(f"  [PASS] filter_links(p={P_MAIN}, threshold_number={T_MAIN}) done [{elapsed()}]")

filter_a_dir = os.path.join(OUTDIR_02, "filter_A_t2000")
os.makedirs(filter_a_dir, exist_ok=True)

for cl in clusters:
    df = links.filtered_links[cl]
    out = os.path.join(filter_a_dir, f"filtered_A_{cl}.csv")
    df.to_csv(out, index=False)
    log(f"    {cl}: {len(df)} edges → {out}")

# MAFB audit for Filter A
log(f"\n  MAFB edge audit (Filter A, t={T_MAIN}):")
mafb_a = {}
for cl in clusters:
    df = links.filtered_links[cl]
    mafb_edges = df[df['source']=='MAFB'] if 'source' in df.columns else pd.DataFrame()
    mafb_a[cl] = mafb_edges
    log(f"    {cl}: MAFB→target edges = {len(mafb_edges)}")
    if len(mafb_edges) > 0:
        top5 = mafb_edges.nlargest(5,'coef_abs')[['source','target','coef_mean','coef_abs']].to_string(index=False)
        log(f"      Top5 by coef_abs:\n{top5}")

# Save MAFB audit A
mafb_a_rows = []
for cl in clusters:
    df = mafb_a[cl].copy()
    df['cluster'] = cl
    df['filter'] = f'A_t{T_MAIN}'
    mafb_a_rows.append(df)
if mafb_a_rows:
    mafb_a_df = pd.concat(mafb_a_rows, ignore_index=True)
    mafb_a_df.to_csv(os.path.join(filter_a_dir, "MAFB_edges_filterA.csv"), index=False)
    log(f"  [PASS] MAFB edges Filter A saved")

# ── 02D: Filter B (threshold_number=5000) ────────────────────
log(f"\n--- 02D-B: Filter B (threshold_number={T_SUPP}, p<{P_MAIN}) ---")
# Reload unfiltered links for Filter B
with open("/tmp/links_unfiltered.pkl","rb") as f:
    links_b = pickle.load(f)
links_b.filter_links(p=P_MAIN, weight="coef_abs", threshold_number=T_SUPP)
log(f"  [PASS] filter_links(p={P_MAIN}, threshold_number={T_SUPP}) done [{elapsed()}]")

filter_b_dir = os.path.join(OUTDIR_02, "filter_B_t5000")
os.makedirs(filter_b_dir, exist_ok=True)

for cl in clusters:
    df = links_b.filtered_links[cl]
    out = os.path.join(filter_b_dir, f"filtered_B_{cl}.csv")
    df.to_csv(out, index=False)
    log(f"    {cl}: {len(df)} edges → {out}")

# MAFB audit for Filter B
log(f"\n  MAFB edge audit (Filter B, t={T_SUPP}):")
mafb_b = {}
for cl in clusters:
    df = links_b.filtered_links[cl]
    mafb_edges = df[df['source']=='MAFB'] if 'source' in df.columns else pd.DataFrame()
    mafb_b[cl] = mafb_edges
    log(f"    {cl}: MAFB→target edges = {len(mafb_edges)}")

# Save MAFB audit B
mafb_b_rows = []
for cl in clusters:
    df = mafb_b[cl].copy()
    df['cluster'] = cl
    df['filter'] = f'B_t{T_SUPP}'
    mafb_b_rows.append(df)
if mafb_b_rows:
    mafb_b_df = pd.concat(mafb_b_rows, ignore_index=True)
    mafb_b_df.to_csv(os.path.join(filter_b_dir, "MAFB_edges_filterB.csv"), index=False)
    log(f"  [PASS] MAFB edges Filter B saved")

# ── 02E: 23-gene 5-class classification ──────────────────────
log(f"\n--- 02E: 23-gene 5-class Classification ---")
# Load HVG list
hvg_path = os.path.join(OUTDIR, "01_preprocessing", "hvg3000_genes.txt")
with open(hvg_path) as f:
    hvg3000 = set(f.read().strip().split("\n"))
log(f"  HVG3000 size: {len(hvg3000)}")

# Load base GRN TFs
grn = pd.read_csv(translate("/mnt/results/GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv"),
                  index_col=0).reset_index()
# TFs are columns (excluding index/target column)
# The GRN has rows=target genes, cols=TF names
grn_raw = pd.read_csv(translate("/mnt/results/GSE182434/celloracle_rerun_20260324_035651/base_GRN_human_promoter.csv"),
                      index_col=0)
tf_set = set(grn_raw.columns.tolist())
log(f"  Base GRN TFs: {len(tf_set)}")

# Get all genes in oracle (after import_TF_data)
with open("/tmp/oracle_after_02C.pkl","rb") as f:
    oracle = pickle.load(f)
oracle_genes = set(oracle.adata.var_names.tolist())
log(f"  Oracle genes (HVG∩GRN): {len(oracle_genes)}")

# 5-class classification
rows = []
for gene in PYSCENIC_23:
    in_hvg = gene in hvg3000
    in_tf  = gene in tf_set
    in_oracle = gene in oracle_genes
    
    # Class assignment
    if in_hvg and in_tf and in_oracle:
        cls = "A_hvg_tf_oracle"
    elif in_hvg and in_tf and not in_oracle:
        cls = "B_hvg_tf_noOracle"
    elif in_hvg and not in_tf:
        cls = "C_hvg_notTF"
    elif not in_hvg and in_tf:
        cls = "D_notHVG_tf"
    else:
        cls = "E_notHVG_notTF"
    
    # MAFB edges in Filter A (as target)
    mafb_target_a = {}
    for cl in clusters:
        df = links.filtered_links[cl]
        if 'target' in df.columns:
            n = len(df[(df['source']=='MAFB') & (df['target']==gene)])
        else:
            n = 0
        mafb_target_a[cl] = n
    
    rows.append({
        'gene': gene,
        'in_HVG3000': in_hvg,
        'in_base_GRN_TF': in_tf,
        'in_oracle': in_oracle,
        'class': cls,
        'MAFB_target_IFN_TAM_filterA': mafb_target_a.get('IFN_TAM',0),
        'MAFB_target_LA_TAM_filterA': mafb_target_a.get('LA_TAM',0),
        'MAFB_target_Mono_filterA': mafb_target_a.get('Mono',0),
    })

gene_class_df = pd.DataFrame(rows)
out_class = os.path.join(OUTDIR_02, "02E_23gene_classification.csv")
gene_class_df.to_csv(out_class, index=False)
log(f"  [PASS] 23-gene classification saved: {out_class}")

# Summary
log(f"\n  Class distribution:")
for cls, grp in gene_class_df.groupby('class'):
    genes = grp['gene'].tolist()
    log(f"    {cls}: {len(genes)} genes — {genes}")

log(f"\n  Genes with MAFB→gene edge in Filter A (any cluster):")
mafb_targets = gene_class_df[
    (gene_class_df['MAFB_target_IFN_TAM_filterA']>0) |
    (gene_class_df['MAFB_target_LA_TAM_filterA']>0) |
    (gene_class_df['MAFB_target_Mono_filterA']>0)
]
log(f"    {mafb_targets['gene'].tolist()}")

# ── Save audit ────────────────────────────────────────────────
total = (datetime.now()-t0).total_seconds()
log(f"\n{'='*60}")
log(f"STEP 02D+02E COMPLETE. Total: {total:.1f}s ({total/60:.1f} min)")
log("="*60)

audit_de = os.path.join(OUTDIR_02, "02DE_filter_audit_log.txt")
with open(audit_de,"w") as f: f.write("\n".join(lines)+"\n")
print("\n".join(lines))
print(f"\n02D+02E COMPLETE — {total:.1f}s")
