import pandas as pd

full = pd.read_csv(r"D:\bulk-download\GSE182434\cellcomm\liana_results_full.csv")
SRC, LIG, REC = "source", "ligand_complex", "receptor_complex"   # ← receptor！
MR = "magnitude_rank"

for thr in [0.05, 0.10]:
    sub = full[full[MR] < thr].sort_values(MR).reset_index(drop=True)
    sub["priority_rank"] = range(1, len(sub)+1)
    print(f"\n=== magnitude_rank < {thr}：共 {len(sub)} 对，前 10 ===")
    print(sub[[SRC, LIG, REC, MR, "priority_rank"]].head(10).to_string(index=False))

    hit = sub[(sub[LIG].str.upper()=="HMGB1") & (sub[REC].str.upper()=="HAVCR2")]
    if len(hit):
        h = hit.iloc[0]
        print(f"\n>>> HMGB1–HAVCR2 在 <{thr} 下: priority_rank = {int(h['priority_rank'])}/{len(sub)}, magnitude_rank = {h[MR]:.4f}")
    else:
        print(f"\n>>> HMGB1–HAVCR2 不在 <{thr} 集合中")

# 存档 top 20（0.1 口径）
relaxed = full[full[MR] < 0.10].sort_values(MR).head(20).reset_index(drop=True)
relaxed["priority_rank_under_0.1"] = range(1, len(relaxed)+1)
relaxed.to_csv(r"D:\bulk-download\GSE182434\cellcomm\liana_top20_under_0.1.csv", index=False)
print("\n已存: liana_top20_under_0.1.csv")

import pandas as pd

full7 = pd.read_csv(r"D:\bulk-download\GSE182434\cellcomm\liana_results_full.csv")
full13 = pd.read_csv(r"D:\bulk-download\GSE182434\cellcomm\liana_results_13subtypes_full.csv")
COLS = ["source","target","ligand_complex","receptor_complex","magnitude_rank","specificity_rank"]

def is_hh(df):
    return (df["ligand_complex"].str.upper()=="HMGB1") & (df["receptor_complex"].str.upper()=="HAVCR2")

# ── 1) HMGB1–HAVCR2 的完整行（两张表都看，含 source/target/specificity_rank）──
print("=== 7 群表中的 HMGB1–HAVCR2 ===")
print(full7[is_hh(full7)][COLS].sort_values("magnitude_rank").to_string(index=False))
print("\n=== 13 亚型表中的 HMGB1–HAVCR2 ===")
print(full13[is_hh(full13)][COLS].sort_values("magnitude_rank").to_string(index=False))

# ── 2) 在 B 系 sender → 巨噬 receiver 子集内的真实排名 ──
B_SOURCES_7 = ["B cells (other)","GC B cells","Naive B cells","Memory B cells",
               "Age-associated B cells","Proliferative GC B cells","Plasma cells"]
for name, df, bsrc in [("7 群", full7, B_SOURCES_7), ("13 亚型", full13, None)]:
    d = df if bsrc is None else df[df["source"].isin(bsrc)]
    if bsrc is None:   # 13 亚型：sender 含 "B" 或 "malignant"，先看有哪些
        print("\n13 亚型表 source 清单：", sorted(d["source"].unique()))
    mac = d[d["target"].str.contains("Macro|Mono", case=False, na=False)]
    for metric in ["magnitude_rank", "specificity_rank"]:
        sub = mac.sort_values(metric).reset_index(drop=True)
        sub["rank_in_subset"] = range(1, len(sub)+1)
        hit = sub[is_hh(sub)]
        if len(hit):
            h = hit.iloc[0]
            print(f"[{name}] B系→Mac 子集内, 按 {metric}: rank {int(h['rank_in_subset'])}/{len(sub)}, "
                  f"值 = {h[metric]:.4f}")
        else:
            print(f"[{name}] B系→Mac 子集({len(mac)}对)内无 HMGB1–HAVCR2 —— 检查 source/target 命名")

import pandas as pd

full = pd.read_csv(r"D:\bulk-download\GSE182434\cellcomm\liana_results_full.csv")
B_SOURCES = ["B cells (other)","GC B cells","Naive B cells","Memory B cells",
             "Age-associated B cells","Proliferative GC B cells","Plasma cells"]

# B→Mac 子集完整排名（回复信数字 86/1584 的来源，作 supplementary 证据）
mac = full[(full["source"].isin(B_SOURCES)) &
           (full["target"]=="Monocytes/Macrophages")].copy()
mac = mac.sort_values("magnitude_rank").reset_index(drop=True)
mac["rank_in_B_to_Mac_axis"] = range(1, len(mac)+1)
mac["passed_0.05"] = mac["magnitude_rank"] < 0.05
mac["retained_under_0.1"] = mac["magnitude_rank"] < 0.10
mac.to_csv(r"D:\bulk-download\GSE182434\cellcomm\TableS24_BtoMac_ranking_sensitivity.csv",
           index=False)
print(f"子集共 {len(mac)} 对；HMGB1–HAVCR2 名次：")
print(mac[mac["ligand_complex"]=="HMGB1"][["source","ligand_complex","receptor_complex",
      "magnitude_rank","rank_in_B_to_Mac_axis"]].to_string(index=False))

import pandas as pd
import numpy as np

df = pd.read_csv(r"D:\bulk-download\GSE182434\cellcomm\TableS24_BtoMac_ranking_sensitivity.csv")

# 保留支撑主张的列；丢弃全空的 specificity_rank 和含 -inf/NaN 的 lr_logfc
keep = ["source","target","ligand_complex","receptor_complex",
        "lr_means","expr_prod","scaled_weight",          # magnitude 证据列
        "magnitude_rank","rank_in_B_to_Mac_axis",
        "passed_0.05","retained_under_0.1"]
out = df[keep].copy()

# 阈值列转为 Yes/No（Excel 可读性）
out["passed_0.05"] = out["passed_0.05"].map({True:"Yes", False:"No"})
out["retained_under_0.1"] = out["retained_under_0.1"].map({True:"Yes", False:"No"})

out.to_csv(r"D:\bulk-download\GSE182434\cellcomm\TableS24_clean.csv", index=False)

# 顺手打几个供回复信/脚注用的统计量
print(f"总对数: {len(out)}")
print(f"passed_0.05: {(df['passed_0.05']).sum()} 对")
print(f"retained_under_0.1: {(df['retained_under_0.1']).sum()} 对")
print(out[out['ligand_complex']=='HMGB1'][['source','receptor_complex',
      'magnitude_rank','rank_in_B_to_Mac_axis']].to_string(index=False))
