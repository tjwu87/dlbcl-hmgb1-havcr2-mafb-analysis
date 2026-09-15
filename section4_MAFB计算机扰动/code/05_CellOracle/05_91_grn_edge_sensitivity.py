# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from config.paths import translate, DATA_ROOT, FIG_DIR, TAB_DIR
# <<< 自动注入结束 >>>

import pandas as pd

base = translate(r"D:\bulk-download\GSE182434\celloracle_paper_final_20260331_085223\02_oracle_grn")
clusters = ["IFN_TAM", "LA_TAM", "Mono"]
rows, rank_rows = [], []

for c in clusters:
    A = pd.read_csv(base + rf"\filter_A_t2000\filtered_A_{c}.csv")
    B = pd.read_csv(base + rf"\filter_B_t5000\filtered_B_{c}.csv")

    ea = set(map(tuple, A[["source","target"]].values))
    eb = set(map(tuple, B[["source","target"]].values))

    mA, mB = A[A["source"]=="MAFB"], B[B["source"]=="MAFB"]
    ta, tb = set(mA["target"]), set(mB["target"])
    topA = set(mA.nlargest(20,"coef_abs")["target"])
    topB = set(mB.nlargest(20,"coef_abs")["target"])

    # MAFB 最强边在全网络中的排位
    b_rank = int((B["coef_abs"] > mB["coef_abs"].max()).sum()) + 1

    rows.append({"cluster": c,
        "edges_t2000": len(ea), "edges_t5000": len(eb),
        "jaccard_t2000_t5000": round(len(ea&eb)/len(ea|eb), 3),
        "MAFB_edges_t2000": len(ta), "MAFB_edges_t5000": len(tb),
        "MAFB_top20_overlap": f"{len(topA&topB)}/20",
        "MAFB_strongest_edge": mB.sort_values("coef_abs", ascending=False).iloc[0]["target"],
        "MAFB_top_edge_rank_in_network": b_rank,
        "t5000_coef_cutoff": round(B["coef_abs"].min(), 4),
        "t5000_max_p": B["p"].max()})

pd.DataFrame(rows).to_csv(base + r"\TableS5f_GRN_edge_cap_concordance.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))

import pandas as pd
base = translate(r"D:\bulk-download\GSE182434\celloracle_paper_final_20260331_085223\02_oracle_grn")
for c in ["IFN_TAM","LA_TAM","Mono"]:
    A = pd.read_csv(base + rf"\filter_A_t2000\filtered_A_{c}.csv")
    B = pd.read_csv(base + rf"\filter_B_t5000\filtered_B_{c}.csv")
    ta = set(A[A["source"]=="MAFB"]["target"])
    topB = set(B[B["source"]=="MAFB"].nlargest(20,"coef_abs")["target"])
    n = min(20, len(ta))
    print(f"{c}: t2000 MAFB targets {len(ta)}, retained in t5000 top-20: {len(ta&topB)}/{n}")
