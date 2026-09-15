"""
一键复现入口（按论文 Results 小节组织；本归档含 section1–6）
================================================================================
    python run_all.py                 # 全部执行（section1 → section6）
    python run_all.py --stage 4       # 只跑 section4
    python run_all.py --list          # 只列出各 section 的脚本
    python run_all.py --check         # 只做输入数据预检，不执行
    python run_all.py --dry-run       # 只打印将要执行的命令

数据位置
--------------------------------------------------------------------------------
config/paths.py 会按以下顺序解析数据根，一般无需手动设置：
    1) 环境变量 DLBCL_DATA_ROOT
    2) <仓库>/00_共用/数据/     ← 归档自带的公用数据（推荐）
    3) <仓库>/data/
若要指向其它位置（例如完整的 D:\\bulk-download）：
    Windows : set DLBCL_DATA_ROOT=D:\\bulk-download
    Linux   : export DLBCL_DATA_ROOT=/path/to/data

依赖
--------------------------------------------------------------------------------
    conda env create -f 00_共用/环境/environment.yml
    R 依赖：Rscript 00_共用/环境/install_R_packages.R
各 section 需要的包见其 README.md 「环境状态」一节。
================================================================================
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))

from config.paths import DATA_ROOT, CELLORACLE, FIG_DIR, GSE182434  # noqa: E402

S1 = "section1_恶性B细胞与髓系通讯"
S2 = "section2_髓系轨迹与LA_TAM终末态"
S3 = "section3_HAVCR2-MAFB调控模块"
S4 = "section4_MAFB计算机扰动"
S5 = "section5_空间转录组验证"
S6 = "section6_MAFB预后评分模型"

# ── 与论文 Results 小节一一对应（本归档含 section1–6）──────────────────────────────────────────
STAGES = [
    ("1", "恶性 B 细胞与髓系通讯（Figure 1 / S1）", [
        ("py", f"{S1}/code/01_scRNA_GSE182434/01_01_qc_integration.py"),
        ("py", f"{S1}/code/01_scRNA_GSE182434/01_03_annotation_concordance.py"),
        ("py", f"{S1}/code/01_scRNA_GSE182434/01_04_supplementary_panels.py"),
        ("py", f"{S1}/code/02_cellcomm_LIANA/02_00_liana_run.py"),
        ("py", f"{S1}/code/02_cellcomm_LIANA/02_01_fig1de_aggregate_circle_heatmap.py"),
        ("py", f"{S1}/code/02_cellcomm_LIANA/02_02_liana_figures_13subtypes_volcano.py"),
        ("py", f"{S1}/code/02_cellcomm_LIANA/02_03_liana_supplementary.py"),
        ("py", f"{S1}/code/10_figures/10_00_umap_cnv_myeloid_panels.py"),
    ]),
    ("2", "髓系轨迹与 LA-TAM 终末态（Figure 2 / S2 / S8）", [
        ("py", f"{S2}/code/01_scRNA_GSE182434/01_02_myeloid_subtyping_and_paga_trajectory.py"),
        ("py", f"{S2}/code/01_scRNA_GSE182434/01_05_fig1b_marker_dotplot.py"),
        ("py", f"{S2}/code/01_scRNA_GSE182434/01_06_mono_mac_5subtypes.py"),
        ("py", f"{S2}/code/03_trajectory_PAGA/03_01_paga_trajectory.py"),
        ("py", f"{S2}/code/03_trajectory_PAGA/03_02_latam_score_comparison.py"),
        ("py", f"{S2}/code/10_figures/FigS2/10_10_FigS2_supplementary.py"),
        ("py", f"{S2}/code/10_figures/FigS8/10_11_FigS8_supplementary.py"),
    ]),
    ("3", "HAVCR2–MAFB 调控模块（Figure 3 / S3 / S9）", [
        ("py", f"{S3}/code/04_pySCENIC/04_01_pyscenic_downstream_stats_figures.py"),
        ("py", f"{S3}/code/10_figures/FigS8/10_12_figS3_S9_receptor_prioritization.py"),
    ]),
    ("4", "MAFB 计算机扰动 / CellOracle（Figure 4）", [
        ("setup_outdir", ""),
        ("py", f"{S4}/code/05_CellOracle/step00A_env.py"),
        ("py", f"{S4}/code/05_CellOracle/step00B_freeze_params.py"),
        ("py", f"{S4}/code/05_CellOracle/step01_preprocess.py"),
        ("py", f"{S4}/code/05_CellOracle/step02AB_oracle_build.py"),
        ("py", f"{S4}/code/05_CellOracle/step02C_get_links.py"),
        ("py", f"{S4}/code/05_CellOracle/step02DE_filter_audit.py"),
        ("py", f"{S4}/code/05_CellOracle/step02F_fit_grn.py"),
        ("py", f"{S4}/code/05_CellOracle/step03_perturbation.py"),
        ("py", f"{S4}/code/05_CellOracle/step03_null.py"),
        ("py", f"{S4}/code/05_CellOracle/step04A_fate_calc.py"),
        ("py", f"{S4}/code/05_CellOracle/step04B_fate_stats.py"),
        ("py", f"{S4}/code/05_CellOracle/step04_5_consistency_check.py"),
        ("py", f"{S4}/code/05_CellOracle/step05_scenic_targets.py"),
        ("py", f"{S4}/code/05_CellOracle/step06_program_score.py"),
        ("py", f"{S4}/code/05_CellOracle/step07_enrichment.py"),
        ("py", f"{S4}/code/05_CellOracle/step07C_ko_response_enrichment.py"),
        ("py", f"{S4}/code/05_CellOracle/step08AB_sensitivity.py"),
        ("py", f"{S4}/code/05_CellOracle/step08C_sensitivity_hvg4000.py"),
        ("py", f"{S4}/code/05_CellOracle/step08D_sensitivity_n30.py"),
        ("py", f"{S4}/code/05_CellOracle/step08E_k_sensitivity.py"),
        ("py", f"{S4}/code/05_CellOracle/step08_summary.py"),
        ("py", f"{S4}/code/05_CellOracle/step09_export_data.py"),
        ("py", f"{S4}/code/05_CellOracle/step09_fig23_deg_ora.py"),
        ("py", f"{S4}/code/05_CellOracle/step04_scoring_figures.py"),
        ("py", f"{S4}/code/05_CellOracle/05_90_reproduce_figures_local.py"),
        ("py", f"{S4}/code/05_CellOracle/05_91_grn_edge_sensitivity.py"),
        ("py", f"{S4}/code/05_CellOracle/05_91_reproduce_bigfont.py"),
    ]),
    ("5", "空间转录组验证 GeoMx（Figure 5 / S5）", [
        ("py", f"{S5}/code/06_spatial_GSE232853/06_01_spatial_analysis.py"),
    ]),
    ("6", "MAFB 预后评分模型（Figure 6 / S6 + 补充表）", [
        ("R", f"{S6}/code/07_bulk_prognosis/07_02_lasso_cv_curves.R"),
        ("R", f"{S6}/code/07_bulk_prognosis/07_01_prognosis_main.R"),
        ("R", f"{S6}/code/07_bulk_prognosis/07_03_immune_deconvolution.R"),
        ("R", f"{S6}/code/07_bulk_prognosis/07_04_supplementary_tables.R"),
        ("R", f"{S6}/code/07_bulk_prognosis/07_05_cibersort.R"),
        ("py", f"{S6}/code/10_figures/supp_tables/tableS23_ARI.py"),
        ("py", f"{S6}/code/10_figures/supp_tables/tableS5C_GRN_5000.py"),
        ("py", f"{S6}/code/10_figures/supp_tables/S24_BtoMac_ranking.py"),
        ("py", f"{S6}/code/10_figures/supp_tables/mono_5gene_vs_16gene_score.py"),
    ]),
]

SKIP_MISSING = True   # 数据/依赖未就绪时跳过而非中断


def setup_outdir() -> None:
    """把 CellOracle 输出目录写入 tmp_manifest/outdir.txt。

    step01–step08 全部通过 outdir.txt 解析 OUTDIR，且期望其中是
    分层结构（00_audit/、01_preprocessing/、02_oracle_grn/…）。
    必须指向 celloracle_paper_final_*，不能指向扁平的 celloracle0331。
    """
    target = REPO / S4 / "code" / "05_CellOracle" / "tmp_manifest"
    target.mkdir(parents=True, exist_ok=True)
    # 优先使用带分层结构的 celloracle_paper_final_*，否则回退到 CELLORACLE。
    _cands = sorted(GSE182434.glob("celloracle_paper_final_*")) \
        if GSE182434.is_dir() else []
    root = _cands[-1] if _cands else CELLORACLE
    (target / "outdir.txt").write_text(str(root), encoding="utf-8")
    print(f"[setup] outdir.txt -> {root}")


def run(kind: str, script: str, dry: bool) -> bool:
    if kind == "setup_outdir":
        print("[setup] 写入 CellOracle 输出根目录")
        if not dry:
            setup_outdir()
        return True

    path = REPO / script
    if not path.exists():
        print(f"[skip] 缺失 {script}")
        return SKIP_MISSING

    cmd = ["python", str(path)] if kind == "py" else ["Rscript", str(path)]
    print(f"\n$ {' '.join(cmd)}")
    if dry:
        return True

    env = os.environ.copy()
    env["DLBCL_REPO_ROOT"] = str(REPO)
    env.setdefault("DLBCL_DATA_ROOT", str(DATA_ROOT))
    ret = subprocess.run(cmd, cwd=str(REPO), env=env).returncode
    if ret != 0:
        print(f"[fail] {script} (exit {ret})")
        return False
    return True


def show_list() -> None:
    for code, desc, steps in STAGES:
        print(f"\nsection{code} — {desc}")
        for kind, script in steps:
            mark = "OK " if (kind == "setup_outdir" or (REPO / script).exists()) else "缺失"
            print(f"   [{mark}] {script or '(setup)'}")


def preflight() -> None:
    """输入数据预检：脚本引用的路径逐个 translate() 后判存在。"""
    import re

    legacy = ("d:/bulk-download", "d:\\bulk-download", "/mnt/results",
              "d:/scrna", "d:\\scrna", "d:/badidunetdiskdownload")
    pat = re.compile(r"""['"]([^'"\n]{4,300})['"]""")
    chk = re.compile(r"\.(h5ad|h5|pkl|npz|csv|tsv|txt|gz|xlsx|rds|rdata|npy|obj|mtx)$", re.I)
    from config.paths import translate

    total = hit = 0
    problems = []
    for code, desc, steps in STAGES:
        for kind, script in steps:
            if not script:
                continue
            p = REPO / script
            if not p.exists():
                problems.append(f"[缺失脚本] {script}")
                continue
            txt = p.read_text(encoding="utf-8", errors="replace")
            for m in pat.finditer(txt):
                s = m.group(1)
                if "{" in s or s.startswith(("http", "ftp")):
                    continue
                sl = s.replace("\\", "/").lower()
                if not (any(sl.startswith(r) for r in legacy) or (sl.count("/") >= 2 and chk.search(sl))):
                    continue
                total += 1
                try:
                    t = translate(s)
                except Exception:
                    continue
                if t.exists():
                    hit += 1
                else:
                    problems.append(f"[输入缺失] {script}\n             {s}\n          -> {t}")

    print(f"\n路径引用检查：{total} 处，命中 {hit}，问题 {len(problems)}")
    for x in problems:
        print("  " + x)
    print("\n（运行期产物如 /tmp/*.pkl、OUTDIR/*.csv 不在检查范围内）")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", help="只运行指定 section，如 4")
    ap.add_argument("--from", dest="from_stage", help="从指定 section 开始")
    ap.add_argument("--host", dest="from_host", help="从指定 section 开始（含）")
    ap.add_argument("--to", dest="to_stage", help="运行到指定 section 为止")
    ap.add_argument("--list", action="store_true", help="只列出脚本")
    ap.add_argument("--check", action="store_true", help="只做数据预检")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.list:
        show_list()
        return

    print("=" * 70)
    print("DLBCL HMGB1–HAVCR2–MAFB 复现流程")
    print(f"  仓库根   : {REPO}")
    print(f"  数据根   : {DATA_ROOT}   {'存在' if DATA_ROOT.is_dir() else '★不存在'}")
    print(f"  结果输出 : {FIG_DIR.parent}")
    print("=" * 70)

    if args.check:
        return preflight()

    start = args.from_stage or args.from_host
    for code, desc, steps in STAGES:
        if args.stage and code != args.stage:
            continue
        if start and code < start:
            continue
        if args.to_stage and code > args.to_stage:
            continue
        print(f"\n{'=' * 70}\nsection{code} — {desc}\n{'=' * 70}")
        for kind, script in steps:
            if not run(kind, script, args.dry_run):
                print(f"\n[中止] section{code} 在 {script} 处失败。")
                sys.exit(1)

    print("\n完成。图片输出于:", FIG_DIR)


if __name__ == "__main__":
    main()
