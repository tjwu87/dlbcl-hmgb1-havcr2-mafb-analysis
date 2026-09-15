# HAVCR2-MAFB调控模块

> **论文对应**：Results 3 — HAVCR2–MAFB 调控模块  
> **图件**：Figure 3（a–g）、Figure S3、Figure S9

## 分析内容

pySCENIC 调控子活性（AUCell）→ RSS 受体优先级排序 → 伪时序相关 → MAFB 靶基因富集

## 目录

```
section3_HAVCR2-MAFB调控模块/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/04_pySCENIC/04_01_pyscenic_downstream_stats_figures.py` | |
| `code/04_pySCENIC/04_00a_setup_scenic_db.sh` | 上游本体（找回件）：下载 pySCENIC 参考数据库（TF 列表 / hg38 rankings / motif 注释），文件头含官方 CLI 三步命令 |
| `code/04_pySCENIC/04_00b_build_loom.py` | 上游本体（找回件）：构建 loom 输入（327 细胞 × 9118 基因） |
| `code/04_pySCENIC/04_00c_pyscenic_run_and_export.py` | 上游本体（找回件）：GRNBoost2 → RcisTarget(NES≥3.0) → AUCell，并导出 `data_*.csv` —— 见 `code/04_pySCENIC/README_recovered.md` |
| `code/10_figures/FigS8/10_12_figS3_S9_receptor_prioritization.py` | |

## 输入数据（公用，位于 `data/`）

- `GSE182434/scenic/data_AUC_matrix.csv`
- `GSE182434/scenic/data_RSS_scores.csv`
- `GSE182434/scenic/data_RcisTarget_regulon_targets.csv`
- `GSE182434/scenic/data_core_TF_AUC_and_receptor_expression.csv`
- `GSE182434/scenic/MAFB_gseapy_*.csv`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 3
```

## 环境状态

- 需要：pyscenic, ctxcore, omnipath, decoupler
- 本机 `scRNA` 环境已具备：（无）
- 尚缺：pyscenic, ctxcore, omnipath, decoupler

## 产出图件

- `figures/论文成图/Fig3.png`
- `figures/论文成图/FigS3.png`
- `figures/论文成图/Figure S3.png`
- `figures/论文成图/Figure3.pdf`
- `figures/论文成图/Figure3.png`
- `figures/重排版/Fig3a.png`
- `figures/重排版/Fig3b.png`
- `figures/面板原件/3A.png`
- `figures/面板原件/Fig1_regulon_activity_heatmap.png`
- `figures/面板原件/Fig1_regulon_activity_heatmap_paper.png`
- `figures/面板原件/Fig2_RSS_scatter.png`
- `figures/面板原件/Fig3(1).png`
- `figures/面板原件/Fig3(2).png`
- `figures/面板原件/FigS3_receptor_prioritization.png`
- `figures/面板原件/Figure 3.1.png`
- `figures/面板原件/Figure3.3.png`
- `figures/面板原件/Figure3.4.png`
- `figures/面板原件/MAFB_enrichment_dotplot.png`
- `figures/面板原件/MAFB_gseapy_enrichment_dotplot.png`
- `figures/面板原件/S3.png`
- `figures/面板原件/lineplot_LA_TAM_TF_receptor_2x4.png`
- `figures/面板原件/receptor_pseudotime_heatmap.png`
- `figures/面板原件/receptor_pseudotime_heatmap_path.png`
- `figures/面板原件/receptor_pseudotime_lineplots.png`
- `figures/面板原件/receptor_pseudotime_lineplots_path.png`
- `figures/面板原件/receptor_pseudotime_lineplots_path_dual.png`
- `figures/面板原件/sankey_7layer_TF_targets_v3.png`
- `figures/面板原件/scenic_tf_dynamics.png`
