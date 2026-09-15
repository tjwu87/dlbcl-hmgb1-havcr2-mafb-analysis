# 髓系轨迹与LA_TAM终末态

> **论文对应**：Results 2 — 髓系细胞轨迹与 LA-TAM 终末态  
> **图件**：Figure 2（a–f）、Figure S2、Figure S8 部分

## 分析内容

髓系再分型 → PAGA 轨迹 → CellRank 命运驱动 → LA-TAM 打分与配对比较

## 目录

```
section2_髓系轨迹与LA_TAM终末态/
├── code/            分析脚本（见下）
├── data/            本部分专有的中间产物（可比对/复用）
├── figures/         论文图件与面板
└── README.md
```

## 脚本清单

| 脚本 | 说明 |
| --- | --- |
| `code/01_scRNA_GSE182434/01_02_myeloid_subtyping_and_paga_trajectory.py` | |
| `code/01_scRNA_GSE182434/01_05_fig1b_marker_dotplot.py` | |
| `code/01_scRNA_GSE182434/01_06_mono_mac_5subtypes.py` | |
| `code/03_trajectory_PAGA/03_01_paga_trajectory.py` | |
| `code/03_trajectory_PAGA/03_02_latam_score_comparison.py` | |
| `code/10_figures/FigS2/10_10_FigS2_supplementary.py` | |
| `code/10_figures/FigS8/10_11_FigS8_supplementary.py` | |

## 输入数据（公用，位于 `data/`）

- `GSE182434/adata_mac_subtyped.h5ad`
- `GSE182434/trajectory/adata_mac_annotated_patched.h5ad`
- `GSE182434/trajectory/cellrank_LA_TAM_fate_drivers_receptors.csv`
- `GSE182434/trajectory/receptor_pseudotime_cell_level.csv`

## 复现命令

```bash
# 数据根无需设置：config/paths.py 默认解析到 data/
cd DLBCL_HMGB1_HAVCR2_MAFB
# 全流程
python run_all.py
# 只跑本部分
python run_all.py --stage 2
```

## 环境状态

- 需要：scanpy, anndata, cellrank, palantir
- 本机 `scRNA` 环境已具备：scanpy, anndata
- 尚缺：cellrank, palantir

## 产出图件

- `figures/论文成图/Fig2.png`
- `figures/论文成图/FigS2.png`
- `figures/论文成图/Figure S2.png`
- `figures/论文成图/Figure2.pdf`
- `figures/论文成图/Figure2.png`
- `figures/重排版/Fig2.png`
- `figures/面板原件/Figure 2.1.png`
- `figures/面板原件/GSE182434_trajectory_monomac_gene_dynamics修.png`
- `figures/面板原件/S2.png`
- `figures/面板原件/enrichment_barplot_IFN_TAM_vs_LA_TAM (2).png`
- `figures/面板原件/enrichment_barplot_IFN_TAM_vs_LA_TAM.png`
- `figures/面板原件/mac_subtype_dotplot.png`
- `figures/面板原件/mac_subtype_dotplot2.png`
- `figures/面板原件/monomac_gene_dynamics.png`
- `figures/面板原件/monomac_paga_trajectory.png`
- `figures/面板原件/monomac_paga_trajectory_rawpt.png`
- `figures/面板原件/monomac_trajectory_directed.png`
- `figures/面板原件/monomac_trajectory_directed_rawpt.png`
- `figures/面板原件/s8.png`
- `figures/面板原件/volcano_IFN_TAM_vs_LA_TAM.png`
- `figures/面板原件/volcano_IFN_TAM_vs_rest.png`
