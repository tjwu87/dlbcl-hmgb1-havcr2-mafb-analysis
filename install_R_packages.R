# =============================================================================
# R 依赖安装脚本
# -----------------------------------------------------------------------------
#   Rscript install_R_packages.R
#
# 用途：bulk 预后建模（LASSO-Cox / 列线图 / 校准曲线 / 时间依赖 ROC）、
#       免疫浸润（CIBERSORT / GSVA / ssGSEA）、富集分析、绘图。
# 建议 R >= 4.2。
# =============================================================================

cran <- c(
  # 基础
  "tidyverse", "data.table", "reshape2", "stringr", "readr",
  # 统计与建模
  "survival", "survminer", "glmnet", "rms", "Hmisc", "caret",
  "e1071", "randomForestSRC", "xgboost", "gbm", "mboost",
  "plsRglm", "ridge", "BART",
  # ROC / 校准
  "timeROC", "pROC", "survcomp", "ROCit",
  # 可视化
  "ggplot2", "ggpubr", "ggrepel", "ggdist", "ggbeeswarm", "cowplot",
  "patchwork", "pheatmap", "ComplexHeatmap", "circlize", "corrplot",
  "RColorBrewer", "paletteer", "scales", "svglite", "gridExtra",
  # 富集与基因集
  "GSVA", "GSEABase", "clusterProfiler", "enrichplot", "msigdbr",
  "org.Hs.eg.db", "GOSemSim", "fgsea",
  # 表达谱处理
  "limma", "Biobase", "SummarizedExperiment", "GEOquery", "impute",
  # 肿瘤基因组
  "TCGAbiolinks", "maftools",
  # 其他
  "Seurat", "harmony", "ConsensusClusterPlus", "forestploter",
  "ezcox", "tinyarray", "rstatix", "openxlsx", "snowfall"
)

bioc <- c(
  "Biobase", "SummarizedExperiment", "GEOquery", "limma", "GSVA",
  "GSEABase", "clusterProfiler", "enrichplot", "org.Hs.eg.db",
  "hgu133plus2.db", "illuminaHumanv4.db", "ComplexHeatmap",
  "TCGAbiolinks", "maftools", "impute", "fgsea"
)

installed <- rownames(installed.packages())

cat("== 安装 CRAN 包 ==\n")
for (p in setdiff(cran, installed)) {
  cat("  ->", p, "\n")
  try(install.packages(p, repos = "https://cloud.r-project.org"), silent = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("== 安装 BiocManager ==\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("== 安装 Bioconductor 包 ==\n")
installed <- rownames(installed.packages())
for (p in setdiff(bioc, installed)) {
  cat("  ->", p, "\n")
  try(BiocManager::install(p, ask = FALSE, update = FALSE), silent = TRUE)
}

# ── 需手动安装的包（不在 CRAN/Bioconductor）───────────────────────────────────
cat("
== 以下包需手动安装（许可或托管限制）==
  1. CIBERSORT
       从 https://ciberx.stanford.edu/ 注册下载 CIBERSORT.R 与
       LM22.txt，放到 data/reference/ 目录。
       注意：CIBERSORT 仅限学术研究使用。
  2. pRRophetic（药物敏感性预测）
       devtools::install_github('paulgeeleher/pRRophetic')
       依赖较大，且需要 GDSC2 训练数据（见 00_download/）。
  3. AUCell / pySCENIC 的 R 接口（如使用）
       BiocManager::install('AUCell')

== 验证 ==
  source('install_R_packages.R') 后运行 check_deps()
")

check_deps <- function() {
  need <- c("survival", "glmnet", "survminer", "timeROC", "rms",
            "ggplot2", "patchwork", "GSVA", "limma", "clusterProfiler")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) cat("缺失:", paste(miss, collapse = ", "), "\n")
  else cat("核心依赖齐全。\n")
}
