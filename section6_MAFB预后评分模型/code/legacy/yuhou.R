# =============================================================================
# 【归档收录说明】本文件取自作者原始脚本 Tcell/yuhou.R（3441 行）。
#
# 为什么收录：section6 的 07_01~07_05 只**读**下列中间表，仓库里没有产出它们的
# 代码；而本文件正是这些表的**唯一产出者**：
#   - LASSO_prognostic_formula.csv  (L289)
#   - GSE10846_risk_scores.csv  (L304)
#   - Cox_univariable_multivariable_results.csv  (L957)
#   - Calibration_data_1_3_5yr.csv  (L1571)
#   - Immune_infiltration_scores.csv / Immune_RiskScore_correlations.csv  (L1731/L1734)
#   - TCGA_DLBC_risk_scores.csv  (L3414)
#   - GSE87371 / GSE11318 / GSE181063_risk_scores.csv  (L3421/L3426/L3431)
#   - All_cohorts_risk_scores.csv  (L3356)
#   - MultiCohort_KM_statistics.csv  (L3404)
#
# 使用前须知：
#  1) 已注入 config/paths.R，并把原文的服务器绝对路径（/mnt/results、/workspace）
#     用 translate_path() 包裹，输出会落到当前数据根下；映射表见
#     code/legacy/README.md。
#  2) 其中若干段是当年的**探索过程**（多次 LASSO 参数尝试、CV 曲线重绘、
#     队列方向性排查、被取代的实现），不是必需步骤 —— 行号清单见
#     code/legacy/README.md，可整段跳过。
#  3) 本文件是**归档件，不进入 run_all.py 流程**：它依赖上游产物
#     （GSE10846 series matrix、LASSO_prognostic_formula.csv 等）。
#  4) 下文中的行号注释均为**原文行号**，收录时仅在文件头部追加了本段说明。
#  5) ★ 本文件是 **R + 内嵌 Python 的混合工作文件**（含 import requests、pd.read_csv
#     等 Python 片段），整文件无法直接 source() 运行；当年应是分块手动执行。
#     收录目的是保留产出逻辑与参数，不是提供一键入口。
# =============================================================================

source(file.path("config", "paths.R"))


library(GEOquery)
library(dplyr)
library(stringr)

# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10846
cat("Parsing GSE10846 series matrix...\n")
gse <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)

cat("Class:", class(gse), "\n")
cat("Dimensions:", nrow(gse), "probes x", ncol(gse), "samples\n")

# Extract phenotype data
pdata <- pData(gse)
cat("\nColumn names in pData:\n")
print(colnames(pdata))

# Inspect characteristics columns
char_cols <- grep("characteristics_ch1", colnames(pdata), value=TRUE)
for(col in char_cols){
  vals <- unique(pdata[[col]])
  cat("\n---", col, "---\n")
  print(head(vals, 5))
}

library(dplyr)
library(stringr)

# Helper to strip prefix
strip <- function(x, prefix) str_trim(str_replace(x, fixed(prefix), ""))

# Build clinical dataframe
clin <- data.frame(
  sample_id   = pdata$geo_accession,
  gender      = strip(pdata$`characteristics_ch1`,   "Gender: "),
  age         = as.numeric(strip(pdata$`characteristics_ch1.1`, "Age: ")),
  status      = strip(pdata$`characteristics_ch1.7`,  "Clinical info: Follow up status: "),
  follow_up_years = as.numeric(strip(pdata$`characteristics_ch1.8`, "Clinical info: Follow up years: ")),
  chemo       = strip(pdata$`characteristics_ch1.9`,  "Clinical info: Chemotherapy: "),
  ecog        = as.numeric(strip(pdata$`characteristics_ch1.10`, "Clinical info: ECOG performance status: ")),
  stage       = as.numeric(strip(pdata$`characteristics_ch1.11`, "Clinical info: Stage: ")),
  ldh_ratio   = as.numeric(strip(pdata$`characteristics_ch1.12`, "Clinical info: LDH ratio: ")),
  extranodal  = as.numeric(strip(pdata$`characteristics_ch1.13`, "Clinical info: Number of extranodal sites: ")),
  subtype     = strip(pdata$`characteristics_ch1.6`,  "Clinical info: Final microarray diagnosis: "),
  stringsAsFactors = FALSE
)

# OS event: DEAD=1, ALIVE=0
clin$OS_event <- ifelse(clin$status == "DEAD", 1, ifelse(clin$status == "ALIVE", 0, NA))
clin$OS_years <- clin$follow_up_years

# IPI score (0-5): age>60, stage III/IV, LDH>1, ECOG>=2, extranodal>1
clin$ipi_age      <- ifelse(!is.na(clin$age)       & clin$age > 60,       1, 0)
clin$ipi_stage    <- ifelse(!is.na(clin$stage)     & clin$stage >= 3,     1, 0)
clin$ipi_ldh      <- ifelse(!is.na(clin$ldh_ratio) & clin$ldh_ratio > 1,  1, 0)
clin$ipi_ecog     <- ifelse(!is.na(clin$ecog)      & clin$ecog >= 2,      1, 0)
clin$ipi_extranod <- ifelse(!is.na(clin$extranodal)& clin$extranodal > 1, 1, 0)
clin$IPI <- clin$ipi_age + clin$ipi_stage + clin$ipi_ldh + clin$ipi_ecog + clin$ipi_extranod

# Keep only R-CHOP treated (more relevant, n=233)
clin_rchop <- clin %>% filter(grepl("R-CHOP", chemo))
cat("R-CHOP patients:", nrow(clin_rchop), "\n")
cat("With complete OS data:", sum(!is.na(clin_rchop$OS_event) & !is.na(clin_rchop$OS_years)), "\n")
cat("Events (deaths):", sum(clin_rchop$OS_event==1, na.rm=TRUE), "\n")
cat("IPI distribution:\n"); print(table(clin_rchop$IPI, useNA="ifany"))
cat("Stage distribution:\n"); print(table(clin_rchop$stage, useNA="ifany"))
cat("Age summary:\n"); print(summary(clin_rchop$age))


library(Biobase)

# Define 25-gene set
gene_set <- c("HMGB1","HAVCR2",
              "CEP192","CHCHD10","NAGK","PPRC1","BASP1","MTR","TIMM17A","AAK1",
              "MGAT1","USP12","FAM20A","APOL3","GSDMD","ELMO1","PICALM","RASSF4",
              "IFIT3","NAE1","NAP1L4","PPP2R5A","MRGBP","PPP1R18","CTSL")
cat("Gene set (n=25):", paste(gene_set, collapse=", "), "\n\n")

# Get expression matrix
expr_mat <- exprs(gse)
cat("Expression matrix:", nrow(expr_mat), "probes x", ncol(expr_mat), "samples\n")

# Download GPL570 annotation to map probe -> gene symbol
cat("Downloading GPL570 annotation...\n")
gpl <- getGEO("GPL570", destdir=translate_path("/workspace/tmp_geo"))
feat <- Table(gpl)
cat("GPL570 annotation columns:", paste(colnames(feat)[1:10], collapse=", "), "\n")
cat("Rows:", nrow(feat), "\n")

library(dplyr)
library(stringr)

# Get Gene Symbol column
cat("GPL570 columns:\n")
print(colnames(feat))

# Gene Symbol column
sym_col <- grep("Gene Symbol", colnames(feat), value=TRUE, ignore.case=TRUE)
cat("Symbol column:", sym_col, "\n")

# Build probe->symbol map
probe_map <- feat %>%
  select(ID, Gene_Symbol = all_of(sym_col)) %>%
  filter(Gene_Symbol != "" & !is.na(Gene_Symbol)) %>%
  # Some probes map to multiple genes (/// separated) - take first
  mutate(Gene_Symbol = str_trim(str_split_fixed(Gene_Symbol, "///", 2)[,1]))

cat("Probes with gene symbol:", nrow(probe_map), "\n")

# Filter to our 25 genes
probe_25 <- probe_map %>% filter(Gene_Symbol %in% gene_set)
cat("Probes matching 25-gene set:", nrow(probe_25), "\n")
cat("Genes found:\n")
print(sort(unique(probe_25$Gene_Symbol)))
cat("Genes NOT found:\n")
print(setdiff(gene_set, unique(probe_25$Gene_Symbol)))

library(dplyr)

# Subset expression to probes of interest
expr_25 <- expr_mat[probe_25$ID, , drop=FALSE]

# For each gene, keep probe with highest mean expression across all samples
probe_25$mean_expr <- rowMeans(expr_mat[probe_25$ID, , drop=FALSE])

best_probes <- probe_25 %>%
  group_by(Gene_Symbol) %>%
  slice_max(mean_expr, n=1, with_ties=FALSE) %>%
  ungroup()

cat("Best probe per gene (n=25):\n")
print(best_probes %>% select(Gene_Symbol, ID, mean_expr) %>% arrange(Gene_Symbol))

# Build gene x sample matrix
expr_gene <- expr_mat[best_probes$ID, , drop=FALSE]
rownames(expr_gene) <- best_probes$Gene_Symbol

# Subset to R-CHOP samples only
rchop_ids <- clin_rchop$sample_id
expr_rchop <- expr_gene[, rchop_ids, drop=FALSE]
cat("\nExpression matrix (R-CHOP): ", nrow(expr_rchop), "genes x", ncol(expr_rchop), "samples\n")

# Transpose: samples x genes (for LASSO)
expr_t <- t(expr_rchop)  # 233 x 25
cat("Transposed: ", nrow(expr_t), "samples x", ncol(expr_t), "genes\n")

# Check alignment
cat("Sample ID match:", all(rownames(expr_t) == clin_rchop$sample_id), "\n")
cat("Any NA in expression:", any(is.na(expr_t)), "\n")
cat("Expression range:", round(range(expr_t), 2), "\n")

library(glmnet)
library(survival)

set.seed(42)

# Prepare survival object
surv_obj <- Surv(clin_rchop$OS_years, clin_rchop$OS_event)

# Scale expression (important for LASSO coefficient interpretation)
X <- scale(expr_t)

# 10-fold CV LASSO Cox (alpha=1)
cv_fit <- cv.glmnet(X, surv_obj, family="cox", alpha=1, nfolds=10,
                    type.measure="deviance", standardize=FALSE)

cat("lambda.min:", cv_fit$lambda.min, "\n")
cat("lambda.1se:", cv_fit$lambda.1se, "\n")
cat("log(lambda.min):", log(cv_fit$lambda.min), "\n")
cat("log(lambda.1se):", log(cv_fit$lambda.1se), "\n")

# Coefficients at lambda.min
coef_min <- coef(cv_fit, s="lambda.min")
nonzero <- coef_min[coef_min[,1] != 0, , drop=FALSE]
cat("\nNon-zero genes at lambda.min (n=", nrow(nonzero), "):\n")
print(round(nonzero, 5))

# Coefficients at lambda.1se
coef_1se <- coef(cv_fit, s="lambda.1se")
nonzero_1se <- coef_1se[coef_1se[,1] != 0, , drop=FALSE]
cat("\nNon-zero genes at lambda.1se (n=", nrow(nonzero_1se), "):\n")
print(round(nonzero_1se, 5))

library(glmnet)
library(survival)

# Check OS_years distribution
cat("OS_years summary:\n"); print(summary(clin_rchop$OS_years))
cat("Zero or negative times:", sum(clin_rchop$OS_years <= 0, na.rm=TRUE), "\n")
cat("NA times:", sum(is.na(clin_rchop$OS_years)), "\n")

# Filter: keep only positive OS_years and non-NA event
keep <- !is.na(clin_rchop$OS_years) & clin_rchop$OS_years > 0 & !is.na(clin_rchop$OS_event)
cat("Samples kept:", sum(keep), "/ removed:", sum(!keep), "\n")

clin_lasso <- clin_rchop[keep, ]
X_lasso    <- scale(expr_t[keep, ])

set.seed(42)
surv_obj <- Surv(clin_lasso$OS_years, clin_lasso$OS_event)
cat("Events:", sum(clin_lasso$OS_event), "/ Total:", nrow(clin_lasso), "\n")

# 10-fold CV LASSO Cox
cv_fit <- cv.glmnet(X_lasso, surv_obj, family="cox", alpha=1, nfolds=10,
                    type.measure="deviance", standardize=FALSE)

cat("\nlambda.min:", cv_fit$lambda.min, "  log:", round(log(cv_fit$lambda.min),3), "\n")
cat("lambda.1se:", cv_fit$lambda.1se,  "  log:", round(log(cv_fit$lambda.1se),3), "\n")

# Coefficients at lambda.min
coef_min <- coef(cv_fit, s="lambda.min")
nonzero  <- coef_min[coef_min[,1] != 0, , drop=FALSE]
cat("\nNon-zero genes at lambda.min (n=", nrow(nonzero), "):\n")
print(round(nonzero, 5))

# Also check lambda.1se
coef_1se    <- coef(cv_fit, s="lambda.1se")
nonzero_1se <- coef_1se[coef_1se[,1] != 0, , drop=FALSE]
cat("\nNon-zero genes at lambda.1se (n=", nrow(nonzero_1se), "):\n")
print(round(nonzero_1se, 5))

library(glmnet)
library(survival)

# Use ALL 420 samples (CHOP + R-CHOP) for more power
clin_all <- clin %>% filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event))
cat("Full cohort: n=", nrow(clin_all), " events=", sum(clin_all$OS_event), "\n")

# Match expression
expr_all_ids <- intersect(clin_all$sample_id, colnames(expr_gene))
clin_all     <- clin_all[clin_all$sample_id %in% expr_all_ids, ]
expr_all_t   <- t(expr_gene[, clin_all$sample_id, drop=FALSE])
cat("After matching: n=", nrow(clin_all), " events=", sum(clin_all$OS_event), "\n")
cat("EPV (events/25 genes):", round(sum(clin_all$OS_event)/25, 1), "\n")

X_all <- scale(expr_all_t)
surv_all <- Surv(clin_all$OS_years, clin_all$OS_event)

set.seed(42)

# Try elastic net (alpha=0.5) for better gene retention
cv_en <- cv.glmnet(X_all, surv_all, family="cox", alpha=0.5, nfolds=10,
                   type.measure="deviance", standardize=FALSE)

coef_en_min <- coef(cv_en, s="lambda.min")
nz_en <- coef_en_min[coef_en_min[,1] != 0, , drop=FALSE]
cat("\nElastic net (alpha=0.5) lambda.min — non-zero genes (n=", nrow(nz_en), "):\n")
print(round(nz_en, 5))

# Also try pure LASSO on full cohort
cv_lasso_all <- cv.glmnet(X_all, surv_all, family="cox", alpha=1, nfolds=10,
                           type.measure="deviance", standardize=FALSE)
coef_lasso_all <- coef(cv_lasso_all, s="lambda.min")
nz_lasso <- coef_lasso_all[coef_lasso_all[,1] != 0, , drop=FALSE]
cat("\nLASSO (alpha=1) full cohort lambda.min — non-zero genes (n=", nrow(nz_lasso), "):\n")
print(round(nz_lasso, 5))

# Also check at a slightly relaxed lambda (lambda at nzero=5)
lam_seq  <- cv_lasso_all$lambda
nzero_seq <- cv_lasso_all$nzero
cat("\nLambda path — nzero genes:\n")
print(data.frame(log_lambda=round(log(lam_seq),3), nzero=nzero_seq) %>% 
      filter(nzero >= 3 & nzero <= 10) %>% head(15))

library(ggplot2)
library(dplyr)
library(stringr)

# ── 1. Final model: LASSO alpha=1, full cohort, lambda.min ──────────────────
lasso_genes <- rownames(nz_lasso)
lasso_coefs <- as.numeric(nz_lasso[,1])
names(lasso_coefs) <- lasso_genes

cat("=== PROGNOSTIC FORMULA ===\n")
cat("RiskScore =")
for(i in seq_along(lasso_coefs)){
  sign_str <- ifelse(lasso_coefs[i] >= 0, " + ", " - ")
  cat(sign_str, round(abs(lasso_coefs[i]),5), "*", lasso_genes[i])
}
cat("\n\n")

# Save formula
formula_df <- data.frame(
  Gene        = lasso_genes,
  Coefficient = lasso_coefs,
  Direction   = ifelse(lasso_coefs > 0, "Risk", "Protective"),
  stringsAsFactors = FALSE
)
write.csv(formula_df, translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"), row.names=FALSE)

# ── 2. Compute risk scores for all 412 patients ──────────────────────────────
# Use UNSCALED expression (raw probe values) with scaled coefficients
# Standard: use the scaled X that was used to fit the model
risk_scores <- as.numeric(X_all[, lasso_genes, drop=FALSE] %*% lasso_coefs)
clin_all$RiskScore <- risk_scores
clin_all$RiskGroup <- ifelse(risk_scores >= median(risk_scores), "High", "Low")

cat("Risk score summary:\n"); print(summary(risk_scores))
cat("High risk:", sum(clin_all$RiskGroup=="High"), "  Low risk:", sum(clin_all$RiskGroup=="Low"), "\n")

# Save risk scores
write.csv(clin_all %>% select(sample_id, RiskScore, RiskGroup, OS_years, OS_event,
                               age, gender, stage, IPI, chemo),
          translate_path("/mnt/results/DLBCL_prognosis/tables/GSE10846_risk_scores.csv"), row.names=FALSE)
cat("Risk scores saved.\n")

library(ggplot2)
library(dplyr)

# ── Extract CV data ──────────────────────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_lasso_all$lambda),
  cvm        = cv_lasso_all$cvm,
  cvsd       = cv_lasso_all$cvsd,
  cvup       = cv_lasso_all$cvup,
  cvlo       = cv_lasso_all$cvlo,
  nzero      = cv_lasso_all$nzero
)

lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)

# ── Plot ─────────────────────────────────────────────────────────────────────
p_cv <- ggplot(cv_df, aes(x=log_lambda, y=cvm)) +
  # Error ribbon
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), fill="#BFCFDF", alpha=0.5) +
  # CV mean line
  geom_line(color="#2166AC", linewidth=0.8) +
  # Points
  geom_point(color="#2166AC", size=1.2, shape=16) +
  # lambda.min line
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.8) +
  # lambda.1se line
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.8) +
  # Annotations
  annotate("text", x=lam_min_log, y=max(cv_df$cvup)*0.98,
           label=paste0("λ.min\n(", cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_min_log))], " genes)"),
           color="#D73027", size=3.2, hjust=1.1, fontface="bold") +
  annotate("text", x=lam_1se_log, y=max(cv_df$cvup)*0.98,
           label=paste0("λ.1se\n(", cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_1se_log))], " genes)"),
           color="#4DAC26", size=3.2, hjust=-0.1, fontface="bold") +
  # Top axis: number of non-zero genes
  scale_x_continuous(
    name = "log(λ)",
    sec.axis = sec_axis(~., name="Number of Non-zero Genes",
                        breaks = cv_df$log_lambda[seq(1, nrow(cv_df), length.out=12)],
                        labels = cv_df$nzero[seq(1, nrow(cv_df), length.out=12)])
  ) +
  labs(
    title    = "LASSO Cox Regression — 10-fold Cross-Validation",
    subtitle = "GSE10846 DLBCL (n=412, events=163)",
    y        = "Partial Likelihood Deviance (±1 SE)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title    = element_text(face="bold", size=13),
    plot.subtitle = element_text(size=10, color="grey40"),
    axis.title    = element_text(size=11),
    axis.text     = element_text(size=10),
    axis.title.x.top = element_text(size=10, color="grey30"),
    axis.text.x.top  = element_text(size=8,  color="grey30"),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1A_LASSO_CV_curve.png"),
       p_cv, width=7, height=5, dpi=300)
ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1A_LASSO_CV_curve.svg"),
       p_cv, width=7, height=5)
cat("CV curve saved.\n")

library(ggplot2)
library(dplyr)

# ── Extract CV data ──────────────────────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_lasso_all$lambda),
  cvm        = cv_lasso_all$cvm,
  cvsd       = cv_lasso_all$cvsd,
  cvup       = cv_lasso_all$cvup,
  cvlo       = cv_lasso_all$cvlo,
  nzero      = cv_lasso_all$nzero
)

lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)

# nzero at each lambda line
nz_min <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_min_log))]
nz_1se <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_1se_log))]

# Top axis: evenly spaced ticks with correct nzero labels
n_ticks <- 10
tick_idx    <- round(seq(1, nrow(cv_df), length.out=n_ticks))
tick_breaks <- cv_df$log_lambda[tick_idx]
tick_labels <- as.character(cv_df$nzero[tick_idx])

# y range for annotation placement
y_top <- max(cv_df$cvup) + 0.02 * diff(range(cv_df$cvup))

p_cv <- ggplot(cv_df, aes(x=log_lambda, y=cvm)) +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), fill="#BFCFDF", alpha=0.5) +
  geom_line(color="#2166AC", linewidth=0.8) +
  geom_point(color="#2166AC", size=1.2, shape=16) +
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.8) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.8) +
  # lambda.min label — placed to the LEFT of the line
  annotate("text", x=lam_min_log - 0.08, y=y_top,
           label=paste0("lambda.min\n(", nz_min, " genes)"),
           color="#D73027", size=3.0, hjust=1.0, fontface="bold") +
  # lambda.1se label — placed to the RIGHT of the line, but ensure it's inside
  annotate("text", x=lam_1se_log + 0.08, y=y_top,
           label=paste0("lambda.1se\n(", nz_1se, " genes)"),
           color="#4DAC26", size=3.0, hjust=0.0, fontface="bold") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(min(cv_df$log_lambda) - 0.1, max(cv_df$log_lambda) + 0.6),
    sec.axis = sec_axis(~., name="Number of Non-zero Genes",
                        breaks = tick_breaks,
                        labels = tick_labels)
  ) +
  labs(
    title    = "LASSO Cox Regression: 10-fold Cross-Validation",
    subtitle = "GSE10846 DLBCL (n=412, events=163)",
    y        = "Partial Likelihood Deviance (+/-1 SE)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.title         = element_text(size=11),
    axis.text          = element_text(size=10),
    axis.title.x.top   = element_text(size=10, color="grey30", margin=margin(b=4)),
    axis.text.x.top    = element_text(size=8,  color="grey30"),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4),
    plot.margin        = margin(t=10, r=20, b=10, l=10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1A_LASSO_CV_curve.png"),
       p_cv, width=8, height=5.5, dpi=300)
cat("CV curve saved.\n")

# Check the actual nzero values across lambda path
cat("Lambda path nzero values (first 30):\n")
print(data.frame(
  idx        = 1:min(30, length(cv_lasso_all$lambda)),
  log_lambda = round(log(cv_lasso_all$lambda[1:30]), 3),
  nzero      = cv_lasso_all$nzero[1:30]
))
cat("\nTotal lambda steps:", length(cv_lasso_all$lambda), "\n")
cat("nzero range:", range(cv_lasso_all$nzero), "\n")

library(ggplot2)
library(dplyr)

# ── CV data ──────────────────────────────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_lasso_all$lambda),
  cvm        = cv_lasso_all$cvm,
  cvup       = cv_lasso_all$cvup,
  cvlo       = cv_lasso_all$cvlo,
  nzero      = cv_lasso_all$nzero
)

lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)
nz_min <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_min_log))]
nz_1se <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_1se_log))]

# ── Top axis: pick ~10 evenly spaced points with UNIQUE nzero transitions ────
# Use points where nzero changes (transitions) + endpoints
transitions <- c(1, which(diff(cv_df$nzero) != 0) + 1, nrow(cv_df))
transitions <- unique(transitions)
# Subsample to ~10 points
if(length(transitions) > 10) {
  transitions <- transitions[round(seq(1, length(transitions), length.out=10))]
}
tick_breaks <- cv_df$log_lambda[transitions]
tick_labels <- as.character(cv_df$nzero[transitions])
cat("Top axis ticks — log_lambda:", round(tick_breaks,2), "\n")
cat("Top axis ticks — nzero:     ", tick_labels, "\n")

# ── y limits with padding ────────────────────────────────────────────────────
y_min <- min(cv_df$cvlo) - 0.05 * diff(range(cv_df$cvm))
y_max <- max(cv_df$cvup) + 0.15 * diff(range(cv_df$cvm))
y_ann <- max(cv_df$cvup) + 0.08 * diff(range(cv_df$cvm))

# ── x limits: extend right to fit lambda.1se label ──────────────────────────
x_min <- min(cv_df$log_lambda) - 0.1
x_max <- max(cv_df$log_lambda) + 0.5   # extra room on right

p_cv <- ggplot(cv_df, aes(x=log_lambda, y=cvm)) +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), fill="#BFCFDF", alpha=0.45) +
  geom_line(color="#2166AC", linewidth=0.9) +
  geom_point(color="#2166AC", size=1.4, shape=16) +
  # Vertical lines
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.85) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.85) +
  # lambda.min annotation — left of line
  annotate("text", x=lam_min_log - 0.05, y=y_ann,
           label=paste0("lambda.min\n(", nz_min, " genes)"),
           color="#D73027", size=3.1, hjust=1.0, fontface="bold") +
  # lambda.1se annotation — right of line, well inside
  annotate("text", x=lam_1se_log + 0.05, y=y_ann,
           label=paste0("lambda.1se\n(", nz_1se, " genes)"),
           color="#4DAC26", size=3.1, hjust=0.0, fontface="bold") +
  # Axes
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(x_min, x_max),
    sec.axis = sec_axis(
      transform = ~.,
      name   = "Number of Non-zero Genes",
      breaks = tick_breaks,
      labels = tick_labels
    )
  ) +
  scale_y_continuous(limits=c(y_min, y_max)) +
  labs(
    title    = "LASSO Cox Regression: 10-fold Cross-Validation",
    subtitle = "GSE10846 DLBCL (n=412, events=163)",
    y        = "Partial Likelihood Deviance (+/- 1 SE)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.title         = element_text(size=11),
    axis.text          = element_text(size=10),
    axis.title.x.top   = element_text(size=10, color="grey35", margin=margin(b=5)),
    axis.text.x.top    = element_text(size=9,  color="grey35"),
    panel.grid.major.y = element_line(color="grey92", linewidth=0.4),
    plot.margin        = margin(t=8, r=15, b=8, l=10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1A_LASSO_CV_curve.png"),
       p_cv, width=8.5, height=5.5, dpi=300)
cat("Saved.\n")

library(ggplot2)
library(dplyr)

# ── CV data ──────────────────────────────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_lasso_all$lambda),
  cvm        = cv_lasso_all$cvm,
  cvup       = cv_lasso_all$cvup,
  cvlo       = cv_lasso_all$cvlo,
  nzero      = cv_lasso_all$nzero
)

lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)
nz_min <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_min_log))]
nz_1se <- cv_df$nzero[which.min(abs(cv_df$log_lambda - lam_1se_log))]
cat("lambda.min: log=", round(lam_min_log,3), " nzero=", nz_min, "\n")
cat("lambda.1se: log=", round(lam_1se_log,3), " nzero=", nz_1se, "\n")

# ── Top axis: pick transition points ─────────────────────────────────────────
transitions <- c(1, which(diff(cv_df$nzero) != 0) + 1, nrow(cv_df))
transitions <- unique(transitions)
if(length(transitions) > 10) {
  transitions <- transitions[round(seq(1, length(transitions), length.out=10))]
}
tick_breaks <- cv_df$log_lambda[transitions]
tick_labels <- as.character(cv_df$nzero[transitions])

# ── y/x limits ───────────────────────────────────────────────────────────────
y_range <- diff(range(cv_df$cvm))
y_min   <- min(cv_df$cvlo) - 0.05 * y_range
y_max   <- max(cv_df$cvup) + 0.18 * y_range
y_ann   <- max(cv_df$cvup) + 0.10 * y_range

# Key: lambda.1se is at the RIGHT (large lambda, few genes)
# Place its label to the LEFT of the line to avoid right-edge clipping
x_min <- min(cv_df$log_lambda) - 0.2
x_max <- max(cv_df$log_lambda) + 0.2   # minimal right extension

p_cv <- ggplot(cv_df, aes(x=log_lambda, y=cvm)) +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), fill="#BFCFDF", alpha=0.45) +
  geom_line(color="#2166AC", linewidth=0.9) +
  geom_point(color="#2166AC", size=1.4, shape=16) +
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.85) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.85) +
  # lambda.min: to the LEFT of line (it's in the middle-left area)
  annotate("text", x=lam_min_log - 0.06, y=y_ann,
           label=paste0("lambda.min\n(", nz_min, " genes)"),
           color="#D73027", size=3.0, hjust=1.0, fontface="bold") +
  # lambda.1se: to the LEFT of line (it's near the right edge)
  annotate("text", x=lam_1se_log - 0.06, y=y_ann,
           label=paste0("lambda.1se\n(", nz_1se, " genes)"),
           color="#4DAC26", size=3.0, hjust=1.0, fontface="bold") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(x_min, x_max),
    sec.axis = sec_axis(
      transform = ~.,
      name   = "Number of Non-zero Genes",
      breaks = tick_breaks,
      labels = tick_labels
    )
  ) +
  scale_y_continuous(limits=c(y_min, y_max)) +
  labs(
    title    = "LASSO Cox Regression: 10-fold Cross-Validation",
    subtitle = "GSE10846 DLBCL (n=412, events=163)",
    y        = "Partial Likelihood Deviance (+/- 1 SE)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.title         = element_text(size=11),
    axis.text          = element_text(size=10),
    axis.title.x.top   = element_text(size=10, color="grey35", margin=margin(b=5)),
    axis.text.x.top    = element_text(size=9,  color="grey35"),
    panel.grid.major.y = element_line(color="grey92", linewidth=0.4),
    plot.margin        = margin(t=8, r=20, b=8, l=10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1A_LASSO_CV_curve.png"),
       p_cv, width=8.5, height=5.5, dpi=300)
cat("Saved.\n")

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# ── Extract full coefficient path ────────────────────────────────────────────
# Fit full LASSO path (not CV)
fit_path <- glmnet(X_all, surv_all, family="cox", alpha=1, standardize=FALSE)

# Build long-format dataframe: gene x lambda
coef_mat  <- as.matrix(coef(fit_path))          # 25 genes x n_lambda
log_lam   <- log(fit_path$lambda)
gene_names <- rownames(coef_mat)

coef_long <- data.frame(
  log_lambda = rep(log_lam, each=nrow(coef_mat)),
  gene       = rep(gene_names, times=ncol(coef_mat)),
  coef       = as.vector(coef_mat)
)

# Identify genes selected at lambda.min (non-zero)
selected_genes <- names(lasso_coefs)   # 13 genes
coef_long$selected <- coef_long$gene %in% selected_genes

# Color palette: selected genes get distinct colors, others grey
n_sel <- length(selected_genes)
pal   <- c(brewer.pal(8,"Set1"), brewer.pal(5,"Dark2"))
sel_colors <- setNames(pal[1:n_sel], selected_genes)

coef_long$color_group <- ifelse(coef_long$selected, coef_long$gene, "Other")

# Vertical lines
lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)

# Label positions: last non-zero value for each selected gene
label_df <- coef_long %>%
  filter(selected, coef != 0) %>%
  group_by(gene) %>%
  slice_min(log_lambda, n=1) %>%   # leftmost (smallest lambda = most genes)
  ungroup()

p_traj <- ggplot() +
  # Grey lines for non-selected genes
  geom_line(data=coef_long %>% filter(!selected),
            aes(x=log_lambda, y=coef, group=gene),
            color="grey75", linewidth=0.4, alpha=0.7) +
  # Colored lines for selected genes
  geom_line(data=coef_long %>% filter(selected),
            aes(x=log_lambda, y=coef, group=gene, color=gene),
            linewidth=0.85) +
  # Vertical lines
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.8) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.8) +
  # Horizontal zero line
  geom_hline(yintercept=0, color="black", linewidth=0.3) +
  # Gene labels at left edge
  geom_text(data=label_df,
            aes(x=log_lambda - 0.05, y=coef, label=gene, color=gene),
            size=2.8, hjust=1.0, fontface="bold", show.legend=FALSE) +
  # Lambda annotations at top
  annotate("text", x=lam_min_log, y=max(abs(coef_mat))*0.95,
           label="lambda.min", color="#D73027", size=3.0, hjust=0.5, fontface="bold") +
  annotate("text", x=lam_1se_log, y=max(abs(coef_mat))*0.95,
           label="lambda.1se", color="#4DAC26", size=3.0, hjust=0.5, fontface="bold") +
  scale_color_manual(values=sel_colors, name="Gene") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(min(log_lam) - 1.5, max(log_lam) + 0.2),
    sec.axis = sec_axis(
      transform = ~.,
      name   = "Number of Non-zero Genes",
      breaks = log(fit_path$lambda[round(seq(1, length(fit_path$lambda), length.out=10))]),
      labels = fit_path$df[round(seq(1, length(fit_path$lambda), length.out=10))]
    )
  ) +
  labs(
    title    = "LASSO Cox Regression: Coefficient Trajectory",
    subtitle = "GSE10846 DLBCL (n=412) — 25-gene set",
    y        = "Coefficient"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.title         = element_text(size=11),
    axis.text          = element_text(size=10),
    axis.title.x.top   = element_text(size=10, color="grey35", margin=margin(b=5)),
    axis.text.x.top    = element_text(size=9,  color="grey35"),
    panel.grid.major.y = element_line(color="grey92", linewidth=0.4),
    legend.position    = "right",
    legend.text        = element_text(size=9),
    legend.title       = element_text(size=10, face="bold"),
    plot.margin        = margin(t=8, r=10, b=8, l=10)
  ) +
  guides(color=guide_legend(ncol=1, keywidth=0.8, keyheight=0.7))

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig1B_LASSO_coef_trajectory.png"),
       p_traj, width=9, height=5.5, dpi=300)
cat("Trajectory plot saved.\n")

library(survival)
library(dplyr)

# ── Prepare clinical data with all covariates ────────────────────────────────
cox_df <- clin_all %>%
  mutate(
    Age60      = ifelse(age >= 60, 1, 0),           # Age >= 60 vs < 60
    Male       = ifelse(gender == "male", 1, 0),    # Male vs Female
    Stage_adv  = ifelse(stage >= 3, 1, 0),          # Stage III-IV vs I-II
    IPI_high   = ifelse(IPI >= 3, 1, 0),            # IPI high (3-5) vs low (0-2)
    RiskScore  = RiskScore                           # continuous
  ) %>%
  filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event))

cat("Cox dataset: n=", nrow(cox_df), " events=", sum(cox_df$OS_event), "\n")
cat("Missing values per covariate:\n")
print(colSums(is.na(cox_df[, c("RiskScore","Age60","Male","Stage_adv","IPI_high")])))

# ── Univariable Cox for each variable ────────────────────────────────────────
vars <- c("RiskScore","Age60","Male","Stage_adv","IPI_high")
uni_results <- lapply(vars, function(v) {
  f  <- as.formula(paste0("Surv(OS_years, OS_event) ~ ", v))
  df <- cox_df[!is.na(cox_df[[v]]), ]
  m  <- coxph(f, data=df)
  s  <- summary(m)
  data.frame(
    Variable = v,
    n        = nrow(df),
    events   = sum(df$OS_event),
    HR       = s$conf.int[1,1],
    CI_low   = s$conf.int[1,3],
    CI_high  = s$conf.int[1,4],
    pval     = s$coefficients[1,5],
    type     = "Univariable",
    stringsAsFactors = FALSE
  )
})
uni_df <- do.call(rbind, uni_results)

# ── Multivariable Cox ────────────────────────────────────────────────────────
# Use complete cases
cox_complete <- cox_df %>%
  filter(!is.na(Age60) & !is.na(Male) & !is.na(Stage_adv) & !is.na(IPI_high))
cat("\nComplete cases for multivariable Cox: n=", nrow(cox_complete),
    " events=", sum(cox_complete$OS_event), "\n")

multi_fit <- coxph(Surv(OS_years, OS_event) ~ RiskScore + Age60 + Male + Stage_adv + IPI_high,
                   data=cox_complete)
s_multi <- summary(multi_fit)
cat("\nMultivariable Cox summary:\n")
print(round(s_multi$conf.int, 4))
print(round(s_multi$coefficients[,"Pr(>|z|)"], 4))

multi_df <- data.frame(
  Variable = rownames(s_multi$conf.int),
  n        = nrow(cox_complete),
  events   = sum(cox_complete$OS_event),
  HR       = s_multi$conf.int[,1],
  CI_low   = s_multi$conf.int[,3],
  CI_high  = s_multi$conf.int[,4],
  pval     = s_multi$coefficients[,"Pr(>|z|)"],
  type     = "Multivariable",
  stringsAsFactors = FALSE
)

# Combine
forest_df <- rbind(uni_df, multi_df)
cat("\nForest data:\n"); print(forest_df)

library(ggplot2)
library(dplyr)

# ── Prepare forest data ──────────────────────────────────────────────────────
var_labels <- c(
  RiskScore = "Risk Score (continuous)",
  Age60     = "Age (>=60 vs <60)",
  Male      = "Sex (Male vs Female)",
  Stage_adv = "Stage (III-IV vs I-II)",
  IPI_high  = "IPI (High vs Low)"
)

# Build combined forest dataframe with section headers
uni_plot <- uni_df %>%
  mutate(
    label    = var_labels[Variable],
    sig      = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns"))),
    pval_str = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    hr_str   = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
  )

multi_plot <- multi_df %>%
  mutate(
    label    = var_labels[Variable],
    sig      = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns"))),
    pval_str = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval)),
    hr_str   = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
  )

# Assign y positions: univariable top, multivariable bottom, with gap
n_vars <- 5
uni_plot$y_pos   <- (n_vars + 2):(3)          # 7,6,5,4,3
multi_plot$y_pos <- (n_vars):(1)              # 5,4,3,2,1

# Significance color
sig_colors <- c("***"="#D73027", "**"="#FC8D59", "*"="#FEE090", "ns"="grey60")

# ── Build plot ───────────────────────────────────────────────────────────────
make_forest <- function(df, title_str, y_header) {
  ggplot(df, aes(y=y_pos, x=HR, xmin=CI_low, xmax=CI_high)) +
    # Reference line
    geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
    # CI bars
    geom_errorbarh(aes(color=sig), height=0.25, linewidth=0.8) +
    # Point estimates
    geom_point(aes(color=sig, size=sig), shape=18) +
    # Variable labels (left)
    geom_text(aes(x=0.05, label=label), hjust=0, size=3.3, color="black") +
    # HR (95% CI) text (right)
    geom_text(aes(x=12, label=hr_str), hjust=0, size=3.0, color="grey20") +
    # p-value text
    geom_text(aes(x=22, label=pval_str, color=sig), hjust=0, size=3.0, fontface="bold") +
    # Section header
    annotate("text", x=0.05, y=y_header, label=title_str,
             hjust=0, size=3.8, fontface="bold", color="grey20") +
    scale_color_manual(values=sig_colors, guide="none") +
    scale_size_manual(values=c("***"=4.5,"**"=4,"*"=3.5,"ns"=3), guide="none") +
    scale_x_log10(limits=c(0.04, 30),
                  breaks=c(0.25,0.5,1,2,4,8),
                  labels=c("0.25","0.5","1","2","4","8")) +
    scale_y_continuous(limits=c(0.5, y_header+0.5)) +
    labs(x=NULL, y=NULL) +
    theme_classic(base_size=11) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x  = element_text(size=9),
      panel.grid.major.x = element_line(color="grey93", linewidth=0.4),
      plot.margin  = margin(2,2,2,2)
    )
}

p_uni   <- make_forest(uni_plot,   "Univariable Analysis",   8)
p_multi <- make_forest(multi_plot, "Multivariable Analysis", 6)

# ── Combined plot with column headers ────────────────────────────────────────
library(gridExtra)
library(grid)

# Header row
header <- ggplot() +
  annotate("text", x=0.05, y=0.5, label="Variable",
           hjust=0, size=3.8, fontface="bold") +
  annotate("text", x=1,    y=0.5, label="HR (95% CI)",
           hjust=0.5, size=3.8, fontface="bold") +
  annotate("text", x=12,   y=0.5, label="HR (95% CI)",
           hjust=0, size=3.8, fontface="bold") +
  annotate("text", x=22,   y=0.5, label="P value",
           hjust=0, size=3.8, fontface="bold") +
  scale_x_log10(limits=c(0.04,30)) +
  scale_y_continuous(limits=c(0,1)) +
  theme_void()

# Full combined forest plot (single panel, both sections)
all_df <- bind_rows(
  uni_plot   %>% mutate(section="Univariable"),
  multi_plot %>% mutate(section="Multivariable")
)

# Reassign y positions: uni at top (rows 13-9), gap, multi (rows 7-3)
all_df$y_pos2 <- c(13,12,11,10,9,  7,6,5,4,3)

p_forest <- ggplot(all_df, aes(y=y_pos2, x=HR, xmin=CI_low, xmax=CI_high)) +
  # Reference line
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  # CI bars
  geom_errorbarh(aes(color=sig), height=0.3, linewidth=0.9) +
  # Points
  geom_point(aes(color=sig, size=sig), shape=18) +
  # Variable labels
  geom_text(aes(x=0.03, label=label), hjust=0, size=3.2, color="black") +
  # HR string
  geom_text(aes(x=11, label=hr_str), hjust=0, size=2.9, color="grey20") +
  # p-value
  geom_text(aes(x=21, label=pval_str, color=sig), hjust=0, size=2.9, fontface="bold") +
  # Section headers
  annotate("text", x=0.03, y=14.2, label="Univariable Analysis",
           hjust=0, size=4.0, fontface="bold", color="#2166AC") +
  annotate("text", x=0.03, y=8.2,  label="Multivariable Analysis",
           hjust=0, size=4.0, fontface="bold", color="#D73027") +
  # Column headers
  annotate("text", x=0.03, y=15.2, label="Variable",
           hjust=0, size=3.5, fontface="bold", color="grey30") +
  annotate("text", x=1.5,  y=15.2, label="Hazard Ratio",
           hjust=0.5, size=3.5, fontface="bold", color="grey30") +
  annotate("text", x=11,   y=15.2, label="HR (95% CI)",
           hjust=0, size=3.5, fontface="bold", color="grey30") +
  annotate("text", x=21,   y=15.2, label="P value",
           hjust=0, size=3.5, fontface="bold", color="grey30") +
  # Separator lines
  geom_hline(yintercept=c(14.7, 8.7), color="grey70", linewidth=0.5) +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_size_manual(values=c("***"=5,"**"=4.5,"*"=4,"ns"=3.5), guide="none") +
  scale_x_log10(
    limits = c(0.025, 35),
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25","0.5","1","2","4","8"),
    name   = "Hazard Ratio (log scale)"
  ) +
  scale_y_continuous(limits=c(2, 16), expand=c(0,0)) +
  labs(
    title    = "Cox Proportional Hazards Regression",
    subtitle = "GSE10846 DLBCL R-CHOP cohort (n=405, events=159)",
    y        = NULL
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.text.x        = element_text(size=10),
    axis.title.x       = element_text(size=11),
    panel.grid.major.x = element_line(color="grey93", linewidth=0.4),
    plot.margin        = margin(t=10, r=20, b=10, l=10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig2_Cox_forest_plot.png"),
       p_forest, width=10, height=7, dpi=300)
cat("Forest plot saved.\n")

# Save Cox results
write.csv(all_df %>% select(Variable, label, type=section, HR, CI_low, CI_high, pval, pval_str, hr_str),
          translate_path("/mnt/results/DLBCL_prognosis/tables/Cox_univariable_multivariable_results.csv"),
          row.names=FALSE)
cat("Cox results saved.\n")

library(ggplot2)
library(dplyr)

# ── Data ─────────────────────────────────────────────────────────────────────
var_labels <- c(
  RiskScore = "Risk Score (continuous)",
  Age60     = "Age (>=60 vs <60)",
  Male      = "Sex (Male vs Female)",
  Stage_adv = "Stage (III-IV vs I-II)",
  IPI_high  = "IPI (High vs Low)"
)
sig_colors <- c("***"="#D73027","**"="#FC8D59","*"="#E08000","ns"="grey55")

all_df <- bind_rows(
  uni_df   %>% mutate(section="Univariable"),
  multi_df %>% mutate(section="Multivariable")
) %>% mutate(
  label    = var_labels[Variable],
  sig      = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
             ifelse(pval < 0.05, "*", "ns"))),
  pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
  hr_str   = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
)
# y positions: uni=13..9, gap at 8, multi=7..3
all_df$y_pos <- c(13,12,11,10,9, 7,6,5,4,3)

# ── x-axis is log-scale for the forest part (0.1 to 10) ─────────────────────
# We use a TRICK: plot the forest on log scale 0.1-10,
# then place text annotations using coord_cartesian to extend the canvas
# Use a wide figure (12 inches) and place text at fixed x positions

# Text x positions (on log10 scale, so we use actual values):
# Variable label: x_label (left, outside plot via negative)
# HR CI text:     x_hr
# p-value text:   x_pv

# We'll use a two-panel approach: left panel = labels, right panel = forest+text
# Simpler: use a single ggplot with xlim extended far right, text in annotation space

x_forest_min <- 0.15
x_forest_max <- 12
x_hr_text    <- 14    # HR (95% CI) column
x_pv_text    <- 26    # p-value column
x_label_text <- 0.03  # variable name (left)

p_forest <- ggplot(all_df, aes(y=y_pos, x=HR, xmin=CI_low, xmax=CI_high)) +
  # Reference line
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  # CI bars
  geom_errorbar(aes(color=sig), width=0.3, linewidth=0.9,
                orientation="y") +
  # Diamond points
  geom_point(aes(color=sig, size=sig), shape=18) +
  # Variable labels (left of forest)
  geom_text(aes(x=x_label_text, label=label),
            hjust=0, size=3.2, color="black") +
  # HR (95% CI) text — right column 1
  geom_text(aes(x=x_hr_text, label=hr_str),
            hjust=0, size=2.9, color="grey20") +
  # p-value text — right column 2
  geom_text(aes(x=x_pv_text, label=pval_str, color=sig),
            hjust=0, size=2.9, fontface="bold") +
  # Section headers
  annotate("text", x=x_label_text, y=14.3,
           label="Univariable Analysis",
           hjust=0, size=4.0, fontface="bold", color="#2166AC") +
  annotate("text", x=x_label_text, y=8.3,
           label="Multivariable Analysis",
           hjust=0, size=4.0, fontface="bold", color="#D73027") +
  # Column headers
  annotate("text", x=x_label_text, y=15.5, label="Variable",
           hjust=0, size=3.5, fontface="bold", color="grey25") +
  annotate("text", x=1.5, y=15.5, label="Hazard Ratio",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  annotate("text", x=x_hr_text, y=15.5, label="HR (95% CI)",
           hjust=0, size=3.5, fontface="bold", color="grey25") +
  annotate("text", x=x_pv_text, y=15.5, label="P value",
           hjust=0, size=3.5, fontface="bold", color="grey25") +
  # Separator lines
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.5) +
  # Scales
  scale_color_manual(values=sig_colors, guide="none") +
  scale_size_manual(values=c("***"=5,"**"=4.5,"*"=4,"ns"=3.5), guide="none") +
  scale_x_log10(
    limits = c(x_label_text, 32),
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25","0.5","1","2","4","8"),
    name   = "Hazard Ratio (log scale)"
  ) +
  coord_cartesian(xlim=c(x_label_text, 32), ylim=c(2, 16.5), clip="off") +
  scale_y_continuous(expand=c(0,0)) +
  labs(
    title    = "Cox Proportional Hazards Regression",
    subtitle = "GSE10846 DLBCL (n=405, events=159)",
    y        = NULL
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title         = element_text(face="bold", size=13),
    plot.subtitle      = element_text(size=10, color="grey40"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.text.x        = element_text(size=10),
    axis.title.x       = element_text(size=11),
    panel.grid.major.x = element_line(color="grey93", linewidth=0.4),
    plot.margin        = margin(t=10, r=80, b=10, l=10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig2_Cox_forest_plot.png"),
       p_forest, width=12, height=7, dpi=300)
cat("Forest plot saved.\n")

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Data ─────────────────────────────────────────────────────────────────────
var_labels <- c(
  RiskScore = "Risk Score",
  Age60     = "Age (>=60 vs <60)",
  Male      = "Sex (Male vs Female)",
  Stage_adv = "Stage (III-IV vs I-II)",
  IPI_high  = "IPI (High vs Low)"
)
sig_colors <- c("***"="#D73027","**"="#FC8D59","*"="#E08000","ns"="grey55")

all_df <- bind_rows(
  uni_df   %>% mutate(section="Univariable"),
  multi_df %>% mutate(section="Multivariable")
) %>% mutate(
  label    = var_labels[Variable],
  sig      = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
             ifelse(pval < 0.05, "*", "ns"))),
  pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
  hr_str   = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
)
all_df$y_pos <- c(13,12,11,10,9, 7,6,5,4,3)

# ── Panel 1: Variable labels (left) ─────────────────────────────────────────
p_left <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=label), x=0, hjust=0, size=3.3, color="black") +
  # Section headers
  annotate("text", x=0, y=14.3, label="Univariable",
           hjust=0, size=3.8, fontface="bold", color="#2166AC") +
  annotate("text", x=0, y=8.3,  label="Multivariable",
           hjust=0, size=3.8, fontface="bold", color="#D73027") +
  annotate("text", x=0, y=15.5, label="Variable",
           hjust=0, size=3.5, fontface="bold", color="grey25") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=0, b=10, l=10))

# ── Panel 2: Forest plot (center) ────────────────────────────────────────────
p_forest <- ggplot(all_df, aes(y=y_pos, x=HR, xmin=CI_low, xmax=CI_high)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  geom_errorbar(aes(color=sig), width=0.3, linewidth=0.9, orientation="y") +
  geom_point(aes(color=sig, size=sig), shape=18) +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=1.5, y=15.5, label="Hazard Ratio",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_size_manual(values=c("***"=5,"**"=4.5,"*"=4,"ns"=3.5), guide="none") +
  scale_x_log10(limits=c(0.2, 12),
                breaks=c(0.25,0.5,1,2,4,8),
                labels=c("0.25","0.5","1","2","4","8"),
                name="Hazard Ratio (log scale)") +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  theme_classic(base_size=11) +
  theme(
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.text.x=element_text(size=9), axis.title.x=element_text(size=10),
    panel.grid.major.x=element_line(color="grey93", linewidth=0.4),
    plot.margin=margin(t=10, r=5, b=10, l=5)
  )

# ── Panel 3: HR (95% CI) text (right-1) ──────────────────────────────────────
p_hr <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=hr_str), x=0.5, hjust=0.5, size=2.9, color="grey20") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=0.5, y=15.5, label="HR (95% CI)",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=5, b=10, l=5))

# ── Panel 4: P-value text (right-2) ──────────────────────────────────────────
p_pv <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=pval_str, color=sig), x=0.5, hjust=0.5,
            size=2.9, fontface="bold") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=0.5, y=15.5, label="P value",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=10, b=10, l=5))

# ── Combine with patchwork ────────────────────────────────────────────────────
combined <- p_left + p_forest + p_hr + p_pv +
  plot_layout(widths=c(2.2, 2.5, 1.8, 1.0)) +
  plot_annotation(
    title    = "Cox Proportional Hazards Regression",
    subtitle = "GSE10846 DLBCL (n=405, events=159)",
    theme=theme(
      plot.title    = element_text(face="bold", size=13, hjust=0),
      plot.subtitle = element_text(size=10, color="grey40", hjust=0),
      plot.margin   = margin(10,10,10,10)
    )
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig2_Cox_forest_plot.png"),
       combined, width=12, height=7, dpi=300)
cat("Forest plot saved.\n")

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Panel 2: Forest plot — add labs(y=NULL) ───────────────────────────────────
p_forest <- ggplot(all_df, aes(y=y_pos, x=HR, xmin=CI_low, xmax=CI_high)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  geom_errorbar(aes(color=sig), width=0.3, linewidth=0.9, orientation="y") +
  geom_point(aes(color=sig, size=sig), shape=18) +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=1.5, y=15.5, label="Hazard Ratio",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_size_manual(values=c("***"=5,"**"=4.5,"*"=4,"ns"=3.5), guide="none") +
  scale_x_log10(limits=c(0.2, 12),
                breaks=c(0.25,0.5,1,2,4,8),
                labels=c("0.25","0.5","1","2","4","8"),
                name="Hazard Ratio (log scale)") +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(y=NULL) +
  theme_classic(base_size=11) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size=9),
    axis.title.x       = element_text(size=10),
    panel.grid.major.x = element_line(color="grey93", linewidth=0.4),
    plot.margin        = margin(t=10, r=5, b=10, l=5)
  )

# Rebuild all other panels (unchanged)
p_left <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=label), x=0, hjust=0, size=3.3, color="black") +
  annotate("text", x=0, y=14.3, label="Univariable",
           hjust=0, size=3.8, fontface="bold", color="#2166AC") +
  annotate("text", x=0, y=8.3,  label="Multivariable",
           hjust=0, size=3.8, fontface="bold", color="#D73027") +
  annotate("text", x=0, y=15.5, label="Variable",
           hjust=0, size=3.5, fontface="bold", color="grey25") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=0, b=10, l=10))

p_hr <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=hr_str), x=0.5, hjust=0.5, size=2.9, color="grey20") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=0.5, y=15.5, label="HR (95% CI)",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=5, b=10, l=5))

p_pv <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=pval_str, color=sig), x=0.5, hjust=0.5,
            size=2.9, fontface="bold") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=0.5, y=15.5, label="P value",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(plot.margin=margin(t=10, r=10, b=10, l=5))

combined <- p_left + p_forest + p_hr + p_pv +
  plot_layout(widths=c(2.2, 2.5, 1.8, 1.0)) +
  plot_annotation(
    title    = "Cox Proportional Hazards Regression",
    subtitle = "GSE10846 DLBCL (n=405, events=159)",
    theme=theme(
      plot.title    = element_text(face="bold", size=13, hjust=0),
      plot.subtitle = element_text(size=10, color="grey40", hjust=0),
      plot.margin   = margin(10,10,10,10)
    )
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig2_Cox_forest_plot.png"),
       combined, width=12, height=7, dpi=300)
cat("Saved.\n")

library(rms)
library(survival)
library(dplyr)

# ── Prepare data for rms ─────────────────────────────────────────────────────
nom_df <- cox_complete %>%
  select(OS_years, OS_event, RiskScore, IPI) %>%
  filter(!is.na(IPI) & !is.na(RiskScore) & OS_years > 0)

cat("Nomogram dataset: n=", nrow(nom_df), " events=", sum(nom_df$OS_event), "\n")
cat("RiskScore range:", round(range(nom_df$RiskScore), 3), "\n")
cat("IPI range:", range(nom_df$IPI), "\n")

# Set datadist for rms
dd <- datadist(nom_df)
options(datadist="dd")

# Fit Cox model with rms::cph
fit_cph <- cph(Surv(OS_years, OS_event) ~ RiskScore + IPI,
               data=nom_df, x=TRUE, y=TRUE, surv=TRUE, time.inc=1)
cat("\nrms cph model summary:\n")
print(fit_cph)

library(rms)

# Survival functions at 1, 3, 5 years
surv_1 <- Survival(fit_cph)(1)
surv_3 <- Survival(fit_cph)(3)
surv_5 <- Survival(fit_cph)(5)

cat("Median survival times used for nomogram:\n")
cat("1-year baseline survival:", round(mean(surv_1), 3), "\n")
cat("3-year baseline survival:", round(mean(surv_3), 3), "\n")
cat("5-year baseline survival:", round(mean(surv_5), 3), "\n")

# Build nomogram
nom <- nomogram(fit_cph,
                fun = list(
                  function(x) 1 - Survival(fit_cph)(1)(x),
                  function(x) 1 - Survival(fit_cph)(3)(x),
                  function(x) 1 - Survival(fit_cph)(5)(x)
                ),
                funlabel = c("1-Year Mortality", "3-Year Mortality", "5-Year Mortality"),
                lp = FALSE)

# Save as PNG
png(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig3A_Nomogram.png"),
    width=3000, height=1800, res=300)
par(mar=c(2, 2, 3, 2), cex=0.9)
plot(nom,
     xfrac=0.35,
     label.every=1,
     col.grid=grey(c(0.8, 0.95)),
     main="Nomogram for 1/3/5-Year Overall Survival\n(GSE10846 DLBCL, n=405)")
dev.off()
cat("Nomogram saved.\n")

library(rms)

# Correct rms syntax: Survival(fit)(times, lp)
# For nomogram fun= argument, use a function of linear predictor
# Get baseline survival at each time point
S1 <- survest(fit_cph, times=1, what="survival")$surv
S3 <- survest(fit_cph, times=3, what="survival")$surv
S5 <- survest(fit_cph, times=5, what="survival")$surv
cat("Baseline S(1):", round(median(S1, na.rm=TRUE), 3),
    " S(3):", round(median(S3, na.rm=TRUE), 3),
    " S(5):", round(median(S5, na.rm=TRUE), 3), "\n")

# Get baseline hazard at each time point using median lp=0
# Use the correct approach: define fun as S0^exp(lp)
# Get S0 at each time
S0_1 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=1)$surv
S0_3 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=3)$surv
S0_5 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=5)$surv
cat("S0(1):", round(S0_1,3), " S0(3):", round(S0_3,3), " S0(5):", round(S0_5,3), "\n")

# Nomogram fun: input is linear predictor (lp), output is 1-S(t)
# S(t|lp) = S0(t)^exp(lp - center)
center <- fit_cph$center
cat("Model center:", round(center, 4), "\n")

f1 <- function(lp) 1 - S0_1^exp(lp - center)
f3 <- function(lp) 1 - S0_3^exp(lp - center)
f5 <- function(lp) 1 - S0_5^exp(lp - center)

# Build nomogram
nom <- nomogram(fit_cph,
                fun      = list(f1, f3, f5),
                funlabel = c("1-Year Mortality", "3-Year Mortality", "5-Year Mortality"),
                fun.at   = list(c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                                c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                                c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)),
                lp = FALSE)

png(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig3A_Nomogram.png"),
    width=3200, height=1800, res=300)
par(mar=c(2, 2, 3, 2), cex=0.95)
plot(nom,
     xfrac    = 0.32,
     label.every = 1,
     col.grid = grey(c(0.8, 0.95)),
     main     = "Nomogram: 1/3/5-Year Overall Survival Prediction\n(IPI + Risk Score, GSE10846 DLBCL n=405)")
dev.off()
cat("Nomogram saved.\n")

library(rms)
library(ggplot2)
library(dplyr)

set.seed(42)

# ── Bootstrap calibration at 1, 3, 5 years ──────────────────────────────────
cal_1 <- calibrate(fit_cph, u=1, B=500, cmethod="KM")
cal_3 <- calibrate(fit_cph, u=3, B=500, cmethod="KM")
cal_5 <- calibrate(fit_cph, u=5, B=500, cmethod="KM")

# Extract calibration data
extract_cal <- function(cal_obj, year_label) {
  d <- as.data.frame(cal_obj)
  # Columns: mean.predicted, KM, std.err, n
  data.frame(
    predicted = d$mean.predicted,
    observed  = d$KM,
    lower     = d$KM - 1.96 * d$std.err,
    upper     = d$KM + 1.96 * d$std.err,
    year      = year_label,
    stringsAsFactors = FALSE
  )
}

cal_df <- bind_rows(
  extract_cal(cal_1, "1-Year"),
  extract_cal(cal_3, "3-Year"),
  extract_cal(cal_5, "5-Year")
)

# Convert mortality to survival (calibrate gives survival probability)
# predicted = predicted survival, observed = KM survival
# We want: x=predicted survival, y=observed survival
cat("Calibration data sample:\n")
print(head(cal_df, 6))
cat("Predicted range:", round(range(cal_df$predicted, na.rm=TRUE), 3), "\n")
cat("Observed range:", round(range(cal_df$observed, na.rm=TRUE), 3), "\n")

library(rms)

cat("OS_years range:", round(range(nom_df$OS_years), 3), "\n")
cat("Quantiles:", round(quantile(nom_df$OS_years, c(0.1,0.25,0.5,0.75,0.9)), 2), "\n")

# Refit with time.inc=1 (years) explicitly
fit_cph2 <- cph(Surv(OS_years, OS_event) ~ RiskScore + IPI,
                data=nom_df, x=TRUE, y=TRUE, surv=TRUE, time.inc=1)

# Test calibrate at 1 year
set.seed(42)
cal_test <- calibrate(fit_cph2, u=1, B=100)
cat("Calibration at 1 year - class:", class(cal_test), "\n")
cat("Dim:", dim(cal_test), "\n")
print(head(as.data.frame(cal_test)))

library(survival)
library(ggplot2)
library(dplyr)

# ── Compute predicted survival at 1, 3, 5 years for each patient ─────────────
# Use standard coxph (not rms) for reliable survival prediction
cox_cal <- coxph(Surv(OS_years, OS_event) ~ RiskScore + IPI,
                 data=nom_df, x=TRUE)

# Predicted survival at each time point
pred_surv <- function(fit, newdata, times) {
  sf <- survfit(fit, newdata=newdata)
  # Extract survival at each requested time
  sapply(times, function(t) {
    idx <- which(sf$time <= t)
    if(length(idx)==0) return(rep(1, nrow(newdata)))
    surv_mat <- summary(sf, times=t, extend=TRUE)$surv
    surv_mat
  })
}

# Get predicted survival for all patients at 1, 3, 5 years
sf_all <- survfit(cox_cal, newdata=nom_df)
get_pred_at_t <- function(sf, t) {
  s <- summary(sf, times=t, extend=TRUE)
  # Returns matrix: rows=times, cols=patients
  if(is.matrix(s$surv)) s$surv[1,] else s$surv
}

pred_1 <- get_pred_at_t(sf_all, 1)
pred_3 <- get_pred_at_t(sf_all, 3)
pred_5 <- get_pred_at_t(sf_all, 5)

cat("Predicted 1-yr survival range:", round(range(pred_1),3), "\n")
cat("Predicted 3-yr survival range:", round(range(pred_3),3), "\n")
cat("Predicted 5-yr survival range:", round(range(pred_5),3), "\n")

# ── Decile-based calibration ─────────────────────────────────────────────────
calibration_decile <- function(pred_surv, time_pt, df, n_groups=10) {
  df2 <- df %>% mutate(pred=pred_surv,
                       grp=ntile(pred, n_groups))
  results <- df2 %>% group_by(grp) %>% do({
    sub <- .
    km  <- survfit(Surv(OS_years, OS_event) ~ 1, data=sub)
    # KM survival at time_pt
    s   <- summary(km, times=time_pt, extend=TRUE)
    obs <- ifelse(length(s$surv)>0, s$surv[1], NA)
    se  <- ifelse(length(s$std.err)>0, s$std.err[1], NA)
    data.frame(
      mean_pred = mean(sub$pred),
      obs_surv  = obs,
      obs_lower = obs - 1.96*se,
      obs_upper = obs + 1.96*se,
      n         = nrow(sub)
    )
  }) %>% ungroup() %>%
    mutate(time=time_pt)
  results
}

cal_1yr <- calibration_decile(pred_1, 1, nom_df)
cal_3yr <- calibration_decile(pred_3, 3, nom_df)
cal_5yr <- calibration_decile(pred_5, 5, nom_df)

cal_all <- bind_rows(
  cal_1yr %>% mutate(label="1-Year OS"),
  cal_3yr %>% mutate(label="3-Year OS"),
  cal_5yr %>% mutate(label="5-Year OS")
) %>% mutate(label=factor(label, levels=c("1-Year OS","3-Year OS","5-Year OS")))

cat("\nCalibration data (1-year, first 5 rows):\n")
print(head(cal_1yr, 5))

library(ggplot2)
library(dplyr)

# Colors for 3 time points
time_colors <- c("1-Year OS"="#2166AC", "3-Year OS"="#4DAC26", "5-Year OS"="#D73027")

# Remove rows with NA observed survival
cal_plot <- cal_all %>% filter(!is.na(obs_surv))

# Compute mean absolute error per time point
mae_df <- cal_plot %>%
  group_by(label) %>%
  summarise(MAE = round(mean(abs(mean_pred - obs_surv)), 3),
            .groups="drop")
cat("Mean Absolute Error per time point:\n"); print(mae_df)

# Add MAE label for annotation
cal_plot <- cal_plot %>% left_join(mae_df, by="label") %>%
  mutate(ann_label = paste0(label, "\nMAE=", MAE))

# Annotation positions (top-left of each facet)
ann_df <- cal_plot %>%
  group_by(label, MAE) %>%
  summarise(.groups="drop") %>%
  mutate(x=0.08, y=0.97,
         txt=paste0("MAE = ", MAE))

p_cal <- ggplot(cal_plot, aes(x=mean_pred, y=obs_surv, color=label)) +
  # Perfect calibration diagonal
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50",
              linewidth=0.7) +
  # Error bars (95% CI)
  geom_errorbar(aes(ymin=pmax(obs_lower,0), ymax=pmin(obs_upper,1)),
                width=0.015, linewidth=0.6, alpha=0.7) +
  # Calibration points
  geom_point(size=3, shape=16) +
  # Smooth loess line
  geom_smooth(method="loess", se=FALSE, linewidth=1.0, span=1.2) +
  # MAE annotation
  geom_text(data=ann_df, aes(x=x, y=y, label=txt, color=label),
            hjust=0, vjust=1, size=3.2, fontface="bold", inherit.aes=FALSE) +
  facet_wrap(~label, ncol=3) +
  scale_color_manual(values=time_colors, guide="none") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2), name="Predicted Survival Probability") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2), name="Observed Survival Probability (KM)") +
  labs(
    title    = "Calibration Curves: Predicted vs Observed Survival",
    subtitle = "Decile-based calibration with 95% CI (GSE10846 DLBCL, n=405, B=10 groups)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title       = element_text(face="bold", size=13),
    plot.subtitle    = element_text(size=9.5, color="grey40"),
    strip.text       = element_text(face="bold", size=11),
    strip.background = element_rect(fill="grey95", color="grey70"),
    axis.title       = element_text(size=11),
    axis.text        = element_text(size=9),
    panel.grid.major = element_line(color="grey93", linewidth=0.4),
    plot.margin      = margin(10,10,10,10)
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig3B_Calibration_curves.png"),
       p_cal, width=11, height=4.5, dpi=300)

# Save calibration data
write.csv(cal_all, translate_path("/mnt/results/DLBCL_prognosis/tables/Calibration_data_1_3_5yr.csv"),
          row.names=FALSE)
cat("Calibration curves saved.\n")

library(dplyr)

# ── Published immune cell gene signatures ─────────────────────────────────────
# Based on Bindea et al. (2013) Immunity and TIMER2 marker sets
immune_sigs <- list(
  M2_Macrophage = c("CD163","MRC1","MSR1","MARCO","CD68","CSF1R","IL10","TGFB1",
                    "CCL18","CCL22","MMP9","VEGFA","PDCD1LG2","CD274","ARG1","FOLR2"),
  M1_Macrophage = c("CD80","CD86","NOS2","IL12A","IL12B","IL23A","TNF","IL6",
                    "CXCL9","CXCL10","CXCL11","IRF5","SOCS1","HLA-DRA","HLA-DRB1"),
  CD8_T_cell    = c("CD8A","CD8B","GZMB","GZMA","PRF1","IFNG","TBX21","EOMES",
                    "CXCR3","CCL5","FASLG","TNFRSF9","LAG3","HAVCR2","PDCD1"),
  CD4_T_cell    = c("CD4","IL2","IL4","IL5","IL13","GATA3","FOXP3","IL17A",
                    "RORC","BCL6","CXCR5","ICOS","CD40LG","CTLA4"),
  NK_cell       = c("NCAM1","NKG7","GNLY","KLRB1","KLRD1","KLRK1","NCR1","NCR3",
                    "FCGR3A","TYROBP","KIR2DL1","KIR3DL1","FGFBP2","CX3CR1"),
  B_cell        = c("CD19","MS4A1","CD79A","CD79B","PAX5","BLK","FCRL5","IGHM",
                    "IGHG1","IGKC","IGLC2","CR2","FCER2","CD22","CD27"),
  Treg          = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TNFRSF4",
                    "ENTPD1","NT5E","TIGIT","LAG3","HAVCR2","PDCD1","CCR8"),
  Neutrophil    = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9","MPO","ELANE",
                    "PRTN3","AZU1","CEACAM8","ITGAM","SELL","CXCL8"),
  Dendritic     = c("ITGAX","CLEC9A","XCR1","SIGLEC6","LILRA4","CLEC4C","IRF7",
                    "IRF8","BATF3","ID2","FLT3","THBD","CD1C","FCER1A")
)

# ── Manual ssGSEA implementation ──────────────────────────────────────────────
ssgsea_score <- function(expr_mat, gene_set) {
  # expr_mat: genes x samples matrix
  # Returns vector of ssGSEA scores per sample
  genes_present <- intersect(gene_set, rownames(expr_mat))
  if(length(genes_present) < 3) return(rep(NA, ncol(expr_mat)))
  
  apply(expr_mat, 2, function(sample_expr) {
    # Rank genes
    ranks <- rank(sample_expr, ties.method="average")
    n     <- length(ranks)
    # Indicator for gene set membership
    in_set <- names(ranks) %in% genes_present
    # Sorted ranks
    sorted_idx <- order(ranks, decreasing=TRUE)
    in_set_sorted <- in_set[sorted_idx]
    ranks_sorted  <- ranks[sorted_idx]
    # Running sum
    n_set <- sum(in_set)
    n_rest <- n - n_set
    step_up   <- ifelse(in_set_sorted, ranks_sorted^0.75 / sum(ranks_sorted[in_set_sorted]^0.75), 0)
    step_down <- ifelse(!in_set_sorted, 1/n_rest, 0)
    running_sum <- cumsum(step_up - step_down)
    # ssGSEA score = area under running sum
    sum(running_sum)
  })
}

# Use full expression matrix (all 420 samples, all genes)
cat("Running ssGSEA on", ncol(expr_mat), "samples...\n")
immune_scores <- sapply(names(immune_sigs), function(cell_type) {
  genes_found <- intersect(immune_sigs[[cell_type]], rownames(expr_mat))
  cat(cell_type, "- genes found:", length(genes_found), "/", length(immune_sigs[[cell_type]]), "\n")
  ssgsea_score(expr_mat, immune_sigs[[cell_type]])
})

immune_df <- as.data.frame(immune_scores)
immune_df$sample_id <- colnames(expr_mat)
cat("\nImmune scores computed for", nrow(immune_df), "samples\n")
cat("Score ranges:\n")
print(round(apply(immune_scores, 2, range, na.rm=TRUE), 2))

library(dplyr)
library(stringr)

# Build full gene-symbol expression matrix (all probes -> best probe per gene)
cat("Building full gene-symbol expression matrix...\n")

# Use probe_map (all probes with gene symbols) from earlier
probe_map_full <- probe_map %>%
  mutate(mean_expr = rowMeans(expr_mat[ID, , drop=FALSE])) %>%
  group_by(Gene_Symbol) %>%
  slice_max(mean_expr, n=1, with_ties=FALSE) %>%
  ungroup()

cat("Total unique genes:", nrow(probe_map_full), "\n")

# Build full gene x sample matrix
expr_full <- expr_mat[probe_map_full$ID, , drop=FALSE]
rownames(expr_full) <- probe_map_full$Gene_Symbol

cat("Full expression matrix:", nrow(expr_full), "genes x", ncol(expr_full), "samples\n")

# Check immune signature gene coverage
for(ct in names(immune_sigs)) {
  found <- intersect(immune_sigs[[ct]], rownames(expr_full))
  cat(ct, ":", length(found), "/", length(immune_sigs[[ct]]), "genes found\n")
}

library(dplyr)

cat("Running ssGSEA on full matrix (22880 genes x 420 samples)...\n")
immune_scores <- sapply(names(immune_sigs), function(cell_type) {
  ssgsea_score(expr_full, immune_sigs[[cell_type]])
})

immune_df <- as.data.frame(immune_scores)
immune_df$sample_id <- colnames(expr_full)

cat("Done. Score ranges:\n")
print(round(apply(immune_scores, 2, range, na.rm=TRUE), 1))

# Merge with risk scores
risk_immune <- clin_all %>%
  select(sample_id, RiskScore, RiskGroup) %>%
  inner_join(immune_df, by="sample_id")

cat("\nMerged dataset: n=", nrow(risk_immune), "\n")

# Spearman correlations with RiskScore
cor_results <- sapply(names(immune_sigs), function(ct) {
  test <- cor.test(risk_immune$RiskScore, risk_immune[[ct]],
                   method="spearman", exact=FALSE)
  c(rho=test$estimate, pval=test$p.value)
})
cor_df <- as.data.frame(t(cor_results))
cor_df$cell_type <- rownames(cor_df)
cor_df$sig <- ifelse(cor_df$pval < 0.001, "***",
              ifelse(cor_df$pval < 0.01, "**",
              ifelse(cor_df$pval < 0.05, "*", "ns")))
cat("\nSpearman correlations with RiskScore:\n")
print(cor_df %>% arrange(rho) %>%
      mutate(rho=round(rho,3), pval=formatC(pval, format="e", digits=2)) %>%
      select(cell_type, rho, pval, sig))

library(dplyr)

# Fix: cor.test returns named vector with "rho.rho" issue
cor_list <- lapply(names(immune_sigs), function(ct) {
  test <- cor.test(risk_immune$RiskScore, risk_immune[[ct]],
                   method="spearman", exact=FALSE)
  data.frame(
    cell_type = ct,
    rho       = as.numeric(test$estimate),
    pval      = test$p.value,
    stringsAsFactors = FALSE
  )
})
cor_df <- do.call(rbind, cor_list) %>%
  mutate(
    sig      = ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
               ifelse(pval < 0.05, "*", "ns"))),
    pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))
  ) %>%
  arrange(rho)

cat("Spearman correlations with RiskScore:\n")
print(cor_df %>% select(cell_type, rho=rho, pval_str, sig) %>%
      mutate(rho=round(rho,3)))

# Save immune scores + risk scores
write.csv(risk_immune,
          translate_path("/mnt/results/DLBCL_prognosis/tables/Immune_infiltration_scores.csv"),
          row.names=FALSE)
write.csv(cor_df,
          translate_path("/mnt/results/DLBCL_prognosis/tables/Immune_RiskScore_correlations.csv"),
          row.names=FALSE)
cat("\nSaved.\n")

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Helper: scatter plot with regression line ─────────────────────────────────
scatter_immune <- function(df, x_col, y_col, color, title_str, rho, pval_str) {
  ggplot(df, aes_string(x=x_col, y=y_col)) +
    geom_point(color=color, alpha=0.45, size=1.5) +
    geom_smooth(method="lm", se=TRUE, color=color, fill=paste0(color,"33"),
                linewidth=1.0) +
    annotate("text", x=Inf, y=Inf,
             label=paste0("Spearman r = ", round(rho,3), "\n", "P ", pval_str),
             hjust=1.1, vjust=1.3, size=3.5, fontface="bold",
             color=ifelse(rho>0,"#D73027","#2166AC")) +
    labs(title=title_str, x="Risk Score", y=paste0(gsub("_"," ",y_col), "\n(ssGSEA score)")) +
    theme_classic(base_size=11) +
    theme(
      plot.title   = element_text(face="bold", size=11),
      axis.title   = element_text(size=10),
      axis.text    = element_text(size=9),
      panel.grid.major = element_line(color="grey93", linewidth=0.3)
    )
}

# ── Key scatter plots ─────────────────────────────────────────────────────────
p_m2 <- scatter_immune(risk_immune, "RiskScore", "M2_Macrophage",
                        "#D73027", "M2 Macrophage vs Risk Score",
                        0.255, "< 0.001")

p_cd8 <- scatter_immune(risk_immune, "RiskScore", "CD8_T_cell",
                         "#2166AC", "CD8+ T Cell vs Risk Score",
                         0.008, "= 0.867")

p_m1 <- scatter_immune(risk_immune, "RiskScore", "M1_Macrophage",
                        "#4DAC26", "M1 Macrophage vs Risk Score",
                        -0.155, "= 0.002")

p_bcell <- scatter_immune(risk_immune, "RiskScore", "B_cell",
                           "#984EA3", "B Cell vs Risk Score",
                           0.236, "< 0.001")

# ── Correlation bar chart (all cell types) ────────────────────────────────────
cor_df_plot <- cor_df %>%
  mutate(
    cell_label = gsub("_", " ", cell_type),
    direction  = ifelse(rho > 0, "Positive", "Negative"),
    sig_label  = paste0(cell_label, " (", sig, ")")
  ) %>%
  arrange(rho) %>%
  mutate(cell_label = factor(cell_label, levels=cell_label))

p_bar <- ggplot(cor_df_plot, aes(x=rho, y=cell_label, fill=direction)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_vline(xintercept=0, color="black", linewidth=0.5) +
  geom_text(aes(label=sig, x=ifelse(rho>=0, rho+0.005, rho-0.005),
                hjust=ifelse(rho>=0, 0, 1)),
            size=3.5, fontface="bold") +
  scale_fill_manual(values=c("Positive"="#D73027","Negative"="#2166AC"),
                    name="Direction") +
  scale_x_continuous(limits=c(-0.25, 0.35), breaks=seq(-0.2,0.3,0.1)) +
  labs(title="Risk Score vs Immune Cell Infiltration",
       subtitle="Spearman correlation (GSE10846, n=412)",
       x="Spearman rho", y=NULL) +
  theme_classic(base_size=11) +
  theme(
    plot.title   = element_text(face="bold", size=12),
    plot.subtitle= element_text(size=9.5, color="grey40"),
    axis.text.y  = element_text(size=10),
    axis.text.x  = element_text(size=9),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color="grey93", linewidth=0.3)
  )

# ── Combine: bar chart + 4 scatter plots ─────────────────────────────────────
top_row    <- p_m2 + p_cd8 + p_m1 + p_bcell + plot_layout(ncol=4)
bottom_row <- p_bar

combined <- top_row / bottom_row +
  plot_layout(heights=c(1.2, 1)) +
  plot_annotation(
    title    = "Immune Infiltration Analysis",
    subtitle = "ssGSEA scores correlated with LASSO Risk Score (GSE10846 DLBCL)",
    theme=theme(
      plot.title    = element_text(face="bold", size=14, hjust=0),
      plot.subtitle = element_text(size=10, color="grey40", hjust=0),
      plot.margin   = margin(10,10,10,10)
    )
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis/figures/Fig4_Immune_infiltration.png"),
       combined, width=14, height=9, dpi=300)
cat("Immune infiltration figure saved.\n")

library(GEOquery)
library(dplyr)
library(stringr)

# ── Helper: parse a GEO series matrix and extract expression + clinical ────────
parse_geo_cohort <- function(filepath, cohort_name) {
  cat("\n=== Parsing", cohort_name, "===\n")
  gse <- getGEO(filename=filepath, GSEMatrix=TRUE, getGPL=FALSE)
  cat("Dimensions:", nrow(gse), "probes x", ncol(gse), "samples\n")
  
  pd <- pData(gse)
  cat("pData columns:\n"); print(colnames(pd))
  
  list(gse=gse, pdata=pd)
}

res87371  <- parse_geo_cohort(translate_path("/workspace/tmp_geo/GSE87371_series_matrix.txt.gz"),  "GSE87371")
res181063 <- parse_geo_cohort(translate_path("/workspace/tmp_geo/GSE181063_series_matrix.txt.gz"), "GSE181063")
res11318  <- parse_geo_cohort(translate_path("/workspace/tmp_geo/GSE11318_series_matrix.txt.gz"),  "GSE11318")

library(Biobase)
library(dplyr)
library(stringr)

# Load LASSO model info
lasso_data <- readRDS(translate_path("/workspace/lasso_results.rds"))
lasso_coefs <- lasso_data$coefs  # named vector: gene -> coefficient
train_means <- lasso_data$train_means
train_sds   <- lasso_data$train_sds
selected_genes <- names(lasso_coefs)

cat("Selected genes:", paste(selected_genes, collapse=", "), "\n")
cat("Coefficients:\n"); print(lasso_coefs)

library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# ── LASSO model parameters ─────────────────────────────────────────────────────
lasso_formula <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs   <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)
cat("Selected genes:", paste(selected_genes, collapse=", "), "\n")

# Training scaling parameters (from GSE10846 - need to recompute from risk scores)
# We'll load the training expression to get means/SDs
# Actually, the risk score was computed as sum(coef * scaled_expr)
# We need the original training means/SDs to scale validation data
# Load GSE10846 to recompute
cat("\nLoading GSE10846 to get training scaling parameters...\n")
gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), 
                   GSEMatrix=TRUE, getGPL=FALSE)

# Load GPL570 annotation
cat("Loading GPL570 annotation...\n")
gpl570 <- read.table(gzfile(translate_path("/workspace/tmp_geo/GPL570.soft.gz")), 
                     sep="\t", header=TRUE, comment.char="!", 
                     quote="", fill=TRUE, stringsAsFactors=FALSE)
# Keep ID and Gene Symbol
gpl570 <- gpl570[, c("ID", "Gene.Symbol")]
gpl570$Gene.Symbol <- str_trim(gpl570$Gene.Symbol)
# Remove probes with empty or multi-gene symbols
gpl570 <- gpl570[gpl570$Gene.Symbol != "" & !grepl("///", gpl570$Gene.Symbol), ]
cat("GPL570 probes with unique gene symbols:", nrow(gpl570), "\n")

library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# LASSO model parameters
lasso_formula <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs   <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

# Check GPL570 columns
gpl570_raw <- read.table(gzfile(translate_path("/workspace/tmp_geo/GPL570.soft.gz")), 
                         sep="\t", header=TRUE, comment.char="!", 
                         quote="", fill=TRUE, stringsAsFactors=FALSE, nrows=5)
cat("GPL570 columns:\n"); print(colnames(gpl570_raw))

library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# LASSO model parameters
lasso_formula <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs   <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

# Load GSE10846 with GPL annotation
cat("Loading GSE10846 with GPL570...\n")
gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"),
                   GSEMatrix=TRUE, getGPL=FALSE)

# Get expression matrix
expr10846 <- exprs(gse10846)
cat("GSE10846 dims:", nrow(expr10846), "x", ncol(expr10846), "\n")

# Use hgu133plus2.db for probe annotation
if (!requireNamespace("hgu133plus2.db", quietly=TRUE)) {
  BiocManager::install("hgu133plus2.db", ask=FALSE)
}
library(hgu133plus2.db)

# Map probes to gene symbols
probe_ids <- rownames(expr10846)
gene_syms <- mapIds(hgu133plus2.db, keys=probe_ids, 
                    column="SYMBOL", keytype="PROBEID", multiVals="first")
cat("Probes mapped:", sum(!is.na(gene_syms)), "/", length(gene_syms), "\n")

# Build probe-gene data frame
probe_gene <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene <- probe_gene[!is.na(probe_gene$gene), ]

# For each selected gene, find best probe (highest mean expression)
get_best_probe <- function(gene, expr_mat, probe_gene_df) {
  probes <- probe_gene_df$probe[probe_gene_df$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NULL)
  if (length(probes) == 1) return(probes)
  means <- rowMeans(expr_mat[probes, , drop=FALSE])
  probes[which.max(means)]
}

# Get best probes for selected genes
best_probes <- sapply(selected_genes, get_best_probe, 
                      expr_mat=expr10846, probe_gene_df=probe_gene)
cat("\nBest probes for selected genes:\n")
print(best_probes)
cat("Genes found:", sum(!is.null(best_probes) & !is.na(best_probes)), "/", length(selected_genes), "\n")

library(GEOquery); library(Biobase); library(dplyr); library(stringr); library(hgu133plus2.db)

# ── Reload LASSO params ────────────────────────────────────────────────────────
lasso_formula  <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs    <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

# ── Training data: GSE10846 ────────────────────────────────────────────────────
gse10846   <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr10846  <- exprs(gse10846)
probe_ids  <- rownames(expr10846)
gene_syms  <- mapIds(hgu133plus2.db, keys=probe_ids, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene <- probe_gene[!is.na(probe_gene$gene), ]

get_best_probe <- function(gene, expr_mat, pg) {
  probes <- pg$probe[pg$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NA)
  if (length(probes) == 1) return(probes)
  probes[which.max(rowMeans(expr_mat[probes, , drop=FALSE]))]
}

best_probes_train <- sapply(selected_genes, get_best_probe, expr_mat=expr10846, pg=probe_gene)

# Extract 13-gene expression matrix from training data
gene_expr_train <- t(expr10846[best_probes_train, ])
colnames(gene_expr_train) <- selected_genes

# Compute training means and SDs (for scaling validation data)
train_means <- colMeans(gene_expr_train)
train_sds   <- apply(gene_expr_train, 2, sd)
cat("Training means and SDs computed for", length(train_means), "genes\n")

# ── Helper: compute risk score for a new cohort ────────────────────────────────
compute_risk_score <- function(expr_mat, probe_gene_df, best_probes=NULL) {
  # If best_probes not provided, find them
  if (is.null(best_probes)) {
    bp <- sapply(selected_genes, get_best_probe, expr_mat=expr_mat, pg=probe_gene_df)
  } else {
    bp <- best_probes
  }
  
  found <- !is.na(bp) & bp %in% rownames(expr_mat)
  cat("  Genes found:", sum(found), "/", length(selected_genes), "\n")
  
  # Extract expression
  gene_expr <- t(expr_mat[bp[found], , drop=FALSE])
  colnames(gene_expr) <- selected_genes[found]
  
  # Scale using TRAINING means/SDs
  gene_expr_scaled <- sweep(gene_expr, 2, train_means[found], "-")
  gene_expr_scaled <- sweep(gene_expr_scaled, 2, train_sds[found], "/")
  
  # Compute risk score
  risk_score <- as.vector(gene_expr_scaled %*% lasso_coefs[found])
  names(risk_score) <- rownames(gene_expr)
  risk_score
}

# ── Process GSE87371 ───────────────────────────────────────────────────────────
cat("\n=== GSE87371 ===\n")
gse87371 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE87371_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr87371 <- exprs(gse87371)
pd87371   <- pData(gse87371)

# Clinical data
clin87371 <- data.frame(
  sample_id  = rownames(pd87371),
  OS_years   = as.numeric(pd87371$`os_time:ch1`),
  OS_event   = as.numeric(pd87371$`cens_os:ch1`),
  stringsAsFactors = FALSE
)
clin87371 <- clin87371[!is.na(clin87371$OS_years) & clin87371$OS_years > 0, ]
cat("  Samples with OS:", nrow(clin87371), "\n")

# Risk scores (GPL570 same platform)
rs87371 <- compute_risk_score(expr87371, probe_gene)
clin87371$RiskScore <- rs87371[clin87371$sample_id]
clin87371 <- clin87371[!is.na(clin87371$RiskScore), ]
clin87371$RiskGroup <- ifelse(clin87371$RiskScore >= median(clin87371$RiskScore), "High", "Low")
cat("  Final n:", nrow(clin87371), "| High:", sum(clin87371$RiskGroup=="High"), "Low:", sum(clin87371$RiskGroup=="Low"), "\n")
cat("  Events:", sum(clin87371$OS_event), "\n")

# ── Process GSE11318 ───────────────────────────────────────────────────────────
cat("\n=== GSE11318 ===\n")
gse11318 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE11318_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr11318 <- exprs(gse11318)
pd11318   <- pData(gse11318)

# Check clinical columns
cat("  Clinical columns:", paste(grep("ch1", colnames(pd11318), value=TRUE), collapse=", "), "\n")
cat("  Sample clinical info:\n")
print(head(pd11318[, grep("ch1|Clinical|Disease|Follow|Survival|OS|time|status", 
                           colnames(pd11318), ignore.case=TRUE)], 3))

library(GEOquery); library(Biobase); library(dplyr); library(stringr); library(hgu133plus2.db)

# ── Reload everything ──────────────────────────────────────────────────────────
lasso_formula  <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs    <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr10846 <- exprs(gse10846)
probe_ids <- rownames(expr10846)
gene_syms <- mapIds(hgu133plus2.db, keys=probe_ids, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene <- probe_gene[!is.na(probe_gene$gene), ]

get_best_probe <- function(gene, expr_mat, pg) {
  probes <- pg$probe[pg$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NA)
  if (length(probes) == 1) return(probes)
  probes[which.max(rowMeans(expr_mat[probes, , drop=FALSE]))]
}

best_probes_train <- sapply(selected_genes, get_best_probe, expr_mat=expr10846, pg=probe_gene)
gene_expr_train   <- t(expr10846[best_probes_train, ])
colnames(gene_expr_train) <- selected_genes
train_means <- colMeans(gene_expr_train)
train_sds   <- apply(gene_expr_train, 2, sd)

compute_risk_score <- function(expr_mat, pg) {
  bp    <- sapply(selected_genes, get_best_probe, expr_mat=expr_mat, pg=pg)
  found <- !is.na(bp) & bp %in% rownames(expr_mat)
  cat("  Genes found:", sum(found), "/", length(selected_genes), "\n")
  gene_expr        <- t(expr_mat[bp[found], , drop=FALSE])
  colnames(gene_expr) <- selected_genes[found]
  gene_expr_scaled <- sweep(gene_expr, 2, train_means[found], "-")
  gene_expr_scaled <- sweep(gene_expr_scaled, 2, train_sds[found], "/")
  rs <- as.vector(gene_expr_scaled %*% lasso_coefs[found])
  names(rs) <- rownames(gene_expr)
  rs
}

# ── GSE11318 ───────────────────────────────────────────────────────────────────
cat("=== GSE11318 ===\n")
gse11318  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE11318_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr11318 <- exprs(gse11318)
pd11318   <- pData(gse11318)

# Parse clinical info from characteristics columns
parse_field <- function(pd, pattern) {
  # Search all characteristics columns
  char_cols <- grep("characteristics_ch1", colnames(pd), value=TRUE)
  result <- rep(NA_character_, nrow(pd))
  for (col in char_cols) {
    vals <- pd[[col]]
    idx  <- grepl(pattern, vals, ignore.case=TRUE)
    result[idx] <- str_extract(vals[idx], paste0("(?<=", pattern, ":\\s?).*"))
  }
  result
}

clin11318 <- data.frame(
  sample_id = rownames(pd11318),
  OS_years  = as.numeric(str_extract(pd11318$`Clinical info:ch1`, "(?<=Follow up years: )[0-9.]+")),
  OS_event  = ifelse(grepl("DEAD", pd11318$`Clinical info:ch1`), 1, 0),
  stringsAsFactors = FALSE
)
clin11318 <- clin11318[!is.na(clin11318$OS_years) & clin11318$OS_years > 0, ]
cat("  Samples with OS:", nrow(clin11318), "| Events:", sum(clin11318$OS_event), "\n")

rs11318 <- compute_risk_score(expr11318, probe_gene)
clin11318$RiskScore <- rs11318[clin11318$sample_id]
clin11318 <- clin11318[!is.na(clin11318$RiskScore), ]
clin11318$RiskGroup <- ifelse(clin11318$RiskScore >= median(clin11318$RiskScore), "High", "Low")
cat("  Final n:", nrow(clin11318), "| High:", sum(clin11318$RiskGroup=="High"), "Low:", sum(clin11318$RiskGroup=="Low"), "\n")

# ── GSE181063 ─────────────────────────────────────────────────────────────────
cat("\n=== GSE181063 ===\n")
gse181063  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE181063_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr181063 <- exprs(gse181063)
pd181063   <- pData(gse181063)

cat("  Platform:", unique(pd181063$platform_id), "\n")
cat("  OS columns:", paste(grep("os_", colnames(pd181063), value=TRUE), collapse=", "), "\n")

clin181063 <- data.frame(
  sample_id = rownames(pd181063),
  OS_years  = as.numeric(pd181063$`os_followup_y:ch1`),
  OS_event  = as.numeric(pd181063$`os_status:ch1`),
  stringsAsFactors = FALSE
)
# Filter: valid OS, DLBCL diagnosis, curative intent, no QC fail
diag_col <- pd181063$`case_diagnostic_subtype_icdo3:ch1`
qc_col   <- pd181063$`qc_fail:ch1`
curative <- pd181063$`curative_intent:ch1`
cat("  Diagnosis values:", paste(unique(diag_col)[1:5], collapse=", "), "\n")
cat("  QC fail values:", paste(unique(qc_col), collapse=", "), "\n")
cat("  Curative values:", paste(unique(curative), collapse=", "), "\n")

keep <- !is.na(clin181063$OS_years) & clin181063$OS_years > 0 &
        !is.na(clin181063$OS_event) &
        (is.na(qc_col) | qc_col != "1")
clin181063 <- clin181063[keep, ]
cat("  Samples after filtering:", nrow(clin181063), "| Events:", sum(clin181063$OS_event, na.rm=TRUE), "\n")

library(GEOquery); library(Biobase); library(dplyr); library(stringr); library(hgu133plus2.db)

# ── Reload core objects ────────────────────────────────────────────────────────
lasso_formula  <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs    <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr10846 <- exprs(gse10846)
probe_ids <- rownames(expr10846)
gene_syms <- mapIds(hgu133plus2.db, keys=probe_ids, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene_affy <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene_affy <- probe_gene_affy[!is.na(probe_gene_affy$gene), ]

get_best_probe <- function(gene, expr_mat, pg) {
  probes <- pg$probe[pg$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NA)
  if (length(probes) == 1) return(probes)
  probes[which.max(rowMeans(expr_mat[probes, , drop=FALSE]))]
}

best_probes_train <- sapply(selected_genes, get_best_probe, expr_mat=expr10846, pg=probe_gene_affy)
gene_expr_train   <- t(expr10846[best_probes_train, ])
colnames(gene_expr_train) <- selected_genes
train_means <- colMeans(gene_expr_train)
train_sds   <- apply(gene_expr_train, 2, sd)

compute_risk_score_affy <- function(expr_mat) {
  bp    <- sapply(selected_genes, get_best_probe, expr_mat=expr_mat, pg=probe_gene_affy)
  found <- !is.na(bp) & bp %in% rownames(expr_mat)
  cat("  Genes found:", sum(found), "/", length(selected_genes), "\n")
  gene_expr        <- t(expr_mat[bp[found], , drop=FALSE])
  colnames(gene_expr) <- selected_genes[found]
  gene_expr_scaled <- sweep(gene_expr, 2, train_means[found], "-")
  gene_expr_scaled <- sweep(gene_expr_scaled, 2, train_sds[found], "/")
  rs <- as.vector(gene_expr_scaled %*% lasso_coefs[found])
  names(rs) <- rownames(gene_expr)
  rs
}

# ── GSE181063: Illumina GPL14951 ───────────────────────────────────────────────
cat("=== GSE181063 (Illumina GPL14951) ===\n")
gse181063  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE181063_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr181063 <- exprs(gse181063)
pd181063   <- pData(gse181063)

# Check probe IDs - Illumina probes often have gene symbol embedded or use ILMN_ IDs
cat("  Sample probe IDs:", paste(rownames(expr181063)[1:5], collapse=", "), "\n")
cat("  Dims:", nrow(expr181063), "x", ncol(expr181063), "\n")

# For Illumina, probe IDs may be gene symbols directly or ILMN IDs
# Check if any selected genes appear directly as row names
direct_match <- selected_genes[selected_genes %in% rownames(expr181063)]
cat("  Direct gene name matches:", length(direct_match), ":", paste(direct_match, collapse=", "), "\n")

# Try illuminaHumanv4.db
if (!requireNamespace("illuminaHumanv4.db", quietly=TRUE)) {
  BiocManager::install("illuminaHumanv4.db", ask=FALSE, quiet=TRUE)
}
library(illuminaHumanv4.db)
ilmn_probes <- rownames(expr181063)
ilmn_syms   <- mapIds(illuminaHumanv4.db, keys=ilmn_probes, 
                      column="SYMBOL", keytype="PROBEID", multiVals="first")
cat("  Illumina probes mapped:", sum(!is.na(ilmn_syms)), "/", length(ilmn_syms), "\n")

probe_gene_ilmn <- data.frame(probe=ilmn_probes, gene=ilmn_syms, stringsAsFactors=FALSE)
probe_gene_ilmn <- probe_gene_ilmn[!is.na(probe_gene_ilmn$gene), ]

# Check selected genes
found_ilmn <- selected_genes[selected_genes %in% probe_gene_ilmn$gene]
cat("  Selected genes found via illuminaHumanv4:", length(found_ilmn), ":", paste(found_ilmn, collapse=", "), "\n")

library(GEOquery); library(Biobase); library(dplyr); library(stringr)
library(hgu133plus2.db); library(illuminaHumanv4.db)

# ── Core objects ───────────────────────────────────────────────────────────────
lasso_formula  <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs    <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

# Training scaling params
gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr10846 <- exprs(gse10846)
probe_ids <- rownames(expr10846)
gene_syms <- mapIds(hgu133plus2.db, keys=probe_ids, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene_affy <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene_affy <- probe_gene_affy[!is.na(probe_gene_affy$gene), ]

get_best_probe <- function(gene, expr_mat, pg) {
  probes <- pg$probe[pg$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NA)
  if (length(probes) == 1) return(probes)
  probes[which.max(rowMeans(expr_mat[probes, , drop=FALSE]))]
}

best_probes_train <- sapply(selected_genes, get_best_probe, expr_mat=expr10846, pg=probe_gene_affy)
gene_expr_train   <- t(expr10846[best_probes_train, ])
colnames(gene_expr_train) <- selected_genes
train_means <- colMeans(gene_expr_train)
train_sds   <- apply(gene_expr_train, 2, sd)

# Generic risk score function
compute_rs <- function(expr_mat, pg) {
  bp    <- sapply(selected_genes, get_best_probe, expr_mat=expr_mat, pg=pg)
  found <- !is.na(bp) & bp %in% rownames(expr_mat)
  cat("  Genes found:", sum(found), "/", length(selected_genes), "\n")
  ge    <- t(expr_mat[bp[found], , drop=FALSE])
  colnames(ge) <- selected_genes[found]
  ge_sc <- sweep(ge, 2, train_means[found], "-")
  ge_sc <- sweep(ge_sc, 2, train_sds[found], "/")
  rs    <- as.vector(ge_sc %*% lasso_coefs[found])
  names(rs) <- rownames(ge)
  rs
}

# ── GSE87371 ───────────────────────────────────────────────────────────────────
cat("=== GSE87371 ===\n")
gse87371  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE87371_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr87371 <- exprs(gse87371); pd87371 <- pData(gse87371)
rs87371   <- compute_rs(expr87371, probe_gene_affy)
clin87371 <- data.frame(
  sample_id = rownames(pd87371),
  OS_years  = as.numeric(pd87371$`os_time:ch1`),
  OS_event  = as.numeric(pd87371$`cens_os:ch1`),
  stringsAsFactors = FALSE
)
clin87371 <- clin87371[!is.na(clin87371$OS_years) & clin87371$OS_years > 0, ]
clin87371$RiskScore <- rs87371[clin87371$sample_id]
clin87371 <- clin87371[!is.na(clin87371$RiskScore), ]
clin87371$RiskGroup <- ifelse(clin87371$RiskScore >= median(clin87371$RiskScore), "High", "Low")
clin87371$Cohort <- "GSE87371"
cat("  n=", nrow(clin87371), "events=", sum(clin87371$OS_event), "\n")

# ── GSE11318 ───────────────────────────────────────────────────────────────────
cat("=== GSE11318 ===\n")
gse11318  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE11318_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr11318 <- exprs(gse11318); pd11318 <- pData(gse11318)
rs11318   <- compute_rs(expr11318, probe_gene_affy)
clin11318 <- data.frame(
  sample_id = rownames(pd11318),
  OS_years  = as.numeric(str_extract(pd11318$`Clinical info:ch1`, "(?<=Follow up years: )[0-9.]+")),
  OS_event  = ifelse(grepl("DEAD", pd11318$`Clinical info:ch1`), 1, 0),
  stringsAsFactors = FALSE
)
clin11318 <- clin11318[!is.na(clin11318$OS_years) & clin11318$OS_years > 0, ]
clin11318$RiskScore <- rs11318[clin11318$sample_id]
clin11318 <- clin11318[!is.na(clin11318$RiskScore), ]
clin11318$RiskGroup <- ifelse(clin11318$RiskScore >= median(clin11318$RiskScore), "High", "Low")
clin11318$Cohort <- "GSE11318"
cat("  n=", nrow(clin11318), "events=", sum(clin11318$OS_event), "\n")

# ── GSE181063 ─────────────────────────────────────────────────────────────────
cat("=== GSE181063 ===\n")
gse181063  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE181063_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr181063 <- exprs(gse181063); pd181063 <- pData(gse181063)

# Illumina probe mapping
ilmn_probes <- rownames(expr181063)
ilmn_syms   <- mapIds(illuminaHumanv4.db, keys=ilmn_probes, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene_ilmn <- data.frame(probe=ilmn_probes, gene=ilmn_syms, stringsAsFactors=FALSE)
probe_gene_ilmn <- probe_gene_ilmn[!is.na(probe_gene_ilmn$gene), ]

rs181063 <- compute_rs(expr181063, probe_gene_ilmn)

# Filter to DLBCL NOS + curative intent + no QC fail
diag_col <- pd181063$`case_diagnostic_subtype_icdo3:ch1`
qc_col   <- pd181063$`qc_fail:ch1`
dlbcl_keep <- grepl("Diffuse large B-cell", diag_col, ignore.case=TRUE) &
              (is.na(qc_col) | qc_col == "0")

clin181063 <- data.frame(
  sample_id = rownames(pd181063),
  OS_years  = as.numeric(pd181063$`os_followup_y:ch1`),
  OS_event  = as.numeric(pd181063$`os_status:ch1`),
  stringsAsFactors = FALSE
)
clin181063 <- clin181063[dlbcl_keep & !is.na(clin181063$OS_years) & clin181063$OS_years > 0 &
                         !is.na(clin181063$OS_event), ]
clin181063$RiskScore <- rs181063[clin181063$sample_id]
clin181063 <- clin181063[!is.na(clin181063$RiskScore), ]
clin181063$RiskGroup <- ifelse(clin181063$RiskScore >= median(clin181063$RiskScore), "High", "Low")
clin181063$Cohort <- "GSE181063"
cat("  n=", nrow(clin181063), "events=", sum(clin181063$OS_event), "\n")

# ── GSE10846 training cohort ───────────────────────────────────────────────────
cat("=== GSE10846 (training) ===\n")
clin10846 <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/GSE10846_risk_scores.csv"))
clin10846 <- clin10846[!is.na(clin10846$OS_years) & clin10846$OS_years > 0, ]
clin10846$Cohort <- "GSE10846"
cat("  n=", nrow(clin10846), "events=", sum(clin10846$OS_event), "\n")

# Save all cohort data
saveRDS(list(
  GSE10846  = clin10846,
  GSE87371  = clin87371,
  GSE11318  = clin11318,
  GSE181063 = clin181063,
  train_means = train_means,
  train_sds   = train_sds
), translate_path("/workspace/cohort_data.rds"))

cat("\nAll GEO cohorts processed and saved.\n")

import requests
import pandas as pd
import numpy as np
import json, os

# ── TCGA-DLBC: fetch clinical data via GDC API ─────────────────────────────────
# Source: https://portal.gdc.cancer.gov/projects/TCGA-DLBC

BASE = "https://api.gdc.cancer.gov"

# Step 1: Get cases with OS data
cases_payload = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "project.project_id", "value": "TCGA-DLBC"}},
            {"op": "=", "content": {"field": "demographic.vital_status", "value": ["Alive", "Dead"]}}
        ]
    },
    "fields": "case_id,submitter_id,demographic.vital_status,demographic.days_to_death,"
              "diagnoses.days_to_last_follow_up,diagnoses.age_at_diagnosis",
    "format": "JSON",
    "size": "500"
}

resp = requests.post(f"{BASE}/cases", json=cases_payload, timeout=60)
data = resp.json()
cases = data["data"]["hits"]
print(f"TCGA-DLBC cases retrieved: {len(cases)}")

# Parse clinical
rows = []
for c in cases:
    demo = c.get("demographic", {})
    diag = c.get("diagnoses", [{}])
    diag = diag[0] if diag else {}
    
    vital = demo.get("vital_status", "")
    days_death = demo.get("days_to_death")
    days_follow = diag.get("days_to_last_follow_up")
    
    # OS time: days_to_death if dead, else days_to_last_follow_up
    if vital == "Dead" and days_death is not None:
        os_days = float(days_death)
        os_event = 1
    elif days_follow is not None:
        os_days = float(days_follow)
        os_event = 0
    else:
        continue
    
    if os_days <= 0:
        continue
        
    rows.append({
        "case_id": c["case_id"],
        "submitter_id": c["submitter_id"],
        "OS_years": os_days / 365.25,
        "OS_event": os_event
    })

clin_tcga = pd.DataFrame(rows)
print(f"Cases with valid OS: {len(clin_tcga)}")
print(f"Events: {clin_tcga['OS_event'].sum()}")
print(clin_tcga.head(3))

import requests
import pandas as pd
import numpy as np
import json, os, io, tarfile, gzip

BASE = "https://api.gdc.cancer.gov"

# Source: https://portal.gdc.cancer.gov/projects/TCGA-DLBC
# Get RNA-seq HTSeq FPKM files for TCGA-DLBC
files_payload = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-DLBC"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
            {"op": "=", "content": {"field": "data_format", "value": "TSV"}}
        ]
    },
    "fields": "file_id,file_name,cases.submitter_id,cases.case_id",
    "format": "JSON",
    "size": "200"
}

resp = requests.post(f"{BASE}/files", json=files_payload, timeout=60)
files_data = resp.json()
files = files_data["data"]["hits"]
print(f"RNA-seq files found: {len(files)}")

# Build file_id -> case mapping
file_case_map = {}
for f in files:
    fid = f["file_id"]
    cases_list = f.get("cases", [])
    if cases_list:
        file_case_map[fid] = cases_list[0]["submitter_id"]

print(f"Files with case mapping: {len(file_case_map)}")
print("Sample:", list(file_case_map.items())[:3])

import requests
import pandas as pd
import numpy as np
import json, os, io, tarfile, gzip

BASE = "https://api.gdc.cancer.gov"

# Selected genes from LASSO model
selected_genes = ["AAK1","BASP1","CEP192","CHCHD10","ELMO1","IFIT3",
                  "MTR","NAE1","PPP2R5A","PPRC1","RASSF4","TIMM17A","USP12"]

# Bulk download via GDC data endpoint
file_ids = list(file_case_map.keys())
print(f"Downloading {len(file_ids)} files...")

# Use manifest-based bulk download
download_payload = {"ids": file_ids}
resp = requests.post(
    f"{BASE}/data",
    json=download_payload,
    headers={"Content-Type": "application/json"},
    stream=True,
    timeout=300
)
print(f"Response status: {resp.status_code}")
print(f"Content-Type: {resp.headers.get('Content-Type','')}")

# Save tar archive
tar_path = translate_path("/workspace/tcga_dlbc_rnaseq.tar.gz")
with open(tar_path, "wb") as f:
    for chunk in resp.iter_content(chunk_size=1024*1024):
        f.write(chunk)

file_size = os.path.getsize(tar_path)
print(f"Downloaded: {file_size/1e6:.1f} MB")

import tarfile, gzip, io, os
import pandas as pd
import numpy as np

selected_genes = ["AAK1","BASP1","CEP192","CHCHD10","ELMO1","IFIT3",
                  "MTR","NAE1","PPP2R5A","PPRC1","RASSF4","TIMM17A","USP12"]

tar_path = translate_path("/workspace/tcga_dlbc_rnaseq.tar.gz")
expr_dict = {}  # file_id -> {gene: tpm}

with tarfile.open(tar_path, "r:gz") as tar:
    members = tar.getmembers()
    print(f"Files in archive: {len(members)}")
    # Show first few
    for m in members[:3]:
        print(" ", m.name)
    
    for member in members:
        if not member.name.endswith(".tsv") and not member.name.endswith(".tsv.gz"):
            continue
        
        # Extract file_id from path (format: file_id/filename.tsv)
        parts = member.name.split("/")
        if len(parts) < 2:
            continue
        file_id = parts[0]
        
        f = tar.extractfile(member)
        if f is None:
            continue
        
        content = f.read()
        if member.name.endswith(".gz"):
            content = gzip.decompress(content)
        
        # Parse TSV - STAR counts format has header lines starting with #
        lines = content.decode("utf-8").split("\n")
        # Find header line
        header_idx = 0
        for i, line in enumerate(lines):
            if line.startswith("gene_id"):
                header_idx = i
                break
        
        df = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])), sep="\t")
        
        # Get gene name column and TPM column
        gene_col = df.columns[0]  # gene_id (Ensembl)
        # Find TPM column
        tpm_cols = [c for c in df.columns if "tpm" in c.lower()]
        if not tpm_cols:
            # Try unstranded/stranded columns
            print(f"  Columns: {list(df.columns)}")
            continue
        tpm_col = tpm_cols[0]
        
        # Get gene_name column if available
        name_cols = [c for c in df.columns if "gene_name" in c.lower()]
        if name_cols:
            df["symbol"] = df[name_cols[0]]
        else:
            # Strip version from Ensembl ID
            df["symbol"] = df[gene_col].str.split(".").str[0]
        
        # Filter to selected genes
        gene_data = df[df["symbol"].isin(selected_genes)][["symbol", tpm_col]]
        gene_data = gene_data.groupby("symbol")[tpm_col].max()  # best probe
        
        expr_dict[file_id] = gene_data.to_dict()

print(f"\nFiles parsed: {len(expr_dict)}")
# Check gene coverage
if expr_dict:
    sample_key = list(expr_dict.keys())[0]
    print(f"Genes found in first file: {list(expr_dict[sample_key].keys())}")

import pandas as pd
import numpy as np

selected_genes = ["AAK1","BASP1","CEP192","CHCHD10","ELMO1","IFIT3",
                  "MTR","NAE1","PPP2R5A","PPRC1","RASSF4","TIMM17A","USP12"]

# LASSO coefficients
lasso_coefs = {
    "AAK1": -0.00286601381553037,
    "BASP1": -0.0920149545567974,
    "CEP192": 0.126325823454384,
    "CHCHD10": 0.0549231088915332,
    "ELMO1": 0.051337888713898,
    "IFIT3": 0.117207661256774,
    "MTR": -0.142201368454816,
    "NAE1": 0.0371959932375997,
    "PPP2R5A": -0.0643732234665509,
    "PPRC1": 0.263151788664124,
    "RASSF4": -0.0281013227812087,
    "TIMM17A": 0.0239819346834286,
    "USP12": -0.118848436868953
}

# Training scaling parameters (from GSE10846 Affymetrix data)
# These were computed in R - need to load them
# Since TCGA uses TPM (RNA-seq) vs microarray, we need to handle cross-platform scaling
# Best approach: z-score within TCGA cohort (cohort-specific scaling)
# This is standard practice for cross-platform validation

# Build expression matrix
expr_df = pd.DataFrame(expr_dict).T  # samples x genes
expr_df.index.name = "file_id"
expr_df = expr_df[selected_genes]  # ensure column order
print(f"Expression matrix: {expr_df.shape}")
print(f"Missing values: {expr_df.isnull().sum().sum()}")

# Log2 transform TPM (add 1 to avoid log(0))
expr_log = np.log2(expr_df + 1)
print(f"\nTPM range before log: {expr_df.values.min():.2f} - {expr_df.values.max():.2f}")
print(f"Log2(TPM+1) range: {expr_log.values.min():.2f} - {expr_log.values.max():.2f}")

# Z-score within TCGA cohort (standard cross-platform approach)
expr_scaled = (expr_log - expr_log.mean()) / expr_log.std()

# Compute risk scores
coef_vec = np.array([lasso_coefs[g] for g in selected_genes])
risk_scores = expr_scaled.values @ coef_vec
print(f"\nRisk score range: {risk_scores.min():.3f} - {risk_scores.max():.3f}")
print(f"Risk score median: {np.median(risk_scores):.3f}")

# Map file_id -> case submitter_id
rs_df = pd.DataFrame({
    "file_id": expr_df.index,
    "RiskScore": risk_scores
})
rs_df["submitter_id"] = rs_df["file_id"].map(file_case_map)
print(f"\nMapped to submitter_id: {rs_df['submitter_id'].notna().sum()}/{len(rs_df)}")

# Merge with clinical data
clin_tcga_merged = clin_tcga.merge(rs_df[["submitter_id","RiskScore"]], on="submitter_id", how="inner")
clin_tcga_merged["RiskGroup"] = np.where(
    clin_tcga_merged["RiskScore"] >= clin_tcga_merged["RiskScore"].median(), "High", "Low"
)
clin_tcga_merged["Cohort"] = "TCGA-DLBC"
print(f"\nTCGA-DLBC final: n={len(clin_tcga_merged)}, events={clin_tcga_merged['OS_event'].sum()}")
print(f"High: {(clin_tcga_merged['RiskGroup']=='High').sum()}, Low: {(clin_tcga_merged['RiskGroup']=='Low').sum()}")

# Save
clin_tcga_merged.to_csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"), index=False)
print("Saved TCGA-DLBC risk scores.")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# ── Load all cohort data ───────────────────────────────────────────────────────
cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

# GSE10846
clin10846 <- cohort_data$GSE10846 %>%
  select(sample_id, RiskScore, RiskGroup, OS_years, OS_event) %>%
  mutate(Cohort = "GSE10846\n(Training, n=412)")

# GSE87371
clin87371 <- cohort_data$GSE87371 %>%
  select(sample_id, RiskScore, RiskGroup, OS_years, OS_event) %>%
  mutate(Cohort = "GSE87371\n(Validation, n=221)")

# GSE11318
clin11318 <- cohort_data$GSE11318 %>%
  select(sample_id, RiskScore, RiskGroup, OS_years, OS_event) %>%
  mutate(Cohort = "GSE11318\n(Validation, n=199)")

# GSE181063
clin181063 <- cohort_data$GSE181063 %>%
  select(sample_id, RiskScore, RiskGroup, OS_years, OS_event) %>%
  mutate(Cohort = "GSE181063\n(Validation, n=882)")

# TCGA-DLBC
clin_tcga <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv")) %>%
  rename(sample_id = submitter_id) %>%
  select(sample_id, RiskScore, RiskGroup, OS_years, OS_event) %>%
  mutate(Cohort = "TCGA-DLBC\n(Validation, n=45)")

cohorts <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
cohort_names <- c("GSE10846\n(Training, n=412)",
                  "GSE87371\n(Validation, n=221)",
                  "GSE11318\n(Validation, n=199)",
                  "GSE181063\n(Validation, n=882)",
                  "TCGA-DLBC\n(Validation, n=45)")

# ── Color palette ──────────────────────────────────────────────────────────────
col_high <- "#E63946"   # red for high risk
col_low  <- "#457B9D"   # blue for low risk

# ── Generate one KM plot per cohort ───────────────────────────────────────────
make_km_plot <- function(df, title) {
  df$RiskGroup <- factor(df$RiskGroup, levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  
  # Log-rank test
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  pval <- 1 - pchisq(lr$chisq, df = 1)
  plab <- ifelse(pval < 0.001, "p < 0.001",
           ifelse(pval < 0.01,  sprintf("p = %.3f", pval),
                                sprintf("p = %.3f", pval)))
  
  # Cox HR
  cx  <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  hr  <- exp(coef(cx))
  ci  <- exp(confint(cx))
  hr_lab <- sprintf("HR = %.2f\n(95%% CI: %.2f–%.2f)", hr, ci[1], ci[2])
  
  # Annotation text
  annot <- paste0(plab, "\n", hr_lab)
  
  p <- ggsurvplot(
    fit,
    data          = df,
    palette       = c(col_low, col_high),
    legend.labs   = c("Low Risk", "High Risk"),
    legend.title  = "",
    legend        = c(0.75, 0.85),
    xlab          = "Time (years)",
    ylab          = "Overall Survival",
    title         = title,
    font.title    = c(11, "bold"),
    font.x        = 9,
    font.y        = 9,
    font.tickslab = 8,
    font.legend   = 8,
    risk.table    = TRUE,
    risk.table.height = 0.22,
    risk.table.fontsize = 2.8,
    risk.table.title = "No. at risk",
    tables.theme  = theme_cleantable(),
    ggtheme       = theme_classic(base_size = 9),
    conf.int      = TRUE,
    conf.int.alpha = 0.12,
    pval          = FALSE,   # we add manually
    surv.median.line = "hv",
    break.time.by = 5,
    xlim          = c(0, max(df$OS_years, na.rm=TRUE) * 1.02)
  )
  
  # Add annotation
  p$plot <- p$plot +
    annotate("text", x = Inf, y = 0.55,
             label = annot, hjust = 1.05, vjust = 0,
             size = 2.8, color = "black",
             fontface = "plain",
             family = "sans")
  
  p
}

# Build all 5 plots
plots <- mapply(make_km_plot, cohorts, cohort_names, SIMPLIFY = FALSE)

cat("KM plots generated for all 5 cohorts.\n")
for (i in seq_along(cohorts)) {
  df <- cohorts[[i]]
  lr <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  pv <- 1 - pchisq(lr$chisq, df=1)
  cx <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  hr <- exp(coef(cx)); ci <- exp(confint(cx))
  cat(sprintf("  %s: n=%d, events=%d, HR=%.2f (%.2f-%.2f), p=%.4f\n",
              gsub("\n.*","",cohort_names[i]),
              nrow(df), sum(df$OS_event), hr, ci[1], ci[2], pv))
}

library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(grid)

# ── Load all cohort data ───────────────────────────────────────────────────────
cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

prep_cohort <- function(df, label) {
  df <- df[, c("sample_id","RiskScore","RiskGroup","OS_years","OS_event")]
  df$Cohort <- label
  df
}

clin10846  <- prep_cohort(cohort_data$GSE10846,  "GSE10846")
clin87371  <- prep_cohort(cohort_data$GSE87371,  "GSE87371")
clin11318  <- prep_cohort(cohort_data$GSE11318,  "GSE11318")
clin181063 <- prep_cohort(cohort_data$GSE181063, "GSE181063")

tcga_raw   <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))
clin_tcga  <- data.frame(
  sample_id  = tcga_raw$submitter_id,
  RiskScore  = tcga_raw$RiskScore,
  RiskGroup  = tcga_raw$RiskGroup,
  OS_years   = tcga_raw$OS_years,
  OS_event   = tcga_raw$OS_event,
  Cohort     = "TCGA-DLBC",
  stringsAsFactors = FALSE
)

cohorts <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
titles  <- c("GSE10846 (Training, n=412)",
             "GSE87371 (Validation, n=221)",
             "GSE11318 (Validation, n=199)",
             "GSE181063 (Validation, n=882)",
             "TCGA-DLBC (Validation, n=45)")

col_high <- "#E63946"
col_low  <- "#457B9D"

# ── KM plot function ───────────────────────────────────────────────────────────
make_km_plot <- function(df, title) {
  df$RiskGroup <- factor(df$RiskGroup, levels = c("Low","High"))
  
  fit <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  pval <- 1 - pchisq(lr$chisq, df = 1)
  plab <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  
  cx  <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  hr  <- exp(coef(cx))
  ci  <- exp(confint(cx))
  hr_lab <- sprintf("HR = %.2f (%.2f\u2013%.2f)", hr, ci[1], ci[2])
  
  annot <- paste0(plab, "\n", hr_lab)
  
  xmax <- ceiling(max(df$OS_years, na.rm=TRUE))
  brks <- seq(0, xmax, by = ifelse(xmax > 15, 5, ifelse(xmax > 8, 4, 3)))
  
  p <- ggsurvplot(
    fit,
    data               = df,
    palette            = c(col_low, col_high),
    legend.labs        = c("Low Risk","High Risk"),
    legend.title       = "",
    legend             = c(0.78, 0.88),
    xlab               = "Time (years)",
    ylab               = "Overall Survival Probability",
    title              = title,
    font.title         = c(10, "bold"),
    font.x             = 8,
    font.y             = 8,
    font.tickslab      = 7.5,
    font.legend        = 8,
    risk.table         = TRUE,
    risk.table.height  = 0.25,
    risk.table.fontsize = 2.6,
    risk.table.title   = "No. at risk",
    tables.theme       = theme_cleantable(),
    ggtheme            = theme_classic(base_size = 9) +
                           theme(plot.title = element_text(hjust=0.5)),
    conf.int           = TRUE,
    conf.int.alpha     = 0.12,
    pval               = FALSE,
    surv.median.line   = "hv",
    break.time.by      = brks[2],
    xlim               = c(0, xmax)
  )
  
  p$plot <- p$plot +
    annotate("text", x = xmax * 0.98, y = 0.52,
             label = annot, hjust = 1, vjust = 0,
             size = 2.8, color = "black")
  p
}

plots <- lapply(seq_along(cohorts), function(i) make_km_plot(cohorts[[i]], titles[[i]]))

# ── Print summary stats ────────────────────────────────────────────────────────
cat("Cohort summary:\n")
for (i in seq_along(cohorts)) {
  df   <- cohorts[[i]]
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pv   <- 1 - pchisq(lr$chisq, df=1)
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr   <- exp(coef(cx)); ci <- exp(confint(cx))
  cat(sprintf("  %-12s n=%4d  events=%3d  HR=%.2f (%.2f-%.2f)  p=%.4f\n",
              cohorts[[i]]$Cohort[1], nrow(df), sum(df$OS_event), hr, ci[1], ci[2], pv))
}
cat("Plots ready.\n")

# Check GSE87371 - why is HR < 1?
df87 <- cohort_data$GSE87371
df87$RiskGroup <- factor(df87$RiskGroup, levels=c("Low","High"))

cat("GSE87371 risk score distribution:\n")
cat("  Range:", range(df87$RiskScore), "\n")
cat("  Median:", median(df87$RiskScore), "\n")
cat("  High group events:", sum(df87$OS_event[df87$RiskGroup=="High"]), "/", sum(df87$RiskGroup=="High"), "\n")
cat("  Low group events:", sum(df87$OS_event[df87$RiskGroup=="Low"]), "/", sum(df87$RiskGroup=="Low"), "\n")

# Median OS by group
cat("\nMedian OS by group:\n")
fit87 <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df87)
print(summary(fit87)$table[, c("records","events","median")])

# Check GSE10846 for comparison
df10 <- cohort_data$GSE10846
df10$RiskGroup <- factor(df10$RiskGroup, levels=c("Low","High"))
cat("\nGSE10846 for comparison:\n")
cat("  High group events:", sum(df10$OS_event[df10$RiskGroup=="High"]), "/", sum(df10$RiskGroup=="High"), "\n")
cat("  Low group events:", sum(df10$OS_event[df10$RiskGroup=="Low"]), "/", sum(df10$RiskGroup=="Low"), "\n")

# Check GSE87371 OS time distribution - maybe it's PFS not OS?
cat("\nGSE87371 OS time range:", range(df87$OS_years, na.rm=TRUE), "\n")
cat("GSE87371 event rate:", mean(df87$OS_event), "\n")

# Very high event rate (168/221 = 76%) suggests this may be a pre-rituximab cohort
# or the cens_os coding may be inverted (0=dead, 1=alive)
# Check: if cens_os=1 means censored (alive), then event=1 means dead - that's correct
# But let's verify by checking if high risk score patients have WORSE survival
cat("\nCorrelation of RiskScore with OS_event:", cor(df87$RiskScore, df87$OS_event), "\n")
cat("Correlation of RiskScore with OS_years:", cor(df87$RiskScore, df87$OS_years), "\n")

library(GEOquery); library(Biobase)

gse87371 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE87371_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
pd87371  <- pData(gse87371)

# Check cens_os values
cat("cens_os values (first 20):\n")
print(table(pd87371$`cens_os:ch1`))

cat("\nos_time range:", range(as.numeric(pd87371$`os_time:ch1`), na.rm=TRUE), "\n")

# Check a few samples: if cens_os=1 means event (dead) or censored
# Look at samples with short follow-up
short_fu <- pd87371[order(as.numeric(pd87371$`os_time:ch1`)), 
                    c("os_time:ch1","cens_os:ch1","pfs_time:ch1","cens_pfs:ch1","treatment:ch1")]
cat("\nShortest follow-up samples:\n")
print(head(short_fu, 10))

# Also check the GEO description for cens_os meaning
cat("\nSample characteristics for first sample:\n")
print(pd87371[1, grep("ch1", colnames(pd87371), value=TRUE)[1:10]])

library(GEOquery); library(Biobase); library(survival); library(survminer)
library(hgu133plus2.db)

# ── Reload core objects ────────────────────────────────────────────────────────
lasso_formula  <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"))
lasso_coefs    <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
selected_genes <- names(lasso_coefs)

gse10846 <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE10846_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr10846 <- exprs(gse10846)
probe_ids <- rownames(expr10846)
gene_syms <- mapIds(hgu133plus2.db, keys=probe_ids, column="SYMBOL", keytype="PROBEID", multiVals="first")
probe_gene_affy <- data.frame(probe=probe_ids, gene=gene_syms, stringsAsFactors=FALSE)
probe_gene_affy <- probe_gene_affy[!is.na(probe_gene_affy$gene), ]

get_best_probe <- function(gene, expr_mat, pg) {
  probes <- pg$probe[pg$gene == gene]
  probes <- probes[probes %in% rownames(expr_mat)]
  if (length(probes) == 0) return(NA)
  if (length(probes) == 1) return(probes)
  probes[which.max(rowMeans(expr_mat[probes, , drop=FALSE]))]
}

best_probes_train <- sapply(selected_genes, get_best_probe, expr_mat=expr10846, pg=probe_gene_affy)
gene_expr_train   <- t(expr10846[best_probes_train, ])
colnames(gene_expr_train) <- selected_genes
train_means <- colMeans(gene_expr_train)
train_sds   <- apply(gene_expr_train, 2, sd)

compute_rs <- function(expr_mat, pg) {
  bp    <- sapply(selected_genes, get_best_probe, expr_mat=expr_mat, pg=pg)
  found <- !is.na(bp) & bp %in% rownames(expr_mat)
  ge    <- t(expr_mat[bp[found], , drop=FALSE])
  colnames(ge) <- selected_genes[found]
  ge_sc <- sweep(ge, 2, train_means[found], "-")
  ge_sc <- sweep(ge_sc, 2, train_sds[found], "/")
  rs    <- as.vector(ge_sc %*% lasso_coefs[found])
  names(rs) <- rownames(ge); rs
}

# ── Fix GSE87371: cens_os=0 means DEAD, cens_os=1 means ALIVE (censored) ──────
cat("=== GSE87371 (corrected event coding) ===\n")
gse87371  <- getGEO(filename=translate_path("/workspace/tmp_geo/GSE87371_series_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE)
expr87371 <- exprs(gse87371); pd87371 <- pData(gse87371)

rs87371 <- compute_rs(expr87371, probe_gene_affy)

clin87371 <- data.frame(
  sample_id = rownames(pd87371),
  OS_years  = as.numeric(pd87371$`os_time:ch1`),
  OS_event  = 1 - as.numeric(pd87371$`cens_os:ch1`),  # FLIP: 0→1 (dead), 1→0 (alive)
  stringsAsFactors = FALSE
)
clin87371 <- clin87371[!is.na(clin87371$OS_years) & clin87371$OS_years > 0 & !is.na(clin87371$OS_event), ]
clin87371$RiskScore <- rs87371[clin87371$sample_id]
clin87371 <- clin87371[!is.na(clin87371$RiskScore), ]
clin87371$RiskGroup <- ifelse(clin87371$RiskScore >= median(clin87371$RiskScore), "High", "Low")
clin87371$Cohort <- "GSE87371"

cat("  n=", nrow(clin87371), "events=", sum(clin87371$OS_event), "\n")
cat("  High events:", sum(clin87371$OS_event[clin87371$RiskGroup=="High"]),
    "/ Low events:", sum(clin87371$OS_event[clin87371$RiskGroup=="Low"]), "\n")

# Verify direction
cx87 <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, 
              data=transform(clin87371, RiskGroup=factor(RiskGroup, levels=c("Low","High"))))
hr87 <- exp(coef(cx87)); ci87 <- exp(confint(cx87))
lr87 <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin87371)
pv87 <- 1 - pchisq(lr87$chisq, df=1)
cat(sprintf("  HR=%.2f (%.2f-%.2f), p=%.4f\n", hr87, ci87[1], ci87[2], pv87))

# Update cohort_data
cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))
cohort_data$GSE87371 <- clin87371
saveRDS(cohort_data, translate_path("/workspace/cohort_data.rds"))
cat("Updated cohort_data saved.\n")

library(survival); library(survminer); library(ggplot2); library(gridExtra); library(grid)

cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

prep_cohort <- function(df, label) {
  df <- df[, c("sample_id","RiskScore","RiskGroup","OS_years","OS_event")]
  df$Cohort <- label; df
}

clin10846  <- prep_cohort(cohort_data$GSE10846,  "GSE10846")
clin87371  <- prep_cohort(cohort_data$GSE87371,  "GSE87371")
clin11318  <- prep_cohort(cohort_data$GSE11318,  "GSE11318")
clin181063 <- prep_cohort(cohort_data$GSE181063, "GSE181063")
tcga_raw   <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))
clin_tcga  <- data.frame(
  sample_id=tcga_raw$submitter_id, RiskScore=tcga_raw$RiskScore,
  RiskGroup=tcga_raw$RiskGroup, OS_years=tcga_raw$OS_years,
  OS_event=tcga_raw$OS_event, Cohort="TCGA-DLBC", stringsAsFactors=FALSE
)

cohorts <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
titles  <- c("GSE10846 (Training, n=412)",
             "GSE87371 (Validation, n=221)",
             "GSE11318 (Validation, n=199)",
             "GSE181063 (Validation, n=882)",
             "TCGA-DLBC (Validation, n=45)")

col_high <- "#E63946"; col_low <- "#457B9D"

make_km_plot <- function(df, title) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval <- 1 - pchisq(lr$chisq, df=1)
  plab <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr   <- exp(coef(cx)); ci <- exp(confint(cx))
  annot <- sprintf("%s\nHR = %.2f (%.2f\u2013%.2f)", plab, hr, ci[1], ci[2])
  xmax  <- ceiling(max(df$OS_years, na.rm=TRUE))
  bby   <- if (xmax > 15) 5L else if (xmax > 8) 4L else 3L

  p <- ggsurvplot(
    fit, data=df,
    palette=c(col_low, col_high),
    legend.labs=c("Low Risk","High Risk"), legend.title="",
    legend=c(0.78, 0.88),
    xlab="Time (years)", ylab="Overall Survival Probability",
    title=title,
    font.title=c(10L,"bold"), font.x=9L, font.y=9L,
    font.tickslab=8L, font.legend=8L,
    risk.table=TRUE, risk.table.height=0.25,
    risk.table.fontsize=3L, risk.table.title="No. at risk",
    tables.theme=theme_cleantable(),
    ggtheme=theme_classic(base_size=9L) +
              theme(plot.title=element_text(hjust=0.5)),
    conf.int=TRUE, conf.int.alpha=0.12,
    pval=FALSE, surv.median.line="hv",
    break.time.by=bby, xlim=c(0, xmax)
  )
  p$plot <- p$plot +
    annotate("text", x=xmax*0.98, y=0.52,
             label=annot, hjust=1, vjust=0, size=3L, color="black")
  p
}

plots <- lapply(seq_along(cohorts), function(i) make_km_plot(cohorts[[i]], titles[[i]]))

# ── Arrange 5 panels: 3 top + 2 bottom centered ───────────────────────────────
outdir <- translate_path("/mnt/results/DLBCL_prognosis/figures")

# Save as combined figure using arrange_ggsurvplots
combined <- arrange_ggsurvplots(plots, nrow=2, ncol=3,
                                 print=FALSE, risk.table.height=0.25)

png(file.path(outdir, "Fig5_MultiCohort_KM_curves.png"),
    width=18, height=12, units="in", res=180)
print(combined)
dev.off()

svg(file.path(outdir, "Fig5_MultiCohort_KM_curves.svg"),
    width=18, height=12)
print(combined)
dev.off()

cat("Saved Fig5_MultiCohort_KM_curves.png and .svg\n")

# Print final summary
cat("\nFinal cohort summary:\n")
cat(sprintf("%-12s %5s %7s %20s %8s\n", "Cohort","n","Events","HR (95% CI)","p-value"))
cat(strrep("-", 60), "\n")
for (i in seq_along(cohorts)) {
  df <- cohorts[[i]]; df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  lr <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pv <- 1 - pchisq(lr$chisq, df=1)
  cx <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr <- exp(coef(cx)); ci <- exp(confint(cx))
  cat(sprintf("%-12s %5d %7d %20s %8s\n",
              cohorts[[i]]$Cohort[1], nrow(df), sum(df$OS_event),
              sprintf("%.2f (%.2f-%.2f)", hr, ci[1], ci[2]),
              ifelse(pv<0.001,"<0.001",sprintf("%.4f",pv))))
}

library(survival); library(survminer); library(ggplot2); library(gridExtra); library(grid)

cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

prep_cohort <- function(df, label) {
  df <- df[, c("sample_id","RiskScore","RiskGroup","OS_years","OS_event")]
  df$Cohort <- label; df
}

clin10846  <- prep_cohort(cohort_data$GSE10846,  "GSE10846")
clin87371  <- prep_cohort(cohort_data$GSE87371,  "GSE87371")
clin11318  <- prep_cohort(cohort_data$GSE11318,  "GSE11318")
clin181063 <- prep_cohort(cohort_data$GSE181063, "GSE181063")
tcga_raw   <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))
clin_tcga  <- data.frame(
  sample_id=tcga_raw$submitter_id, RiskScore=tcga_raw$RiskScore,
  RiskGroup=tcga_raw$RiskGroup, OS_years=tcga_raw$OS_years,
  OS_event=tcga_raw$OS_event, Cohort="TCGA-DLBC", stringsAsFactors=FALSE
)

cohorts <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
titles  <- c("GSE10846 (Training, n=412)",
             "GSE87371 (Validation, n=221)",
             "GSE11318 (Validation, n=199)",
             "GSE181063 (Validation, n=882)",
             "TCGA-DLBC (Validation, n=45)")

# Custom x-axis limits per cohort (cap at clinically relevant range)
xlims <- c(22, 20, 22, 15, 18)  # years

col_high <- "#E63946"; col_low <- "#457B9D"

make_km_plot <- function(df, title, xmax_cap) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  # Cap OS time at xmax for display
  df_plot <- df
  df_plot$OS_event[df_plot$OS_years > xmax_cap] <- 0
  df_plot$OS_years <- pmin(df_plot$OS_years, xmax_cap)
  
  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df_plot)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)  # use original for stats
  pval <- 1 - pchisq(lr$chisq, df=1)
  plab <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr   <- exp(coef(cx)); ci <- exp(confint(cx))
  annot <- sprintf("%s\nHR = %.2f (%.2f\u2013%.2f)", plab, hr, ci[1], ci[2])
  bby   <- if (xmax_cap > 15) 5L else if (xmax_cap > 8) 4L else 3L

  p <- ggsurvplot(
    fit, data=df_plot,
    palette=c(col_low, col_high),
    legend.labs=c("Low Risk","High Risk"), legend.title="",
    legend=c(0.78, 0.88),
    xlab="Time (years)", ylab="Overall Survival Probability",
    title=title,
    font.title=c(10L,"bold"), font.x=9L, font.y=9L,
    font.tickslab=8L, font.legend=8L,
    risk.table=TRUE, risk.table.height=0.25,
    risk.table.fontsize=3L, risk.table.title="No. at risk",
    tables.theme=theme_cleantable(),
    ggtheme=theme_classic(base_size=9L) +
              theme(plot.title=element_text(hjust=0.5)),
    conf.int=TRUE, conf.int.alpha=0.12,
    pval=FALSE, surv.median.line="hv",
    break.time.by=bby, xlim=c(0, xmax_cap)
  )
  # Place annotation at left side to avoid right-edge clipping
  p$plot <- p$plot +
    annotate("text", x=xmax_cap * 0.05, y=0.12,
             label=annot, hjust=0, vjust=0, size=3L, color="black")
  p
}

plots <- mapply(make_km_plot, cohorts, titles, xlims, SIMPLIFY=FALSE)

# ── Save combined figure ───────────────────────────────────────────────────────
outdir <- translate_path("/mnt/results/DLBCL_prognosis/figures")
combined <- arrange_ggsurvplots(plots, nrow=2, ncol=3, print=FALSE, risk.table.height=0.25)

png(file.path(outdir, "Fig5_MultiCohort_KM_curves.png"),
    width=18, height=12, units="in", res=180)
print(combined)
dev.off()

svg(file.path(outdir, "Fig5_MultiCohort_KM_curves.svg"),
    width=18, height=12)
print(combined)
dev.off()

cat("Saved updated Fig5_MultiCohort_KM_curves\n")

library(survival); library(survminer); library(ggplot2); library(gridExtra); library(grid)

cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

prep_cohort <- function(df, label) {
  df <- df[, c("sample_id","RiskScore","RiskGroup","OS_years","OS_event")]
  df$Cohort <- label; df
}

clin10846  <- prep_cohort(cohort_data$GSE10846,  "GSE10846")
clin87371  <- prep_cohort(cohort_data$GSE87371,  "GSE87371")
clin11318  <- prep_cohort(cohort_data$GSE11318,  "GSE11318")
clin181063 <- prep_cohort(cohort_data$GSE181063, "GSE181063")
tcga_raw   <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))
clin_tcga  <- data.frame(
  sample_id=tcga_raw$submitter_id, RiskScore=tcga_raw$RiskScore,
  RiskGroup=tcga_raw$RiskGroup, OS_years=tcga_raw$OS_years,
  OS_event=tcga_raw$OS_event, Cohort="TCGA-DLBC", stringsAsFactors=FALSE
)

cohorts <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
titles  <- c("GSE10846 (Training, n=412)",
             "GSE87371 (Validation, n=221)",
             "GSE11318 (Validation, n=199)",
             "GSE181063 (Validation, n=882)",
             "TCGA-DLBC (Validation, n=45)")
xlims   <- c(22, 20, 22, 15, 18)

col_high <- "#E63946"; col_low <- "#457B9D"

make_km_plot <- function(df, title, xmax_cap) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  df_plot <- df
  df_plot$OS_event[df_plot$OS_years > xmax_cap] <- 0
  df_plot$OS_years <- pmin(df_plot$OS_years, xmax_cap)

  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df_plot)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval <- 1 - pchisq(lr$chisq, df=1)
  plab <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr   <- exp(coef(cx)); ci <- exp(confint(cx))
  annot <- sprintf("%s\nHR = %.2f (%.2f\u2013%.2f)", plab, hr, ci[1], ci[2])
  bby   <- if (xmax_cap > 15) 5L else if (xmax_cap > 8) 4L else 3L

  p <- ggsurvplot(
    fit, data=df_plot,
    palette=c(col_low, col_high),
    legend.labs=c("Low Risk","High Risk"), legend.title="",
    legend=c(0.78, 0.88),
    xlab="Time (years)", ylab="Overall Survival Probability",
    title=title,
    font.title=c(10L,"bold"), font.x=9L, font.y=9L,
    font.tickslab=8L, font.legend=8L,
    risk.table=TRUE, risk.table.height=0.25,
    risk.table.fontsize=3L, risk.table.title="No. at risk",
    tables.theme=theme_cleantable(),
    ggtheme=theme_classic(base_size=9L) +
              theme(plot.title=element_text(hjust=0.5)),
    conf.int=TRUE, conf.int.alpha=0.12,
    pval=FALSE, surv.median.line="hv",
    break.time.by=bby, xlim=c(0, xmax_cap)
  )
  # Place annotation at 20% from left, bottom area — avoids both edges
  p$plot <- p$plot +
    annotate("text", x=xmax_cap * 0.20, y=0.10,
             label=annot, hjust=0, vjust=0, size=3L, color="black")
  p
}

plots <- mapply(make_km_plot, cohorts, titles, xlims, SIMPLIFY=FALSE)

# ── Manually arrange: 3 top + 2 bottom centered using cowplot ─────────────────
# Extract grobs from each survplot
extract_grob <- function(p) {
  ggsurvplot_combined <- cowplot::plot_grid(
    p$plot, p$table, ncol=1, rel_heights=c(0.75, 0.25)
  )
  ggsurvplot_combined
}

if (!requireNamespace("cowplot", quietly=TRUE)) install.packages("cowplot", quiet=TRUE)
library(cowplot)

grobs <- lapply(plots, extract_grob)

# Top row: 3 panels
top_row    <- plot_grid(grobs[[1]], grobs[[2]], grobs[[3]], nrow=1, ncol=3)
# Bottom row: 2 panels centered (add empty spacer on each side)
spacer     <- ggplot() + theme_void()
bottom_row <- plot_grid(spacer, grobs[[4]], grobs[[5]], spacer,
                        nrow=1, ncol=4, rel_widths=c(0.5, 1, 1, 0.5))

final_fig  <- plot_grid(top_row, bottom_row, nrow=2, ncol=1)

outdir <- translate_path("/mnt/results/DLBCL_prognosis/figures")
png(file.path(outdir, "Fig5_MultiCohort_KM_curves.png"),
    width=18, height=12, units="in", res=180)
print(final_fig)
dev.off()

svg(file.path(outdir, "Fig5_MultiCohort_KM_curves.svg"),
    width=18, height=12)
print(final_fig)
dev.off()

cat("Final Fig5 saved.\n")

library(survival)

cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))
tcga_raw    <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))

prep_cohort <- function(df, label) {
  df <- df[, c("sample_id","RiskScore","RiskGroup","OS_years","OS_event")]
  df$Cohort <- label; df
}

clin10846  <- prep_cohort(cohort_data$GSE10846,  "GSE10846")
clin87371  <- prep_cohort(cohort_data$GSE87371,  "GSE87371")
clin11318  <- prep_cohort(cohort_data$GSE11318,  "GSE11318")
clin181063 <- prep_cohort(cohort_data$GSE181063, "GSE181063")
clin_tcga  <- data.frame(
  sample_id=tcga_raw$submitter_id, RiskScore=tcga_raw$RiskScore,
  RiskGroup=tcga_raw$RiskGroup, OS_years=tcga_raw$OS_years,
  OS_event=tcga_raw$OS_event, Cohort="TCGA-DLBC", stringsAsFactors=FALSE
)

cohorts      <- list(clin10846, clin87371, clin11318, clin181063, clin_tcga)
cohort_names <- c("GSE10846","GSE87371","GSE11318","GSE181063","TCGA-DLBC")
xlims        <- c(22, 20, 22, 15, 18)

# ── 1. Combined risk scores CSV ────────────────────────────────────────────────
all_scores <- do.call(rbind, cohorts)
write.csv(all_scores, 
          translate_path("/mnt/results/DLBCL_prognosis/tables/All_cohorts_risk_scores.csv"),
          row.names=FALSE)
cat("Saved All_cohorts_risk_scores.csv:", nrow(all_scores), "samples\n")

# ── 2. KM statistics summary CSV ──────────────────────────────────────────────
km_stats <- do.call(rbind, lapply(seq_along(cohorts), function(i) {
  df   <- cohorts[[i]]
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  
  # Cap for stats consistency
  df_cap <- df
  df_cap$OS_event[df_cap$OS_years > xlims[i]] <- 0
  df_cap$OS_years <- pmin(df_cap$OS_years, xlims[i])
  
  # Log-rank (on original data)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval <- 1 - pchisq(lr$chisq, df=1)
  
  # Cox HR
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  hr   <- exp(coef(cx))
  ci   <- exp(confint(cx))
  
  # Median OS per group
  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df_cap)
  med  <- summary(fit)$table[, "median"]
  
  data.frame(
    Cohort         = cohort_names[i],
    Type           = ifelse(i==1, "Training", "Validation"),
    N              = nrow(df),
    N_High         = sum(df$RiskGroup=="High"),
    N_Low          = sum(df$RiskGroup=="Low"),
    Events_total   = sum(df$OS_event),
    Events_High    = sum(df$OS_event[df$RiskGroup=="High"]),
    Events_Low     = sum(df$OS_event[df$RiskGroup=="Low"]),
    HR             = round(hr, 3),
    HR_CI_lower    = round(ci[1], 3),
    HR_CI_upper    = round(ci[2], 3),
    LogRank_pvalue = signif(pval, 4),
    Median_OS_High = round(med["RiskGroup=High"], 2),
    Median_OS_Low  = round(med["RiskGroup=Low"], 2),
    stringsAsFactors = FALSE
  )
}))

write.csv(km_stats,
          translate_path("/mnt/results/DLBCL_prognosis/tables/MultiCohort_KM_statistics.csv"),
          row.names=FALSE)
cat("\nSaved MultiCohort_KM_statistics.csv\n")
cat("\nKM Statistics Summary:\n")
print(km_stats[, c("Cohort","Type","N","Events_total","HR","HR_CI_lower","HR_CI_upper","LogRank_pvalue")])

outdir <- translate_path("/mnt/results/DLBCL_prognosis/tables")

# ── 1. TCGA-DLBC risk scores ───────────────────────────────────────────────────
tcga <- read.csv(translate_path("/workspace/TCGA_DLBC_risk_scores.csv"))
write.csv(tcga, file.path(outdir, "TCGA_DLBC_risk_scores.csv"), row.names=FALSE)
cat("Saved TCGA_DLBC_risk_scores.csv:", nrow(tcga), "samples\n")

# ── 2. Individual GEO validation cohort risk scores ────────────────────────────
cohort_data <- readRDS(translate_path("/workspace/cohort_data.rds"))

# GSE87371
write.csv(cohort_data$GSE87371,
          file.path(outdir, "GSE87371_risk_scores.csv"), row.names=FALSE)
cat("Saved GSE87371_risk_scores.csv:", nrow(cohort_data$GSE87371), "samples\n")

# GSE11318
write.csv(cohort_data$GSE11318,
          file.path(outdir, "GSE11318_risk_scores.csv"), row.names=FALSE)
cat("Saved GSE11318_risk_scores.csv:", nrow(cohort_data$GSE11318), "samples\n")

# GSE181063
write.csv(cohort_data$GSE181063,
          file.path(outdir, "GSE181063_risk_scores.csv"), row.names=FALSE)
cat("Saved GSE181063_risk_scores.csv:", nrow(cohort_data$GSE181063), "samples\n")

# ── 3. Verify all tables now saved ────────────────────────────────────────────
cat("\n=== All CSV files in tables/ ===\n")
files <- list.files(outdir, pattern="\\.csv$", full.names=FALSE)
for (f in sort(files)) {
  df <- read.csv(file.path(outdir, f))
  cat(sprintf("  %-45s  %d rows x %d cols\n", f, nrow(df), ncol(df)))
}
