# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
# 注意：本仓库路径含中文，R 在 Windows 下对含中文的**绝对路径**做编码
# 转换会失败（unable to translate to UTF-8），因此统一使用相对路径，
# 并确保在仓库根目录运行。
source(file.path("config", "paths.R"))
# <<< 自动注入结束 >>>

# =============================================================================
# 0. 环境准备
# =============================================================================
# rm(list = ls())   # 迁移工具已注释：会清除注入的路径配置
options(stringsAsFactors = FALSE)

# ---- 如未安装请先取消注释 ----
# install.packages(c("dplyr", "stringr", "ggplot2", "RColorBrewer"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "Biobase", "glmnet", "survival",
#                        "AnnotationDbi", "hgu133plus2.db"))

library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)
library(glmnet)
library(survival)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(hgu133plus2.db)

# =============================================================================
# 1. 可调参数区（只需改这里）
# =============================================================================

gse10846_file <- translate_path("D:/BadiduNetdiskDownload/R/GSE10846survive/GSE10846_series_matrix.txt.gz")
outdir_fig    <- translate_path("D:/bulk-download/DLBCL_prognosis/figures")
outdir_tab    <- translate_path("D:/bulk-download/DLBCL_prognosis/tables")
lasso_seed    <- 42

gene_set <- c(
              "CEP192","CHCHD10","NAGK","PPRC1","BASP1","MTR","TIMM17A","AAK1",
              "MGAT1","USP12","FAM20A","APOL3","GSDMD","ELMO1","PICALM","RASSF4",
              "IFIT3","NAE1","NAP1L4","PPP2R5A","MRGBP","PPP1R18","CTSL")

file_formula_csv <- file.path(outdir_tab, "LASSO_prognostic_formula.csv")
file_risk_csv    <- file.path(outdir_tab, "GSE10846_risk_scores.csv")
file_cv_png      <- file.path(outdir_fig, "Fig1A_LASSO_CV_curve.png")
file_cv_svg      <- file.path(outdir_fig, "Fig1A_LASSO_CV_curve.svg")
file_traj_png    <- file.path(outdir_fig, "Fig1B_LASSO_coef_trajectory.png")
file_traj_svg    <- file.path(outdir_fig, "Fig1B_LASSO_coef_trajectory.svg")

dir.create(outdir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_tab, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. 读取 GSE10846 + 整理临床信息
# =============================================================================

cat("Parsing GSE10846 series matrix...\n")
gse <- getGEO(filename = gse10846_file, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]

cat("Class:", class(gse), "\n")
cat("Dimensions:", nrow(gse), "probes x", ncol(gse), "samples\n")

pdata <- pData(gse)
cat("\nColumn names in pData:\n"); print(colnames(pdata))

strip <- function(x, prefix) str_trim(str_replace(x, fixed(prefix), ""))

clin <- data.frame(
  sample_id       = pdata$geo_accession,
  gender          = strip(pdata$`characteristics_ch1`,    "Gender: "),
  age             = as.numeric(strip(pdata$`characteristics_ch1.1`,  "Age: ")),
  subtype         = strip(pdata$`characteristics_ch1.6`,  "Clinical info: Final microarray diagnosis: "),
  status          = strip(pdata$`characteristics_ch1.7`,  "Clinical info: Follow up status: "),
  follow_up_years = as.numeric(strip(pdata$`characteristics_ch1.8`,  "Clinical info: Follow up years: ")),
  chemo           = strip(pdata$`characteristics_ch1.9`,  "Clinical info: Chemotherapy: "),
  ecog            = as.numeric(strip(pdata$`characteristics_ch1.10`, "Clinical info: ECOG performance status: ")),
  stage           = as.numeric(strip(pdata$`characteristics_ch1.11`, "Clinical info: Stage: ")),
  ldh_ratio       = as.numeric(strip(pdata$`characteristics_ch1.12`, "Clinical info: LDH ratio: ")),
  extranodal      = as.numeric(strip(pdata$`characteristics_ch1.13`, "Clinical info: Number of extranodal sites: ")),
  stringsAsFactors = FALSE
)

clin$OS_event <- ifelse(clin$status == "DEAD", 1,
                        ifelse(clin$status == "ALIVE", 0, NA))
clin$OS_years <- clin$follow_up_years

clin$ipi_age      <- ifelse(!is.na(clin$age)       & clin$age > 60,       1, 0)
clin$ipi_stage    <- ifelse(!is.na(clin$stage)     & clin$stage >= 3,     1, 0)
clin$ipi_ldh      <- ifelse(!is.na(clin$ldh_ratio) & clin$ldh_ratio > 1,  1, 0)
clin$ipi_ecog     <- ifelse(!is.na(clin$ecog)      & clin$ecog >= 2,      1, 0)
clin$ipi_extranod <- ifelse(!is.na(clin$extranodal)& clin$extranodal > 1, 1, 0)
clin$IPI          <- clin$ipi_age + clin$ipi_stage + clin$ipi_ldh +
                     clin$ipi_ecog + clin$ipi_extranod

# =============================================================================
# 3. Probe 注释（统一使用 hgu133plus2.db）
# =============================================================================

expr_mat  <- exprs(gse)
probe_ids <- rownames(expr_mat)

cat("\nMapping probes to gene symbols via hgu133plus2.db...\n")
gene_syms <- mapIds(hgu133plus2.db,
                    keys      = probe_ids,
                    column    = "SYMBOL",
                    keytype   = "PROBEID",
                    multiVals = "first")

probe_gene <- data.frame(
  probe = probe_ids,
  gene  = as.character(gene_syms),
  stringsAsFactors = FALSE
) %>% filter(!is.na(gene) & gene != "")

cat("Probes mapped to symbols:", nrow(probe_gene), "\n")

probe_25 <- probe_gene %>% filter(gene %in% gene_set)
cat("Probes matching 25-gene set:", nrow(probe_25), "\n")
cat("Genes found:\n");     print(sort(unique(probe_25$gene)))
cat("Genes NOT found:\n"); print(setdiff(gene_set, unique(probe_25$gene)))

probe_25$mean_expr <- rowMeans(expr_mat[probe_25$probe, , drop = FALSE])

best_probes <- probe_25 %>%
  group_by(gene) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\nBest probe per gene:\n")
print(best_probes %>% dplyr::select(gene, probe, mean_expr) %>% arrange(gene))

expr_gene <- expr_mat[best_probes$probe, , drop = FALSE]
rownames(expr_gene) <- best_probes$gene
expr_gene <- expr_gene[intersect(gene_set, rownames(expr_gene)), , drop = FALSE]



# =============================================================================
# 5. 全队列重新建模（412 个样本，最终模型）
# =============================================================================

cat("\n================ Full cohort remodeling ================\n")

clin_all <- clin %>%
  filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event))
cat("Full cohort before expression match: n =", nrow(clin_all),
    " events =", sum(clin_all$OS_event), "\n")

expr_all_ids <- intersect(clin_all$sample_id, colnames(expr_gene))
clin_all     <- clin_all[clin_all$sample_id %in% expr_all_ids, , drop = FALSE]
clin_all     <- clin_all[match(expr_all_ids, clin_all$sample_id), , drop = FALSE]
expr_all_t   <- t(expr_gene[, clin_all$sample_id, drop = FALSE])

cat("After matching: n =", nrow(clin_all),
    " events =", sum(clin_all$OS_event), "\n")
cat("EPV (events/25 genes):",
    round(sum(clin_all$OS_event) / ncol(expr_all_t), 1), "\n")

X_all    <- scale(expr_all_t)
surv_all <- Surv(clin_all$OS_years, clin_all$OS_event)

set.seed(lasso_seed)

# Elastic net (alpha=0.5) — 与原代码保持一致
cv_en <- cv.glmnet(
  X_all, surv_all,
  family = "cox", alpha = 0.5, nfolds = 10,
  type.measure = "deviance", standardize = FALSE
)
coef_en <- coef(cv_en, s = "lambda.min")
nz_en   <- coef_en[coef_en[,1] != 0, , drop = FALSE]
cat("\nElastic net (alpha=0.5) lambda.min — non-zero genes (n=",
    nrow(nz_en), "):\n", sep = "")
print(round(nz_en, 5))

# Full cohort LASSO (alpha=1) — 最终模型
cv_lasso_all <- cv.glmnet(
  X_all, surv_all,
  family = "cox", alpha = 1, nfolds = 10,
  type.measure = "deviance", standardize = FALSE
)
coef_lasso_all <- coef(cv_lasso_all, s = "lambda.min")
nz_lasso       <- coef_lasso_all[coef_lasso_all[,1] != 0, , drop = FALSE]
cat("\nLASSO (alpha=1) full cohort lambda.min — non-zero genes (n=",
    nrow(nz_lasso), "):\n", sep = "")
print(round(nz_lasso, 5))

# Lambda 路径参考
lam_seq   <- cv_lasso_all$lambda
nzero_seq <- cv_lasso_all$nzero
cat("\nLambda path — nzero genes (3~10 range):\n")
print(
  data.frame(log_lambda = round(log(lam_seq), 3), nzero = nzero_seq) %>%
    filter(nzero >= 3 & nzero <= 10) %>% head(15)
)

# =============================================================================
# 6. 导出公式 + 风险评分
# =============================================================================

lasso_genes <- rownames(nz_lasso)
lasso_coefs <- as.numeric(nz_lasso[, 1])
names(lasso_coefs) <- lasso_genes

cat("\n=== PROGNOSTIC FORMULA ===\n")
cat("RiskScore =")
for (i in seq_along(lasso_coefs)) {
  sign_str <- ifelse(lasso_coefs[i] >= 0, " + ", " - ")
  cat(sign_str, round(abs(lasso_coefs[i]), 5), "*", lasso_genes[i])
}
cat("\n\n")

formula_df <- data.frame(
  Gene        = lasso_genes,
  Coefficient = lasso_coefs,
  Direction   = ifelse(lasso_coefs > 0, "Risk", "Protective"),
  stringsAsFactors = FALSE
)
write.csv(formula_df, file_formula_csv, row.names = FALSE)

risk_scores <- as.numeric(X_all[, lasso_genes, drop = FALSE] %*% lasso_coefs)
clin_all$RiskScore <- risk_scores
clin_all$RiskGroup <- ifelse(risk_scores >= median(risk_scores), "High", "Low")

cat("Risk score summary:\n"); print(summary(risk_scores))
cat("High risk:", sum(clin_all$RiskGroup == "High"),
    " Low risk:", sum(clin_all$RiskGroup == "Low"), "\n")

write.csv(
  clin_all %>% dplyr::select(sample_id, RiskScore, RiskGroup, OS_years, OS_event,
                      age, gender, stage, IPI, chemo),
  file_risk_csv,
  row.names = FALSE
)
cat("Risk scores saved.\n")

# =============================================================================
# 7. Fig1A：LASSO CV 曲线
# =============================================================================

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

cat("lambda.min: log =", round(lam_min_log, 3), " nzero =", nz_min, "\n")
cat("lambda.1se: log =", round(lam_1se_log, 3), " nzero =", nz_1se, "\n")

# 顶轴刻度：取 nzero 发生变化的转折点
transitions <- c(1, which(diff(cv_df$nzero) != 0) + 1, nrow(cv_df))
transitions <- unique(transitions)
if (length(transitions) > 10) {
  transitions <- transitions[round(seq(1, length(transitions), length.out = 10))]
}
tick_breaks <- cv_df$log_lambda[transitions]
tick_labels <- as.character(cv_df$nzero[transitions])

y_range <- diff(range(cv_df$cvm))
y_min   <- min(cv_df$cvlo) - 0.05 * y_range
y_max   <- max(cv_df$cvup) + 0.18 * y_range
y_ann   <- max(cv_df$cvup) + 0.10 * y_range
x_min   <- min(cv_df$log_lambda) - 0.2
x_max   <- max(cv_df$log_lambda) + 0.2

p_cv <- ggplot(cv_df, aes(x = log_lambda, y = cvm)) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup),
              fill = "#BFCFDF", alpha = 0.45) +
  geom_line(color = "#2166AC", linewidth = 0.9) +
  geom_point(color = "#2166AC", size = 1.4, shape = 16) +
  geom_vline(xintercept = lam_min_log,
             linetype = "dashed", color = "#D73027", linewidth = 0.85) +
  geom_vline(xintercept = lam_1se_log,
             linetype = "dashed", color = "#4DAC26", linewidth = 0.85) +
  annotate("text", x = lam_min_log - 0.06, y = y_ann,
           label = paste0("lambda.min\n(", nz_min, " genes)"),
           color = "#D73027", size = 3.0, hjust = 1.0, fontface = "bold") +
  annotate("text", x = lam_1se_log - 0.06, y = y_ann,
           label = paste0("lambda.1se\n(", nz_1se, " genes)"),
           color = "#4DAC26", size = 3.0, hjust = 1.0, fontface = "bold") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(x_min, x_max),
    sec.axis = sec_axis(
      transform = ~ .,
      name   = "Number of Non-zero Genes",
      breaks = tick_breaks,
      labels = tick_labels
    )
  ) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  labs(
    title    = "LASSO Cox Regression: 10-fold Cross-Validation",
    subtitle = paste0("GSE10846 DLBCL (n=", nrow(clin_all),
                      ", events=", sum(clin_all$OS_event), ")"),
    y        = "Partial Likelihood Deviance (+/- 1 SE)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    axis.title.x.top   = element_text(size = 10, color = "grey35",
                                      margin = margin(b = 5)),
    axis.text.x.top    = element_text(size = 9, color = "grey35"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    plot.margin        = margin(t = 8, r = 20, b = 8, l = 10)
  )

ggsave(file_cv_png, p_cv, width = 8.5, height = 5.5, dpi = 300)
ggsave(file_cv_svg, p_cv, width = 8.5, height = 5.5)
cat("Saved CV curve.\n")

# =============================================================================
# 8. Fig1B：LASSO 系数轨迹图
# =============================================================================

fit_path <- glmnet(X_all, surv_all,
                   family = "cox", alpha = 1, standardize = FALSE)

coef_mat   <- as.matrix(coef(fit_path))
log_lam    <- log(fit_path$lambda)
gene_names <- rownames(coef_mat)

coef_long <- data.frame(
  log_lambda = rep(log_lam, each = nrow(coef_mat)),
  gene       = rep(gene_names, times = ncol(coef_mat)),
  coef       = as.vector(coef_mat)
)

selected_genes <- names(lasso_coefs)
coef_long$selected <- coef_long$gene %in% selected_genes

n_sel      <- length(selected_genes)
pal        <- c(brewer.pal(8, "Set1"), brewer.pal(5, "Dark2"))
sel_colors <- setNames(pal[1:n_sel], selected_genes)

label_df <- coef_long %>%
  filter(selected, coef != 0) %>%
  group_by(gene) %>%
  slice_min(log_lambda, n = 1) %>%
  ungroup()

p_traj <- ggplot() +
  geom_line(
    data = coef_long %>% filter(!selected),
    aes(x = log_lambda, y = coef, group = gene),
    color = "grey75", linewidth = 0.4, alpha = 0.7
  ) +
  geom_line(
    data = coef_long %>% filter(selected),
    aes(x = log_lambda, y = coef, group = gene, color = gene),
    linewidth = 0.85
  ) +
  geom_vline(xintercept = lam_min_log,
             linetype = "dashed", color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = lam_1se_log,
             linetype = "dashed", color = "#4DAC26", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(
    data = label_df,
    aes(x = log_lambda - 0.05, y = coef, label = gene, color = gene),
    size = 2.8, hjust = 1.0, fontface = "bold", show.legend = FALSE
  ) +
  annotate("text",
           x = lam_min_log,
           y = max(abs(coef_mat)) * 0.95,
           label = "lambda.min",
           color = "#D73027", size = 3.0, hjust = 0.5, fontface = "bold") +
  annotate("text",
           x = lam_1se_log,
           y = max(abs(coef_mat)) * 0.95,
           label = "lambda.1se",
           color = "#4DAC26", size = 3.0, hjust = 0.5, fontface = "bold") +
  scale_color_manual(values = sel_colors, name = "Gene") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(min(log_lam) - 1.5, max(log_lam) + 0.2),
    sec.axis = sec_axis(
      transform = ~ .,
      name   = "Number of Non-zero Genes",
      breaks = log(fit_path$lambda[
        round(seq(1, length(fit_path$lambda), length.out = 10))
      ]),
      labels = fit_path$df[
        round(seq(1, length(fit_path$lambda), length.out = 10))
      ]
    )
  ) +
  labs(
    title    = "LASSO Cox Regression: Coefficient Trajectory",
    subtitle = paste0("GSE10846 DLBCL (n=", nrow(clin_all),
                      ") — 25-gene set"),
    y        = "Coefficient"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey40"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    axis.title.x.top   = element_text(size = 10, color = "grey35",
                                      margin = margin(b = 5)),
    axis.text.x.top    = element_text(size = 9, color = "grey35"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    legend.position    = "right",
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 10, face = "bold"),
    plot.margin        = margin(t = 8, r = 10, b = 8, l = 10)
  ) +
  guides(color = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

ggsave(file_traj_png, p_traj, width = 9, height = 5.5, dpi = 300)
ggsave(file_traj_svg, p_traj, width = 9, height = 5.5)
cat("Trajectory plot saved.\n")

cat("\nAll done.\n")
