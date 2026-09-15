# =============================================================================
# 【归档收录说明】本文件取自作者原始脚本 Tcell/yuhou2.R（2394 行）。
#
# 为什么收录：section6 的 07_01~07_05 只**读**下列中间表，仓库里没有产出它们的
# 代码；而本文件正是这些表的**唯一产出者**：
#   - LASSO_prognostic_formula.csv（23 基因 = 07_02 的候选集）  (L80)
#   - GSE10846_risk_scores.csv  (L96)
#   - Cox_univariable_multivariable_results.csv  (L411)
#   - Calibration_data_1_3_5yr.csv  (L528)
#   - MultiCohort_KM_statistics.csv（6 队列，含 NCICCR）  (L1995 起)
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


library(dplyr)
library(survival)
library(glmnet)

# ── 23-gene set (HMGB1 and HAVCR2 removed) ───────────────────────────────────
gene_set_23 <- c("CEP192","CHCHD10","NAGK","PPRC1","BASP1","MTR","TIMM17A","AAK1",
                 "MGAT1","USP12","FAM20A","APOL3","GSDMD","ELMO1","PICALM","RASSF4",
                 "IFIT3","NAE1","NAP1L4","PPP2R5A","MRGBP","PPP1R18","CTSL")
cat("23-gene set:", paste(gene_set_23, collapse=", "), "\n")
cat("n =", length(gene_set_23), "\n\n")

# ── Load cached full expression matrix (22880 genes x 420 samples) ───────────
expr_full <- as.matrix(read.csv(translate_path("/mnt/results/GSE10846_expression_gene_level.csv"),
                                 row.names=1, check.names=FALSE))
cat("Full expression matrix:", nrow(expr_full), "genes x", ncol(expr_full), "samples\n")

# ── Check which 23 genes are present ─────────────────────────────────────────
found    <- intersect(gene_set_23, rownames(expr_full))
missing  <- setdiff(gene_set_23, rownames(expr_full))
cat("Genes found:", length(found), "/", length(gene_set_23), "\n")
cat("Found:", paste(sort(found), collapse=", "), "\n")
if (length(missing) > 0) cat("MISSING:", paste(missing, collapse=", "), "\n")

# ── Load GSE10846 clinical data ───────────────────────────────────────────────
clin_all <- read.csv(translate_path("/mnt/results/DLBCL_prognosis/tables/GSE10846_risk_scores.csv"),
                     stringsAsFactors=FALSE)
cat("\nGSE10846 clinical: n=", nrow(clin_all), " events=", sum(clin_all$OS_event), "\n")

# ── Build 23-gene expression matrix for training samples ─────────────────────
train_ids <- clin_all$sample_id
expr_23   <- t(expr_full[found, train_ids, drop=FALSE])  # 412 x 23
colnames(expr_23) <- found
cat("Training expression matrix:", nrow(expr_23), "x", ncol(expr_23), "\n")

# ── Scale ─────────────────────────────────────────────────────────────────────
X_all <- scale(expr_23)
train_means <- attr(X_all, "scaled:center")
train_sds   <- attr(X_all, "scaled:scale")

surv_all <- Surv(clin_all$OS_years, clin_all$OS_event)
cat("Survival object: n=", nrow(clin_all), " events=", sum(clin_all$OS_event), "\n")
cat("EPV (events/23 genes):", round(sum(clin_all$OS_event)/23, 1), "\n")

set.seed(42)

# ── 10-fold CV LASSO Cox ──────────────────────────────────────────────────────
cv_lasso_all <- cv.glmnet(X_all, surv_all, family="cox", alpha=1, nfolds=10,
                           type.measure="deviance", standardize=FALSE)

cat("lambda.min:", cv_lasso_all$lambda.min, "  log:", round(log(cv_lasso_all$lambda.min),3), "\n")
cat("lambda.1se:", cv_lasso_all$lambda.1se, "  log:", round(log(cv_lasso_all$lambda.1se),3), "\n")

# Coefficients at lambda.min
coef_min  <- coef(cv_lasso_all, s="lambda.min")
nz_lasso  <- coef_min[coef_min[,1] != 0, , drop=FALSE]
cat("\nNon-zero genes at lambda.min (n=", nrow(nz_lasso), "):\n")
print(round(nz_lasso, 6))

lasso_genes <- rownames(nz_lasso)
lasso_coefs <- setNames(as.numeric(nz_lasso[,1]), lasso_genes)

# ── Prognostic formula ────────────────────────────────────────────────────────
cat("\n=== PROGNOSTIC FORMULA (23-gene LASSO) ===\n")
cat("RiskScore =")
for(i in seq_along(lasso_coefs)){
  sign_str <- ifelse(lasso_coefs[i] >= 0, " + ", " - ")
  cat(sign_str, round(abs(lasso_coefs[i]),5), "*", lasso_genes[i])
}
cat("\n")

# ── Save formula ──────────────────────────────────────────────────────────────
formula_df <- data.frame(
  Gene        = lasso_genes,
  Coefficient = lasso_coefs,
  Direction   = ifelse(lasso_coefs > 0, "Risk", "Protective"),
  stringsAsFactors = FALSE
)
write.csv(formula_df,
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/LASSO_prognostic_formula.csv"),
          row.names=FALSE)
cat("\nFormula saved.\n")

# ── Compute risk scores for all 412 patients ──────────────────────────────────
risk_scores <- as.numeric(X_all[, lasso_genes, drop=FALSE] %*% lasso_coefs)
clin_all$RiskScore <- risk_scores
clin_all$RiskGroup <- ifelse(risk_scores >= median(risk_scores), "High", "Low")

cat("Risk score summary:\n"); print(summary(risk_scores))
cat("High risk:", sum(clin_all$RiskGroup=="High"),
    "  Low risk:", sum(clin_all$RiskGroup=="Low"), "\n")

# ── Save risk scores ──────────────────────────────────────────────────────────
write.csv(clin_all %>% select(sample_id, RiskScore, RiskGroup, OS_years, OS_event,
                               age, gender, stage, IPI, chemo),
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/GSE10846_risk_scores.csv"),
          row.names=FALSE)
cat("Risk scores saved.\n")

# ── Verify KM direction ───────────────────────────────────────────────────────
clin_all$RiskGroup <- factor(clin_all$RiskGroup, levels=c("Low","High"))
cx_check <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=clin_all)
hr_check  <- exp(coef(cx_check)); ci_check <- exp(confint(cx_check))
lr_check  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin_all)
pv_check  <- 1 - pchisq(lr_check$chisq, df=1)
cat(sprintf("\nGSE10846 KM check: HR=%.2f (%.2f-%.2f), p=%.4f\n",
            hr_check, ci_check[1], ci_check[2], pv_check))

library(ggplot2)
library(dplyr)

# ── CV data ───────────────────────────────────────────────────────────────────
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

# Top axis: transition points
transitions <- c(1, which(diff(cv_df$nzero) != 0) + 1, nrow(cv_df))
transitions <- unique(transitions)
if(length(transitions) > 10) {
  transitions <- transitions[round(seq(1, length(transitions), length.out=10))]
}
tick_breaks <- cv_df$log_lambda[transitions]
tick_labels <- as.character(cv_df$nzero[transitions])

# y/x limits
y_range <- diff(range(cv_df$cvm))
y_min   <- min(cv_df$cvlo) - 0.05 * y_range
y_max   <- max(cv_df$cvup) + 0.18 * y_range
y_ann   <- max(cv_df$cvup) + 0.10 * y_range
x_min   <- min(cv_df$log_lambda) - 0.2
x_max   <- max(cv_df$log_lambda) + 0.2

p_cv <- ggplot(cv_df, aes(x=log_lambda, y=cvm)) +
  geom_ribbon(aes(ymin=cvlo, ymax=cvup), fill="#BFCFDF", alpha=0.45) +
  geom_line(color="#2166AC", linewidth=0.9) +
  geom_point(color="#2166AC", size=1.4, shape=16) +
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.85) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.85) +
  annotate("text", x=lam_min_log - 0.06, y=y_ann,
           label=paste0("lambda.min\n(", nz_min, " genes)"),
           color="#D73027", size=3.0, hjust=1.0, fontface="bold") +
  annotate("text", x=lam_1se_log - 0.06, y=y_ann,
           label=paste0("lambda.1se\n(", nz_1se, " genes)"),
           color="#4DAC26", size=3.0, hjust=1.0, fontface="bold") +
  scale_x_continuous(
    name   = "log(Lambda)",
    limits = c(x_min, x_max),
    sec.axis = sec_axis(transform=~., name="Number of Non-zero Genes",
                        breaks=tick_breaks, labels=tick_labels)
  ) +
  scale_y_continuous(limits=c(y_min, y_max)) +
  labs(
    title    = "LASSO Cox Regression: 10-fold Cross-Validation",
    subtitle = "GSE10846 DLBCL (n=412, events=163) — 23-gene set",
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

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig1A_LASSO_CV_curve.png"),
       p_cv, width=8.5, height=5.5, dpi=300)
cat("Fig1A saved.\n")

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Full LASSO path
fit_path  <- glmnet(X_all, surv_all, family="cox", alpha=1, standardize=FALSE)
coef_mat  <- as.matrix(coef(fit_path))
log_lam   <- log(fit_path$lambda)
gene_names <- rownames(coef_mat)

coef_long <- data.frame(
  log_lambda = rep(log_lam, each=nrow(coef_mat)),
  gene       = rep(gene_names, times=ncol(coef_mat)),
  coef       = as.vector(coef_mat)
)

selected_genes <- names(lasso_coefs)
coef_long$selected <- coef_long$gene %in% selected_genes

n_sel      <- length(selected_genes)
pal        <- c(brewer.pal(8,"Set1"), brewer.pal(5,"Dark2"))
sel_colors <- setNames(pal[1:n_sel], selected_genes)

lam_min_log <- log(cv_lasso_all$lambda.min)
lam_1se_log <- log(cv_lasso_all$lambda.1se)

# Label positions: leftmost non-zero point per selected gene
label_df <- coef_long %>%
  filter(selected, coef != 0) %>%
  group_by(gene) %>%
  slice_min(log_lambda, n=1) %>%
  ungroup()

p_traj <- ggplot() +
  geom_line(data=coef_long %>% filter(!selected),
            aes(x=log_lambda, y=coef, group=gene),
            color="grey75", linewidth=0.4, alpha=0.7) +
  geom_line(data=coef_long %>% filter(selected),
            aes(x=log_lambda, y=coef, group=gene, color=gene),
            linewidth=0.85) +
  geom_vline(xintercept=lam_min_log, linetype="dashed", color="#D73027", linewidth=0.8) +
  geom_vline(xintercept=lam_1se_log, linetype="dashed", color="#4DAC26", linewidth=0.8) +
  geom_hline(yintercept=0, color="black", linewidth=0.3) +
  geom_text(data=label_df,
            aes(x=log_lambda - 0.05, y=coef, label=gene, color=gene),
            size=2.8, hjust=1.0, fontface="bold", show.legend=FALSE) +
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
    subtitle = "GSE10846 DLBCL (n=412) — 23-gene set",
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

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig1B_LASSO_coef_trajectory.png"),
       p_traj, width=9, height=5.5, dpi=300)
cat("Fig1B saved.\n")

library(survival)
library(dplyr)
library(ggplot2)
library(patchwork)

# ── Prepare Cox dataset ───────────────────────────────────────────────────────
cox_df <- clin_all %>%
  mutate(
    Age60     = ifelse(age >= 60, 1, 0),
    Male      = ifelse(gender == "male", 1, 0),
    Stage_adv = ifelse(stage >= 3, 1, 0),
    IPI_high  = ifelse(IPI >= 3, 1, 0)
  ) %>%
  filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event))

cat("Cox dataset: n=", nrow(cox_df), " events=", sum(cox_df$OS_event), "\n")
cat("Missing:\n")
print(colSums(is.na(cox_df[, c("RiskScore","Age60","Male","Stage_adv","IPI_high")])))

# ── Univariable Cox ───────────────────────────────────────────────────────────
vars <- c("RiskScore","Age60","Male","Stage_adv","IPI_high")
uni_results <- lapply(vars, function(v) {
  df <- cox_df[!is.na(cox_df[[v]]), ]
  m  <- coxph(as.formula(paste0("Surv(OS_years, OS_event) ~ ", v)), data=df)
  s  <- summary(m)
  data.frame(Variable=v, n=nrow(df), events=sum(df$OS_event),
             HR=s$conf.int[1,1], CI_low=s$conf.int[1,3], CI_high=s$conf.int[1,4],
             pval=s$coefficients[1,5], type="Univariable", stringsAsFactors=FALSE)
})
uni_df <- do.call(rbind, uni_results)

# ── Multivariable Cox ─────────────────────────────────────────────────────────
cox_complete <- cox_df %>%
  filter(!is.na(Age60) & !is.na(Male) & !is.na(Stage_adv) & !is.na(IPI_high))
cat("\nComplete cases: n=", nrow(cox_complete), " events=", sum(cox_complete$OS_event), "\n")

multi_fit <- coxph(Surv(OS_years, OS_event) ~ RiskScore + Age60 + Male + Stage_adv + IPI_high,
                   data=cox_complete)
s_multi <- summary(multi_fit)
cat("\nMultivariable Cox:\n")
print(round(s_multi$conf.int, 4))

multi_df <- data.frame(
  Variable = rownames(s_multi$conf.int),
  n=nrow(cox_complete), events=sum(cox_complete$OS_event),
  HR=s_multi$conf.int[,1], CI_low=s_multi$conf.int[,3], CI_high=s_multi$conf.int[,4],
  pval=s_multi$coefficients[,"Pr(>|z|)"], type="Multivariable", stringsAsFactors=FALSE
)

# ── Build forest plot data ────────────────────────────────────────────────────
var_labels <- c(RiskScore="Risk Score", Age60="Age (>=60 vs <60)",
                Male="Sex (Male vs Female)", Stage_adv="Stage (III-IV vs I-II)",
                IPI_high="IPI (High vs Low)")
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

# ── 4-panel patchwork forest plot ────────────────────────────────────────────
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
  labs(x=NULL, y=NULL) + theme_void() +
  theme(plot.margin=margin(t=10, r=0, b=10, l=10))

p_forest <- ggplot(all_df, aes(y=y_pos, x=HR, xmin=CI_low, xmax=CI_high)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.6) +
  geom_errorbar(aes(color=sig), width=0.3, linewidth=0.9, orientation="y") +
  geom_point(aes(color=sig, size=sig), shape=18) +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=1.5, y=15.5, label="Hazard Ratio",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_color_manual(values=sig_colors, guide="none") +
  scale_size_manual(values=c("***"=5,"**"=4.5,"*"=4,"ns"=3.5), guide="none") +
  scale_x_log10(limits=c(0.2, 12), breaks=c(0.25,0.5,1,2,4,8),
                labels=c("0.25","0.5","1","2","4","8"),
                name="Hazard Ratio (log scale)") +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(y=NULL) + theme_classic(base_size=11) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.line.y=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_text(size=9), axis.title.x=element_text(size=10),
        panel.grid.major.x=element_line(color="grey93", linewidth=0.4),
        plot.margin=margin(t=10, r=5, b=10, l=5))

p_hr <- ggplot(all_df, aes(y=y_pos, x=1)) +
  geom_text(aes(label=hr_str), x=0.5, hjust=0.5, size=2.9, color="grey20") +
  geom_hline(yintercept=c(15.0, 8.7), color="grey65", linewidth=0.4) +
  annotate("text", x=0.5, y=15.5, label="HR (95% CI)",
           hjust=0.5, size=3.5, fontface="bold", color="grey25") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(2, 16.5), expand=c(0,0)) +
  labs(x=NULL, y=NULL) + theme_void() +
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
  labs(x=NULL, y=NULL) + theme_void() +
  theme(plot.margin=margin(t=10, r=10, b=10, l=5))

combined <- p_left + p_forest + p_hr + p_pv +
  plot_layout(widths=c(2.2, 2.5, 1.8, 1.0)) +
  plot_annotation(
    title    = "Cox Proportional Hazards Regression",
    subtitle = "GSE10846 DLBCL (n=405, events=159)",
    theme=theme(plot.title=element_text(face="bold", size=13, hjust=0),
                plot.subtitle=element_text(size=10, color="grey40", hjust=0),
                plot.margin=margin(10,10,10,10))
  )

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig2_Cox_forest_plot.png"),
       combined, width=12, height=7, dpi=300)
cat("Fig2 saved.\n")

# Save Cox results
write.csv(all_df %>% select(Variable, label, type=section, HR, CI_low, CI_high, pval, pval_str, hr_str),
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/Cox_univariable_multivariable_results.csv"),
          row.names=FALSE)

library(rms)
library(survival)
library(ggplot2)
library(dplyr)

# ── Nomogram ──────────────────────────────────────────────────────────────────
nom_df <- cox_complete %>%
  select(OS_years, OS_event, RiskScore, IPI) %>%
  filter(!is.na(IPI) & !is.na(RiskScore) & OS_years > 0)

cat("Nomogram dataset: n=", nrow(nom_df), " events=", sum(nom_df$OS_event), "\n")

dd <- datadist(nom_df)
options(datadist="dd")

fit_cph <- cph(Surv(OS_years, OS_event) ~ RiskScore + IPI,
               data=nom_df, x=TRUE, y=TRUE, surv=TRUE, time.inc=1)

# Baseline survival
S0_1 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=1)$surv
S0_3 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=3)$surv
S0_5 <- survest(fit_cph, newdata=data.frame(RiskScore=0, IPI=0), times=5)$surv
center <- fit_cph$center
cat("S0(1):", round(S0_1,3), " S0(3):", round(S0_3,3), " S0(5):", round(S0_5,3), "\n")

f1 <- function(lp) 1 - S0_1^exp(lp - center)
f3 <- function(lp) 1 - S0_3^exp(lp - center)
f5 <- function(lp) 1 - S0_5^exp(lp - center)

nom <- nomogram(fit_cph,
                fun=list(f1, f3, f5),
                funlabel=c("1-Year Mortality","3-Year Mortality","5-Year Mortality"),
                fun.at=list(c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                            c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                            c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)),
                lp=FALSE)

png(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig3A_Nomogram.png"),
    width=3200, height=1800, res=300)
par(mar=c(2, 2, 3, 2), cex=0.95)
plot(nom, xfrac=0.32, label.every=1, col.grid=grey(c(0.8, 0.95)),
     main="Nomogram: 1/3/5-Year Overall Survival Prediction\n(IPI + Risk Score, GSE10846 DLBCL n=405)")
dev.off()
cat("Fig3A saved.\n")

# ── Calibration curves ────────────────────────────────────────────────────────
cox_cal <- coxph(Surv(OS_years, OS_event) ~ RiskScore + IPI, data=nom_df, x=TRUE)
sf_all  <- survfit(cox_cal, newdata=nom_df)

get_pred_at_t <- function(sf, t) {
  s <- summary(sf, times=t, extend=TRUE)
  if(is.matrix(s$surv)) s$surv[1,] else s$surv
}

pred_1 <- get_pred_at_t(sf_all, 1)
pred_3 <- get_pred_at_t(sf_all, 3)
pred_5 <- get_pred_at_t(sf_all, 5)

calibration_decile <- function(pred_surv, time_pt, df, n_groups=10) {
  df2 <- df %>% mutate(pred=pred_surv, grp=ntile(pred, n_groups))
  df2 %>% group_by(grp) %>% do({
    sub <- .
    km  <- survfit(Surv(OS_years, OS_event) ~ 1, data=sub)
    s   <- summary(km, times=time_pt, extend=TRUE)
    obs <- ifelse(length(s$surv)>0, s$surv[1], NA)
    se  <- ifelse(length(s$std.err)>0, s$std.err[1], NA)
    data.frame(mean_pred=mean(sub$pred), obs_surv=obs,
               obs_lower=obs-1.96*se, obs_upper=obs+1.96*se, n=nrow(sub))
  }) %>% ungroup() %>% mutate(time=time_pt)
}

cal_all <- bind_rows(
  calibration_decile(pred_1, 1, nom_df) %>% mutate(label="1-Year OS"),
  calibration_decile(pred_3, 3, nom_df) %>% mutate(label="3-Year OS"),
  calibration_decile(pred_5, 5, nom_df) %>% mutate(label="5-Year OS")
) %>% mutate(label=factor(label, levels=c("1-Year OS","3-Year OS","5-Year OS")))

time_colors <- c("1-Year OS"="#2166AC","3-Year OS"="#4DAC26","5-Year OS"="#D73027")
cal_plot    <- cal_all %>% filter(!is.na(obs_surv))

mae_df <- cal_plot %>%
  group_by(label) %>%
  summarise(MAE=round(mean(abs(mean_pred - obs_surv)), 3), .groups="drop")
cat("MAE per time point:\n"); print(mae_df)

cal_plot <- cal_plot %>% left_join(mae_df, by="label")
ann_df   <- cal_plot %>% group_by(label, MAE) %>% summarise(.groups="drop") %>%
  mutate(x=0.08, y=0.97, txt=paste0("MAE = ", MAE))

p_cal <- ggplot(cal_plot, aes(x=mean_pred, y=obs_surv, color=label)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50", linewidth=0.7) +
  geom_errorbar(aes(ymin=pmax(obs_lower,0), ymax=pmin(obs_upper,1)),
                width=0.015, linewidth=0.6, alpha=0.7) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="loess", se=FALSE, linewidth=1.0, span=1.2) +
  geom_text(data=ann_df, aes(x=x, y=y, label=txt, color=label),
            hjust=0, vjust=1, size=3.2, fontface="bold", inherit.aes=FALSE) +
  facet_wrap(~label, ncol=3) +
  scale_color_manual(values=time_colors, guide="none") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2), name="Predicted Survival Probability") +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), name="Observed Survival Probability (KM)") +
  labs(title="Calibration Curves: Predicted vs Observed Survival",
       subtitle="Decile-based calibration with 95% CI (GSE10846 DLBCL, n=405, B=10 groups)") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9.5, color="grey40"),
        strip.text=element_text(face="bold", size=11),
        strip.background=element_rect(fill="grey95", color="grey70"),
        axis.title=element_text(size=11), axis.text=element_text(size=9),
        panel.grid.major=element_line(color="grey93", linewidth=0.4),
        plot.margin=margin(10,10,10,10))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig3B_Calibration_curves.png"),
       p_cal, width=11, height=4.5, dpi=300)
write.csv(cal_all, translate_path("/mnt/results/DLBCL_prognosis_v2/tables/Calibration_data_1_3_5yr.csv"), row.names=FALSE)
cat("Fig3B saved.\n")

library(GSVA)
library(GSEABase)
library(dplyr)
library(ggplot2)
library(reshape2)
library(msigdbr)

# ── 1. Get immune cell gene sets from MSigDB C7 / curated ────────────────────
# Use LM22 / CIBERSORT-like gene sets from MSigDB
# We'll use curated immune cell signatures from published sources
immune_sigs <- list(
  B_cells       = c("CD19","MS4A1","CD79A","CD79B","BANK1","FCRL5","IGHM","IGHD"),
  CD8_T_cells   = c("CD8A","CD8B","GZMB","PRF1","IFNG","TBX21","EOMES","GZMK"),
  CD4_T_cells   = c("CD4","IL7R","TCF7","LEF1","CCR7","SELL","IL2RA","FOXP3"),
  NK_cells      = c("NCAM1","KLRB1","KLRD1","NKG7","GNLY","FCGR3A","KLRC1","KLRC2"),
  Macrophages   = c("CD68","CD163","MRC1","MSR1","MARCO","VSIG4","TREM2","C1QA"),
  Monocytes     = c("CD14","LYZ","S100A8","S100A9","VCAN","FCN1","SELL","CD300E"),
  Dendritic     = c("ITGAX","CLEC9A","XCR1","SIGLEC6","LILRA4","CLEC4C","IRF7","IRF8"),
  Neutrophils   = c("FCGR3B","CSF3R","CXCR2","S100A12","PGLYRP1","OLFM4","MMP8","MMP9"),
  Mast_cells    = c("KIT","TPSAB1","TPSB2","CPA3","HPGDS","MS4A2","FCER1A","HDC")
)

# ── 2. Subset expression matrix to genes present ──────────────────────────────
# expr_full is 22880 x 420 (genes x samples), rownames = gene symbols
cat("expr_full dim:", dim(expr_full), "\n")
cat("clin_all dim:", dim(clin_all), "\n")

# Keep only samples in clin_all
common_samps <- intersect(colnames(expr_full), clin_all$geo_accession)
cat("Common samples:", length(common_samps), "\n")
expr_sub <- expr_full[, common_samps]

# Check gene coverage
for(ct in names(immune_sigs)) {
  found <- sum(immune_sigs[[ct]] %in% rownames(expr_sub))
  cat(ct, ":", found, "/", length(immune_sigs[[ct]]), "genes found\n")
}

library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

# ── Prepare long format for violin plots ─────────────────────────────────────
long_df <- merged_df %>%
  select(sample_id, RiskGroup, all_of(cell_types)) %>%
  melt(id.vars=c("sample_id","RiskGroup"), variable.name="CellType", value.name="Score") %>%
  mutate(
    CellType  = gsub("_", " ", as.character(CellType)),
    RiskGroup = factor(RiskGroup, levels=c("Low","High"))
  )

# Wilcoxon test per cell type
wilcox_res <- long_df %>%
  group_by(CellType) %>%
  do({
    sub <- .
    lo  <- sub$Score[sub$RiskGroup=="Low"]
    hi  <- sub$Score[sub$RiskGroup=="High"]
    wt  <- wilcox.test(lo, hi)
    data.frame(pval=wt$p.value)
  }) %>%
  ungroup() %>%
  mutate(padj=p.adjust(pval, method="BH"),
         sig=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,"*","ns"))))

cat("Wilcoxon results:\n")
print(wilcox_res %>% arrange(pval))

# ── Violin plot ───────────────────────────────────────────────────────────────
# Compute y positions for significance brackets
y_max_df <- long_df %>%
  group_by(CellType) %>%
  summarise(y_max=max(Score, na.rm=TRUE), .groups="drop")
sig_df <- wilcox_res %>%
  left_join(y_max_df, by="CellType") %>%
  mutate(y_pos=y_max * 1.08)

risk_colors <- c("Low"="#2166AC","High"="#D73027")

p_violin <- ggplot(long_df, aes(x=RiskGroup, y=Score, fill=RiskGroup)) +
  geom_violin(alpha=0.7, trim=TRUE, linewidth=0.4) +
  geom_boxplot(width=0.15, outlier.size=0.5, outlier.alpha=0.4,
               fill="white", linewidth=0.5) +
  geom_text(data=sig_df, aes(x=1.5, y=y_pos, label=sig),
            inherit.aes=FALSE, size=4.5, fontface="bold", color="grey20") +
  facet_wrap(~CellType, nrow=3, scales="free_y") +
  scale_fill_manual(values=risk_colors, name="Risk Group") +
  labs(x=NULL, y="ssGSEA Enrichment Score",
       title="Immune Cell Infiltration by Risk Group",
       subtitle="ssGSEA scores, Wilcoxon test (BH-adjusted)") +
  theme_classic(base_size=11) +
  theme(
    plot.title=element_text(face="bold", size=13),
    plot.subtitle=element_text(size=9.5, color="grey40"),
    strip.text=element_text(face="bold", size=9.5),
    strip.background=element_rect(fill="grey95", color="grey70"),
    axis.text.x=element_text(size=9, face="bold"),
    axis.text.y=element_text(size=8),
    legend.position="bottom",
    panel.grid.major.y=element_line(color="grey93", linewidth=0.4),
    plot.margin=margin(10,10,10,10)
  )

# ── Correlation bar plot ──────────────────────────────────────────────────────
cor_plot_df <- cor_df %>%
  mutate(
    label    = gsub("_", " ", CellType),
    label    = factor(label, levels=label[order(rho)]),
    bar_col  = ifelse(rho > 0, "#D73027", "#2166AC"),
    sig_star = sig
  )

p_cor <- ggplot(cor_plot_df, aes(x=label, y=rho, fill=bar_col)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_text(aes(label=sig_star,
                y=ifelse(rho>=0, rho+0.005, rho-0.005),
                vjust=ifelse(rho>=0, 0, 1)),
            size=4.5, fontface="bold", color="grey20") +
  geom_hline(yintercept=0, linewidth=0.5, color="grey40") +
  scale_fill_identity() +
  scale_y_continuous(limits=c(-0.12, 0.22), breaks=seq(-0.1,0.2,0.05),
                     name="Spearman ρ (vs Risk Score)") +
  labs(x=NULL, title="Spearman Correlation: Risk Score vs Immune Infiltration",
       subtitle="BH-adjusted p-values: * <0.05, ** <0.01, *** <0.001") +
  coord_flip() +
  theme_classic(base_size=11) +
  theme(
    plot.title=element_text(face="bold", size=12),
    plot.subtitle=element_text(size=9.5, color="grey40"),
    axis.text.y=element_text(size=10),
    axis.text.x=element_text(size=9),
    panel.grid.major.x=element_line(color="grey93", linewidth=0.4),
    plot.margin=margin(10,10,10,10)
  )

# ── Combine ───────────────────────────────────────────────────────────────────
fig4 <- p_violin / p_cor + plot_layout(heights=c(2.5, 1)) +
  plot_annotation(tag_levels="A",
    theme=theme(plot.margin=margin(10,10,10,10)))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig4_Immune_infiltration.png"),
       fig4, width=13, height=14, dpi=300)

# Save data
write.csv(cor_df, translate_path("/mnt/results/DLBCL_prognosis_v2/tables/ssGSEA_RiskScore_correlations.csv"), row.names=FALSE)
write.csv(ssgsea_mat %>% as.data.frame() %>% tibble::rownames_to_column("sample_id"),
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/ssGSEA_immune_scores.csv"), row.names=FALSE)
cat("Fig4 saved.\n")

import requests
import pandas as pd
import json

# ── GDC API: fetch NCICCR-DLBCL cases ────────────────────────────────────────
# Project: NCICCR-DLBCL
GDC_CASES = "https://api.gdc.cancer.gov/cases"

payload = {
    "filters": {
        "op": "=",
        "content": {"field": "project.project_id", "value": "NCICCR-DLBCL"}
    },
    "fields": "case_id,submitter_id,demographic.vital_status,demographic.days_to_death,"
              "diagnoses.days_to_last_follow_up,diagnoses.age_at_diagnosis,"
              "diagnoses.tumor_stage",
    "format": "JSON",
    "size": 500
}

resp = requests.post(GDC_CASES, json=payload, timeout=60)
print("Status:", resp.status_code)
data = resp.json()
hits = data["data"]["hits"]
print("Total cases:", data["data"]["pagination"]["total"])
print("Returned:", len(hits))

import pandas as pd
import numpy as np

records = []
for h in hits:
    demo = h.get("demographic", {}) or {}
    diag_list = h.get("diagnoses", []) or []
    diag = diag_list[0] if diag_list else {}
    
    records.append({
        "case_id":           h.get("case_id"),
        "submitter_id":      h.get("submitter_id"),
        "vital_status":      demo.get("vital_status"),
        "days_to_death":     demo.get("days_to_death"),
        "days_to_last_fu":   diag.get("days_to_last_follow_up"),
        "age_at_diagnosis":  diag.get("age_at_diagnosis"),
        "tumor_stage":       diag.get("tumor_stage"),
    })

clin_nci = pd.DataFrame(records)
print("Shape:", clin_nci.shape)
print("\n=== vital_status unique values ===")
print(clin_nci["vital_status"].value_counts(dropna=False))
print("\n=== days_to_death sample ===")
print(clin_nci["days_to_death"].describe())
print("\n=== days_to_last_fu sample ===")
print(clin_nci["days_to_last_fu"].describe())
print("\n=== Missing summary ===")
print(clin_nci.isnull().sum())

# Inspect a Dead case to see where survival time is stored
dead_hits = [h for h in hits if h.get("demographic", {}).get("vital_status") == "Dead"]
print(f"Dead cases: {len(dead_hits)}")
print("\nFirst dead case raw structure:")
import json
print(json.dumps(dead_hits[0], indent=2))

# Build clinical table: only cases with OS time available
records_clean = []
for h in hits:
    demo = h.get("demographic", {}) or {}
    diag_list = h.get("diagnoses", []) or []
    diag = diag_list[0] if diag_list else {}
    
    vs  = demo.get("vital_status")
    fu  = diag.get("days_to_last_follow_up")
    dtd = demo.get("days_to_death")
    
    # OS time: prefer days_to_death if available, else days_to_last_follow_up
    os_days = dtd if dtd is not None else fu
    
    if vs not in ("Dead","Alive") or os_days is None:
        continue
    
    records_clean.append({
        "case_id":          h.get("case_id"),
        "submitter_id":     h.get("submitter_id"),
        "vital_status":     vs,
        "OS_event":         1 if vs == "Dead" else 0,
        "OS_days":          float(os_days),
        "OS_years":         float(os_days) / 365.25,
        "age_at_diagnosis": diag.get("age_at_diagnosis"),  # in days
    })

clin_nci_clean = pd.DataFrame(records_clean)
print("Usable NCICCR-DLBCL cases:", len(clin_nci_clean))
print("Events (Dead):", clin_nci_clean["OS_event"].sum())
print("Censored (Alive):", (clin_nci_clean["OS_event"]==0).sum())
print("\nOS_years stats:")
print(clin_nci_clean["OS_years"].describe().round(2))
print("\nOS_days > 0:", (clin_nci_clean["OS_days"] > 0).sum())
print("OS_days == 0:", (clin_nci_clean["OS_days"] == 0).sum())

import requests, json, time

GDC_FILES = "https://api.gdc.cancer.gov/files"

# Query for RNA-seq gene expression files for NCICCR-DLBCL
payload = {
    "filters": {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "NCICCR-DLBCL"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}}
        ]
    },
    "fields": "file_id,file_name,cases.case_id,cases.submitter_id,data_type,analysis.workflow_type",
    "format": "JSON",
    "size": 600
}

resp = requests.post(GDC_FILES, json=payload, timeout=60)
print("Status:", resp.status_code)
data = resp.json()
print("Total files:", data["data"]["pagination"]["total"])
print("Returned:", len(data["data"]["hits"]))

# Show sample
if data["data"]["hits"]:
    print("\nSample file entry:")
    print(json.dumps(data["data"]["hits"][0], indent=2))

import pandas as pd

# Build file → case mapping
file_hits = data["data"]["hits"]
file_map = []
for f in file_hits:
    cases = f.get("cases", [])
    if cases:
        file_map.append({
            "file_id":      f["file_id"],
            "file_name":    f["file_name"],
            "submitter_id": cases[0]["submitter_id"],
            "case_id":      cases[0]["case_id"]
        })

file_df = pd.DataFrame(file_map)
print("Total files with case mapping:", len(file_df))
print("Unique submitter_ids:", file_df["submitter_id"].nunique())

# Usable cases (have survival time)
usable_ids = set(clin_nci_clean["submitter_id"])
print("Usable cases (have OS time):", len(usable_ids))

# Filter files to usable cases
file_df_usable = file_df[file_df["submitter_id"].isin(usable_ids)].copy()
print("Files for usable cases:", len(file_df_usable))

# Check for duplicates
dup = file_df_usable["submitter_id"].duplicated().sum()
print("Duplicate submitter_ids:", dup)
if dup > 0:
    print(file_df_usable[file_df_usable["submitter_id"].duplicated(keep=False)])

import requests, os, tarfile, io, time
import pandas as pd
import numpy as np

GDC_DOWNLOAD = "https://api.gdc.cancer.gov/data"

file_ids = file_df_usable["file_id"].tolist()
print(f"Downloading {len(file_ids)} files via GDC bulk API...")

# GDC bulk download: POST list of file IDs, returns a tar archive
payload = {"ids": file_ids}
resp = requests.post(
    GDC_DOWNLOAD,
    json=payload,
    headers={"Content-Type": "application/json"},
    stream=True,
    timeout=300
)
print("Status:", resp.status_code)
print("Content-Type:", resp.headers.get("Content-Type",""))

# Save to disk
os.makedirs("/tmp/nciccr_expr", exist_ok=True)
tar_path = "/tmp/nciccr_expr/nciccr_expr.tar.gz"
total = 0
with open(tar_path, "wb") as f:
    for chunk in resp.iter_content(chunk_size=1024*1024):
        if chunk:
            f.write(chunk)
            total += len(chunk)

print(f"Downloaded: {total/1e6:.1f} MB")

import tarfile, os

tar_path = "/tmp/nciccr_expr/nciccr_expr.tar.gz"

with tarfile.open(tar_path, "r:gz") as tar:
    members = tar.getmembers()
    print(f"Total members in archive: {len(members)}")
    # Show first few
    for m in members[:6]:
        print(f"  {m.name}  ({m.size/1e3:.1f} KB)")
    print("  ...")
    for m in members[-3:]:
        print(f"  {m.name}  ({m.size/1e3:.1f} KB)")

import tarfile, io, pandas as pd

tar_path = "/tmp/nciccr_expr/nciccr_expr.tar.gz"

with tarfile.open(tar_path, "r:gz") as tar:
    members = [m for m in tar.getmembers() if m.name.endswith(".tsv")]
    # Read first file
    f = tar.extractfile(members[0])
    df_sample = pd.read_csv(f, sep="\t", comment=None)

print("Shape:", df_sample.shape)
print("\nFirst 8 rows:")
print(df_sample.head(8).to_string())
print("\nColumns:", df_sample.columns.tolist())

import tarfile, io, pandas as pd

tar_path = "/tmp/nciccr_expr/nciccr_expr.tar.gz"

with tarfile.open(tar_path, "r:gz") as tar:
    members = [m for m in tar.getmembers() if m.name.endswith(".tsv")]
    f = tar.extractfile(members[0])
    df_sample = pd.read_csv(f, sep="\t", comment="#", header=0)

print("Shape:", df_sample.shape)
print("\nFirst 8 rows:")
print(df_sample.head(8).to_string())
print("\nColumns:", df_sample.columns.tolist())

import tarfile, io, pandas as pd, numpy as np

# The 11 LASSO model genes
lasso_genes = ["CEP192","CHCHD10","PPRC1","ELMO1","IFIT3","NAE1",
               "BASP1","MTR","USP12","RASSF4","NAP1L4"]

# file_id → submitter_id lookup
fid2sub = dict(zip(file_df_usable["file_id"], file_df_usable["submitter_id"]))

tar_path = "/tmp/nciccr_expr/nciccr_expr.tar.gz"

expr_records = {}
missing_genes_report = {}

with tarfile.open(tar_path, "r:gz") as tar:
    members = [m for m in tar.getmembers() if m.name.endswith(".tsv")]
    print(f"Processing {len(members)} TSV files...")
    
    for i, m in enumerate(members):
        # Extract file_id from path: "<file_id>/<filename>.tsv"
        file_id = m.name.split("/")[0]
        submitter_id = fid2sub.get(file_id)
        if submitter_id is None:
            continue  # not in usable set
        
        f = tar.extractfile(m)
        df = pd.read_csv(f, sep="\t", comment="#", header=0,
                         usecols=["gene_name","tpm_unstranded"])
        # Drop summary rows (N_unmapped etc.)
        df = df.dropna(subset=["gene_name"])
        # Keep model genes only
        df_genes = df[df["gene_name"].isin(lasso_genes)].copy()
        # Handle duplicate gene names: take max TPM
        df_genes = df_genes.groupby("gene_name")["tpm_unstranded"].max().reset_index()
        
        row = {}
        for _, r in df_genes.iterrows():
            row[r["gene_name"]] = np.log2(r["tpm_unstranded"] + 1)
        expr_records[submitter_id] = row
        
        if i % 50 == 0:
            print(f"  {i+1}/{len(members)} done...")

print(f"\nParsed {len(expr_records)} samples")

# Build matrix
expr_nci = pd.DataFrame(expr_records).T  # samples x genes
expr_nci.index.name = "submitter_id"
print("Expression matrix shape:", expr_nci.shape)
print("\nGenes found:", sorted(expr_nci.columns.tolist()))
missing = set(lasso_genes) - set(expr_nci.columns)
print("Missing genes:", missing if missing else "None")
print("\nSample stats (log2 TPM+1):")
print(expr_nci.describe().round(3))

library(dplyr)
library(survival)

# ── Load LASSO formula ────────────────────────────────────────────────────────
formula_df <- read.csv(translate_path("/mnt/results/DLBCL_prognosis_v2/tables/LASSO_prognostic_formula.csv"))
cat("LASSO formula:\n")
print(formula_df)

# ── Load NCICCR expression matrix (written from Python) ──────────────────────
# Pass via Python → R via reticulate or write to CSV first
# We'll write from Python side and read here

# Save expression matrix
expr_nci.reset_index().to_csv("/tmp/nciccr_expr_matrix.csv", index=False)
print("Expression matrix saved:", expr_nci.shape)

# Save clinical data
clin_nci_clean.to_csv("/tmp/nciccr_clinical.csv", index=False)
print("Clinical data saved:", clin_nci_clean.shape)
print("\nClinical columns:", clin_nci_clean.columns.tolist())
print("\nFirst 3 rows:")
print(clin_nci_clean.head(3))

library(dplyr)
library(survival)
library(survminer)

# ── Load data ─────────────────────────────────────────────────────────────────
expr_nci  <- read.csv("/tmp/nciccr_expr_matrix.csv", row.names=1, check.names=FALSE)
clin_nci  <- read.csv("/tmp/nciccr_clinical.csv")

cat("Expression matrix:", dim(expr_nci), "\n")
cat("Clinical data:", dim(clin_nci), "\n")

# ── Apply training-set scaling (train_means / train_sds from GSE10846) ────────
cat("\ntrain_means genes:", names(train_means), "\n")
cat("train_sds genes:",   names(train_sds), "\n")

lasso_genes_ordered <- formula_df$Gene
cat("LASSO genes:", lasso_genes_ordered, "\n")

# Scale each gene using training parameters
expr_scaled <- expr_nci[, lasso_genes_ordered, drop=FALSE]
for (g in lasso_genes_ordered) {
  expr_scaled[[g]] <- (expr_nci[[g]] - train_means[g]) / train_sds[g]
}
cat("\nScaled expression (first 3 rows):\n")
print(round(head(expr_scaled, 3), 3))

# ── Compute risk scores ───────────────────────────────────────────────────────
coefs <- setNames(formula_df$Coefficient, formula_df$Gene)
risk_scores <- as.numeric(as.matrix(expr_scaled[, names(coefs)]) %*% coefs)
cat("\nRisk score stats:\n")
print(summary(risk_scores))

# ── Merge with clinical data ──────────────────────────────────────────────────
# expr_nci rownames = submitter_id
nci_df <- data.frame(
  submitter_id = rownames(expr_nci),
  RiskScore    = risk_scores,
  stringsAsFactors = FALSE
)

# Merge with clinical
nci_merged <- clin_nci %>%
  inner_join(nci_df, by="submitter_id") %>%
  filter(OS_years > 0)

cat("\nMerged NCICCR dataset: n=", nrow(nci_merged),
    " events=", sum(nci_merged$OS_event), "\n")

# ── Optimal cutpoint (median of training set risk scores) ────────────────────
# Use training-set median as cutpoint for consistency
train_median <- median(clin_all$RiskScore)
cat("Training-set median RiskScore:", round(train_median, 4), "\n")

nci_merged <- nci_merged %>%
  mutate(RiskGroup = ifelse(RiskScore >= train_median, "High", "Low"),
         RiskGroup = factor(RiskGroup, levels=c("Low","High")))

cat("Risk group distribution:\n")
print(table(nci_merged$RiskGroup))

# ── KM test ───────────────────────────────────────────────────────────────────
km_fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=nci_merged)
cox_fit <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=nci_merged)
s_cox   <- summary(cox_fit)
HR      <- round(s_cox$conf.int[1,1], 2)
CI_lo   <- round(s_cox$conf.int[1,3], 2)
CI_hi   <- round(s_cox$conf.int[1,4], 2)
pval    <- s_cox$coefficients[1,5]
cat(sprintf("\nNCI-CCR KM: HR=%.2f (%.2f-%.2f), p=%.4f\n", HR, CI_lo, CI_hi, pval))

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(nci_merged %>% select(submitter_id, case_id, vital_status,
                                 OS_days, OS_years, OS_event, RiskScore, RiskGroup),
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/NCICCR_DLBCL_risk_scores.csv"),
          row.names=FALSE)
cat("NCICCR risk scores saved.\n")

library(dplyr)
library(survival)

# ── Within-cohort z-score scaling (standard for cross-platform validation) ────
# Each gene scaled to mean=0, sd=1 within the NCICCR cohort
expr_nci_r <- read.csv("/tmp/nciccr_expr_matrix.csv", row.names=1, check.names=FALSE)
lasso_genes_ordered <- formula_df$Gene

expr_sub <- expr_nci_r[, lasso_genes_ordered, drop=FALSE]

# Z-score within cohort
expr_z <- as.data.frame(scale(expr_sub))  # scale() centers and scales columns
cat("Within-cohort z-score stats (should be ~mean=0, sd=1):\n")
print(round(apply(expr_z, 2, function(x) c(mean=mean(x), sd=sd(x))), 3))

# ── Compute risk scores ───────────────────────────────────────────────────────
coefs <- setNames(formula_df$Coefficient, formula_df$Gene)
risk_scores_z <- as.numeric(as.matrix(expr_z[, names(coefs)]) %*% coefs)
cat("\nRisk score stats (within-cohort z-scored):\n")
print(summary(risk_scores_z))

# ── Merge with clinical ───────────────────────────────────────────────────────
nci_df2 <- data.frame(
  submitter_id = rownames(expr_nci_r),
  RiskScore    = risk_scores_z,
  stringsAsFactors = FALSE
)

nci_merged2 <- clin_nci %>%
  inner_join(nci_df2, by="submitter_id") %>%
  filter(OS_years > 0)

cat("\nMerged: n=", nrow(nci_merged2), " events=", sum(nci_merged2$OS_event), "\n")

# ── Cutpoint: median of THIS cohort's risk scores ─────────────────────────────
nci_median <- median(nci_merged2$RiskScore)
cat("NCICCR median RiskScore:", round(nci_median, 4), "\n")

nci_merged2 <- nci_merged2 %>%
  mutate(RiskGroup = ifelse(RiskScore >= nci_median, "High", "Low"),
         RiskGroup = factor(RiskGroup, levels=c("Low","High")))

cat("Risk group distribution:\n")
print(table(nci_merged2$RiskGroup))

# ── KM + Cox ──────────────────────────────────────────────────────────────────
km_nci   <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=nci_merged2)
cox_nci  <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=nci_merged2)
s_nci    <- summary(cox_nci)
HR_nci   <- round(s_nci$conf.int[1,1], 2)
CI_lo    <- round(s_nci$conf.int[1,3], 2)
CI_hi    <- round(s_nci$conf.int[1,4], 2)
pval_nci <- s_nci$coefficients[1,5]
cat(sprintf("\nNCICCR Cox: HR=%.2f (%.2f-%.2f), p=%.4f\n", HR_nci, CI_lo, CI_hi, pval_nci))

# Log-rank p-value
lr_test <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=nci_merged2)
lr_p    <- 1 - pchisq(lr_test$chisq, df=1)
cat(sprintf("Log-rank p=%.4f\n", lr_p))

# ── Save updated file ─────────────────────────────────────────────────────────
write.csv(nci_merged2 %>% select(submitter_id, case_id, vital_status,
                                  OS_days, OS_years, OS_event, RiskScore, RiskGroup),
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/NCICCR_DLBCL_risk_scores.csv"),
          row.names=FALSE)
cat("NCICCR risk scores (within-cohort z-scored) saved.\n")

library(GEOparse)
library(dplyr)
library(survival)

lasso_genes <- formula_df$Gene
coefs       <- setNames(formula_df$Coefficient, formula_df$Gene)

compute_risk <- function(expr_mat, coefs) {
  # expr_mat: samples x genes data.frame, genes = names(coefs)
  genes_present <- intersect(names(coefs), colnames(expr_mat))
  missing <- setdiff(names(coefs), colnames(expr_mat))
  if (length(missing) > 0) cat("  WARNING: missing genes:", paste(missing, collapse=", "), "\n")
  
  expr_sub <- expr_mat[, genes_present, drop=FALSE]
  # Within-cohort z-score
  expr_z <- as.data.frame(scale(expr_sub))
  # Risk score
  as.numeric(as.matrix(expr_z[, names(coefs)[names(coefs) %in% genes_present]]) %*%
               coefs[names(coefs) %in% genes_present])
}

# ── GSE87371 ──────────────────────────────────────────────────────────────────
cat("=== GSE87371 ===\n")
gse87371 <- GEOparse::getGEO("GSE87371", destdir="/tmp/", GSEMatrix=TRUE, AnnotGPL=FALSE)
eset87   <- gse87371[[1]]

# Expression matrix
expr87   <- exprs(eset87)
cat("Expression dim:", dim(expr87), "\n")

# Platform annotation
gpl87    <- annotation(eset87)
cat("Platform:", gpl87, "\n")

# Get feature data for gene symbol mapping
fdata87  <- fData(eset87)
cat("Feature data columns:", names(fdata87), "\n")
cat("First few rows:\n")
print(head(fdata87[, 1:min(5, ncol(fdata87))]))

library(GEOquery)
library(dplyr)
library(survival)

lasso_genes <- formula_df$Gene
coefs       <- setNames(formula_df$Coefficient, formula_df$Gene)

# Helper: within-cohort z-score → risk score
compute_risk <- function(expr_mat, coefs) {
  # expr_mat: samples x genes
  genes_ok <- intersect(names(coefs), colnames(expr_mat))
  missing  <- setdiff(names(coefs), colnames(expr_mat))
  if (length(missing) > 0) cat("  Missing genes:", paste(missing, collapse=", "), "\n")
  expr_z <- as.data.frame(scale(expr_mat[, genes_ok, drop=FALSE]))
  as.numeric(as.matrix(expr_z) %*% coefs[genes_ok])
}

# Helper: map probe → gene symbol, keep max-mean probe per gene
probe_to_gene <- function(expr_mat, fdata, symbol_col) {
  # expr_mat: probes x samples; fdata has symbol_col
  fdata2 <- fdata[rownames(expr_mat), , drop=FALSE]
  fdata2$symbol <- as.character(fdata2[[symbol_col]])
  fdata2$probe  <- rownames(fdata2)
  # Remove probes with no symbol or multiple mappings
  fdata2 <- fdata2[!is.na(fdata2$symbol) & fdata2$symbol != "" & 
                   !grepl("///", fdata2$symbol), ]
  expr_sub <- expr_mat[fdata2$probe, , drop=FALSE]
  rownames(expr_sub) <- fdata2$symbol
  # For duplicate gene symbols: keep probe with highest mean expression
  means <- rowMeans(expr_sub, na.rm=TRUE)
  keep  <- tapply(seq_along(means), rownames(expr_sub), function(idx) idx[which.max(means[idx])])
  expr_sub[unlist(keep), , drop=FALSE]
}

# ── GSE87371 ──────────────────────────────────────────────────────────────────
cat("=== Fetching GSE87371 ===\n")
gse87 <- getGEO("GSE87371", destdir="/tmp/", GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
eset87 <- gse87[[1]]
cat("Samples:", ncol(eset87), "\n")
cat("Platform:", annotation(eset87), "\n")
cat("Feature data cols:", names(fData(eset87)), "\n")
cat("pData cols:", names(pData(eset87))[1:10], "\n")

library(hgu133plus2.db)
library(dplyr)
library(survival)

# ── GPL570 probe → gene symbol map (reuse from training) ─────────────────────
probe_syms <- AnnotationDbi::select(hgu133plus2.db,
                                     keys=rownames(exprs(eset87)),
                                     columns="SYMBOL", keytype="PROBEID")
probe_syms <- probe_syms[!is.na(probe_syms$SYMBOL) & probe_syms$SYMBOL != "", ]
cat("Probe-symbol mappings:", nrow(probe_syms), "\n")

# Keep only probes mapping to LASSO genes
probe_lasso <- probe_syms[probe_syms$SYMBOL %in% lasso_genes, ]
cat("Probes for LASSO genes:", nrow(probe_lasso), "\n")
print(table(probe_lasso$SYMBOL))

# Expression matrix: probes x samples
expr87_mat <- exprs(eset87)

# For each gene: keep probe with highest mean expression
gene_expr87 <- lapply(lasso_genes, function(g) {
  probes <- probe_lasso$PROBEID[probe_lasso$SYMBOL == g]
  probes <- intersect(probes, rownames(expr87_mat))
  if (length(probes) == 0) return(NULL)
  sub <- expr87_mat[probes, , drop=FALSE]
  best <- probes[which.max(rowMeans(sub))]
  data.frame(t(sub[best, , drop=FALSE]), check.names=FALSE)
})
names(gene_expr87) <- lasso_genes
missing87 <- sapply(gene_expr87, is.null)
cat("Missing genes:", paste(lasso_genes[missing87], collapse=", "), 
    if(!any(missing87)) "None" else "", "\n")

# Combine into samples x genes matrix
expr87_df <- do.call(cbind, gene_expr87[!missing87])
colnames(expr87_df) <- lasso_genes[!missing87]
cat("Expression matrix:", dim(expr87_df), "\n")

# ── Clinical data ─────────────────────────────────────────────────────────────
pd87 <- pData(eset87)
cat("\npData columns:\n")
print(names(pd87))
cat("\nFirst characteristics_ch1 values:\n")
print(head(pd87$characteristics_ch1, 5))
# Show all characteristics columns
char_cols <- grep("characteristics", names(pd87), value=TRUE)
cat("\nAll characteristics columns:", char_cols, "\n")
for (cc in char_cols) {
  cat("\n---", cc, "---\n")
  print(head(pd87[[cc]], 3))
}

BiocManager::install("hgu133plus2.db", ask=FALSE, update=FALSE)
library(hgu133plus2.db)
cat("hgu133plus2.db loaded.\n")

library(hgu133plus2.db)
library(GEOquery)
library(dplyr)
library(survival)

# ── Build GPL570 probe → gene symbol map ──────────────────────────────────────
all_probes <- keys(hgu133plus2.db, keytype="PROBEID")
probe_map  <- AnnotationDbi::select(hgu133plus2.db,
                                     keys=all_probes,
                                     columns="SYMBOL",
                                     keytype="PROBEID")
probe_map  <- probe_map[!is.na(probe_map$SYMBOL) & probe_map$SYMBOL != "", ]
# Keep only probes for LASSO genes
probe_lasso <- probe_map[probe_map$SYMBOL %in% lasso_genes, ]
cat("Probes for LASSO genes:", nrow(probe_lasso), "\n")
print(table(probe_lasso$SYMBOL))

# ── Helper: GPL570 expr matrix → samples x lasso_genes ───────────────────────
gpl570_to_gene_expr <- function(expr_mat) {
  # expr_mat: probes x samples
  result <- lapply(lasso_genes, function(g) {
    probes <- probe_lasso$PROBEID[probe_lasso$SYMBOL == g]
    probes <- intersect(probes, rownames(expr_mat))
    if (length(probes) == 0) return(NULL)
    sub  <- expr_mat[probes, , drop=FALSE]
    best <- probes[which.max(rowMeans(sub, na.rm=TRUE))]
    as.numeric(sub[best, ])
  })
  names(result) <- lasso_genes
  missing <- sapply(result, is.null)
  if (any(missing)) cat("  Missing genes:", paste(lasso_genes[missing], collapse=", "), "\n")
  df <- as.data.frame(do.call(cbind, result[!missing]))
  colnames(df) <- lasso_genes[!missing]
  rownames(df) <- colnames(expr_mat)
  df
}

# ── Helper: within-cohort z-score → risk score ───────────────────────────────
compute_risk <- function(expr_df, coefs) {
  genes_ok <- intersect(names(coefs), colnames(expr_df))
  expr_z   <- as.data.frame(scale(expr_df[, genes_ok, drop=FALSE]))
  as.numeric(as.matrix(expr_z) %*% coefs[genes_ok])
}

# ═══════════════════════════════════════════════════════════════════════════════
# GSE87371
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== GSE87371 ===\n")
pd87   <- pData(eset87)
expr87 <- gpl570_to_gene_expr(exprs(eset87))
cat("Expression:", dim(expr87), "\n")

# Inspect all characteristics columns
char_cols <- grep("characteristics", names(pd87), value=TRUE)
cat("Characteristics columns:", char_cols, "\n")
for (cc in char_cols[1:min(4, length(char_cols))]) {
  cat("\n", cc, "- unique values (first 5):\n")
  print(head(unique(pd87[[cc]]), 5))
}

# Show remaining characteristics columns
for (cc in char_cols[5:length(char_cols)]) {
  cat("\n", cc, "- unique values (first 5):\n")
  print(head(unique(pd87[[cc]]), 5))
}

# ── Parse GSE87371 clinical ───────────────────────────────────────────────────
parse_char <- function(pd, col) {
  # Extract value after ": " from characteristics column
  as.character(gsub("^.*: ", "", pd[[col]]))
}

clin87 <- data.frame(
  sample_id  = rownames(pd87),
  os_time    = as.numeric(parse_char(pd87, "characteristics_ch1.10")),  # months
  cens_os    = as.numeric(parse_char(pd87, "characteristics_ch1.11")),
  stringsAsFactors = FALSE
) %>%
  mutate(
    OS_event = 1 - cens_os,   # flip: cens_os=1 means censored (alive)
    OS_years = os_time / 12
  ) %>%
  filter(!is.na(os_time) & !is.na(cens_os) & os_time > 0)

cat("GSE87371 clinical: n=", nrow(clin87), " events=", sum(clin87$OS_event), "\n")
cat("cens_os unique values:", unique(clin87$cens_os), "\n")
cat("OS_event distribution:", table(clin87$OS_event), "\n")

# ── Risk scores ───────────────────────────────────────────────────────────────
expr87_sub <- expr87[clin87$sample_id, , drop=FALSE]
clin87$RiskScore <- compute_risk(expr87_sub, coefs)
clin87$RiskGroup <- ifelse(clin87$RiskScore >= median(clin87$RiskScore), "High", "Low")
clin87$RiskGroup <- factor(clin87$RiskGroup, levels=c("Low","High"))
clin87$cohort    <- "GSE87371"

# ── KM + Cox ──────────────────────────────────────────────────────────────────
cox87 <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=clin87)
s87   <- summary(cox87)
lr87  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin87)
cat(sprintf("GSE87371: HR=%.2f (%.2f-%.2f), Cox p=%.4f, LR p=%.4f\n",
    s87$conf.int[1,1], s87$conf.int[1,3], s87$conf.int[1,4],
    s87$coefficients[1,5],
    1 - pchisq(lr87$chisq, df=1)))

# ═══════════════════════════════════════════════════════════════════════════════
# GSE11318
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== Fetching GSE11318 ===\n")
gse113 <- getGEO("GSE11318", destdir="/tmp/", GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
cat("Number of esets:", length(gse113), "\n")
eset113 <- gse113[[1]]
cat("Samples:", ncol(eset113), "  Platform:", annotation(eset113), "\n")

pd113 <- pData(eset113)
char_cols113 <- grep("characteristics", names(pd113), value=TRUE)
cat("Characteristics columns:", length(char_cols113), "\n")
for (cc in char_cols113) {
  cat("\n", cc, ":\n")
  print(head(unique(pd113[[cc]]), 4))
}

# ── GSE11318 clinical ─────────────────────────────────────────────────────────
expr113 <- gpl570_to_gene_expr(exprs(eset113))
cat("GSE11318 expression:", dim(expr113), "\n")

clin113 <- data.frame(
  sample_id  = rownames(pd113),
  status_raw = parse_char(pd113, "characteristics_ch1.7"),
  fu_years   = as.numeric(parse_char(pd113, "characteristics_ch1.8")),
  stringsAsFactors = FALSE
) %>%
  mutate(
    OS_event = ifelse(status_raw == "DEAD", 1, ifelse(status_raw == "ALIVE", 0, NA)),
    OS_years = fu_years
  ) %>%
  filter(!is.na(OS_event) & !is.na(OS_years) & OS_years > 0)

cat("GSE11318: n=", nrow(clin113), " events=", sum(clin113$OS_event), "\n")

expr113_sub      <- expr113[clin113$sample_id, , drop=FALSE]
clin113$RiskScore <- compute_risk(expr113_sub, coefs)
clin113$RiskGroup <- factor(ifelse(clin113$RiskScore >= median(clin113$RiskScore), "High", "Low"),
                             levels=c("Low","High"))
clin113$cohort   <- "GSE11318"

cox113 <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=clin113)
s113   <- summary(cox113)
lr113  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin113)
cat(sprintf("GSE11318: HR=%.2f (%.2f-%.2f), Cox p=%.4f, LR p=%.4f\n",
    s113$conf.int[1,1], s113$conf.int[1,3], s113$conf.int[1,4],
    s113$coefficients[1,5], 1 - pchisq(lr113$chisq, df=1)))

# ═══════════════════════════════════════════════════════════════════════════════
# GSE181063
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Fetching GSE181063 ===\n")
gse181 <- getGEO("GSE181063", destdir="/tmp/", GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
cat("Number of esets:", length(gse181), "\n")
eset181 <- gse181[[1]]
cat("Samples:", ncol(eset181), "  Platform:", annotation(eset181), "\n")

pd181 <- pData(eset181)
char_cols181 <- grep("characteristics", names(pd181), value=TRUE)
cat("Characteristics columns:", length(char_cols181), "\n")
for (cc in char_cols181) {
  cat("\n", cc, ":\n")
  print(head(unique(pd181[[cc]]), 4))
}

# ── GSE181063: GPL14951 (Illumina HT-12 v4) ──────────────────────────────────
# Feature data may have gene symbols embedded
fdata181 <- fData(eset181)
cat("Feature data columns:", names(fdata181), "\n")
cat("Rows:", nrow(fdata181), "\n")
if (ncol(fdata181) > 0) print(head(fdata181[, 1:min(4, ncol(fdata181))]))

BiocManager::install("illuminaHumanv4.db", ask=FALSE, update=FALSE)
library(illuminaHumanv4.db)
cat("illuminaHumanv4.db loaded.\n")

# Build probe → gene symbol map for LASSO genes
ilmn_probes <- keys(illuminaHumanv4.db, keytype="PROBEID")
ilmn_map <- AnnotationDbi::select(illuminaHumanv4.db,
                                   keys=ilmn_probes,
                                   columns="SYMBOL",
                                   keytype="PROBEID")
ilmn_map <- ilmn_map[!is.na(ilmn_map$SYMBOL) & ilmn_map$SYMBOL != "", ]

# Keep only probes for LASSO genes
ilmn_lasso <- ilmn_map[ilmn_map$SYMBOL %in% lasso_genes, ]
cat("Illumina probes for LASSO genes:", nrow(ilmn_lasso), "\n")
print(table(ilmn_lasso$SYMBOL))

# Check overlap with expression matrix probes
expr181_mat <- exprs(eset181)
cat("Expression matrix probes:", nrow(expr181_mat), "\n")
ilmn_lasso_present <- ilmn_lasso[ilmn_lasso$PROBEID %in% rownames(expr181_mat), ]
cat("Probes present in matrix:", nrow(ilmn_lasso_present), "\n")
print(table(ilmn_lasso_present$SYMBOL))

# ── Build Illumina gene expression matrix ─────────────────────────────────────
ilmn_to_gene_expr <- function(expr_mat, probe_lasso_df) {
  result <- lapply(lasso_genes, function(g) {
    probes <- probe_lasso_df$PROBEID[probe_lasso_df$SYMBOL == g]
    probes <- intersect(probes, rownames(expr_mat))
    if (length(probes) == 0) return(NULL)
    sub  <- expr_mat[probes, , drop=FALSE]
    best <- probes[which.max(rowMeans(sub, na.rm=TRUE))]
    as.numeric(sub[best, ])
  })
  names(result) <- lasso_genes
  missing <- sapply(result, is.null)
  if (any(missing)) cat("  Missing:", paste(lasso_genes[missing], collapse=", "), "\n")
  df <- as.data.frame(do.call(cbind, result[!missing]))
  colnames(df) <- lasso_genes[!missing]
  rownames(df) <- colnames(expr_mat)
  df
}

expr181 <- ilmn_to_gene_expr(expr181_mat, ilmn_lasso_present)
cat("GSE181063 expression:", dim(expr181), "\n")

# ── Filter to DLBCL samples only ─────────────────────────────────────────────
diag181 <- parse_char(pd181, "characteristics_ch1")   # diagnostic_group
cat("Diagnostic groups:\n"); print(table(diag181))

dlbcl_idx <- diag181 == "DLBCL"
cat("DLBCL samples:", sum(dlbcl_idx), "\n")

# ── Clinical data ─────────────────────────────────────────────────────────────
clin181 <- data.frame(
  sample_id = rownames(pd181),
  diag      = diag181,
  os_status = as.numeric(parse_char(pd181, "characteristics_ch1.14")),
  os_years  = as.numeric(parse_char(pd181, "characteristics_ch1.15")),
  stringsAsFactors = FALSE
) %>%
  filter(diag == "DLBCL" & !is.na(os_status) & !is.na(os_years) & os_years > 0)

cat("GSE181063 DLBCL with OS: n=", nrow(clin181), " events=", sum(clin181$os_status), "\n")

# ── Risk scores ───────────────────────────────────────────────────────────────
expr181_sub       <- expr181[clin181$sample_id, , drop=FALSE]
clin181$RiskScore <- compute_risk(expr181_sub, coefs)
clin181$OS_event  <- clin181$os_status
clin181$OS_years  <- clin181$os_years
clin181$RiskGroup <- factor(ifelse(clin181$RiskScore >= median(clin181$RiskScore),
                                    "High", "Low"), levels=c("Low","High"))
clin181$cohort    <- "GSE181063"

cox181 <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=clin181)
s181   <- summary(cox181)
lr181  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin181)
cat(sprintf("GSE181063: HR=%.2f (%.2f-%.2f), Cox p=%.4f, LR p=%.4f\n",
    s181$conf.int[1,1], s181$conf.int[1,3], s181$conf.int[1,4],
    s181$coefficients[1,5], 1 - pchisq(lr181$chisq, df=1)))

# Try TCGAbiolinks first
if (!requireNamespace("TCGAbiolinks", quietly=TRUE)) {
  BiocManager::install("TCGAbiolinks", ask=FALSE, update=FALSE)
}
library(TCGAbiolinks)
cat("TCGAbiolinks loaded.\n")

# Query TCGA-DLBC RNA-seq (STAR counts)
query_dlbc <- GDCquery(
  project      = "TCGA-DLBC",
  data.category = "Transcriptome Profiling",
  data.type    = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
cat("Files found:", nrow(getResults(query_dlbc)), "\n")

# Download expression data
GDCdownload(query_dlbc, directory="/tmp/TCGA_DLBC/", method="api")
cat("Download complete.\n")

# Prepare SummarizedExperiment
se_dlbc <- GDCprepare(query_dlbc, directory="/tmp/TCGA_DLBC/")
cat("SE dimensions:", dim(se_dlbc), "\n")
cat("Assay names:", assayNames(se_dlbc), "\n")

cat("SE class:", class(se_dlbc), "\n")
cat("SE dimensions:", dim(se_dlbc), "\n")
cat("Assay names:", assayNames(se_dlbc), "\n")
cat("rowData names:", names(rowData(se_dlbc)), "\n")
cat("colData names:", names(colData(se_dlbc)), "\n")

# Check gene names
cat("\nFirst 5 rowData gene symbols:\n")
print(head(rowData(se_dlbc)$gene_name, 5))

# Check clinical columns relevant to OS
cd <- as.data.frame(colData(se_dlbc))
os_cols <- grep("vital|days_to|os|survival|death|follow", names(cd), 
                value=TRUE, ignore.case=TRUE)
cat("\nOS-related columns:", os_cols, "\n")
for (cc in os_cols[1:min(6, length(os_cols))]) {
  cat("\n", cc, ":\n")
  print(head(unique(cd[[cc]]), 4))
}

cat("se_dlbc exists:", exists("se_dlbc"), "\n")
if (exists("se_dlbc")) {
  cat("Class:", class(se_dlbc), "\n")
  cat("Dim:", dim(se_dlbc), "\n")
  cat("Assays:", assayNames(se_dlbc), "\n")
  cat("rowData cols:", names(rowData(se_dlbc)), "\n")
  cat("colData cols (first 20):", names(colData(se_dlbc))[1:20], "\n")
  
  # Check gene names
  cat("\nFirst 5 gene names:\n")
  print(head(rowData(se_dlbc)$gene_name, 5))
  
  # OS-related columns
  cd <- as.data.frame(colData(se_dlbc))
  os_cols <- grep("vital|days_to|os_|survival|death|follow", names(cd), 
                  value=TRUE, ignore.case=TRUE)
  cat("\nOS-related columns:", os_cols, "\n")
  for (cc in os_cols[1:min(6, length(os_cols))]) {
    cat("\n", cc, ":\n")
    print(head(unique(cd[[cc]]), 4))
  }
} else {
  cat("se_dlbc not found — need to re-run GDCprepare\n")
}

library(SummarizedExperiment)

cat("Dim:", dim(se_dlbc), "\n")
cat("Assays:", assayNames(se_dlbc), "\n")
cat("rowData cols:", names(rowData(se_dlbc)), "\n")

# Gene names
gene_names <- rowData(se_dlbc)$gene_name
cat("First 5 genes:", head(gene_names, 5), "\n")
cat("LASSO genes present:", sum(lasso_genes %in% gene_names), "/", length(lasso_genes), "\n")

# Clinical data
cd <- as.data.frame(colData(se_dlbc))
os_cols <- grep("vital|days_to|os_|survival|death|follow", names(cd), 
                value=TRUE, ignore.case=TRUE)
cat("\nOS-related columns:", os_cols, "\n")
for (cc in os_cols[1:min(8, length(os_cols))]) {
  cat("\n", cc, ":\n")
  print(head(unique(cd[[cc]]), 4))
}

library(SummarizedExperiment)
library(dplyr)
library(survival)

# ── Extract log2(TPM+1) for LASSO genes ──────────────────────────────────────
tpm_mat   <- assay(se_dlbc, "tpm_unstrand")   # genes x samples
gene_names <- rowData(se_dlbc)$gene_name

# For each LASSO gene: pick row with highest mean TPM (handles duplicates)
expr_tcga <- lapply(lasso_genes, function(g) {
  idx <- which(gene_names == g)
  if (length(idx) == 0) return(NULL)
  sub <- tpm_mat[idx, , drop=FALSE]
  best <- idx[which.max(rowMeans(sub, na.rm=TRUE))]
  log2(as.numeric(tpm_mat[best, ]) + 1)
})
names(expr_tcga) <- lasso_genes
missing_tcga <- sapply(expr_tcga, is.null)
cat("Missing genes:", if(any(missing_tcga)) paste(lasso_genes[missing_tcga], collapse=", ") else "None", "\n")

expr_tcga_df <- as.data.frame(do.call(cbind, expr_tcga[!missing_tcga]))
colnames(expr_tcga_df) <- lasso_genes[!missing_tcga]
rownames(expr_tcga_df) <- colnames(se_dlbc)
cat("TCGA expression matrix:", dim(expr_tcga_df), "\n")

# ── Clinical data ─────────────────────────────────────────────────────────────
cd <- as.data.frame(colData(se_dlbc))
clin_tcga <- data.frame(
  sample_id    = rownames(cd),
  vital_status = cd$vital_status,
  days_to_death = as.numeric(cd$days_to_death),
  days_to_last_fu = as.numeric(cd$days_to_last_follow_up),
  stringsAsFactors = FALSE
) %>%
  mutate(
    OS_event = ifelse(vital_status == "Dead", 1, ifelse(vital_status == "Alive", 0, NA)),
    OS_days  = ifelse(!is.na(days_to_death), days_to_death, days_to_last_fu),
    OS_years = OS_days / 365.25
  ) %>%
  filter(!is.na(OS_event) & !is.na(OS_days) & OS_days > 0)

cat("TCGA-DLBC clinical: n=", nrow(clin_tcga), " events=", sum(clin_tcga$OS_event), "\n")
cat("vital_status values:", table(clin_tcga$vital_status), "\n")

# ── Risk scores ───────────────────────────────────────────────────────────────
expr_tcga_sub      <- expr_tcga_df[clin_tcga$sample_id, , drop=FALSE]
clin_tcga$RiskScore <- compute_risk(expr_tcga_sub, coefs)
clin_tcga$RiskGroup <- factor(
  ifelse(clin_tcga$RiskScore >= median(clin_tcga$RiskScore), "High", "Low"),
  levels=c("Low","High"))
clin_tcga$cohort <- "TCGA-DLBC"

# ── KM + Cox ──────────────────────────────────────────────────────────────────
cox_tcga <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=clin_tcga)
s_tcga   <- summary(cox_tcga)
lr_tcga  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=clin_tcga)
cat(sprintf("TCGA-DLBC: HR=%.2f (%.2f-%.2f), Cox p=%.4f, LR p=%.4f\n",
    s_tcga$conf.int[1,1], s_tcga$conf.int[1,3], s_tcga$conf.int[1,4],
    s_tcga$coefficients[1,5], 1 - pchisq(lr_tcga$chisq, df=1)))

library(survival)
library(survminer)
library(cowplot)
library(dplyr)
library(ggplot2)

# ── Prepare all 6 cohort datasets ────────────────────────────────────────────
# 1. GSE10846 (training)
df1 <- clin_all %>%
  select(sample_id, OS_years, OS_event, RiskGroup) %>%
  filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event)) %>%
  mutate(cohort="GSE10846\n(Training, n=412)")

# 2. GSE87371
df2 <- clin87 %>%
  select(sample_id, OS_years, OS_event, RiskGroup) %>%
  mutate(cohort="GSE87371\n(Validation, n=221)")

# 3. GSE11318
df3 <- clin113 %>%
  select(sample_id, OS_years, OS_event, RiskGroup) %>%
  mutate(cohort="GSE11318\n(Validation, n=199)")

# 4. GSE181063
df4 <- clin181 %>%
  select(sample_id, OS_years, OS_event, RiskGroup) %>%
  mutate(cohort="GSE181063\n(Validation, n=1144)")

# 5. NCICCR-DLBCL
nci_data <- read.csv(translate_path("/mnt/results/DLBCL_prognosis_v2/tables/NCICCR_DLBCL_risk_scores.csv"))
df5 <- nci_data %>%
  select(sample_id=submitter_id, OS_years, OS_event, RiskGroup) %>%
  mutate(RiskGroup=factor(RiskGroup, levels=c("Low","High")),
         cohort="NCICCR-DLBCL\n(Validation, n=234)")

# 6. TCGA-DLBC
df6 <- clin_tcga %>%
  select(sample_id, OS_years, OS_event, RiskGroup) %>%
  mutate(cohort="TCGA-DLBC\n(Validation, n=45)")

# ── Helper: draw one KM panel ─────────────────────────────────────────────────
risk_pal <- c("Low"="#2166AC", "High"="#D73027")

make_km <- function(df, title_str, show_table=TRUE) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  
  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  cox  <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc   <- summary(cox)
  HR   <- round(sc$conf.int[1,1], 2)
  CIlo <- round(sc$conf.int[1,3], 2)
  CIhi <- round(sc$conf.int[1,4], 2)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval <- 1 - pchisq(lr$chisq, df=1)
  pstr <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  
  subtitle_str <- sprintf("HR = %.2f (%.2f–%.2f), %s", HR, CIlo, CIhi, pstr)
  
  p <- ggsurvplot(
    fit, data=df,
    palette        = risk_pal,
    legend.labs    = c("Low Risk","High Risk"),
    legend.title   = "",
    legend         = c(0.75, 0.85),
    xlab           = "Time (years)",
    ylab           = "Overall Survival",
    title          = title_str,
    subtitle       = subtitle_str,
    font.title     = c(11, "bold"),
    font.subtitle  = c(8.5, "plain", "grey30"),
    font.x         = c(10, "plain"),
    font.y         = c(10, "plain"),
    font.tickslab  = c(9, "plain"),
    font.legend    = c(9, "plain"),
    risk.table     = show_table,
    risk.table.height = 0.22,
    risk.table.fontsize = 3.0,
    risk.table.title = "No. at risk",
    tables.theme   = theme_cleantable(),
    conf.int       = TRUE,
    conf.int.alpha = 0.12,
    linetype       = "solid",
    size           = 0.8,
    ggtheme        = theme_classic(base_size=10) +
      theme(plot.title    = element_text(face="bold", size=11, hjust=0),
            plot.subtitle = element_text(size=8.5, color="grey30", hjust=0),
            legend.background = element_rect(fill=alpha("white",0.7), color=NA),
            panel.grid.major.y = element_line(color="grey93", linewidth=0.3))
  )
  p
}

# ── Draw all 6 panels ─────────────────────────────────────────────────────────
p1 <- make_km(df1, "GSE10846 (Training)")
p2 <- make_km(df2, "GSE87371")
p3 <- make_km(df3, "GSE11318")
p4 <- make_km(df4, "GSE181063")
p5 <- make_km(df5, "NCICCR-DLBCL")
p6 <- make_km(df6, "TCGA-DLBC")

cat("All 6 KM plots created.\n")

# ── Combine with cowplot ──────────────────────────────────────────────────────
plots_list <- list(p1$plot, p2$plot, p3$plot, p4$plot, p5$plot, p6$plot)
# Add risk tables below each
combined <- plot_grid(
  plot_grid(p1$plot, p1$table, ncol=1, rel_heights=c(3,0.8)),
  plot_grid(p2$plot, p2$table, ncol=1, rel_heights=c(3,0.8)),
  plot_grid(p3$plot, p3$table, ncol=1, rel_heights=c(3,0.8)),
  plot_grid(p4$plot, p4$table, ncol=1, rel_heights=c(3,0.8)),
  plot_grid(p5$plot, p5$table, ncol=1, rel_heights=c(3,0.8)),
  plot_grid(p6$plot, p6$table, ncol=1, rel_heights=c(3,0.8)),
  ncol=3, nrow=2,
  labels=c("A","B","C","D","E","F"),
  label_size=14, label_fontface="bold"
)

title_grob <- ggdraw() +
  draw_label("11-Gene Risk Score: Kaplan-Meier Survival Analysis Across 6 DLBCL Cohorts",
             fontface="bold", size=13, x=0.5, hjust=0.5)

final_fig5 <- plot_grid(title_grob, combined, ncol=1, rel_heights=c(0.04, 1))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.png"),
       final_fig5, width=21, height=14, dpi=300)
ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.svg"),
       final_fig5, width=21, height=14)
cat("Fig5 saved.\n")

library(survival)
library(survminer)
library(cowplot)
library(ggplot2)

# Explicitly use dplyr::select to avoid AnnotationDbi masking
sel <- dplyr::select

# ── Prepare all 6 cohort datasets ────────────────────────────────────────────
df1 <- clin_all %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event)) %>%
  dplyr::mutate(cohort = "GSE10846")

df2 <- clin87 %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::mutate(cohort = "GSE87371")

df3 <- clin113 %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::mutate(cohort = "GSE11318")

df4 <- clin181 %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::mutate(cohort = "GSE181063")

nci_data <- read.csv(translate_path("/mnt/results/DLBCL_prognosis_v2/tables/NCICCR_DLBCL_risk_scores.csv"))
df5 <- nci_data %>%
  dplyr::rename(sample_id = submitter_id) %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::mutate(RiskGroup = factor(RiskGroup, levels=c("Low","High")),
                cohort = "NCICCR-DLBCL")

df6 <- clin_tcga %>%
  sel(sample_id, OS_years, OS_event, RiskGroup) %>%
  dplyr::mutate(cohort = "TCGA-DLBC")

cat("Dataset sizes:\n")
for (df in list(df1,df2,df3,df4,df5,df6)) {
  cat(" n=", nrow(df), " events=", sum(df$OS_event), "\n")
}

# ── Helper: one KM panel ──────────────────────────────────────────────────────
risk_pal <- c("Low"="#2166AC","High"="#D73027")

make_km <- function(df, title_str) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  cox  <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc   <- summary(cox)
  HR   <- round(sc$conf.int[1,1], 2)
  CIlo <- round(sc$conf.int[1,3], 2)
  CIhi <- round(sc$conf.int[1,4], 2)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval <- 1 - pchisq(lr$chisq, df=1)
  pstr <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  subtitle_str <- sprintf("HR = %.2f (%.2f\u2013%.2f), %s", HR, CIlo, CIhi, pstr)
  n_lo <- sum(df$RiskGroup=="Low");  n_hi <- sum(df$RiskGroup=="High")
  
  ggsurvplot(
    fit, data=df,
    palette        = risk_pal,
    legend.labs    = c(sprintf("Low (n=%d)", n_lo), sprintf("High (n=%d)", n_hi)),
    legend.title   = "",
    legend         = c(0.78, 0.88),
    xlab           = "Time (years)",
    ylab           = "Overall Survival Probability",
    title          = title_str,
    subtitle       = subtitle_str,
    font.title     = c(11, "bold"),
    font.subtitle  = c(8.5, "plain", "grey30"),
    font.x         = c(9.5, "plain"),
    font.y         = c(9.5, "plain"),
    font.tickslab  = c(8.5, "plain"),
    font.legend    = c(8.5, "plain"),
    risk.table     = TRUE,
    risk.table.height = 0.20,
    risk.table.fontsize = 2.8,
    risk.table.title = "No. at risk",
    tables.theme   = theme_cleantable(),
    conf.int       = TRUE,
    conf.int.alpha = 0.12,
    size           = 0.85,
    ggtheme        = theme_classic(base_size=10) +
      theme(plot.title    = element_text(face="bold", size=11, hjust=0),
            plot.subtitle = element_text(size=8.5, color="grey35", hjust=0),
            legend.background = element_rect(fill=alpha("white",0.75), color=NA),
            panel.grid.major.y = element_line(color="grey93", linewidth=0.3))
  )
}

# ── Draw all 6 ────────────────────────────────────────────────────────────────
km_list <- list(
  make_km(df1, "GSE10846 (Training)"),
  make_km(df2, "GSE87371"),
  make_km(df3, "GSE11318"),
  make_km(df4, "GSE181063"),
  make_km(df5, "NCICCR-DLBCL"),
  make_km(df6, "TCGA-DLBC")
)
cat("All 6 KM plots created.\n")

# ── Combine: each panel = KM plot + risk table stacked ────────────────────────
panels <- lapply(km_list, function(km) {
  plot_grid(km$plot, km$table, ncol=1, rel_heights=c(3.2, 0.75),
            align="v", axis="lr")
})

combined <- plot_grid(
  panels[[1]], panels[[2]], panels[[3]],
  panels[[4]], panels[[5]], panels[[6]],
  ncol=3, nrow=2,
  labels=c("A","B","C","D","E","F"),
  label_size=14, label_fontface="bold",
  label_x=0.01, label_y=0.99
)

title_grob <- ggdraw() +
  draw_label(
    "11-Gene Prognostic Risk Score: Kaplan-Meier Analysis Across 6 DLBCL Cohorts",
    fontface="bold", size=13, x=0.5, hjust=0.5, color="grey10"
  )

final_fig5 <- plot_grid(title_grob, combined, ncol=1, rel_heights=c(0.035, 1))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.png"),
       final_fig5, width=21, height=14, dpi=300)
ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.svg"),
       final_fig5, width=21, height=14)
cat("Fig5 saved.\n")

install.packages("svglite", quiet=TRUE)
library(svglite)

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.svg"),
       final_fig5, width=21, height=14, device=svglite)
cat("SVG saved.\n")

# ── Save MultiCohort KM statistics CSV ───────────────────────────────────────
cohort_list <- list(
  list(name="GSE10846 (Training)",  df=df1, platform="GPL570 microarray",    n_total=412),
  list(name="GSE87371",             df=df2, platform="GPL570 microarray",    n_total=221),
  list(name="GSE11318",             df=df3, platform="GPL570 microarray",    n_total=199),
  list(name="GSE181063",            df=df4, platform="GPL14951 microarray",  n_total=1144),
  list(name="NCICCR-DLBCL",        df=df5, platform="RNA-seq (STAR counts)", n_total=234),
  list(name="TCGA-DLBC",           df=df6, platform="RNA-seq (STAR counts)", n_total=45)
)

stats_rows <- lapply(cohort_list, function(x) {
  df  <- x$df
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  cox <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc  <- summary(cox)
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval_lr  <- 1 - pchisq(lr$chisq, df=1)
  pval_cox <- sc$coefficients[1,5]
  data.frame(
    Cohort      = x$name,
    Platform    = x$platform,
    N_total     = nrow(df),
    N_events    = sum(df$OS_event),
    N_Low       = sum(df$RiskGroup=="Low"),
    N_High      = sum(df$RiskGroup=="High"),
    HR          = round(sc$conf.int[1,1], 3),
    CI_low      = round(sc$conf.int[1,3], 3),
    CI_high     = round(sc$conf.int[1,4], 3),
    Cox_pval    = signif(pval_cox, 3),
    LogRank_pval = signif(pval_lr, 3),
    stringsAsFactors = FALSE
  )
})

stats_df <- do.call(rbind, stats_rows)
cat("\nMulti-cohort KM statistics:\n")
print(stats_df)

write.csv(stats_df,
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_KM_statistics.csv"),
          row.names=FALSE)
cat("\nStatistics CSV saved.\n")

cohort_list <- list(
  list(name="GSE10846 (Training)",  df=df1, platform="GPL570 microarray"),
  list(name="GSE87371",             df=df2, platform="GPL570 microarray"),
  list(name="GSE11318",             df=df3, platform="GPL570 microarray"),
  list(name="GSE181063",            df=df4, platform="GPL14951 microarray"),
  list(name="NCICCR-DLBCL",        df=df5, platform="RNA-seq (STAR counts)"),
  list(name="TCGA-DLBC",           df=df6, platform="RNA-seq (STAR counts)")
)

stats_rows <- lapply(cohort_list, function(x) {
  df  <- x$df
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  cox <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc  <- summary(cox)
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  data.frame(
    Cohort       = x$name,
    Platform     = x$platform,
    N_total      = nrow(df),
    N_events     = sum(df$OS_event),
    N_Low        = sum(df$RiskGroup=="Low"),
    N_High       = sum(df$RiskGroup=="High"),
    HR           = round(sc$conf.int[1,1], 3),
    CI_low       = round(sc$conf.int[1,3], 3),
    CI_high      = round(sc$conf.int[1,4], 3),
    Cox_pval     = signif(sc$coefficients[1,5], 3),
    LogRank_pval = signif(1 - pchisq(lr$chisq, df=1), 3),
    stringsAsFactors = FALSE
  )
})

stats_df <- do.call(rbind, stats_rows)
cat("Multi-cohort KM statistics:\n")
print(stats_df[, c("Cohort","N_total","N_events","HR","CI_low","CI_high","Cox_pval","LogRank_pval")])

write.csv(stats_df,
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_KM_statistics.csv"),
          row.names=FALSE)
cat("\nStatistics CSV saved.\n")

library(survival); library(survminer); library(cowplot); library(ggplot2)

sel <- dplyr::select

make_km_v2 <- function(df, title_str) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  n_lo <- sum(df$RiskGroup=="Low");  n_hi <- sum(df$RiskGroup=="High")
  
  fit <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  cox <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc  <- summary(cox)
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  HR   <- round(sc$conf.int[1,1], 2)
  CIlo <- round(sc$conf.int[1,3], 2)
  CIhi <- round(sc$conf.int[1,4], 2)
  pval <- 1 - pchisq(lr$chisq, df=1)
  pstr <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  subtitle_str <- sprintf("HR = %.2f (%.2f\u2013%.2f), %s", HR, CIlo, CIhi, pstr)
  
  # strata names as they appear in survfit object
  strata_names <- names(fit$strata)  # e.g. "RiskGroup=Low", "RiskGroup=High"
  # positional palette: Low=blue, High=red
  pal <- c("#2166AC", "#D73027")
  names(pal) <- strata_names
  
  lab_lo <- sprintf("Low (n=%d)", n_lo)
  lab_hi <- sprintf("High (n=%d)", n_hi)
  
  p <- ggsurvplot(
    fit, data=df,
    palette        = pal,
    legend.labs    = c(lab_lo, lab_hi),
    legend.title   = "",
    legend         = c(0.78, 0.88),
    xlab           = "Time (years)",
    ylab           = "Overall Survival",
    title          = title_str,
    subtitle       = subtitle_str,
    font.title     = c(11, "bold"),
    font.subtitle  = c(8.5, "plain", "grey30"),
    font.x         = c(9.5, "plain"),
    font.y         = c(9.5, "plain"),
    font.tickslab  = c(8.5, "plain"),
    font.legend    = c(8.5, "plain"),
    risk.table     = TRUE,
    risk.table.height = 0.20,
    risk.table.fontsize = 2.8,
    risk.table.title = "No. at risk",
    tables.theme   = theme_cleantable(),
    conf.int       = TRUE,
    conf.int.alpha = 0.12,
    size           = 0.85,
    ggtheme        = theme_classic(base_size=10) +
      theme(plot.title    = element_text(face="bold", size=11, hjust=0),
            plot.subtitle = element_text(size=8.5, color="grey35", hjust=0),
            legend.background = element_rect(fill=alpha("white",0.75), color=NA),
            panel.grid.major.y = element_line(color="grey93", linewidth=0.3))
  )
  p
}

# Redraw all 6
km2 <- list(
  make_km_v2(df1, "GSE10846 (Training)"),
  make_km_v2(df2, "GSE87371"),
  make_km_v2(df3, "GSE11318"),
  make_km_v2(df4, "GSE181063"),
  make_km_v2(df5, "NCICCR-DLBCL"),
  make_km_v2(df6, "TCGA-DLBC")
)
cat("Plots created.\n")

panels2 <- lapply(km2, function(km) {
  plot_grid(km$plot, km$table, ncol=1, rel_heights=c(3.2, 0.75), align="v", axis="lr")
})

combined2 <- plot_grid(
  panels2[[1]], panels2[[2]], panels2[[3]],
  panels2[[4]], panels2[[5]], panels2[[6]],
  ncol=3, nrow=2,
  labels=c("A","B","C","D","E","F"),
  label_size=14, label_fontface="bold",
  label_x=0.01, label_y=0.99
)

title_grob2 <- ggdraw() +
  draw_label("11-Gene Prognostic Risk Score: Kaplan-Meier Analysis Across 6 DLBCL Cohorts",
             fontface="bold", size=13, x=0.5, hjust=0.5, color="grey10")

final_fig5v2 <- plot_grid(title_grob2, combined2, ncol=1, rel_heights=c(0.035, 1))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.png"),
       final_fig5v2, width=21, height=14, dpi=300)
cat("Fig5 v2 saved.\n")

library(survival); library(survminer); library(cowplot); library(ggplot2)

sel <- dplyr::select

make_km_v3 <- function(df, title_str) {
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  n_lo <- sum(df$RiskGroup=="Low");  n_hi <- sum(df$RiskGroup=="High")
  
  fit <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  cox <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc  <- summary(cox)
  lr  <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  HR   <- round(sc$conf.int[1,1], 2)
  CIlo <- round(sc$conf.int[1,3], 2)
  CIhi <- round(sc$conf.int[1,4], 2)
  pval <- 1 - pchisq(lr$chisq, df=1)
  pstr <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  subtitle_str <- sprintf("HR = %.2f (%.2f\u2013%.2f), %s", HR, CIlo, CIhi, pstr)
  
  # Use unnamed vector — ggsurvplot assigns positionally (Low first, High second)
  p <- ggsurvplot(
    fit, data=df,
    palette        = c("#2166AC", "#D73027"),   # positional: Low=blue, High=red
    legend.labs    = c(sprintf("Low (n=%d)", n_lo), sprintf("High (n=%d)", n_hi)),
    legend.title   = "",
    legend         = c(0.78, 0.88),
    xlab           = "Time (years)",
    ylab           = "Overall Survival",
    title          = title_str,
    subtitle       = subtitle_str,
    font.title     = c(11, "bold"),
    font.subtitle  = c(8.5, "plain", "grey30"),
    font.x         = c(9.5, "plain"),
    font.y         = c(9.5, "plain"),
    font.tickslab  = c(8.5, "plain"),
    font.legend    = c(8.5, "plain"),
    risk.table     = TRUE,
    risk.table.height = 0.20,
    risk.table.fontsize = 2.8,
    risk.table.title = "No. at risk",
    tables.theme   = theme_cleantable(),
    conf.int       = TRUE,
    conf.int.alpha = 0.12,
    size           = 0.85,
    ggtheme        = theme_classic(base_size=10) +
      theme(plot.title    = element_text(face="bold", size=11, hjust=0),
            plot.subtitle = element_text(size=8.5, color="grey35", hjust=0),
            legend.background = element_rect(fill=alpha("white",0.75), color=NA),
            panel.grid.major.y = element_line(color="grey93", linewidth=0.3))
  )
  
  # Force colours on the plot object directly
  p$plot <- p$plot +
    scale_color_manual(values=c("#2166AC","#D73027"),
                       labels=c(sprintf("Low (n=%d)",n_lo), sprintf("High (n=%d)",n_hi))) +
    scale_fill_manual(values=c("#2166AC","#D73027"),
                      labels=c(sprintf("Low (n=%d)",n_lo), sprintf("High (n=%d)",n_hi)))
  p
}

# Redraw all 6
km3 <- list(
  make_km_v3(df1, "GSE10846 (Training)"),
  make_km_v3(df2, "GSE87371"),
  make_km_v3(df3, "GSE11318"),
  make_km_v3(df4, "GSE181063"),
  make_km_v3(df5, "NCICCR-DLBCL"),
  make_km_v3(df6, "TCGA-DLBC")
)
cat("Plots created.\n")

panels3 <- lapply(km3, function(km) {
  plot_grid(km$plot, km$table, ncol=1, rel_heights=c(3.2, 0.75), align="v", axis="lr")
})

combined3 <- plot_grid(
  panels3[[1]], panels3[[2]], panels3[[3]],
  panels3[[4]], panels3[[5]], panels3[[6]],
  ncol=3, nrow=2,
  labels=c("A","B","C","D","E","F"),
  label_size=14, label_fontface="bold",
  label_x=0.01, label_y=0.99
)

title_grob3 <- ggdraw() +
  draw_label("11-Gene Prognostic Risk Score: Kaplan-Meier Analysis Across 6 DLBCL Cohorts",
             fontface="bold", size=13, x=0.5, hjust=0.5, color="grey10")

final_fig5v3 <- plot_grid(title_grob3, combined3, ncol=1, rel_heights=c(0.035, 1))

ggsave(translate_path("/mnt/results/DLBCL_prognosis_v2/figures/Fig5_MultiCohort_KM.png"),
       final_fig5v3, width=21, height=14, dpi=300)
cat("Fig5 v3 saved.\n")

library(survival); library(dplyr)

cohort_list <- list(
  list(name="GSE10846 (Training)",  df=df1, platform="GPL570 microarray",     type="Training"),
  list(name="GSE87371",             df=df2, platform="GPL570 microarray",     type="Validation"),
  list(name="GSE11318",             df=df3, platform="GPL570 microarray",     type="Validation"),
  list(name="GSE181063",            df=df4, platform="GPL14951 microarray",   type="Validation"),
  list(name="NCICCR-DLBCL",        df=df5, platform="RNA-seq (STAR counts)", type="Validation"),
  list(name="TCGA-DLBC",           df=df6, platform="RNA-seq (STAR counts)", type="Validation")
)

stats_rows <- lapply(cohort_list, function(x) {
  df  <- x$df
  df$RiskGroup <- factor(df$RiskGroup, levels=c("Low","High"))
  
  # Binary RiskGroup Cox (for KM plot)
  cox_bin <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc_bin  <- summary(cox_bin)
  
  # Continuous RiskScore Cox
  cox_cont <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=df)
  sc_cont  <- summary(cox_cont)
  
  # Log-rank
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval_lr <- 1 - pchisq(lr$chisq, df=1)
  
  data.frame(
    Cohort            = x$name,
    Type              = x$type,
    Platform          = x$platform,
    Scaling           = "Within-cohort z-score (cohort-specific median cutpoint)",
    N_total           = nrow(df),
    N_events          = sum(df$OS_event),
    N_Low             = sum(df$RiskGroup=="Low"),
    N_High            = sum(df$RiskGroup=="High"),
    HR_binary         = round(sc_bin$conf.int[1,1], 3),
    CI_low_binary     = round(sc_bin$conf.int[1,3], 3),
    CI_high_binary    = round(sc_bin$conf.int[1,4], 3),
    Pval_binary_Cox   = signif(sc_bin$coefficients[1,5], 3),
    HR_continuous     = round(sc_cont$conf.int[1,1], 3),
    CI_low_continuous = round(sc_cont$conf.int[1,3], 3),
    CI_high_continuous= round(sc_cont$conf.int[1,4], 3),
    Pval_continuous_Cox = signif(sc_cont$coefficients[1,5], 3),
    LogRank_pval      = signif(pval_lr, 3),
    stringsAsFactors  = FALSE
  )
})

stats_df2 <- do.call(rbind, stats_rows)

cat("Final multi-cohort statistics (binary HR for KM, continuous HR for Cox):\n")
print(stats_df2[, c("Cohort","Type","N_total","N_events",
                     "HR_binary","CI_low_binary","CI_high_binary","Pval_binary_Cox",
                     "HR_continuous","Pval_continuous_Cox","LogRank_pval")])

write.csv(stats_df2,
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_KM_statistics.csv"),
          row.names=FALSE)
cat("\nUpdated statistics CSV saved.\n")

# ── Methodology note ──────────────────────────────────────────────────────────
note <- paste(
  "METHODOLOGY NOTES:",
  "1. Expression scaling: Within-cohort z-score normalization applied to each validation cohort",
  "   independently. This is standard practice for cross-platform validation (microarray vs RNA-seq).",
  "   Consequence: cohort-specific median cutpoints guarantee balanced Low/High groups (50/50 split).",
  "   This should be acknowledged as a limitation vs. fixed-cutpoint external validation.",
  "2. NCICCR-DLBCL survival time: days_to_last_follow_up used for both Dead and Alive cases.",
  "   For Dead cases, GDC convention stores date of death in days_to_last_follow_up when",
  "   days_to_death is missing (confirmed by GDC data model documentation).",
  "3. HR_binary: Cox regression on binary RiskGroup (High vs Low), used for KM plot annotations.",
  "   HR_continuous: Cox regression on continuous RiskScore, used for model performance assessment.",
  "4. TCGA-DLBC: n=45 with only 9 events — severely underpowered; non-significant result expected.",
  sep="\n"
)
writeLines(note, translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_methodology_notes.txt"))
cat("Methodology notes saved.\n")

library(survival); library(dplyr)

# Use full clinical objects (which have RiskScore column)
nci_full <- read.csv(translate_path("/mnt/results/DLBCL_prognosis_v2/tables/NCICCR_DLBCL_risk_scores.csv")) %>%
  dplyr::rename(sample_id=submitter_id) %>%
  dplyr::mutate(RiskGroup=factor(RiskGroup, levels=c("Low","High")))

cohort_list2 <- list(
  list(name="GSE10846 (Training)", df=clin_all,   platform="GPL570 microarray",     type="Training"),
  list(name="GSE87371",            df=clin87,      platform="GPL570 microarray",     type="Validation"),
  list(name="GSE11318",            df=clin113,     platform="GPL570 microarray",     type="Validation"),
  list(name="GSE181063",           df=clin181,     platform="GPL14951 microarray",   type="Validation"),
  list(name="NCICCR-DLBCL",       df=nci_full,    platform="RNA-seq (STAR counts)", type="Validation"),
  list(name="TCGA-DLBC",          df=clin_tcga,   platform="RNA-seq (STAR counts)", type="Validation")
)

stats_rows2 <- lapply(cohort_list2, function(x) {
  df <- x$df %>%
    dplyr::filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event) &
                  !is.na(RiskScore) & !is.na(RiskGroup)) %>%
    dplyr::mutate(RiskGroup=factor(RiskGroup, levels=c("Low","High")))
  
  # Binary RiskGroup Cox
  cox_bin <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  sc_bin  <- summary(cox_bin)
  
  # Continuous RiskScore Cox
  cox_cont <- coxph(Surv(OS_years, OS_event) ~ RiskScore, data=df)
  sc_cont  <- summary(cox_cont)
  
  # Log-rank
  lr      <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data=df)
  pval_lr <- 1 - pchisq(lr$chisq, df=1)
  
  data.frame(
    Cohort              = x$name,
    Type                = x$type,
    Platform            = x$platform,
    Scaling             = "Within-cohort z-score; cohort-specific median cutpoint",
    N_total             = nrow(df),
    N_events            = sum(df$OS_event),
    N_Low               = sum(df$RiskGroup=="Low"),
    N_High              = sum(df$RiskGroup=="High"),
    HR_binary           = round(sc_bin$conf.int[1,1], 3),
    CI_low_binary       = round(sc_bin$conf.int[1,3], 3),
    CI_high_binary      = round(sc_bin$conf.int[1,4], 3),
    Pval_binary_Cox     = signif(sc_bin$coefficients[1,5], 3),
    HR_continuous       = round(sc_cont$conf.int[1,1], 3),
    CI_low_continuous   = round(sc_cont$conf.int[1,3], 3),
    CI_high_continuous  = round(sc_cont$conf.int[1,4], 3),
    Pval_continuous_Cox = signif(sc_cont$coefficients[1,5], 3),
    LogRank_pval        = signif(pval_lr, 3),
    stringsAsFactors    = FALSE
  )
})

stats_df2 <- do.call(rbind, stats_rows2)

cat("Final multi-cohort statistics:\n")
print(stats_df2[, c("Cohort","Type","N_total","N_events",
                     "HR_binary","CI_low_binary","CI_high_binary","Pval_binary_Cox",
                     "HR_continuous","Pval_continuous_Cox","LogRank_pval")])

write.csv(stats_df2,
          translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_KM_statistics.csv"),
          row.names=FALSE)

# Methodology notes
note <- c(
  "METHODOLOGY NOTES",
  "==================",
  "1. Expression scaling: Within-cohort z-score normalization applied independently to each cohort.",
  "   Required due to cross-platform incompatibility (microarray vs RNA-seq).",
  "   Consequence: cohort-specific median cutpoints guarantee balanced Low/High groups (~50/50).",
  "   Limitation: not equivalent to fixed-cutpoint external validation.",
  "",
  "2. NCICCR-DLBCL survival time: days_to_last_follow_up used for both Dead and Alive cases.",
  "   GDC convention: for Dead cases, days_to_last_follow_up records date of death when",
  "   days_to_death is missing (standard GDC data model).",
  "",
  "3. HR_binary: Cox on binary RiskGroup (High vs Low) — used for KM plot HR annotations.",
  "   HR_continuous: Cox on continuous RiskScore — used for model discrimination assessment.",
  "",
  "4. TCGA-DLBC: n=45, only 9 events — severely underpowered; non-significant result expected",
  "   and does not indicate model failure."
)
writeLines(note, translate_path("/mnt/results/DLBCL_prognosis_v2/tables/MultiCohort_methodology_notes.txt"))
cat("\nAll files saved.\n")
