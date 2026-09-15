# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
# 注意：本仓库路径含中文，R 在 Windows 下对含中文的**绝对路径**做编码
# 转换会失败（unable to translate to UTF-8），因此统一使用相对路径，
# 并确保在仓库根目录运行。
source(file.path("config", "paths.R"))
# <<< 自动注入结束 >>>

# ════════════════════════════════════════════════════════════════════════════════
# DLBCL 预后分析 — 完整一键运行脚本
# 数据来源：Agent 生成的 CSV（D:/bulk-download/DLBCL_prognosis/tables/）
# 图表：Part A（Agent 风格）+ Part B（你原来的风格）
# ════════════════════════════════════════════════════════════════════════════════

# rm(list = ls())   # 迁移工具已注释：会清除注入的路径配置
# ── 0. 加载包 ─────────────────────────────────────────────────────────────────
pkgs <- c("dplyr", "stringr", "survival", "survminer", "ggplot2", "patchwork",
          "RColorBrewer", "rms", "timeROC", "Hmisc", "cowplot", "tinyarray")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ── 1. 路径设置 ───────────────────────────────────────────────────────────────
data_dir <- translate_path("D:/bulk-download/DLBCL_prognosis/tables")
fig_dir  <- translate_path("D:/bulk-download/DLBCL_prognosis/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 2. 读入所有 CSV ───────────────────────────────────────────────────────────
lasso_formula <- read.csv(file.path(data_dir, "LASSO_prognostic_formula.csv"), stringsAsFactors = FALSE)
lasso_coefs   <- setNames(lasso_formula$Coefficient, lasso_formula$Gene)
lasso_genes   <- names(lasso_coefs)
cat("LASSO 基因 (n=", length(lasso_genes), "):", paste(lasso_genes, collapse = ", "), "\n")

cox_raw    <- read.csv(file.path(data_dir, "Cox_univariable_multivariable_results.csv"), stringsAsFactors = FALSE)
cal_raw    <- read.csv(file.path(data_dir, "Calibration_data_1_3_5yr.csv"),              stringsAsFactors = FALSE)
immune_raw <- read.csv(file.path(data_dir, "Immune_infiltration_scores.csv"),             stringsAsFactors = FALSE)
cor_raw    <- read.csv(file.path(data_dir, "Immune_RiskScore_correlations.csv"),          stringsAsFactors = FALSE)

# ── 3. 读入各队列数据并整理 ───────────────────────────────────────────────────
prep_cohort <- function(path, cohort_name) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  # 统一列名
  if ("submitter_id" %in% colnames(df)) df <- rename(df, sample_id = submitter_id)
  df$fp        <- df$RiskScore
  df$Risk      <- factor(ifelse(df$RiskScore > median(df$RiskScore, na.rm = TRUE),
                                "high_risk", "low_risk"), levels = c("low_risk", "high_risk"))
  df$RiskGroup <- ifelse(df$RiskScore >= median(df$RiskScore, na.rm = TRUE), "High", "Low")
  df$time      <- df$OS_years
  df$event     <- df$OS_event
  df$Cohort    <- cohort_name
  df <- df[!is.na(df$time) & df$time > 0 & !is.na(df$event), ]
  df
}

dat    <- prep_cohort(file.path(data_dir, "GSE10846_risk_scores.csv"),  "GSE10846")
dat87  <- prep_cohort(file.path(data_dir, "GSE87371_risk_scores.csv"),  "GSE87371")
dat11  <- prep_cohort(file.path(data_dir, "GSE11318_risk_scores.csv"),  "GSE11318")
dat181 <- prep_cohort(file.path(data_dir, "GSE181063_risk_scores.csv"), "GSE181063")
dat_tcga <- prep_cohort(file.path(data_dir, "TCGA_DLBC_risk_scores.csv"), "TCGA-DLBC")

cat("队列样本量：\n")
for (d in list(dat, dat87, dat11, dat181, dat_tcga)) {
  cat(sprintf("  %-12s n=%-5d events=%d\n", d$Cohort[1], nrow(d), sum(d$event)))
}

# ════════════════════════════════════════════════════════════════════════════════
# PART A：Agent 风格图
# ════════════════════════════════════════════════════════════════════════════════

cat("\n====== PART A: Agent 风格图 ======\n")

# ── A1. Cox 森林图（单因素 + 多因素）────────────────────────────────────────
cat("A1: Cox 森林图...\n")

var_labels <- c(
  RiskScore = "Risk Score",
  Age60     = "Age (>=60 vs <60)",
  Male      = "Sex (Male vs Female)",
  Stage_adv = "Stage (III-IV vs I-II)",
  IPI_high  = "IPI (High vs Low)"
)
sig_colors <- c("***​" = "#D73027", "​**" = "#FC8D59", "*" = "#E08000", "ns" = "grey55")

forest_df <- cox_raw %>%
  mutate(
    label    = var_labels[Variable],
    sig      = ifelse(pval < 0.001, "***​", ifelse(pval < 0.01, "​**",
               ifelse(pval < 0.05, "*", "ns"))),
    pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
    hr_str   = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
    y_pos    = c(13, 12, 11, 10, 9, 7, 6, 5, 4, 3)
  )

p_left <- ggplot(forest_df, aes(y = y_pos)) +
  geom_text(aes(label = label), x = 0, hjust = 0, size = 4.6, color = "black") +
  annotate("text", x = 0, y = 14.3, label = "Univariable",
           hjust = 0, size = 5.0, fontface = "bold", color = "#2166AC") +
  annotate("text", x = 0, y = 8.3,  label = "Multivariable",
           hjust = 0, size = 5.0, fontface = "bold", color = "#D73027") +
  annotate("text", x = 0, y = 15.5, label = "Variable",
           hjust = 0, size = 4.8, fontface = "bold", color = "grey25") +
  geom_hline(yintercept = c(15.0, 8.7), color = "grey65", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(2, 16.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_void() +
  theme(plot.margin = margin(t = 10, r = 0, b = 10, l = 10))

p_forest_mid <- ggplot(forest_df, aes(y = y_pos, x = HR, xmin = CI_low, xmax = CI_high)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbar(aes(color = sig), width = 0.3, linewidth = 0.9, orientation = "y") +
  geom_point(aes(color = sig, size = sig), shape = 18) +
  geom_hline(yintercept = c(15.0, 8.7), color = "grey65", linewidth = 0.4) +
  annotate("text", x = 1.5, y = 15.5, label = "Hazard Ratio",
           hjust = 0.5, size = 4.8, fontface = "bold", color = "grey25") +
  scale_color_manual(values = sig_colors, guide = "none") +
  scale_size_manual(values = c("***​" = 5, "​**" = 4.5, "*" = 4, "ns" = 3.5), guide = "none") +
  scale_x_log10(limits = c(0.2, 12), breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8"),
                name = "Hazard Ratio (log scale)") +
  scale_y_continuous(limits = c(2, 16.5), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_classic(base_size = 15) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 14),
        panel.grid.major.x = element_line(color = "grey93", linewidth = 0.4),
        plot.margin = margin(t = 10, r = 5, b = 10, l = 5))

p_hr_col <- ggplot(forest_df, aes(y = y_pos)) +
  geom_text(aes(label = hr_str), x = 0.5, hjust = 0.5, size = 4.2, color = "grey20") +
  geom_hline(yintercept = c(15.0, 8.7), color = "grey65", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 15.5, label = "HR (95% CI)",
           hjust = 0.5, size = 4.8, fontface = "bold", color = "grey25") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(2, 16.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_void() +
  theme(plot.margin = margin(t = 10, r = 5, b = 10, l = 5))

p_pv_col <- ggplot(forest_df, aes(y = y_pos)) +
  geom_text(aes(label = pval_str, color = sig), x = 0.5, hjust = 0.5,
            size = 4.2, fontface = "bold") +
  geom_hline(yintercept = c(15.0, 8.7), color = "grey65", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 15.5, label = "P value",
           hjust = 0.5, size = 4.8, fontface = "bold", color = "grey25") +
  scale_color_manual(values = sig_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(2, 16.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_void() +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 5))

p_A1 <- p_left + p_forest_mid + p_hr_col + p_pv_col +
  plot_layout(widths = c(2.2, 2.5, 1.8, 1.0)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))  # ★ 标题去除

ggsave(file.path(fig_dir, "A1_Cox_forest.png"), p_A1, width = 6.5, height = 4.6, dpi = 300)  # ★ 紧凑版：85mm 栏宽下正文字号约 7-9pt
cat("  已保存 A1_Cox_forest.png\n")

# ── A2. 诺莫图（rms 重跑）────────────────────────────────────────────────────
cat("A2: 诺莫图...\n")

# 准备数据（时间转为天，rms 要求）
nom_df <- dat %>%
  filter(!is.na(IPI) & !is.na(RiskScore) & time > 0) %>%
  mutate(time_day = time * 365.25) %>%
  select(time_day, event, RiskScore, IPI)

dd <- datadist(nom_df)
options(datadist = "dd")

fit_cph <- cph(Surv(time_day, event) ~ RiskScore + IPI,
               data = nom_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)

S0_1 <- survest(fit_cph, newdata = data.frame(RiskScore = 0, IPI = 0), times = 365)$surv
S0_3 <- survest(fit_cph, newdata = data.frame(RiskScore = 0, IPI = 0), times = 1095)$surv
S0_5 <- survest(fit_cph, newdata = data.frame(RiskScore = 0, IPI = 0), times = 1825)$surv
center_cph <- fit_cph$center

f1 <- function(lp) 1 - S0_1^exp(lp - center_cph)
f3 <- function(lp) 1 - S0_3^exp(lp - center_cph)
f5 <- function(lp) 1 - S0_5^exp(lp - center_cph)

nom <- nomogram(fit_cph,
                fun      = list(f1, f3, f5),
                funlabel = c("1-Year Mortality", "3-Year Mortality", "5-Year Mortality"),
                fun.at   = list(seq(0.05, 0.9, 0.1), seq(0.05, 0.9, 0.1), seq(0.05, 0.9, 0.1)),
                lp = FALSE)

png(file.path(fig_dir, "A2_Nomogram.png"), width = 3200, height = 1800, res = 300)
par(mar = c(2, 2, 3, 2), cex = 0.95)
plot(nom, xfrac = 0.32, label.every = 1, col.grid = grey(c(0.8, 0.95)),
     main = "Nomogram: 1/3/5-Year Overall Survival\n(IPI + Risk Score, GSE10846 DLBCL n=405)")
dev.off()
cat("  已保存 A2_Nomogram.png\n")

# ── A3. 1/3/5 年校准曲线（ggplot2 版）────────────────────────────────────────
cat("A3: 校准曲线...\n")

# 用 coxph 预测生存概率，再做十分位组校准
cox_cal <- coxph(Surv(time_day, event) ~ RiskScore + IPI, data = nom_df, x = TRUE)
sf_all  <- survfit(cox_cal, newdata = nom_df)

get_pred_t <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)
  if (is.matrix(s$surv)) s$surv[1, ] else s$surv
}

pred_1 <- get_pred_t(sf_all, 365)
pred_3 <- get_pred_t(sf_all, 1095)
pred_5 <- get_pred_t(sf_all, 1825)

cal_decile <- function(pred_vec, t_days, df, n_grp = 10) {
  df2 <- df %>% mutate(pred = pred_vec, grp = ntile(pred, n_grp))
  df2 %>% group_by(grp) %>% do({
    sub <- .
    km  <- survfit(Surv(time_day, event) ~ 1, data = sub)
    s   <- summary(km, times = t_days, extend = TRUE)
    obs <- ifelse(length(s$surv) > 0,    s$surv[1],    NA)
    se  <- ifelse(length(s$std.err) > 0, s$std.err[1], NA)
    data.frame(mean_pred = mean(sub$pred), obs_surv = obs,
               obs_lo = obs - 1.96 * se, obs_hi = obs + 1.96 * se, n = nrow(sub))
  }) %>% ungroup()
}

cal_all <- bind_rows(
  cal_decile(pred_1, 365,  nom_df) %>% mutate(label = "1-Year OS"),
  cal_decile(pred_3, 1095, nom_df) %>% mutate(label = "3-Year OS"),
  cal_decile(pred_5, 1825, nom_df) %>% mutate(label = "5-Year OS")
) %>% mutate(label = factor(label, levels = c("1-Year OS", "3-Year OS", "5-Year OS")))

time_colors <- c("1-Year OS" = "#2166AC", "3-Year OS" = "#4DAC26", "5-Year OS" = "#D73027")
cal_plot <- cal_all %>% filter(!is.na(obs_surv))
mae_df   <- cal_plot %>% group_by(label) %>%
  summarise(MAE = round(mean(abs(mean_pred - obs_surv)), 3), .groups = "drop")
ann_df   <- mae_df %>% mutate(x = 0.08, y = 0.97, txt = paste0("MAE = ", MAE))

p_A3 <- ggplot(cal_plot, aes(x = mean_pred, y = obs_surv, color = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_errorbar(aes(ymin = pmax(obs_lo, 0), ymax = pmin(obs_hi, 1)),
                width = 0.015, linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 3, shape = 16) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.0, span = 1.2) +
  geom_text(data = ann_df, aes(x = x, y = y, label = txt, color = label),
            hjust = 0, vjust = 1, size = 3.2, fontface = "bold", inherit.aes = FALSE) +
  facet_wrap(~label, ncol = 3) +
  scale_color_manual(values = time_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Predicted Survival") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), name = "Observed Survival (KM)") +
  labs(title    = "Calibration Curves: Predicted vs Observed Survival",
       subtitle = "Decile-based calibration (GSE10846 DLBCL, n=405)") +
  theme_classic(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9.5, color = "grey40"),
        strip.text       = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.4),
        plot.margin      = margin(10, 10, 10, 10))

ggsave(file.path(fig_dir, "A3_Calibration_curves.png"), p_A3, width = 11, height = 4.5, dpi = 300)
cat("  已保存 A3_Calibration_curves.png\n")

# ════════════════════════════════════════════════════════════════════════════════
# A4 修正版：条形图 + 4个固定颜色散点图（完全对照 yuhou.txt）
# ════════════════════════════════════════════════════════════════════════════════

cat("A4: 免疫浸润图（修正版）...
")

risk_immune <- dat %>%
  select(sample_id, RiskScore, RiskGroup) %>%
  inner_join(
    immune_raw %>% select(-any_of(c("RiskScore", "RiskGroup"))),
    by = "sample_id"
  )
cell_types <- names(risk_immune)[!names(risk_immune) %in%
                                   c("sample_id", "RiskScore", "RiskGroup")]

cor_df <- do.call(rbind, lapply(cell_types, function(ct) {
  test <- cor.test(risk_immune$RiskScore, risk_immune[[ct]],
                   method = "spearman", exact = FALSE)
  data.frame(cell_type = ct,
             rho  = as.numeric(test$estimate),
             pval = test$p.value,
             stringsAsFactors = FALSE)
})) %>%
  mutate(
    sig      = ifelse(pval < 0.001, "***​", ifelse(pval < 0.01, "​**",
                                                   ifelse(pval < 0.05, "*", "ns"))),
    pval_str = ifelse(pval < 0.001, "< 0.001", sprintf("= %.3f", pval))
  ) %>%
  arrange(rho)

# ── 散点图辅助函数（完全对照 yuhou.txt scatter_immune）────────────────────────
scatter_immune <- function(df, y_col, color, title_str, rho_val, pval_str) {
  ggplot(df, aes_string(x = "RiskScore", y = y_col)) +
    geom_point(color = color, alpha = 0.45, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = color,
                fill = paste0(color, "33"), linewidth = 1.0) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Spearman r = ", round(rho_val, 3),
                            "\nP ", pval_str),
             hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "bold",
             color = ifelse(rho_val > 0, "#D73027", "#2166AC")) +
    labs(title = title_str,
         x     = "Risk Score",
         y     = paste0(gsub("_", " ", y_col), "\n(ssGSEA score)")) +
    theme_classic(base_size = 15) +
    theme(plot.title       = element_text(face = "bold", size = 11),
          panel.grid.major = element_line(color = "grey93", linewidth = 0.3))
}

# ── 4 个固定散点图，颜色和顺序完全对照原代码 ─────────────────────────────────
get_cor <- function(ct) {
  row <- cor_df[cor_df$cell_type == ct, ]
  list(rho = row$rho, pval_str = row$pval_str)
}

p_m2    <- scatter_immune(risk_immune, "M2_Macrophage", "#D73027",
                          "M2 Macrophage vs Risk Score",
                          get_cor("M2_Macrophage")$rho,
                          get_cor("M2_Macrophage")$pval_str)

p_cd8   <- scatter_immune(risk_immune, "CD8_T_cell",    "#2166AC",
                          "CD8+ T Cell vs Risk Score",
                          get_cor("CD8_T_cell")$rho,
                          get_cor("CD8_T_cell")$pval_str)

p_m1    <- scatter_immune(risk_immune, "M1_Macrophage", "#4DAC26",
                          "M1 Macrophage vs Risk Score",
                          get_cor("M1_Macrophage")$rho,
                          get_cor("M1_Macrophage")$pval_str)

p_bcell <- scatter_immune(risk_immune, "B_cell",        "#984EA3",
                          "B Cell vs Risk Score",
                          get_cor("B_cell")$rho,
                          get_cor("B_cell")$pval_str)

# ── 条形图（全部细胞类型）────────────────────────────────────────────────────
cor_df_plot <- cor_df %>%
  mutate(
    cell_label = gsub("_", " ", cell_type),
    direction  = ifelse(rho > 0, "Positive", "Negative")
  ) %>%
  arrange(rho) %>%
  mutate(cell_label = factor(cell_label, levels = cell_label))

p_bar <- ggplot(cor_df_plot, aes(x = rho, y = cell_label, fill = direction)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sig,
                x     = ifelse(rho >= 0, rho + 0.005, rho - 0.005),
                hjust = ifelse(rho >= 0, 0, 1)),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Positive" = "#D73027", "Negative" = "#2166AC"),
                    name = "Direction") +
  scale_x_continuous(limits = c(-0.25, 0.35), breaks = seq(-0.2, 0.3, 0.1)) +
  labs(title    = "Risk Score vs Immune Cell Infiltration",
       subtitle = "Spearman correlation (GSE10846, n=412)",
       x = "Spearman rho", y = NULL) +
  theme_classic(base_size = 15) +
  theme(plot.title         = element_text(face = "bold", size = 12),
        plot.subtitle      = element_text(size = 9.5, color = "grey40"),
        axis.text.y        = element_text(size = 10),
        panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3))

# ── 组合（上4散点 + 下条形）─────────────────────────────────────────────────
top_row_A4 <- p_m2 + p_cd8 + p_m1 + p_bcell + plot_layout(ncol = 4)
p_A4 <- top_row_A4 / p_bar +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title    = "Immune Infiltration Analysis",
    subtitle = "ssGSEA scores correlated with LASSO Risk Score (GSE10846 DLBCL)",
    theme = theme(plot.title    = element_text(face = "bold", size = 14, hjust = 0),
                  plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0),
                  plot.margin   = margin(10, 10, 10, 10))
  )

ggsave(file.path(fig_dir, "A4_Immune_infiltration.png"), p_A4,
       width = 14, height = 9, dpi = 300)
cat("  已保存 A4_Immune_infiltration.png
")




# ── A5. 多队列 KM 曲线（5 个队列）────────────────────────────────────────────
cat("A5: 多队列 KM 曲线...\n")

cohort_list  <- list(dat, dat87, dat11, dat181, dat_tcga)
cohort_title <- c("GSE10846 (Training, n=412)",
                  "GSE87371 (Validation, n=221)",
                  "GSE11318 (Validation, n=199)",
                  "GSE181063 (Validation, n=882)",
                  "TCGA-DLBC (Validation, n=45)")
xlims        <- c(22, 20, 22, 15, 18)

col_high <- "#E63946"
col_low  <- "#457B9D"

make_km <- function(df, title, xmax_cap) {
  df$RiskGroup <- factor(df$RiskGroup, levels = c("Low", "High"))
  df_plot <- df
  df_plot$OS_event[df_plot$OS_years > xmax_cap] <- 0
  df_plot$OS_years <- pmin(df_plot$OS_years, xmax_cap)

  fit  <- survfit(Surv(OS_years, OS_event) ~ RiskGroup, data = df_plot)
  lr   <- survdiff(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  pval <- 1 - pchisq(lr$chisq, df = 1)
  plab <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
  cx   <- coxph(Surv(OS_years, OS_event) ~ RiskGroup, data = df)
  hr   <- exp(coef(cx)); ci <- exp(confint(cx))
  annot <- sprintf("%s\nHR = %.2f (%.2f\u2013%.2f)", plab, hr, ci[1], ci[2])
  bby   <- if (xmax_cap > 15) 5L else if (xmax_cap > 8) 4L else 3L

  p <- ggsurvplot(
    fit, data = df_plot,
    palette = c(col_low, col_high),
    legend.labs = c("Low Risk", "High Risk"), legend.title = "",
    legend = c(0.78, 0.88),
    xlab = "Time (years)", ylab = "Overall Survival Probability",
    title = title,
    font.title = c(10L, "bold"), font.x = 9L, font.y = 9L,
    font.tickslab = 8L, font.legend = 8L,
    risk.table = TRUE, risk.table.height = 0.25,
    risk.table.fontsize = 3L, risk.table.title = "No. at risk",
    tables.theme = theme_cleantable(),
    ggtheme = theme_classic(base_size = 9L) +
              theme(plot.title = element_text(hjust = 0.5)),
    conf.int = TRUE, conf.int.alpha = 0.12,
    pval = FALSE, surv.median.line = "hv",
    break.time.by = bby, xlim = c(0, xmax_cap)
  )
  p$plot <- p$plot +
    annotate("text", x = xmax_cap * 0.20, y = 0.10,
             label = annot, hjust = 0, vjust = 0, size = 3L, color = "black")
  p
}

plots_km <- mapply(make_km, cohort_list, cohort_title, xlims, SIMPLIFY = FALSE)

extract_grob <- function(p) {
  cowplot::plot_grid(p$plot, p$table, ncol = 1, rel_heights = c(0.75, 0.25))
}
grobs_km <- lapply(plots_km, extract_grob)

top_row_km    <- plot_grid(grobs_km[[1]], grobs_km[[2]], grobs_km[[3]], nrow = 1, ncol = 3)
spacer        <- ggplot() + theme_void()
bottom_row_km <- plot_grid(spacer, grobs_km[[4]], grobs_km[[5]], spacer,
                           nrow = 1, ncol = 4, rel_widths = c(0.5, 1, 1, 0.5))
final_km      <- plot_grid(top_row_km, bottom_row_km, nrow = 2, ncol = 1)

ggsave(file.path(fig_dir, "A5_MultiCohort_KM.png"), final_km,
       width = 18, height = 12, dpi = 180)
cat("  已保存 A5_MultiCohort_KM.png\n")

# ════════════════════════════════════════════════════════════════════════════════
# PART B：你原来的风格图（基于 Agent CSV 数据）
# ════════════════════════════════════════════════════════════════════════════════

cat("\n====== PART B: 你原来风格的图 ======\n")

# ════════════════════════════════════════════════════════════════════════════════
# B1 修正版：从本地 Rdata 读取基因表达矩阵，完全对照原代码
# ════════════════════════════════════════════════════════════════════════════════

cat("B1: 风险评分三联图（修正版）...
")

library(tinyarray)

# ── 从本地 Rdata 读取原始数据（你原来的文件）────────────────────────────────
rdata_path <- translate_path("D:/bulk-download/DLBCL_prognosis/GSE10846_sur_model.Rdata")
# 如果路径不对，改成你本地实际路径
if (!file.exists(rdata_path)) {
  stop("找不到 GSE10846_sur_model.Rdata，请修改 rdata_path 为正确路径")
}
load(rdata_path)   # 加载进来应该有 exprSet 和 meta

# LASSO 13 基因
gene_cols <- intersect(lasso_genes, rownames(exprSet))
cat("找到 LASSO 基因:", length(gene_cols), "个:", paste(gene_cols, collapse = ", "), "
")

exprSet_hub <- exprSet[gene_cols, ]   # genes × samples

# ── 风险评分（用 dat 里已有的 RiskScore）─────────────────────────────────────
riskscore <- dat$RiskScore
names(riskscore) <- dat$sample_id

# 确保样本顺序一致
common_s  <- intersect(names(riskscore), colnames(exprSet_hub))
riskscore <- riskscore[common_s]
exprSet_hub <- exprSet_hub[, common_s]
meta_sub    <- meta[match(common_s, rownames(meta)), ]

cut_val <- median(riskscore)
color   <- c("#3288BD", "#D73027")

# ── 构建三个子图数据 ──────────────────────────────────────────────────────────
fp_dat <- data.frame(
  patientid = seq_len(length(riskscore)),
  riskscore = as.numeric(sort(riskscore)),
  ri        = factor(ifelse(sort(riskscore) > cut_val, "high", "low"),
                     levels = c("low", "high"))
)

sur_dat <- data.frame(
  patientid = seq_len(length(riskscore)),
  time      = meta_sub[order(riskscore), "time"],
  event     = meta_sub[order(riskscore), "event"]
)
sur_dat$event <- ifelse(sur_dat$event == 0, "alive", "death")

exp_dat <- scale(t(exprSet_hub[, order(riskscore)]))
n_cutoff <- 3
exp_dat[exp_dat >  n_cutoff] <-  n_cutoff
exp_dat[exp_dat < -n_cutoff] <- -n_cutoff

risk_order <- factor(ifelse(sort(riskscore) > cut_val, "high", "low"),
                     levels = c("low", "high"))

# ── p1: 风险评分点图 ──────────────────────────────────────────────────────────
p1 <- ggplot(fp_dat, aes(x = patientid, y = riskscore, color = ri)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = color) +
  geom_vline(xintercept = sum(riskscore < cut_val), lty = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.line        = element_line(color = "black"),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "                                     GSE10846", x = "", y = "Riskscore")

# ── p2: 生存时间点图（alive/death 颜色，Y 轴 0-120）─────────────────────────
p2 <- ggplot(sur_dat, aes(x = patientid, y = time)) +
  geom_point(aes(col = event), size = 0.8) +
  geom_vline(xintercept = sum(riskscore < cut_val), lty = 2) +
  scale_color_manual(values = color) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120)) +
  theme_bw() +
  labs(color = "Event", y = "Time", x = "") +
  theme(axis.line        = element_line(color = "black"),
        plot.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ── p3: 热图 ─────────────────────────────────────────────────────────────────
p3 <- ggheat(exp_dat,
             risk_order,
             show_rownames = TRUE,
             color         = c(color[1], "white", color[2]),
             legend_color  = color,
             groupname     = "Risk",
             expname       = "Expression")

n_height <- ifelse(length(gene_cols) < 15, 3,
                   ifelse(length(gene_cols) < 20, 4, 5))

p_B1 <- wrap_plots(p1, p2, p3, ncol = 1, heights = c(1, 1, n_height))

ggsave(file.path(fig_dir, "B1_RiskScore_tripanel.png"), p_B1,
       width = 10, height = 4 + n_height * 1.5, dpi = 300)
cat("  已保存 B1_RiskScore_tripanel.png
")



# ── B2. KM 生存曲线（训练队列）────────────────────────────────────────────────
cat("B2: KM 生存曲线（训练队列）...\n")

dat$Risk <- factor(ifelse(dat$RiskScore >= median(dat$RiskScore), "High", "Low"),
                   levels = c("Low", "High"))
fit_km <- survfit(Surv(time, event) ~ Risk, data = dat)
lr_km  <- survdiff(Surv(time, event) ~ Risk, data = dat)
pval_km <- 1 - pchisq(lr_km$chisq, df = 1)
plab_km <- ifelse(pval_km < 0.001, "p < 0.001", sprintf("p = %.3f", pval_km))
cx_km   <- coxph(Surv(time, event) ~ Risk, data = dat)
hr_km   <- exp(coef(cx_km)); ci_km <- exp(confint(cx_km))

p_B2 <- ggsurvplot(
  fit_km, data = dat,
  palette = c("#457B9D", "#E63946"),
  legend.labs  = c("Low Risk", "High Risk"), legend.title = "",
  legend       = c(0.80, 0.90),
  xlab = "Time (years)", ylab = "Overall Survival Probability",
  title = "Kaplan-Meier Survival Curves\nGSE10846 DLBCL (Training Cohort, n=412)",
  font.title   = c(12, "bold"), font.x = 11, font.y = 11,
  font.tickslab = 10, font.legend = 10,
  risk.table   = TRUE, risk.table.height = 0.25,
  risk.table.fontsize = 3.5, risk.table.title = "No. at risk",
  tables.theme = theme_cleantable(),
  ggtheme      = theme_classic(base_size = 15) +
                 theme(plot.title = element_text(hjust = 0.5, face = "bold")),
  conf.int     = TRUE, conf.int.alpha = 0.15,
  pval         = FALSE, surv.median.line = "hv",
  break.time.by = 5
)
p_B2$plot <- p_B2$plot +
  annotate("text", x = 1, y = 0.12,
           label = sprintf("%s\nHR = %.2f (%.2f\u2013%.2f)", plab_km, hr_km, ci_km[1], ci_km[2]),
           hjust = 0, vjust = 0, size = 3.8, color = "black")

png(file.path(fig_dir, "B2_KM_training.png"), width = 2800, height = 2200, res = 300)
print(p_B2)
dev.off()
cat("  已保存 B2_KM_training.png\n")

# ── B3. Time-dependent ROC 曲线（1/3/5 年）───────────────────────────────────
cat("B3: Time-ROC 曲线...\n")

library(timeROC)
roc_obj <- timeROC(
  T      = dat$time,
  delta  = dat$event,
  marker = dat$RiskScore,
  cause  = 1,
  times  = c(1, 3, 5),
  iid    = TRUE
)

auc_vals <- roc_obj$AUC
ci_vals  <- confint(roc_obj, level = 0.95)$CI_AUC

roc_df <- do.call(rbind, lapply(1:3, function(i) {
  data.frame(
    FPR   = roc_obj$FP[, i],
    TPR   = roc_obj$TP[, i],
    time  = c(1, 3, 5)[i],
    label = sprintf("%d-Year AUC = %.3f (%.3f\u2013%.3f)",
                    c(1, 3, 5)[i], auc_vals[i], ci_vals[i, 1], ci_vals[i, 2])
  )
})) %>% mutate(label = factor(label, levels = unique(label)))

p_B3 <- ggplot(roc_df, aes(x = FPR, y = TPR, color = label)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#2166AC", "#4DAC26", "#D73027"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "1 - Specificity (FPR)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     name = "Sensitivity (TPR)") +
  labs(title    = "Time-Dependent ROC Curves",
       subtitle = "GSE10846 DLBCL Training Cohort (n=412)") +
  theme_classic(base_size = 12) +
  theme(plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle   = element_text(size = 10, color = "grey40", hjust = 0.5),
        legend.position = c(0.72, 0.18),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text     = element_text(size = 9),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3))

ggsave(file.path(fig_dir, "B3_TimeROC.png"), p_B3, width = 6.5, height = 6, dpi = 300)
cat("  已保存 B3_TimeROC.png\n")

# ── B4. 诺莫图（你原来风格，rms）─────────────────────────────────────────────
# （已在 A2 中生成，此处用相同数据生成一个 base R 版本存档）
cat("B4: 诺莫图（base R 版）— 已在 A2 中生成，跳过重复保存\n")

# ════════════════════════════════════════════════════════════════════════════════
# B5 修正版：时间单位统一为天（* 365.25），calibrate u= 对应天数
# ════════════════════════════════════════════════════════════════════════════════

cat("B5: 校准曲线（修正版）...
")

library(rms)

# time 是年，* 365.25 转为天，和 u=365/1095/1825 对应
cal_df <- dat %>%
  filter(!is.na(RiskScore) & time > 0) %>%
  mutate(time_day = time * 365.25)

# 检查是否有临床变量
has_clin <- all(c("age", "gender", "subtype", "stage") %in% colnames(cal_df))

if (has_clin) {
  cal_df$gender  <- factor(cal_df$gender)
  cal_df$stage   <- factor(as.character(cal_df$stage))
  cal_df$subtype <- factor(cal_df$subtype)
  dd_b5 <- datadist(cal_df); options(datadist = "dd_b5")
  formula_b5 <- Surv(time_day, event) ~ age + gender + subtype + stage + RiskScore
  cat("  使用完整临床公式：age + gender + subtype + stage + RiskScore
")
} else {
  dd_b5 <- datadist(cal_df); options(datadist = "dd_b5")
  formula_b5 <- Surv(time_day, event) ~ RiskScore
  cat("  未找到临床变量，仅用 RiskScore
")
}

# 三个时间点分别拟合（u 单位必须和 time_day 一致，即天）
f1_b5 <- cph(formula_b5, data = cal_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
f3_b5 <- cph(formula_b5, data = cal_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095)
f5_b5 <- cph(formula_b5, data = cal_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1825)

cat("  运行 calibrate（B=200，约需 1-2 分钟）...
")
cal1_b5 <- calibrate(f1_b5, cmethod = "KM", method = "boot", u = 365,  m = 40, B = 200)
cal3_b5 <- calibrate(f3_b5, cmethod = "KM", method = "boot", u = 1095, m = 40, B = 200)
cal5_b5 <- calibrate(f5_b5, cmethod = "KM", method = "boot", u = 1825, m = 40, B = 200)

data1_b5 <- data.frame(Time      = "1-year",
                       Predicted = cal1_b5[, "mean.predicted"],
                       Observed  = cal1_b5[, "KM"])
data3_b5 <- data.frame(Time      = "3-year",
                       Predicted = cal3_b5[, "mean.predicted"],
                       Observed  = cal3_b5[, "KM"])
data5_b5 <- data.frame(Time      = "5-year",
                       Predicted = cal5_b5[, "mean.predicted"],
                       Observed  = cal5_b5[, "KM"])

all_cal <- rbind(data1_b5, data3_b5, data5_b5)
all_cal$Time <- factor(all_cal$Time, levels = c("1-year", "3-year", "5-year"))

p_B5 <- ggplot(all_cal, aes(x = Predicted, y = Observed)) +
  geom_abline(intercept = 0, slope = 1,
              size = 2, linetype = 3, color = "darkgrey") +
  geom_line(aes(color = Time)) +
  scale_color_manual(values = c("1-year" = "#e6b707",
                                "3-year" = "#3288BD",
                                "5-year" = "#D73027")) +
  geom_point(aes(shape = Time, fill = Time),
             color = "white", size = 5, stroke = 2) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("1-year" = "#e6b707",
                               "3-year" = "#3288BD",
                               "5-year" = "#D73027")) +
  labs(title = "                                          GSE10846",
       x     = "Predicted Probability",
       y     = "Observed Probability",
       color = "Time") +
  xlim(0, 1) + ylim(0, 1) +
  theme_classic() +
  coord_fixed() +
  theme(legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12))

ggsave(file.path(fig_dir, "B5_Calibration_ggplot.png"), p_B5,
       width = 6, height = 6, dpi = 300)
cat("  已保存 B5_Calibration_ggplot.png
")


# ════════════════════════════════════════════════════════════════════════════════
# 完成汇总
# ════════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════\n")
cat("  全部图表生成完毕！保存路径：\n")
cat(paste0("  ", fig_dir, "\n"))
cat("════════════════════════════════════════════════════\n")

fig_list <- list(
  "A1_Cox_forest.png"          = "Part A | Cox 单/多因素森林图",
  "A2_Nomogram.png"            = "Part A | 诺莫图（1/3/5 年，rms）",
  "A3_Calibration_curves.png"  = "Part A | 校准曲线（十分位 ggplot2）",
  "A4_Immune_infiltration.png" = "Part A | 免疫浸润相关图",
  "A5_MultiCohort_KM.png"      = "Part A | 5 队列 KM 曲线",
  "B1_RiskScore_tripanel.png"  = "Part B | 风险评分三联图",
  "B2_KM_training.png"         = "Part B | KM 生存曲线（训练队列）",
  "B3_TimeROC.png"             = "Part B | Time-ROC 曲线（1/3/5 年）",
  "B5_Calibration_rms.png"     = "Part B | 校准曲线（rms bootstrap）"
)

cat("\n  图表清单：\n")
for (fn in names(fig_list)) {
  status <- ifelse(file.exists(file.path(fig_dir, fn)), "[✓]", "[✗]")
  cat(sprintf("  %s %-35s %s\n", status, fn, fig_list[[fn]]))
}

# ════════════════════════════════════════════════════════════════════════════════
# LASSO 图：CV 曲线 + 系数轨迹图（Agent 风格 ggplot2）
# 数据来源：你本地的 GSE10846_sur_model.Rdata
# ════════════════════════════════════════════════════════════════════════════════

cat("LASSO 图...
")

library(glmnet)
library(ggplot2)
library(patchwork)

# ── 1. 读入原始数据 ───────────────────────────────────────────────────────────
rdata_path <- translate_path("D:/bulk-download/DLBCL_prognosis/GSE10846_sur_model.Rdata")
if (!file.exists(rdata_path)) stop("找不到 Rdata，请修改 rdata_path")
load(rdata_path)   # 加载 exprSet, meta

# ── 2. 准备 x / y（和你原来代码完全一致）────────────────────────────────────
# 用 Agent 确定的 13 个 LASSO 基因
lasso_genes <- c("AAK1","BASP1","CEP192","CHCHD10","ELMO1",
                 "IFIT3","MTR","NAE1","PPP2R5A","PPRC1",
                 "RASSF4","TIMM17A","USP12","NAGK", "MGAT1",
                 'FAM20A',"APOL3", "GSDMD","PICALM","NAP1L4",
                 "MRGBP","PPP1R18","CTSL","HMGB1","HAVCR2"
)

gene_use <- intersect(lasso_genes, rownames(exprSet))
cat("使用基因:", length(gene_use), "个
")

s <- intersect(rownames(meta), colnames(exprSet))
x <- t(exprSet[gene_use, s])
y <- ifelse(meta[s, "event"] == 0, 0, 1)

# ── 3. 重跑 LASSO（固定随机种子保证可重复）──────────────────────────────────
set.seed(1006)
cv_fit <- cv.glmnet(x = x, y = y, alpha = 1, family = "binomial",
                    nfolds = 10, type.measure = "deviance")
fit    <- glmnet(x = x, y = y, alpha = 1, family = "binomial")

cat("lambda.min =", cv_fit$lambda.min, "
")
cat("lambda.1se =", cv_fit$lambda.1se, "
")

# ── 4. CV 曲线（ggplot2 Agent 风格）─────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_fit$lambda),
  cvm        = cv_fit$cvm,
  cvup       = cv_fit$cvup,
  cvlo       = cv_fit$cvlo,
  nzero      = cv_fit$nzero
)

lmin_log <- log(cv_fit$lambda.min)
l1se_log <- log(cv_fit$lambda.1se)

p_cv <- ggplot(cv_df, aes(x = log_lambda, y = cvm)) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup),
              fill = "#AEC6E8", alpha = 0.4) +
  geom_line(color = "#2166AC", linewidth = 0.9) +
  geom_point(color = "#D73027", size = 1.8, shape = 16) +
  geom_vline(xintercept = lmin_log, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = l1se_log, linetype = "dashed",
             color = "#4DAC26", linewidth = 0.8) +
  annotate("text", x = lmin_log, y = max(cv_df$cvup),
           label = "lambda.min", hjust = -0.1, vjust = 1.2,
           size = 3.5, color = "#D73027", fontface = "bold") +
  annotate("text", x = l1se_log, y = max(cv_df$cvup),
           label = "lambda.1se", hjust = -0.1, vjust = 1.2,
           size = 3.5, color = "#4DAC26", fontface = "bold") +
  # 顶部显示非零系数个数
  scale_x_continuous(
    sec.axis = sec_axis(
      ~ .,
      breaks = log(cv_fit$lambda)[seq(1, length(cv_fit$lambda), length.out = 8)],
      labels = cv_df$nzero[seq(1, nrow(cv_df), length.out = 8)],
      name   = "Number of Non-zero Coefficients"
    )
  ) +
  labs(title    = "LASSO Cross-Validation",
       subtitle = paste0("10-fold CV  |  lambda.min = ",
                         round(cv_fit$lambda.min, 4),
                         "  |  Genes selected = ",
                         sum(coef(cv_fit, s = "lambda.min")[-1] != 0)),
       x = expression(log(lambda)),
       y = "Binomial Deviance") +
  theme_classic(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9.5, color = "grey40"),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
        axis.title.x.top = element_text(size = 10, color = "grey30"),
        axis.text.x.top  = element_text(size = 8,  color = "grey30"))

# ── 5. 系数轨迹图（ggplot2 Agent 风格）──────────────────────────────────────
# 提取所有 lambda 下的系数矩阵
coef_mat   <- as.matrix(coef(fit))[-1, ]   # 去掉截距行，genes × lambdas
coef_df    <- as.data.frame(t(coef_mat))
coef_df$log_lambda <- log(fit$lambda)

coef_long <- tidyr::pivot_longer(coef_df,
                                 cols      = -log_lambda,
                                 names_to  = "gene",
                                 values_to = "coef")

# 标记最终入选基因（lambda.min 下非零）
final_genes <- rownames(coef_mat)[coef_mat[, which.min(abs(fit$lambda - cv_fit$lambda.min))] != 0]

# 为入选基因准备标签位置（取最右侧 lambda 处的系数值）
label_df <- coef_long %>%
  filter(gene %in% final_genes) %>%
  group_by(gene) %>%
  slice_min(log_lambda, n = 1) %>%
  ungroup()

# 颜色：入选基因高亮，其余灰色
n_final  <- length(final_genes)
pal_cols <- colorRampPalette(c("#E63946","#F4A261","#2A9D8F",
                               "#457B9D","#9B2226","#6A4C93",
                               "#3A86FF","#FB5607","#8AC926",
                               "#FFBE0B","#FF006E","#3D405B",
                               "#06D6A0"))(n_final)
names(pal_cols) <- final_genes

coef_long$highlight <- ifelse(coef_long$gene %in% final_genes,
                              coef_long$gene, "other")

p_coef <- ggplot(coef_long, aes(x = log_lambda, y = coef,
                                group = gene, color = highlight)) +
  geom_line(data = subset(coef_long, highlight == "other"),
            color = "grey80", linewidth = 0.4, alpha = 0.6) +
  geom_line(data = subset(coef_long, highlight != "other"),
            aes(color = highlight), linewidth = 0.9) +
  geom_vline(xintercept = lmin_log, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = l1se_log, linetype = "dashed",
             color = "#4DAC26", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data        = label_df,
    aes(label   = gene, color = gene),
    size        = 3.2,
    fontface    = "bold",
    nudge_x     = -0.3,
    direction   = "y",
    segment.size = 0.3,
    segment.color = "grey60",
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c(pal_cols, "other" = "grey80"),
    guide  = "none"
  ) +
  labs(title    = "LASSO Coefficient Paths",
       subtitle = paste0("Regularization path  |  ",
                         n_final, " genes selected at lambda.min"),
       x = expression(log(lambda)),
       y = "Coefficient") +
  theme_classic(base_size = 12) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9.5, color = "grey40"),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.3))

# ── 6. 组合并保存 ─────────────────────────────────────────────────────────────
p_lasso <- p_cv + p_coef +
  plot_annotation(
    title    = "LASSO Penalized Logistic Regression — GSE10846 DLBCL",
    subtitle = paste0("Input: ", length(gene_use), " candidate genes  |  ",
                      "n = ", nrow(x), " samples"),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(fig_dir, "A0_LASSO_CV_coef.png"), p_lasso,
       width = 14, height = 6, dpi = 300)

# 也单独保存两张
ggsave(file.path(fig_dir, "A0a_LASSO_CV.png"),   p_cv,   width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "A0b_LASSO_coef.png"), p_coef, width = 7, height = 5.5, dpi = 300)

cat("  已保存 A0_LASSO_CV_coef.png / A0a / A0b
")
