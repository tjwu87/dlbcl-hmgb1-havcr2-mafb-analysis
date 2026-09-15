# >>> 自动注入：统一路径配置（由 tools/migrate_scripts.py 添加）<<<
# 注意：本仓库路径含中文，R 在 Windows 下对含中文的**绝对路径**做编码
# 转换会失败（unable to translate to UTF-8），因此统一使用相对路径，
# 并确保在仓库根目录运行。
source(file.path("config", "paths.R"))
# <<< 自动注入结束 >>>

# ==============================================================================
# 1. 环境准备与包加载
# ==============================================================================
# rm(list = ls())   # 迁移工具已注释：会清除注入的路径配置
options(timeout = 100000, scipen = 20)

library(GEOquery)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(Biobase)
library(hgu133plus2.db)
library(glmnet)
library(survival)
library(RColorBrewer)
library(patchwork)
library(fmsb)
library(ggpubr)
library(ggunchained)
library(tinyarray)

# 【全局配色：低风险=蓝，高风险=红，所有图统一】
my_colors <- c("#3288BD", "#D73027")

# ==============================================================================
# 2. GEO 数据下载与临床信息提取
# ==============================================================================
cat("Parsing GSE10846 series matrix...\n")
setwd(translate_path("D:/BadiduNetdiskDownload/R/Tcell"))
gse <- getGEO(filename=translate_path("D:/BadiduNetdiskDownload/R/GSE10846survive/GSE10846_series_matrix.txt.gz"),
              GSEMatrix=TRUE, getGPL=FALSE)
pdata <- pData(gse)

strip <- function(x, prefix) str_trim(str_replace(x, fixed(prefix), ""))

clin <- data.frame(
  sample_id       = pdata$geo_accession,
  age             = as.numeric(strip(pdata$`characteristics_ch1.1`, "Age: ")),
  status          = strip(pdata$`characteristics_ch1.7`,  "Clinical info: Follow up status: "),
  follow_up_years = as.numeric(strip(pdata$`characteristics_ch1.8`, "Clinical info: Follow up years: ")),
  stringsAsFactors = FALSE
)
clin$OS_event <- ifelse(clin$status == "DEAD", 1, ifelse(clin$status == "ALIVE", 0, NA))
clin$OS_years <- clin$follow_up_years

# ==============================================================================
# 3. 读取 Agent LASSO 系数，计算风险评分
#    【核心改动】：不重新跑 LASSO，直接读取已保存的13基因系数
# ==============================================================================
cat("读取 Agent LASSO 模型系数 (13基因)...\n")

formula_df  <- read.csv(translate_path("D:/bulk-download/DLBCL_prognosis/tables/LASSO_prognostic_formula.csv"),
                         stringsAsFactors = FALSE)
lasso_genes <- formula_df$Gene
lasso_coefs <- formula_df$Coefficient
cat("模型基因:", paste(lasso_genes, collapse=", "), "\n")

# 用 hgu133plus2.db 找每个基因的最优 probe（与 Agent 训练时一致）
expr_mat  <- exprs(gse)
probe_ids <- rownames(expr_mat)
gene_syms <- mapIds(hgu133plus2.db,
                    keys      = probe_ids,
                    column    = "SYMBOL",
                    keytype   = "PROBEID",
                    multiVals = "first")
probe_gene <- data.frame(probe = probe_ids,
                          gene  = as.character(gene_syms),
                          stringsAsFactors = FALSE)
probe_gene <- probe_gene[!is.na(probe_gene$gene), ]

get_best_probe <- function(g, em, pg) {
  pr <- pg$probe[pg$gene == g]
  pr <- pr[pr %in% rownames(em)]
  if (length(pr) == 0) return(NA)
  if (length(pr) == 1) return(pr)
  pr[which.max(rowMeans(em[pr, , drop=FALSE]))]
}

best_probe_vec <- sapply(lasso_genes, get_best_probe, em=expr_mat, pg=probe_gene)
found          <- !is.na(best_probe_vec)
cat("找到 probe 的基因:", sum(found), "/", length(lasso_genes), "\n")

# 构建 13基因 × 样本 表达矩阵
expr_13           <- expr_mat[best_probe_vec[found], , drop=FALSE]
rownames(expr_13) <- lasso_genes[found]

# 筛选有生存数据的样本
clin_all     <- clin %>% filter(!is.na(OS_years) & OS_years > 0 & !is.na(OS_event))
expr_all_ids <- intersect(clin_all$sample_id, colnames(expr_13))
clin_all     <- clin_all[clin_all$sample_id %in% expr_all_ids, ]
expr_all_t   <- t(expr_13[, clin_all$sample_id, drop=FALSE])  # samples × genes

# scale 全队列（与 Agent 训练时一致，保证风险评分可比）
X_all <- scale(expr_all_t)

# 计算风险评分并分组
risk_scores        <- as.numeric(X_all %*% lasso_coefs[found])
clin_all$RiskScore <- risk_scores
clin_all$RiskGroup <- factor(
  ifelse(risk_scores >= median(risk_scores), "High", "Low"),
  levels = c("Low", "High")
)
cat("High:", sum(clin_all$RiskGroup=="High"),
    " Low:", sum(clin_all$RiskGroup=="Low"), "\n")

# ==============================================================================
# 4. 免疫浸润 (CIBERSORT) 分析
# ==============================================================================

# ── 加载 CIBERSORT 结果 ──────────────────────────────────────
f = "GSE10846_ciber.Rdata"
if(!file.exists(f)){
  library(CIBERSORT)
  lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
  TME.results = cibersort(lm22f, "exp.txt", perm = 1000, QN = T)
  save(TME.results, file = f)
}
load(f)

re <- TME.results[, -(23:25)]   # 去掉最后3列统计列
common_samples <- intersect(rownames(re), clin_all$sample_id)
re_sub   <- re[common_samples, ]
clin_sub <- clin_all[clin_all$sample_id %in% common_samples, ]

# 整理长格式画图数据（按 RiskGroup 排序，使堆叠图分区清晰）
dat3 <- re_sub %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(group = clin_sub$RiskGroup) %>%
  gather(key = Cell_type, value = Proportion, -Sample, -group) %>%
  arrange(group)
dat3$Sample <- factor(dat3$Sample, ordered = T, levels = unique(dat3$Sample))

# 顶部颜色条数据（按分组排序）
dat4 <- data.frame(
  a     = 1:nrow(clin_sub),
  b     = 1,
  group = sort(clin_sub$RiskGroup)
)

# ══════════════════════════════════════════
# 4.1 堆叠柱状图 + 顶部颜色条 — 可调参数区
# ══════════════════════════════════════════
bar_legend_text_size  <- 12    # 图例文字大小
bar_legend_title_size <- 14    # 图例标题大小
bar_heights_ratio     <- c(1, 10)  # 顶部颜色条与主图的高度比例；颜色条太窄改大第1个数
bar_x_text_size       <- 20     # X轴样本名文字大小；设为0则隐藏
bar_x_text_angle      <- 60    # X轴文字旋转角度；常用45或90
bar_x_text_hjust      <- 1     # 水平对齐；角度90时设1，角度45时设1，角度0时设0.5

# ══════════════════════════════════════════
# 4.1 堆叠柱状图 + 顶部颜色条 — 主代码区
# ══════════════════════════════════════════
mypalette <- colorRampPalette(brewer.pal(8, "Set1"))

p1 <- ggplot(dat4, aes(x = a, y = b)) +
  geom_tile(aes(fill = group)) +
  scale_fill_manual(values = my_colors, labels = c("Low Risk", "High Risk")) +
  theme_void() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")

p2 <- ggplot(dat3, aes(x = Sample, y = Proportion, fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type", x = "", y = "Estimated Proportion") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size  = bar_x_text_size,
                                angle = bar_x_text_angle,
                                hjust = bar_x_text_hjust),
    axis.ticks.x = element_line()
  ) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = mypalette(22))

patch_plot <- p1 / p2 +
  plot_layout(heights = bar_heights_ratio, guides = "collect") &
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = bar_legend_text_size),
        legend.title     = element_text(size = bar_legend_title_size, face = "bold"))

print(patch_plot)


# ══════════════════════════════════════════
# 4.2 箱线图 — 可调参数区
# ══════════════════════════════════════════

# ── X 轴 ──────────────────────────────────
box_x_text_size  <- 60    # X 轴标签文字大小
box_x_angle      <- 70    # X 轴标签旋转角度（标签太长可改为 90）
box_x_vjust      <- 1     # X 轴标签垂直对齐
box_x_hjust      <- 1     # X 轴标签水平对齐
box_x_title      <- "Cell Type"          # X 轴标题文字
box_x_title_size <- 40    # X 轴标题文字大小

# ── Y 轴 ──────────────────────────────────
box_y_text_size  <- 40    # Y 轴刻度数字大小
box_y_title      <- "Estimated Proportion"  # Y 轴标题文字
box_y_title_size <- 40    # Y 轴标题文字大小

# ── 箱体 ──────────────────────────────────
box_line_width    <- 2    # 箱体边框线粗细
box_whisker_width <- 2    # 须线粗细（通常与 box_line_width 保持一致）
box_width         <- 0.6  # 箱体宽度（0~1，越大越宽；原来设的 2 会超出范围）

# ── 离群点 ────────────────────────────────
box_outlier_size  <- 2    # 离群点大小
box_outlier_shape <- 16   # 离群点形状（16=实心圆，4=叉，NA=隐藏）

# ── 显著性标注 ────────────────────────────
box_signif_size   <- 20    # 星号文字大小（默认约 3.88）

# ── 颜色 ──────────────────────────────────
box_colors <- c("#3288BD", "#D73027")  # 低风险、高风险填充色

# ── 图例 ──────────────────────────────────
box_legend_title      <- "Risk Group"              # 图例标题
box_legend_labels     <- c("Low Risk", "High Risk") # 图例标签（顺序与 box_colors 对应）
box_legend_title_size <- 40    # 图例标题文字大小
box_legend_text_size  <- 40    # 图例标签文字大小
box_legend_position   <- "top" # 图例位置："top","bottom","left","right","none"

# ══════════════════════════════════════════
# 4.2 箱线图 — 主代码区
# ══════════════════════════════════════════

# 去掉全为 0 的细胞类型列（无意义）
k     <- colSums(re_sub) > 0
re    <- re_sub[, k]
Group <- clin_sub$RiskGroup   # draw_boxplot 需要此变量名

# 生成基础图对象
p_box <- draw_boxplot(t(re), factor(Group),
                      drop  = TRUE,
                      color = box_colors)

# 查看所有图层类型（调试用，确认图层编号后可注释掉）
for (i in seq_along(p_box$layers)) {
  cat(i, ":", class(p_box$layers[[i]]$geom)[1], "\n")
}

# 直接修改已有图层参数，避免叠加第二层箱线图
for (i in seq_along(p_box$layers)) {

  # 箱体图层
  if (inherits(p_box$layers[[i]]$geom, "GeomBoxplot")) {
    p_box$layers[[i]]$aes_params$linewidth      <- box_line_width
    p_box$layers[[i]]$aes_params$width          <- box_width
    p_box$layers[[i]]$aes_params$outlier.size   <- box_outlier_size
    p_box$layers[[i]]$geom_params$outlier.shape <- box_outlier_shape
  }

  # 显著性标注图层
  if (inherits(p_box$layers[[i]]$geom, "GeomText")) {
    p_box$layers[[i]]$aes_params$textsize <- box_signif_size
    p_box$layers[[i]]$aes_params$size     <- box_signif_size  # 兼容不同版本 ggsignif
  }

}

# 叠加样式调整
p_box <- p_box +
  scale_color_manual(values = box_colors,
                     labels = box_legend_labels) +
  scale_fill_manual(values  = box_colors,
                    labels  = box_legend_labels) +
  labs(x     = box_x_title,
       y     = box_y_title,
       color = box_legend_title,
       fill  = box_legend_title) +
  theme(
    # X 轴
    axis.text.x  = element_text(size  = box_x_text_size,
                                angle = box_x_angle,
                                vjust = box_x_vjust,
                                hjust = box_x_hjust),
    axis.title.x = element_text(size  = box_x_title_size,
                                face  = "bold"),
    # Y 轴
    axis.text.y  = element_text(size  = box_y_text_size),
    axis.title.y = element_text(size  = box_y_title_size,
                                face  = "bold"),
    # 图例
    legend.position = box_legend_position,
    legend.title    = element_text(size = box_legend_title_size,
                                   face = "bold"),
    legend.text     = element_text(size = box_legend_text_size)
  )

print(p_box)


# ══════════════════════════════════════════
# 4.3 雷达图 — 可调参数区
# ══════════════════════════════════════════

# 1. 雷达图范围控制
radar_expand <- 1.1     # 最大值放大倍数：图太挤改大(1.3/1.4)，太空改小(1.1)

# 2. 颜色与透明度
fill_alpha   <- 0.3     # 填充透明度：0.1-0.2 更透明，0.4-0.5 更实
line_width   <- 5       # 轮廓线宽：2 细，4-5 粗
grid_color   <- "grey80"  # 网格线颜色
grid_width   <- 2.5     # 网格线宽
label_size   <- 4       # 细胞类型标签字号：重叠时调小(1.5)，图大时调大(2.5)

# 3. 图例位置与大小
legend_x         <- 0.5   # 图例 x 位置：越大越往右
legend_y         <- -0.8  # 图例 y 位置：越大越往上
legend_cex       <- 5     # 图例文字大小
legend_pt_cex    <- 6     # 图例点大小
legend_spacing   <- 0.2   # 图例行距：越大越疏
legend_x_intersp <- 0.3   # 图例点与文字的水平间距

# 4. 显著性星号阈值
sig_high <- 0.001   # *** 阈值
sig_mid  <- 0.01    # **  阈值
sig_low  <- 0.05    # *   阈值

# 5. 图例文字
radar_legend_labels <- c("Low Risk", "High Risk")  # 可直接改组名

# ══════════════════════════════════════════
# 4.3 雷达图 — 主代码区
# ══════════════════════════════════════════

p_vals <- dat3 %>%
  group_by(Cell_type) %>%
  dplyr::summarise(
    p_value = wilcox.test(Proportion ~ group)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    sig_star = case_when(
      p_value < sig_high ~ "***",
      p_value < sig_mid  ~ "**",
      p_value < sig_low  ~ "*",
      TRUE ~ ""
    ),
    New_Cell_type = paste0(Cell_type, sig_star)
  )

su_mean <- dat3 %>%
  group_by(group, Cell_type) %>%
  dplyr::summarise(Proportion = mean(Proportion), .groups = 'drop') %>%
  left_join(p_vals %>% dplyr::select(Cell_type, New_Cell_type), by = "Cell_type")

dat_radar <- su_mean %>%
  dplyr::select(group, New_Cell_type, Proportion) %>%
  pivot_wider(names_from = New_Cell_type, values_from = Proportion) %>%
  column_to_rownames("group")

dat_radar <- rbind(
  Max = apply(dat_radar, 2, max) * radar_expand,
  Min = rep(0, ncol(dat_radar)),
  dat_radar
)

colors_fill <- c(alpha(my_colors[1], fill_alpha),
                 alpha(my_colors[2], fill_alpha))

radarchart(dat_radar,
           pcol = my_colors, pfcol = colors_fill, plwd = line_width,
           cglcol = grid_color, cglty = 1, cglwd = grid_width,
           axislabcol = "grey20", vlcex = label_size)

legend(x = legend_x, y = legend_y,
       legend = radar_legend_labels,
       bty = "n", pch = 20, col = my_colors, text.col = "black",
       cex = legend_cex, pt.cex = legend_pt_cex,
       y.intersp = legend_spacing, x.intersp = legend_x_intersp)

# ══════════════════════════════════════════
# 4.4 CIBERSORT 分割小提琴图 — 可调参数区
# ══════════════════════════════════════════

ciber_point_size      <- 5     # 均值点大小：越大越明显
ciber_errorbar_width  <- 0     # 误差线横杠宽度：0 = 无横杠，0.2 = 有横杠
ciber_errorbar_lwd    <- 1     # 误差线粗细
ciber_dodge_width     <- 0.5   # 两组点/误差线左右分开距离
ciber_stat_size       <- 15    # 显著性星号字体大小：调大改这里
ciber_x_angle         <- 70    # x 轴标签旋转角度
ciber_x_hjust         <- 1     # x 轴标签水平对齐
ciber_x_size          <- 25    # x 轴标签字体大小
ciber_legend_text     <- 25    # 图例文字大小
ciber_legend_title    <- 27    # 图例标题大小
ciber_legend_key      <- 2     # 图例色块大小（单位 cm）
ciber_legend_title_str <- "Group"  # 图例标题文字

# ══════════════════════════════════════════
# 4.4 CIBERSORT 分割小提琴图 — 主代码区
# ══════════════════════════════════════════

su_ciber <- dat3 %>%
  group_by(group, Cell_type) %>%
  summarise(s = sd(Proportion), Proportion = mean(Proportion), .groups = "drop")

p_ciber_violin <- ggplot(dat3, aes(x = Cell_type, y = Proportion, fill = group)) +
  geom_split_violin(color = NA, scale = "width") +
  scale_fill_manual(values = my_colors, name = ciber_legend_title_str) +
  theme_bw() +
  geom_point(data = su_ciber, aes(x = Cell_type, y = Proportion),
             pch = 19, position = position_dodge(ciber_dodge_width),
             size = ciber_point_size) +
  geom_errorbar(data = su_ciber,
                mapping = aes(ymin = Proportion - s, ymax = Proportion + s),
                width = ciber_errorbar_width,
                position = position_dodge(ciber_dodge_width),
                color = "black", linewidth = ciber_errorbar_lwd) +
  stat_compare_means(aes(group = group), label = "p.signif",
                     size = ciber_stat_size) +
  theme(
    axis.text.x     = element_text(angle = ciber_x_angle,
                                   hjust = ciber_x_hjust,
                                   size  = ciber_x_size),
    legend.text      = element_text(size = ciber_legend_text),
    legend.title     = element_text(size = ciber_legend_title),
    legend.key.size  = unit(ciber_legend_key, "cm")
  )
print(p_ciber_violin)

# ==============================================================================
# 5. 免疫检查点 (ICP) 分析
# ==============================================================================

# ══════════════════════════════════════════
# 5.1 ICP 分割小提琴图 — 可调参数区
# ══════════════════════════════════════════

# 1. 显著性筛选阈值
icp_sig_cutoff <- 0.05    # 只展示 p < 此值的基因；改大会显示更多基因

# 2. 分组显示名称
icp_legend_title  <- "Group"
icp_legend_labels <- c("Low Risk", "High Risk")

# 3. 汇总点和误差线参数
icp_point_size    <- 5    # 均值点大小
icp_dodge_width   <- 1    # 两组左右分开距离：改大分得更开
icp_errorbar_width <- 0   # 误差线横杠宽度：0 = 无横杠
icp_errorbar_lwd  <- 2    # 误差线粗细

# 4. 统计检验字体参数
icp_stat_size <- 24       # 显著性星号字体大小：调大改这里

# 5. 坐标轴字体参数
icp_x_angle  <- 45        # x 轴标签旋转角度
icp_x_hjust  <- 1         # x 轴标签水平对齐
icp_x_face   <- "italic"  # x 轴标签字体样式：italic=斜体，plain=正常，bold=加粗
icp_x_color  <- "black"   # x 轴标签颜色
icp_x_size   <- 60        # x 轴标签字体大小：图小时调小（如 20-30）
icp_y_color  <- "black"   # y 轴刻度颜色
icp_y_size   <- 60        # y 轴刻度字体大小
icp_axis_title_size  <- 60   # 坐标轴标题字体大小
icp_axis_title_color <- "black"

# 6. 图例字体参数
icp_legend_title_size <- 48  # 图例标题字体大小
icp_legend_text_size  <- 44  # 图例文字字体大小

# 7. 坐标轴标题文字
icp_x_label <- "Immune Checkpoints"
icp_y_label <- "Expression Level"

# ══════════════════════════════════════════
# 5.1 ICP 分割小提琴图 — 主代码区
# ══════════════════════════════════════════

tmp       <- readLines("immu_check_point.txt")
icp_genes <- str_split(tmp, ", ")[[1]]

exprSet   <- trans_array(expr_mat, toTable(hgu133plus2SYMBOL))
icp_genes <- icp_genes[icp_genes %in% rownames(exprSet)]

icp_expr        <- as.data.frame(t(exprSet[icp_genes, common_samples]))
icp_expr$Sample <- rownames(icp_expr)
icp_expr$group  <- clin_sub$RiskGroup

icp_long <- icp_expr %>%
  gather(key = Gene, value = Expression, -Sample, -group)

# 筛选有显著差异的基因
p_vals_icp <- icp_long %>%
  group_by(Gene) %>%
  dplyr::summarise(p_value = wilcox.test(Expression ~ group)$p.value,
                   .groups = "drop")
sig_genes    <- p_vals_icp %>% filter(p_value < icp_sig_cutoff) %>% pull(Gene)
icp_long_sig <- icp_long %>% filter(Gene %in% sig_genes)

su_icp_sig <- icp_long_sig %>%
  group_by(group, Gene) %>%
  dplyr::summarise(s = sd(Expression, na.rm=TRUE),
                   Expression = mean(Expression, na.rm=TRUE),
                   .groups = "drop")

p_icp <- ggplot(icp_long_sig, aes(x = Gene, y = Expression, fill = group)) +
  geom_split_violin(color = NA, scale = "width") +
  guides(fill = guide_legend(title = icp_legend_title)) +
  scale_fill_manual(values = my_colors, labels = icp_legend_labels) +
  theme_bw() +
  geom_point(data = su_icp_sig, aes(x = Gene, y = Expression),
             pch = 19, position = position_dodge(icp_dodge_width),
             size = icp_point_size) +
  geom_errorbar(data = su_icp_sig,
                mapping = aes(ymin = Expression - s, ymax = Expression + s),
                width = icp_errorbar_width,
                position = position_dodge(icp_dodge_width),
                color = "black", linewidth = icp_errorbar_lwd) +
  stat_compare_means(aes(group = group), label = "p.signif",
                     size = icp_stat_size) +
  labs(x = icp_x_label, y = icp_y_label) +
  theme(
    axis.text.x  = element_text(angle = icp_x_angle, hjust = icp_x_hjust,
                                face  = icp_x_face,  color = icp_x_color,
                                size  = icp_x_size),
    axis.text.y  = element_text(color = icp_y_color, size = icp_y_size),
    axis.title   = element_text(size = icp_axis_title_size,
                                color = icp_axis_title_color),
    legend.title = element_text(size = icp_legend_title_size),
    legend.text  = element_text(size = icp_legend_text_size)
  )
print(p_icp)

# ══════════════════════════════════════════
# 5.2 ICP 箱线图 — 可调参数区
# ══════════════════════════════════════════
# 注意：此图展示所有 ICP 基因（不做显著性筛选）
icp_box_drop <- TRUE   # TRUE = 去掉全为0的基因；FALSE = 保留所有基因

# ══════════════════════════════════════════
# 5.2 ICP 箱线图 — 主代码区
# ══════════════════════════════════════════
draw_boxplot(exprSet[icp_genes, common_samples],
             clin_sub$RiskGroup,
             drop  = icp_box_drop,
             color = my_colors) +
  xlab("Immune checkpoint")

# ==============================================================================
# 6. TIDE 免疫治疗响应预测分析
# ==============================================================================

# ══════════════════════════════════════════
# TIDE 全局 — 可调参数区
# ══════════════════════════════════════════

# 1. 数据文件路径
tide_file <- "DLBCL_MRGs_Cluster_PD1.csv"

# 2. 配色
tide_risk_colors <- c("Low Risk"  = "#A6C2DE",
                      "High Risk" = "#EFACA9")
tide_resp_colors <- c("False" = "#D1AFD6",
                      "True"  = "#AFDA9D")

# 3. 图例字体（三张图共用）
tide_legend_title_size <- 40
tide_legend_text_size  <- 40

# ══════════════════════════════════════════
# 图A：ICB 响应率堆叠柱状图 — 可调参数区
# ══════════════════════════════════════════

bar_width_tide        <- 0.6    # 柱子宽度
bar_label_size_tide   <- 30     # 柱内百分比字体大小
bar_axis_text_size    <- 40     # 坐标轴刻度字体大小
bar_axis_title_size_A <- 40     # 坐标轴标题字体大小（X轴和Y轴共用）
bar_title_A           <- "ICB Response Rate"  # 图A标题文字
bar_title_size_A      <- 40     # 图A标题字体大小
bar_x_title_A         <- "Risk Group"         # X轴标题文字
bar_y_title_A         <- "Percentage"         # Y轴标题文字

# ══════════════════════════════════════════
# 图B：TIDE 得分小提琴图 — 可调参数区
# ══════════════════════════════════════════

tide_violin_alpha     <- 0.7    # 小提琴填充透明度
tide_boxplot_width    <- 0.15   # 内嵌箱线图宽度
tide_stat_size_B      <- 20     # 显著性星号字体大小
tide_stat_label_x     <- 1.4    # 显著性标签 x 位置
tide_strip_size       <- 40     # 分面标题字体大小
tide_x_text_size      <- 40     # X轴刻度字体大小
tide_title_B          <- "TIDE Component Scores"  # 图B标题文字
tide_title_size_B     <- 40     # 图B标题字体大小
tide_y_title_B        <- "Score"                  # 图B Y轴标题文字
tide_y_title_size_B   <- 40     # 图B Y轴标题字体大小

# ══════════════════════════════════════════
# 图C：TIDE 瀑布图 — 可调参数区
# ══════════════════════════════════════════

waterfall_bar_width    <- 1     # 柱子宽度
waterfall_hline_lwd    <- 0.6   # 零线粗细
waterfall_x_title      <- "Patients (Ordered by Risk Group and TIDE Score)"  # X轴标题文字
waterfall_x_title_size <- 40    # X轴标题字体大小
waterfall_y_title      <- "TIDE Score"  # Y轴标题文字
waterfall_y_text_size  <- 40    # Y轴刻度字体大小
waterfall_y_title_size <- 40    # Y轴标题字体大小

# ══════════════════════════════════════════
# 拼图比例 — 可调参数区
# ══════════════════════════════════════════

tide_layout_heights <- c(1, 1.5)
tide_strip_ratio    <- c(0.05, 1)

# ══════════════════════════════════════════
# TIDE 分析 — 主代码区（通常不需要改）
# ══════════════════════════════════════════

tide_res <- read.csv(tide_file, header=TRUE, row.names=1,
                     check.names=FALSE, stringsAsFactors=FALSE)

common_tide <- intersect(rownames(tide_res), clin_all$sample_id)
tide_res    <- tide_res[common_tide, ]
temp_risk   <- clin_all$RiskGroup[match(common_tide, clin_all$sample_id)]
tide_res$RiskGroup <- factor(temp_risk,
                             levels = c("Low",  "High"),
                             labels = c("Low Risk", "High Risk"))

my_colors <- tide_risk_colors

# --- 图A：响应率堆叠柱状图 ---
resp_summary <- tide_res %>%
  group_by(RiskGroup, Responder) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(RiskGroup) %>%
  mutate(Freq  = Count / sum(Count),
         Label = scales::percent(Freq, accuracy = 1))

p_response_bar <- ggplot(resp_summary,
                         aes(x = RiskGroup, y = Freq, fill = Responder)) +
  geom_bar(stat = "identity", position = "fill", width = bar_width_tide) +
  geom_text(aes(label = Label),
            position = position_fill(vjust = 0.5),
            size = bar_label_size_tide, color = "white", fontface = "bold") +
  scale_fill_manual(values = tide_resp_colors, name = "Response") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(title = bar_title_A,
       x     = bar_x_title_A,
       y     = bar_y_title_A) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5,
                                size = bar_title_size_A),
    axis.text    = element_text(color = "black", size = bar_axis_text_size),
    axis.title   = element_text(face = "bold", size = bar_axis_title_size_A),
    legend.title = element_text(size = tide_legend_title_size),
    legend.text  = element_text(size = tide_legend_text_size)
  )

# --- 图B：TIDE 三大核心得分小提琴图 ---
tide_melt <- tide_res %>%
  dplyr::select(RiskGroup, TIDE, Dysfunction, Exclusion) %>%
  pivot_longer(cols = c(TIDE, Dysfunction, Exclusion),
               names_to = "ScoreType", values_to = "Score")

tide_melt$ScoreType <- factor(tide_melt$ScoreType,
                              levels = c("TIDE", "Dysfunction", "Exclusion"))

p_scores_violin <- ggplot(tide_melt,
                          aes(x = RiskGroup, y = Score, fill = RiskGroup)) +
  geom_violin(alpha = tide_violin_alpha, trim = FALSE, color = NA) +
  geom_boxplot(width = tide_boxplot_width, fill = "white",
               outlier.shape = NA, color = "grey30") +
  facet_wrap(~ScoreType, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x = tide_stat_label_x, size = tide_stat_size_B) +
  scale_fill_manual(values = my_colors, name = "Risk Group") +
  theme_bw() +
  labs(title = tide_title_B,
       x     = NULL,
       y     = tide_y_title_B) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5,
                                    size = tide_title_size_B),
    strip.background = element_rect(fill = "grey90", color = "white"),
    strip.text       = element_text(size = tide_strip_size, face = "bold"),
    axis.text.x      = element_text(size = tide_x_text_size, color = "black"),
    axis.title.y     = element_text(face = "bold", size = tide_y_title_size_B),
    legend.title     = element_text(size = tide_legend_title_size),
    legend.text      = element_text(size = tide_legend_text_size)
  )

# --- 图C：瀑布图 ---
plot_tide <- tide_res %>%
  arrange(RiskGroup, desc(TIDE)) %>%
  mutate(x_order = row_number())

p_top <- ggplot(plot_tide, aes(x = x_order, y = 1, fill = RiskGroup)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = my_colors, guide = "none") +
  theme_void() +
  scale_x_continuous(expand = c(0.01, 0.01))

p_main <- ggplot(plot_tide, aes(x = x_order, y = TIDE, fill = Responder)) +
  geom_bar(stat = "identity", width = waterfall_bar_width) +
  geom_hline(yintercept = 0, color = "grey30", linetype = "dashed",
             linewidth = waterfall_hline_lwd) +
  scale_fill_manual(values = tide_resp_colors, guide = "none") +
  labs(x = waterfall_x_title,
       y = waterfall_y_title) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  theme_classic() +
  theme(
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = waterfall_x_title_size, face = "bold",
                                margin = margin(t = 10)),
    axis.text.y  = element_text(size = waterfall_y_text_size, color = "black"),
    axis.title.y = element_text(size = waterfall_y_title_size, face = "bold")
  )

# --- 拼图 ---
layout_design <- (
  (p_response_bar | p_scores_violin) /
  (p_top / p_main + plot_layout(heights = tide_strip_ratio))
) +
  plot_layout(heights = tide_layout_heights, guides = "collect") &
  theme(
    legend.position = "right",
    legend.title    = element_text(size = tide_legend_title_size),
    legend.text     = element_text(size = tide_legend_text_size)
  )

print(layout_design)

