# =============================================================================
# 统一出图规范（R 侧，ggplot2）—— 与 config/plot_style.py 规格一致
# -----------------------------------------------------------------------------
# 版面：\textwidth = 160 mm = 6.30 in，\textheight = 216 mm = 8.50 in
# 字号：以 PDF 中的真实 pt 为准（8 pt 起，不小于 6 pt）
#
# 用法：
#   source("config/plot_style.R")
#   p <- ggplot(...) + theme_paper()
#   save_figure(p, file.path(FIG_DIR, "Fig6A_KM"), width_in = 3.1, height_in = 2.6)
#
# 旧脚本的问题与 Python 侧相同：画布开得过大（13×13 in 等）而字号用默认值，
# 排版缩放后字号只剩 1~2 pt。这里一律按最终版面尺寸出图。
# =============================================================================

TEXTWIDTH_IN  <- 160 / 25.4   # 6.30
TEXTHEIGHT_IN <- 216 / 25.4   # 8.50
SAFE_HEIGHT_IN <- 7.60
BASE_FONT_PT  <- 8
MIN_FONT_PT   <- 6
DEFAULT_DPI   <- 600

theme_paper <- function(base_size = BASE_FONT_PT, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base_size + 1, hjust = 0.5),
      axis.title       = ggplot2::element_text(size = base_size + 0.5),
      axis.text        = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(linewidth = 0.4, colour = "black"),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 0.5),
      legend.background = ggplot2::element_blank(),
      legend.key       = ggplot2::element_blank(),
      legend.key.size  = grid::unit(8, "pt"),
      strip.text       = ggplot2::element_text(size = base_size),
      strip.background = ggplot2::element_blank(),
      plot.margin      = grid::unit(c(2, 2, 2, 2), "pt"),
      complete = TRUE
    )
}

# 按 panel 网格推算画布尺寸（与 Python panel_figsize 同逻辑）
panel_size <- function(cols, rows, panel_aspect = 1.0, width_frac = 1.0) {
  total_w <- TEXTWIDTH_IN * width_frac
  gap_x <- (BASE_FONT_PT * 0.55) / 72 * (cols - 1)
  gap_y <- (BASE_FONT_PT * 0.65) / 72 * (rows - 1)
  panel_w <- (total_w - gap_x) / cols
  panel_h <- panel_w / panel_aspect
  c(width = round(total_w, 3),
    height = round(panel_h * rows + gap_y + 0.25, 3))
}

# 统一导出：PDF（矢量，投稿用）+ PNG（600 dpi，单独上传用）
save_figure <- function(plot, filename, width_in, height_in,
                        dpi = DEFAULT_DPI, formats = c("pdf", "png"),
                        verbose = TRUE) {
  stem <- sub("\\.(png|pdf|svg|tiff?|eps)$", "", filename)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  for (fmt in formats) {
    out <- paste0(stem, ".", fmt)
    if (fmt == "pdf") {
      ggplot2::ggsave(out, plot, width = width_in, height = height_in,
                      units = "in", device = grDevices::cairo_pdf,
                      useDingbats = FALSE)
    } else {
      ggplot2::ggsave(out, plot, width = width_in, height = height_in,
                      units = "in", dpi = dpi, type = "cairo", bg = "white")
    }
    if (verbose) message(sprintf("[saved] %s (%.2f x %.2f in)", basename(out), width_in, height_in))
  }
  invisible(NULL)
}

# 出图前自检
check_size <- function(width_in, height_in, name = "") {
  msgs <- character(0)
  if (width_in > TEXTWIDTH_IN * 1.02)
    msgs <- c(msgs, sprintf("宽度 %.2fin 超过 textwidth(%.2fin)", width_in, TEXTWIDTH_IN))
  if (height_in > SAFE_HEIGHT_IN)
    msgs <- c(msgs, sprintf("高度 %.2fin 超过安全高度 %.2fin，建议拆页", height_in, SAFE_HEIGHT_IN))
  if (length(msgs)) message("[plot_style] ", name, " ", paste(msgs, collapse = "; "))
  invisible(msgs)
}
