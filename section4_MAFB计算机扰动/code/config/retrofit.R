# =============================================================================
# 出图改造层（R 侧，路线 B·修订版）—— 与 config/retrofit.py 同原理
# -----------------------------------------------------------------------------
# 作者提出的三点：
#   1. 改用矢量图
#   2. 按最终物理尺寸设计字号
#   3. 字与图的比例同步放大/缩小（字号与图形元素放大相同倍数，不失调）
#
# 做法：**画布保持脚本原始尺寸不缩放**，仅按「最终显示宽度」反推字号放大倍数：
#
#   s      = 最终显示宽 / 画布宽            （显示缩放）
#   factor = (目标字号 / 基准字号 10) / s   （字号与元素同步放大倍数）
#
# 旧版会把画布缩到 3.2 in 并把字号放大数倍，导致"字非常大、图很小"，已废弃。
#
# 用法：
#   Rscript tools/run_with_retrofit.R 07_bulk_prognosis/07_01_prognosis_main.R 6.3 7
# =============================================================================

if (!exists("RETRO_W"))  RETRO_W  <- 6.3    # 最终显示宽度（= \textwidth）
if (!exists("RETRO_H"))  RETRO_H  <- 7.6    # 最终显示高度上限
if (!exists("RETRO_PT")) RETRO_PT <- 7      # 排版后希望显示的字号
if (!exists("RETRO_OUT")) RETRO_OUT <- ""
if (!exists("RETRO_VECTOR_ONLY")) RETRO_VECTOR_ONLY <- TRUE   # 只出矢量 PDF

BASE_FONT_PT <- 10

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggsave <- function(filename, plot = ggplot2::last_plot(),
                     width = NA, height = NA, dpi = 300, ...) {
    # 画布尺寸保持原样，只算显示缩放
    s <- 1
    if (!is.na(width) && !is.na(height) && is.numeric(width) && is.numeric(height)) {
      s <- min(RETRO_W / width, RETRO_H / height, 1)
    }
    if (s < 1) {
      factor <- (RETRO_PT / BASE_FONT_PT) / s
      plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = BASE_FONT_PT * factor))
    }
    dpi <- max(dpi, 600, na.rm = TRUE)
    if (nzchar(RETRO_OUT)) {
      filename <- file.path(RETRO_OUT, basename(filename))
      dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    }
    # ① 矢量 PDF —— 投稿用，缩放不失真（首选）
    stem <- sub("\\.(png|tiff?|jpg|jpeg)$", "", filename, ignore.case = TRUE)
    pdf_file <- paste0(stem, ".pdf")
    message(sprintf("[retrofit] %s  画布 %.2f x %.2f in，显示缩放 %.3f，字号 x%.2f",
                    basename(pdf_file), width, height, s,
                    if (s < 1) (RETRO_PT / BASE_FONT_PT) / s else 1))
    ggplot2::ggsave(pdf_file, plot, width = width, height = height,
                    device = grDevices::cairo_pdf, bg = "white")
    # ② 原格式（多为 PNG）仅作快速预览，体积大且缩放会发糊，不用于投稿
    if (!isTRUE(RETRO_VECTOR_ONLY)) {
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
    }
    invisible(NULL)
  }
}

# base graphics 的 png()：改走矢量 PDF 设备。
# 这些图（如列线图、KM 曲线）本身就是矢量绘制的，只是被写成了位图；
# 换成 cairo_pdf 后同样清晰，且缩放不失真。
png <- function(filename = "Rplot%03d.png", width = 480, height = 480,
                res = 72, ...) {
  if (nzchar(RETRO_OUT)) {
    filename <- file.path(RETRO_OUT, basename(filename))
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  }
  pdf_file <- sub("\\.(png|tiff?|jpg|jpeg)$", ".pdf", filename, ignore.case = TRUE)
  if (isTRUE(RETRO_VECTOR_ONLY)) {
    message(sprintf("[retrofit] %s  %.2f x %.2f in (矢量)",
                    basename(pdf_file), width / res, height / res))
    return(grDevices::cairo_pdf(pdf_file, width = width / res, height = height / res))
  }
  res <- max(res, 600, na.rm = TRUE)
  grDevices::png(filename, width = width, height = height, res = res, ...)
}
