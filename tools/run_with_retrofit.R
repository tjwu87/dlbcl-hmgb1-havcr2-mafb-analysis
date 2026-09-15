#!/usr/bin/env Rscript
# =============================================================================
# 以出图改造层运行既有 R 画图脚本（路线 B 的 R 侧入口）
#
#   Rscript tools/run_with_retrofit.R <script.R> [width_in] [font_pt]
#
# 示例：
#   Rscript tools/run_with_retrofit.R 07_bulk_prognosis/07_01_prognosis_main.R 3.2 7
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# 关键：本仓库路径含中文。R 在 Windows 下默认的 native encoding 无法转换
# 中文路径（unable to translate to native encoding / file name conversion
# problem），必须先切到 UTF-8 locale，否则读写数据文件都会失败。
invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))

if (length(args) < 1) {
  cat("用法: Rscript tools/run_with_retrofit.R <script.R> [width_in] [font_pt]\n")
  quit(status = 1)
}

# 解析 --out <目录>（位置参数顺序：script [width] [font]）
out_flag <- NULL
pos <- character(0)
i <- 1
while (i <= length(args)) {
  if (args[i] == "--out" && i < length(args)) {
    out_flag <- args[i + 1]
    i <- i + 2
  } else {
    pos <- c(pos, args[i])
    i <- i + 1
  }
}
args <- pos

script <- args[1]
if (!file.exists(script)) { cat("脚本不存在:", script, "\n"); quit(status = 1) }

# 注意：本仓库路径含中文（...\论文修改\...），R 在 Windows 下对含中文的
# **绝对路径** 做编码转换会失败（unable to translate to UTF-8）。
# 因此这里一律使用相对路径，并确保在仓库根目录运行。
RETRO_W   <- if (length(args) >= 2) as.numeric(args[2]) else 3.2
RETRO_PT  <- if (length(args) >= 3) as.numeric(args[3]) else 7
RETRO_OUT <- if (!is.null(out_flag)) out_flag else
  file.path("results", "figures",
            tools::file_path_sans_ext(basename(script)))

# 先定义变量，再 source retrofit.R（其中的函数会引用这些变量）
source(file.path("config", "retrofit.R"), local = FALSE)

cat(sprintf("[retrofit] 目标 panel 宽 %.2f in，排版后字号 %.0f pt\n", RETRO_W, RETRO_PT))
cat(sprintf("[retrofit] 输出重定向到 %s\n", RETRO_OUT))
cat(sprintf("[retrofit] 运行 %s\n", script))

source(script, local = FALSE)
cat("[retrofit] 完成\n")
