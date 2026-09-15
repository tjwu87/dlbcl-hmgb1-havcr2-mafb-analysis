# =============================================================================
# 统一路径配置（R 侧）—— 与 config/paths.py 保持同一套目录约定
# -----------------------------------------------------------------------------
# 使用方法（脚本头部已由 tools/migrate_scripts.py 自动注入）：
#   if (Sys.getenv("DLBCL_REPO_ROOT") == "") Sys.setenv(DLBCL_REPO_ROOT = getwd())
#   source(file.path(Sys.getenv("DLBCL_REPO_ROOT"), "config/paths.R"))
#
# 切换数据位置（不改动任何代码）：
#   Sys.setenv(DLBCL_DATA_ROOT = "D:/bulk-download")
#   export DLBCL_DATA_ROOT=/mnt/results        # Linux
# =============================================================================

REPO_ROOT <- Sys.getenv("DLBCL_REPO_ROOT")
if (!nzchar(REPO_ROOT)) REPO_ROOT <- getwd()

DATA_ROOT <- Sys.getenv("DLBCL_DATA_ROOT")
if (!nzchar(DATA_ROOT)) DATA_ROOT <- file.path(REPO_ROOT, "data")

# ── 数据集目录（与论文结果部分一一对应）──────────────────────────────────────
GSE182434   <- file.path(DATA_ROOT, "GSE182434")
CELLORACLE  <- file.path(DATA_ROOT, "celloracle0331")
GSE232853   <- file.path(DATA_ROOT, "GSE232853")
GSE232853_V2 <- file.path(DATA_ROOT, "GSE232853_v2")
PROGNOSIS   <- file.path(DATA_ROOT, "DLBCL_prognosis")
GDSC2       <- file.path(DATA_ROOT, "GDSC2")
SPATIAL     <- file.path(DATA_ROOT, "spatial_analysis")
CNV_DIR     <- file.path(DATA_ROOT, "CNV")

# ── 输出目录 ──────────────────────────────────────────────────────────────────
RESULTS <- Sys.getenv("DLBCL_OUT_ROOT")
if (!nzchar(RESULTS)) RESULTS <- file.path(REPO_ROOT, "results")
FIG_DIR <- file.path(RESULTS, "figures")
TAB_DIR <- file.path(RESULTS, "tables")

# ── 历史路径 → 新路径映射（与 paths.py 的 LEGACY_MAP 一致）────────────────────
LEGACY_MAP <- list(
  "d:/bulk-download/gse182434"        = "GSE182434",
  "d:/bulk-download/celloracle0331"   = "CELLORACLE",
  "d:/bulk-download/gse232853_v2"     = "GSE232853_V2",
  "d:/bulk-download/gse232853"        = "GSE232853",
  "d:/bulk-download/dlbcl_prognosis"  = "PROGNOSIS",
  # 需注册 / 手工获取的外部文件（放数据根的 external/ 下，见 data/README.md）
  "d:/badidunetdiskdownload/r/gse10846survive" = "external/GSE10846survive",
  "d:/badidunetdiskdownload/r/tcell"           = "external/Tcell",
  # 历史服务器暂存目录（原始母脚本用；见 section6/code/legacy/README.md）
  "/workspace"                                 = "external/workspace",
  "d:/bulk-download"                  = ".",
  "/mnt/results/dlbcl_prognosis_v2"   = "PROGNOSIS",
  "/mnt/results/dlbcl_prognosis"      = "PROGNOSIS",
  "/mnt/results/gse182434"            = "GSE182434",
  "/mnt/results"                      = ".",
  "/data"                             = "CELLORACLE"
)

NAME_MAP <- c(
  "GSE182434"   = GSE182434,
  "CELLORACLE"  = CELLORACLE,
  "GSE232853_V2" = GSE232853_V2,
  "GSE232853"   = GSE232853,
  "PROGNOSIS"   = PROGNOSIS,
  "CNV"         = CNV_DIR
)

#' 把历史硬编码绝对路径翻译到当前数据根下
translate_path <- function(path) {
  s <- gsub("\\\\", "/", path)
  low <- tolower(sub("/$", "", s))
  for (prefix in names(LEGACY_MAP)) {
    if (low == prefix) {
      target <- LEGACY_MAP[[prefix]]
      return(if (target %in% names(NAME_MAP)) unname(NAME_MAP[target]) else file.path(DATA_ROOT, target))
    }
    if (startsWith(low, paste0(prefix, "/"))) {
      rest <- substring(s, nchar(prefix) + 2)
      target <- LEGACY_MAP[[prefix]]
      base <- if (target %in% names(NAME_MAP)) unname(NAME_MAP[target]) else file.path(DATA_ROOT, target)
      return(file.path(base, rest))
    }
  }
  path
}

ensure_dir <- function(...) {
  for (d in list(...)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(NULL)
}
