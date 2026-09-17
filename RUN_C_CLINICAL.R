# Stage D (POST HOC clinical robustness) 를 실행합니다.
# 18 -> 19 -> 20 을 먼저 돌린 뒤 23_clinical.R 로 이어집니다.
# -----------------------------------------------------------------------------
.run_dir <- function() {
  fr <- sys.frames()
  if (length(fr)) for (i in rev(seq_along(fr))) {
    o <- fr[[i]]$ofile
    if (!is.null(o) && is.character(o) && length(o) == 1L && nzchar(o))
      return(dirname(normalizePath(o)))
  }
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)
    if (!inherits(p, "try-error") && length(p) == 1L && nzchar(p))
      return(dirname(normalizePath(p)))
  }
  getwd()
}
setwd(.run_dir())
if (!file.exists(file.path("analysis", "18_subtype.R")))
  stop("analysis/18_subtype.R not found. 현재 작업 디렉터리: ", getwd())
message("[RUN] working directory: ", getwd())

for (.f in c("18_subtype.R", "19_addendum.R", "20_validation.R", "23_clinical.R"))
  source(file.path("analysis", .f))
