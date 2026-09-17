# 전체 파이프라인을 한 번에 실행합니다 (Stage 1-3, PH, 그림, 임상 민감도, 참고인 분석).
# 순서는 기존 러너들과 같습니다: RUN_C_FIGURES.R -> 24_figure5.R -> RUN_C_CLINICAL.R -> 25_referee.R.
# 각 스크립트는 앞 단계가 메모리에 만들어 둔 객체를 그대로 씁니다.
# -----------------------------------------------------------------------------
# 이 파일이 있는 폴더를 작업 디렉터리로 삼는다.
#  · source("경로/RUN_*.R") 로 실행해도, RStudio 의 Source 버튼으로 실행해도 같다.
#  · 예전 판은 rstudioapi 로 '편집기에 열려 있는 파일'을 물어봤는데, 콘솔에서
#    경로로 source() 하면 엉뚱한 파일(그때 열려 있던 탭)의 폴더로 이동해
#    다른 사본이 조용히 실행됐다. ofile 을 먼저 보도록 고쳤다.
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

## 24_figure5.R must run AFTER 23_clinical.R: block D5 of 23 draws Figure 5
## without a numbers-at-risk table, and 24 redraws it with one.
for (.f in c("18_subtype.R", "19_addendum.R", "20_validation.R", "21_ph.R",
             "22_figures.R", "23_clinical.R", "25_referee.R", "24_figure5.R"))
  source(file.path("analysis", .f))

message("\n[RUN] all stages complete.")
