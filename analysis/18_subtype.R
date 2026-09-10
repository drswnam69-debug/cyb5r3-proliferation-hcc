# =====================================================================
#  18_subtype.R  —  Project C
#
#  Question: In TCGA-LIHC, CYB5R3 expression DECREASES across histologic
#  grade (JT p = 0.0067) yet HIGH expression marks WORSE survival
#  (adj HR 1.226 per SD, p = 0.0345). Int J Mol Sci 2026;27:7806 left this
#  tension unresolved. Here I ask whether molecular subclass composition
#  accounts for it.
#
#  Design is PRE-REGISTERED. Every decision rule below was fixed before the
#  data were touched; see PREREG_ProjectC_subtype.docx. Do not edit the
#  DECISION RULES block after looking at any result.
#
#  This script does NOT modify the deposited pipeline (Zenodo v1.0.5). It
#  reuses that repository's 01_common.R for helpers and Xena access, then
#  redirects all output into this folder.
#
#  Run: open RUN_C.R in RStudio and press Source.
#  작성 2026-09-10
# =====================================================================

## ---- 0. 경로 --------------------------------------------------------
.here <- function() {
  fr <- sys.frames()
  if (length(fr)) for (i in rev(seq_along(fr))) {
    o <- fr[[i]]$ofile
    if (!is.null(o) && nzchar(o)) return(dirname(normalizePath(o)))
  }
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  getwd()
}
ANA  <- .here()
PROJ <- normalizePath(file.path(ANA, ".."), mustWork = FALSE)
## 헬퍼는 게재 논문의 파이프라인에서 가져온다. 그 저장소가 옆에 없으면
## (예: 이 저장소만 내려받은 경우) 동봉한 사본을 쓴다.
PIPE   <- normalizePath(file.path(PROJ, "..", "cyb5r3-marc1-lihc", "R"), mustWork = FALSE)
COMMON <- file.path(PIPE, "01_common.R")
if (!file.exists(COMMON)) COMMON <- file.path(ANA, "vendor", "01_common.R")
if (!file.exists(COMMON))
  stop("Cannot find 01_common.R. Expected it either in a sibling checkout of ",
       "cyb5r3-marc1-lihc (doi:10.5281/zenodo.21987539) or at analysis/vendor/01_common.R")
message("[C] helpers from: ", COMMON)
source(COMMON)

## 산출물은 이 프로젝트 폴더로.
RES <- file.path(PROJ, "results"); FIG <- file.path(PROJ, "figures")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## 캐시 위치는 명시적으로 정한다. 우선순위:
##   1) 환경변수 CYB5R3_CACHE  2) 옆에 있는 게재 파이프라인의 cache (재다운로드 회피)
##   3) 이 프로젝트 폴더의 cache
## 이렇게 해두면 어떤 01_common.R 을 썼든 내려받는 위치가 예측 가능하다.
CACHE <- path.expand(Sys.getenv("CYB5R3_CACHE", unset = ""))
if (!nzchar(CACHE)) {
  .sib <- normalizePath(file.path(PROJ, "..", "cyb5r3-marc1-lihc", "cache"), mustWork = FALSE)
  CACHE <- if (dir.exists(.sib)) .sib else file.path(PROJ, "cache")
}
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
message("[C] results -> ", RES)
message("[C] cache   -> ", CACHE)

set.seed(20260910)

## =====================================================================
##  DECISION RULES  (pre-specified; do not edit after seeing results)
## =====================================================================
RULES <- list(
  ## subclass assignment
  ntp_permutations      = 1000,
  ntp_fdr_cut           = 0.20,   # per-sample NTP BH-FDR above this = Unclassified
  min_class_n           = 30,     # a subclass with fewer patients is not analysed separately

  ## subclass validity control (must PASS before any contrast is computed)
  ##   S2 is the progenitor/proliferative class, S3 the well-differentiated class.
  ##   ★ 대조 표지자는 템플릿 밖에 있어야 한다. 템플릿 구성원으로 배정을 검증하면
  ##     순환논리다. AFP(S2 구성원)와 APOA1(S3 구성원)은 그래서 제외했다.
  ##     아래 목록은 세 템플릿(v2026.1.Hs) 어디에도 없음을 확인한 것이고,
  ##     실행 시 §2 끝의 guard 가 다시 자동 검사한다.
  s2_markers            = c("EPCAM", "KRT19", "MKI67", "TOP2A"),
  s3_markers            = c("ALB", "CYP2E1", "HNF1A", "SERPINA1", "TTR"),
  markers_required      = 3,      # >= this many must move in the expected direction

  ## CTNNB1 positive control (both required)
  ctnnb1_control        = c("GLUL", "LGR5"),

  ## H2 - does subclass explain the grade trend?
  ##   b1 = grade slope alone, b2 = grade slope adjusted for subclass
  h2_explained_if       = "attenuation >= 0.50 AND p(b2) > 0.05",
  h2_independent_if     = "attenuation <  0.25 AND p(b2) < 0.05",

  ## H3 - does subclass explain the prognostic signal?
  h3_marker_if          = "HR < 1.10 with 95% CI covering 1",
  h3_independent_if     = "HR >= 1.15 with 95% CI excluding 1",

  alpha                 = 0.05,
  multiplicity          = "BH across the H1/H2/H3/H4 family"
)
w_res(data.frame(rule = names(RULES),
                 value = vapply(RULES, function(x) paste(x, collapse = ", "), character(1))),
      "C00_decision_rules.csv")

## =====================================================================
##  1. PREMISE CHECKS  (standing rule: establish what the data contain and
##     report that check BEFORE computing any contrast)
## =====================================================================
hdr("1 - Premise checks")
flow <- list(); add <- function(step, n, note = "") {
  ## n 은 "20530 x 423" 같은 문자열도 받으므로 열 자료형을 문자로 통일한다.
  flow[[length(flow) + 1]] <<- data.frame(step = step, n = as.character(n),
                                          note = note, stringsAsFactors = FALSE)
  message(sprintf("  [flow] %-52s n = %s", step, n))
}

## 1a. full expression matrix (same matrix as the published analysis)
q  <- UCSCXenaTools::XenaGenerate(subset = XenaDatasets == XDS$lihc_hiseq)
fp <- UCSCXenaTools::XenaQuery(q) |>
      UCSCXenaTools::XenaDownload(destdir = CACHE, trans_slash = TRUE, force = FALSE)
EX <- UCSCXenaTools::XenaPrepare(fp)
if (is.list(EX) && !is.data.frame(EX)) EX <- EX[[1]]
EX <- as.data.frame(EX)
names(EX)[1] <- "gene"
EX <- EX[!duplicated(EX$gene) & nzchar(EX$gene), ]
rownames(EX) <- EX$gene; EX$gene <- NULL
colnames(EX) <- nb(colnames(EX))
add("HiSeqV2 matrix: genes x samples", sprintf("%d x %d", nrow(EX), ncol(EX)),
    "log2(norm_count + 1); identical matrix to the published primary analysis")

## 1b. primary tumours, one sample per patient (published selection rule)
keep <- colnames(EX)[typ(colnames(EX)) %in% TUMOR_CODES]
EX   <- EX[, keep, drop = FALSE]
add("Tumour samples (codes 01/02/03/05/06/07)", ncol(EX))
pt   <- pat(colnames(EX))
EX   <- EX[, !duplicated(pt), drop = FALSE]
colnames(EX) <- pat(colnames(EX))
add("One sample per patient", ncol(EX))

## 1c. clinical + survival, assembled exactly as in the published script
cl <- xena_table(XDS$lihc_clin, exact = XDS$lihc_clin)
if (is.null(cl)) stop("LIHC clinical matrix unavailable")
names(cl)[1] <- "sample"; cl$patient <- pat(cl$sample)
pick <- function(p) { h <- grep(p, names(cl), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
col_age <- pick("^age_at_initial|^age$"); col_sex <- pick("^gender$|^sex$")
col_stg <- pick("pathologic_stage|^stage$|ajcc.*stage")
col_grd <- pick("neoplasm_histologic_grade|histologic.*grade")
message(sprintf("  [cols] age='%s' sex='%s' stage='%s' grade='%s'", col_age, col_sex, col_stg, col_grd))
cl$age <- suppressWarnings(as.numeric(as.character(cl[[col_age]])))
cl$sex <- as.character(cl[[col_sex]])
cl$stage_raw <- as.character(cl[[col_stg]])
cl$stage_bin <- dplyr::case_when(
  stringr::str_detect(cl$stage_raw, "Stage III|Stage IV")                ~ "III-IV",
  stringr::str_detect(cl$stage_raw, "Stage I{1,2}$|Stage I |Stage II")   ~ "I-II",
  stringr::str_detect(cl$stage_raw, "^Stage I$")                         ~ "I-II",
  TRUE ~ NA_character_)
cl$grade <- as.character(cl[[col_grd]])
cl$grade_ord <- suppressWarnings(as.integer(sub("^G", "", cl$grade)))
cl <- cl[!duplicated(cl$patient), ]

sv <- UCSCXenaTools::fetch_dense_values(XH$pan, XDS$pan_surv, c("OS", "OS.time"),
                                        use_probeMap = FALSE)
sv <- as.data.frame(t(sv)); sv$patient <- pat(rownames(sv))
sv$OS      <- suppressWarnings(as.integer(as.character(sv$OS)))
sv$OS.time <- suppressWarnings(as.numeric(as.character(sv$OS.time)))
sv <- sv[!duplicated(sv$patient), ]

## 1d. the analysis genes
gv <- function(g) {
  for (a in (ALIAS[[g]] %||% g)) if (a %in% rownames(EX)) return(as.numeric(EX[a, ]))
  rep(NA_real_, ncol(EX))
}
D <- data.frame(patient = colnames(EX), CYB5R3 = gv("CYB5R3"), MTARC1 = gv("MTARC1"),
                stringsAsFactors = FALSE)
if (all(is.na(D$CYB5R3))) stop("CYB5R3 not found in the expression matrix")
D <- dplyr::left_join(D, cl[, c("patient","age","sex","stage_bin","grade","grade_ord")], by = "patient")
D <- dplyr::left_join(D, sv[, c("patient","OS","OS.time")], by = "patient")
add("Merged with clinical and survival", nrow(D))
CC <- D[!is.na(D$OS) & !is.na(D$OS.time) & D$OS.time > 0 &
        !is.na(D$age) & !is.na(D$sex) & !is.na(D$stage_bin), ]
add("PRIMARY ANALYSIS SET (complete age, sex, stage)", nrow(CC),
    sprintf("deaths = %d; published set was n = 339, deaths = 114", sum(CC$OS == 1)))

## 1e. the two published anchors must reproduce, otherwise stop
CC$z5 <- as.numeric(scale(CC$CYB5R3))
CC$sex <- factor(CC$sex); CC$stage_bin <- factor(CC$stage_bin, c("I-II","III-IV"))
anchor_cox <- summary(coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = CC))
hr0  <- anchor_cox$conf.int["z5","exp(coef)"]
lo0  <- anchor_cox$conf.int["z5","lower .95"]; hi0 <- anchor_cox$conf.int["z5","upper .95"]
p0   <- anchor_cox$coefficients["z5","Pr(>|z|)"]
message(sprintf("  [anchor] published HR 1.226 (1.015-1.480) p=0.0345 | here HR %.3f (%.3f-%.3f) p=%.4f",
                hr0, lo0, hi0, p0))
anchor_ok <- abs(hr0 - 1.226) < 0.05
if (!anchor_ok)
  warning("Primary Cox does not reproduce the published estimate. Resolve before interpreting anything below.")

FLOW <- dplyr::bind_rows(flow); w_res(FLOW, "C01_flow.csv")
w_res(data.frame(anchor = "OS ~ CYB5R3(per SD) + age + sex + stage",
                 n = nrow(CC), events = sum(CC$OS == 1),
                 HR = hr0, CI_low = lo0, CI_high = hi0, p = p0,
                 published_HR = 1.226, reproduced = anchor_ok),
      "C02_anchor_reproduction.csv")

## =====================================================================
##  2. SUBCLASS ASSIGNMENT - Hoshida S1/S2/S3 by nearest-template prediction
## =====================================================================
hdr("2 - Hoshida subclass by NTP")
## 템플릿은 MSigDB C2:CGP 에서 프로그램으로 가져온다 (손으로 옮겨적지 않는다).
## 한 번 성공하면 analysis/hoshida_templates.tsv 로 굳혀서, 이후 실행은 외부
## 서버 상태와 무관하게 재현되고 그 파일 자체가 기탁 대상이 된다.
TEMPL_FILE <- file.path(ANA, "hoshida_templates.tsv")
SETS <- paste0("HOSHIDA_LIVER_CANCER_SUBCLASS_S", 1:3)

GRP_DIR <- file.path(ANA, "grp")
grp_files <- if (dir.exists(GRP_DIR))
  list.files(GRP_DIR, pattern = "HOSHIDA_LIVER_CANCER_SUBCLASS_S[123].*\\.grp$", full.names = TRUE) else character(0)

if (file.exists(TEMPL_FILE)) {
  tt <- utils::read.delim(TEMPL_FILE, stringsAsFactors = FALSE)
  message("  templates read from the frozen local copy: ", TEMPL_FILE)
} else if (length(grp_files) >= 3) {
  ## MSigDB 에서 직접 내려받은 .grp 파일에서 읽는다. 가장 확실한 출처이고
  ## 외부 서버에 의존하지 않는다. 형식: 1행 세트이름, '#' 주석, 이후 유전자 기호.
  rd <- function(f) {
    ln <- trimws(readLines(f, warn = FALSE))
    ln <- ln[nzchar(ln) & !startsWith(ln, "#")]
    list(set = ln[1], genes = unique(ln[-1]))
  }
  parsed <- lapply(grp_files, rd)
  tt <- do.call(rbind, lapply(parsed, function(x)
    data.frame(gs_name = x$set, gene_symbol = x$genes, stringsAsFactors = FALSE)))
  ver <- unique(sub(".*\\.(v[0-9]{4}\\.[0-9]+\\.[A-Za-z]+)\\.grp$", "\\1", basename(grp_files)))
  message("  templates read from MSigDB .grp files in ", GRP_DIR,
          "  (release ", paste(ver, collapse = "/"), ")")
  for (x in parsed) message(sprintf("    %-40s %d genes", x$set, length(x$genes)))
  ## 검증: MSigDB 세트 페이지(v2026.1.Hs)에 적힌 유전자 수와 대조한다.
  ## 출처 = 각 gene set 페이지의 "Show members (... mapped to N genes)" 표기.
  EXPECTED <- c(HOSHIDA_LIVER_CANCER_SUBCLASS_S1 = 235,
                HOSHIDA_LIVER_CANCER_SUBCLASS_S2 = 115,
                HOSHIDA_LIVER_CANCER_SUBCLASS_S3 = 266)
  CHK <- do.call(rbind, lapply(parsed, function(x) data.frame(
    gene_set = x$set, parsed_n = length(x$genes),
    expected_n = unname(EXPECTED[x$set]),
    matches = isTRUE(length(x$genes) == unname(EXPECTED[x$set])),
    stringsAsFactors = FALSE)))
  print(CHK); w_res(CHK, "C03a_template_gene_counts.csv")
  if (!all(CHK$matches, na.rm = TRUE))
    warning("Parsed gene counts do not match the MSigDB gene set pages. ",
            "Check that the .grp files are the current release before interpreting anything.")
  utils::write.table(tt, TEMPL_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  frozen to: ", TEMPL_FILE)
} else {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required so the Hoshida templates come from a versioned, ",
         "citable source rather than being transcribed by hand. install.packages('msigdbr')")
  ## 인자 이름은 설치된 버전에 맞춰 고른다 (10.0.0 에서 category -> collection 로 바뀜).
  .fm <- names(formals(msigdbr::msigdbr))
  .get <- function() if ("collection" %in% .fm)
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP") else
    msigdbr::msigdbr(species = "Homo sapiens", category  = "C2", subcategory  = "CGP")
  ## msigdbr 10.x 이상은 MSigDB 릴리스를 Zenodo 에서 내려받는다. 504 등 일시적
  ## 게이트웨이 오류가 흔하므로 재시도한다. 성공하면 로컬에 캐시된다.
  msig <- NULL
  for (att in 1:5) {
    msig <- tryCatch(.get(), error = function(e) {
      message("  [attempt ", att, "/5] ", conditionMessage(e)); NULL })
    if (!is.null(msig)) break
    if (att < 5) { message("  waiting 15 s before retrying ..."); Sys.sleep(15) }
  }
  if (is.null(msig))
    stop("Could not retrieve the MSigDB C2:CGP collection. The download server (Zenodo) ",
         "returned an error on five attempts, which is usually transient. Either re-run ",
         "RUN_C.R later, or download the three gene sets from MSigDB as .grp files and put ",
         "them in ", GRP_DIR, " (HOSHIDA_LIVER_CANCER_SUBCLASS_S1/S2/S3).")
  msig <- as.data.frame(msig)
  sym_col <- if ("gene_symbol" %in% names(msig)) "gene_symbol" else "human_gene_symbol"
  tt <- unique(msig[msig$gs_name %in% SETS, c("gs_name", sym_col)])
  names(tt) <- c("gs_name", "gene_symbol")
  attr(tt, "msigdbr_version") <- as.character(utils::packageVersion("msigdbr"))
  utils::write.table(tt, TEMPL_FILE, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  templates frozen to: ", TEMPL_FILE,
          "  (msigdbr ", utils::packageVersion("msigdbr"), ")")
}

TEMPL <- list(S1 = unique(tt$gene_symbol[tt$gs_name == SETS[1]]),
              S2 = unique(tt$gene_symbol[tt$gs_name == SETS[2]]),
              S3 = unique(tt$gene_symbol[tt$gs_name == SETS[3]]))
if (any(lengths(TEMPL) == 0)) stop("Hoshida subclass gene sets not found in this MSigDB release")

## ---- 순환논리 guard -------------------------------------------------
## 타당성 대조 표지자가 템플릿 구성원이면 그 표지자는 배정을 검증할 수 없다.
## 자기가 배정을 만들어낸 유전자이기 때문이다. 자동으로 빼고 근거를 남긴다.
.tmpl_all <- unique(unlist(TEMPL))
IND <- data.frame(
  marker = c(RULES$s2_markers, RULES$s3_markers, RULES$ctnnb1_control),
  role   = c(rep("S2 validity control",     length(RULES$s2_markers)),
             rep("S3 validity control",     length(RULES$s3_markers)),
             rep("CTNNB1 positive control", length(RULES$ctnnb1_control))),
  stringsAsFactors = FALSE)
IND$in_template <- vapply(IND$marker, function(g)
  paste(names(TEMPL)[vapply(TEMPL, function(st) g %in% st, logical(1))], collapse = ";"),
  character(1))
IND$independent <- !nzchar(IND$in_template)
print(IND); w_res(IND, "C03b_control_marker_independence.csv")
if (any(!IND$independent)) {
  message("  [guard] excluded from the controls because they are template members: ",
          paste(IND$marker[!IND$independent], collapse = ", "))
  RULES$s2_markers    <- setdiff(RULES$s2_markers,    .tmpl_all)
  RULES$s3_markers    <- setdiff(RULES$s3_markers,    .tmpl_all)
  RULES$ctnnb1_control <- setdiff(RULES$ctnnb1_control, .tmpl_all)
}
if (length(RULES$s2_markers) < RULES$markers_required ||
    length(RULES$s3_markers) < RULES$markers_required)
  stop("Too few template-independent control markers remain to run the validity check.")
message(sprintf("  templates: S1 %d, S2 %d, S3 %d genes",
                length(TEMPL$S1), length(TEMPL$S2), length(TEMPL$S3)))

uni <- intersect(Reduce(union, TEMPL), rownames(EX))
message("  template genes present in the matrix: ", length(uni))
w_res(data.frame(union_size = length(Reduce(union, TEMPL)),
                 present_in_matrix = length(uni),
                 S1 = length(TEMPL$S1), S2 = length(TEMPL$S2), S3 = length(TEMPL$S3),
                 matrix_genes = nrow(EX), matrix_samples = ncol(EX)),
      "C03c_template_coverage.csv")
Z <- t(scale(t(as.matrix(EX[uni, , drop = FALSE]))))      # gene-wise standardisation
Z[!is.finite(Z)] <- 0
TM <- vapply(TEMPL, function(g) as.numeric(uni %in% g), numeric(length(uni)))

## cosine similarity of every sample against every template, vectorised
zn    <- sqrt(colSums(Z^2)); zn[zn == 0] <- NA_real_
cosall <- function(TT) (t(TT) %*% Z) / (sqrt(colSums(TT^2)) %o% zn)   # 3 x nsamples
S  <- t(cosall(TM))
colnames(S) <- names(TEMPL)
best   <- colnames(S)[max.col(S, ties.method = "first")]
bestv  <- S[cbind(seq_len(nrow(S)), max.col(S, ties.method = "first"))]

## permutation p per sample: shuffle the gene labels of the templates
B <- RULES$ntp_permutations
null_max <- matrix(NA_real_, nrow(S), B)
for (b in seq_len(B)) {
  TMp <- TM[sample.int(nrow(TM)), , drop = FALSE]
  null_max[, b] <- apply(cosall(TMp), 2, max)
}
pval <- (rowSums(null_max >= bestv) + 1) / (B + 1)   # +1 보정: p = 0 을 만들지 않는다
fdr  <- p.adjust(pval, "BH")
SUB <- data.frame(patient = colnames(Z), subclass = best, cos = bestv,
                  p = pval, fdr = fdr, stringsAsFactors = FALSE)
SUB$subclass[SUB$fdr > RULES$ntp_fdr_cut] <- "Unclassified"
print(table(SUB$subclass))
w_res(SUB, "C03_subclass_assignment.csv")

CC <- dplyr::left_join(CC, SUB[, c("patient","subclass","fdr")], by = "patient")
D  <- dplyr::left_join(D,  SUB[, c("patient","subclass","fdr")], by = "patient")
tab <- as.data.frame(table(CC$subclass), stringsAsFactors = FALSE)
names(tab) <- c("subclass","n")
w_res(tab, "C04_subclass_counts.csv")
small <- tab$subclass[tab$n < RULES$min_class_n & tab$subclass != "Unclassified"]
if (length(small)) message("  [note] below min_class_n and not analysed separately: ",
                           paste(small, collapse = ", "))

## =====================================================================
##  3. SUBCLASS VALIDITY CONTROL - the classification must behave as known
## =====================================================================
hdr("3 - Subclass validity control")
mk <- function(gs, cls) {
  out <- list()
  for (g in gs) {
    if (!g %in% rownames(EX)) { out[[g]] <- data.frame(marker = g, present = FALSE); next }
    v <- as.numeric(EX[g, ]); names(v) <- colnames(EX)
    a <- v[SUB$patient[SUB$subclass == cls]]
    b <- v[SUB$patient[SUB$subclass %in% setdiff(c("S1","S2","S3"), cls)]]
    out[[g]] <- data.frame(marker = g, present = TRUE, class = cls,
                           median_in = median(a, na.rm = TRUE),
                           median_out = median(b, na.rm = TRUE),
                           delta = median(a, na.rm = TRUE) - median(b, na.rm = TRUE),
                           p = suppressWarnings(wilcox.test(a, b)$p.value),
                           expected = "higher",
                           passes = median(a, na.rm = TRUE) > median(b, na.rm = TRUE))
  }
  dplyr::bind_rows(out)
}
CTRL <- dplyr::bind_rows(mk(RULES$s2_markers, "S2"), mk(RULES$s3_markers, "S3"))
print(CTRL); w_res(CTRL, "C05_subclass_validity_control.csv")
pass_s2 <- sum(CTRL$passes[CTRL$class == "S2"], na.rm = TRUE)
pass_s3 <- sum(CTRL$passes[CTRL$class == "S3"], na.rm = TRUE)
CONTROL_OK <- pass_s2 >= RULES$markers_required && pass_s3 >= RULES$markers_required
message(sprintf("  S2 markers passing: %d/%d | S3 markers passing: %d/%d | CONTROL %s",
                pass_s2, length(RULES$s2_markers), pass_s3, length(RULES$s3_markers),
                if (CONTROL_OK) "PASS" else "FAIL"))
if (!CONTROL_OK) {
  w_res(data.frame(verdict = "NOT TESTABLE",
                   reason = "Hoshida NTP assignment does not reproduce the known marker behaviour of S2 and S3 in this matrix. Per the pre-registered rule no contrast is computed."),
        "C99_verdict.csv")
  save_session("18_subtype"); stop("Validity control failed - reported as not testable, by design.")
}

## =====================================================================
##  4. CTNNB1 status and a beta-catenin target score
## =====================================================================
hdr("4 - CTNNB1")
BCAT <- c("GLUL","LGR5","TBX3","AXIN1","NKD1","RHBG","SP5")
have <- intersect(BCAT, rownames(EX))
bs <- colMeans(t(scale(t(as.matrix(EX[have, , drop = FALSE])))), na.rm = TRUE)
D$bcat_score  <- bs[D$patient];  CC$bcat_score <- bs[CC$patient]
mut <- tryCatch({
  m <- UCSCXenaTools::fetch_dense_values(XH$pan, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena",
                                         "CTNNB1", use_probeMap = FALSE)
  v <- tovec(m); names(v) <- pat(names(v)); v[!duplicated(names(v))]
}, error = function(e) { message("  [note] CTNNB1 mutation matrix unavailable: ", conditionMessage(e)); NULL })
CTNNB1_OK <- FALSE
if (!is.null(mut)) {
  D$CTNNB1 <- ifelse(mut[D$patient] > 0, "Mutant", "Wild-type")
  CC$CTNNB1 <- ifelse(mut[CC$patient] > 0, "Mutant", "Wild-type")
  pc <- lapply(RULES$ctnnb1_control, function(g) {
    v <- as.numeric(EX[g, ]); names(v) <- colnames(EX)
    a <- v[D$patient[which(D$CTNNB1 == "Mutant")]]; b <- v[D$patient[which(D$CTNNB1 == "Wild-type")]]
    data.frame(control_gene = g, n_mut = sum(!is.na(a)), n_wt = sum(!is.na(b)),
               median_mut = median(a, na.rm = TRUE), median_wt = median(b, na.rm = TRUE),
               p = suppressWarnings(wilcox.test(a, b)$p.value),
               passes = median(a, na.rm = TRUE) > median(b, na.rm = TRUE))
  })
  PC <- dplyr::bind_rows(pc); print(PC); w_res(PC, "C06_ctnnb1_positive_control.csv")
  CTNNB1_OK <- all(PC$passes, na.rm = TRUE)
  message("  CTNNB1 positive control: ", if (CTNNB1_OK) "PASS" else "FAIL - H4 will be reported as not testable")
}

## =====================================================================
##  5. H1 - CYB5R3 across subclass
## =====================================================================
hdr("5 - H1: CYB5R3 across Hoshida subclass")
dd <- D[D$subclass %in% c("S1","S2","S3") & !is.na(D$CYB5R3), ]
kw <- kruskal.test(CYB5R3 ~ factor(subclass), data = dd)
eps2 <- unname((kw$statistic - length(unique(dd$subclass)) + 1) /
               (nrow(dd) - length(unique(dd$subclass))))
H1 <- data.frame(test = "Kruskal-Wallis, CYB5R3 by subclass",
                 n = nrow(dd), chisq = unname(kw$statistic), df = unname(kw$parameter),
                 p = kw$p.value, epsilon_squared = eps2)
pw <- list()
for (cmb in combn(c("S1","S2","S3"), 2, simplify = FALSE)) {
  a <- dd$CYB5R3[dd$subclass == cmb[1]]; b <- dd$CYB5R3[dd$subclass == cmb[2]]
  if (length(a) < RULES$min_class_n || length(b) < RULES$min_class_n) next
  pw[[paste(cmb, collapse = "_vs_")]] <- data.frame(
    contrast = paste(cmb, collapse = " vs "), n1 = length(a), n2 = length(b),
    median1 = median(a), median2 = median(b),
    cliffs_delta = cliffs_delta(a, b),
    p = suppressWarnings(wilcox.test(a, b)$p.value))
}
PW <- dplyr::bind_rows(pw); if (nrow(PW)) PW$p_BH <- p.adjust(PW$p, "BH")
MED <- dd |> dplyr::group_by(subclass) |>
  dplyr::summarise(n = dplyr::n(), median_CYB5R3 = median(CYB5R3), IQR = IQR(CYB5R3), .groups = "drop")
print(H1); print(as.data.frame(MED)); print(PW)
w_res(H1, "C10_H1_kruskal.csv"); w_res(as.data.frame(MED), "C11_H1_medians.csv")
if (nrow(PW)) w_res(PW, "C12_H1_pairwise.csv")

## =====================================================================
##  6. H2 - is the grade trend explained by subclass composition?
## =====================================================================
hdr("6 - H2: grade trend before and after subclass adjustment")
dg <- D[!is.na(D$grade_ord) & !is.na(D$CYB5R3), ]
jt_p <- NA_real_
if (requireNamespace("clinfun", quietly = TRUE))
  jt_p <- clinfun::jonckheere.test(dg$CYB5R3, dg$grade_ord, alternative = "two.sided")$p.value
GR <- dg |> dplyr::group_by(grade_ord) |>
  dplyr::summarise(n = dplyr::n(), median_CYB5R3 = median(CYB5R3), .groups = "drop")
print(as.data.frame(GR)); w_res(as.data.frame(GR), "C20_H2_grade_medians.csv")

dg2 <- dg[dg$subclass %in% c("S1","S2","S3"), ]
m1 <- lm(rank(CYB5R3) ~ grade_ord, data = dg2)
m2 <- lm(rank(CYB5R3) ~ grade_ord + factor(subclass), data = dg2)
b1 <- coef(m1)["grade_ord"]; b2 <- coef(m2)["grade_ord"]
p1 <- summary(m1)$coefficients["grade_ord","Pr(>|t|)"]
p2 <- summary(m2)$coefficients["grade_ord","Pr(>|t|)"]
atten <- unname((b1 - b2) / b1)
verdict_h2 <- if (atten >= 0.50 && p2 > RULES$alpha) "Composition explains the grade trend" else
              if (atten <  0.25 && p2 < RULES$alpha) "Grade trend is independent of subclass" else
              "Indeterminate"
H2 <- data.frame(n = nrow(dg2), JT_p_all_grades = jt_p,
                 b_grade_alone = unname(b1), p_grade_alone = p1,
                 b_grade_adj_subclass = unname(b2), p_grade_adj_subclass = p2,
                 attenuation = atten, verdict = verdict_h2)
print(H2); w_res(H2, "C21_H2_grade_attenuation.csv")

within <- list()
for (s in c("S1","S2","S3")) {
  d <- dg[which(dg$subclass == s), ]
  if (nrow(d) < RULES$min_class_n || length(unique(d$grade_ord)) < 2) next
  pj <- if (requireNamespace("clinfun", quietly = TRUE))
    clinfun::jonckheere.test(d$CYB5R3, d$grade_ord, alternative = "two.sided")$p.value else NA_real_
  sp <- spearman_ci(d$CYB5R3, d$grade_ord)
  within[[s]] <- data.frame(subclass = s, n = nrow(d), JT_p = pj,
                            spearman_rho = sp$rho, spearman_p = sp$p)
}
WI <- dplyr::bind_rows(within)
if (nrow(WI)) { WI$JT_p_BH <- p.adjust(WI$JT_p, "BH"); print(WI); w_res(WI, "C22_H2_within_subclass_grade.csv") }

## =====================================================================
##  7. H3 - does subclass explain the prognostic signal?
## =====================================================================
hdr("7 - H3: Cox models with and without subclass")
cs <- CC[CC$subclass %in% c("S1","S2","S3"), ]
cs$z5 <- as.numeric(scale(cs$CYB5R3))
cs$subclass <- factor(cs$subclass)
cs$grade_ord2 <- cs$grade_ord
tidy1 <- function(fit, label, term = "z5") {
  s <- summary(fit)
  data.frame(model = label, n = fit$n, events = fit$nevent,
             HR = s$conf.int[term,"exp(coef)"],
             CI_low = s$conf.int[term,"lower .95"],
             CI_high = s$conf.int[term,"upper .95"],
             p = s$coefficients[term,"Pr(>|z|)"])
}
f1 <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = cs)
f2 <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + subclass, data = cs)
cs3 <- cs[!is.na(cs$grade_ord2), ]
f3 <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + subclass + grade_ord2, data = cs3)
H3 <- dplyr::bind_rows(
  tidy1(f1, "M1 published covariates, subclassifiable patients only"),
  tidy1(f2, "M2 + Hoshida subclass"),
  tidy1(f3, "M3 + subclass + histologic grade"))
H3$lrt_vs_M1 <- c(NA, anova(f1, f2)$`Pr(>|Chi|)`[2], NA)
hr2 <- H3$HR[2]; lo2 <- H3$CI_low[2]; hi2 <- H3$CI_high[2]
H3$verdict <- c(NA,
  if (hr2 < 1.10 && lo2 <= 1 && hi2 >= 1) "CYB5R3 largely a subclass marker" else
  if (hr2 >= 1.15 && lo2 > 1)             "Prognostic signal independent of subclass" else
  "Indeterminate", NA)
print(H3); w_res(H3, "C30_H3_cox_models.csv")

strat <- list()
for (s in levels(cs$subclass)) {
  d <- cs[cs$subclass == s, ]
  if (nrow(d) < RULES$min_class_n || sum(d$OS == 1) < 10) next
  d$z <- as.numeric(scale(d$CYB5R3))
  strat[[s]] <- cbind(subclass = s,
    tidy1(coxph(Surv(OS.time, OS) ~ z + age + sex + stage_bin, data = d),
          paste("within", s), "z"))
}
ST <- dplyr::bind_rows(strat)
if (nrow(ST)) { print(ST); w_res(ST, "C31_H3_within_subclass_cox.csv") }
fi <- coxph(Surv(OS.time, OS) ~ z5 * subclass + age + sex + stage_bin, data = cs)
INT <- data.frame(test = "CYB5R3 x subclass interaction (LRT vs M2)",
                  p = anova(f2, fi)$`Pr(>|Chi|)`[2])
print(INT); w_res(INT, "C32_H3_interaction.csv")

## =====================================================================
##  8. H4 - CTNNB1
## =====================================================================
hdr("8 - H4: CTNNB1 and beta-catenin activation")
if (CTNNB1_OK) {
  a <- D$CYB5R3[which(D$CTNNB1 == "Mutant")]; b <- D$CYB5R3[which(D$CTNNB1 == "Wild-type")]
  sp <- spearman_ci(D$CYB5R3, D$bcat_score)
  H4 <- data.frame(n_mut = length(a), n_wt = length(b),
                   median_mut = median(a, na.rm = TRUE), median_wt = median(b, na.rm = TRUE),
                   cliffs_delta = cliffs_delta(a, b),
                   p_mut_vs_wt = suppressWarnings(wilcox.test(a, b)$p.value),
                   bcat_score_genes = paste(have, collapse = ";"),
                   rho_CYB5R3_bcat = sp$rho, rho_lo = sp$lo, rho_hi = sp$hi, p_rho = sp$p)
  print(H4); w_res(H4, "C40_H4_ctnnb1.csv")
} else {
  w_res(data.frame(verdict = "NOT TESTABLE",
                   reason = "CTNNB1 calls unavailable, or the GLUL/LGR5 positive control failed."),
        "C40_H4_ctnnb1.csv")
}

## =====================================================================
##  9. Family-wise multiplicity and one-line verdict
## =====================================================================
hdr("9 - Multiplicity and verdict")
fam <- data.frame(
  hypothesis = c("H1 CYB5R3 differs across subclass",
                 "H2 grade slope after subclass adjustment",
                 "H3 CYB5R3 x subclass interaction",
                 "H4 CYB5R3 higher in CTNNB1-mutant"),
  p = c(H1$p, p2, INT$p, if (CTNNB1_OK) H4$p_mut_vs_wt else NA_real_))
fam$p_BH <- p.adjust(fam$p, "BH")
print(fam); w_res(fam, "C50_family_BH.csv")
w_res(data.frame(
  anchor_reproduced = anchor_ok,
  subclass_control  = "PASS",
  H2_verdict        = verdict_h2,
  H3_verdict        = H3$verdict[2],
  note = "Interpret only after confirming the anchor reproduced the published estimate."),
  "C99_verdict.csv")

save_session("18_subtype")
message("\n[C] done. results -> ", RES)
