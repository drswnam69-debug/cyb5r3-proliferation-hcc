# =====================================================================
#  23_clinical.R  -  Project C, Stage D: clinical robustness
#
#  ★★ POST HOC. Designed on 2026-09-11, after the manuscript had been
#     drafted and audited from a reviewer's perspective. Nothing here is
#     pre-registered. It is reported as post hoc sensitivity analysis and
#     never alongside H1-H4 or the pre-registered replication.
#
#  Why it exists. Four gaps were identified in review:
#    D1  Overall survival mixes tumor death with liver-related death in a
#        resected cirrhotic population. The claim of a proliferation
#        independent association needs disease-specific and progression
#        endpoints, which the survival table already loaded contains.
#    D2  Neither model carries liver function or etiology, although both
#        cohorts record them.
#    D3  The change in the hazard ratio, which is the central quantity of
#        the article, has no confidence interval.
#    D4  There is no Kaplan-Meier curve and no absolute survival estimate.
#
#  Standing rule: establish what the data contain and report that check
#  BEFORE computing any contrast. Section D0 does that and sets the
#  coverage threshold below which a variable is described but not modeled.
#
#  Run: open RUN_C_CLINICAL.R and Source. It runs 18, 19, 20 first.
# =====================================================================
hdr("D - POST HOC clinical robustness")
message("  script version: 2026-09-17g")
message("  NOTE: every result below is post hoc and must be labeled as such.")

## Joining a column that is already present makes dplyr rename both copies to
## .x and .y, so a second run of this script in the same session finds neither
## and fails. Drop the incoming columns from the target first, which makes each
## join idempotent.
join_fresh <- function(target, incoming, by) {
  drop <- setdiff(names(incoming), by)
  target <- target[, setdiff(names(target), drop), drop = FALSE]
  dplyr::left_join(target, incoming, by = by)
}


RULES_D <- list(
  designed_on            = "2026-09-11, after the manuscript was drafted",
  status                 = "POST HOC, sensitivity only; no confirmatory claim rests on it",
  min_coverage_to_model  = 0.60,
  coverage_rule          = "a variable available in < 60% of a model set is described, not modeled; meeting the threshold makes a variable eligible for modeling and does not mean it entered any model reported here",
  endpoints_tcga         = "OS (primary, as published), DSS, PFI",
  endpoint_gse           = "OS (primary), recurrence-free survival",
  delta_bootstrap_B      = 2000,
  delta_rule             = "paired resamples; the same bootstrap sample fits both models",
  km_groups              = "tertiles of CYB5R3 within each cohort",
  five_year_mark         = "1825 days in TCGA-LIHC, 60 months in GSE14520"
)
w_res(data.frame(rule = names(RULES_D),
                 value = vapply(RULES_D, function(x) paste(x, collapse = ", "), character(1))),
      "CD00_clinical_rules.csv")

zsafe <- function(x) as.numeric(scale(x))
cov_of <- function(x) mean(!is.na(x))

## =====================================================================
##  D0. PREMISE CHECK - what do these cohorts actually contain?
## =====================================================================
hdr("D0 - Premise check: coverage of every clinical variable")

## ---- D0a. TCGA alternative endpoints ---------------------------------
sv2 <- tryCatch(
  UCSCXenaTools::fetch_dense_values(XH$pan, XDS$pan_surv,
                                    c("DSS", "DSS.time", "PFI", "PFI.time"),
                                    use_probeMap = FALSE),
  error = function(e) { message("  [!] could not fetch DSS/PFI: ", conditionMessage(e)); NULL })

have_alt <- !is.null(sv2)
if (have_alt) {
  sv2 <- as.data.frame(t(sv2)); sv2$patient <- pat(rownames(sv2))
  for (v in c("DSS", "PFI"))
    sv2[[v]] <- suppressWarnings(as.integer(as.character(sv2[[v]])))
  for (v in c("DSS.time", "PFI.time"))
    sv2[[v]] <- suppressWarnings(as.numeric(as.character(sv2[[v]])))
  sv2 <- sv2[!duplicated(sv2$patient), ]
  CC <- join_fresh(CC, sv2[, c("patient","DSS","DSS.time","PFI","PFI.time")], "patient")
}

## ---- D0b. TCGA liver function and etiology ---------------------------
## The clinical matrix column names are not stable across releases, so the
## columns are located by pattern and the match is printed for the record.
find_col <- function(p) { h <- grep(p, names(cl), ignore.case = TRUE, value = TRUE)
                          if (length(h)) h[1] else NA_character_ }
## NOTE (2026-09-11): the first version of this script matched "albumin" and
## "bilirubin" to the first column carrying those words, which are the assay
## REFERENCE LIMITS, not the patient's measurement. Adjusting for a laboratory
## reference bound is meaningless. The patterns below take the specified value
## where one exists. For bilirubin no such column exists in this matrix: only
## lower/upper limit fields, whose contents cannot be shown to be patient
## values, so bilirubin is declared not testable and is not modeled.
LIVER <- list(
  child_pugh = "^child_pugh_classification_grade$",
  ishak      = "^fibrosis_ishak_score$",
  albumin    = "^albumin_result_specified_value$",
  platelet   = "^platelet_result_count$",
  ptime      = "^prothrombin_time_result_value$",
  serology   = "^viral_hepatitis_serology$",
  residual   = "^residual_tumor$"
)
lv_cols <- vapply(LIVER, find_col, character(1))
for (nm in names(lv_cols))
  message(sprintf("  [col] %-11s -> %s", nm, ifelse(is.na(lv_cols[[nm]]), "NOT FOUND", lv_cols[[nm]])))

for (nm in names(lv_cols)) {
  v <- lv_cols[[nm]]
  cl[[paste0("lv_", nm)]] <- if (is.na(v)) NA else as.character(cl[[v]])
}
## Ishak is recorded as text such as "3,4 - Fibrous Septa"; take the first number.
cl$lv_ishak_n <- suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", cl$lv_ishak)))

## Units are not consistent across submitting centres, so each numeric field is
## brought to one scale by an explicit, reported rule rather than used raw.
unit_fix <- function(x, lo, hi, factor, label) {
  x <- suppressWarnings(as.numeric(x))
  n_conv <- sum(!is.na(x) & x > hi)
  x[!is.na(x) & x > hi] <- x[!is.na(x) & x > hi] / factor
  n_drop <- sum(!is.na(x) & (x < lo | x > hi))
  x[!is.na(x) & (x < lo | x > hi)] <- NA
  message(sprintf("  [units] %-10s converted %d value(s) by /%g, set %d implausible value(s) to missing",
                  label, n_conv, factor, n_drop))
  x
}
## albumin: plausible 1 to 6 g/dL; values recorded in g/L are divided by 10
cl$lv_albumin_n  <- unit_fix(cl$lv_albumin, 1, 6, 10, "albumin")
## platelet: plausible 10 to 1000 x10^9/L; values recorded per microlitre are divided by 1000
cl$lv_platelet_n <- unit_fix(cl$lv_platelet, 10, 1000, 1000, "platelet")
## prothrombin time reported as INR in most records; plausible 0.7 to 4
cl$lv_ptime_n    <- unit_fix(cl$lv_ptime, 0.7, 4, 1, "INR")
cl$lv_cp_class <- ifelse(is.na(cl$lv_child_pugh) | !nzchar(cl$lv_child_pugh), NA_character_,
                         toupper(substr(trimws(cl$lv_child_pugh), 1, 1)))
cl$lv_cp_class[!cl$lv_cp_class %in% c("A","B","C")] <- NA_character_

## Etiology. The free-text risk factor field is populated for very few
## patients, so serology is used instead. A serology record is treated as
## informative only where the field is non-empty.
se <- toupper(ifelse(is.na(cl$lv_serology), "", cl$lv_serology))
cl$et_any <- ifelse(nzchar(se), TRUE, NA)
cl$et_hbv <- ifelse(nzchar(se), grepl("HEPATITIS B SURFACE ANTIGEN", se), NA)
cl$et_hcv <- ifelse(nzchar(se), grepl("HEPATITIS +C ANTIBODY", se), NA)

lv_keep <- c("lv_albumin_n", "lv_platelet_n", "lv_ptime_n", "lv_ishak_n",
             "lv_cp_class", "lv_residual", "et_any", "et_hbv", "et_hcv")
CC <- join_fresh(CC, cl[, c("patient", lv_keep)], "patient")

## ---- D0c. GSE14520 liver variables and recurrence --------------------
sup2 <- file.path(CACHE, "GSE14520_Extra_Supplement.txt.gz")
if (!file.exists(sup2)) sup2 <- file.path(CACHE, "GSE14520_Extra_Supplement.txt")
if (!file.exists(sup2) && exists("EXTERNAL_GEO_CACHE"))
  sup2 <- list.files(dirname(EXTERNAL_GEO_CACHE), pattern = "GSE14520_Extra_Supplement",
                     recursive = TRUE, full.names = TRUE)[1]
have_gse_extra <- !is.na(sup2) && file.exists(sup2)
if (have_gse_extra) {
  s14 <- data.table::fread(sup2)
  g2 <- function(p) grep(p, names(s14), ignore.case = TRUE, value = TRUE)[1]
  ex14 <- data.frame(
    geo       = as.character(s14[[g2("Affy.?GSM")]]),
    cirrhosis = as.character(s14[[g2("^Cirrhosis")]]),
    bclc      = as.character(s14[[g2("BCLC")]]),
    clip      = as.character(s14[[g2("CLIP")]]),
    hbv       = as.character(s14[[g2("HBV")]]),
    afp       = as.character(s14[[g2("^AFP")]]),
    size      = as.character(s14[[g2("Main Tumor Size")]]),
    multinod  = as.character(s14[[g2("Multinodular")]]),
    rec_stat  = suppressWarnings(as.integer(s14[[g2("Recurr.*status")]])),
    rec_time  = suppressWarnings(as.numeric(s14[[g2("Recurr.*month")]])),
    stringsAsFactors = FALSE)
  ex14 <- ex14[!duplicated(ex14$geo), ]
  d14 <- join_fresh(d14, ex14, "geo")
  d14$cirr_bin <- ifelse(grepl("^Y|^1", toupper(trimws(d14$cirrhosis))), "yes",
                  ifelse(grepl("^N|^0", toupper(trimws(d14$cirrhosis))), "no", NA))
  d14$bclc_bin <- ifelse(grepl("^0|^A", toupper(trimws(d14$bclc))), "0-A",
                  ifelse(grepl("^B|^C", toupper(trimws(d14$bclc))), "B-C", NA))
}

## ---- D0d. report coverage BEFORE any contrast ------------------------
cov_rows <- list()
addcov <- function(cohort, variable, x, note = "") {
  cov_rows[[length(cov_rows) + 1]] <<- data.frame(
    cohort = cohort, variable = variable, n_available = sum(!is.na(x)),
    n_set = length(x), coverage = round(mean(!is.na(x)), 3),
    meets_coverage_threshold = mean(!is.na(x)) >= RULES_D$min_coverage_to_model,
    note = note, stringsAsFactors = FALSE)
}
if (have_alt) {
  addcov("TCGA-LIHC", "DSS",  CC$DSS,  "disease-specific survival, Liu et al CDR")
  addcov("TCGA-LIHC", "PFI",  CC$PFI,  "progression-free interval, Liu et al CDR")
}
addcov("TCGA-LIHC", "Child-Pugh class", CC$lv_cp_class)
addcov("TCGA-LIHC", "Ishak fibrosis",   CC$lv_ishak_n)
addcov("TCGA-LIHC", "albumin",          CC$lv_albumin_n, "specified value, unit corrected")
addcov("TCGA-LIHC", "platelet",         CC$lv_platelet_n, "unit corrected")
addcov("TCGA-LIHC", "INR",              CC$lv_ptime_n)
addcov("TCGA-LIHC", "serology recorded", CC$et_any, "etiology; the free-text risk field is near empty")
cov_rows[[length(cov_rows) + 1]] <- data.frame(
  cohort = "TCGA-LIHC", variable = "bilirubin", n_available = NA_integer_,
  n_set = nrow(CC), coverage = NA_real_, meets_coverage_threshold = FALSE,
  note = "NOT TESTABLE: the matrix carries assay limit fields only, with no column that can be shown to hold the patient value",
  stringsAsFactors = FALSE)
if (have_gse_extra) {
  addcov("GSE14520", "cirrhosis",  d14$cirr_bin)
  addcov("GSE14520", "BCLC stage", d14$bclc_bin)
  addcov("GSE14520", "recurrence", d14$rec_stat)
  addcov("GSE14520", "HBV status", d14$hbv)
}
COV <- dplyr::bind_rows(cov_rows)
print(COV); w_res(COV, "CD01_variable_coverage.csv")
message("  Variables below ", 100 * RULES_D$min_coverage_to_model,
        "% coverage are described only, by the rule fixed above.")

## ---- D0e. cohort description (what a clinical reader expects) --------
desc <- function(cohort, variable, level, n, N)
  data.frame(cohort = cohort, variable = variable, level = level,
             n = n, N = N, percent = round(100 * n / N, 1), stringsAsFactors = FALSE)
dl <- list()
N1 <- nrow(CC)
Nse <- sum(CC$et_any %in% TRUE)
dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "serology recorded", "any", Nse, N1)
if (Nse > 0) for (e in c("et_hbv", "et_hcv")) {
  lab <- c(et_hbv = "hepatitis B surface antigen", et_hcv = "hepatitis C antibody")[[e]]
  dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "serology, among those with a record", lab,
                               sum(CC[[e]] %in% TRUE), Nse)
}
if (any(!is.na(CC$lv_cp_class)))
  for (k in c("A","B","C"))
    dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "Child-Pugh class", k,
                                 sum(CC$lv_cp_class %in% k), N1)
if (any(!is.na(CC$lv_ishak_n))) {
  Nis <- sum(!is.na(CC$lv_ishak_n))
  dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "Ishak recorded", "any", Nis, N1)
  dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "Ishak 5-6, among those with a score", "5-6",
                               sum(CC$lv_ishak_n %in% c(5, 6)), Nis)
  dl[[length(dl) + 1]] <- desc("TCGA-LIHC", "Ishak 0, among those with a score", "0",
                               sum(CC$lv_ishak_n %in% 0), Nis)
}
if (have_gse_extra) {
  N2 <- nrow(d14)
  dl[[length(dl) + 1]] <- desc("GSE14520", "cirrhosis", "yes", sum(d14$cirr_bin %in% "yes"), N2)
  dl[[length(dl) + 1]] <- desc("GSE14520", "BCLC stage", "0-A", sum(d14$bclc_bin %in% "0-A"), N2)
  dl[[length(dl) + 1]] <- desc("GSE14520", "BCLC stage", "B-C", sum(d14$bclc_bin %in% "B-C"), N2)
  dl[[length(dl) + 1]] <- desc("GSE14520", "HBV positive", "yes",
                               sum(grepl("^AVR-CC|^CC|^Y|^1|POSITIVE", toupper(trimws(d14$hbv)))), N2)
}
DESC <- dplyr::bind_rows(dl)
print(DESC); w_res(DESC, "CD02_cohort_description.csv")

## =====================================================================
##  D1. Alternative endpoints in TCGA-LIHC
## =====================================================================
hdr("D1 - Does the proliferation effect appear on DSS and PFI as well?")
tidD <- function(fit, label, endpoint, term = "z5") {
  s <- summary(fit)
  data.frame(endpoint = endpoint, model = label, term = term,
             n = fit$n, events = fit$nevent,
             HR = s$conf.int[term,"exp(coef)"], CI_low = s$conf.int[term,"lower .95"],
             CI_high = s$conf.int[term,"upper .95"], p = s$coefficients[term,"Pr(>|z|)"],
             stringsAsFactors = FALSE)
}
ep_rows <- list()
ep_pairs <- list(c("OS.time","OS"))
if (have_alt) ep_pairs <- c(ep_pairs, list(c("DSS.time","DSS"), c("PFI.time","PFI")))
for (pr in ep_pairs) {
  tv <- pr[1]; ev <- pr[2]
  dd <- CC[complete.cases(CC[, c("CYB5R3","prolif","grade_ord","age","sex","stage_bin")]), ]
  dd <- dd[!is.na(dd[[tv]]) & !is.na(dd[[ev]]) & dd[[tv]] > 0, ]
  if (nrow(dd) < 50 || sum(dd[[ev]] == 1) < 20) {
    message("  [skip] ", ev, ": too few events"); next
  }
  dd$z5 <- zsafe(dd$CYB5R3); dd$zp <- zsafe(dd$prolif)
  f0 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ z5 + age + sex + stage_bin + grade_ord", tv, ev)), data = dd)
  f1 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ z5 + age + sex + stage_bin + grade_ord + zp", tv, ev)), data = dd)
  ep_rows[[length(ep_rows) + 1]] <- tidD(f0, "without proliferation", ev)
  ep_rows[[length(ep_rows) + 1]] <- tidD(f1, "with proliferation", ev)
  ep_rows[[length(ep_rows) + 1]] <- tidD(f1, "proliferation term itself", ev, term = "zp")
}
EP <- dplyr::bind_rows(ep_rows)
print(EP); w_res(EP, "CD03_tcga_alternative_endpoints.csv")

## =====================================================================
##  D2. Recurrence-free survival in GSE14520
## =====================================================================
if (have_gse_extra && sum(!is.na(d14$rec_stat)) >= 50) {
  hdr("D2 - GSE14520 recurrence-free survival")
  dr <- d14[!is.na(d14$rec_stat) & !is.na(d14$rec_time) & d14$rec_time > 0, ]
  r0 <- coxph(Surv(rec_time, rec_stat) ~ z5 + agev + sexv + tnm_bin, data = dr)
  r1 <- coxph(Surv(rec_time, rec_stat) ~ z5 + agev + sexv + tnm_bin + zp, data = dr)
  RFS <- dplyr::bind_rows(
    tidD(r0, "GSE14520 RFS: without proliferation", "recurrence"),
    tidD(r1, "GSE14520 RFS: with proliferation",    "recurrence"),
    tidD(r1, "GSE14520 RFS: proliferation term",    "recurrence", term = "zp"))
  print(RFS); w_res(RFS, "CD04_gse14520_recurrence.csv")
} else {
  message("  [skip] GSE14520 recurrence not available")
}

## =====================================================================
##  D3. Liver function as a covariate
##
##  Every liver variable is recorded for fewer patients than the full model
##  set, so adding one changes the sample as well as the adjustment. The
##  covariate row alone cannot separate the two. Each block therefore fits the
##  SAME model without the covariate on the SAME patients first, and the
##  attenuation that matters is the one between those two rows.
## =====================================================================
hdr("D3 - Does the result survive adjustment for liver function?")
base_terms <- "age + sex + stage_bin + grade_ord + zp"
dd <- CC[complete.cases(CC[, c("CYB5R3","prolif","grade_ord","age","sex","stage_bin","OS","OS.time")]), ]
dd$z5 <- zsafe(dd$CYB5R3); dd$zp <- zsafe(dd$prolif)
dd$cp2  <- factor(ifelse(is.na(dd$lv_cp_class), NA,
                         ifelse(dd$lv_cp_class == "A", "A", "B-C")), c("A","B-C"))
dd$cirr <- factor(ifelse(is.na(dd$lv_ishak_n), NA,
                         ifelse(dd$lv_ishak_n >= 5, "Ishak 5-6", "Ishak 0-4")),
                  c("Ishak 0-4","Ishak 5-6"))

pair_block <- function(label, expr, keepvar, within_rule = TRUE) {
  keep <- which(!is.na(dd[[keepvar]]))
  if (length(keep) < 80) { message("  [skip] ", label, ": too few patients"); return(NULL) }
  d2 <- dd[keep, , drop = FALSE]
  if (sum(d2$OS == 1) < 30) { message("  [skip] ", label, ": too few events"); return(NULL) }
  d2$z5 <- zsafe(d2$CYB5R3); d2$zp <- zsafe(d2$prolif)
  f_ref <- try(coxph(as.formula(paste("Surv(OS.time, OS) ~ z5 +", base_terms)), data = d2), silent = TRUE)
  f_adj <- try(coxph(as.formula(paste("Surv(OS.time, OS) ~ z5 +", base_terms, "+", expr)), data = d2), silent = TRUE)
  if (inherits(f_ref, "try-error") || inherits(f_adj, "try-error")) {
    message("  [skip] ", label, ": model failed"); return(NULL) }
  out <- dplyr::bind_rows(
    tidD(f_ref, paste0("same patients, WITHOUT ", label), "OS"),
    tidD(f_adj, paste0("same patients, WITH ", label),    "OS"))
  out$block <- label
  out$within_coverage_rule <- within_rule
  out
}
SENS <- list(cbind(tidD(coxph(as.formula(paste("Surv(OS.time, OS) ~ z5 +", base_terms)), data = dd),
                        "full model set, reference", "OS"),
                   block = "reference", within_coverage_rule = TRUE))
if (cov_of(dd$cp2) >= RULES_D$min_coverage_to_model)
  SENS[[length(SENS) + 1]] <- pair_block("Child-Pugh (A vs B-C)", "cp2", "cp2")
if (cov_of(dd$lv_albumin_n) >= RULES_D$min_coverage_to_model)
  SENS[[length(SENS) + 1]] <- pair_block("serum albumin", "lv_albumin_n", "lv_albumin_n")
if (cov_of(dd$lv_platelet_n) >= RULES_D$min_coverage_to_model)
  SENS[[length(SENS) + 1]] <- pair_block("platelet count", "lv_platelet_n", "lv_platelet_n")
if (cov_of(dd$lv_ptime_n) >= RULES_D$min_coverage_to_model)
  SENS[[length(SENS) + 1]] <- pair_block("INR", "lv_ptime_n", "lv_ptime_n")
## Ishak fibrosis falls below the pre-fixed coverage threshold. The rule stands,
## so it is not part of the sensitivity set; it is shown separately and labeled,
## because cirrhosis is the variable a clinical reader will ask about first.
if (cov_of(dd$cirr) < RULES_D$min_coverage_to_model) {
  extra <- pair_block("Ishak fibrosis 5-6", "cirr", "cirr", within_rule = FALSE)
  if (!is.null(extra)) SENS[[length(SENS) + 1]] <- extra
}
SENS <- dplyr::bind_rows(Filter(Negate(is.null), SENS))
print(SENS); w_res(SENS, "CD05_liver_function_sensitivity.csv")
message("  Read each block as a pair. A drop from the WITHOUT row to the WITH row")
message("  is attributable to the covariate; a drop from the full-set reference to")
message("  the WITHOUT row is attributable to losing patients.")

if (have_gse_extra) {
  sg <- list(cbind(tidD(coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp, data = d14),
                        "GSE14520 reference", "OS"), block = "reference"))
  if (cov_of(d14$cirr_bin) >= RULES_D$min_coverage_to_model) {
    dg <- d14[!is.na(d14$cirr_bin), ]; dg$cf <- factor(dg$cirr_bin, c("no","yes"))
    sg[[length(sg) + 1]] <- cbind(dplyr::bind_rows(
      tidD(coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp, data = dg),
           "same patients, WITHOUT cirrhosis", "OS"),
      tidD(coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp + cf, data = dg),
           "same patients, WITH cirrhosis", "OS")), block = "cirrhosis")
  }
  if (cov_of(d14$bclc_bin) >= RULES_D$min_coverage_to_model) {
    dg <- d14[!is.na(d14$bclc_bin), ]; dg$bf <- factor(dg$bclc_bin, c("0-A","B-C"))
    sg[[length(sg) + 1]] <- cbind(
      tidD(coxph(Surv(time, status) ~ z5 + agev + sexv + bf + zp, data = dg),
           "BCLC in place of TNM", "OS"), block = "BCLC")
  }
  SG <- dplyr::bind_rows(sg); print(SG); w_res(SG, "CD06_gse14520_liver_sensitivity.csv")
}

## =====================================================================
##  D4. Confidence interval for the CHANGE in the hazard ratio
## =====================================================================
hdr("D4 - Bootstrap interval for the change in the hazard ratio")
delta_boot <- function(d, tv, ev, rhs0, rhs1, B = RULES_D$delta_bootstrap_B, seed = 20260911) {
  f0 <- as.formula(sprintf("Surv(%s, %s) ~ %s", tv, ev, rhs0))
  f1 <- as.formula(sprintf("Surv(%s, %s) ~ %s", tv, ev, rhs1))
  obs0 <- coef(coxph(f0, data = d))[["z5"]]
  obs1 <- coef(coxph(f1, data = d))[["z5"]]
  set.seed(seed)
  n <- nrow(d); out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    db <- d[idx, , drop = FALSE]
    v <- try({
      c0 <- coef(coxph(f0, data = db))[["z5"]]
      c1 <- coef(coxph(f1, data = db))[["z5"]]
      c1 - c0
    }, silent = TRUE)
    if (!inherits(v, "try-error") && is.finite(v)) out[b] <- v
  }
  out <- out[is.finite(out)]
  ci <- as.numeric(quantile(out, c(0.025, 0.975)))
  data.frame(HR_without = exp(obs0), HR_with = exp(obs1),
             ratio_of_HR = exp(obs1 - obs0),
             delta_logHR = obs1 - obs0,
             CI_low = exp(ci[1]), CI_high = exp(ci[2]),
             boot_n = length(out), B = B,
             excludes_no_change = ci[1] > 0,
             stringsAsFactors = FALSE)
}
DB <- delta_boot(dd, "OS.time", "OS",
                 "z5 + age + sex + stage_bin + grade_ord",
                 "z5 + age + sex + stage_bin + grade_ord + zp")
DB$cohort <- "TCGA-LIHC"
DB2 <- delta_boot(d14, "time", "status",
                  "z5 + agev + sexv + tnm_bin",
                  "z5 + agev + sexv + tnm_bin + zp")
DB2$cohort <- "GSE14520"
DELTA <- dplyr::bind_rows(DB, DB2)[, c("cohort","HR_without","HR_with","ratio_of_HR",
                                       "delta_logHR","CI_low","CI_high","boot_n","B",
                                       "excludes_no_change")]
DELTA$reading <- ifelse(DELTA$excludes_no_change,
  "The increase in the hazard ratio is larger than resampling noise.",
  "The increase is not distinguishable from resampling noise; report the point estimate only.")
print(DELTA); w_res(DELTA, "CD07_delta_HR_bootstrap.csv")

## =====================================================================
##  D5. Kaplan-Meier by CYB5R3 tertile, with absolute survival
## =====================================================================
hdr("D5 - Kaplan-Meier and absolute survival")
km_tab <- function(d, tv, ev, cohort, mark) {
  q <- quantile(d$CYB5R3, c(1/3, 2/3), na.rm = TRUE)
  d$tert <- factor(ifelse(d$CYB5R3 <= q[1], "low",
                   ifelse(d$CYB5R3 <= q[2], "middle", "high")),
                   c("low","middle","high"))
  fit <- survfit(as.formula(sprintf("Surv(%s, %s) ~ tert", tv, ev)), data = d)
  sm <- summary(fit)$table
  lab <- sub("^tert=", "", rownames(sm))
  at <- summary(fit, times = mark, extend = TRUE)
  five <- tapply(at$surv, sub("^tert=", "", as.character(at$strata)), function(x) x[1])
  lr <- survdiff(as.formula(sprintf("Surv(%s, %s) ~ tert", tv, ev)), data = d)
  p  <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  list(tab = data.frame(cohort = cohort, tertile = lab,
                        n = as.numeric(sm[,"records"]), events = as.numeric(sm[,"events"]),
                        median = as.numeric(sm[,"median"]),
                        survival_at_mark = round(as.numeric(five[lab]), 3),
                        mark = mark, logrank_p = p, stringsAsFactors = FALSE),
       fit = fit, data = d)
}
K1 <- km_tab(CC, "OS.time", "OS", "TCGA-LIHC", 1825)
K2 <- km_tab(d14, "time", "status", "GSE14520", 60)
KM <- dplyr::bind_rows(K1$tab, K2$tab)
print(KM); w_res(KM, "CD08_km_absolute_survival.csv")

## curves, drawn from the fitted objects so the figure cannot drift from the table
km_df <- function(K, cohort, xlab_unit) {
  f <- K$fit
  data.frame(time = f$time, surv = f$surv,
             tert = factor(rep(sub("^tert=", "", names(f$strata)), f$strata),
                           c("low","middle","high")),
             cohort = cohort, stringsAsFactors = FALSE)
}
KD <- rbind(km_df(K1, "TCGA-LIHC (days)"), km_df(K2, "GSE14520 (months)"))
KD$cohort <- factor(KD$cohort, c("TCGA-LIHC (days)", "GSE14520 (months)"))

if (!exists("th")) th <- ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold", size = 8.5),
                 plot.subtitle = ggplot2::element_text(size = 7))

## Figure saving is self-contained so that this stage does not depend on
## 22_figures.R having been sourced first. The TIFF compression tag is read
## back from the written file, because on macOS without cairo the tiff device
## falls back to quartz and writes an uncompressed file while reporting lzw.
if (!exists("save3")) {
  tiff_compression <- function(fn) {
    con <- file(fn, "rb"); on.exit(close(con), add = TRUE)
    h <- readBin(con, "raw", 8)
    e <- if (rawToChar(h[1:2]) == "II") "little" else "big"
    off <- readBin(h[5:8], "integer", size = 4, endian = e)
    seek(con, off)
    n <- readBin(con, "integer", size = 2, n = 1, signed = FALSE, endian = e)
    for (i in seq_len(n)) {
      ent <- readBin(con, "raw", 12)
      tag <- readBin(ent[1:2], "integer", size = 2, signed = FALSE, endian = e)
      if (tag == 259L) return(readBin(ent[9:10], "integer", size = 2, signed = FALSE, endian = e))
    }
    NA_integer_
  }
  save3 <- function(p, name, w = 175, h = 126) {
    base <- file.path(FIG, name); win <- w / 25.4; hin <- h / 25.4
    ggplot2::ggsave(paste0(base, ".png"), p, width = win, height = hin, units = "in", dpi = 300)
    tif <- paste0(base, ".tiff"); ok <- FALSE
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_tiff(tif, width = win, height = hin, units = "in", res = 300,
                     compression = "lzw", background = "white")
      print(p); grDevices::dev.off()
      ok <- isTRUE(tiff_compression(tif) == 5L)
    }
    if (!ok) {
      ggplot2::ggsave(tif, p, width = win, height = hin, units = "in", dpi = 300)
      message("  [!] ", basename(tif), " is UNCOMPRESSED. install.packages(\"ragg\") for LZW.")
    }
    svg <- tryCatch({ suppressMessages(ggplot2::ggsave(paste0(base, ".svg"), p,
                        width = win, height = hin, units = "in")); TRUE },
                    error = function(e) FALSE)
    eps <- paste0(base, ".eps"); epsok <- FALSE
    if (file.exists(eps)) unlink(eps)   # cairo_ps() only warns; see 22_figures.R
    if (isTRUE(capabilities("cairo"))) {
      epsok <- tryCatch({
        grDevices::cairo_ps(eps, width = win, height = hin, fallback_resolution = 800,
                            onefile = FALSE, bg = "white")
        print(p); grDevices::dev.off(); file.exists(eps) && file.size(eps) > 5000
      }, error = function(e) FALSE)
    }
    if (!epsok) {
      epsok <- tryCatch({
        grDevices::postscript(eps, width = win, height = hin, onefile = FALSE,
                              horizontal = FALSE, paper = "special", bg = "white")
        print(p); grDevices::dev.off(); file.exists(eps) && file.size(eps) > 5000
      }, error = function(e) FALSE)
    }
    message("  -> ", basename(base), ": png, tiff (",
            if (ok) "LZW verified" else "uncompressed", ")", if (svg) ", svg" else "",
            if (epsok) ", eps" else " [!] EPS FAILED")
  }
}

pk <- ggplot2::ggplot(KD, ggplot2::aes(time, surv, colour = tert)) +
  ggplot2::geom_step(linewidth = 0.5) +
  ggplot2::facet_wrap(~ cohort, scales = "free_x") +
  ggplot2::scale_colour_manual(values = c(low = "#2E5496", middle = "#7F7F7F", high = "#C00000")) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = "Time since diagnosis", y = "Overall survival",
                colour = "CYB5R3 tertile",
                title = "Overall survival by CYB5R3 tertile",
                subtitle = "Unadjusted. The crude curves are the point: without adjustment the groups barely separate.") +
  th + ggplot2::theme(legend.position = "bottom")
save3(pk, "Figure5_km", w = 175, h = 87)

hdr("D - done")
message("  All Stage D output is POST HOC. Files: CD00 to CD08.")
if (exists("save_session")) save_session("23_clinical")
