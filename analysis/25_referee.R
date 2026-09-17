# =====================================================================
#  25_referee.R  -  Project C, Stage E: analyses added in referee review
#
#  ★★ POST HOC. Designed on 2026-09-17, after the manuscript had been
#     drafted, audited and rewritten for a second journal. Nothing here is
#     pre-registered. Every result is reported as post hoc.
#
#  Why it exists. A reviewer-perspective audit raised five gaps:
#    E1  The Cox hazard ratio is NOT collapsible. Adding a strongly
#        prognostic covariate raises the conditional hazard ratio even
#        when that covariate is independent of the exposure and no
#        confounding exists. Part of 1.241 -> 1.371 may therefore be
#        arithmetic rather than confounding. E1 measures how much.
#    E2  No baseline characteristics table, no median follow-up, no
#        censoring account. REMARK and STROBE both require them.
#    E3  The grade trend rests on 12 grade-4 tumors. E3 removes them.
#    E4  Tissue composition was carried by PC1 of marker genes rather
#        than by a published purity estimate. E4 uses ESTIMATE.
#    E5  The proliferation score is eight hand-picked genes. E5 repeats
#        the central contrast with HALLMARK_G2M_CHECKPOINT instead.
#
#  Run: open RUN_C_CLINICAL.R and Source it first (it runs 18, 19, 20,
#  23 and leaves CC, d14 and EX in the workspace), then Source this file.
# =====================================================================
if (!exists("CC") || !exists("d14"))
  stop("Run RUN_C_CLINICAL.R first; CC and d14 must be in the workspace.")

suppressPackageStartupMessages({ library(survival) })

## CC carries CYB5R3 and prolif but not the z-scores that the models use;
## 23_clinical.R builds them inside each block. Do the same, once.
zf <- if (exists("zsafe")) zsafe else function(x) as.numeric(scale(x))
CCE <- CC[stats::complete.cases(
            CC[, c("CYB5R3", "prolif", "grade_ord", "age", "sex", "stage_bin")]), ]
CCE <- CCE[!is.na(CCE$OS.time) & !is.na(CCE$OS) & CCE$OS.time > 0, ]
CCE$z5 <- zf(CCE$CYB5R3); CCE$zp <- zf(CCE$prolif)
message("  TCGA analysis frame: ", nrow(CCE), " patients, ", sum(CCE$OS), " deaths")
if (!nrow(CCE)) stop("the TCGA analysis frame is empty; check that CC has CYB5R3, prolif, grade_ord, age, sex and stage_bin")
hdr("E - POST HOC analyses added in referee review")
message("  script version: 2026-09-17e")
message("  NOTE: every result below is post hoc and must be labeled as such.")

RULES_E <- list(
  designed_on   = "2026-09-17, after the manuscript was rewritten",
  status        = "POST HOC, sensitivity and reporting only",
  e1_question   = "how much of the rise in the hazard ratio is non-collapsibility",
  e1_design     = paste("simulate survival from the fitted adjusted model with the",
                        "proliferation score permuted, so that it keeps its marginal",
                        "distribution and its coefficient but is independent of CYB5R3;",
                        "the ratio of hazard ratios that remains is non-collapsibility alone"),
  e1_B          = 1000,
  e1_seed       = 20260917,
  e3_question   = "does the grade trend survive removal of the 12 grade-4 tumors",
  e4_tool       = "ESTIMATE stromal and immune scores; on RNA-seq the calibrated purity conversion is not defined, so the ESTIMATE score itself is used as the covariate",
  e5_signature  = "HALLMARK_G2M_CHECKPOINT via msigdbr, if available"
)
w_res(data.frame(rule = names(RULES_E), value = unlist(RULES_E)), "CE00_referee_rules.csv")

## =====================================================================
##  E1 - How much of the rise is non-collapsibility?
## =====================================================================
hdr("E1 - non-collapsibility of the Cox hazard ratio")

## Simulate survival times from a fitted Cox model, by inverting the
## fitted baseline cumulative hazard. Returns time and status.
sim_cox <- function(fit, lp, tmax, cens_fit) {
  bh <- survival::basehaz(fit, centered = FALSE)      # H0(t)
  bh <- bh[!duplicated(bh$hazard), , drop = FALSE]    # approxfun needs distinct x
  inv <- stats::approxfun(bh$hazard, bh$time, method = "linear",
                          yleft = 0, rule = 2)
  n <- length(lp)
  target <- -log(stats::runif(n)) / exp(lp)
  t_ev <- inv(target)
  t_ev[target > max(bh$hazard)] <- tmax               # beyond follow-up
  ## draw censoring from the observed censoring distribution
  cs <- summary(cens_fit, times = sort(unique(c(0, bh$time))), extend = TRUE)
  ci <- stats::approxfun(cs$time, cs$surv, method = "constant", yleft = 1, rule = 2)
  u  <- stats::runif(n)
  grid <- sort(unique(c(0, bh$time, tmax)))
  sc   <- ci(grid)
  t_ce <- vapply(u, function(uu) {
    k <- which(sc <= uu)[1]
    if (is.na(k)) tmax else grid[k]
  }, numeric(1))
  data.frame(time = pmin(t_ev, t_ce), status = as.integer(t_ev <= t_ce))
}

noncollapse <- function(d, tv, ev, rhs_base, expo = "z5", zvar = "zp",
                        B = RULES_E$e1_B, label = "") {
  f_base <- stats::as.formula(sprintf("Surv(%s, %s) ~ %s", tv, ev, rhs_base))
  f_full <- stats::as.formula(sprintf("Surv(%s, %s) ~ %s + %s", tv, ev, rhs_base, zvar))
  d <- d[stats::complete.cases(d[, all.vars(f_full), drop = FALSE]), ]
  d <- d[d[[tv]] > 0, ]
  if (!nrow(d)) stop("no complete observations for ", label)
  m_base <- survival::coxph(f_base, data = d)
  m_full <- survival::coxph(f_full, data = d)
  obs <- exp(unname(coef(m_full)[expo]) - unname(coef(m_base)[expo]))

  ## the truth we simulate under: the fitted conditional model, but with
  ## the proliferation score permuted so it is independent of CYB5R3
  X    <- stats::model.matrix(f_full, data = d)
  if ("(Intercept)" %in% colnames(X)) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  cf   <- coef(m_full)
  if (!all(names(cf) %in% colnames(X)))
    stop("design matrix does not match the fitted coefficients for ", label)
  X    <- X[, names(cf), drop = FALSE]
  tmax <- max(d[[tv]], na.rm = TRUE)
  cens <- survival::survfit(stats::as.formula(sprintf("Surv(%s, 1 - %s) ~ 1", tv, ev)),
                            data = d)
  set.seed(RULES_E$e1_seed)
  out <- numeric(B)
  for (b in seq_len(B)) {
    Xb <- X
    Xb[, zvar] <- sample(X[, zvar])          # independent of everything else
    lp <- as.vector(Xb %*% cf)
    sdat <- sim_cox(m_full, lp, tmax, cens)
    db <- cbind(d[, setdiff(names(d), c(tv, ev, "time", "status")), drop = FALSE], sdat)
    db[[zvar]] <- Xb[, zvar]
    fb <- try(survival::coxph(stats::as.formula(
                sprintf("Surv(time, status) ~ %s", rhs_base)), data = db), silent = TRUE)
    ff <- try(survival::coxph(stats::as.formula(
                sprintf("Surv(time, status) ~ %s + %s", rhs_base, zvar)), data = db), silent = TRUE)
    out[b] <- if (inherits(fb, "try-error") || inherits(ff, "try-error")) NA_real_
              else exp(unname(coef(ff)[expo]) - unname(coef(fb)[expo]))
  }
  q <- stats::quantile(out, c(0.025, 0.5, 0.975), na.rm = TRUE)
  ## Share of the OBSERVED RISE that non-collapsibility explains, on the log
  ## scale. An earlier version divided one ratio by the other, which is not a
  ## share of anything; that column was wrong and has been replaced.
  data.frame(cohort = label,
             observed_ratio_of_HR = obs,
             observed_rise_percent = 100 * (obs - 1),
             noncollapsibility_median = unname(q[2]),
             noncollapsibility_rise_percent = 100 * (unname(q[2]) - 1),
             noncollapsibility_lo = unname(q[1]),
             noncollapsibility_hi = unname(q[3]),
             share_of_rise_from_noncollapsibility = log(unname(q[2])) / log(obs),
             p_one_sided_obs_vs_noncollapsibility = mean(out >= obs, na.rm = TRUE),
             B = sum(!is.na(out)))
}

NC <- rbind(
  noncollapse(CCE, "OS.time", "OS", "z5 + age + sex + stage_bin + grade_ord",
              label = "TCGA-LIHC"),
  noncollapse(d14, "time", "status",
              "z5 + agev + sexv + tnm_bin", label = "GSE14520"))
print(NC); w_res(NC, "CE01_noncollapsibility.csv")
message("  Reading: 'share_of_rise_from_noncollapsibility' is the fraction of the",
        " observed rise that non-collapsibility alone would produce.")

## =====================================================================
##  E2 - Baseline characteristics, follow-up and censoring
## =====================================================================
hdr("E2 - baseline characteristics and follow-up")

med_fu <- function(d, tv, ev) {                     # reverse Kaplan-Meier
  f <- try(survival::survfit(stats::as.formula(sprintf("Surv(%s, 1 - %s) ~ 1", tv, ev)),
                             data = d), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  tb <- summary(f)$table
  v <- if (is.matrix(tb)) tb[1, "median"] else tb["median"]
  unname(v)
}
num_row <- function(lab, x, unit) data.frame(
  characteristic = lab, value = sprintf("%.0f (%.0f to %.0f)", stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE), stats::quantile(x, 0.75, na.rm = TRUE)),
  unit = unit, n_available = sum(!is.na(x)))
cat_rows <- function(lab, x) {
  t <- table(x, useNA = "no")
  data.frame(characteristic = paste0(lab, ": ", names(t)),
             value = sprintf("%d (%.1f%%)", as.integer(t), 100 * as.integer(t) / sum(t)),
             unit = "n (%)", n_available = sum(t))
}
B1 <- rbind(
  data.frame(characteristic = "Patients", value = as.character(nrow(CC)), unit = "n",
             n_available = nrow(CC)),
  data.frame(characteristic = "Deaths", value = as.character(sum(CC$OS, na.rm = TRUE)),
             unit = "n", n_available = nrow(CC)),
  data.frame(characteristic = "Censored", value = as.character(sum(CC$OS == 0, na.rm = TRUE)),
             unit = "n", n_available = nrow(CC)),
  num_row("Age", CC$age, "years"),
  cat_rows("Sex", CC$sex),
  cat_rows("Stage", CC$stage_bin),
  cat_rows("Grade", CC$grade_ord),
  data.frame(characteristic = "Median follow-up (reverse Kaplan-Meier)",
             value = sprintf("%.0f", med_fu(CC, "OS.time", "OS")), unit = "days",
             n_available = nrow(CC)))
B1$cohort <- "TCGA-LIHC"
B2 <- rbind(
  data.frame(characteristic = "Patients", value = as.character(nrow(d14)), unit = "n",
             n_available = nrow(d14)),
  data.frame(characteristic = "Deaths", value = as.character(sum(d14$status, na.rm = TRUE)),
             unit = "n", n_available = nrow(d14)),
  data.frame(characteristic = "Censored", value = as.character(sum(d14$status == 0, na.rm = TRUE)),
             unit = "n", n_available = nrow(d14)),
  num_row("Age", d14$agev, "years"),
  cat_rows("Sex", d14$sexv),
  cat_rows("TNM stage", d14$tnm_bin),
  data.frame(characteristic = "Median follow-up (reverse Kaplan-Meier)",
             value = sprintf("%.1f", med_fu(d14, "time", "status")), unit = "months",
             n_available = nrow(d14)))
if (!is.null(d14$cirr_bin)) B2 <- rbind(B2, cat_rows("Cirrhosis", d14$cirr_bin))
if (!is.null(d14$bclc_bin)) B2 <- rbind(B2, cat_rows("BCLC stage", d14$bclc_bin))
B2$cohort <- "GSE14520"
BASE <- rbind(B1, B2)[, c("cohort", "characteristic", "value", "unit", "n_available")]
print(BASE, row.names = FALSE); w_res(BASE, "CE02_baseline_characteristics.csv")

## =====================================================================
##  E3 - Does the grade trend survive removal of grade 4?
## =====================================================================
hdr("E3 - grade trend without the 12 grade-4 tumors")

grade_block <- function(d, lab) {
  d <- d[!is.na(d$grade_ord) & !is.na(d$z5), ]
  jt <- try(clinfun::jonckheere.test(d$z5, d$grade_ord, alternative = "two.sided",
                                     nperm = 10000), silent = TRUE)
  rho <- suppressWarnings(stats::cor.test(d$z5, as.numeric(d$grade_ord),
                                          method = "spearman"))
  sl_none <- stats::lm(rank(z5) ~ grade_ord + age + sex + stage_bin, data = d)
  sl_prol <- stats::lm(rank(z5) ~ grade_ord + age + sex + stage_bin + zp, data = d)
  b0 <- unname(coef(sl_none)["grade_ord"]); b1 <- unname(coef(sl_prol)["grade_ord"])
  data.frame(set = lab, n = nrow(d),
             JT_p = if (inherits(jt, "try-error")) NA_real_ else unname(jt$p.value),
             spearman_rho = unname(rho$estimate), spearman_p = unname(rho$p.value),
             slope_unadjusted = b0, p_unadjusted = summary(sl_none)$coefficients["grade_ord", 4],
             slope_with_proliferation = b1,
             p_with_proliferation = summary(sl_prol)$coefficients["grade_ord", 4],
             attenuation = 1 - abs(b1) / abs(b0))
}
CG <- CCE
CG34 <- CG; CG34$grade_ord <- pmin(as.numeric(CG34$grade_ord), 3)
GR <- rbind(grade_block(CG, "all grades"),
            grade_block(CG[as.numeric(CG$grade_ord) < 4, ], "grade 4 excluded"),
            grade_block(CG34, "grades 3 and 4 collapsed"))
print(GR); w_res(GR, "CE03_grade_trend_sensitivity.csv")

PROG <- do.call(rbind, lapply(list(
  list(d = CG, lab = "all grades"),
  list(d = CG[as.numeric(CG$grade_ord) < 4, ], lab = "grade 4 excluded")), function(s) {
    m <- survival::coxph(Surv(OS.time, OS) ~ grade_ord + age + sex + stage_bin, data = s$d)
    ci <- summary(m)$conf.int["grade_ord", ]
    data.frame(set = s$lab, n = nrow(s$d), events = sum(s$d$OS),
               HR_grade = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
               p = summary(m)$coefficients["grade_ord", 5])
  }))
print(PROG); w_res(PROG, "CE04_grade_prognostic_sensitivity.csv")

## =====================================================================
##  E4 - Tumor purity by ESTIMATE, in place of the PC1 proxy
## =====================================================================
hdr("E4 - ESTIMATE stromal, immune and combined scores")

has_est <- requireNamespace("estimate", quietly = TRUE)
if (!has_est) {
  message("  installing 'estimate' from R-Forge ...")
  try(utils::install.packages("estimate", repos = "http://r-forge.r-project.org",
                              dependencies = TRUE), silent = TRUE)
  has_est <- requireNamespace("estimate", quietly = TRUE)
}
if (!has_est || !exists("EX")) {
  message("  [!] ESTIMATE not available (or EX not in the workspace). E4 SKIPPED.")
  w_res(data.frame(note = "ESTIMATE unavailable; E4 not run"), "CE05_estimate_purity.csv")
} else { e4 <- try({
  library(estimate)   # attaching loads common_genes and SI_geneset
  tf <- tempfile(fileext = ".txt"); gf <- tempfile(fileext = ".gct")
  sf <- tempfile(fileext = ".gct"); of <- tempfile(fileext = ".gct")
  utils::write.table(cbind(NAME = rownames(EX), EX), tf, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = tf, output.f = gf, id = "GeneSymbol")
  estimateScore(gf, sf, platform = "illumina")
  es <- utils::read.delim(sf, skip = 2, row.names = 1, check.names = FALSE)
  es <- as.data.frame(t(es[, -1, drop = FALSE]))
  es$patient <- substr(gsub("\\.", "-", rownames(es)), 1, 12)
  CP <- merge(CCE, es[, c("patient", "StromalScore", "ImmuneScore", "ESTIMATEScore")],
              by = "patient")
  CP$zpur <- zsc(CP$ESTIMATEScore)
  CP <- CP[!is.na(CP$grade_ord) & !is.na(CP$zp) & !is.na(CP$zpur), ]
  rows <- list(
    c("without the ESTIMATE score", "z5 + age + sex + stage_bin + grade_ord + zp"),
    c("with the ESTIMATE score", "z5 + age + sex + stage_bin + grade_ord + zp + zpur"),
    c("ESTIMATE score but no proliferation", "z5 + age + sex + stage_bin + grade_ord + zpur"))
  PUR <- do.call(rbind, lapply(rows, function(r) {
    m <- survival::coxph(stats::as.formula(paste("Surv(OS.time, OS) ~", r[2])), data = CP)
    ci <- summary(m)$conf.int["z5", ]
    data.frame(model = r[1], n = nrow(CP), events = sum(CP$OS),
               HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
               p = summary(m)$coefficients["z5", 5])
  }))
  ## the quantity the article cites: is the composition score measuring the
  ## same thing as the proliferation score, or something different?
  PUR$rho_ESTIMATE_vs_proliferation <- suppressWarnings(
    stats::cor(CP$zpur, CP$zp, method = "spearman", use = "complete.obs"))
  print(PUR); w_res(PUR, "CE05_estimate_purity.csv")
}, silent = TRUE)
  if (inherits(e4, "try-error")) {
    message("  [!] E4 failed, continuing: ", as.character(e4))
    w_res(data.frame(note = paste("E4 failed:", as.character(e4))), "CE05_estimate_purity.csv")
  }
}

## =====================================================================
##  E5 - An independently defined proliferation signature
## =====================================================================
hdr("E5 - HALLMARK_G2M_CHECKPOINT in place of the eight-gene score")

has_ms <- requireNamespace("msigdbr", quietly = TRUE)
if (!has_ms) {
  message("  installing 'msigdbr' ...")
  try(utils::install.packages("msigdbr"), silent = TRUE)
  has_ms <- requireNamespace("msigdbr", quietly = TRUE)
}
if (!has_ms || !exists("EX")) {
  message("  [!] msigdbr not available (or EX not in the workspace). E5 SKIPPED.")
  w_res(data.frame(note = "msigdbr unavailable; E5 not run"), "CE06_hallmark_g2m.csv")
} else { e5 <- try({
  h <- try(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"), silent = TRUE)
  if (inherits(h, "try-error") || !nrow(h))
    h <- try(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), silent = TRUE)
  if (inherits(h, "try-error")) stop("msigdbr returned no collection")
  gs <- unique(h$gene_symbol[h$gs_name == "HALLMARK_G2M_CHECKPOINT"])
  tmpl <- if (exists("TEMPL")) unique(unlist(TEMPL)) else character(0)
  gs <- setdiff(intersect(gs, rownames(EX)), tmpl)                    # template independent
  message("  HALLMARK_G2M_CHECKPOINT genes used: ", length(gs))
  g2m <- colMeans(t(scale(t(EX[gs, , drop = FALSE]))), na.rm = TRUE)
  gd <- data.frame(patient = substr(gsub("\\.", "-", names(g2m)), 1, 12),
                   g2m = as.numeric(g2m))
  CH <- merge(CCE, gd, by = "patient")
  CH$zg <- zsc(CH$g2m)
  CH <- CH[!is.na(CH$grade_ord) & !is.na(CH$zg) & !is.na(CH$z5), ]
  rows <- list(
    c("without any proliferation term", "z5 + age + sex + stage_bin + grade_ord"),
    c("with the eight-gene score",      "z5 + age + sex + stage_bin + grade_ord + zp"),
    c("with HALLMARK G2M",              "z5 + age + sex + stage_bin + grade_ord + zg"))
  G2 <- do.call(rbind, lapply(rows, function(r) {
    m <- survival::coxph(stats::as.formula(paste("Surv(OS.time, OS) ~", r[2])), data = CH)
    ci <- summary(m)$conf.int["z5", ]
    data.frame(model = r[1], n = nrow(CH), events = sum(CH$OS),
               HR = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]),
               p = summary(m)$coefficients["z5", 5])
  }))
  G2$n_genes_G2M <- length(gs)
  G2$rho_G2M_vs_eight_gene <- suppressWarnings(
    stats::cor(CH$zg, CH$zp, method = "spearman", use = "complete.obs"))
  print(G2); w_res(G2, "CE06_hallmark_g2m.csv")
}, silent = TRUE)
  if (inherits(e5, "try-error")) {
    message("  [!] E5 failed, continuing: ", as.character(e5))
    w_res(data.frame(note = paste("E5 failed:", as.character(e5))), "CE06_hallmark_g2m.csv")
  }
}

hdr("E - done")
message("  All Stage E output is POST HOC. Files: CE00 to CE06.")
if (exists("save_session")) save_session("25_referee")
