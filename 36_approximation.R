# =====================================================================
#  36_approximation.R  -  Stage O: where the first-order omitted
#                          variable expression fails in real Cox data
#
#  ★★ POST HOC. Designed on 2026-09-29, after the registered test had
#  been seen and after a referee objected that the landscape slope is
#  close to an algebraic consequence of the definitions. The rules
#  below were fixed before this script was first run and are written
#  to CO00 on every run.
#
#  The objection is correct as far as it goes. In a linear model,
#  omitting Z from y ~ X b + Z g gives
#
#        E[b_omit] = b + delta * g
#
#  with delta the coefficient of Z regressed on X. Applied to one gene
#  at a time, with the shift defined as the adjusted coefficient minus
#  the unadjusted one, this predicts
#
#        shift_g  ~=  - gamma * delta_g
#
#  where gamma is the log hazard ratio of the proliferation score and
#  delta_g is the coefficient of the standardized gene in an auxiliary
#  regression of the score on that gene and the other covariates.
#
#  That expression is exact for linear models and approximate for Cox.
#  What this script measures is the approximation error: how much of
#  the observed shift the expression accounts for, how large what is
#  left over is, and whether what is left over is the non-collapsibility
#  term that Stage G estimated by a different route. The error, not the
#  slope, is the part of this that the definitions do not give away.
#
#  Requires: GW, CCF, Mt, Zb from 26_genomewide.R (and, if present,
#  GG, dd, Mg, Zg for the replication cohort).
#
#  Writes CO00 to CO06.
# =====================================================================

if (!exists("hdr") || !exists("w_res"))
  stop("Run RUN_C_ALL.R first; hdr and w_res must exist.")
if (!exists("GW") || !exists("CCF") || !exists("Mt") || !exists("Zb"))
  stop("Run analysis/26_genomewide.R first; GW, CCF, Mt and Zb must exist.")

hdr("O - POST HOC: where the first-order approximation fails")
message("  script version: 2026-09-29a")
message("  NOTE: every result below is post hoc and must be labeled as such.")

RULES <- data.frame(
  rule = c("status",
           "designed_on",
           "question",
           "gamma",
           "delta",
           "prediction",
           "residual",
           "primary_reading",
           "secondary_reading",
           "failure_reading",
           "bins",
           "cohorts"),
  value = c(
    "POST HOC, designed 2026-09-29 after the registered test was seen",
    "2026-09-29",
    paste("how much of the per-gene shift does the first-order omitted",
          "variable expression account for, and what is left over"),
    paste("the log hazard ratio of the standardized proliferation score in",
          "a Cox model containing the score and the covariates but no gene"),
    paste("the coefficient of the standardized gene in an ordinary least",
          "squares regression of the standardized score on that gene and",
          "the same covariates, obtained by the Frisch-Waugh-Lovell",
          "residual form so that it is exact rather than iterative"),
    "predicted log shift = -gamma * delta",
    "residual = observed log shift - predicted log shift",
    paste("if the predicted value explains at least 0.95 of the variance of",
          "the observed shift and the regression of observed on predicted",
          "has a slope within 0.05 of 1, the first-order expression is an",
          "adequate description and the departure is second order"),
    paste("if the residual is proportional to the gene's own log hazard",
          "ratio with a coefficient close to the non-collapsibility",
          "inflation estimated independently in Stage G (0.064), then what",
          "the expression omits is non-collapsibility and the two analyses",
          "agree by different routes"),
    paste("if neither holds, the departure is not explained by either",
          "account and is reported as unexplained"),
    "deciles of the predicted shift and of the gene's own log hazard ratio",
    "TCGA-LIHC, and GSE14520 where the objects from Stage F are present"),
  stringsAsFactors = FALSE)
print(RULES); w_res(RULES, "CO00_approximation_rules.csv")

## ---------------------------------------------------------------
##  a reusable worker, so the two cohorts are treated identically
## ---------------------------------------------------------------
approx_cohort <- function(M, frame, ycol, ecol, Zbase, zp, sweep, label) {
  ## gamma: the score's own coefficient, with no gene in the model
  dg <- data.frame(zp = zp, Zbase)
  fg <- survival::coxph(survival::Surv(frame[[ycol]], frame[[ecol]]) ~ .,
                        data = dg)
  gamma <- unname(stats::coef(fg)["zp"])

  ## delta, for every gene at once, by Frisch-Waugh-Lovell:
  ## project the covariates out of both sides, then take the ratio
  X   <- stats::model.matrix(~ ., data = Zbase)
  qrX <- qr(X)
  G   <- scale(t(as.matrix(M)))          # patients x genes, standardized
  Eg  <- qr.resid(qrX, G)
  ez  <- as.numeric(qr.resid(qrX, matrix(zp, ncol = 1)))
  ss  <- colSums(Eg^2)
  delta <- as.numeric(crossprod(Eg, ez)) / ss
  names(delta) <- colnames(G)

  D <- data.frame(gene = names(delta), delta = delta,
                  stringsAsFactors = FALSE)
  D <- merge(D, sweep[, c("gene", "HR_without", "log_shift")], by = "gene")
  D$coef_without <- log(D$HR_without)
  D$predicted    <- -gamma * D$delta
  D$observed     <- D$log_shift
  D$residual     <- D$observed - D$predicted
  D <- D[is.finite(D$predicted) & is.finite(D$observed) &
         is.finite(D$coef_without) & is.finite(D$delta), ]
  attr(D, "gamma") <- gamma
  attr(D, "label") <- label
  D
}

## ---------------------------------------------------------------
##  O1 - the discovery cohort
## ---------------------------------------------------------------
hdr("O1 - the approximation in TCGA-LIHC")
DO <- approx_cohort(Mt, CCF, "OS.time", "OS", Zb, CCF$zp, GW, "TCGA-LIHC")
message("  genes matched: ", nrow(DO),
        " ; gamma = ", round(attr(DO, "gamma"), 4))
w_res(DO[order(-abs(DO$residual)),
         c("gene", "delta", "coef_without", "predicted", "observed",
           "residual")],
      "CO01_tcga_approximation.csv")

summarise_fit <- function(D) {
  m1 <- stats::lm(observed ~ predicted, data = D)
  m2 <- stats::lm(observed ~ predicted + coef_without, data = D)
  mr <- stats::lm(residual ~ 0 + coef_without, data = D)
  cir <- stats::confint(mr)
  data.frame(
    quantity = c("genes",
                 "gamma",
                 "correlation of observed with predicted",
                 "variance of observed explained by predicted",
                 "slope of observed on predicted",
                 "intercept of observed on predicted",
                 "median residual",
                 "interquartile range of the residual",
                 "median absolute residual",
                 "median absolute observed shift",
                 "median absolute residual as a share of the observed",
                 "variance explained after adding the gene's own coefficient",
                 "residual per unit of the gene's own log hazard ratio",
                 "lower confidence limit",
                 "upper confidence limit",
                 "Stage G non-collapsibility inflation minus one"),
    value = c(nrow(D),
              round(attr(D, "gamma"), 4),
              round(stats::cor(D$observed, D$predicted), 4),
              round(summary(m1)$r.squared, 4),
              round(unname(stats::coef(m1)[2]), 4),
              round(unname(stats::coef(m1)[1]), 5),
              round(stats::median(D$residual), 5),
              round(stats::IQR(D$residual), 5),
              round(stats::median(abs(D$residual)), 5),
              round(stats::median(abs(D$observed)), 5),
              round(stats::median(abs(D$residual)) /
                    stats::median(abs(D$observed)), 4),
              round(summary(m2)$r.squared, 4),
              round(unname(stats::coef(mr)[1]), 4),
              round(unname(cir[1, 1]), 4),
              round(unname(cir[1, 2]), 4),
              0.0643),
    stringsAsFactors = FALSE)
}
SO <- summarise_fit(DO)
print(SO); w_res(SO, "CO02_tcga_approximation_summary.csv")

## ---------------------------------------------------------------
##  O2 - where it breaks down
## ---------------------------------------------------------------
hdr("O2 - where the approximation breaks down")
by_decile <- function(D, var, lab) {
  b <- stats::quantile(D[[var]], probs = seq(0, 1, 0.1), na.rm = TRUE)
  b[1] <- b[1] - 1e-9; b[length(b)] <- b[length(b)] + 1e-9
  D$bin <- cut(D[[var]], breaks = unique(b), include.lowest = TRUE)
  do.call(rbind, lapply(split(D, D$bin), function(d) data.frame(
    grouping = lab,
    bin = as.character(d$bin[1]),
    n = nrow(d),
    median_value = round(stats::median(d[[var]]), 4),
    median_abs_residual = round(stats::median(abs(d$residual)), 5),
    median_abs_observed = round(stats::median(abs(d$observed)), 5),
    share = round(stats::median(abs(d$residual)) /
                  stats::median(abs(d$observed)), 4),
    stringsAsFactors = FALSE)))
}
BD <- rbind(by_decile(DO, "predicted", "decile of the predicted shift"),
            by_decile(DO, "coef_without", "decile of the gene's own coefficient"))
print(BD); w_res(BD, "CO03_breakdown_by_decile.csv")

w_res(utils::head(DO[order(-abs(DO$residual)),
                     c("gene", "delta", "coef_without", "predicted",
                       "observed", "residual")], 50),
      "CO04_largest_residuals.csv")

## ---------------------------------------------------------------
##  O3 - the replication cohort, if Stage F left it behind
## ---------------------------------------------------------------
if (exists("GG") && exists("dd") && exists("Mg") && exists("Zg") &&
    is.data.frame(GG) && nrow(GG) > 0) {
  hdr("O3 - the approximation in GSE14520")
  DG <- approx_cohort(Mg, dd, "time", "status", Zg, dd$zp, GG, "GSE14520")
  message("  genes matched: ", nrow(DG),
          " ; gamma = ", round(attr(DG, "gamma"), 4))
  w_res(DG[order(-abs(DG$residual)),
           c("gene", "delta", "coef_without", "predicted", "observed",
             "residual")],
        "CO05_gse_approximation.csv")
  SG <- summarise_fit(DG)
  print(SG); w_res(SG, "CO06_gse_approximation_summary.csv")
} else {
  message("  [!] GSE14520 objects not in the workspace; O3 skipped.")
  w_res(data.frame(note = "GSE14520 objects absent; O3 not run"),
        "CO06_gse_approximation_summary.csv")
}

hdr("O - done")
message("  All Stage O output is POST HOC. Files: CO00 to CO06.")
if (exists("save_session")) save_session("36_approximation")
