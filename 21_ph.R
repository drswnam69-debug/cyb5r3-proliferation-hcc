# =====================================================================
#  21_ph.R  —  Project C, Addendum C
#  PRE-SPECIFIED handling of the proportional-hazards violation found in B6.
#
#  What B6 showed: in the TCGA model, the proliferation term violates PH
#  (Schoenfeld p = 0.0023; global p = 0.0094). CYB5R3 itself is fine
#  (p = 0.808). So the proliferation effect is time-dependent, and a Cox
#  model that forces it to be constant is misspecified.
#
#  The question this script settles, fixed before running:
#    Does the CYB5R3 estimate survive once proliferation is handled without
#    a proportional-hazards assumption?
#
#  Pre-specified reading (do not edit after seeing output):
#    ROBUST   if CYB5R3 HR >= 1.25 with a 95% CI excluding 1 under BOTH
#             stratification and the time-varying specification
#    WEAKENED if the CI covers 1 under either
#    Anything else is reported as indeterminate.
# =====================================================================
hdr("C - PRE-SPECIFIED handling of the PH violation")

RULES_C <- list(
  robust_if   = "HR >= 1.25 and 95% CI excludes 1 under BOTH strata() and tt()",
  weakened_if = "95% CI covers 1 under either specification",
  strata_def  = "proliferation score tertiles; no PH assumption is placed on proliferation",
  tt_def      = "time-varying coefficient for the proliferation score, zp * log(time)",
  landmark    = "period split at 24 months")
w_res(data.frame(rule = names(RULES_C),
                 value = vapply(RULES_C, function(x) paste(x, collapse = ", "), character(1))),
      "CC00_ph_rules.csv")

dph <- CC[complete.cases(CC[, c("CYB5R3","prolif","grade_ord","age","sex","stage_bin","OS","OS.time")]), ]
dph$z5 <- as.numeric(scale(dph$CYB5R3))
dph$zp <- as.numeric(scale(dph$prolif))
dph$ptert <- cut(dph$zp, breaks = quantile(dph$zp, c(0, 1/3, 2/3, 1)),
                 include.lowest = TRUE, labels = c("low","mid","high"))
message("  n = ", nrow(dph), " ; deaths = ", sum(dph$OS == 1))
print(table(dph$ptert))

tidC <- function(fit, label, term = "z5") { s <- summary(fit)
  data.frame(model = label, term = term, n = fit$n, events = fit$nevent,
             HR = s$conf.int[term,"exp(coef)"], CI_low = s$conf.int[term,"lower .95"],
             CI_high = s$conf.int[term,"upper .95"], p = s$coefficients[term,"Pr(>|z|)"]) }

## C1. reference: the misspecified model from A5/B6
f_ref <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + zp, data = dph)

## C2. stratify on proliferation tertiles - no PH assumption on proliferation
f_str <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + strata(ptert), data = dph)

## C3. time-varying coefficient for proliferation
f_tt  <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + zp + tt(zp),
               data = dph, tt = function(x, t, ...) x * log(t + 1))

C1 <- dplyr::bind_rows(
  tidC(f_ref, "reference: proliferation as a constant term (violates PH)"),
  tidC(f_str, "PRE-SPECIFIED: stratified on proliferation tertiles"),
  tidC(f_tt,  "PRE-SPECIFIED: time-varying coefficient for proliferation"))
print(C1); w_res(C1, "CC01_ph_corrected_models.csv")

## C4. PH re-test on the stratified model
z2 <- cox.zph(f_str)
C2 <- data.frame(term = rownames(z2$table), chisq = z2$table[,"chisq"],
                 df = z2$table[,"df"], p = z2$table[,"p"], row.names = NULL)
print(C2); w_res(C2, "CC02_ph_retest_stratified.csv")

## C5. landmark split at 24 months
cut24 <- 24 * 30.44
early <- survival::survSplit(Surv(OS.time, OS) ~ ., data = dph, cut = cut24, episode = "period")
C3 <- dplyr::bind_rows(lapply(sort(unique(early$period)), function(k) {
  d <- early[early$period == k, ]
  if (nrow(d) < 60 || sum(d$OS == 1) < 15) return(NULL)
  cbind(period = if (k == 1) "0-24 months" else ">24 months",
        tidC(coxph(Surv(tstart, OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + zp, data = d),
             "landmark period"))
}))
if (!is.null(C3) && nrow(C3)) { print(C3); w_res(C3, "CC03_landmark_periods.csv") }

## C6. PH for the GSE14520 replication model
zg <- cox.zph(coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp, data = d14))
C4 <- data.frame(cohort = "GSE14520", term = rownames(zg$table), chisq = zg$table[,"chisq"],
                 df = zg$table[,"df"], p = zg$table[,"p"], row.names = NULL)
print(C4); w_res(C4, "CC04_ph_gse14520.csv")

## verdict
s_ok <- C1$HR[2] >= 1.25 && C1$CI_low[2] > 1
t_ok <- C1$HR[3] >= 1.25 && C1$CI_low[3] > 1
verdict <- if (s_ok && t_ok) "ROBUST to the PH violation" else
           if (C1$CI_low[2] <= 1 || C1$CI_low[3] <= 1) "WEAKENED once PH is handled properly" else
           "Indeterminate"
message("  VERDICT: ", verdict)
w_res(data.frame(HR_stratified = C1$HR[2], CI_low_stratified = C1$CI_low[2],
                 HR_timevarying = C1$HR[3], CI_low_timevarying = C1$CI_low[3],
                 verdict = verdict), "CC99_ph_verdict.csv")
save_session("21_ph")
message("\n[C PH addendum] done.")
