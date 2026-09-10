# =====================================================================
#  20_validation.R  —  Project C, Addendum B
#  PRE-REGISTERED external validation of the post hoc finding from 19.
#
#  Hypothesis (fixed before this script was run, and before GSE14520 was
#  touched for this purpose):
#    In TCGA-LIHC, adjusting for tumour proliferation RAISED the CYB5R3
#    hazard ratio (1.231 -> 1.371) and strengthened it (p 0.031 -> 0.0013),
#    because CYB5R3 is negatively correlated with proliferation while
#    proliferation is adverse. That is negative confounding, or suppression.
#    If real, the same direction must appear in an independent cohort.
#
#  This is confirmatory, so the reading is directional and fixed below.
#  Do not edit the RULES_B block after seeing any output.
# =====================================================================
hdr("B - PRE-REGISTERED external validation (GSE14520)")

RULES_B <- list(
  primary = "GSE14520: OS ~ CYB5R3(per SD) + age + sex + TNM, with and without the proliferation score",
  replication_declared_if = "HR increases when proliferation is added AND the adjusted 95% CI excludes 1",
  anchor_tolerance_HR = 0.05,   # published GSE14520 adjusted HR was 1.244
  min_genes_mapped    = 6,      # of the 8 proliferation genes
  prolif_genes = c("MKI67","TOP2A","CCNB1","PCNA","BUB1","CCNA2","AURKA","RRM2"),
  probe_rule   = "mean of all probes mapped to the gene symbol (same rule as the published pipeline)",
  multiplicity = "primary is one directional test; secondary analyses carry BH adjustment")
w_res(data.frame(rule = names(RULES_B),
                 value = vapply(RULES_B, function(x) paste(x, collapse = ", "), character(1))),
      "CB00_validation_rules.csv")

## ---- B1. GSE14520 assembly (mirrors the deposited pipeline) ----------
hdr("B1 - GSE14520 assembly")
suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })
gl14 <- getGEO("GSE14520", GSEMatrix = TRUE, destdir = CACHE)
## 19 에서 정한 기질/면역 유전자 목록도 함께 가져온다 (B8 탐색적 분석용).
WANT <- unique(c("CYB5R3", RULES_B$prolif_genes,
                 if (exists("strom_ok")) strom_ok else character(0)))
rows <- list(); plog <- list()
for (i in seq_along(gl14)) {
  es <- gl14[[i]]; ex <- exprs(es); fd <- fData(es)
  scol <- grep("Gene.?Symbol|GENE_SYMBOL|Symbol", names(fd), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(scol)) next
  syms <- toupper(as.character(fd[[scol]])); plat <- annotation(es)
  r <- data.frame(geo = colnames(ex), stringsAsFactors = FALSE)
  for (g in WANT) {
    idx <- which(syms %in% toupper(ALIAS[[g]] %||% g))
    if (!length(idx)) next
    r[[g]] <- as.numeric(colMeans(ex[idx, , drop = FALSE], na.rm = TRUE))
    plog[[length(plog)+1]] <- data.frame(platform = plat, gene = g, n_probes = length(idx),
      probes = paste(rownames(ex)[idx], collapse = ";"), rule = RULES_B$probe_rule)
  }
  if (ncol(r) > 1) rows[[length(rows)+1]] <- r
}
g14 <- dplyr::bind_rows(rows)
PL <- dplyr::bind_rows(plog); w_res(PL, "CB01_gse14520_probe_selection.csv")
mapped <- intersect(RULES_B$prolif_genes, names(g14))
message("  proliferation genes mapped: ", length(mapped), "/", length(RULES_B$prolif_genes),
        " (", paste(mapped, collapse = ", "), ")")
if (length(mapped) < RULES_B$min_genes_mapped) {
  w_res(data.frame(verdict = "NOT TESTABLE",
                   reason = "Too few proliferation genes map to probes on this platform."),
        "CB99_validation_verdict.csv")
  stop("Validation not testable, by the pre-registered rule.")
}

sup <- file.path(CACHE, "GSE14520_Extra_Supplement.txt.gz")
if (!file.exists(sup)) sup <- file.path(CACHE, "GSE14520_Extra_Supplement.txt")
cl14 <- data.table::fread(sup)
gg <- function(p) grep(p, names(cl14), ignore.case = TRUE, value = TRUE)[1]
key <- data.frame(geo = as.character(cl14[[gg("Affy.?GSM")]]),
                  time = suppressWarnings(as.numeric(cl14[[gg("Survival.*month")]])),
                  status = suppressWarnings(as.integer(cl14[[gg("Survival.*status")]])),
                  sexv = as.character(cl14[[gg("^Gender")]]),
                  agev = suppressWarnings(as.numeric(cl14[[gg("^Age")]])),
                  tnm  = as.character(cl14[[gg("TNM")]]), stringsAsFactors = FALSE)
m14 <- dplyr::inner_join(g14, key, by = "geo")
m14 <- m14[!is.na(m14$time) & m14$time > 0 & !is.na(m14$status), ]
tn <- toupper(trimws(m14$tnm))
m14$tnm_bin <- ifelse(grepl("^(III|IV|3|4)", tn), "III-IV",
               ifelse(grepl("^(I|II|1|2)",  tn), "I-II", NA))
message("  merged n = ", nrow(m14), " ; deaths = ", sum(m14$status == 1))

## proliferation score: same construction as TCGA (mean of z-scores)
zsc <- function(v) as.numeric(scale(v))
m14$prolif <- rowMeans(as.data.frame(lapply(m14[mapped], zsc)), na.rm = TRUE)

## ---- B2. Anchor: reproduce the published GSE14520 adjusted estimate ---
hdr("B2 - Anchor")
d14 <- m14[complete.cases(m14[, c("CYB5R3","agev","sexv","tnm_bin","prolif")]), ]
d14$sexv <- factor(d14$sexv); d14$tnm_bin <- factor(d14$tnm_bin, c("I-II","III-IV"))
d14$z5 <- zsc(d14$CYB5R3); d14$zp <- zsc(d14$prolif)
w_res(data.frame(
  step = c("Merged with survival", "Complete covariates (age, sex, TNM, proliferation)"),
  n = c(nrow(m14), nrow(d14)),
  deaths = c(sum(m14$status == 1), sum(d14$status == 1))),
  "CB01b_gse14520_flow.csv")
tidB <- function(fit, label, term) { s <- summary(fit)
  data.frame(model = label, term = term, n = fit$n, events = fit$nevent,
             HR = s$conf.int[term,"exp(coef)"], CI_low = s$conf.int[term,"lower .95"],
             CI_high = s$conf.int[term,"upper .95"], p = s$coefficients[term,"Pr(>|z|)"]) }
fU <- coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin, data = d14)
AN <- tidB(fU, "GSE14520 ADJUSTED (age, sex, TNM)", "z5")
AN$published_HR <- 1.244
AN$within_tolerance <- abs(AN$HR - 1.244) < RULES_B$anchor_tolerance_HR
print(AN); w_res(AN, "CB02_gse14520_anchor.csv")
if (!AN$within_tolerance)
  warning("GSE14520 anchor differs from the published estimate. Resolve before interpreting B3.")

## ---- B3. PRIMARY pre-registered test ---------------------------------
hdr("B3 - PRIMARY: does adjusting for proliferation raise the CYB5R3 HR?")
fP <- coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp, data = d14)
B3 <- dplyr::bind_rows(
  tidB(fU, "GSE14520: without proliferation", "z5"),
  tidB(fP, "GSE14520: with proliferation",    "z5"))
rose <- B3$HR[2] > B3$HR[1]
excl <- B3$CI_low[2] > 1
B3$replicated <- c(NA, rose && excl)
B3$reading <- c(NA, if (rose && excl)
  "REPLICATED. Proliferation suppresses the CYB5R3 hazard ratio in an independent cohort, as in TCGA-LIHC." else
  if (rose) "PARTIAL. The hazard ratio moves in the predicted direction but the interval still covers 1." else
  "NOT REPLICATED. The predicted direction did not appear.")
print(B3); w_res(B3, "CB03_PRIMARY_replication.csv")

## ---- B4. Is proliferation itself adverse, and is it negatively
##          correlated with CYB5R3, in both cohorts? --------------------
hdr("B4 - The two ingredients of suppression, in both cohorts")
pr14 <- tidB(coxph(Surv(time, status) ~ zp + agev + sexv + tnm_bin, data = d14),
             "GSE14520: proliferation score", "zp")
cj <- CC[complete.cases(CC[, c("CYB5R3","prolif","age","sex","stage_bin","OS","OS.time")]), ]
cj$zp <- zsc(cj$prolif)
prT <- tidB(coxph(Surv(OS.time, OS) ~ zp + age + sex + stage_bin, data = cj),
            "TCGA-LIHC: proliferation score", "zp")
B4a <- dplyr::bind_rows(prT, pr14)
B4b <- dplyr::bind_rows(
  cbind(cohort = "TCGA-LIHC", spearman_ci(cj$CYB5R3, cj$prolif)),
  cbind(cohort = "GSE14520",  spearman_ci(d14$CYB5R3, d14$prolif)))
B4b$expected <- "negative"
B4b$as_expected <- B4b$rho < 0
print(B4a); print(B4b)
w_res(B4a, "CB04_proliferation_prognostic.csv")
w_res(B4b, "CB05_cyb5r3_vs_proliferation.csv")

## ---- B5. TCGA: does CYB5R3 now add discrimination? -------------------
hdr("B5 - TCGA delta C-index with proliferation in the base model")
cg <- CC[complete.cases(CC[, c("CYB5R3","prolif","grade_ord","age","sex","stage_bin","OS","OS.time")]), ]
cg$z5 <- zsc(cg$CYB5R3); cg$zp <- zsc(cg$prolif)
B5 <- c_index_delta(cg, "OS.time", "OS",
                    "age + sex + stage_bin + grade_ord + zp",
                    "age + sex + stage_bin + grade_ord + zp + z5", B = 1000)
B5$comparison <- "adding CYB5R3 to a model that already contains proliferation and grade"
B5$published_dC_without_proliferation <- 0.0047
print(B5); w_res(B5, "CB06_tcga_delta_cindex.csv")

## ---- B6. Proportional hazards for the proliferation-adjusted model ---
hdr("B6 - PH assumption")
fT <- coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + zp, data = cg)
zph <- cox.zph(fT)
B6 <- data.frame(term = rownames(zph$table), chisq = zph$table[,"chisq"],
                 df = zph$table[,"df"], p = zph$table[,"p"], row.names = NULL)
print(B6); w_res(B6, "CB07_ph_test.csv")

## ---- B7. Sensitivity: does the result hang on one gene? --------------
hdr("B7 - Sensitivity of the proliferation score")
sens <- list()
base_rhs <- "age + sex + stage_bin + grade_ord"
one <- function(v, label) {
  d <- cg; d$zz <- zsc(v[d$patient])
  f <- coxph(as.formula(paste("Surv(OS.time, OS) ~ z5 +", base_rhs, "+ zz")), data = d)
  cbind(definition = label, tidB(f, "TCGA: CYB5R3 with this proliferation definition", "z5"))
}
allp <- intersect(RULES_B$prolif_genes, rownames(EX))
zmat <- t(scale(t(as.matrix(EX[allp, , drop = FALSE]))))
sens[["all8"]] <- one(colMeans(zmat, na.rm = TRUE), paste0("all ", length(allp), " genes"))
for (g in allp) {
  keep <- setdiff(allp, g)
  sens[[paste0("drop_", g)]] <- one(colMeans(zmat[keep, , drop = FALSE], na.rm = TRUE),
                                    paste0("leave out ", g))
}
for (g in c("MKI67","TOP2A")) if (g %in% allp)
  sens[[paste0("only_", g)]] <- one(zmat[g, ], paste0(g, " alone"))
B7 <- dplyr::bind_rows(sens)
print(B7); w_res(B7, "CB08_sensitivity_proliferation_definition.csv")

## ---- B8. EXPLORATORY: tissue composition in GSE14520 -----------------
##  ★ 사전지정 아님. TCGA-LIHC 의 조성 결과를 본 뒤에 추가한 탐색적 분석이다.
##    원고에서 반드시 exploratory 로 표시하고, 사전등록된 증식 결과와 같은
##    지위로 제시하지 않는다. 판정 문구는 실행 전에 아래에 고정했다.
hdr("B8 - EXPLORATORY: does tissue composition behave the same way here?")
RULES_B8 <- list(
  status = "EXPLORATORY, not pre-specified; added after the TCGA-LIHC composition result was seen",
  consistent_if = "adding composition on top of proliferation raises the CYB5R3 HR AND the 95% CI excludes 1",
  min_genes = "at least 6 of the stromal and immune markers must map to probes",
  sign_note = "the sign of a principal component is arbitrary and need not agree between cohorts; it is used as a covariate only")
w_res(data.frame(rule = names(RULES_B8), value = unlist(RULES_B8, use.names = FALSE)),
      "CB09_exploratory_composition_rules.csv")

strom_col <- if (exists("strom_ok")) intersect(strom_ok, names(m14)) else character(0)
## GSE14520 은 두 플랫폼(GPL3921, GPL571)을 합친 자료라, 한쪽에만 있는 유전자는
## 반대쪽 플랫폼 검체에서 결측이 된다. 그런 유전자를 남겨 두면 PCA 가 실패한다.
## 검체를 잃지 않기 위해, 두 플랫폼 모두에 있는 유전자만 사용한다.
strom14 <- if (length(strom_col))
  strom_col[vapply(strom_col, function(g) all(is.finite(m14[[g]])), logical(1))] else character(0)
drop14 <- setdiff(strom_col, strom14)
message("  stromal/immune markers: ", length(strom14), " usable",
        if (exists("strom_ok")) paste0(" of ", length(strom_ok), " used in TCGA-LIHC") else "")
if (length(drop14))
  message("  dropped, present on only one of the two platforms: ", paste(drop14, collapse = ", "))
w_res(data.frame(gene = c(strom14, drop14),
                 used = c(rep(TRUE, length(strom14)), rep(FALSE, length(drop14))),
                 reason = c(rep("present on both platforms", length(strom14)),
                            rep("present on only one platform; would force missing values", length(drop14)))),
      "CB09b_exploratory_composition_genes.csv")

if (length(strom14) >= 6) {
  pc14 <- prcomp(as.matrix(m14[, strom14, drop = FALSE]), center = TRUE, scale. = TRUE)
  m14$comp <- pc14$x[, 1]
  vex <- 100 * pc14$sdev[1]^2 / sum(pc14$sdev^2)
  message("  PC1 explains ", round(vex, 1), "% of variance (TCGA-LIHC value was 47.6%)")

  d8 <- m14[complete.cases(m14[, c("CYB5R3","agev","sexv","tnm_bin","prolif","comp")]), ]
  d8$sexv <- factor(d8$sexv); d8$tnm_bin <- factor(d8$tnm_bin, c("I-II","III-IV"))
  d8$z5 <- zsc(d8$CYB5R3); d8$zp <- zsc(d8$prolif); d8$zc <- zsc(d8$comp)

  e1 <- coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin, data = d8)
  e2 <- coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp, data = d8)
  e3 <- coxph(Surv(time, status) ~ z5 + agev + sexv + tnm_bin + zp + zc, data = d8)
  B8 <- dplyr::bind_rows(
    tidB(e1, "GSE14520 EXPLORATORY: age, sex, TNM", "z5"),
    tidB(e2, "GSE14520 EXPLORATORY: + proliferation", "z5"),
    tidB(e3, "GSE14520 EXPLORATORY: + proliferation + composition", "z5"))
  rose <- B8$HR[3] > B8$HR[2]; excl <- B8$CI_low[3] > 1
  B8$pc1_variance_percent <- c(NA, NA, round(vex, 1))
  B8$composition_genes <- c(NA, NA, paste(strom14, collapse = ";"))
  B8$reading <- c(NA, NA,
    if (rose && excl) "CONSISTENT with TCGA-LIHC: composition masks part of the association here too" else
    if (rose) "PARTIAL: the direction matches but the interval covers 1" else
    "NOT CONSISTENT: the direction seen in TCGA-LIHC does not appear here")
  print(B8); w_res(B8, "CB10_exploratory_composition.csv")

  ## 조성 자체의 예후, 그리고 CYB5R3 와의 관계
  B8b <- dplyr::bind_rows(
    cbind(quantity = "composition score, own hazard ratio",
          tidB(coxph(Surv(time, status) ~ zc + agev + sexv + tnm_bin, data = d8),
               "GSE14520 EXPLORATORY", "zc")[, c("n","events","HR","CI_low","CI_high","p")]))
  cc8 <- cbind(quantity = "CYB5R3 vs composition (Spearman)",
               spearman_ci(d8$CYB5R3, d8$comp))
  print(B8b); print(cc8)
  w_res(B8b, "CB11_exploratory_composition_prognostic.csv")
  w_res(cc8, "CB12_exploratory_cyb5r3_vs_composition.csv")
} else {
  w_res(data.frame(verdict = "NOT TESTABLE",
                   reason = "Too few stromal and immune markers map to probes on this platform."),
        "CB10_exploratory_composition.csv")
  message("  not testable on this platform")
}

## ---- verdict ---------------------------------------------------------
hdr("B - verdict")
w_res(data.frame(
  gse14520_anchor_ok = AN$within_tolerance,
  HR_without_proliferation = B3$HR[1], HR_with_proliferation = B3$HR[2],
  replicated = B3$replicated[2], reading = B3$reading[2],
  proliferation_adverse_TCGA = prT$HR > 1, proliferation_adverse_GSE14520 = pr14$HR > 1,
  negative_correlation_both = all(B4b$as_expected)),
  "CB99_validation_verdict.csv")
save_session("20_validation")
message("\n[C validation] done.")
