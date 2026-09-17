# =====================================================================
#  19_addendum.R  —  Project C, POST HOC
#
#  ★★ 이 스크립트의 모든 분석은 POST HOC 이다. 18_subtype.R 의 결과를 본 뒤에
#     설계했다. 사전등록된 H1-H4 와 같은 지위로 보고해서는 안 되며, 원고에서
#     반드시 "post hoc, hypothesis-generating" 으로 표시한다.
#     Everything here was designed AFTER seeing the pre-registered results and
#     must be reported as post hoc, never alongside H1-H4 as if planned.
#
#  왜 필요한가: 아형 조성 가설이 기각되면서 역설이 그대로 남았다. 그런데 아직
#  가장 값싸고 결정적인 질문을 하지 않았다.
#     "이 코호트에서 조직등급이 애초에 예후인자인가?"
#  등급이 생존과 무관하다면 '등급이 오를수록 CYB5R3 감소'와 'CYB5R3 높으면
#  생존 불량'은 서로 모순이 아니다. 역설 자체가 해소된다.
#
#  Run: RUN_C_ADDENDUM.R 을 열고 Source (18 을 먼저 돌린 뒤 이어서 실행된다)
# =====================================================================
hdr("A - POST HOC addendum")
message("  NOTE: every result below is post hoc and must be labeled as such.")

PROLIF <- c("MKI67","TOP2A","CCNB1","PCNA","BUB1","CCNA2","AURKA","RRM2")
prolif_ok <- setdiff(intersect(PROLIF, rownames(EX)), unique(unlist(TEMPL)))
message("  proliferation genes used (template-independent): ",
        paste(prolif_ok, collapse = ", "))
ps <- colMeans(t(scale(t(as.matrix(EX[prolif_ok, , drop = FALSE])))), na.rm = TRUE)

strom_ok <- setdiff(intersect(STROMAL_IMMUNE, rownames(EX)), unique(unlist(TEMPL)))
## 표본이 관측치, 유전자가 변수가 되도록 전치한 뒤 PCA 를 수행한다.
## (이전 판은 유전자를 관측치로 두어 표본 수가 변수 개수가 되는 퇴화된 PCA 였다.)
pcm <- prcomp(t(as.matrix(EX[strom_ok, , drop = FALSE])), center = TRUE, scale. = TRUE)
cs_proxy <- pcm$x[, 1]; names(cs_proxy) <- colnames(EX)
message("  composition proxy: PC1 of ", length(strom_ok),
        " stromal/immune markers (samples as observations), ",
        round(100 * pcm$sdev[1]^2 / sum(pcm$sdev^2), 1), "% of variance")
## 주성분의 부호는 임의이다. 공변량으로만 쓰고, 상관계수의 부호는 해석하지 않는다.
message("  note: the sign of a principal component is arbitrary; it is used as a covariate only")
w_res(data.frame(
  score = c(rep("proliferation", length(prolif_ok)), rep("composition", length(strom_ok))),
  gene  = c(prolif_ok, strom_ok),
  note  = c(rep("mean of gene-wise z scores", length(prolif_ok)),
            rep(sprintf("PC1, %.1f%% of variance, samples as observations",
                        100 * pcm$sdev[1]^2 / sum(pcm$sdev^2)), length(strom_ok)))),
  "CA0_score_genes.csv")

for (nm in c("D", "CC")) {
  x <- get(nm); x$prolif <- ps[x$patient]; x$comp <- cs_proxy[x$patient]; assign(nm, x)
}

## ---- A1. 이 코호트에서 등급이 예후인자인가? (핵심 질문) ----------------
hdr("A1 - Is histologic grade prognostic here at all?")
g <- CC[!is.na(CC$grade_ord), ]
g$grade_f <- factor(g$grade_ord)
tid <- function(fit, label, term) {
  s <- summary(fit)
  data.frame(model = label, term = term, n = fit$n, events = fit$nevent,
             HR = s$conf.int[term,"exp(coef)"], CI_low = s$conf.int[term,"lower .95"],
             CI_high = s$conf.int[term,"upper .95"], p = s$coefficients[term,"Pr(>|z|)"])
}
A1 <- dplyr::bind_rows(
  tid(coxph(Surv(OS.time, OS) ~ grade_ord, data = g),
      "POST HOC: OS ~ grade (unadjusted)", "grade_ord"),
  tid(coxph(Surv(OS.time, OS) ~ grade_ord + age + sex + stage_bin, data = g),
      "POST HOC: OS ~ grade + age + sex + stage", "grade_ord"))
f_g0 <- coxph(Surv(OS.time, OS) ~ age + sex + stage_bin, data = g)
f_g1 <- coxph(Surv(OS.time, OS) ~ age + sex + stage_bin + grade_f, data = g)
A1$lrt_grade_factor <- c(NA, anova(f_g0, f_g1)$`Pr(>|Chi|)`[2])
A1$reading <- c(NA, if (A1$p[2] > 0.05)
  "Grade is NOT prognostic in this cohort. The discordance dissolves: grade and survival are decoupled here." else
  "Grade IS prognostic, so the discordance is real and remains unexplained.")
print(A1); w_res(A1, "CA1_posthoc_grade_prognostic.csv")

## ---- A2. JT 순열 p (정규근사 경고 제거) -------------------------------
hdr("A2 - Jonckheere-Terpstra by permutation")
dg <- D[!is.na(D$grade_ord) & !is.na(D$CYB5R3), ]
set.seed(20260910)
jt_perm <- clinfun::jonckheere.test(dg$CYB5R3, dg$grade_ord,
                                    alternative = "two.sided", nperm = 10000)
A2 <- data.frame(n = nrow(dg), JT_statistic = unname(jt_perm$statistic),
                 p_permutation = jt_perm$p.value, nperm = 10000,
                 note = "Replaces the normal-approximation p value; ties are present.")
print(A2); w_res(A2, "CA2_posthoc_JT_permutation.csv")

## ---- A3/A4. 증식·조성이 등급 추세를 설명하는가 -------------------------
hdr("A3/A4 - Does proliferation or tissue composition explain the grade trend?")
dg2 <- dg[!is.na(dg$prolif) & !is.na(dg$comp), ]
b <- function(f) { m <- lm(f, data = dg2)
  c(slope = unname(coef(m)["grade_ord"]),
    p = unname(summary(m)$coefficients["grade_ord","Pr(>|t|)"])) }
m0 <- b(rank(CYB5R3) ~ grade_ord)
mp <- b(rank(CYB5R3) ~ grade_ord + prolif)
mc <- b(rank(CYB5R3) ~ grade_ord + comp)
mb <- b(rank(CYB5R3) ~ grade_ord + prolif + comp)
A3 <- data.frame(
  adjustment = c("none", "+ proliferation score", "+ composition proxy", "+ both"),
  grade_slope = c(m0["slope"], mp["slope"], mc["slope"], mb["slope"]),
  p = c(m0["p"], mp["p"], mc["p"], mb["p"]))
A3$attenuation_vs_none <- (A3$grade_slope[1] - A3$grade_slope) / A3$grade_slope[1]
print(A3); w_res(A3, "CA3_posthoc_grade_slope_adjustments.csv")

A4 <- dplyr::bind_rows(
  cbind(pair = "CYB5R3 ~ proliferation", spearman_ci(dg2$CYB5R3, dg2$prolif)),
  cbind(pair = "CYB5R3 ~ composition proxy", spearman_ci(dg2$CYB5R3, dg2$comp)),
  cbind(pair = "grade ~ proliferation", spearman_ci(dg2$grade_ord, dg2$prolif)))
print(A4); w_res(A4, "CA4_posthoc_correlations.csv")

## ---- A5. 셋을 함께 넣은 예후 모형 --------------------------------------
hdr("A5 - Prognostic model with grade, proliferation and composition")
cj <- CC[!is.na(CC$grade_ord) & !is.na(CC$prolif) & !is.na(CC$comp), ]
cj$z5 <- as.numeric(scale(cj$CYB5R3))
A5 <- dplyr::bind_rows(
  tid(coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin, data = cj),
      "POST HOC: base", "z5"),
  tid(coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord, data = cj),
      "POST HOC: + grade", "z5"),
  tid(coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + prolif, data = cj),
      "POST HOC: + grade + proliferation", "z5"),
  tid(coxph(Surv(OS.time, OS) ~ z5 + age + sex + stage_bin + grade_ord + prolif + comp, data = cj),
      "POST HOC: + grade + proliferation + composition", "z5"))
print(A5); w_res(A5, "CA5_posthoc_cox_full.csv")

save_session("19_addendum")
message("\n[C addendum] done. All results above are POST HOC.")
