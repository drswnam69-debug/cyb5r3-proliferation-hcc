# =====================================================================
#  32_figure_prereg.R  -  Figure 2 for the manuscript
#
#  ★★ POST HOC. This figure displays a relationship that was formed after
#  the registered test had been seen. The panel label says so, because a
#  reader should not have to consult the Methods to learn it.
#
#  (A) The registered prediction and what happened. The slope of each
#      cohort's landscape, ordered, with the line the registration
#      predicted every cohort would fall below.
#  (B) The structure in the failure. Each cohort's slope against the log
#      hazard ratio of its own proliferation score.
#
#  Inputs: CJ03_confirmatory_summary.csv (the ten registered cohorts) and
#  CI02_pancancer_summary.csv (the three exploratory cohorts). The two are
#  kept visually distinct throughout; they are never pooled into one
#  undifferentiated cloud.
#
#  Run: after 30_pancancer.R and 31_pancancer_confirm.R.
# =====================================================================
if (!exists("save3") || !exists("hdr"))
  stop("Run RUN_C_ALL.R first; save3 and hdr must exist.")

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })
hdr("Figure 2 - the registered test and the structure in its failure")
message("  script version: 2026-09-21c")
message("  NOTE: panel B is POST HOC.")

rd <- function(f) {
  p <- file.path(RES, f)
  if (!file.exists(p)) stop(f, " is missing from the results folder.")
  utils::read.csv(p, stringsAsFactors = FALSE)
}
CONF <- rd("CJ03_confirmatory_summary.csv")
PILOT <- rd("CI02_pancancer_summary.csv")

## the pilot file carries "TCGA-XXXX"; the confirmatory file carries "XXXX"
PILOT$cohort <- sub("^TCGA-", "", PILOT$cohort)
keep <- c("cohort", "slope", "r_squared", "proliferation_HR",
          "proliferation_lo", "proliferation_hi", "patients", "deaths")
D <- rbind(
  data.frame(CONF[, keep], set = "registered", stringsAsFactors = FALSE),
  data.frame(PILOT[PILOT$cohort %in% c("LIHC", "KIRC", "LUAD"), keep],
             set = "exploratory", stringsAsFactors = FALSE))
D$set <- factor(D$set, levels = c("registered", "exploratory"))
D$L <- log(D$proliferation_HR)
message("  cohorts plotted: ", nrow(D),
        " (registered ", sum(D$set == "registered"),
        ", exploratory ", sum(D$set == "exploratory"), ")")

OK <- "#2E5496"; BAD <- "#C00000"; GREY <- "#7F7F7F"
th2 <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 7),
        legend.margin = margin(t = -4))

## ---- (A) the registered prediction, and the outcome -----------------
A <- D[order(D$slope), ]
A$cohort <- factor(A$cohort, levels = A$cohort)
f2a <- ggplot(A, aes(cohort, slope, fill = set)) +
  geom_hline(yintercept = 0, colour = BAD, linewidth = 0.5) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = c(registered = OK, exploratory = GREY)) +
  coord_flip() +
  labs(title = "A  The registered prediction",
       subtitle = "all ten registered cohorts were predicted to be negative",
       x = NULL, y = "Slope of the landscape") + th2

## ---- (B) what the failure is made of (POST HOC) ---------------------
fit <- stats::lm(slope ~ L, data = D[D$set == "registered", ])
cf  <- stats::coef(fit)
r   <- stats::cor(D$slope[D$set == "registered"], D$L[D$set == "registered"])
sub <- sprintf("registered cohorts: r %.3f, slope %.3f, intercept %.3f",
               r, cf[2], cf[1])

f2b <- ggplot(D, aes(L, slope)) +
  geom_hline(yintercept = 0, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_abline(intercept = cf[1], slope = cf[2], colour = BAD, linewidth = 0.5) +
  geom_point(aes(shape = set, fill = set), size = 2.1, stroke = 0.4, colour = "black") +
  scale_shape_manual(values = c(registered = 21, exploratory = 24)) +
  scale_fill_manual(values = c(registered = OK, exploratory = "white")) +
  ## the five cohorts near the origin sit close together, so the labels need
  ## room to move; without this they land on top of their own points
  ggrepel::geom_text_repel(aes(label = cohort), size = 2.1, seed = 1,
                           min.segment.length = 0.15, box.padding = 0.45,
                           point.padding = 0.30, segment.size = 0.22,
                           segment.colour = "grey45", max.overlaps = Inf) +
  labs(title = "B  The structure in the failure (post hoc)",
       subtitle = sub,
       x = "Log hazard ratio of the proliferation score",
       y = "Slope of the landscape") + th2

f2 <- ggarrange(f2a, f2b, ncol = 2, widths = c(0.85, 1),
                common.legend = TRUE, legend = "bottom")
save3(f2, "Figure2_prereg", w = 170, h = 95)   # 170 mm = BMC full column width
w_res(D[order(D$L), c("cohort", "set", "patients", "deaths",
                      "proliferation_HR", "proliferation_lo",
                      "proliferation_hi", "L", "slope", "r_squared")],
      "CK01_figure2_data.csv")

## the numbers the manuscript quotes, recomputed here so the figure and the
## text cannot drift apart
fit_all <- stats::lm(slope ~ L, data = D)
S <- data.frame(
  quantity = c("cohorts, registered", "Pearson r, registered",
               "coefficient, registered", "intercept, registered",
               "cohorts, all", "Pearson r, all",
               "coefficient, all", "intercept, all",
               "sign concordance, all"),
  value = c(sum(D$set == "registered"), round(r, 4),
            round(cf[2], 4), round(cf[1], 4),
            nrow(D),
            round(stats::cor(D$slope, D$L), 4),
            round(stats::coef(fit_all)[2], 4),
            round(stats::coef(fit_all)[1], 4),
            paste0(sum((D$proliferation_HR > 1) == (D$slope < 0)), " of ", nrow(D))),
  stringsAsFactors = FALSE)
print(S); w_res(S, "CK02_figure2_summary.csv")
message("  Figure 2 written.")
if (exists("save_session")) save_session("32_figure_prereg")
