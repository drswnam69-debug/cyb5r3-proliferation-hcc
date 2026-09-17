# =====================================================================
#  24_figure5.R  -  Figure 5 with a numbers-at-risk table
#
#  Optional. Reviewers of a clinical journal expect a risk table under a
#  Kaplan-Meier plot; this redraws Figure 5 with one. It needs only CC
#  (TCGA-LIHC) and d14 (GSE14520), which are already in the workspace after
#  RUN_C_CLINICAL.R or RUN_C_FIGURES.R, so it runs in seconds.
#
#  Run: open this file and Source it, after either of those wrappers.
# =====================================================================
if (!exists("CC") || !exists("d14"))
  stop("Run RUN_C_CLINICAL.R (or RUN_C_FIGURES.R) first; CC and d14 must be in the workspace.")

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(survival) })
hdr("Figure 5 with numbers at risk")

if (!exists("th")) th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 7))
COLS <- c(low = "#2E5496", middle = "#7F7F7F", high = "#C00000")

tert <- function(x) {
  q <- quantile(x, c(1/3, 2/3), na.rm = TRUE)
  factor(ifelse(x <= q[1], "low", ifelse(x <= q[2], "middle", "high")),
         c("low", "middle", "high"))
}

pfmt <- function(p) if (p < 0.1) sprintf("log-rank p = %.3f", p) else sprintf("log-rank p = %.2f", p)

panel <- function(d, tv, ev, title, xlab, breaks) {
  d$tert <- tert(d$CYB5R3)
  ## the axis covers the whole of follow-up, so no event is hidden by the display
  xlim <- c(0, max(d[[tv]], na.rm = TRUE))
  f <- survfit(as.formula(sprintf("Surv(%s, %s) ~ tert", tv, ev)), data = d)
  lr <- survdiff(as.formula(sprintf("Surv(%s, %s) ~ tert", tv, ev)), data = d)
  p  <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  kd <- data.frame(time = f$time, surv = f$surv,
                   tert = factor(rep(sub("^tert=", "", names(f$strata)), f$strata),
                                 c("low", "middle", "high")))
  ## start every curve at 1 so the step plot does not begin mid-air
  kd <- rbind(data.frame(time = 0, surv = 1,
                         tert = factor(c("low","middle","high"), c("low","middle","high"))), kd)
  g <- ggplot(kd, aes(time, surv, colour = tert)) +
    geom_step(linewidth = 0.5) +
    scale_colour_manual(values = COLS) +
    scale_x_continuous(limits = xlim, breaks = breaks) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Overall survival", title = title, subtitle = pfmt(p)) +
    th + theme(legend.position = "none",
               axis.text.x = element_blank(), axis.ticks.x = element_blank(),
               plot.margin = margin(4, 6, 0, 4))

  at <- summary(f, times = breaks, extend = TRUE)
  rt <- data.frame(time = at$time,
                   n = at$n.risk,
                   tert = factor(sub("^tert=", "", as.character(at$strata)),
                                 c("low", "middle", "high")))
  r <- ggplot(rt, aes(time, tert, label = n, colour = tert)) +
    geom_text(size = 2.1) +
    scale_colour_manual(values = COLS) +
    scale_x_continuous(limits = xlim, breaks = breaks) +
    scale_y_discrete(limits = rev(c("low", "middle", "high"))) +
    labs(x = xlab, y = NULL, title = "Number at risk") +
    th + theme(legend.position = "none",
               panel.grid = element_blank(), panel.border = element_blank(),
               axis.ticks.y = element_blank(),
               plot.title = element_text(face = "plain", size = 7),
               plot.margin = margin(0, 6, 4, 4))
  ggarrange(g, r, ncol = 1, heights = c(3.1, 1))
}

p1 <- panel(CC,  "OS.time", "OS",     "A  TCGA-LIHC", "Days since diagnosis",   seq(0, 3000, 1000))
p2 <- panel(d14, "time",    "status", "B  GSE14520",  "Months since diagnosis", seq(0, 60, 20))

leg <- get_legend(
  ggplot(data.frame(x = 1:3, tert = factor(c("low","middle","high"),
                                           c("low","middle","high"))),
         aes(x, x, colour = tert)) + geom_line() +
    scale_colour_manual(values = COLS) +
    labs(colour = "CYB5R3 tertile") + th + theme(legend.position = "bottom"))

f5 <- ggarrange(ggarrange(p1, p2, ncol = 2), leg, ncol = 1, heights = c(10, 1))
save3(f5, "Figure5_km", w = 175, h = 97)
message("  Figure 5 redrawn with numbers at risk.")
if (exists("save_session")) save_session("24_figure5")
