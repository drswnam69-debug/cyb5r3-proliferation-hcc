# =====================================================================
#  28_figure6.R  -  Figure 6: the genome-wide landscape
#
#  ★★ POST HOC, designed on 2026-09-20 with Stage F and Stage G.
#
#  (A) TCGA-LIHC. Every expressed gene, its correlation with the
#      proliferation score against the shift its hazard ratio takes when
#      proliferation enters the model. Density rather than points,
#      because there are sixteen thousand of them. The fitted line and
#      the decile medians are drawn on top, and CYB5R3 is marked.
#  (B) The same in GSE14520.
#  (C) What the shift costs. The number of genes reaching each threshold
#      before and after proliferation is adjusted for.
#
#  Run: after 26_genomewide.R and 27_sweep_refine.R.
# =====================================================================
if (!exists("GW") || !exists("save3"))
  stop("Run 26_genomewide.R and RUN_C_ALL.R first; GW and save3 must exist.")

suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })
hdr("Figure 6 - the genome-wide landscape")
message("  script version: 2026-09-20b")

OK <- "#2E5496"; BAD <- "#C00000"; GREY <- "#7F7F7F"
th6 <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 7.5),
        plot.subtitle = element_text(size = 6.3, colour = "grey30"),
        legend.key.size = unit(3, "mm"), legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6.2))

sweep_panel <- function(D, title, subtitle, mark_gene = NULL) {
  D <- D[is.finite(D$rho_with_proliferation) & is.finite(D$log_shift), ]
  ## trim the display to the central 99.5% so a handful of extreme genes do
  ## not set the scale; nothing is removed from any fit or any count
  ylim <- stats::quantile(D$log_shift, c(0.0025, 0.9975), na.rm = TRUE)
  dec <- stats::quantile(D$rho_with_proliferation, probs = seq(0, 1, 0.1), na.rm = TRUE)
  D$bin <- cut(D$rho_with_proliferation, breaks = dec, include.lowest = TRUE)
  M <- do.call(rbind, lapply(split(D, D$bin), function(d) data.frame(
    x = stats::median(d$rho_with_proliferation), y = stats::median(d$log_shift))))
  p <- ggplot(D, aes(rho_with_proliferation, log_shift)) +
    geom_hline(yintercept = 0, colour = GREY, linetype = 2, linewidth = 0.3) +
    geom_bin2d(bins = 90) +
    ## the shading is gene density on a log scale; it carries no number the
    ## reader needs, and two colour bars cost more width than they return,
    ## so the guide is suppressed and the scale is named in the legend text
    scale_fill_gradient(low = "#E8EEF6", high = "#1F3864",
                        trans = "log10", guide = "none") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = BAD, linewidth = 0.5) +
    geom_point(data = M, aes(x, y), inherit.aes = FALSE,
               shape = 21, fill = "white", colour = "black",
               size = 1.3, stroke = 0.35) +
    coord_cartesian(ylim = ylim) +
    labs(title = title, subtitle = subtitle,
         x = "Correlation with proliferation (Spearman)",
         y = "Log shift in the hazard ratio") + th6
  if (!is.null(mark_gene) && mark_gene %in% D$gene) {
    g <- D[D$gene == mark_gene, ][1, ]
    p <- p +
      geom_point(data = g, aes(rho_with_proliferation, log_shift), inherit.aes = FALSE,
                 shape = 21, fill = BAD, colour = "black", size = 1.8, stroke = 0.4) +
      geom_text(data = g, aes(rho_with_proliferation, log_shift, label = mark_gene),
                inherit.aes = FALSE, hjust = 1.25, vjust = -0.6,
                size = 2.1, colour = BAD, fontface = "bold")
  }
  p
}

slope_txt <- function(D) {
  m <- summary(stats::lm(log_shift ~ rho_with_proliferation, data = D))
  sprintf("slope %.3f, R² %.2f, %s genes",
          m$coefficients["rho_with_proliferation", "Estimate"],
          m$r.squared, format(nrow(D), big.mark = ","))
}

f6a <- sweep_panel(GW, "A  TCGA-LIHC", slope_txt(GW), "CYB5R3")

GGok <- exists("GG") && is.data.frame(GG) && nrow(GG) > 0
f6b <- if (GGok) sweep_panel(GG, "B  GSE14520", slope_txt(GG), "CYB5R3") else
  ggplot() + labs(title = "B  GSE14520", subtitle = "not available") + th6

## ---- (C) how many genes reach each threshold, before and after ------
cnt <- data.frame(
  threshold = factor(rep(c("p < 0.05\n(uncorrected)", "BH q < 0.10"), each = 2),
                     levels = c("p < 0.05\n(uncorrected)", "BH q < 0.10")),
  model = factor(rep(c("without proliferation", "with proliferation"), 2),
                 levels = c("without proliferation", "with proliferation")),
  n = c(sum(GW$p_without < 0.05, na.rm = TRUE),
        sum(GW$p_with    < 0.05, na.rm = TRUE),
        sum(GW$p_without_BH < 0.10, na.rm = TRUE),
        sum(GW$p_with_BH    < 0.10, na.rm = TRUE)))
print(cnt); w_res(cnt, "CF10_threshold_counts.csv")

f6c <- ggplot(cnt, aes(threshold, n, fill = model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = format(n, big.mark = ",")),
            position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.1) +
  scale_fill_manual(values = c("without proliferation" = GREY,
                               "with proliferation" = OK), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "C  Genes passing each threshold",
       subtitle = "TCGA-LIHC, the same 16,902 genes",
       x = NULL, y = "Genes") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  th6 + theme(legend.position = "bottom",
              legend.margin = margin(t = -4))

f6 <- ggarrange(f6a, f6b, f6c, ncol = 3, widths = c(1, 1, 0.9),
                common.legend = FALSE)
save3(f6, "Figure6_genomewide", w = 170, h = 80)   # 170 mm = BMC full column width
message("  Figure 6 written.")
if (exists("save_session")) save_session("28_figure6")
