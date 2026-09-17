# =====================================================================
#  22_figures.R  —  Project C, manuscript figures
#  Figures 1 to 4 as PNG (300 dpi), TIFF (LZW) and SVG.
#  Run after 18 to 21 so every value comes from the objects that
#  produced the numbers in the text.
# =====================================================================
hdr("Figures")
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })

OK <- "#2E5496"; BAD <- "#C00000"; GREY <- "#7F7F7F"
th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 7),
        axis.title    = element_text(size = 8),
        axis.text     = element_text(size = 7.5),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 4))

## TIFF 의 실제 압축 방식을 파일에서 읽어 확인한다 (라벨을 믿지 않는다).
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
  ggsave(paste0(base, ".png"), p, width = win, height = hin, units = "in", dpi = 300)

  ## LZW 압축 TIFF. macOS 에는 cairo 가 없는 경우가 많으므로 ragg 를 쓴다.
  ## (grDevices::tiff(type="cairo") 는 cairo 가 없으면 경고만 내고 quartz 로
  ##  조용히 넘어가 무압축 파일을 만든다. 그래서 결과를 파일에서 검증한다.)
  tif <- paste0(base, ".tiff"); ok <- FALSE
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(tif, width = win, height = hin, units = "in", res = 300,
                   compression = "lzw", background = "white")
    print(p); grDevices::dev.off()
    ok <- isTRUE(tiff_compression(tif) == 5L)
  }
  if (!ok) {
    ggsave(tif, p, width = win, height = hin, units = "in", dpi = 300)
    message("  [!] ", basename(tif), " is UNCOMPRESSED. Run install.packages(\"ragg\") for LZW.")
  }
  svg <- tryCatch({ suppressMessages(
    ggsave(paste0(base, ".svg"), p, width = win, height = hin, units = "in")); TRUE },
    error = function(e) FALSE)

  ## Vector EPS for journals that require vector line art (Scientific Reports).
  ## Remove any earlier file first. cairo_ps() only WARNS when the cairo DLL
  ## cannot be loaded, so without this the existence-and-size check below
  ## succeeds on the previous run's file and the console reports an EPS that
  ## this run never wrote.
  eps <- paste0(base, ".eps"); epsok <- FALSE
  if (file.exists(eps)) unlink(eps)
  if (isTRUE(capabilities("cairo"))) {
    epsok <- tryCatch({
      grDevices::cairo_ps(eps, width = win, height = hin, fallback_resolution = 800,
                          onefile = FALSE, bg = "white")
      print(p); grDevices::dev.off(); file.exists(eps) && file.size(eps) > 5000
    }, error = function(e) FALSE)
  }
  epsfb <- FALSE
  if (!epsok) {
    ## No cairo on this machine. postscript() writes a valid EPS but uses
    ## Helvetica AFM metrics, which are wider than the metrics used for the
    ## PNG/SVG, so panel titles can be clipped. Such a file is a placeholder:
    ## the submission EPS is converted from the SVG (see figures/README_EPS.txt).
    epsok <- tryCatch({
      grDevices::postscript(eps, width = win, height = hin, onefile = FALSE,
                            horizontal = FALSE, paper = "special", bg = "white")
      print(p); grDevices::dev.off(); file.exists(eps) && file.size(eps) > 5000
    }, error = function(e) FALSE)
    epsfb <- epsok
  }
  epstag <- if (!epsok) " [!] EPS FAILED" else if (epsfb)
    ", eps [!] postscript fallback, convert from SVG" else ", eps"
  message("  -> ", basename(base), ": png, tiff (",
          if (ok) "LZW verified" else "uncompressed", ")", if (svg) ", svg" else "",
          epstag)
}

## ---------------------------------------------------------------- Fig 1
B <- data.frame(
  x = c(1, 1, 4.3, 4.3, 7.6, 7.6),
  y = c(3, 1.85, 3, 1.85, 3, 1.85),
  w = 2.7, h = c(1.0, 0.62, 1.0, 0.62, 1.0, 0.62),
  fill = c("#DEEAF6", "#F2F2F2", "#FFF2CC", "#F2F2F2", "#DEEAF6", "#F2F2F2"),
  lab = c("PRE-REGISTERED\nDoes molecular subclass explain\nthe discordance? (H1 to H4)",
          "All four fail.\nSubclass explains nothing.",
          "POST HOC\nIs grade prognostic at all?\nWhat is the grade trend made of?",
          "Grade is not prognostic.\nThe trend is proliferation.",
          "RULE FIXED IN ADVANCE\nDoes adjusting for proliferation\nraise the hazard ratio?",
          "Replicated in GSE14520.\nRobust to the PH violation."),
  stringsAsFactors = FALSE)
f1 <- ggplot() +
  geom_rect(data = B, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey45", linewidth = 0.3) +
  scale_fill_identity() +
  geom_text(data = B, aes(x, y, label = lab), size = 2.3, lineheight = 1.15) +
  geom_segment(data = data.frame(x = c(2.45, 5.75)), aes(x = x, xend = x + 0.5, y = 3, yend = 3),
               arrow = arrow(length = unit(1.8, "mm")), linewidth = 0.4, colour = "grey35") +
  geom_segment(data = B[c(1, 3, 5), ], aes(x = x, xend = x, y = y - 0.5, yend = y - 0.75),
               arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.3, colour = "grey55") +
  annotate("text", x = c(1, 4.3, 7.6), y = 1.35, size = 2.1, colour = "grey35",
           label = c("rules fixed before any data were seen",
                     "designed after the results above",
                     "rule fixed before this cohort was used")) +
  coord_cartesian(xlim = c(-0.5, 9.1), ylim = c(1.15, 3.65)) +
  theme_void(base_size = 8)
save3(f1, "Figure1_design", 175, 60)

## ---------------------------------------------------------------- Fig 2
cnt <- as.data.frame(table(SUB$subclass), stringsAsFactors = FALSE)
names(cnt) <- c("subclass", "n")
f2a <- ggplot(cnt, aes(subclass, n)) +
  geom_col(fill = OK, width = 0.62) +
  geom_text(aes(label = n), vjust = -0.45, size = 2.4) +
  expand_limits(y = max(cnt$n) * 1.18) +
  labs(title = "A  Subclass assignment", x = NULL, y = "Tumors") + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

cd <- CTRL[!is.na(CTRL$class), ]
cd$marker <- factor(cd$marker, levels = rev(cd$marker))
f2b <- ggplot(cd, aes(delta, marker, colour = class)) +
  geom_vline(xintercept = 0, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c(S2 = BAD, S3 = OK)) +
  expand_limits(x = c(0, max(cd$delta) * 1.1)) +
  labs(title = "B  Validity controls", x = "Median difference, in minus out (log2)", y = NULL) +
  th + theme(legend.position = "bottom", legend.title = element_blank(),
             legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7),
             legend.margin = margin(0, 0, 0, 0))

dd2 <- D[D$subclass %in% c("S1", "S2", "S3") & !is.na(D$CYB5R3), ]
f2c <- ggplot(dd2, aes(subclass, CYB5R3)) +
  geom_boxplot(outlier.size = 0.4, width = 0.52, fill = "#EDEDED", linewidth = 0.3) +
  labs(title = "C  CYB5R3 by subclass",
       subtitle = sprintf("Kruskal-Wallis p = %.3f", H1$p),
       x = NULL, y = "CYB5R3 (log2)") + th
f2 <- ggarrange(f2a, f2b, f2c, ncol = 3, widths = c(0.85, 1.5, 0.95))
save3(f2, "Figure2_subclass", 175, 70)

## ---------------------------------------------------------------- Fig 3
dg3 <- D[!is.na(D$grade_ord) & !is.na(D$CYB5R3), ]
dg3$g <- factor(paste0("G", dg3$grade_ord))
f3a <- ggplot(dg3, aes(g, CYB5R3)) +
  geom_boxplot(outlier.size = 0.4, width = 0.55, fill = "#EDEDED", linewidth = 0.3) +
  labs(title = "A  CYB5R3 by grade",
       subtitle = sprintf("permutation p = %.3f", A2$p_permutation),
       x = NULL, y = "CYB5R3 (log2)") + th

A3b <- A3
A3b$adjustment <- factor(A3b$adjustment, levels = rev(A3b$adjustment))
A3b$sig <- A3b$p < 0.05
xr <- range(c(A3b$grade_slope, 0))
f3b <- ggplot(A3b, aes(grade_slope, adjustment, colour = sig)) +
  geom_vline(xintercept = 0, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_point(size = 2.1) +
  geom_text(aes(x = 1.5, label = sprintf("p = %.3f", p)), hjust = 0, size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = BAD, `FALSE` = GREY)) +
  coord_cartesian(xlim = c(xr[1] * 1.08, 12)) +
  labs(title = "B  Grade slope",
       subtitle = "red, p < 0.05", x = "Slope on ranked CYB5R3", y = NULL) + th

A1b <- A1; A1b$lab <- factor(c("Unadjusted", "Age, sex, stage"),
                             levels = c("Age, sex, stage", "Unadjusted"))
f3c <- ggplot(A1b, aes(HR, lab)) +
  geom_vline(xintercept = 1, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.1, colour = GREY, linewidth = 0.4) +
  geom_point(size = 2.1, colour = GREY) +
  geom_text(aes(x = 1.58, label = sprintf("%.2f (%.2f-%.2f)\np = %.2f", HR, CI_low, CI_high, p)),
            hjust = 0, size = 2.1, lineheight = 1.1) +
  coord_cartesian(xlim = c(0.82, 2.25)) +
  labs(title = "C  Grade is not prognostic", x = "Hazard ratio for overall survival", y = NULL) + th
f3 <- ggarrange(f3a, f3b, f3c, ncol = 3, widths = c(0.85, 1.25, 1.2))
save3(f3, "Figure3_grade_proliferation", 175, 70)

## ---------------------------------------------------------------- Fig 4
gr <- function(d, i, coh, lab) data.frame(cohort = coh, lab = lab,
  HR = d$HR[i], CI_low = d$CI_low[i], CI_high = d$CI_high[i], p = d$p[i])
FA <- rbind(
  gr(A5, 1, "TCGA-LIHC", "Age, sex, stage"),
  gr(A5, 2, "TCGA-LIHC", "+ grade"),
  gr(A5, 3, "TCGA-LIHC", "+ grade + proliferation"),
  gr(C1, 2, "TCGA-LIHC", "Stratified on proliferation tertiles"),
  gr(C1, 3, "TCGA-LIHC", "Time-varying proliferation coefficient"),
  gr(B3, 1, "GSE14520",  "Age, sex, TNM"),
  gr(B3, 2, "GSE14520",  "+ proliferation"))
## 조성 보정은 탐색적이고 재현되지 않았으므로 본문 그림에서 제외하고
## Supporting Information 에만 싣는다 (원고 Table 3 과 일치시킨다).
FA$cohort <- factor(FA$cohort, levels = c("TCGA-LIHC", "GSE14520"))
FA$lab <- factor(FA$lab, levels = rev(FA$lab))
FA$adj <- grepl("proliferation", FA$lab, ignore.case = TRUE)
xtxt <- max(FA$CI_high) * 1.04
f4a <- ggplot(FA, aes(HR, lab, colour = adj)) +
  geom_vline(xintercept = 1, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.1, linewidth = 0.4) +
  geom_point(size = 1.9) +
  geom_text(aes(x = xtxt, label = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)),
            hjust = 0, size = 2.1) +
  scale_colour_manual(values = c(`FALSE` = GREY, `TRUE` = BAD)) +
  coord_cartesian(xlim = c(0.92, xtxt * 1.40)) +
  facet_grid(cohort ~ ., scales = "free_y", space = "free_y") +
  labs(title = "A  Hazard ratio by adjustment set",
       subtitle = "red, models containing proliferation",
       x = "Hazard ratio per standard deviation", y = NULL) +
  th + theme(strip.background = element_rect(fill = "#EDEDED"),
             strip.text = element_text(size = 7.5))

B7b <- B7; B7b$definition <- factor(B7b$definition, levels = rev(B7b$definition))
f4b <- ggplot(B7b, aes(HR, definition)) +
  geom_vline(xintercept = 1, colour = GREY, linetype = 2, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.1, colour = BAD, linewidth = 0.4) +
  geom_point(size = 1.8, colour = BAD) +
  coord_cartesian(xlim = c(0.96, 1.78)) +
  labs(title = "B  Sensitivity of the score",
       subtitle = "TCGA-LIHC, eleven definitions",
       x = "Hazard ratio per standard deviation", y = NULL) + th
f4 <- ggarrange(f4a, f4b, ncol = 2, widths = c(1.5, 1))
save3(f4, "Figure4_suppression", 175, 89)

message("\n[figures] done -> ", FIG)
