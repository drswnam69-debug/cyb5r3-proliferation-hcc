Vector figure files (EPS and PDF)
=================================

Each figure is archived in four raster and vector formats at a fixed width of
175 mm: PNG (300 dpi), TIFF (300 dpi, LZW), SVG, EPS and PDF.

PNG, TIFF and SVG are written directly by analysis/22_figures.R and
analysis/24_figure5.R through save3().

EPS and PDF are derived from the archived SVG, not written by R. The Mac used
for this analysis has no cairo device (X11/XQuartz absent), so
grDevices::cairo_ps() is unavailable and save3() falls back to
grDevices::postscript(). That fallback produces a valid EPS but lays text out
with Helvetica AFM metrics, which are wider than the metrics used for the PNG
and SVG; panel titles are clipped as a result. Files written by the fallback
are placeholders only and are flagged in the console output of save3().

The submission EPS and PDF were converted from the archived SVG, which carries
explicit textLength and lengthAdjust attributes for every text element and
therefore reproduces the original geometry exactly:

  rsvg-convert -f pdf -o FigureN.pdf FigureN.svg
  gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=eps2write \
     -sOutputFile=FigureN.eps FigureN.pdf

librsvg 2.x and Ghostscript 10.x. Fonts are embedded in both outputs. The
bounding box of Figure1_design.eps was set to the full 496.06 x 170.08 pt page,
because eps2write crops to the inked area and that panel has white margins.

Each converted file was rendered back to PNG and compared against the PNG
written by R to confirm that no text, symbol or line was lost.
