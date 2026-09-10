# Changelog

## v1.0.0 — 2026-09-10

First release. Accompanies the manuscript "Tumor proliferation masks the prognostic
association of CYB5R3 in hepatocellular carcinoma: a pre-registered two-cohort analysis".

Contents at this release:

- Pre-registered decision rules, fixed in code before the data were examined, and
  written to `results/C00_decision_rules.csv` on every run.
- Analysis scripts `analysis/18_subtype.R` through `analysis/22_figures.R`.
- The three Hoshida subclass gene sets as downloaded from MSigDB v2026.1.Hs, in
  `analysis/grp/`, frozen to `analysis/hoshida_templates.tsv` on first use.
- All result files, run logs and session information under `results/`.
- Manuscript figures in PNG (300 dpi), TIFF (LZW) and SVG under `figures/`.
- The pre-registration document, `PREREG_ProjectC_subtype.docx`.

No results were changed after the pre-registration was fixed. Analyses added after
the pre-registered tests were seen are labeled POST HOC in both the code and the
output files.
