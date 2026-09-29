# Changelog

## v2.1.0 - 2026-09-29

A figure-format release. No analysis, no decision rule and no reported number
changed. Every result file in this release, including the three the two figure
scripts rewrite on each run, is byte-identical to v2.0.0.

Changed:

- `analysis/28_figure6.R` and `analysis/32_figure_prereg.R` now write their
  figures at 170 mm rather than 175 mm, which is the full text width of the
  journal the manuscript is now submitted to. The figures in `figures/` were
  regenerated at that width and the exported EPS files were rebuilt from the
  SVG through `rsvg-convert` and Ghostscript, because `cairo_ps()` is not
  available on the analysis machine and the `postscript()` fallback clips text.

Unchanged and reconfirmed on this run: `CF10` (2308, 1333, 338, 1) and `CK02`
(registered Pearson r -0.9876, coefficient -0.9559, intercept -0.0077; all
thirteen cohorts r -0.9916; sign concordance 13 of 13).

## v2.0.0 — 2026-09-23

The study changed scope. What began as a single-gene analysis in hepatocellular
carcinoma is now a methodological study of what omitting one covariate does to
prognostic gene expression estimates, measured across every expressed gene and
across thirteen cancer cohorts. The version is a major one because the earlier
releases cannot reproduce anything in the current manuscript.

Added, with each stage labeled in the code and in its output files:

- `analysis/26_genomewide.R`, Stage F, exploratory. The genome-wide sweep in
  TCGA-LIHC, its replication in GSE14520, and a repeat with
  HALLMARK_G2M_CHECKPOINT in place of the eight-gene score. Results `CF00` to
  `CF10`.
- `analysis/27_sweep_refine.R`, Stage G, post hoc. Non-collapsibility read off
  the genes that carry no confounding, reversals counted among genes that carry
  an effect, and the consequences under Benjamini-Hochberg control. Results
  `CG00` to `CG06`.
- `analysis/28_figure6.R`. The landscape figure.
- `analysis/29_litsurvey.R`, Stage H. A survey of 50 recent prognostic gene
  expression papers drawn under a fixed seed from a frame of 424 records. The
  design was fixed before any paper was read and is written to `CH00` on every
  run. Results `CH00` to `CH05`.
- `analysis/30_pancancer.R`, Stage I, exploratory pilot. The sweep in TCGA-KIRC
  and TCGA-LUAD, with TCGA-LIHC refitted under a common covariate base so the
  three are comparable. Results `CI00` to `CI02`.
- `analysis/31_pancancer_confirm.R`, Stage J, PREREGISTERED. The confirmatory
  pan-cancer sweep. The plan, `PREREG_PanCancer_20260920.docx`, was registered
  as doi:10.17605/OSF.IO/X5DCF on 2026-09-20; the script ran on 2026-09-21 and
  its run log records that. Ten cancer types met the registered inclusion rule.
  **Both registered hypotheses were refuted**, and `CJ05_decision.csv` records
  that verdict. Results `CJ00` to `CJ05`.
- `analysis/32_figure_prereg.R`, Stage K. The Stage J figure, which recomputes
  the numbers the manuscript quotes into `CK02` so that figure and text cannot
  drift apart.
- `analysis/33_score_validity.R`, Stage L, post hoc. Whether the eight-gene
  proliferation score behaves as an instrument outside the liver. Each cohort's
  analysis set is rebuilt by the same rule and checked against that sweep's own
  patient count. Results `CL00` to `CL03`.
- `analysis/34_attenuation.R`, Stage M, post hoc. Whether measurement error can
  account for the relationship the refutation revealed. Results `CM01` to `CM03`.
- `PREREG_PanCancer_20260920.docx`, the Stage J preregistration.
- `PREREG2_ConfoundingLaw_20260920.docx`, doi:10.17605/OSF.IO/N8DW6, registering
  the test of the post hoc relationship in cohorts that have not been examined.
  No analysis covered by it has been run, and none of it is in this release.
- Figures now carry EPS and PDF alongside PNG, TIFF and SVG, with
  `figures/README_EPS.txt` recording how the vector files were produced and why.

Changed:

- `README.md` rewritten for the current scope, with a table mapping every result
  prefix to its stage and to its status as preregistered, post hoc or
  exploratory.
- `analysis/23_clinical.R` made idempotent, so that re-sourcing it no longer
  produces duplicated join columns.
- `RUN_C_ALL.R` runs `24_figure5.R` last, so that Figure 5 keeps its numbers at
  risk table.

Nothing in the preregistered stages was changed after the fact. The refuted
Stage J hypotheses are reported as refuted, in the manuscript and here.

## v1.1.0 — 2026-09-17

Stage D clinical robustness (`analysis/23_clinical.R`, results `CD00` to `CD08`)
and Stage E referee analyses (`analysis/25_referee.R`, results `CE01` onward),
plus the completed vector figure package. Both stages are post hoc and labeled
as such.

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
