# Changelog

## v2.3.0 - 2026-09-30

One analysis was added, in answer to an objection raised during pre-submission
review: that the landscape is a re-expression of the fitted model rather than a
measurement of anything. Nothing else changed. No model, no decision rule and no
estimate from an earlier stage was altered, and every number the earlier releases
hold is reproduced by this one.

Added:

- `analysis/36_approximation.R`, Stage O, post hoc. How much of the per-gene
  shift the first-order omitted variable expression accounts for, and what is
  left over. For each gene the expression predicts the log shift as minus the log
  hazard ratio of the proliferation score times the coefficient of the gene in a
  regression of the score on that gene and the same covariates, obtained in
  Frisch-Waugh-Lovell residual form so that it is exact for all 16,902 genes at
  once rather than iterative. Three readings were fixed before the fits were seen
  and are written to `CO00` on every run: the expression is an adequate
  description if it explains at least 0.95 of the variance of the observed shift
  with a slope within 0.05 of one; what it omits is non-collapsibility if the
  remainder is proportional to the gene's own log hazard ratio with a coefficient
  close to the independent Stage G estimate of 0.064; if neither holds, the
  departure is reported as unexplained. Results `CO00` to `CO06`.

Both pre-specified readings failed. In TCGA-LIHC the prediction explains 0.8386
of the variance of the observed shift with a slope of 0.8968, short of both
adequacy thresholds, and the remainder rises with the gene's own log hazard ratio
at 0.0977 per unit (0.0945 to 0.1010), which is not the 0.0643 that
non-collapsibility alone would give. GSE14520 behaves the same way: 0.8201 of the
variance, a slope of 0.8669, and 0.0531 per unit (0.0509 to 0.0554). Adding the
gene's own coefficient to the prediction raises the variance explained only to
0.8594 and 0.8375. The departure is therefore reported as unexplained by either
account, which is the informative outcome: the objection is right about the slope
of the landscape and wrong about the per-gene shift, which carries something
neither the first-order expression nor non-collapsibility accounts for.

The ratio of the remainder coefficient to the proliferation log hazard ratio is
0.260 (0.252 to 0.269) in TCGA-LIHC and 0.245 (0.235 to 0.256) in GSE14520. That
regularity is recorded as an observation for a later test and not as a result.

Analyses ran under R 4.6.1 with survival 3.8-6.

## v2.2.0 - 2026-09-29

A reporting and provenance release, prepared during a full pre-submission audit
of the manuscript. No model, no decision rule and no estimate changed. What
changed is which cohorts the primary analysis uses, and how much of what the
manuscript quotes can be checked against a result file.

The confirmatory set. The registered rule was applied to the cohorts as The
Cancer Genome Atlas distributes them, and that distribution contains aggregate
cohorts alongside single-organ ones. Ten cohorts met the rule, but COADREAD is
COAD together with rectal adenocarcinoma and LUNG is LUSC together with LUAD, so
the ten covered only eight disjoint groups of patients. Rectal adenocarcinoma had
been excluded on its own for having 84 patients and lung adenocarcinoma had been
excluded by name, and each re-entered inside an aggregate. The primary analysis
is now the eight cohorts in which no patient appears twice, with the original ten
reported as a sensitivity analysis. Neither registered hypothesis changes verdict:
the slope is negative in four of eight and five of ten, against a refutation
threshold of 0.75, and the median pairwise agreement is -0.268 and -0.293 against
a refutation threshold of 0.10. This departure from the registered set was
identified after the analysis had been run and is labeled post hoc.

Changed:

- `analysis/26_genomewide.R` now records the analysis set of each sweep. The
  discovery sweep adjusts for histologic grade, so it is fitted in the 339
  patients with a recorded grade, among whom there are 114 deaths, and not in the
  341 of the cohort. `CF03` and `CF08` carry those counts, and `CF08` now also
  carries the slope standard error and the variance explained, so the GSE14520
  R-squared of 0.7959 that the manuscript quotes no longer exists only inside a
  figure subtitle.
- `analysis/27_sweep_refine.R` records the size of the Hoshida-restricted
  HALLMARK_G2M_CHECKPOINT score, 172 genes, in `CG06`.
- `analysis/33_score_validity.R` records the size of the hallmark comparison set,
  195 genes, in the new `CL04`.
- `analysis/28_figure6.R` writes the slope, its standard error and the variance
  explained for both landscape panels to the new `CF11`, and draws at 170 mm.
- `analysis/32_figure_prereg.R` marks the two aggregate cohorts apart from the
  eight primary ones in both panels, fits panel B to the primary eight, and
  reports the regression for four sets in `CK02`: the primary eight, the
  registered ten, the eight with the three exploratory cohorts, and all thirteen.
  It draws at 170 mm.
- Both figures were regenerated and their EPS files rebuilt from the SVG through
  `rsvg-convert` and Ghostscript, because `cairo_ps()` is unavailable on the
  analysis machine and the `postscript()` fallback clips text.

Reproduced unchanged on this run: the discovery slope -0.2782 with R-squared
0.7851 over 16,902 genes, the GSE14520 slope -0.1798 over 10,861 genes, the
threshold counts 2308, 1333, 338 and 1, the reversal fractions 18.24 and 3.74
percent, the null band of 2963 genes with an inflation of 0.0643 per unit
coefficient, and every cohort-level estimate in `CJ03` and `CI02`.

Analyses ran under R 4.6.1 with survival 3.8-6.

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
