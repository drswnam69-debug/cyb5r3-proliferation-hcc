# Changelog

## v1.1.0 — 2026-09-17

Adds the two stages written after the first release, and retargets the accompanying
manuscript from Liver International, which desk rejected it, to Scientific Reports.

- `analysis/23_clinical.R`, Stage D: post hoc clinical robustness. Alternative endpoints,
  liver function sensitivity, a bootstrap interval for the change in hazard ratio and
  absolute survival by tertile. Results `CD00` to `CD08`.
- `analysis/24_figure5.R`: Figure 5, Kaplan-Meier by CYB5R3 tertile with numbers at risk.
- `analysis/25_referee.R`, Stage E: analyses added in referee review. The Cox hazard ratio
  is not collapsible, so Stage E measures how much of the rise seen when proliferation enters
  the model is arithmetic rather than causal, by simulating from the fitted model with the
  proliferation score permuted. It also adds a baseline characteristics table, median
  follow-up by reverse Kaplan-Meier, a grade trend sensitivity check without the twelve
  grade 4 tumors, ESTIMATE stromal and immune scores and a HALLMARK_G2M_CHECKPOINT score.
  Results `CE00` to `CE06`.
- `RUN_C_ALL.R`: runs every stage in one pass.
- Figure width changed from 180 mm to 175 mm throughout, to fit a two-column page. EPS and
  PDF added to the figure formats. `figures/README_EPS.txt` records how they were made.
- Two labels in the result files were tightened. `CC01_ph_corrected_models.csv` no longer
  calls the stratified and time-varying models PRE-SPECIFIED, because the rule that governs
  them was written after the proportional hazards violation was detected, though before those
  models were fitted; they are now labeled RULE FIXED BEFORE FITTING. In
  `CD01_variable_coverage.csv` the column `modelled` is renamed `meets_coverage_threshold`,
  because passing the 60 percent coverage threshold makes a variable eligible for modeling
  and does not mean it entered any model reported here.
- British spellings in comments, console messages and result strings changed to American,
  matching the manuscript. Function and argument names from ggplot2 and dplyr are unchanged.

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
