# Does adjusting for tumor proliferation bias prognostic gene expression estimates?

Analysis code, preregistered decision rules, result files, run logs and figures for a study of
what one omitted covariate does to prognostic gene expression estimates, measured across every
expressed gene in a cohort and across thirteen cancer cohorts.

The work began as a single-gene question in hepatocellular carcinoma and became a methodological
one. It follows up Nam SW, *Divergence of the CYB5R3-mARC1 Redox Axis Across the Human
MASLD-to-Hepatocellular-Carcinoma Continuum*, Int J Mol Sci 2026;27:7806
(doi:10.3390/ijms27177806), whose code is deposited separately at doi:10.5281/zenodo.21987539.

## What the study did, in the order it did it

Each stage is labeled below as **preregistered**, **post hoc**, or **exploratory**, and the same
label appears in the code and in the output files. The distinction is the point of this
repository, not decoration.

**Stage 1, preregistered.** Does Hoshida molecular subclass explain the discordance between
histologic grade and survival for CYB5R3? Hypotheses H1 to H4, the decision rules, the subclass
assignment method and the validity controls were fixed in code before the data were examined.
All four hypotheses failed.

**Stage 2, post hoc.** Grade turns out not to be prognostic in TCGA-LIHC at all, so the
discordance was apparent rather than real. The grade trend tracks proliferation, and
proliferation negatively confounds the CYB5R3 association: adjusting for it raises the hazard
ratio rather than lowering it.

**Stage 3, preregistered.** That finding was tested in GSE14520 under a rule fixed before the
cohort was used for the purpose. It replicated. A proportional hazards violation found along the
way was handled under a further rule fixed in advance.

**Stage D, post hoc.** Clinical robustness: alternative endpoints, liver function, a paired
bootstrap for the change in hazard ratio, and absolute survival by tertile.

**Stage E, post hoc.** The Cox hazard ratio is not collapsible, so part of any rise on adjustment
is arithmetic. A simulation with proliferation made independent of every other covariate
measures how much.

**Stage F, exploratory.** The genome-wide sweep. For every expressed gene in TCGA-LIHC, two Cox
models are fitted, one with and one without a proliferation score. The ratio of the two hazard
ratios is the shift. The shift is regressed on the gene's Spearman correlation with the score.
Replicated in GSE14520 and repeated with HALLMARK_G2M_CHECKPOINT in place of the score.

**Stage G, post hoc.** Refinements to the sweep: non-collapsibility read off the genes that carry
no confounding, reversals counted among genes that carry an effect, and the consequences under
Benjamini-Hochberg control.

**Stage H, post hoc in timing, with its design fixed before any paper was read.** A survey of
50 recent prognostic gene expression papers in hepatocellular carcinoma, drawn under a fixed
seed from a frame of 424 records, coded for whether the survival model adjusted for
proliferation.

**Stage I, exploratory pilot.** The same sweep in TCGA-KIRC and TCGA-LUAD, with TCGA-LIHC
refitted under a common covariate base so the three are comparable.

**Stage J, PREREGISTERED.** The confirmatory pan-cancer sweep. The plan is
`PREREG_PanCancer_20260920.docx`, registered as doi:10.17605/OSF.IO/X5DCF on 2026-09-20, before
`analysis/31_pancancer_confirm.R` was run on 2026-09-21. Ten cancer types met the registered
inclusion rule. **Both registered hypotheses were refuted**, and they are reported as refuted.

**Stage K.** The figure for Stage J, which recomputes the numbers the manuscript quotes so that
the figure and the text cannot drift apart.

**Stage L, post hoc.** Is the eight-gene proliferation score a valid instrument in every cohort,
or only in the liver? Internal consistency, first principal component share, and agreement with
a 195-gene hallmark proliferation signature, computed on each sweep's own analysis set and
checked against that sweep's patient count.

**Stage M, post hoc.** Could measurement error account for the relationship that Stage J's
refutation revealed? Attenuation-corrected estimates and a test of whether reliability adds
anything once the proliferation hazard ratio is in the model.

**Stage O, post hoc.** Is the landscape a measurement or only a re-expression of the fitted
model? For each gene the first-order omitted variable expression predicts the shift from the
proliferation log hazard ratio and the gene's partial association with the score, obtained in
Frisch-Waugh-Lovell residual form. Two readings were fixed before the fits were seen, one for the
expression being an adequate description and one for the remainder being non-collapsibility.
**Neither holds**, and the departure is reported as unexplained by either account.

A second plan, `PREREG2_ConfoundingLaw_20260920.docx`, doi:10.17605/OSF.IO/N8DW6, registers the
test of the post hoc relationship in cohorts that have not been examined. **No analysis covered
by that registration has been run, and none of it is in this repository.**

## Requirements

- R 4.x with `survival`, `dplyr`, `data.table`, `stringr`, `ggplot2`, `ggpubr`, `ggrepel`,
  `UCSCXenaTools`, `GEOquery`, `Biobase`, `msigdbr`
- `clinfun` for the Jonckheere-Terpstra tests, `ragg` for LZW-compressed TIFF output
- Internet access on first run: TCGA matrices are fetched from UCSC Xena and GSE14520 from GEO,
  then cached. Downloads land in `cache/` beside this README unless the environment variable
  `CYB5R3_CACHE` names another directory, or a checkout of the earlier pipeline sits beside this
  repository, in which case its cache is reused. The console prints the cache in use on every run.
- Stage J downloads an expression matrix for every TCGA cancer type it screens, roughly one
  gigabyte in total. Downloads are cached, so a second run is fitting only.

## How to run

Open one of these in RStudio and press Source. Each one runs everything before it.

| File | Runs |
|---|---|
| `RUN_C.R` | Stage 1 only |
| `RUN_C_ADDENDUM.R` | + Stage 2 |
| `RUN_C_VALIDATION.R` | + Stage 3 |
| `RUN_C_PH.R` | + proportional hazards handling |
| `RUN_C_FIGURES.R` | + figures |
| `RUN_C_ALL.R` | Stages 1 to E and all figures. Run this first |

Stages F onward are sourced individually, in this order, after `RUN_C_ALL.R`:

```
analysis/26_genomewide.R        Stage F
analysis/27_sweep_refine.R      Stage G
analysis/28_figure6.R           the landscape figure
analysis/29_litsurvey.R         Stage H   (run twice; see below)
analysis/30_pancancer.R         Stage I
analysis/31_pancancer_confirm.R Stage J   PREREGISTERED
analysis/32_figure_prereg.R     Stage K
analysis/33_score_validity.R    Stage L
analysis/34_attenuation.R       Stage M
analysis/36_approximation.R     Stage O
```

`29_litsurvey.R` draws the sample and stops, because the coding sheet does not yet exist. The
coded sheet is `results/CH02_litsurvey_coded.csv`; with it present, a second run tabulates.

## Layout

```
analysis/       18 to 25   Stages 1 to E
                26 to 36   Stages F to O
                grp/       the three Hoshida gene sets as downloaded from MSigDB v2026.1.Hs
                hoshida_templates.tsv  those gene sets frozen on first use
                vendor/01_common.R     helper file from the published pipeline, MIT, same author
results/        every result file, run log and session record
figures/        PNG (300 dpi), TIFF (LZW), SVG, EPS and PDF
PREREG_ProjectC_subtype.docx      Stage 1 preregistration
PREREG_PanCancer_20260920.docx    Stage J preregistration, doi:10.17605/OSF.IO/X5DCF
PREREG2_ConfoundingLaw_20260920.docx  the next study, doi:10.17605/OSF.IO/N8DW6
```

## Where the numbers come from

| Prefix | Holds | Status |
|---|---|---|
| `C00` to `C99` | Stage 1, subclass assignment, validity controls, H1 to H4, verdict | preregistered |
| `CA0` to `CA5` | Stage 2 | post hoc |
| `CB00` to `CB08` | Stage 3 replication | preregistered |
| `CB09` to `CB12` | tissue composition, which did not replicate | exploratory |
| `CC00` to `CC99` | proportional hazards handling and its verdict | pre-specified |
| `CD00` to `CD08` | clinical robustness | post hoc |
| `CE01` | the non-collapsibility simulation | post hoc |
| `CF00` to `CF10` | the genome-wide sweep, its replication and its consequences | exploratory |
| `CG00` to `CG06` | null band, reversals, Benjamini-Hochberg consequences, G2M check | post hoc |
| `CH00` to `CH05` | the literature survey, its frame, draw, coding sheet and result | design fixed in advance |
| `CI00` to `CI02` | the two-cohort pilot and the comparable LIHC refit | exploratory |
| `CJ00` to `CJ05` | the confirmatory sweep, exclusions, per-cohort results, and the verdict | PREREGISTERED |
| `CK01`, `CK02` | the Stage J figure and the numbers it recomputes | post hoc |
| `CL00` to `CL03` | instrument validity by cohort and by slope sign | post hoc |
| `CM01` to `CM03` | attenuation correction and whether reliability adds anything | post hoc |
| `CO00` to `CO06` | how far the first-order omitted variable expression accounts for the shift | post hoc |

`CJ05_decision.csv` records the preregistered verdict. Both hypotheses read REFUTED.

## Reproducibility notes

- Seeds are fixed. A run log and session record are written for every execution.
- Every script writes its own decision rules to a result file on each run, so an edit made after
  the fact is visible in the deposited output rather than only in the code history.
- Stage L rebuilds each sweep's analysis set by the same rule and checks the patient count
  against that sweep's own output, so a validity statistic computed on a different set would be
  reported rather than passed over.
- The subclass templates are frozen to `analysis/hoshida_templates.tsv` on first use. The
  original `.grp` downloads are included so the freeze can be checked.
- On macOS without XQuartz, `grDevices::tiff(type = "cairo")` falls back to an uncompressed
  device with only a warning, so TIFF output is written with `ragg` and verified by reading its
  compression tag. For the same reason `grDevices::cairo_ps()` is unavailable and the EPS written
  by the `postscript()` fallback is a placeholder; the deposited EPS and PDF files were converted
  from the archived SVG. See `figures/README_EPS.txt`.

## Data

All data are public. TCGA expression, clinical and survival tables come from UCSC Xena; GSE14520
comes from GEO; the CTNNB1 mutation calls come from the PanCanAtlas MC3 matrix; gene sets come
from MSigDB. No new human or animal data were generated and no identifiable data are included.

## Citation

See `CITATION.cff`. Cite the Zenodo DOI for the version you used alongside the article.

## License

MIT. See `LICENSE`.
