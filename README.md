# Tumor proliferation masks the prognostic association of CYB5R3 in hepatocellular carcinoma

Pre-registered analysis code, decision rules, results and figures for a two-cohort study of
CYB5R3 in hepatocellular carcinoma (HCC), using TCGA-LIHC and GSE14520.

This work follows up Nam SW, *Divergence of the CYB5R3-mARC1 Redox Axis Across the Human
MASLD-to-Hepatocellular-Carcinoma Continuum*, Int J Mol Sci 2026;27:7806
(doi:10.3390/ijms27177806), whose code is deposited separately at doi:10.5281/zenodo.21987539.
That article reported that CYB5R3 expression fell across histologic grade while higher
expression marked worse survival, and stated that the two directions could not be reconciled.
This repository contains the analysis that resolves it.

## What the study did, in the order it did it

**Stage 1, pre-registered.** Does Hoshida molecular subclass explain the discordance?
Hypotheses H1 to H4, the decision rules, the subclass assignment method and the validity
controls were fixed in code before the data were examined. All four hypotheses failed.

**Stage 2, post hoc.** Designed after those results were seen, and labeled POST HOC
throughout. Histologic grade turns out not to be prognostic in TCGA-LIHC at all, so the
discordance was apparent rather than real. The grade trend is a proliferation trend, and
proliferation negatively confounds the association between CYB5R3 and survival: adjusting
for it raises the hazard ratio rather than lowering it.

**Stage 3, pre-registered.** That finding was then tested in GSE14520 under a rule fixed
before the cohort was used for the purpose. It replicated. A proportional hazards violation
found along the way was handled under a further rule fixed in advance.

An exploratory tissue composition analysis, added after Stage 2, did **not** replicate in
GSE14520 and is reported as such.

## Requirements

- R 4.x with `survival`, `dplyr`, `data.table`, `stringr`, `ggplot2`, `ggpubr`,
  `UCSCXenaTools`, `GEOquery`, `Biobase`
- `clinfun` for the Jonckheere-Terpstra tests, `ragg` for LZW-compressed TIFF output
- Internet access on first run: TCGA matrices are fetched from UCSC Xena and GSE14520 from GEO,
  then cached. Downloads land in `cache/` beside this README unless the environment variable
  `CYB5R3_CACHE` names another directory, or a checkout of the earlier pipeline sits beside this
  repository, in which case its cache is reused. The console prints the cache in use on every run.

## How to run

Open one of these in RStudio and press Source. Each one runs everything before it.

| File | Runs |
|---|---|
| `RUN_C.R` | Stage 1 only (pre-registered subclass test) |
| `RUN_C_ADDENDUM.R` | + Stage 2 (post hoc) |
| `RUN_C_VALIDATION.R` | + Stage 3 (GSE14520 replication) |
| `RUN_C_PH.R` | + proportional hazards handling |
| `RUN_C_FIGURES.R` | + figures. This is the full pipeline |
| `RUN_C_FIGURES_ONLY.R` | figures only, if the objects are already in the session |

## Layout

```
analysis/       18_subtype.R  Stage 1, pre-registered
                19_addendum.R Stage 2, post hoc
                20_validation.R Stage 3, pre-registered replication + exploratory composition
                21_ph.R       proportional hazards, pre-specified
                22_figures.R  manuscript figures
                grp/          the three Hoshida gene sets as downloaded from MSigDB v2026.1.Hs
                hoshida_templates.tsv  those gene sets frozen on first use
                vendor/01_common.R     helper file from the published pipeline, MIT, same author
results/        every result file, run logs and session information
figures/        Figure 1 to 4 as PNG (300 dpi), TIFF (LZW) and SVG
PREREG_ProjectC_subtype.docx   the pre-registration
```

## Where the numbers come from

| File | Holds |
|---|---|
| `C00_decision_rules.csv` | the pre-registered rules, rewritten on every run |
| `C01_flow.csv`, `C02_anchor_reproduction.csv` | patient flow, and reproduction of the published estimate |
| `C03*`, `C04`, `C05`, `C06` | subclass assignment, template coverage, validity and positive controls |
| `C10` to `C50` | H1 to H4 and the multiplicity adjustment |
| `C99_verdict.csv` | Stage 1 verdict |
| `CA0` to `CA5` | Stage 2, all POST HOC |
| `CB00` to `CB08` | Stage 3, the pre-registered replication |
| `CB09` to `CB12` | the exploratory composition analysis, which did not replicate |
| `CC00` to `CC99` | proportional hazards handling and its verdict |

## Reproducibility notes

- The seed is fixed. Session information and a run log are written for every execution.
- The subclass templates are frozen to `analysis/hoshida_templates.tsv` on first use, so later
  runs do not depend on the availability of an external gene set server. The original `.grp`
  downloads are included so the freeze can be checked.
- Gene counts parsed from those files are compared against the counts printed on the MSigDB
  gene set pages, and the comparison is written to `results/C03a_template_gene_counts.csv`.
- Control markers are checked against the templates on every run, and any that belong to a
  template are excluded automatically. This is why AFP and APOA1 are not among the controls.
- On macOS without XQuartz, `grDevices::tiff(type = "cairo")` falls back to an uncompressed
  device with only a warning. `analysis/22_figures.R` therefore uses `ragg` and verifies the
  written TIFF by reading its compression tag.

## Data

All data are public. TCGA-LIHC expression, clinical and survival tables come from UCSC Xena;
GSE14520 comes from GEO; the CTNNB1 mutation calls come from the PanCanAtlas MC3 matrix; the
subclass templates come from MSigDB. No new human or animal data were generated and no
identifiable data are included here.

## Citation

See `CITATION.cff`. The DOI assigned by Zenodo on publication of this record should be cited
alongside the article.

## License

MIT. See `LICENSE`.
