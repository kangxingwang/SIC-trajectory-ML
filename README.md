# SIC trajectory machine-learning analysis

Reproducible feature engineering, nested machine learning, temporal/external
validation, and SHAP interpretation for classifying the persistently high SIC
trajectory (`Class 2`) versus `Class 1/3`.

MIMIC-IV is the development database. Admissions from 2008–2019 form the
historical development cohort; 2020–2022 is held out for temporal validation.
MIMIC-III, eICU, and NWICU are locked external validation cohorts.

## Analysis contract

- This is a **trajectory-classification model**, not a mortality model.
- Predictors must exist in all four databases, have `<20%` missingness in every
  database, and be available within 24 h of ICU admission.
- Treatment, outcome, and follow-up variables are excluded from predictors.
- `Plt_S`, `INR_S`, and `SOFA_S` are prespecified SIC components. The total SIC
  score is excluded to avoid deterministic redundancy.
- Boruta and LASSO selection are repeated inside each outer training fold.
  Their stable intersection is combined with the three prespecified components.
- Preprocessing is fold-specific: training-fold medians for numeric variables,
  explicit unknown/novel levels for categorical variables, one-hot encoding,
  and model-appropriate normalization.
- Five candidates are compared: elastic net, random forest, XGBoost, radial
  SVM, and KNN. The manuscript analysis uses 10 outer × 10 inner folds.
- The winning model is selected using nested out-of-fold AUROC and locked before
  temporal/external validation. The classification threshold is also frozen
  from development OOF predictions.
- Evaluation includes AUROC/AUPRC, bootstrap confidence intervals,
  calibration intercept/slope, Brier score, ECE/MCE, threshold metrics,
  DeLong tests with BH-FDR correction, decision-curve analysis, and SHAP.

The completed manuscript run selected **XGBoost** with 10 final predictors:
`Plt_S`, `INR_S`, `SOFA_S`, `MAP`, `Na`, `SOFA`, `TBil`, `Cre`,
`Temperature`, and `GI_Infection`.

## Project structure

```text
R/00_config.R              packages, paths, seeds, helpers
R/01_harmonize_data.R      import, rename, harmonize, missingness audit
R/02_feature_engineering.R nested Boruta/LASSO selection and Figure ML1
R/03_nested_modeling.R     nested tuning and locked external prediction
R/04_model_evaluation.R    metrics, ROC, calibration, and DCA
R/05_shap_analysis.R       global/stability/local SHAP and calculator bundle
R/06_reproducibility.R     hashes, package versions, session information, QA
scripts/build_workbooks.mjs formatted aggregate Excel workbooks
app/app.R                  local Shiny calculator for the locked model
run_all.R                  one-command execution
results/figures/           publication PNGs only
results/tables/            canonical machine-readable aggregate CSVs
results/workbooks/         three formatted aggregate Excel workbooks
private/                   patient-level objects; ignored by Git
```

Each R module is organized as: packages/configuration → data import → analysis
→ plotting → table export.

## Data and environment variables

Raw ICU data are not distributed. Copy `.Renviron.example` to `.Renviron` and
set these local paths:

```text
SIC_DATA_MIMICIV
SIC_DATA_MIMICIII
SIC_DATA_EICU
SIC_DATA_NWICU
```

The expected variable names are documented in
`results/tables/variable_mapping.csv`. `.Renviron` is ignored by Git, so local
private paths are not uploaded.

## Reproduce the analysis

Install R 4.4 or later, then restore the recorded environment:

```r
install.packages("renv")
renv::restore()
```

Run a quick end-to-end smoke test:

```r
Sys.setenv(SIC_FAST_TEST = "true", SIC_RESUME = "false")
source("run_all.R")
```

Run the manuscript analysis in a fresh session:

```r
Sys.setenv(SIC_FAST_TEST = "false", SIC_RESUME = "true")
source("run_all.R")
```

The full nested analysis is computationally intensive. With `SIC_RESUME=true`,
completed fold-level checkpoints are reused after interruption. The fixed seed
is `19971117` unless overridden in `.Renviron`.

Aggregate CSV files are the canonical tables. Formatted Excel copies are
included for manuscript editing. In a Codex runtime with `@oai/artifact-tool`,
set `SIC_BUILD_WORKBOOKS=true` to regenerate the three workbooks.

## Local research calculator

After the full pipeline creates `private/model_bundle.rds`:

```r
shiny::runApp("app")
```

The calculator is for research demonstration only. Its output is the
probability of `Class 2`; it is not a mortality risk, causal estimate, or
treatment recommendation.

## Privacy and public release

- Never commit raw data, `.Renviron`, patient-level predictions, fitted RDS
  objects, SHAP background rows, or calculator bundles.
- Public workbooks contain aggregate results only. Representative local SHAP
  cases use descriptive case labels and contain no `stay_id`.
- Before release, replace the placeholder author/email in `DESCRIPTION`, choose
  a final license, and add the article DOI to the repository citation metadata.
- SHAP values describe predictive associations and must not be interpreted
  causally.

## Current public outputs

- `Feature_Engineering.xlsx`: one-page nested Boruta/LASSO summary.
- `Model_Performance.xlsx`: model ranking, hyperparameters, performance with
  confidence intervals, DeLong-FDR comparisons, and frozen thresholds.
- `SHAP_Interpretation.xlsx`: cross-cohort importance, stability, and
  de-identified representative cases.
- Nine final PNG figures: feature selection; four cohort-specific ROC panels;
  XGBoost DCA; XGBoost calibration; SHAP importance; and SHAP beeswarm.

