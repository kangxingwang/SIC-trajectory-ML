# Contributing

Please open an issue before changing the statistical design. Pull requests
should preserve the no-leakage contract: all preprocessing, feature selection,
and tuning must be estimated inside training folds, and external cohorts must
never be used for model selection or refitting.

Before submitting a change:

1. Parse all R files and run the fast smoke test.
2. Confirm no patient-level data, local paths, `.Renviron`, RDS files, or model
   bundles are staged by Git.
3. Regenerate aggregate tables and figures with the fixed seed.
4. Document any change to candidate predictors, resampling, metrics, or model
   selection in the pull-request description.

