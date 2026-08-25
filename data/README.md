# Data access and expected files

The source cohorts are **not distributed in this repository**. MIMIC and eICU
data are credentialed/restricted resources, and the merged study files contain
patient-level records. Keep all source data outside Git and provide their local
paths through `.Renviron`.

Expected inputs:

| Environment variable | Cohort | Required special field |
|---|---|---|
| `SIC_DATA_MIMICIV` | MIMIC-IV | `admission_era` for the temporal split |
| `SIC_DATA_MIMICIII` | MIMIC-III | none beyond the common data dictionary |
| `SIC_DATA_EICU` | eICU | none beyond the common data dictionary |
| `SIC_DATA_NWICU` | NWICU | none beyond the common data dictionary |

Each row must represent one unique ICU stay (`stay_id`). The target column is
`class`, coded 1, 2, or 3. The analysis predicts Class 2 versus Classes 1/3.

The exact original-to-analysis variable mapping is defined and exported by
`R/01_harmonize_data.R`.
