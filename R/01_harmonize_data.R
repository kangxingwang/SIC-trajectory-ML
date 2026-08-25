## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset = getwd()), "R", "00_config.R"), encoding = "UTF-8")
assert_input_files(SIC_CONFIG$data_files)

## ============================================================================
## 2. Import data and harmonize variable names
## ============================================================================
RENAME_MAP <- c(
  stay_id="stay_id", class="Class", admission_era="Admission_Era", age="Age", gender="Gender", race="Race", icu_type="ICU_Type",
  sofa="SOFA", oasis="OASIS", charlson="Charlson", icd2_hypertension="Hypertension", icd2_diabetes="Diabetes", icd2_copd="COPD",
  icd2_hf="HF", icd2_stroke="Stroke", icd_malignancy="Malignancy", lung_infection="Lung_Infection", gi_infection="GI_Infection",
  genito_urinary_infection="GU_Infection", other_infection="Other_Infection", vs_24h_map_min="MAP", vs_24h_heart_rate_max="Heart_Rate",
  vs_24h_resp_rate_max="Resp_Rate", vs_24h_temp_max="Temperature", lab_in_icu_hemoglobin_first="Hb", lab_in_icu_wbc_first="WBC",
  lab_in_icu_platelet_first="Plt", lab_in_icu_creatinine_first="Cre", lab_in_icu_ast_first="AST", lab_in_icu_tbil_first="TBil",
  lab_in_icu_na_first="Na", lab_in_icu_k_first="K", lab_in_icu_ca_first="Ca", lab_in_icu_cl_first="Cl",
  baseline_score_plt_sic="Plt_S", baseline_score_inr_sic="INR_S", baseline_sofa_sic_points="SOFA_S",
  day28_los="Day28_LOS", day28_outcome="Day28_Outcome"
)
FEATURES <- c("Age","Gender","Race","ICU_Type","SOFA","OASIS","Charlson","Hypertension","Diabetes","COPD","HF","Stroke","Malignancy",
              "Lung_Infection","GI_Infection","GU_Infection","Other_Infection","MAP","Heart_Rate","Resp_Rate","Temperature","Hb","WBC","Plt",
              "Cre","AST","TBil","Na","K","Ca","Cl","Plt_S","INR_S","SOFA_S")
CAT_VARS <- c("Gender","Race","ICU_Type","Hypertension","Diabetes","COPD","HF","Stroke","Malignancy","Lung_Infection","GI_Infection","GU_Infection","Other_Infection")
NUM_VARS <- setdiff(FEATURES, CAT_VARS); FORCED_FEATURES <- c("Plt_S","INR_S","SOFA_S")

harmonize_one <- function(path, cohort) {
  x <- data.table::fread(path, data.table = FALSE, na.strings = c("","NA","NULL","NaN"), encoding = "UTF-8")
  if (!"stay_id" %in% names(x)) stop(cohort, ": stay_id is absent.")
  if (anyDuplicated(x$stay_id)) stop(cohort, ": stay_id is not unique.")
  needed_new <- c("stay_id","Class",FEATURES,"Day28_LOS","Day28_Outcome")
  needed_original <- names(RENAME_MAP)[match(needed_new, unname(RENAME_MAP))]
  if (cohort == "MIMIC-IV") needed_original <- unique(c(needed_original, "admission_era"))
  missing <- setdiff(needed_original, names(x)); if (length(missing)) stop(cohort, " missing fields: ", paste(missing, collapse = ", "))
  x <- x[, unique(needed_original), drop = FALSE]
  old <- intersect(names(RENAME_MAP), names(x)); names(x)[match(old,names(x))] <- unname(RENAME_MAP[old])
  x$Cohort <- cohort; x$Class <- suppressWarnings(as.integer(gsub("[^123]", "", as.character(x$Class))))
  if (any(!x$Class %in% 1:3 | is.na(x$Class))) stop(cohort, ": invalid Class values.")
  x$Target <- factor(ifelse(x$Class == 2, "Class2", "Other"), levels = c("Class2","Other"))
  x$Gender <- clean_gender(x$Gender); x$Race <- clean_race(x$Race); x$ICU_Type <- clean_icu(x$ICU_Type)
  for (v in setdiff(CAT_VARS,c("Gender","Race","ICU_Type"))) x[[v]] <- clean_binary(x[[v]])
  for (v in c(NUM_VARS,"Day28_LOS","Day28_Outcome")) x[[v]] <- suppressWarnings(as.numeric(as.character(x[[v]])))
  if (median(x$Day28_LOS, na.rm = TRUE) > 60) x$Day28_LOS <- x$Day28_LOS / 24
  if (cohort == "MIMIC-IV") {
    x$Admission_Era <- factor(trimws(as.character(x$Admission_Era)), levels = c("2008-2010","2011-2013","2014-2016","2017-2019","2020-2022"))
    if (anyNA(x$Admission_Era)) stop("MIMIC-IV: invalid or missing Admission_Era.")
  }
  x[, c("stay_id","Cohort","Class","Target",if(cohort=="MIMIC-IV") "Admission_Era",FEATURES,"Day28_LOS","Day28_Outcome"), drop = FALSE]
}
cohorts <- lapply(names(SIC_CONFIG$data_files), function(nm) harmonize_one(SIC_CONFIG$data_files[[nm]], nm)); names(cohorts) <- names(SIC_CONFIG$data_files)

## ============================================================================
## 3. Cohort and missingness audit
## ============================================================================
cohort_report <- bind_rows(lapply(cohorts, function(x) x %>% count(Cohort, Class, name = "N") %>% group_by(Cohort) %>% mutate(Percent = 100*N/sum(N)) %>% ungroup()))
missing_report <- bind_rows(lapply(cohorts, function(x) tibble(Cohort = unique(x$Cohort), Variable = FEATURES,
  Missing_N = vapply(x[FEATURES], function(z) sum(is.na(z)), integer(1)), Missing_Pct = 100*vapply(x[FEATURES], function(z) mean(is.na(z)), numeric(1)))))
bad <- missing_report %>% group_by(Variable) %>% summarise(Max_Missing_Pct = max(Missing_Pct), .groups = "drop") %>% filter(Max_Missing_Pct >= 20)
if (nrow(bad)) stop("Candidate predictors with >=20% missingness in at least one cohort: ", paste(bad$Variable, collapse = ", "))
stopifnot(!any(c("Day28_LOS","Day28_Outcome","Class","Target") %in% FEATURES), all(FORCED_FEATURES %in% FEATURES))

## ============================================================================
## 4. Missingness plot
## ============================================================================
pal_cohort <- c("MIMIC-IV"="#377EB8","MIMIC-III"="#E41A1C","eICU"="#FF7F00","NWICU"="#4DAF4A")
p_missing <- missing_report %>% mutate(Variable = factor(Variable, levels = rev(FEATURES))) %>%
  ggplot(aes(Missing_Pct, Variable, color = Cohort, shape = Cohort)) + geom_vline(xintercept = 20, linetype = 2, color = "#B2182B", linewidth = .6) +
  geom_point(size = 2.1, alpha = .9) + scale_color_manual(values = pal_cohort) + labs(title = "Cross-cohort predictor missingness", x = "Missingness (%)", y = NULL, color = NULL, shape = NULL) +
  theme_manuscript(10) + theme(legend.position = "top", axis.text.y = element_text(size = 7))
save_png(p_missing, "Supplement_Missingness.png", 8.3, 8.8)

## ============================================================================
## 5. Export aggregate tables and private harmonized object
## ============================================================================
fwrite(cohort_report, file.path(SIC_CONFIG$dirs$tables,"cohort_distribution.csv"), bom = TRUE)
fwrite(missing_report, file.path(SIC_CONFIG$dirs$tables,"missingness.csv"), bom = TRUE)
fwrite(data.frame(Original_Name = names(RENAME_MAP), Unified_Name = unname(RENAME_MAP)), file.path(SIC_CONFIG$dirs$tables,"variable_mapping.csv"), bom = TRUE)
fwrite(data.frame(Feature = FEATURES, Type = ifelse(FEATURES %in% CAT_VARS,"Categorical","Numeric"), Forced = FEATURES %in% FORCED_FEATURES), file.path(SIC_CONFIG$dirs$tables,"candidate_predictors.csv"), bom = TRUE)
saveRDS(list(cohorts=cohorts,features=FEATURES,cat_vars=CAT_VARS,num_vars=NUM_VARS,forced=FORCED_FEATURES,rename_map=RENAME_MAP), file.path(SIC_CONFIG$dirs$private,"harmonized_cohorts.rds"))
cat("Harmonization complete | cohorts:", paste(names(cohorts), collapse = ", "), "| candidate predictors:", length(FEATURES), "\n")

