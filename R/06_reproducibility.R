## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset=getwd()),"R","00_config.R"),encoding="UTF-8")

## ============================================================================
## 2. Import final aggregate outputs for QA
## ============================================================================
required_public <- c(
  "results/tables/feature_engineering_summary.csv","results/tables/nested_oof_model_ranking.csv",
  "results/tables/performance_summary.csv","results/tables/shap_global_importance.csv",
  "results/figures/Figure_ML1_Feature_Selection.png","results/figures/Figure_ML2A_ROC_MIMICIV_Temporal.png",
  "results/figures/Figure_ML2E_XGBoost_DCA.png","results/figures/Figure_ML2F_XGBoost_Calibration.png",
  "results/figures/Figure_ML2G_XGBoost_SHAP_Importance.png","results/figures/Figure_ML2H_XGBoost_SHAP_Beeswarm.png"
)
missing_public <- required_public[!file.exists(file.path(SIC_CONFIG$project_root,required_public))]
if(length(missing_public)) stop("Reproducibility QA failed; missing outputs:\n",paste(missing_public,collapse="\n"))

## ============================================================================
## 3. Reproducibility metadata and source-file fingerprints
## ============================================================================
core <- readRDS(file.path(SIC_CONFIG$dirs$private,"modeling_core.rds"))
settings <- tibble(Setting=c("Seed","Outer folds","Inner folds","Tuning grid size","Boruta max runs","Bootstrap replicates","SHAP simulations","SHAP sample per cohort","Workers","Fast test","Locked model","Final features"),
  Value=c(SIC_CONFIG$seed,SIC_CONFIG$outer_v,SIC_CONFIG$inner_v,SIC_CONFIG$grid_n,SIC_CONFIG$boruta_runs,SIC_CONFIG$bootstrap_n,SIC_CONFIG$shap_nsim,SIC_CONFIG$shap_sample,SIC_CONFIG$workers,SIC_CONFIG$fast_test,core$best_model,paste(core$final_features,collapse="; ")))
source_fingerprints <- tibble(Cohort=names(SIC_CONFIG$data_files),File_Name=basename(SIC_CONFIG$data_files),Size_Bytes=as.numeric(file.info(SIC_CONFIG$data_files)$size),
  SHA256=vapply(SIC_CONFIG$data_files,digest::digest,character(1),file=TRUE,algo="sha256"))
package_versions <- tibble(Package=required_packages,Version=vapply(required_packages,function(p)as.character(utils::packageVersion(p)),character(1)))
fwrite(settings,file.path(SIC_CONFIG$dirs$tables,"analysis_settings.csv"),bom=TRUE)
fwrite(source_fingerprints,file.path(SIC_CONFIG$dirs$tables,"source_file_fingerprints.csv"),bom=TRUE)
fwrite(package_versions,file.path(SIC_CONFIG$dirs$tables,"package_versions.csv"),bom=TRUE)
writeLines(c(capture.output(sessionInfo()),"",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z"))),file.path(SIC_CONFIG$dirs$logs,"sessionInfo.txt"),useBytes=TRUE)

## ============================================================================
## 4. Public-output inventory
## ============================================================================
all_files <- list.files(SIC_CONFIG$project_root,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)
info <- file.info(all_files); files <- all_files[!info$isdir]; rel <- gsub("\\\\","/",substring(files,nchar(SIC_CONFIG$project_root)+2L))
keep <- !grepl("^(private/|data/raw/|renv/library/|\\.git/|node_modules/|\\.Renviron$)",rel); files <- files[keep]; rel <- rel[keep]; info <- file.info(files)
manifest <- tibble(Relative_Path=rel,Size_Bytes=as.numeric(info$size),Modified=as.character(info$mtime),SHA256=vapply(files,digest::digest,character(1),file=TRUE,algo="sha256")) %>%
  filter(Relative_Path!="results/tables/public_file_manifest_sha256.csv") %>% arrange(Relative_Path)
fwrite(manifest,file.path(SIC_CONFIG$dirs$tables,"public_file_manifest_sha256.csv"),bom=TRUE)

## ============================================================================
## 5. Final QA report
## ============================================================================
qa <- c("SIC trajectory ML reproducibility report",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),paste0("Locked model: ",core$best_model),
  paste0("Final features: ",paste(core$final_features,collapse=", ")),paste0("Public files inventoried: ",nrow(manifest)),
  "Patient-level data and fitted model objects are restricted to private/ and ignored by Git.")
writeLines(qa,file.path(SIC_CONFIG$dirs$logs,"reproducibility_report.txt"),useBytes=TRUE)
cat("Reproducibility QA passed | public files inventoried:",nrow(manifest),"\n")
