## ============================================================================
## 0. Packages and reproducible configuration
## ============================================================================
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "tibble", "forcats", "stringr",
  "ggplot2", "scales", "patchwork", "tidymodels", "Boruta", "glmnet",
  "ranger", "xgboost", "kernlab", "kknn", "future", "pROC", "PRROC",
  "fastshap", "shapviz", "digest", "viridisLite"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop(
  "Missing packages: ", paste(missing_packages, collapse = ", "),
  "\nRun renv::restore() from the project root before analysis."
)
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(scales); library(tidymodels)
})
tidymodels::tidymodels_prefer()

`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "")) y else x
env_int <- function(name, default) as.integer(Sys.getenv(name, unset = as.character(default)))
env_bool <- function(name, default = FALSE) tolower(Sys.getenv(name, unset = ifelse(default, "true", "false"))) %in% c("1", "true", "yes", "y")
project_root <- normalizePath(Sys.getenv("SIC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = FALSE)
fast_test <- env_bool("SIC_FAST_TEST", FALSE)

SIC_CONFIG <- list(
  project_root = project_root,
  data_files = c(
    `MIMIC-IV`  = Sys.getenv("SIC_DATA_MIMICIV"),
    `MIMIC-III` = Sys.getenv("SIC_DATA_MIMICIII"),
    eICU        = Sys.getenv("SIC_DATA_EICU"),
    NWICU       = Sys.getenv("SIC_DATA_NWICU")
  ),
  seed = env_int("SIC_SEED", 19971117L),
  outer_v = if (fast_test) 2L else env_int("SIC_OUTER_FOLDS", 10L),
  inner_v = if (fast_test) 2L else env_int("SIC_INNER_FOLDS", 10L),
  grid_n = if (fast_test) 2L else env_int("SIC_GRID_SIZE", 12L),
  boruta_runs = if (fast_test) 20L else env_int("SIC_BORUTA_RUNS", 100L),
  bootstrap_n = if (fast_test) 100L else env_int("SIC_BOOTSTRAP_N", 1000L),
  shap_nsim = if (fast_test) 8L else env_int("SIC_SHAP_NSIM", 64L),
  shap_sample = if (fast_test) 80L else env_int("SIC_SHAP_SAMPLE", 300L),
  workers = max(1L, env_int("SIC_WORKERS", min(6L, max(1L, parallel::detectCores() - 2L)))),
  resume = env_bool("SIC_RESUME", TRUE), fast_test = fast_test
)
SIC_CONFIG$dirs <- list(
  results = file.path(project_root, "results"), figures = file.path(project_root, "results", "figures"),
  tables = file.path(project_root, "results", "tables"), workbooks = file.path(project_root, "results", "workbooks"),
  logs = file.path(project_root, "results", "logs"), private = file.path(project_root, "private")
)
invisible(lapply(SIC_CONFIG$dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

assert_input_files <- function(paths) {
  absent_env <- names(paths)[!nzchar(paths)]
  if (length(absent_env)) stop("Unset data environment variables for: ", paste(absent_env, collapse = ", "), "\nCopy .Renviron.example to .Renviron and fill in local paths.")
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files)) stop("Input files not found:\n", paste(missing_files, collapse = "\n"))
}
mode_value <- function(x) { z <- na.omit(as.character(x)); if (!length(z)) return(NA_character_); names(sort(table(z), decreasing = TRUE))[1] }
clean_binary <- function(x) {
  z <- trimws(tolower(as.character(x)))
  factor(ifelse(z %in% c("1","yes","y","true","t"), "Yes", ifelse(z %in% c("0","no","n","false","f"), "No", NA_character_)), levels = c("No","Yes"))
}
clean_gender <- function(x) {
  z <- trimws(toupper(as.character(x)))
  factor(ifelse(z %in% c("M","MALE","1"), "Male", ifelse(z %in% c("F","FEMALE","0"), "Female", NA_character_)), levels = c("Female","Male"))
}
clean_race <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- case_when(grepl("white",z) ~ "White", grepl("black|afric",z) ~ "Black", grepl("asian",z) ~ "Asian", is.na(z) | z == "" ~ NA_character_, TRUE ~ "Other")
  factor(out, levels = c("White","Black","Asian","Other"))
}
clean_icu <- function(x) {
  z <- trimws(tolower(as.character(x)))
  out <- case_when(grepl("coronary|cardiac|ccu",z) ~ "CCU", grepl("neuro",z) ~ "Neuro ICU", grepl("surg|sicu|trauma",z) ~ "Surgical ICU", grepl("med|micu",z) ~ "Medical ICU", is.na(z) | z == "" ~ NA_character_, TRUE ~ "Other ICU")
  factor(out, levels = c("Medical ICU","Surgical ICU","CCU","Neuro ICU","Other ICU"))
}
save_png <- function(plot, filename, width, height, dpi = 600) ggplot2::ggsave(file.path(SIC_CONFIG$dirs$figures, filename), plot, width = width, height = height, dpi = dpi, bg = "white")
theme_manuscript <- function(base_size = 10) theme_classic(base_size = base_size) + theme(axis.text = element_text(color = "black"), plot.title = element_text(face = "bold", hjust = .5))

set.seed(SIC_CONFIG$seed)
public_config <- SIC_CONFIG
public_config$project_root <- "."
public_config$data_files <- setNames(paste0(names(SIC_CONFIG$data_files), ".csv"), names(SIC_CONFIG$data_files))
public_config$dirs <- lapply(SIC_CONFIG$dirs, function(x) gsub("^.*?/(results|private)(/.*)?$", "\\1\\2", x))
writeLines(capture.output(str(public_config)), file.path(SIC_CONFIG$dirs$logs, "analysis_configuration.txt"), useBytes = TRUE)
cat("Configuration loaded | seed:", SIC_CONFIG$seed, "| outer/inner folds:", SIC_CONFIG$outer_v, "/", SIC_CONFIG$inner_v, "| fast test:", SIC_CONFIG$fast_test, "\n")
