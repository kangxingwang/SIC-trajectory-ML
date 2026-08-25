options(stringsAsFactors = FALSE, scipen = 999)

command <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", command[grep("^--file=", command)])
project_root <- if (length(file_arg)) dirname(normalizePath(file_arg[1], winslash = "/")) else normalizePath(getwd(), winslash = "/")
Sys.setenv(SIC_PROJECT_ROOT = project_root)

scripts <- c(
  "R/00_config.R",
  "R/01_harmonize_data.R",
  "R/02_feature_engineering.R",
  "R/03_nested_modeling.R",
  "R/04_model_evaluation.R",
  "R/05_shap_analysis.R",
  "R/06_reproducibility.R"
)

started <- Sys.time()
for (script in scripts) {
  cat(sprintf("\n===== %s =====\n", script))
  source(file.path(project_root, script), encoding = "UTF-8", echo = FALSE)
}
if (tolower(Sys.getenv("SIC_BUILD_WORKBOOKS", unset = "false")) %in% c("1","true","yes")) {
  node <- Sys.getenv("SIC_NODE", unset = Sys.which("node"))
  if (!nzchar(node)) warning("SIC_BUILD_WORKBOOKS=true but Node.js was not found; CSV tables are complete and remain the canonical outputs.")
  else system2(node, file.path(project_root,"scripts","build_workbooks.mjs"), env = paste0("SIC_PROJECT_ROOT=",project_root))
}
cat(sprintf("\nPipeline completed in %.1f minutes.\n", as.numeric(difftime(Sys.time(), started, units = "mins"))))
