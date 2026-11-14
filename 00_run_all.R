# Run from the repository root: Rscript 00_run_all.R
# Configuration: DP_DATA_DIR, DP_OUTPUT_DIR, DP_TASKS, DP_ITERATIONS.
# Example: DP_TASKS=binary,continuous Rscript 00_run_all.R
source("replication_helpers.R", local = TRUE)
analyses <- c(binary = "binary_simulations.R", continuous = "continuous_simulations.R",
              villages = "karnataka_villages_analysis.R", twitch = "twitch_analysis.R",
              amazon = "amazon_analysis.R")
requested <- strsplit(Sys.getenv("DP_TASKS", paste(names(analyses), collapse = ",")), ",", fixed = TRUE)[[1]]
if (!length(requested) || any(!requested %in% names(analyses))) {
  stop("DP_TASKS must be a comma-separated subset of binary,continuous,villages,twitch,amazon.")
}
require_packages(c("data.table", "ggplot2", "igraph", "Matrix", "haven", "dplyr", "scales", "ggrepel"))
# Detect missing empirical inputs before spending hours on the simulation grids.
if ("villages" %in% requested) check_data_files(village_input_files(config$data))
if ("twitch" %in% requested) check_data_files(file.path(config$data,
  c("large_twitch_features.csv", "large_twitch_edges.csv")))
if ("amazon" %in% requested) check_data_files(file.path(config$data,
  c("amazon-meta.txt", "com-amazon.ungraph.txt")))
for (name in requested) {
  # A new environment prevents similarly named analysis objects from leaking.
  sys.source(analyses[[name]], envir = new.env(parent = environment()))
}
cat("All requested analyses completed.\n")
