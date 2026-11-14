# Twitch language subnetworks. All plots and private runs use the same within-
# language ranks, induced graphs, and zero-outcome convention for isolates.
# These ranks, languages, and the vertex universe are fixed preprocessing inputs.
source("replication_helpers.R", local = TRUE)
require_packages(c("Matrix", "data.table", "ggplot2", "ggrepel", "scales"))
library(Matrix)
library(data.table)
library(ggplot2)
# Language names used in plots generated from this run.
twitch_language_names <- c(EN="English", FR="French", KO="Korean", JA="Japanese", RU="Russian",
  PL="Polish", DE="German", ES="Spanish", IT="Italian", PT="Portuguese", TR="Turkish",
  ZH="Chinese", SV="Swedish", NL="Dutch", TH="Thai", CS="Czech", DA="Danish", HU="Hungarian",
  FI="Finnish", NO="Norwegian", OTHER="Other")

plot_twitch_english <- function(bins, alpha, beta, correlation) {
  p <- ggplot(bins, aes(x, y)) + geom_point(color = "dodgerblue4", size = 2) +
    geom_abline(intercept = alpha, slope = beta, color = "firebrick", linetype = "dashed") +
    theme_classic() + scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Own View Rank (Percentile)", y = "Average Friend View Rank",
         subtitle = sprintf("English Sub-network (Correlation = %.2f)", correlation))
  save_plot(p, "Twitch_EN_Binscatter", png = TRUE)
}

plot_twitch_privacy <- function(draws, correlations) {
  p <- ggplot(correlations[is.finite(Correlation)], aes(Correlation, fill = factor(Epsilon))) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5,
                   position = "identity", color = "white") +
    scale_fill_manual(values = c("4" = "red", "6" = "green3", "8" = "dodgerblue4")) +
    theme_classic() + labs(x = "Correlation Coefficient (True vs Private Slope across Languages)",
                           y = "Frequency Density", fill = "Epsilon") +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 11), legend.position = "top")
  save_plot(p, "Twitch_Correlation_Histogram_Overlaid_Eps4_6_8", png = TRUE)
  scatter <- copy(draws[epsilon == 8 & iteration == 1 & status == "ok" &
                          is.finite(true_beta) & is.finite(beta)])
  scatter[, language_name := unname(twitch_language_names[language])]
  # Separate axis scales keep the plot legible when private errors are much larger
  # than the variation in true slopes. The dashed reference is y=x, not a 45-degree
  # line on the page; neither axis clips the private observations.
  p <- ggplot(scatter, aes(true_beta, beta)) +
    geom_point(aes(size = n_nodes), color = "dodgerblue4", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(aes(label = language_name), size = 3, seed = 20261107,
                            max.overlaps = Inf, max.time = Inf, max.iter = 10000) +
    theme_classic() + scale_size_continuous(range = c(1,5), labels = scales::comma) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    theme(legend.position = "bottom", plot.subtitle = element_text(hjust = 0.5)) +
    labs(x = "True Slope Coefficient", y = "Private Slope Estimate (Epsilon = 8)", size = "Users",
         subtitle = sprintf("Correlation Coefficient = %.2f",
                            complete_correlation(scatter$true_beta, scatter$beta)))
  save_plot(p, "Twitch_Scatter_Eps8_Run1", png = TRUE)
}

start_experiment("twitch", 2026110740L)

meta_file <- file.path(config$data, "large_twitch_features.csv")
graph_file <- file.path(config$data, "large_twitch_edges.csv")
check_data_files(c(meta_file, graph_file))
features <- fread(meta_file, select = c("numeric_id", "language", "views"))
edges <- fread(graph_file)
stopifnot(!anyDuplicated(features$numeric_id), all(is.finite(features$views)), !anyNA(features$language))
features[, x_raw := rank(views, ties.method = "average") / .N, by = language]
setorder(features, numeric_id)
from_index <- match(edges$numeric_id_1, features$numeric_id)
to_index <- match(edges$numeric_id_2, features$numeric_id)
if (anyNA(c(from_index, to_index))) stop("Twitch edges have endpoints missing from metadata.")
edge_language <- features$language[from_index]
within_language <- edge_language == features$language[to_index]
edges <- edges[within_language]
edges[, language := edge_language[within_language]]
rm(from_index, to_index, edge_language, within_language)

records <- list()
language_ids <- sort(unique(features$language))
# Assign distinct seeds using the sorted language set, independently of execution order.
language_seeds <- setNames(2026110740L + 100L * seq_along(language_ids), language_ids)
for (language_id in language_ids) {
  local_nodes <- features[language == language_id]
  local_edges <- edges[language == language_id]
  weights <- weights_from_edges(local_nodes$numeric_id, local_edges$numeric_id_1, local_edges$numeric_id_2)
  x <- local_nodes$x_raw
  y <- as.numeric(weights %*% x)
  truth <- fit_line(x, y)
  cat("Twitch language:", language_id, "nodes:", length(x), "isolates:", sum(rowSums(weights) == 0), "\n")

  if (language_id == "EN") {
    en_data <- data.table(x = x, y = y)
    en_data[, bin := pmin(20L, floor(x * 20) + 1L)]
    bins <- en_data[, .(x = mean(x), y = mean(y), n = .N), by = bin]
    plot_twitch_english(bins, truth["alpha"], truth["beta"], cor(x, y))
  }

  for (epsilon in c(4, 6, 8)) {
    set.seed(language_seeds[[language_id]] + epsilon)
    draws <- simulate_continuous(weights, x, epsilon / 2, epsilon / 2, 1e-3, iteration_count(1000L))
    draws[, `:=`(language = language_id, true_beta = truth["beta"], n_nodes = length(x))]
    records[[paste(language_id, epsilon)]] <- draws
  }
}
draws <- rbindlist(records)
summary <- draws[, accuracy_summary(beta, true_beta, status),
                 by = .(language, n_nodes, true_beta, epsilon, epsilon_label, epsilon_edge, delta_label)]
fwrite(summary, output_path("DerivedData", "twitch_language_summary.csv"))
# Correlations are conditional on successful, finite estimates within each draw.
# The retained language set can vary across draws; record its size and exclusions.
correlations <- draws[, {
  valid <- status == "ok" & is.finite(true_beta) & is.finite(beta)
  .(Correlation = complete_correlation(true_beta[valid], beta[valid]),
    n_languages = .N, n_languages_used = sum(valid), n_failed = sum(!valid))
}, by = .(epsilon, epsilon_label, epsilon_edge, delta_label, iteration)]
setnames(correlations, "epsilon", "Epsilon")
fwrite(correlations, output_path("DerivedData", "twitch_dp_results.csv"))
print(correlations[, .(mean_correlation = mean(Correlation, na.rm = TRUE),
                       n_valid = sum(is.finite(Correlation)), n_attempted = .N,
                       mean_languages_used = mean(n_languages_used),
                       min_languages_used = min(n_languages_used)), by = Epsilon])
plot_twitch_privacy(draws, correlations)
cat("Twitch analysis complete.\n")
