# Binary-label experiments for Algorithm 1.
# All five experiments estimate cross-type connectedness.
source("replication_helpers.R", local = TRUE)
require_packages(c("igraph", "Matrix", "data.table", "ggplot2", "scales"))
library(Matrix)
library(data.table)
library(ggplot2)
start_experiment("binary", 2026110701L)

# Each graph is reused for several privacy draws; budgets share random numbers.
# Save the graph ID so uncertainty can be calculated across independent graphs.
simulate_binary_graph <- function(weights, labels, budgets, repetitions, graph_id) {
  truth <- true_connectedness(weights, labels)
  records <- vector("list", repetitions * nrow(budgets))
  k <- 0L
  for (iteration in seq_len(repetitions)) {
    flips <- runif(length(labels))
    laplace <- r_laplace(1)
    for (j in seq_len(nrow(budgets))) {
      result <- private_connectedness(weights, labels, budgets$epsilon_label[j],
                                       budgets$epsilon_edge[j], flips, laplace, details = TRUE)
      k <- k + 1L
      records[[k]] <- cbind(budgets[j], data.table(
        graph_id = graph_id, iteration = iteration, true_ec = truth,
        estimate = result$estimate, squared_error = (result$estimate - truth)^2,
        stabilized = result$stabilized))
    }
  }
  rbindlist(records)
}

equal_budgets <- function(epsilon) data.table(
  epsilon = epsilon, epsilon_label = epsilon / 2, epsilon_edge = epsilon / 2,
  delta_label = 0)

summarize_binary <- function(draws, groups) {
  summary <- draws[, .(mse = mean(squared_error), n_attempted = .N,
                       n_valid = sum(is.finite(estimate)),
                       n_stabilized = sum(stabilized)), by = groups]
  # Graph-level means, not individual draws, are independent sampling units.
  graph_means <- draws[, .(graph_mse = mean(squared_error)), by = c(groups, "graph_id")]
  uncertainty <- graph_means[, .(mse_se = if (.N > 1) sd(graph_mse) / sqrt(.N) else NA_real_,
                                 n_graphs = .N), by = groups]
  result <- merge(summary, uncertainty, by = groups, sort = FALSE)
  result[, `:=`(ci_lower = pmax(0, mse - 1.96 * mse_se), ci_upper = mse + 1.96 * mse_se)]
  result
}

# 1. Homophily at approximately degree 80; N=2000, as implied by probabilities.
set.seed(2026110701L)
N <- 2000L
labels <- rep(c(0, 1), each = N / 2)
graph_count <- iteration_count(75L)
noise_count <- 15L
records <- list()
for (within in seq(0.04, 0.08, by = 0.005)) {
  cat("Binary homophily:", within, "\n")
  for (graph_id in seq_len(graph_count)) {
    weights <- generate_sbm_network(N / 2, N / 2, within, 0.08 - within, within)
    draws <- simulate_binary_graph(weights, labels, equal_budgets(c(1, 2, 4)), noise_count, graph_id)
    draws[, p_within := within]
    records[[length(records) + 1L]] <- draws
  }
}
draws <- rbindlist(records)
results <- summarize_binary(draws, c("p_within", "epsilon", "epsilon_label", "epsilon_edge"))
fwrite(results, output_path("DerivedData", "MSE_vs_Homophily.csv"))
p <- ggplot(results, aes(p_within, mse, color = factor(epsilon), group = epsilon)) +
  geom_line() + geom_point() + scale_y_log10() + theme_classic() +
  labs(x = "Homophily (Within-Group Probability)", y = "Mean Squared Error (Log Scale)", color = "Epsilon")
save_plot(p, "MSE_vs_Homophily")

# 2. Privacy budget sweep, with uncertainty over graph-level mean errors.
set.seed(2026110702L)
scenarios <- c("No Homophily" = 0.04, "Low Homophily" = 0.06, "High Homophily" = 0.08)
records <- list()
for (scenario in names(scenarios)) {
  cat("Binary privacy sweep:", scenario, "\n")
  within <- scenarios[[scenario]]
  for (graph_id in seq_len(iteration_count(50L))) {
    weights <- generate_sbm_network(N / 2, N / 2, within, 0.08 - within, within)
    draws <- simulate_binary_graph(weights, labels, equal_budgets(1:8), 10L, graph_id)
    draws[, Homophily := scenario]
    records[[length(records) + 1L]] <- draws
  }
}
draws <- rbindlist(records)
results <- summarize_binary(draws, c("Homophily", "epsilon", "epsilon_label", "epsilon_edge"))
fwrite(results, output_path("DerivedData", "accuracy_results_Amin_VarianceReduced.csv"))
p <- ggplot(results, aes(epsilon, mse, color = Homophily, fill = Homophily)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line() + scale_x_continuous(breaks = 1:8) + theme_classic() +
  labs(x = "Epsilon (Privacy Parameter)", y = "Mean Squared Error", color = "Network Type", fill = "Network Type") +
  guides(fill = "none") +
  theme(legend.position = "bottom", legend.key.width = grid::unit(0.8, "lines"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9),
        legend.spacing.x = grid::unit(0.15, "lines"))
save_plot(p, "accuracyVsEpsilon_Amin_VarianceReduced")

# 3. Cell composition; construct the second block as N - n_low to avoid rounding drift.
set.seed(2026110703L)
records <- list()
for (proportion in seq(0.1, 0.9, by = 0.1)) {
  cat("Binary composition:", proportion, "\n")
  n_low <- round(N * (1 - proportion))
  labels <- c(rep(0, n_low), rep(1, N - n_low))
  for (graph_id in seq_len(graph_count)) {
    weights <- generate_sbm_network(n_low, N - n_low, 0.04, 0.04, 0.04)
    draws <- simulate_binary_graph(weights, labels, equal_budgets(c(1, 2, 4)), noise_count, graph_id)
    draws[, prop_high_ses := proportion]
    records[[length(records) + 1L]] <- draws
  }
}
draws <- rbindlist(records)
results <- summarize_binary(draws, c("prop_high_ses", "epsilon", "epsilon_label", "epsilon_edge"))
fwrite(results, output_path("DerivedData", "MSE_vs_CellComposition.csv"))
p <- ggplot(results, aes(prop_high_ses, mse, color = factor(epsilon), group = epsilon)) +
  geom_line() + geom_point() + scale_y_log10() + theme_classic() +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(x = "Proportion High-SES", y = "Mean Squared Error (Log Scale)", color = "Epsilon")
save_plot(p, "MSE_vs_CellComposition")

# 4. Budget allocation for connectedness, conditional on one graph per size.
records <- list()
budgets <- CJ(epsilon = c(0.5, 1, 2, 4), epsilon_label_proportion = seq(0.1, 0.9, by = 0.1))
budgets[, `:=`(epsilon_label = epsilon * epsilon_label_proportion,
               epsilon_edge = epsilon * (1 - epsilon_label_proportion), delta_label = 0)]
for (N in c(2000L, 20000L)) {
  set.seed(2026110704L + N)
  cat("Binary allocation: N =", N, "\n")
  scale <- min(2000 / N, 1) # Maintain expected degree near 80 across network sizes.
  weights <- generate_sbm_network(N / 2, N / 2, 0.06 * scale, 0.02 * scale, 0.06 * scale)
  draws <- simulate_binary_graph(weights, rep(c(0, 1), each = N / 2), budgets,
                                  iteration_count(1000L), 1L)
  draws[, n_nodes := N]
  records[[as.character(N)]] <- draws
}
draws <- rbindlist(records)
results <- summarize_binary(draws, c("n_nodes", "epsilon", "epsilon_label_proportion",
                                    "epsilon_label", "epsilon_edge"))
fwrite(results, output_path("DerivedData", "sbm_dp_mse_results_combined.csv"))
for (N in unique(results$n_nodes)) {
  p <- ggplot(results[n_nodes == N], aes(epsilon_label_proportion, mse,
                                        color = factor(epsilon), group = epsilon)) +
    geom_line() + geom_point() + scale_y_log10() + theme_classic() +
    labs(x = expression(Ratio ~ epsilon[l] / (epsilon[e] + epsilon[l])), y = "Mean Squared Error",
         color = expression(epsilon[e] + epsilon[l]))
  save_plot(p, paste0("sbm_dp_mse_plot_", N))
}

# 5. Network size at fixed expected degree. Every replicate draws a new graph.
records <- list()
sizes <- c(500L, 1000L, 5000L, 10000L)
for (N in sizes) {
  set.seed(2026110705L + N)
  cat("Binary size sweep: N =", N, "\n")
  count <- iteration_count(10000L)
  local_records <- vector("list", count)
  for (graph_id in seq_len(count)) {
    graph <- igraph::sample_gnp(N, min(20 / (N - 1), 1), directed = FALSE, loops = FALSE)
    labels <- sample(rep(c(0, 1), each = N / 2))
    local_records[[graph_id]] <- simulate_binary_graph(graph_weights(graph), labels,
                                                       equal_budgets(c(0.5, 1, 2, 4, 8)), 1L, graph_id)
    if (graph_id %% 1000 == 0) cat("  graph", graph_id, "of", count, "\n")
  }
  draws <- rbindlist(local_records)
  draws[, n_nodes := N]
  records[[as.character(N)]] <- draws
}
draws <- rbindlist(records)
fwrite(draws, output_path("DerivedData", "mse_simulation_results_ec.csv"))
results <- summarize_binary(draws, c("n_nodes", "epsilon", "epsilon_label", "epsilon_edge"))
p <- ggplot(results, aes(n_nodes, mse, color = factor(epsilon), group = epsilon)) +
  geom_line() + geom_point() + scale_x_log10(breaks = sizes) + scale_y_log10() +
  theme_classic() + labs(x = "Number of Nodes (Log Scale)", y = "Mean Squared Error (Log Scale)", color = "Epsilon")
save_plot(p, "mse_vs_size_constant_degree_ec")
cat("Binary experiments complete.\n")
