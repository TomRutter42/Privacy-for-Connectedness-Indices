# Continuous-label graphon experiments for manuscript Algorithm 2.
# Conditional covariance/mean noise; no variance perturbation.
source("replication_helpers.R", local = TRUE)
require_packages(c("igraph", "Matrix", "data.table", "ggplot2", "scales"))
library(Matrix)
library(data.table)
library(ggplot2)
start_experiment("continuous", 2026110710L)

# 1. Allocation at fixed total epsilon 8 and delta 0.001.
# All noise replicates for each N condition on one generated graph and label vector.
records <- list()
for (N in c(20000L, 50000L, 100000L, 500000L)) {
  set.seed(2026110710L + N)
  network <- generate_graphon_network(N, 20, 0.8)
  weights <- graph_weights(network$g)
  truth <- fit_line(network$x_raw, as.numeric(weights %*% network$x_raw))
  for (edge_fraction in seq(0.1, 0.8, by = 0.1)) {
    cat("Continuous allocation: N =", N, "edge fraction =", edge_fraction, "\n")
    draws <- simulate_continuous(weights, network$x_raw, 8 * (1 - edge_fraction),
                                 8 * edge_fraction, 1e-3, iteration_count(3000L))
    draws[, `:=`(N = N, frac_edge_dp = edge_fraction, true_slope = truth["beta"])]
    records[[length(records) + 1L]] <- draws
  }
}
draws <- rbindlist(records)
results <- draws[, accuracy_summary(beta, true_slope, status),
                 by = .(N, frac_edge_dp, epsilon, epsilon_label, epsilon_edge, delta_label, true_slope)]
results[, frac_local_dp := 1 - frac_edge_dp]
p <- ggplot(results[is.finite(mse) & mse > 0], aes(frac_local_dp, mse, color = factor(N), group = N)) +
  geom_line() + geom_point() + scale_y_log10() + theme_classic() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0.2, 0.9, by = 0.1)) +
  scale_color_discrete(labels = function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)) +
  labs(x = expression(Ratio ~ epsilon[l] / (epsilon[l] + epsilon[e])),
       y = "Mean Squared Error (Log Scale)", color = "Number of Nodes")
save_plot(p, "graphon_allocation_mse_plot", png = TRUE)

# 2 and 3. The same draws estimate slope and bottom-quartile mean friend rank.
# Mean friend rank on [0,1/4] is alpha + beta/8.
records <- list()
for (N in c(10000L, 20000L, 50000L, 100000L)) {
  for (homophily in seq(0, 1, by = 0.2)) {
    set.seed(2026110711L + N + as.integer(round(homophily * 100)))
    cat("Continuous homophily: N =", N, "h =", homophily, "\n")
    network <- generate_graphon_network(N, 20, homophily)
    weights <- graph_weights(network$g)
    truth <- fit_line(network$x_raw, as.numeric(weights %*% network$x_raw))
    draws <- simulate_continuous(weights, network$x_raw, 2, 2, 1e-3, iteration_count(2500L))
    draws[, `:=`(N = N, homophily = homophily, true_slope = truth["beta"],
                 true_alpha = truth["alpha"], true_r_squared = truth["r_squared"],
                 true_mafr = mean_friend_rank(truth["alpha"], truth["beta"]))]
    records[[length(records) + 1L]] <- draws
  }
}
draws <- rbindlist(records)
groups <- c("N", "homophily", "epsilon", "epsilon_label", "epsilon_edge", "delta_label",
            "true_slope", "true_alpha", "true_r_squared", "true_mafr")
regression <- draws[, accuracy_summary(beta, true_slope, status), by = groups]
mafr <- draws[, accuracy_summary(mafr, true_mafr, status), by = groups]
fwrite(regression, output_path("DerivedData", "graphon_regression_mse_results.csv"))
for (quantity in c("slope", "mafr")) {
  result <- if (quantity == "slope") regression else mafr
  p <- ggplot(result[is.finite(mse) & mse > 0], aes(homophily, mse, color = factor(N), group = N)) +
    geom_line() + geom_point() + scale_y_log10() + theme_classic() +
    scale_color_discrete(labels = function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(x = "Homophily Parameter",
         y = if (quantity == "slope") "Mean Squared Error" else "Mean Squared Error (Log Scale)",
         color = "Number of Nodes")
  # graphon_integral_mse_plot reports the MSE of bottom-quartile mean friend rank.
  save_plot(p, if (quantity == "slope") "graphon_regression_mse_plot" else "graphon_integral_mse_plot", png = TRUE)
}
p <- ggplot(regression[N == max(N)], aes(homophily, true_r_squared)) +
  geom_line(color = "dodgerblue4") + geom_point(color = "dodgerblue4") + theme_classic() +
  scale_y_continuous(labels = scales::percent) + labs(x = "Homophily Parameter", y = expression(R^2))
save_plot(p, "graphon_homophily_rsquared_plot", png = TRUE)
cat("Continuous experiments complete.\n")
