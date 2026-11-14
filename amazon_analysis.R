# Amazon co-purchase application. Ranks, retained product set, and metadata
# availability are treated as fixed inputs to the privacy experiment.
source("replication_helpers.R", local = TRUE)
require_packages(c("Matrix", "data.table", "ggplot2"))
library(Matrix)
library(data.table)
library(ggplot2)
start_experiment("amazon", 2026110730L)

meta_file <- file.path(config$data, "amazon-meta.txt")
graph_file <- file.path(config$data, "com-amazon.ungraph.txt")
check_data_files(c(meta_file, graph_file))
cat("Reading Amazon metadata...\n")
nodes <- read_amazon_metadata(meta_file)
# Compute ranks over all metadata products with positive sales rank.
nodes[, x_raw := 1 - rank(salesrank, ties.method = "average") / .N]
edges <- fread(graph_file, skip = "# FromNodeId", header = TRUE)
setnames(edges, names(edges), c("from", "to"))
edges <- edges[from %in% nodes$id & to %in% nodes$id]
# Freeze this preprocessed vertex set throughout the privacy experiment.
# It is not selected again after a neighboring edge is changed.
vertices <- sort(unique(c(edges$from, edges$to)))
nodes <- nodes[match(vertices, id)]
weights <- weights_from_edges(vertices, edges$from, edges$to)
x <- nodes$x_raw
y <- as.numeric(weights %*% x)
truth <- fit_line(x, y)
cat("Amazon nodes:", length(x), "edges:", nrow(edges), "true slope:", truth["beta"], "\n")

# Plot the product-level OLS line. Fitting bin means would change the estimand.
plot_data <- data.table(x = x, y = y)
plot_data[, bin := cut(x, breaks = 50, include.lowest = TRUE, labels = FALSE)]
bins <- plot_data[, .(x = mean(x), y = mean(y), n = .N), by = bin]
p <- ggplot(bins, aes(x, y)) + geom_point(color = "dodgerblue4") +
  geom_abline(intercept = truth["alpha"], slope = truth["beta"], color = "firebrick", linetype = "dashed") +
  theme_classic() + labs(x = "Product sales percentile", y = "Average co-purchased product percentile")
save_plot(p, "amazon_salesrank_regression", png = TRUE)

records <- list()
for (epsilon in c(4, 6, 8)) {
  set.seed(2026110730L + epsilon)
  cat("Amazon total epsilon:", epsilon, "\n")
  draws <- simulate_continuous(weights, x, epsilon / 2, epsilon / 2, 1e-3, iteration_count(1000L))
  draws[, true_slope := truth["beta"]]
  records[[as.character(epsilon)]] <- draws
}
results <- rbindlist(records)
setnames(results, "beta", "slope")
fwrite(results, output_path("DerivedData", "dp_simulation_amazon_multi_eps.csv"))
summary <- results[, accuracy_summary(slope, true_slope, status),
                   by = .(epsilon, epsilon_label, epsilon_edge, delta_label)]
print(summary)
p <- ggplot(results[status == "ok"], aes(slope, y = after_stat(density),
                                          fill = factor(epsilon), color = factor(epsilon))) +
  geom_histogram(bins = 50, alpha = 0.4, position = "identity") +
  geom_vline(xintercept = truth["beta"], color = "black", linetype = "dashed") +
  theme_classic() + labs(x = "Private slope estimate", y = "Density among successful estimates",
                         fill = "Total epsilon", color = "Total epsilon")
save_plot(p, "dp_coefficient_histogram_overlaid", png = TRUE)
cat("Amazon analysis complete.\n")
