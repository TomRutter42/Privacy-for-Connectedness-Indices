# Amazon Co-purchase Analysis
# Regresses average sales rank of neighbors on own sales rank, then runs
# DP simulation of the slope estimator.

library(data.table)
library(dplyr)
library(ggplot2)
library(VGAM)
library(stringr)

# --- Configuration ---
base_dir  <- ext_data_dir
meta_file <- file.path(base_dir, "amazon-meta.txt")
graph_file <- file.path(base_dir, "com-amazon.ungraph.txt")

# Simulation parameters
NUM_ITERATIONS  <- 1000
epsilons_to_test <- c(4, 6, 8)
delta_ell       <- 10^(-3)
bounds <- list(a = 0, b = 1, a_prime = 0, b_prime = 1)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading Amazon metadata from:", meta_file, "\n")
if (!file.exists(meta_file)) stop(paste("Meta file not found:", meta_file))
lines <- readLines(meta_file)

id_lines       <- grep("^Id:\\s+\\d+", lines)
salesrank_lines <- grep("^\\s+salesrank:\\s+\\d+", lines)
cat("Found", length(id_lines), "products.\n")

raw_ids  <- lines[id_lines]
ids      <- as.integer(str_extract(raw_ids, "\\d+"))
id_df    <- data.table(line_idx = id_lines, id = ids, type = "id")

rank_vals <- as.integer(str_extract(lines[salesrank_lines], "\\d+"))
rank_df   <- data.table(line_idx = salesrank_lines, salesrank = rank_vals, type = "rank")

combined <- rbind(
  id_df[, .(line_idx, val = id, type = "id")],
  rank_df[, .(line_idx, val = salesrank, type = "rank")]
)
setkey(combined, line_idx)
combined[, current_id := NA_integer_]
combined[type == "id", current_id := val]
combined[, current_id := nafill(current_id, type = "locf")]

node_data <- combined[type == "rank", .(id = current_id, salesrank = val)]
node_data <- node_data[salesrank > 0]
node_data[, x_raw := 1 - (rank(salesrank, ties.method = "average") / .N)]
cat("Extracted sales rank for", nrow(node_data), "products.\n")

cat("Loading graph from:", graph_file, "\n")
if (!file.exists(graph_file)) stop(paste("Graph file not found:", graph_file))
edges     <- fread(graph_file, skip = 4, header = FALSE, col.names = c("from", "to"))
edges_rev <- edges[, .(from = to, to = from)]
all_edges <- unique(rbind(edges, edges_rev))
cat("Symmetrized edge count:", nrow(all_edges), "\n")

# Filter to nodes with attributes and map IDs to 1..N
valid_ids   <- unique(node_data$id)
all_edges   <- all_edges[from %in% valid_ids & to %in% valid_ids]
unique_nodes <- unique(c(all_edges$from, all_edges$to))
node_map    <- data.table(original_id = unique_nodes, mapped_id = seq_along(unique_nodes))
setkey(node_map, original_id)

setkey(node_data, id)
node_data <- merge(node_data, node_map, by.x = "id", by.y = "original_id")
setkey(node_data, mapped_id)
x_raw    <- node_data$x_raw
n_nodes  <- length(x_raw)

all_edges <- merge(all_edges, node_map, by.x = "from", by.y = "original_id")
setnames(all_edges, "mapped_id", "from_mapped")
all_edges <- merge(all_edges, node_map, by.x = "to", by.y = "original_id")
setnames(all_edges, "mapped_id", "to_mapped")

edge_from <- all_edges$from_mapped
edge_to   <- all_edges$to_mapped

degree_counts <- table(edge_from)
degrees <- numeric(n_nodes)
degrees[as.integer(names(degree_counts))] <- as.integer(degree_counts)

cat("Data loaded: N =", n_nodes, "nodes,", length(edge_from), "edges.\n")

# ==============================================================================
# 2. BINSCATTER: Average neighbor rank vs own rank
# ==============================================================================

cat("\nComputing true average neighbor sales rank...\n")
neighbor_vals_true <- x_raw[edge_to]
node_sums_true     <- rowsum(neighbor_vals_true, edge_from)
y_true             <- numeric(n_nodes)
active_ids_true    <- as.integer(rownames(node_sums_true))
y_true[active_ids_true] <- node_sums_true[, 1] / degrees[active_ids_true]
y_true[is.na(y_true) | is.infinite(y_true)] <- 0

final_data <- data.table(own_rank = x_raw, avg_neighbor_rank = y_true)[y_true > 0]
cat("Regression dataset size:", nrow(final_data), "\n")

model <- lm(avg_neighbor_rank ~ own_rank, data = final_data)
true_slope <- coef(model)[2]
cat(sprintf("True slope: %.5f\n", true_slope))
print(summary(model))

final_data[, bin_cut := cut(own_rank, breaks = 50, include.lowest = TRUE, labels = FALSE)]
binned_data <- final_data[, .(
  mean_own_rank      = mean(own_rank),
  mean_neighbor_rank = mean(avg_neighbor_rank)
), by = bin_cut]

p_binscatter <- ggplot(binned_data, aes(x = mean_own_rank, y = mean_neighbor_rank)) +
  geom_point(size = 2, color = "dodgerblue4") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +
  labs(
    x = "Product Sales Rank (Percentile)",
    y = "Average Sales Rank of Co-Purchased Products (Percentile)"
  ) +
  theme_classic()

ggsave("Figures/amazon_salesrank_regression.pdf", p_binscatter, width = 5, height = 5)
ggsave("Figures/amazon_salesrank_regression.png", p_binscatter, width = 5, height = 5)
cat("Binscatter saved to Figures/amazon_salesrank_regression.pdf and .png\n")

# ==============================================================================
# 3. HELPER FUNCTIONS FOR DP SIMULATION
# ==============================================================================

get_truncated_laplace_variance <- function(lambda, A) {
  term1     <- 2 * lambda^2
  numerator <- exp(-A / lambda) * (A^2 + 2 * lambda * A)
  denominator <- 1 - exp(-A / lambda)
  term1 - (numerator / denominator)
}

r_truncated_laplace <- function(n, lambda, A) {
  samples <- numeric(n)
  count   <- 0
  while (count < n) {
    raw   <- rlaplace((n - count) * 2, location = 0, scale = lambda)
    valid <- raw[abs(raw) <= A]
    take  <- min(length(valid), n - count)
    if (take > 0) {
      samples[(count + 1):(count + take)] <- valid[1:take]
      count <- count + take
    }
  }
  samples
}

algorithm_3_noisystats <- function(x, y, epsilon, bounds) {
  n <- length(x)
  if (n < 2) return(list(beta = NA, alpha = NA))

  a <- bounds$a; b <- bounds$b
  a_prime <- bounds$a_prime; b_prime <- bounds$b_prime

  bar_x  <- mean(x); bar_y <- mean(y)
  nvar_x <- sum((x - bar_x)^2)
  ncov_xy <- sum((x - bar_x) * (y - bar_y))

  Delta_1 <- (1 - 1 / n) * (b - a)^2
  Delta_2 <- (1 - 1 / n) * (b - a) * (b_prime - a_prime)

  L_1 <- rlaplace(1, 0, 3 * Delta_1 / epsilon)
  L_2 <- rlaplace(1, 0, 3 * Delta_2 / epsilon)

  denom <- nvar_x + L_2
  num   <- ncov_xy + L_1

  if (denom > 1e-10) {
    beta_tilde <- num / denom
    Delta_3    <- (b_prime - a_prime) / n + abs(beta_tilde) * (b - a) / n
    L_3        <- rlaplace(1, 0, 3 * Delta_3 / epsilon)
    alpha_tilde <- (bar_y - beta_tilde * bar_x) + L_3
    list(beta = beta_tilde, alpha = alpha_tilde)
  } else {
    list(beta = NA, alpha = NA)
  }
}

# ==============================================================================
# 4. DP SIMULATION
# ==============================================================================

cat("\nStarting DP simulation (", NUM_ITERATIONS, "runs, epsilons:",
    paste(epsilons_to_test, collapse = ", "), ")...\n")

all_results <- list()

for (eps_val in epsilons_to_test) {
  epsilon_ell <- eps_val
  epsilon_e   <- eps_val

  Delta    <- 1.0
  lambda   <- Delta / epsilon_ell
  A        <- lambda * log(1 + (exp(epsilon_ell) - 1) / (2 * delta_ell))
  sigma_sq <- get_truncated_laplace_variance(lambda, A)

  cat("\n-------------------------------------------------\n")
  cat("Running for epsilon =", eps_val, "\n")
  cat("  lambda:", lambda, " A:", A, " Sigma^2:", sigma_sq, "\n")

  private_slopes <- numeric(NUM_ITERATIONS)

  for (iter in 1:NUM_ITERATIONS) {
    if (iter %% 100 == 0) cat("  Iteration", iter, "...\n")

    z     <- r_truncated_laplace(n_nodes, lambda, A)
    x_hat <- x_raw + z

    neighbor_vals  <- x_hat[edge_to]
    node_sums_mat  <- rowsum(neighbor_vals, edge_from)
    y_hat          <- numeric(n_nodes)
    active_ids     <- as.integer(rownames(node_sums_mat))
    y_hat[active_ids] <- node_sums_mat[, 1] / degrees[active_ids]
    y_hat[is.na(y_hat) | is.infinite(y_hat)] <- 0

    res      <- algorithm_3_noisystats(x_hat, y_hat, epsilon_e, bounds)
    hat_beta <- res$beta

    if (!is.na(hat_beta)) {
      s2_x <- var(x_hat)
      if ((s2_x - sigma_sq) > 1e-6) {
        private_slopes[iter] <- hat_beta * s2_x / (s2_x - sigma_sq)
      } else {
        private_slopes[iter] <- hat_beta
      }
    } else {
      private_slopes[iter] <- NA
    }
  }

  valid_slopes <- private_slopes[!is.na(private_slopes)]
  all_results[[as.character(eps_val)]] <- data.table(
    iteration = seq_along(valid_slopes),
    slope     = valid_slopes,
    epsilon   = as.character(eps_val)
  )
}

results <- rbindlist(all_results)
results[, epsilon := factor(epsilon, levels = as.character(epsilons_to_test))]

fwrite(results, "DerivedData/dp_simulation_amazon_multi_eps.csv")

p_hist <- ggplot(results, aes(x = slope, y = after_stat(density), fill = epsilon, color = epsilon)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  geom_vline(xintercept = true_slope, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    x     = "Private Slope Estimate",
    y     = "Density",
    fill  = "Epsilon",
    color = "Epsilon"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("Figures/dp_coefficient_histogram_overlaid.pdf", p_hist, width = 5, height = 5)
ggsave("Figures/dp_coefficient_histogram_overlaid.png", p_hist, width = 5, height = 5)
cat("Histogram saved to Figures/dp_coefficient_histogram_overlaid.pdf and .png\n")
