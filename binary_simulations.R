# Binary Label Simulations
# Five simulations studying DP estimation of the exposure coefficient (EC)
# under SBM and Erdos-Renyi networks with binary node labels.
#
# 1. MSE vs Homophily         — sweeps within-group edge probability
# 2. MSE vs Epsilon by Scenario — three homophily scenarios across epsilon 1-8
# 3. MSE vs Cell Composition  — sweeps proportion of group-1 nodes
# 4. Optimal epsilon split     — sweeps label/edge budget split in SBM
# 5. MSE vs Network Size       — sweeps N under constant-degree Erdos-Renyi

library(VGAM)
library(ggplot2)
library(igraph)
library(Matrix)
library(data.table)
library(dplyr)
library(scales)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

generate_symmetric_matrix <- function(n, prob) {
  if (n <= 0) {
    return(matrix(nrow = 0, ncol = 0))
  }
  mat <- matrix(
    sample(c(0, 1), size = n^2, prob = c(1 - prob, prob), replace = TRUE),
    ncol = n, nrow = n
  )
  mat[lower.tri(mat)] <- t(mat[upper.tri(mat)])
  diag(mat) <- 0
  return(mat)
}

# Returns row-normalised adjacency matrix for a two-group SBM.
# Used by simulations 1 and 2 (fixed 50/50 composition).
generate_sbm_network <- function(num_low, num_high, p_ll, p_lh, p_hh) {
  ll <- generate_symmetric_matrix(num_low, p_ll)
  hh <- generate_symmetric_matrix(num_high, p_hh)
  lh <- matrix(
    sample(c(0, 1), size = num_low * num_high, prob = c(1 - p_lh, p_lh), replace = TRUE),
    ncol = num_low, nrow = num_high
  )
  A <- cbind(rbind(ll, lh), rbind(t(lh), hh))
  row_sums <- rowSums(A)
  rA <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  has_friends <- row_sums > 0
  rA[has_friends, ] <- A[has_friends, ] / row_sums[has_friends]
  return(rA)
}

# Returns list(rA, ses_vector) for arbitrary composition.
# Used by simulation 3.
generate_composition_network <- function(num_nodes, prop_high, p_ll, p_lh, p_hh) {
  num_low <- round(num_nodes * (1 - prop_high))
  num_high <- round(num_nodes * prop_high)

  if (num_low == 0) {
    A <- generate_symmetric_matrix(num_high, p_hh)
    rs <- rowSums(A)
    rA <- matrix(0, num_high, num_high)
    rA[rs > 0, ] <- A[rs > 0, ] / rs[rs > 0]
    return(list(rA = rA, ses_vector = rep(1, num_high)))
  }
  if (num_high == 0) {
    A <- generate_symmetric_matrix(num_low, p_ll)
    rs <- rowSums(A)
    rA <- matrix(0, num_low, num_low)
    rA[rs > 0, ] <- A[rs > 0, ] / rs[rs > 0]
    return(list(rA = rA, ses_vector = rep(0, num_low)))
  }

  ses_vector <- c(rep(0, num_low), rep(1, num_high))
  rA <- generate_sbm_network(num_low, num_high, p_ll, p_lh, p_hh)
  return(list(rA = rA, ses_vector = ses_vector))
}

# True exposure coefficient: mean fraction of group-1 connections among group-0 nodes.
compute_true_iec <- function(rA, ses_vector) {
  low_mask <- ses_vector == 0
  if (sum(low_mask) == 0) {
    return(NA)
  }
  mean(as.numeric(rA %*% ses_vector)[low_mask], na.rm = TRUE)
}

# DP estimate using coupled randomised response + Laplace noise.
compute_noisy_iec_coupled <- function(rA, ses_vector, epsilon_s, epsilon_e, noise_seed) {
  p <- 1 / (1 + exp(epsilon_s))
  hat_ses <- ifelse(noise_seed$U_flip[seq_along(ses_vector)] < p, 1 - ses_vector, ses_vector)

  rho_hat <- as.numeric(rA %*% hat_ses)
  w <- (as.numeric(hat_ses == 0) - p) / (1 - 2 * p)
  rho_tilde <- (rho_hat - p) / (1 - 2 * p)

  S0 <- sum(w)
  S1 <- sum(w * rho_tilde) + noise_seed$U_laplace * (2 / ((1 - 2 * p) * epsilon_e))

  if (abs(S0) > 1e-10) S1 / S0 else NA
}

# Row-normalise a sparse adjacency matrix (used by simulation 5).
row_normalize_sparse <- function(A) {
  inv_d <- ifelse(rowSums(A) > 0, 1 / rowSums(A), 0)
  Diagonal(x = inv_d) %*% A
}

# ==============================================================================
# SIMULATION 1: MSE vs Homophily (constant average degree)
# ==============================================================================
cat("\n=== Simulation 1: MSE vs Homophily ===\n")

NUM_NODES_1 <- 2000
NUM_LOW_1 <- 1000
NUM_HIGH_1 <- 1000
ses_vector_1 <- c(rep(0, NUM_LOW_1), rep(1, NUM_HIGH_1))

p_within_values <- seq(0.04, 0.08, by = 0.005)
epsilon_values_1 <- c(1, 2, 4)
num_fixed_graphs_1 <- 75
num_coupled_seeds_1 <- 15

results_list_1 <- list()

for (p_w in p_within_values) {
  p_b <- 0.08 - p_w
  cat(sprintf("  p_within = %.3f, p_between = %.3f\n", p_w, p_b))

  fixed_graphs <- vector("list", num_fixed_graphs_1)
  true_iecs <- numeric(num_fixed_graphs_1)
  for (g in seq_len(num_fixed_graphs_1)) {
    fixed_graphs[[g]] <- generate_sbm_network(NUM_LOW_1, NUM_HIGH_1, p_w, p_b, p_w)
    true_iecs[g] <- compute_true_iec(fixed_graphs[[g]], ses_vector_1)
  }

  total_samples <- num_fixed_graphs_1 * num_coupled_seeds_1
  sq_errors_mat <- matrix(NA, nrow = length(epsilon_values_1), ncol = total_samples)
  sample_idx <- 0

  for (g in seq_len(num_fixed_graphs_1)) {
    rA <- fixed_graphs[[g]]
    true_ec <- true_iecs[g]
    for (s in seq_len(num_coupled_seeds_1)) {
      sample_idx <- sample_idx + 1
      noise_seed <- list(U_flip = runif(NUM_NODES_1), U_laplace = rlaplace(1, 0, 1))
      for (i in seq_along(epsilon_values_1)) {
        eps <- epsilon_values_1[i]
        est <- compute_noisy_iec_coupled(rA, ses_vector_1, eps / 2, eps / 2, noise_seed)
        sq_err <- (est - true_ec)^2
        if (is.finite(sq_err)) sq_errors_mat[i, sample_idx] <- sq_err
      }
    }
  }

  results_list_1[[as.character(p_w)]] <- data.frame(
    p_within = p_w,
    epsilon  = epsilon_values_1,
    mse      = rowMeans(sq_errors_mat, na.rm = TRUE)
  )
}

results_1 <- do.call(rbind, results_list_1)
rownames(results_1) <- NULL
results_1$epsilon_factor <- as.factor(results_1$epsilon)

p1 <- ggplot(results_1, aes(x = p_within, y = mse, color = epsilon_factor)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  labs(
    x     = "Homophily (Within-Group Probability)",
    y     = "Mean Squared Error (Log Scale)",
    color = "Epsilon"
  ) +
  theme_classic() +
  theme(legend.position = "right")

ggsave("Figures/MSE_vs_Homophily.pdf", p1, width = 5, height = 5)
write.csv(results_1, "DerivedData/MSE_vs_Homophily.csv", row.names = FALSE)
cat("Saved: Figures/MSE_vs_Homophily.pdf\n")

# ==============================================================================
# SIMULATION 2: MSE vs Epsilon by Homophily
# ==============================================================================
cat("\n=== Simulation 2: MSE vs Epsilon by Homophily ===\n")

scenarios <- list(
  list(name = "No Homophily", p_ll = 0.04, p_lh = 0.04, p_hh = 0.04),
  list(name = "Low Homophily", p_ll = 0.06, p_lh = 0.02, p_hh = 0.06),
  list(name = "High Homophily", p_ll = 0.08, p_lh = 0, p_hh = 0.08)
)

epsilon_values_2 <- 1:8
num_fixed_graphs_2 <- 50
num_coupled_seeds_2 <- 10

results_list_2 <- list()

for (scenario in scenarios) {
  cat("  Scenario:", scenario$name, "\n")

  fixed_graphs <- vector("list", num_fixed_graphs_2)
  true_iecs <- numeric(num_fixed_graphs_2)
  for (g in seq_len(num_fixed_graphs_2)) {
    fixed_graphs[[g]] <- generate_sbm_network(
      NUM_LOW_1, NUM_HIGH_1,
      scenario$p_ll, scenario$p_lh, scenario$p_hh
    )
    true_iecs[g] <- compute_true_iec(fixed_graphs[[g]], ses_vector_1)
  }

  total_samples <- num_fixed_graphs_2 * num_coupled_seeds_2
  sq_errors_mat <- matrix(NA, nrow = length(epsilon_values_2), ncol = total_samples)
  sample_idx <- 0

  for (g in seq_len(num_fixed_graphs_2)) {
    rA <- fixed_graphs[[g]]
    true_ec <- true_iecs[g]
    for (s in seq_len(num_coupled_seeds_2)) {
      sample_idx <- sample_idx + 1
      noise_seed <- list(U_flip = runif(NUM_NODES_1), U_laplace = rlaplace(1, 0, 1))
      for (i in seq_along(epsilon_values_2)) {
        eps <- epsilon_values_2[i]
        est <- compute_noisy_iec_coupled(rA, ses_vector_1, eps / 2, eps / 2, noise_seed)
        sq_err <- (est - true_ec)^2
        if (is.finite(sq_err)) sq_errors_mat[i, sample_idx] <- sq_err
      }
    }
  }

  mse_vals <- rowMeans(sq_errors_mat, na.rm = TRUE)
  mse_se <- apply(sq_errors_mat, 1, sd, na.rm = TRUE) / sqrt(total_samples)

  results_list_2[[scenario$name]] <- data.frame(
    epsilon = epsilon_values_2,
    mse = mse_vals,
    ci_lower = mse_vals - 1.96 * mse_se,
    ci_upper = mse_vals + 1.96 * mse_se,
    Homophily = scenario$name
  )
}

results_2 <- do.call(rbind, results_list_2)
rownames(results_2) <- NULL

p2 <- ggplot(results_2, aes(x = epsilon, y = mse, color = Homophily, fill = Homophily)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_x_continuous(breaks = 1:8, limits = c(1, 8)) +
  labs(
    x     = "Epsilon (Privacy Parameter)",
    y     = "Mean Squared Error",
    color = "Network Type"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("Figures/accuracyVsEpsilon_Amin_VarianceReduced.pdf", p2, width = 5, height = 5)
write.csv(results_2, "DerivedData/accuracy_results_Amin_VarianceReduced.csv", row.names = FALSE)
cat("Saved: Figures/accuracyVsEpsilon_Amin_VarianceReduced.pdf\n")

# ==============================================================================
# SIMULATION 3: MSE vs Cell Composition
# ==============================================================================
cat("\n=== Simulation 3: MSE vs Cell Composition ===\n")

prop_high_values <- seq(0.1, 0.9, by = 0.1)
epsilon_values_3 <- c(1, 2, 4)
num_fixed_graphs_3 <- 75
num_coupled_seeds_3 <- 15
P_EDGE <- 0.04

results_list_3 <- list()

for (prop_high in prop_high_values) {
  cat(sprintf("  prop_high_ses = %.1f\n", prop_high))

  fixed_graphs <- vector("list", num_fixed_graphs_3)
  true_iecs <- numeric(num_fixed_graphs_3)
  for (g in seq_len(num_fixed_graphs_3)) {
    net <- generate_composition_network(2000, prop_high, P_EDGE, P_EDGE, P_EDGE)
    fixed_graphs[[g]] <- net
    true_iecs[g] <- compute_true_iec(net$rA, net$ses_vector)
  }

  total_samples <- num_fixed_graphs_3 * num_coupled_seeds_3
  sq_errors_mat <- matrix(NA, nrow = length(epsilon_values_3), ncol = total_samples)
  sample_idx <- 0

  for (g in seq_len(num_fixed_graphs_3)) {
    net <- fixed_graphs[[g]]
    true_ec <- true_iecs[g]
    if (is.na(true_ec)) next

    for (s in seq_len(num_coupled_seeds_3)) {
      sample_idx <- sample_idx + 1
      noise_seed <- list(U_flip = runif(length(net$ses_vector)), U_laplace = rlaplace(1, 0, 1))
      for (i in seq_along(epsilon_values_3)) {
        eps <- epsilon_values_3[i]
        est <- compute_noisy_iec_coupled(net$rA, net$ses_vector, eps / 2, eps / 2, noise_seed)
        sq_err <- (est - true_ec)^2
        if (is.finite(sq_err)) sq_errors_mat[i, sample_idx] <- sq_err
      }
    }
  }

  results_list_3[[as.character(prop_high)]] <- data.frame(
    prop_high_ses = prop_high,
    epsilon       = epsilon_values_3,
    mse           = rowMeans(sq_errors_mat, na.rm = TRUE)
  )
}

results_3 <- do.call(rbind, results_list_3)
rownames(results_3) <- NULL
results_3$epsilon_factor <- as.factor(results_3$epsilon)

p3 <- ggplot(results_3, aes(x = prop_high_ses, y = mse, color = epsilon_factor)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(
    x     = "Proportion High-SES",
    y     = "Mean Squared Error (Log Scale)",
    color = "Epsilon"
  ) +
  theme_classic() +
  theme(legend.position = "right")

ggsave("Figures/MSE_vs_CellComposition.pdf", p3, width = 5, height = 5)
write.csv(results_3, "DerivedData/MSE_vs_CellComposition.csv", row.names = FALSE)
cat("Saved: Figures/MSE_vs_CellComposition.pdf\n")

# ==============================================================================
# SIMULATION 4: Optimal epsilon split in SBM (binary regression)
# ==============================================================================
cat("\n=== Simulation 4: Optimal Epsilon Split in SBM ===\n")

N_values_4 <- c(2000, 20000)
NUM_ITER_4 <- 1000
eps_sums_4 <- c(0.5, 1, 2, 4)
props_4 <- seq(0.1, 0.9, by = 0.1)

results_list_4 <- list()

for (N in N_values_4) {
  p_in <- 0.06 * (2000 / N)
  p_out <- 0.02 * (2000 / N)
  prob_matrix <- matrix(c(p_in, p_out, p_out, p_in), nrow = 2)

  set.seed(42)
  cat(sprintf("  Generating SBM for N = %d...\n", N))
  g <- sample_sbm(N, pref.matrix = prob_matrix, block.sizes = c(N / 2, N / 2), directed = FALSE)

  x_raw <- c(rep(1, N / 2), rep(0, N / 2))

  edges <- as_edgelist(g, names = FALSE)
  edge_from <- c(edges[, 1], edges[, 2])
  edge_to <- c(edges[, 2], edges[, 1])
  degrees <- degree(g)

  y_true <- numeric(N)
  node_sums_true <- rowsum(x_raw[edge_to], edge_from)
  active <- as.integer(rownames(node_sums_true))
  y_true[active] <- node_sums_true[, 1] / degrees[active]
  y_true[is.na(y_true) | is.infinite(y_true)] <- 0

  true_slope <- coef(lm(y_true ~ x_raw))[2]
  cat(sprintf("  True slope: %.4f\n", true_slope))

  for (e_sum in eps_sums_4) {
    for (prop in props_4) {
      eps_ell <- e_sum * prop
      eps_e <- e_sum * (1 - prop)
      p_flip <- 1 / (1 + exp(eps_ell))
      sigma_sq <- p_flip * (1 - p_flip) / (1 - 2 * p_flip)^2

      min_val <- -p_flip / (1 - 2 * p_flip)
      max_val <- (1 - p_flip) / (1 - 2 * p_flip)

      private_slopes <- numeric(NUM_ITER_4)

      for (iter in seq_len(NUM_ITER_4)) {
        flip <- runif(N) < p_flip
        x_hat <- ((ifelse(flip, 1 - x_raw, x_raw)) - p_flip) / (1 - 2 * p_flip)

        node_sums_mat <- rowsum(x_hat[edge_to], edge_from)
        y_hat <- numeric(N)
        active_ids <- as.integer(rownames(node_sums_mat))
        y_hat[active_ids] <- ifelse(degrees[active_ids] > 0,
          node_sums_mat[, 1] / degrees[active_ids], 0
        )
        y_hat[is.na(y_hat) | is.infinite(y_hat)] <- 0

        bar_x <- mean(x_hat)
        bar_y <- mean(y_hat)
        nvar_x <- sum((x_hat - bar_x)^2)
        ncov_xy <- sum((x_hat - bar_x) * (y_hat - bar_y))

        Delta <- (1 - 1 / N) * (max_val - min_val)^2
        denom <- nvar_x + rlaplace(1, 0, 2 * Delta / eps_e)
        num <- ncov_xy + rlaplace(1, 0, 2 * Delta / eps_e)

        if (denom > 1e-10) {
          hat_beta <- num / denom
          s2_x <- var(x_hat)
          private_slopes[iter] <- if ((s2_x - sigma_sq) > 1e-6) {
            hat_beta * s2_x / (s2_x - sigma_sq)
          } else {
            hat_beta
          }
        } else {
          private_slopes[iter] <- NA
        }
      }

      valid <- private_slopes[!is.na(private_slopes)]
      results_list_4[[length(results_list_4) + 1]] <- data.table(
        N            = N,
        eps_sum      = e_sum,
        eps_ell_prop = prop,
        eps_ell      = eps_ell,
        eps_e        = eps_e,
        mse          = mean((valid - true_slope)^2),
        n_valid      = length(valid)
      )
    }
  }
}

results_4 <- rbindlist(results_list_4)
fwrite(results_4, "DerivedData/sbm_dp_mse_results_combined.csv")

p4_base <- ggplot(
  results_4[!is.na(mse)],
  aes(
    x = eps_ell_prop, y = mse,
    color = as.factor(eps_sum), group = eps_sum
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  labs(
    x     = expression(Ratio ~ epsilon[l] / (epsilon[e] + epsilon[l])),
    y     = "Mean Squared Error",
    color = expression(epsilon[e] + epsilon[l])
  ) +
  theme_minimal()

for (n_val in unique(results_4$N)) {
  ggsave(paste0("Figures/sbm_dp_mse_plot_", n_val, ".pdf"),
    p4_base %+% results_4[!is.na(mse) & N == n_val],
    width = 5, height = 5
  )
}
cat("Saved: Figures/sbm_dp_mse_plot_2000.pdf and sbm_dp_mse_plot_20000.pdf\n")

# ==============================================================================
# SIMULATION 5: MSE vs Network Size (Erdos-Renyi, constant degree)
# ==============================================================================
cat("\n=== Simulation 5: MSE vs Network Size ===\n")

NETWORK_SIZES <- c(500, 1000, 5000, 10000)
EPSILON_VALUES <- c(0.5, 1, 2, 4, 8)
NUM_ITER_5 <- 10000
DEG_CONST <- 20

results_list_5 <- list()
start_time <- Sys.time()

for (n in NETWORK_SIZES) {
  p <- min(DEG_CONST / (n - 1), 1)
  cat(sprintf("  N = %d (p = %.5f)\n", n, p))

  for (iter in seq_len(NUM_ITER_5)) {
    g <- sample_gnp(n, p, directed = FALSE, loops = FALSE)
    adj <- as_adjacency_matrix(g, sparse = TRUE)
    rA <- row_normalize_sparse(adj)

    n0 <- floor(n / 2)
    caste_vector <- sample(c(rep(0, n0), rep(1, n - n0)))

    true_ec <- compute_true_iec(rA, caste_vector)
    if (is.na(true_ec)) next

    noise_seed <- list(U_flip = runif(n), U_laplace = rlaplace(1, 0, 1))

    for (eps in EPSILON_VALUES) {
      ec_priv <- compute_noisy_iec_coupled(rA, caste_vector, eps / 2, eps / 2, noise_seed)
      if (is.na(ec_priv)) ec_priv <- 0

      results_list_5[[length(results_list_5) + 1]] <- list(
        n_nodes       = n,
        epsilon       = eps,
        squared_error = (ec_priv - true_ec)^2
      )
    }

    if (iter %% 1000 == 0) cat(sprintf("    iteration %d\n", iter))
  }
}

cat(sprintf(
  "\nSimulation 5 finished in %.1f minutes.\n",
  as.numeric(difftime(Sys.time(), start_time, units = "mins"))
))

results_5 <- rbindlist(results_list_5)
fwrite(results_5, "DerivedData/mse_simulation_results_ec.csv")

summary_5 <- results_5[, .(mse = mean(squared_error, na.rm = TRUE)), by = .(n_nodes, epsilon)]
log_min <- floor(log10(min(summary_5$mse)))
log_max <- ceiling(log10(max(summary_5$mse)))
y_breaks <- 10^(seq(log_min, log_max))

p5 <- ggplot(summary_5, aes(
  x = n_nodes, y = mse,
  color = factor(epsilon), group = epsilon
)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10(breaks = NETWORK_SIZES) +
  scale_y_log10(
    limits = c(min(summary_5$mse), max(summary_5$mse)),
    breaks = y_breaks,
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  labs(
    x     = "Number of Nodes (Log Scale)",
    y     = "Mean Squared Error (Log Scale)",
    color = "Epsilon"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("Figures/mse_vs_size_constant_degree_ec.pdf", p5, width = 5, height = 5)
cat("Saved: Figures/mse_vs_size_constant_degree_ec.pdf\n")
