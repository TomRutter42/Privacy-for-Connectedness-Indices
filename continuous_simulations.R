# Continuous Label Simulations
# Three simulations studying DP estimation of regression slopes on graphon networks.
#
# 1. Epsilon allocation sweep  — MSE vs label/edge budget split, N = 20k–500k
# 2. Regression slope          — MSE of slope vs homophily, N = 10k–100k
# 3. Integral estimation       — MSE of integral under regression line and R² vs homophily, N = 10k–100k
#
# Simulations 2 and 3 share the same network generation loop.

library(data.table)
library(ggplot2)
library(VGAM)
library(igraph)
library(Matrix)
library(parallel)
library(scales)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

get_truncated_laplace_variance <- function(lambda, A) {
  term1 <- 2 * lambda^2
  numerator <- exp(-A / lambda) * (A^2 + 2 * lambda * A)
  denominator <- 1 - exp(-A / lambda)
  term1 - (numerator / denominator)
}

r_truncated_laplace <- function(n, lambda, A) {
  u <- runif(n)
  max_cdf <- 1 - exp(-A / lambda)
  abs_val <- -lambda * log(1 - u * max_cdf)
  sgn <- sample(c(-1, 1), size = n, replace = TRUE)
  abs_val * sgn
}

algorithm_3_noisystats <- function(x, y, epsilon, bounds) {
  n <- length(x)
  if (n < 2) {
    return(list(beta = NA, alpha = NA))
  }

  bar_x <- mean(x)
  bar_y <- mean(y)
  nvar_x <- sum((x - bar_x)^2)
  ncov_xy <- sum((x - bar_x) * (y - bar_y))

  Delta_1 <- (1 - 1 / n) * (bounds$b - bounds$a)^2
  Delta_2 <- (1 - 1 / n) * (bounds$b - bounds$a) * (bounds$b_prime - bounds$a_prime)

  L_1 <- rlaplace(1, 0, 3 * Delta_1 / epsilon)
  L_2 <- rlaplace(1, 0, 3 * Delta_2 / epsilon)

  denom <- nvar_x + L_2
  num <- ncov_xy + L_1

  if (denom > 1e-10) {
    beta_tilde <- num / denom
    Delta_3 <- (bounds$b_prime - bounds$a_prime) / n + abs(beta_tilde) * (bounds$b - bounds$a) / n
    L_3 <- rlaplace(1, 0, 3 * Delta_3 / epsilon)
    alpha_tilde <- (bar_y - beta_tilde * bar_x) + L_3
    list(beta = beta_tilde, alpha = alpha_tilde)
  } else {
    list(beta = NA, alpha = NA)
  }
}

generate_graphon_network <- function(N, expected_degree, homophily) {
  x_raw <- runif(N)
  mu_h <- if (abs(homophily) < 1e-6) {
    1.0
  } else {
    2 / homophily - 2 / (homophily^2) * (1 - exp(-homophily))
  }
  p_base <- min(expected_degree / ((N - 1) * mu_h), 1.0)

  g_er <- sample_gnp(N, p_base, directed = FALSE)
  edges <- as_edgelist(g_er, names = FALSE)

  if (nrow(edges) > 0) {
    p_keep <- exp(-homophily * abs(x_raw[edges[, 1]] - x_raw[edges[, 2]]))
    edges_final <- edges[runif(length(p_keep)) < p_keep, , drop = FALSE]
    g_final <- graph_from_edgelist(edges_final, directed = FALSE)
  } else {
    g_final <- make_empty_graph(n = N, directed = FALSE)
  }
  list(g = g_final, x_raw = x_raw)
}

# ==============================================================================
# SIMULATION 1: Epsilon allocation sweep (parallelised)
# ==============================================================================
cat("\n=== Simulation 1: Epsilon Allocation Sweep ===\n")

N_values_1 <- c(20000, 50000, 100000, 500000)
NUM_ITER_1 <- 3000
total_epsilon_1 <- 8.0
delta_ell_1 <- 1e-3
expected_deg_1 <- 20
homophily_1 <- 0.8
allocations_1 <- seq(0.1, 0.8, by = 0.1)
CHUNK_SIZE <- 50

process_N_alloc <- function(N) {
  set.seed(123 + N)
  cat(sprintf("  Generating network for N = %d...\n", N))
  gd <- generate_graphon_network(N, expected_deg_1, homophily_1)
  x_raw <- gd$x_raw
  edges <- as_edgelist(gd$g, names = FALSE)
  if (nrow(edges) == 0) {
    return(NULL)
  }

  degrees <- degree(gd$g, v = seq_len(N))
  edge_from <- c(edges[, 1], edges[, 2])
  edge_to <- c(edges[, 2], edges[, 1])

  A_sparse <- sparseMatrix(
    i = edge_from, j = edge_to,
    x = 1.0 / degrees[edge_from],
    dims = c(N, N)
  )
  y_true <- as.numeric(A_sparse %*% x_raw)
  y_true[degrees == 0] <- 0
  true_slope <- coef(lm(y_true ~ x_raw))[2]
  cat(sprintf("  N = %d  True slope: %.4f\n", N, true_slope))

  local_results <- list()
  n_chunks <- ceiling(NUM_ITER_1 / CHUNK_SIZE)

  for (frac_e in allocations_1) {
    eps_e <- total_epsilon_1 * frac_e
    eps_l <- total_epsilon_1 * (1 - frac_e)
    lambda <- 1.0 / eps_l
    A_t <- lambda * log(1 + (exp(eps_l) - 1) / (2 * delta_ell_1))
    sigma_sq <- get_truncated_laplace_variance(lambda, A_t)
    rng <- 1 + 2 * A_t
    Delta <- (1 - 1 / N) * rng^2

    private_slopes <- numeric(NUM_ITER_1)

    for (chunk in seq_len(n_chunks)) {
      idx <- ((chunk - 1) * CHUNK_SIZE + 1):min(chunk * CHUNK_SIZE, NUM_ITER_1)
      m <- length(idx)
      z_mat <- matrix(r_truncated_laplace(N * m, lambda, A_t), nrow = N)
      x_hat_mat <- x_raw + z_mat
      y_hat_mat <- as.matrix(A_sparse %*% x_hat_mat)
      y_hat_mat[degrees == 0, ] <- 0

      bar_x <- colMeans(x_hat_mat)
      bar_y <- colMeans(y_hat_mat)
      x_c <- sweep(x_hat_mat, 2, bar_x, "-")
      y_c <- sweep(y_hat_mat, 2, bar_y, "-")
      nvar_x <- colSums(x_c^2)
      ncov_xy <- colSums(x_c * y_c)

      L1 <- rlaplace(m, 0, 3 * Delta / eps_e)
      L2 <- rlaplace(m, 0, 3 * Delta / eps_e)
      denom <- nvar_x + L2
      hat_beta <- ifelse(denom > 1e-10, (ncov_xy + L1) / denom, NA_real_)

      s2_x <- nvar_x / (N - 1)
      corrected <- ifelse(!is.na(hat_beta) & (s2_x - sigma_sq) > 1e-6,
        hat_beta * s2_x / (s2_x - sigma_sq), hat_beta
      )
      private_slopes[idx] <- corrected
    }

    valid <- private_slopes[!is.na(private_slopes)]
    local_results[[length(local_results) + 1]] <- data.table(
      N            = N,
      frac_edge_dp = frac_e,
      true_slope   = true_slope,
      mse          = mean((valid - true_slope)^2)
    )
  }
  local_results
}

n_cores <- if (.Platform$OS.type == "windows") 1L else min(detectCores() - 1, length(N_values_1))
cat(sprintf("Using %d cores\n", n_cores))
all_results_1 <- mclapply(N_values_1, process_N_alloc, mc.cores = n_cores)
results_1 <- rbindlist(unlist(all_results_1, recursive = FALSE))
results_1[, frac_local_dp := 1 - frac_edge_dp]

p1 <- ggplot(results_1, aes(x = frac_local_dp, y = mse, color = as.factor(N), group = N)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = 1 - allocations_1) +
  scale_color_discrete(labels = function(x) {
    format(as.numeric(x),
      big.mark = ",",
      scientific = FALSE, trim = TRUE
    )
  }) +
  labs(
    x     = expression(Ratio ~ epsilon[l] / (epsilon[l] + epsilon[e])),
    y     = "Mean Squared Error (Log Scale)",
    color = "Number of Nodes"
  ) +
  theme_classic()

ggsave("Figures/graphon_allocation_mse_plot.pdf", p1, width = 5, height = 5)
ggsave("Figures/graphon_allocation_mse_plot.png", p1, width = 5, height = 5)
cat("Saved: Figures/graphon_allocation_mse_plot.pdf\n")

# ==============================================================================
# SIMULATIONS 2 & 3: Graphon regression MSE and integral MSE vs homophily
# (share the same network generation loop)
# ==============================================================================
cat("\n=== Simulations 2 & 3: Regression and Integral vs Homophily ===\n")

N_values_23 <- c(10000, 20000, 50000, 100000)
homophily_values <- seq(0, 1, by = 0.2)
NUM_ITER_reg <- 2500
NUM_ITER_integral <- 2500
epsilon_ell_23 <- 4.0
epsilon_e_23 <- 4.0
delta_ell_23 <- 1e-3
expected_deg_23 <- 20

results_reg <- list()
results_integral <- list()

for (N in N_values_23) {
  for (homophily in homophily_values) {
    set.seed(42 + N + homophily * 100)
    cat(sprintf("  N = %d, homophily = %.1f\n", N, homophily))

    gd <- generate_graphon_network(N, expected_deg_23, homophily)
    x_raw <- gd$x_raw
    edges <- as_edgelist(gd$g, names = FALSE)
    if (nrow(edges) == 0) {
      cat("  Empty graph, skipping.\n")
      next
    }

    degrees <- degree(gd$g, v = seq_len(N))
    edge_from <- c(edges[, 1], edges[, 2])
    edge_to <- c(edges[, 2], edges[, 1])

    node_sums <- rowsum(x_raw[edge_to], edge_from)
    y_true <- numeric(N)
    active <- as.integer(rownames(node_sums))
    y_true[active] <- node_sums[, 1] / degrees[active]
    y_true[is.na(y_true) | is.infinite(y_true)] <- 0

    true_model <- lm(y_true ~ x_raw)
    true_slope <- coef(true_model)[2]
    true_alpha <- coef(true_model)[1]
    true_r2 <- summary(true_model)$r.squared
    true_integral <- 0.25 * true_alpha + 0.03125 * true_slope

    lambda <- 1.0 / epsilon_ell_23
    A <- lambda * log(1 + (exp(epsilon_ell_23) - 1) / (2 * delta_ell_23))
    sigma_sq <- get_truncated_laplace_variance(lambda, A)
    bounds <- list(a = -A, b = 1 + A, a_prime = -A, b_prime = 1 + A)

    NUM_ITER_23 <- max(NUM_ITER_reg, NUM_ITER_integral)
    private_slopes <- numeric(NUM_ITER_23)
    private_integrals <- numeric(NUM_ITER_23)

    for (iter in seq_len(NUM_ITER_23)) {
      z <- r_truncated_laplace(N, lambda, A)
      x_hat <- x_raw + z

      ns <- rowsum(x_hat[edge_to], edge_from)
      y_hat <- numeric(N)
      y_hat[active] <- ns[, 1] / degrees[active]
      y_hat[is.na(y_hat) | is.infinite(y_hat)] <- 0

      res <- algorithm_3_noisystats(x_hat, y_hat, epsilon_e_23, bounds)
      hat_beta <- res$beta
      hat_alpha <- res$alpha

      if (!is.na(hat_beta)) {
        s2_x <- var(x_hat)
        if ((s2_x - sigma_sq) > 1e-6) {
          tilde_beta <- hat_beta * s2_x / (s2_x - sigma_sq)
          private_slopes[iter] <- tilde_beta
          if (!is.na(hat_alpha)) {
            tilde_alpha <- hat_alpha - (tilde_beta - hat_beta) * mean(x_hat)
            private_integrals[iter] <- 0.25 * tilde_alpha + 0.03125 * tilde_beta
          } else {
            private_integrals[iter] <- NA
          }
        } else {
          private_slopes[iter] <- hat_beta
          private_integrals[iter] <- if (!is.na(hat_alpha)) 0.25 * hat_alpha + 0.03125 * hat_beta else NA
        }
      } else {
        private_slopes[iter] <- NA
        private_integrals[iter] <- NA
      }
    }

    valid_slopes <- private_slopes[seq_len(NUM_ITER_reg)]
    valid_slopes <- valid_slopes[!is.na(valid_slopes)]
    results_reg[[length(results_reg) + 1]] <- data.table(
      N          = N,
      homophily  = homophily,
      true_slope = true_slope,
      mse        = mean((valid_slopes - true_slope)^2),
      n_valid    = length(valid_slopes)
    )

    valid_integrals <- private_integrals[seq_len(NUM_ITER_integral)]
    valid_integrals <- valid_integrals[!is.na(valid_integrals)]
    results_integral[[length(results_integral) + 1]] <- data.table(
      N              = N,
      homophily      = homophily,
      true_r_squared = true_r2,
      true_integral  = true_integral,
      mse_integral   = mean((valid_integrals - true_integral)^2)
    )
  }
}

# --- Simulation 2: Regression MSE ---
results_reg_dt <- rbindlist(results_reg)
fwrite(results_reg_dt, "DerivedData/graphon_regression_mse_results.csv")

p2 <- ggplot(
  results_reg_dt[!is.na(mse)],
  aes(x = homophily, y = mse, color = as.factor(N), group = N)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_color_discrete(labels = function(x) {
    format(as.numeric(x),
      big.mark = ",",
      scientific = FALSE, trim = TRUE
    )
  }) +
  labs(
    x     = "Homophily Parameter",
    y     = "Mean Squared Error",
    color = "Number of Nodes"
  ) +
  theme_minimal()

ggsave("Figures/graphon_regression_mse_plot.pdf", p2, width = 5, height = 5)
ggsave("Figures/graphon_regression_mse_plot.png", p2, width = 5, height = 5)
cat("Saved: Figures/graphon_regression_mse_plot.pdf\n")

# --- Simulation 3: Integral MSE and R² ---
results_integral_dt <- rbindlist(results_integral)

p3_r2 <- ggplot(
  results_integral_dt[N == max(N_values_23)],
  aes(x = homophily, y = true_r_squared)
) +
  geom_line(color = "dodgerblue4", linewidth = 1) +
  geom_point(size = 3, color = "dodgerblue4") +
  labs(x = "Homophily Parameter", y = expression(R^2)) +
  theme_classic() +
  scale_y_continuous(labels = percent_format(accuracy = 1))

ggsave("Figures/graphon_homophily_rsquared_plot.pdf", p3_r2, width = 5, height = 5)
ggsave("Figures/graphon_homophily_rsquared_plot.png", p3_r2, width = 5, height = 5)

p3_mse <- ggplot(
  results_integral_dt[!is.na(mse_integral)],
  aes(x = homophily, y = mse_integral, color = as.factor(N), group = N)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_color_discrete(labels = function(x) {
    format(as.numeric(x),
      big.mark = ",",
      scientific = FALSE, trim = TRUE
    )
  }) +
  labs(
    x     = "Homophily Parameter",
    y     = "Mean Squared Error of Integral (Log Scale)",
    color = "Number of Nodes"
  ) +
  theme_minimal()

ggsave("Figures/graphon_integral_mse_plot.pdf", p3_mse, width = 5, height = 5)
ggsave("Figures/graphon_integral_mse_plot.png", p3_mse, width = 5, height = 5)
cat("Saved: Figures/graphon_integral_mse_plot.pdf and graphon_homophily_rsquared_plot.pdf\n")
