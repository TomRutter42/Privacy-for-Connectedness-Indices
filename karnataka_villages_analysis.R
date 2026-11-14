# Karnataka Villages Analysis
# Loads village networks with caste data, runs DP analysis, and generates figures.
# Caste coding: 0 = SC/ST/Minority, 1 = General/OBC, NA = blank

source("replication_helpers.R", local = TRUE)
require_packages(c("haven", "dplyr", "ggplot2", "Matrix", "data.table"))
library(Matrix)
library(data.table)
start_experiment("villages", 2026110720L)
ext_data_dir <- config$data

library(haven)
library(dplyr)
library(ggplot2)

# --- Configuration ---
data_dir <- file.path(ext_data_dir, "datav4.0/datav4.0/Data")
adj_dir <- file.path(data_dir, "1. Network Data/Adjacency Matrices")
key_dir <- file.path(data_dir, "1. Network Data/Adjacency Matrix Keys")
demo_file <- file.path(data_dir, "2. Demographics and Outcomes/household_characteristics.dta")

# --- Helper Functions ---

clean_village_network <- function(village_data) {
  caste <- village_data$caste_vector
  A <- village_data$adj_matrix

  keep_idx <- which(!is.na(caste))

  if (length(keep_idx) < 2) {
    return(NULL)
  }

  A_clean <- A[keep_idx, keep_idx, drop = FALSE]
  caste_clean <- caste[keep_idx]

  list(
    adj_matrix = A_clean,
    caste_vector = caste_clean,
    n_households = length(caste_clean),
    n_caste_0 = sum(caste_clean == 0),
    n_caste_1 = sum(caste_clean == 1)
  )
}

load_village_network <- function(village_num, adj_dir, key_dir, hh_data) {
  adj_file <- file.path(adj_dir, paste0("adj_allVillageRelationships_HH_vilno_", village_num, ".csv"))
  key_file <- file.path(key_dir, paste0("key_HH_vilno_", village_num, ".csv"))

  if (!file.exists(adj_file) || !file.exists(key_file)) {
    warning(paste("Missing files for village", village_num))
    return(NULL)
  }

  adj_matrix <- as.matrix(read.csv(adj_file, header = FALSE))
  key_df <- read.csv(key_file, header = FALSE)
  hh_keys <- key_df$V1

  n <- length(hh_keys)
  if (nrow(adj_matrix) != n || ncol(adj_matrix) != n) {
    warning(paste(
      "Dimension mismatch for village", village_num,
      "- matrix:", nrow(adj_matrix), "x", ncol(adj_matrix),
      ", keys:", n
    ))
    return(NULL)
  }

  village_hh <- hh_data %>%
    filter(village == village_num) %>%
    select(adjmatrix_key, caste_binary, castesubcaste)

  caste_vector <- sapply(hh_keys, function(key) {
    match_row <- village_hh %>% filter(adjmatrix_key == key)
    if (nrow(match_row) == 0) {
      return(NA_real_)
    } else {
      return(match_row$caste_binary[1])
    }
  })

  list(
    village = village_num,
    n_households = n,
    adj_matrix = adj_matrix,
    caste_vector = caste_vector,
    n_edges = sum(adj_matrix) / 2,
    n_caste_0 = sum(caste_vector == 0, na.rm = TRUE),
    n_caste_1 = sum(caste_vector == 1, na.rm = TRUE),
    n_caste_NA = sum(is.na(caste_vector))
  )
}

# --- Main Script ---

# 1. Load household characteristics
check_data_files(village_input_files(config$data))
cat("Loading household characteristics...\n")
hh_data <- read_dta(demo_file)

hh_data <- hh_data %>%
  mutate(
    caste_binary = case_when(
      castesubcaste %in% c("SCHEDULE CASTE", "SCHEDULE TRIBE", "MINORITY") ~ 0,
      castesubcaste %in% c("GENERAL", "OBC") ~ 1,
      TRUE ~ NA_real_
    )
  )

# 2. Get list of available villages
adj_files <- list.files(adj_dir, pattern = "adj_allVillageRelationships_HH_vilno_\\d+\\.csv")
village_nums <- as.numeric(gsub("adj_allVillageRelationships_HH_vilno_(\\d+)\\.csv", "\\1", adj_files))
village_nums <- sort(village_nums)

cat("\nFound", length(village_nums), "villages with adjacency matrices\n")

# 3. Load all villages
cat("\nLoading village networks...\n")
village_networks <- list()

for (v in village_nums) {
  result <- load_village_network(v, adj_dir, key_dir, hh_data)
  if (!is.null(result)) {
    as_char_v <- as.character(v)
    village_networks[[as_char_v]] <- result
  }
  if (v %% 10 == 0) cat("Processed village", v, "\n")
}

cat("\nSuccessfully loaded", length(village_networks), "village networks\n")

# 4. Filter villages with usable caste data (>50% coverage and >6 historically
#    disadvantaged households, following the sample restrictions in the paper)
village_networks_with_caste <- village_networks[sapply(village_networks, function(x) {
  coverage <- (x$n_caste_0 + x$n_caste_1) / x$n_households
  coverage > 0.5 && x$n_caste_0 > 6
})]

cat("\nVillages meeting sample restrictions:", length(village_networks_with_caste), "\n")

# --- Village Inclusion/Exclusion Diagnostic ---
cat("\n--- Village Inclusion/Exclusion Diagnostic ---\n")
village_diagnostic <- do.call(rbind, lapply(village_networks, function(x) {
  coverage <- (x$n_caste_0 + x$n_caste_1) / x$n_households
  status <- "Included"
  reason <- "N/A"

  if (coverage <= 0.5) {
    status <- "Excluded"
    reason <- if (coverage == 0) "Zero Coverage" else paste0("Low Coverage (", round(100 * coverage, 1), "%)")
  } else if (x$n_caste_0 <= 6) {
    status <- "Excluded"
    reason <- paste0("Too few historically disadvantaged households (n = ", x$n_caste_0, ")")
  }

  data.frame(
    village = x$village, n_hh = x$n_households, n_caste0 = x$n_caste_0,
    n_caste1 = x$n_caste_1, coverage_pct = round(100 * coverage, 1),
    status = status, reason = reason
  )
}))
print(village_diagnostic, row.names = FALSE)
cat("\nSummary of Coverage:\n")
print(table(village_diagnostic$status, village_diagnostic$reason))

# --- DP Parameters ---
epsilon_values <- c(0.5, 1, 2, 4, 6, 8)
num_iterations <- iteration_count(1000L)

cat("\n")
cat("============================================================\n")
cat("DIFFERENTIAL PRIVACY ANALYSIS FOR CASTE-BASED CONNECTEDNESS\n")
cat("============================================================\n")

villages_to_process <- names(village_networks_with_caste)

# Storage for results
results_list <- list()
raw_data_list <- list()
plot_data_list <- list()
long_results <- list()

for (v_name in villages_to_process) {
  set.seed(2026110720L + as.integer(v_name))
  v_data <- village_networks_with_caste[[v_name]]

  v_clean <- clean_village_network(v_data)

  if (is.null(v_clean) || v_clean$n_caste_0 < 2 || v_clean$n_caste_1 < 2) {
    next
  }

  validate_adjacency(v_clean$adj_matrix)
  rA <- row_normalize(Matrix(v_clean$adj_matrix, sparse = TRUE))
  caste <- v_clean$caste_vector
  n <- v_clean$n_households

  true_ec <- true_connectedness(rA, caste)

  if (is.na(true_ec)) {
    next
  }

  sq_errors <- matrix(NA, nrow = length(epsilon_values), ncol = num_iterations)
  noisy_ec_vals <- matrix(NA, nrow = length(epsilon_values), ncol = num_iterations)
  first_run_eps_8 <- NA
  stabilized <- matrix(FALSE, nrow = length(epsilon_values), ncol = num_iterations)

  for (iter in 1:num_iterations) {
    noise_seed <- list(
      U_flip = runif(n),
      U_laplace = r_laplace(1)
    )

    for (i in seq_along(epsilon_values)) {
      eps <- epsilon_values[i]
      eps_s <- eps / 2
      eps_e <- eps / 2

      private <- private_connectedness(rA, caste, eps_s, eps_e,
                                       noise_seed$U_flip, noise_seed$U_laplace, details = TRUE)
      noisy_ec <- private$estimate
      stabilized[i, iter] <- private$stabilized

      if (eps == 8 && iter == 1) {
        first_run_eps_8 <- noisy_ec
      }

      if (!is.na(noisy_ec) && is.finite(noisy_ec)) {
        sq_errors[i, iter] <- (noisy_ec - true_ec)^2
        noisy_ec_vals[i, iter] <- noisy_ec
      }
    }
  }

  mse_values <- rowMeans(sq_errors, na.rm = TRUE)
  sd_values <- apply(noisy_ec_vals, 1, sd, na.rm = TRUE)

  avg_degree <- sum(v_clean$adj_matrix) / v_clean$n_households

  results_list[[v_name]] <- data.frame(
    village = as.integer(v_name),
    n_households = v_clean$n_households,
    n_caste_0 = v_clean$n_caste_0,
    n_caste_1 = v_clean$n_caste_1,
    avg_degree = avg_degree,
    true_ec = true_ec,
    mse_eps_8 = mse_values[match(8, epsilon_values)],
    sd_eps_8 = sd_values[match(8, epsilon_values)]
  )

  raw_data_list[[v_name]] <- noisy_ec_vals
  long_results[[v_name]] <- data.frame(
    village = as.integer(v_name), epsilon = epsilon_values,
    epsilon_label = epsilon_values / 2, epsilon_edge = epsilon_values / 2,
    delta_label = 0, true_ec = true_ec, mse = mse_values, sd = sd_values,
    bias = rowMeans(noisy_ec_vals) - true_ec,
    n_attempted = num_iterations, n_stabilized = rowSums(stabilized))

  if (!is.na(first_run_eps_8)) {
    plot_data_list[[v_name]] <- data.frame(
      village = v_name,
      true_ec = true_ec,
      noisy_ec = first_run_eps_8,
      n_households = v_clean$n_households
    )
  }

  if (as.integer(v_name) %% 10 == 0) {
    cat("Processed village", v_name, "\n")
  }
}

# Combine results
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

cat("\n--- Results Summary ---\n")
cat("Villages processed:", nrow(results_df), "\n\n")

# --- Signal vs Noise Variance Analysis ---
cat("--- Signal vs Noise Variance ---\n")

signal_var <- var(results_df$true_ec, na.rm = TRUE)
cat("Signal Variance (Var of True EC across villages):", round(signal_var, 6), "\n\n")

long_summary <- rbindlist(long_results)
signal_noise <- long_summary[, .(
  signal_variance = signal_var, mean_noise_variance = mean(sd^2),
  mean_mse = mean(mse), mean_noise_sd = mean(sd),
  variance_ratio = signal_var / mean(sd^2), mse_ratio = signal_var / mean(mse)), by = epsilon]
print(signal_noise)

# --- Summary Table ---
cat("\nExporting summary table to LaTeX...\n")
tex_table <- results_df[, c("village", "n_households", "n_caste_0", "n_caste_1", "avg_degree", "true_ec", "sd_eps_8")]
tex_table$avg_degree <- sprintf("%.1f", tex_table$avg_degree)
tex_table$true_ec    <- sprintf("%.2f", tex_table$true_ec)
tex_table$sd_eps_8   <- sprintf("%.2f", tex_table$sd_eps_8)

headers <- c("Village \\#", "\\# HH", "\\# Hist.\ Dis.", "\\# Hist.\ Non-Dis.", "Avg.\ Deg.", "Cross-Caste Conn.", "SD ($\\varepsilon$=8)")

tex_con <- file(output_path("Tables", "village_summary.tex"), "w")
cat("\\begin{tabular}{lrrrrrr}\n", file = tex_con)
cat("\\hline\n", file = tex_con)
cat(paste(headers, collapse = " & "), "\\\\\n", file = tex_con)
cat("\\hline\n", file = tex_con)
for (i in 1:nrow(tex_table)) {
  cat(paste(tex_table[i, ], collapse = " & "), "\\\\\n", file = tex_con)
}
cat("\\hline\n", file = tex_con)
cat("\\end{tabular}\n", file = tex_con)
close(tex_con)
cat("Table saved to: Tables/village_summary.tex\n")

write.csv(results_df, output_path("DerivedData", "village_caste_dp_results.csv"), row.names = FALSE)
cat("Results saved to: DerivedData/village_caste_dp_results.csv\n")

# --- True vs Noisy EC Scatter for the Analysis Sample ---
cat("\nGenerating chart: True vs Noisy EC (Epsilon = 8)...\n")
plot_data <- do.call(rbind, plot_data_list)

if (!is.null(plot_data) && nrow(plot_data) > 0) {
  cor_val <- cor(plot_data$true_ec, plot_data$noisy_ec, use = "complete.obs")

  p_scatter <- ggplot(plot_data, aes(x = true_ec, y = noisy_ec)) +
    geom_point(aes(size = n_households), color = "dodgerblue4", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "twodash", color = "grey50", linewidth = 0.8) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(
      subtitle = sprintf("Correlation Coefficient = %.2f", cor_val),
      x = "True Across-Caste Connectedness",
      y = "Private Across-Caste Connectedness (Epsilon = 8)",
      size = "Households"
    ) +
    theme_classic() +
    coord_fixed(ratio = 1) +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(hjust = 0.5)
    )

  ggsave(output_path("Figures", "True_vs_Noisy_EC_Village_No_56.pdf"), p_scatter, width = 5, height = 5)
  ggsave(output_path("Figures", "True_vs_Noisy_EC_Village_No_56.png"), p_scatter, width = 5, height = 5)
  cat("Chart saved to: Figures/True_vs_Noisy_EC_Village_No_56.pdf and .png\n")
} else {
  cat("Warning: No data for chart (plot_data is empty)\n")
}

# --- Overlaid Correlation Histogram (Epsilon = 4, 6, 8) ---
cat("\nGenerating overlaid histogram of correlations (Epsilon = 4, 6, 8)...\n")

true_ec_vec <- sapply(results_list, function(x) x$true_ec)

eps4_idx <- which(epsilon_values == 4)
eps6_idx <- which(epsilon_values == 6)
eps8_idx <- which(epsilon_values == 8)

if (length(eps4_idx) > 0 && length(eps6_idx) > 0 && length(eps8_idx) > 0) {
  get_correlations <- function(eps_idx) {
    cor_vals <- numeric(num_iterations)
    n_eps <- length(epsilon_values)
    noisy_ec_run <- numeric(length(raw_data_list))

    for (iter in 1:num_iterations) {
      for (v_idx in seq_along(raw_data_list)) {
        lin_idx <- eps_idx + (iter - 1) * n_eps
        val <- raw_data_list[[v_idx]][lin_idx]
        if (length(val) == 0 || is.na(val)) val <- NA
        noisy_ec_run[v_idx] <- as.numeric(val)
      }
      cor_vals[iter] <- complete_correlation(true_ec_vec, noisy_ec_run)
    }
    return(cor_vals)
  }

  cor_eps4 <- get_correlations(eps4_idx)
  cor_eps6 <- get_correlations(eps6_idx)
  cor_eps8 <- get_correlations(eps8_idx)

  plot_df <- rbind(
    data.frame(Correlation = cor_eps4, Epsilon = "4"),
    data.frame(Correlation = cor_eps6, Epsilon = "6"),
    data.frame(Correlation = cor_eps8, Epsilon = "8")
  )

  mean_4 <- mean(cor_eps4, na.rm = TRUE)
  mean_6 <- mean(cor_eps6, na.rm = TRUE)
  mean_8 <- mean(cor_eps8, na.rm = TRUE)
  cat(sprintf("Means: Eps 4 = %.2f, Eps 6 = %.2f, Eps 8 = %.2f\n", mean_4, mean_6, mean_8))

  min_val <- min(plot_df$Correlation, na.rm = TRUE)
  max_val <- max(plot_df$Correlation, na.rm = TRUE)
  padding <- (max_val - min_val) * 0.15

  p_overlaid <- ggplot(plot_df, aes(x = Correlation, fill = Epsilon)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5, position = "identity", color = "white") +
    scale_fill_manual(values = c("4" = "red", "6" = "green3", "8" = "dodgerblue4")) +
    labs(
      x = "Correlation Coefficient (True vs Private Connectedness)",
      y = "Frequency Density"
    ) +
    coord_cartesian(xlim = c(min_val - padding, max_val + padding)) +
    theme_classic() +
    theme(
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "top"
    )

  ggsave(output_path("Figures", "Correlation_Histogram_Overlaid_Eps4_6_8.pdf"), plot = p_overlaid, width = 5, height = 5)
  ggsave(output_path("Figures", "Correlation_Histogram_Overlaid_Eps4_6_8.png"), plot = p_overlaid, width = 5, height = 5)
  cat("Chart saved to: Figures/Correlation_Histogram_Overlaid_Eps4_6_8.pdf and .png\n")
} else {
  cat("Warning: One or more epsilons (4, 6, 8) missing from simulation values.\n")
}

# --- Boxplot: Largest Village Connectivity ---
cat("\nGenerating boxplot for largest village...\n")

largest_village <- NULL
max_hh <- 0
for (v_data in village_networks) {
  coverage <- (v_data$n_caste_0 + v_data$n_caste_1) / v_data$n_households
  if (coverage > 0.5 && v_data$n_households > max_hh) {
    max_hh <- v_data$n_households
    largest_village <- v_data
  }
}

cat(sprintf("Largest village is %s with %d households.\n", largest_village$village, max_hh))

v_clean_box <- clean_village_network(largest_village)
rA_box <- row_normalize(v_clean_box$adj_matrix)
exposure <- as.numeric(rA_box %*% v_clean_box$caste_vector)

boxplot_data <- data.frame(
  exposure = exposure,
  caste = factor(v_clean_box$caste_vector,
    levels = c(0, 1),
    labels = c("Historically\nDisadvantaged", "Historically\nNon-Disadvantaged")
  )
)

fill_colors <- c("Historically\nDisadvantaged" = "#E69F00", "Historically\nNon-Disadvantaged" = "#56B4E9")
point_colors <- c("Historically\nDisadvantaged" = "#D55E00", "Historically\nNon-Disadvantaged" = "#0072B2")

set.seed(2026110720L)
p_box <- ggplot(boxplot_data, aes(x = exposure, y = caste, fill = caste, color = caste)) +
  geom_jitter(height = 0.15, alpha = 0.4, size = 1.2) +
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.shape = NA, color = "black") +
  labs(
    x = "Proportion of connections to\nhistorically non-disadvantaged\nhouseholds",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linetype = "dashed"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 12), color = "black"),
    plot.margin = margin(15, 20, 15, 10)
  )

ggsave(output_path("Figures", "Largest_Village_Connectivity_Boxplot.pdf"), p_box, width = 5, height = 5)
ggsave(output_path("Figures", "Largest_Village_Connectivity_Boxplot.png"), p_box, width = 5, height = 5)
cat("Village figures and table complete.\n")
