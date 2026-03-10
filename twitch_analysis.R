# Twitch Gamers Analysis
# 1. Binscatter of average friend view rank vs own view rank (English sub-network)
# 2. DP simulation of slope estimator across all language sub-networks

library(data.table)
library(igraph)
library(ggplot2)
library(ggrepel)
library(VGAM)
library(dplyr)
library(scales)

options(scipen = 999)

# --- Configuration ---
meta_file  <- file.path(ext_data_dir, "large_twitch_features.csv")
graph_file <- file.path(ext_data_dir, "large_twitch_edges.csv")

NUM_ITERATIONS   <- 1000
epsilons_to_test <- c(4, 6, 8)
delta_ell        <- 10^(-3)
bounds <- list(a = 0, b = 1, a_prime = 0, b_prime = 1)

lang_map <- c(
  "EN" = "English",  "FR" = "French",  "KO" = "Korean",   "JA" = "Japanese",
  "RU" = "Russian",  "PL" = "Polish",  "DE" = "German",   "ES" = "Spanish",
  "IT" = "Italian",  "PT" = "Portuguese", "TR" = "Turkish", "ZH" = "Chinese",
  "SV" = "Swedish",  "NL" = "Dutch",   "TH" = "Thai",     "CS" = "Czech",
  "DA" = "Danish",   "HU" = "Hungarian", "FI" = "Finnish", "NO" = "Norwegian",
  "OTHER" = "Other"
)

# --- Helper Functions ---
get_truncated_laplace_variance <- function(lambda, A) {
  term1       <- 2 * lambda^2
  numerator   <- exp(-A / lambda) * (A^2 + 2 * lambda * A)
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
  if (n < 2) return(NA_real_)

  bar_x   <- mean(x); bar_y <- mean(y)
  nvar_x  <- sum((x - bar_x)^2)
  ncov_xy <- sum((x - bar_x) * (y - bar_y))

  Delta_1 <- (1 - 1 / n) * (bounds$b - bounds$a)^2
  Delta_2 <- (1 - 1 / n) * (bounds$b - bounds$a) * (bounds$b_prime - bounds$a_prime)

  L_1 <- rlaplace(1, 0, 3 * Delta_1 / epsilon)
  L_2 <- rlaplace(1, 0, 3 * Delta_2 / epsilon)

  denom <- nvar_x + L_2
  num   <- ncov_xy + L_1

  if (denom > 1e-10) num / denom else NA_real_
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading features from:", meta_file, "\n")
features <- fread(meta_file)
cat("Loading edges from:", graph_file, "\n")
edges <- fread(graph_file)

cat("Loaded", nrow(features), "users and", nrow(edges), "edges.\n")
cat("Analyzing", length(unique(features$language)), "unique languages.\n")

# ==============================================================================
# 2. BINSCATTER: English sub-network
# ==============================================================================

cat("\nBuilding English sub-network for binscatter...\n")
setkey(features, numeric_id)

en_users <- features[language == "EN"]
en_ids   <- en_users$numeric_id
en_edges <- edges[numeric_id_1 %in% en_ids & numeric_id_2 %in% en_ids]

g       <- graph_from_data_frame(en_edges, directed = FALSE, vertices = en_users[, .(numeric_id)])
degrees_g <- degree(g)
en_users$deg <- degrees_g[match(en_users$numeric_id, V(g)$name)]
en_users[is.na(deg), deg := 0]

en_users$view_rank <- (rank(en_users$views, ties.method = "average") - 1) / (nrow(en_users) - 1)

val_lookup  <- setNames(en_users$view_rank, as.character(en_users$numeric_id))
edge_df     <- as.data.frame(as_edgelist(g))
edge_to_en  <- c(edge_df$V2, edge_df$V1)
edge_from_en <- c(edge_df$V1, edge_df$V2)

neighbor_vals <- val_lookup[edge_to_en]
node_sums_en  <- rowsum(neighbor_vals, edge_from_en)

y_true_en <- numeric(nrow(en_users))
names(y_true_en) <- as.character(en_users$numeric_id)
active_ids <- rownames(node_sums_en)
y_true_en[active_ids] <- node_sums_en[, 1] / en_users$deg[match(as.integer(active_ids), en_users$numeric_id)]
y_true_en[is.na(y_true_en) | is.infinite(y_true_en)] <- 0
en_users$y_true <- y_true_en[as.character(en_users$numeric_id)]

valid_users <- en_users[deg > 0]
valid_users[, bin := ntile(view_rank, 20)]
bin_data <- valid_users[, .(
  mean_x = mean(view_rank),
  mean_y = mean(y_true),
  se_y   = sd(y_true) / sqrt(.N)
), by = bin]
setorder(bin_data, bin)

cor_val <- cor(valid_users$view_rank, valid_users$y_true)

p_binscatter <- ggplot(bin_data, aes(x = mean_x, y = mean_y)) +
  geom_point(size = 3, color = "dodgerblue4") +
  geom_errorbar(aes(ymin = mean_y - 1.96 * se_y, ymax = mean_y + 1.96 * se_y),
                width = 0.01, color = "dodgerblue4") +
  geom_smooth(data = valid_users, aes(x = view_rank, y = y_true),
              method = "lm", color = "firebrick", se = FALSE, linetype = "dashed") +
  labs(
    x        = "Own View Rank (Percentile)",
    y        = "Average Friend View Rank",
    subtitle = sprintf("Correlation = %.2f", cor_val)
  ) +
  theme_classic() +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(labels = percent)

ggsave("Figures/Twitch_EN_Binscatter.pdf", p_binscatter, width = 5, height = 5)
ggsave("Figures/Twitch_EN_Binscatter.png", p_binscatter, width = 5, height = 5)
cat("Binscatter saved to Figures/Twitch_EN_Binscatter.pdf and .png\n")

# ==============================================================================
# 3. DP SIMULATION: All language sub-networks
# ==============================================================================

cat("\nPreparing all-language data for DP simulation...\n")

features <- features[, x_raw := rank(views, ties.method = "average") / .N, keyby = language]

nodes        <- unique(features$numeric_id)
edges_filt   <- edges[numeric_id_1 %in% nodes & numeric_id_2 %in% nodes]
edges_rev    <- data.table(numeric_id_1 = edges_filt$numeric_id_2, numeric_id_2 = edges_filt$numeric_id_1)
all_edges_dt <- unique(rbind(edges_filt, edges_rev))

setkey(features, numeric_id)
all_edges_lang <- merge(all_edges_dt, features[, .(numeric_id, lang1 = language)], by.x = "numeric_id_1", by.y = "numeric_id")
all_edges_lang <- merge(all_edges_lang, features[, .(numeric_id, lang2 = language)], by.x = "numeric_id_2", by.y = "numeric_id")
all_edges_sub  <- all_edges_lang[lang1 == lang2]
all_edges_sub[, language := lang1]

cat("Computing true slopes by language...\n")
neighbor_info <- all_edges_sub[, .(from = numeric_id_1, neighbor = numeric_id_2)]
neighbor_info <- merge(neighbor_info, features[, .(numeric_id, n_trait = x_raw)], by.x = "neighbor", by.y = "numeric_id")
node_sums_all <- neighbor_info[, .(y_true = sum(n_trait)), by = from]

degree_counts <- all_edges_sub[, .N, by = numeric_id_1]
setnames(degree_counts, "N", "degree")

features_indexed <- copy(features)
features_indexed <- merge(features_indexed, degree_counts, by.x = "numeric_id", by.y = "numeric_id_1", all.x = TRUE)
features_indexed <- merge(features_indexed, node_sums_all, by.x = "numeric_id", by.y = "from", all.x = TRUE)
features_indexed[is.na(degree), degree := 0]
features_indexed[is.na(y_true),  y_true  := 0]
features_indexed[degree > 0, y_true := y_true / degree]

true_slopes <- features_indexed[, .(
  true_beta = if (.N > 1 && var(x_raw) > 1e-10) coef(lm(y_true ~ x_raw))[2] else 0
), by = language]
cat("True slopes computed for", nrow(true_slopes), "languages.\n")

# --- DP simulation loop ---
all_eps_results  <- list()
plot_data_scatter <- NULL

for (eps in epsilons_to_test) {
  lambda   <- 1.0 / eps
  A        <- lambda * log(1 + (exp(eps) - 1) / (2 * delta_ell))
  sigma_sq <- get_truncated_laplace_variance(lambda, A)

  cat("\nRunning DP simulation for epsilon =", eps, "...\n")
  correlations <- numeric(NUM_ITERATIONS)

  for (iter in 1:NUM_ITERATIONS) {
    if (iter %% 100 == 0) cat("  Iter", iter, "/", NUM_ITERATIONS, "\n")

    noise_vec       <- r_truncated_laplace(nrow(features_indexed), lambda, A)
    features_iter   <- copy(features_indexed)
    features_iter[, x_hat := x_raw + noise_vec]

    ni_iter  <- merge(neighbor_info[, .(from, neighbor)],
                      features_iter[, .(numeric_id, n_trait = x_hat)],
                      by.x = "neighbor", by.y = "numeric_id")
    ns_iter  <- ni_iter[, .(y_hat_sum = sum(n_trait)), by = from]
    features_iter <- merge(features_iter, ns_iter, by.x = "numeric_id", by.y = "from", all.x = TRUE)
    features_iter[is.na(y_hat_sum), y_hat_sum := 0]
    features_iter[, y_hat := ifelse(degree > 0, y_hat_sum / degree, 0)]

    slopes_lang <- features_iter[, .(
      hat_beta = algorithm_3_noisystats(x_hat, y_hat, eps, bounds),
      s2_x     = var(x_hat)
    ), by = language]

    slopes_lang[, tilde_beta := ifelse(
      !is.na(hat_beta) & (s2_x - sigma_sq) > 1e-6,
      hat_beta * (s2_x / (s2_x - sigma_sq)),
      hat_beta
    )]

    slopes_bound <- merge(true_slopes, slopes_lang, by = "language")
    valid        <- slopes_bound[!is.na(tilde_beta) & !is.na(true_beta)]

    if (eps == 8 && iter == 1) plot_data_scatter <- copy(valid)

    correlations[iter] <- if (nrow(valid) > 2) cor(valid$true_beta, valid$tilde_beta) else NA
  }

  all_eps_results[[as.character(eps)]] <- data.table(Correlation = correlations, Epsilon = factor(eps))
}

results_df <- rbindlist(all_eps_results)
fwrite(results_df, "DerivedData/twitch_dp_results.csv")

# Correlation histogram
cor_sums <- results_df[, .(mean_cor = mean(Correlation, na.rm = TRUE)), by = Epsilon]
cat("\nMean correlations:", paste(cor_sums$Epsilon, "=", round(cor_sums$mean_cor, 2), collapse = ", "), "\n")

p_hist <- ggplot(results_df, aes(x = Correlation, fill = Epsilon)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.5,
                 position = "identity", color = "white") +
  scale_fill_manual(values = c("4" = "red", "6" = "green3", "8" = "dodgerblue4")) +
  labs(
    x = "Correlation Coefficient (True vs Private Slope across Languages)",
    y = "Frequency Density"
  ) +
  theme_classic() +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 11), legend.position = "top")

ggsave("Figures/Twitch_Correlation_Histogram_Overlaid_Eps4_6_8.pdf", p_hist, width = 5, height = 5)
ggsave("Figures/Twitch_Correlation_Histogram_Overlaid_Eps4_6_8.png", p_hist, width = 5, height = 5)
cat("Histogram saved to Figures/Twitch_Correlation_Histogram_Overlaid_Eps4_6_8.pdf and .png\n")

# Scatter plot: single run at epsilon = 8
cat("\nGenerating scatter plot for single run (Eps = 8)...\n")
sizes <- features_indexed[, .(N = .N), by = language]
plot_data_scatter <- merge(plot_data_scatter, sizes, by = "language", all.x = TRUE)
plot_data_scatter[, language_full := lang_map[language]]

cor_scatter <- cor(plot_data_scatter$true_beta, plot_data_scatter$tilde_beta)

p_scatter <- ggplot(plot_data_scatter, aes(x = true_beta, y = tilde_beta)) +
  geom_point(aes(size = N), color = "dodgerblue4", alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "twodash", color = "grey50", linewidth = 0.8) +
  geom_text_repel(aes(label = language_full), size = 3, max.overlaps = 15) +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    x        = "True Slope Coefficient",
    y        = "Private Slope Estimate (Epsilon = 8)",
    size     = "Users",
    subtitle = sprintf("Correlation Coefficient = %.2f", cor_scatter)
  ) +
  theme_classic() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "right", plot.subtitle = element_text(hjust = 0.5))

ggsave("Figures/Twitch_Scatter_Eps8_Run1.pdf", p_scatter, width = 5, height = 5)
ggsave("Figures/Twitch_Scatter_Eps8_Run1.png", p_scatter, width = 5, height = 5)
cat("Scatter plot saved to Figures/Twitch_Scatter_Eps8_Run1.pdf and .png\n")
