# Shared configuration and output handling. Run scripts from the repository root.
# Defaults use the paper's network sizes and repetition counts.

replication_config <- function() {
  output <- Sys.getenv("DP_OUTPUT_DIR", ".")
  data <- Sys.getenv("DP_DATA_DIR", ".data")
  iteration_input <- Sys.getenv("DP_ITERATIONS", "")
  iterations <- NULL
  if (nzchar(iteration_input)) {
    iterations <- suppressWarnings(as.numeric(iteration_input))
    if (!is.finite(iterations) || iterations < 2 || iterations != floor(iterations) ||
        iterations > .Machine$integer.max) {
      stop("DP_ITERATIONS must be an integer between 2 and ", .Machine$integer.max, ".")
    }
    iterations <- as.integer(iterations)
  }
  list(output = output, data = data,
       iterations = iterations, batch_size = 10L)
}
config <- replication_config()

iteration_count <- function(default) {
  if (!is.null(config$iterations)) return(config$iterations)
  default
}

output_path <- function(directory, filename) {
  path <- file.path(config$output, directory)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  file.path(path, filename)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". See the package installation command in README.md.")
}

start_experiment <- function(name, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  # Avoid all-core defaults and large allocations on laptops.
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(2)
  cat(sprintf("\n%s: output=%s, seed=%d\n", name, config$output, seed))
}

check_data_files <- function(paths) {
  bad <- paths[!file.exists(paths) | is.na(file.info(paths)$size) | file.info(paths)$size == 0]
  if (length(bad)) stop("Missing or empty data files:\n", paste(bad, collapse = "\n"),
                       "\nSet DP_DATA_DIR to the downloaded data directory.")
}

save_plot <- function(plot, filename, width = 5, height = 5, png = FALSE) {
  ggplot2::ggsave(output_path("Figures", paste0(filename, ".pdf")), plot,
                 width = width, height = height)
  if (png) ggplot2::ggsave(output_path("Figures", paste0(filename, ".png")), plot,
                          width = width, height = height, dpi = 160)
}

# MSE is explicitly conditional on success when a regression can fail.
accuracy_summary <- function(estimate, truth, status = rep("ok", length(estimate))) {
  valid <- status == "ok" & is.finite(estimate)
  errors <- estimate[valid] - rep_len(truth, length(estimate))[valid]
  data.frame(n_attempted = length(estimate), n_valid = sum(valid),
             n_failed = sum(!valid), failure_rate = mean(!valid),
             mse = if (any(valid)) mean(errors^2) else NA_real_,
             bias = if (any(valid)) mean(errors) else NA_real_,
             sd = if (sum(valid) > 1) sd(estimate[valid]) else NA_real_)
}

# Correlate the paired cells supplied by the caller; this helper does not omit
# nonfinite values. Callers define the analysis sample: villages supply the full
# sample, while Twitch selects successful, finite estimates within each draw.
complete_correlation <- function(x, y) {
  if (length(x) < 3 || any(!is.finite(c(x, y))) || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}

# Privacy mechanisms used by the replication scripts.
# Graphs are undirected, nonnegative, loop-free, and have a fixed vertex set.
# Labels, ranks, cell membership, and sample selection enter as specified inputs;
# these functions do not privatize preprocessing decisions. See README.md.

positive_scalar <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(name, " must be a finite positive number.", call. = FALSE)
  }
}

# Inverse-CDF sampler: scale is the Laplace scale, not its standard deviation.
r_laplace <- function(n, scale = 1) {
  if (any(!is.finite(scale)) || any(scale < 0)) stop("Invalid Laplace scale.")
  u <- runif(n, -0.5, 0.5)
  -scale * sign(u) * log1p(-2 * abs(u))
}

validate_adjacency <- function(adjacency) {
  if (nrow(adjacency) != ncol(adjacency) || nrow(adjacency) < 1L ||
      any(!is.finite(adjacency)) || any(adjacency < 0) ||
      any(diag(adjacency) != 0) || !isSymmetric(unname(adjacency))) {
    stop("Adjacency must be finite, nonnegative, symmetric, square, and loop-free.")
  }
  invisible(adjacency)
}

# Zero rows stay zero: isolates contribute a zero exposure and remain in the sample.
row_normalize <- function(adjacency) {
  degree <- rowSums(adjacency)
  inverse_degree <- numeric(length(degree))
  inverse_degree[degree > 0] <- 1 / degree[degree > 0]
  if (inherits(adjacency, "Matrix")) {
    Matrix::Diagonal(x = inverse_degree) %*% adjacency
  } else {
    adjacency * inverse_degree
  }
}

true_connectedness <- function(weights, labels) {
  if (length(labels) != nrow(weights) || any(!labels %in% c(0, 1))) {
    stop("Supply one binary label per vertex.")
  }
  if (!any(labels == 0)) return(NA_real_)
  mean(as.numeric(weights %*% labels)[labels == 0])
}

# Algorithm 1, including the isolate correction and B_n floor.
# uniform_flips and standard_laplace support common random numbers in simulations.
# Do not publish these random draws alongside releases from confidential inputs.
# Optional diagnostics depend only on sanitized labels and public parameters;
# the unnoised edge-dependent numerator must remain internal.
private_connectedness <- function(weights, labels, epsilon_label, epsilon_edge,
                                  uniform_flips = runif(length(labels)),
                                  standard_laplace = r_laplace(1), details = FALSE) {
  positive_scalar(epsilon_label, "epsilon_label")
  positive_scalar(epsilon_edge, "epsilon_edge")
  n <- length(labels)
  if (nrow(weights) != n || ncol(weights) != n || any(!labels %in% c(0, 1)) ||
      length(uniform_flips) != n || any(!is.finite(uniform_flips)) ||
      any(uniform_flips < 0 | uniform_flips > 1) ||
      length(standard_laplace) != 1 || !is.finite(standard_laplace)) {
    stop("Invalid binary mechanism inputs.")
  }
  p <- plogis(-epsilon_label)
  q <- tanh(epsilon_label / 2)  # Stable evaluation of 1 - 2p.
  private_labels <- ifelse(uniform_flips < p, 1 - labels, labels)
  membership <- ((private_labels == 0) - p) / q
  # Debias inside the matrix product: this is exactly (rho_hat - p*1{d>0})/q.
  exposure <- as.numeric(weights %*% ((private_labels - p) / q))
  group_size <- sum(membership)
  numerator <- sum(membership * exposure)
  denominator <- max(group_size, 1)
  sensitivity <- 2 * (1 - p) / q^2
  estimate <- (numerator + standard_laplace * sensitivity / epsilon_edge) / denominator
  if (!details) return(estimate)
  list(estimate = estimate, stabilized = group_size < 1, group_size = group_size,
       sensitivity = sensitivity)
}

truncated_laplace_parameters <- function(epsilon_label, delta_label) {
  positive_scalar(epsilon_label, "epsilon_label")
  if (length(delta_label) != 1 || !is.finite(delta_label) ||
      delta_label <= 0 || delta_label >= 0.5) stop("delta_label must be in (0, 0.5).")
  lambda <- 1 / epsilon_label
  # log(1 + (exp(epsilon)-1)/(2*delta)), without overflowing exp(epsilon).
  log_term <- if (epsilon_label < 50) {
    log1p(expm1(epsilon_label) / (2 * delta_label))
  } else {
    epsilon_label + log1p((2 * delta_label - 1) * exp(-epsilon_label)) - log(2 * delta_label)
  }
  A <- lambda * log_term
  list(lambda = lambda, A = A, variance = truncated_laplace_variance(lambda, A),
       bounds = list(a = -A, b = 1 + A, a_prime = -A, b_prime = 1 + A))
}

truncated_laplace_variance <- function(lambda, A) {
  positive_scalar(lambda, "lambda")
  positive_scalar(A, "A")
  t <- A / lambda
  # The series avoids cancellation when truncation is very tight.
  if (t < 1e-3) return(A^2 * (1 / 3 - t / 12 + t^2 / 360 + t^3 / 720))
  lambda^2 * (2 - exp(-t) * (t^2 + 2 * t) / (-expm1(-t)))
}

r_truncated_laplace <- function(n, lambda, A) {
  positive_scalar(lambda, "lambda")
  positive_scalar(A, "A")
  magnitude <- -lambda * log1p(-runif(n) * (-expm1(-A / lambda)))
  magnitude * ifelse(runif(n) < 0.5, -1, 1)
}

# Mean of the fitted line on a publicly chosen interval.
mean_friend_rank <- function(alpha, beta, lower = 0, upper = 0.25) {
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) stop("Invalid rank interval.")
  alpha + beta * (lower + upper) / 2
}

# Manuscript Algorithm 2 edge stage for simple, unweighted networks.
# The labels supplied here must be privatized independently of edges.
# Main simulations and empirical analyses use the full-regression mode below.

# Fix sanitized labels h. One edge changes only two friend averages, each by
# at most K=max(abs(h)), including transitions into/out of isolation. Since
# ncov=sum((h-mean(h))*y), K times the two largest absolute centered labels
# bounds its edge sensitivity. The mean outcome changes by at most 2*K/n.
# These bounds use no unprotected graph information.
conditional_network_sensitivities <- function(sanitized_labels) {
  h <- sanitized_labels
  if (length(h) < 2L || any(!is.finite(h))) {
    stop("Supply at least two finite sanitized labels.", call. = FALSE)
  }
  centered <- abs(h - mean(h))
  largest_index <- which.max(centered)
  largest <- centered[largest_index]
  centered[largest_index] <- 0
  second_largest <- max(centered)
  maximum_label <- max(abs(h))
  c(covariance = maximum_label * (largest + second_largest),
    mean_y = 2 * maximum_label / length(h))
}

# Release a slope adjusted for label-noise variance, optionally with an intercept.
# Conditional on sanitized labels, nvar and the correction denominator have zero
# edge sensitivity. Failure depends only on the sanitized labels as well.
#
# For a full regression, release covariance and mean_y using epsilon_edge/2
# each, then form beta and alpha by post-processing. A slope-only release uses
# all epsilon_edge for covariance and returns alpha=NA intentionally.
#
# A successful release has corrected variance greater than variance_floor, so
# its slope denominator exceeds (n - 1) * variance_floor.
#
# Optional standard noise supports paired public-data experiments and tests;
# never disclose these draws with releases of confidential inputs.
private_conditional_regression_stats <- function(
    n, mean_x, mean_y, nvar, ncov, epsilon_edge,
    covariance_sensitivity, mean_y_sensitivity, noise_variance = 0,
    variance_floor = 1e-8, output = c("regression", "slope"),
    standard_covariance_noise = r_laplace(length(nvar)),
    standard_mean_noise = r_laplace(length(nvar))) {
  output <- match.arg(output)
  positive_scalar(epsilon_edge, "epsilon_edge")
  positive_scalar(variance_floor, "variance_floor")
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != floor(n)) {
    stop("n must be an integer of at least two.", call. = FALSE)
  }
  if (length(noise_variance) != 1L || !is.finite(noise_variance) || noise_variance < 0) {
    stop("Noise variance must be finite and nonnegative.", call. = FALSE)
  }
  values <- list(mean_x, mean_y, nvar, ncov, covariance_sensitivity,
                 mean_y_sensitivity, standard_covariance_noise)
  if (output == "regression") values <- c(values, list(standard_mean_noise))
  if (!length(nvar) || length(unique(lengths(values))) != 1L ||
      any(!is.finite(unlist(values, use.names = FALSE))) ||
      any(nvar < 0) || any(covariance_sensitivity < 0) || any(mean_y_sensitivity < 0)) {
    stop("Supply finite vectors of equal length with nonnegative variances and sensitivities.",
         call. = FALSE)
  }

  corrected_variance <- nvar / (n - 1) - noise_variance
  status <- ifelse(corrected_variance > variance_floor, "ok", "nonpositive_corrected_variance")
  ok <- status == "ok"
  beta <- alpha <- rep(NA_real_, length(nvar))
  epsilon_covariance <- if (output == "regression") epsilon_edge / 2 else epsilon_edge
  noisy_covariance <- ncov + standard_covariance_noise * covariance_sensitivity / epsilon_covariance
  beta[ok] <- noisy_covariance[ok] / ((n - 1) * corrected_variance[ok])
  if (output == "regression") {
    noisy_mean <- mean_y + standard_mean_noise * mean_y_sensitivity / (epsilon_edge / 2)
    alpha[ok] <- noisy_mean[ok] - beta[ok] * mean_x[ok]
  }
  invalid <- ok & (!is.finite(beta) | (output == "regression" & !is.finite(alpha)))
  status[invalid] <- "nonfinite_estimate"
  beta[status != "ok"] <- alpha[status != "ok"] <- NA_real_
  data.frame(beta = beta, alpha = alpha, status = status)
}

# Sparse graph construction keeps large Monte Carlo runs within laptop memory.
graph_weights <- function(graph) {
  row_normalize(igraph::as_adjacency_matrix(graph, sparse = TRUE))
}

generate_sbm_network <- function(num_low, num_high, p_ll, p_lh, p_hh) {
  sizes <- c(num_low, num_high)
  probabilities <- matrix(c(p_ll, p_lh, p_lh, p_hh), 2, 2)
  present <- sizes > 0
  graph <- igraph::sample_sbm(sum(sizes), probabilities[present, present, drop = FALSE],
                              sizes[present], directed = FALSE, loops = FALSE)
  graph_weights(graph)
}

generate_graphon_network <- function(N, expected_degree, homophily) {
  if (N < 2 || homophily < 0 || expected_degree < 0) stop("Invalid graphon parameters.")
  labels <- runif(N)
  # Integral of exp(-h*|x-y|) over the unit square; use a series near zero.
  h <- homophily
  kernel_mean <- if (h < 1e-4) 1 - h / 3 + h^2 / 12 - h^3 / 60 else
    2 / h + 2 * expm1(-h) / h^2
  probability <- expected_degree / ((N - 1) * kernel_mean)
  if (probability > 1) stop("Requested degree cannot use this rejection sampler: base probability > 1.")
  candidate <- igraph::sample_gnp(N, probability, directed = FALSE, loops = FALSE)
  edges <- igraph::as_edgelist(candidate, names = FALSE)
  graph <- igraph::make_empty_graph(N, directed = FALSE)
  if (nrow(edges)) {
    keep <- runif(nrow(edges)) < exp(-h * abs(labels[edges[, 1]] - labels[edges[, 2]]))
    if (any(keep)) graph <- igraph::add_edges(graph, as.vector(t(edges[keep, , drop = FALSE])))
  }
  # Explicit N preserves trailing isolates even when no edge survives thinning.
  list(g = graph, x_raw = labels)
}

# Shared continuous-label simulation. The graph and raw labels are fixed within
# an experiment; each row of the returned table is one complete private release.
fit_line <- function(x, y) {
  centered <- x - mean(x)
  denominator <- sum(centered^2)
  if (denominator <= 0) return(c(alpha = NA_real_, beta = NA_real_, r_squared = NA_real_))
  beta <- sum(centered * (y - mean(y))) / denominator
  total <- sum((y - mean(y))^2)
  c(alpha = mean(y) - beta * mean(x), beta = beta,
    r_squared = if (total > 0) beta^2 * denominator / total else NA_real_)
}

simulate_continuous <- function(weights, x, epsilon_label, epsilon_edge, delta_label,
                                repetitions, batch_size = config$batch_size) {
  if (any(!is.finite(x)) || any(x < 0 | x > 1) || nrow(weights) != length(x)) {
    stop("Continuous labels must be in [0,1], one per vertex.")
  }
  parameters <- truncated_laplace_parameters(epsilon_label, delta_label)
  n <- length(x)
  records <- list()
  for (first in seq.int(1L, repetitions, by = batch_size)) {
    count <- min(batch_size, repetitions - first + 1L)
    x_private <- x + matrix(r_truncated_laplace(n * count, parameters$lambda, parameters$A), nrow = n)
    y_private <- as.matrix(weights %*% x_private)
    x_mean <- colMeans(x_private)
    y_mean <- colMeans(y_private)
    xc <- sweep(x_private, 2, x_mean, "-")
    # Centering y is unnecessary in this cross product because colSums(xc)=0.
    nvar <- colSums(xc^2)
    ncov <- colSums(xc * y_private)
    # Each sensitivity is a function only of the already sanitized label vector.
    sensitivities <- vapply(seq_len(count), function(j) {
      conditional_network_sensitivities(x_private[, j])
    }, c(covariance = 0, mean_y = 0))
    result <- private_conditional_regression_stats(
      n, x_mean, y_mean, nvar, ncov, epsilon_edge,
      sensitivities["covariance", ], sensitivities["mean_y", ], parameters$variance)
    result$iteration <- seq.int(first, length.out = count)
    result$mafr <- mean_friend_rank(result$alpha, result$beta)
    records[[length(records) + 1L]] <- result
  }
  draws <- data.table::as.data.table(do.call(rbind, records))
  draws[, `:=`(epsilon_label = epsilon_label, epsilon_edge = epsilon_edge,
               epsilon = epsilon_label + epsilon_edge, delta_label = delta_label)]
  draws
}

# Discover paired village inputs before starting a potentially expensive run.
village_input_files <- function(data_directory) {
  base <- file.path(data_directory, "datav4.0/datav4.0/Data")
  adjacency_directory <- file.path(base, "1. Network Data/Adjacency Matrices")
  adjacency_files <- list.files(adjacency_directory,
    pattern = "^adj_allVillageRelationships_HH_vilno_[0-9]+[.]csv$", full.names = TRUE)
  if (!length(adjacency_files)) stop("No village adjacency files found in: ", adjacency_directory)
  key_names <- sub("^adj_allVillageRelationships_HH_", "key_HH_", basename(adjacency_files))
  c(file.path(base, "2. Demographics and Outcomes/household_characteristics.dta"),
    adjacency_files, file.path(base, "1. Network Data/Adjacency Matrix Keys", key_names))
}

# Streaming parser for SNAP's Amazon metadata. Product records may cross chunks.
read_amazon_metadata <- function(path, chunk_size = 50000L) {
  connection <- file(path, open = "r")
  on.exit(close(connection))
  current_id <- NA_integer_
  records <- list()
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    id_rows <- grep("^Id:[[:space:]]+[0-9]+", lines)
    rank_rows <- grep("^[[:space:]]+salesrank:[[:space:]]+-?[0-9]+", lines)
    ids <- as.integer(sub("^Id:[[:space:]]+", "", lines[id_rows]))
    if (length(rank_rows)) {
      rank_ids <- c(current_id, ids)[findInterval(rank_rows, id_rows) + 1L]
      ranks <- as.integer(sub("^[[:space:]]+salesrank:[[:space:]]+", "", lines[rank_rows]))
      records[[length(records) + 1L]] <- data.table::data.table(id = rank_ids, salesrank = ranks)
    }
    if (length(ids)) current_id <- tail(ids, 1L)
  }
  result <- data.table::rbindlist(records)
  if (!nrow(result) || anyNA(result$id) || anyDuplicated(result$id)) stop("Malformed Amazon product metadata.")
  result[salesrank > 0]
}

# Build a sparse row-normalized graph on an explicitly supplied vertex universe.
# Duplicate undirected edges are collapsed and isolates are retained.
weights_from_edges <- function(vertices, from, to) {
  if (anyDuplicated(vertices)) stop("Vertex IDs must be unique.")
  i <- match(from, vertices)
  j <- match(to, vertices)
  if (anyNA(c(i, j))) stop("Edge endpoint absent from the specified vertex set.")
  if (any(i == j)) stop("Self-loops are outside this paper's graph model.")
  A <- Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1,
                            dims = c(length(vertices), length(vertices)))
  # Sparse logical replacement can overflow R's 32-bit matrix indices on large
  # vertex sets. sign() collapses duplicates without constructing linear indices.
  A <- sign(A)
  row_normalize(A)
}
