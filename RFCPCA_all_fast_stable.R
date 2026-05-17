###############################################################################
# RFCPCA ALL-IN-ONE: FCPCA + RFCPCA-E + RFCPCA-N + RFCPCA-T
# Stable, optimized, and faster-converging implementation
###############################################################################

suppressPackageStartupMessages({
  library(stats)
})

.get_parallel_map <- function(parallel = TRUE, seed = TRUE) {
  if (parallel && requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE)) {
    future::plan(future::multisession)
    return(function(x, f) furrr::future_map(x, f, .options = furrr::furrr_options(seed = seed)))
  }
  if (parallel) warning("Packages 'future' and/or 'furrr' unavailable; using sequential execution.")
  function(x, f) lapply(x, f)
}

standardize_mts_fast <- function(mts_data) {
  lapply(mts_data, function(x) {
    x <- as.matrix(x)
    mu <- colMeans(x, na.rm = TRUE)
    sig <- apply(x, 2, sd, na.rm = TRUE)
    sig[!is.finite(sig) | sig == 0] <- 1
    sweep(sweep(x, 2, mu, "-"), 2, sig, "/")
  })
}

cross_covariance_block <- function(x, lag = 1L) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n <= lag) stop("Time series must have more rows than lag.")
  gamma0 <- stats::cov(x)
  x_t <- x[seq_len(n - lag), , drop = FALSE]
  x_lag <- x[(lag + 1L):n, , drop = FALSE]
  gamma_lag <- stats::cov(x_t, x_lag)
  block <- rbind(cbind(gamma0, gamma_lag), cbind(t(gamma_lag), gamma0))
  combined <- cbind(x_t, x_lag)
  list(block = block, combined = combined, norm2 = sum(combined * combined))
}

prepare_fcpca_data <- function(ts) {
  ts <- standardize_mts_fast(ts)
  lag1 <- lapply(ts, cross_covariance_block, lag = 1L)
  lag2 <- lapply(ts, cross_covariance_block, lag = 2L)
  list(
    ts = ts,
    sigma1 = lapply(lag1, `[[`, "block"),
    sigma2 = lapply(lag2, `[[`, "block"),
    combined1 = lapply(lag1, `[[`, "combined"),
    combined2 = lapply(lag2, `[[`, "combined"),
    norm1 = vapply(lag1, `[[`, numeric(1), "norm2"),
    norm2 = vapply(lag2, `[[`, numeric(1), "norm2")
  )
}

weighted_covariances_fast <- function(U, m, sigma1, sigma2) {
  W <- U^m
  totals <- colSums(W)
  if (any(!is.finite(totals)) || any(totals <= 1e-12)) {
    stop("At least one cluster has nearly zero total membership weight.")
  }
  d <- nrow(sigma1[[1]])
  flat1 <- vapply(sigma1, as.vector, numeric(d * d))
  flat2 <- vapply(sigma2, as.vector, numeric(d * d))
  cov1 <- flat1 %*% W
  cov2 <- flat2 %*% W
  cov1 <- sweep(cov1, 2, totals, "/")
  cov2 <- sweep(cov2, 2, totals, "/")
  list(
    weighted_cross_cov_lag0_lag1 = lapply(seq_len(ncol(U)), function(j) matrix(cov1[, j], d, d)),
    weighted_cross_cov_lag0_lag2 = lapply(seq_len(ncol(U)), function(j) matrix(cov2[, j], d, d))
  )
}

choose_svd_axes <- function(mat, variance_threshold = 0.95, n_components = NULL) {
  sv <- svd(mat)
  if (is.null(n_components)) {
    denom <- sum(sv$d)
    if (!is.finite(denom) || denom <= 0) {
      n_components <- 1L
    } else {
      n_components <- which(cumsum(sv$d) / denom >= variance_threshold)[1]
      if (is.na(n_components)) n_components <- length(sv$d)
    }
  }
  n_components <- max(1L, min(as.integer(n_components), ncol(sv$u)))
  sv$u[, seq_len(n_components), drop = FALSE]
}

projection_axes_fast <- function(weighted_covariances, variance_threshold = 0.95,
                                 fixed_k1 = NULL, fixed_k2 = NULL) {
  k <- length(weighted_covariances$weighted_cross_cov_lag0_lag1)
  axes <- vector("list", k)
  for (cluster in seq_len(k)) {
    lag01 <- choose_svd_axes(weighted_covariances$weighted_cross_cov_lag0_lag1[[cluster]],
                             variance_threshold = variance_threshold, n_components = fixed_k1)
    lag02 <- choose_svd_axes(weighted_covariances$weighted_cross_cov_lag0_lag2[[cluster]],
                             variance_threshold = variance_threshold, n_components = fixed_k2)
    axes[[cluster]] <- list(lag01 = lag01, lag02 = lag02, k1 = ncol(lag01), k2 = ncol(lag02))
  }
  axes
}

reconstruction_error_fast <- function(prep, projection_axes) {
  n <- length(prep$combined1)
  k <- length(projection_axes)
  err <- matrix(0, nrow = n, ncol = k)
  for (cluster in seq_len(k)) {
    C1 <- projection_axes[[cluster]]$lag01
    C2 <- projection_axes[[cluster]]$lag02
    err1 <- vapply(seq_len(n), function(i) {
      z <- prep$combined1[[i]] %*% C1
      max(prep$norm1[i] - sum(z * z), 0)
    }, numeric(1))
    err2 <- vapply(seq_len(n), function(i) {
      z <- prep$combined2[[i]] %*% C2
      max(prep$norm2[i] - sum(z * z), 0)
    }, numeric(1))
    err[, cluster] <- err1 + err2
  }
  err
}

update_membership_fast <- function(distance, m, eps = 1e-12) {
  n <- nrow(distance); k <- ncol(distance)
  U <- matrix(0, nrow = n, ncol = k)
  zero_rows <- apply(distance <= eps, 1, any)
  if (any(zero_rows)) {
    for (i in which(zero_rows)) {
      z <- which(distance[i, ] <= eps)
      U[i, z] <- 1 / length(z)
    }
  }
  regular <- !zero_rows
  if (any(regular)) {
    power <- -1 / (m - 1)
    D <- pmax(distance[regular, , drop = FALSE], eps)
    numerator <- D^power
    U[regular, ] <- numerator / rowSums(numerator)
  }
  U
}

damp_and_normalize <- function(U_old, U_raw, damping = 0.8) {
  U_new <- damping * U_raw + (1 - damping) * U_old
  rs <- rowSums(U_new)
  rs[!is.finite(rs) | rs <= 0] <- 1
  U_new / rs
}

random_membership <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  U <- matrix(stats::runif(n * k), nrow = n)
  U / rowSums(U)
}

relative_change <- function(old, new) {
  if (!is.finite(old)) return(Inf)
  abs(old - new) / (abs(old) + 1e-8)
}

compute_S_cov_fast <- function(reconstruction_error, weighted_covariances, n, k) {
  cov_dist <- function(a, b) {
    d1 <- sqrt(sum((weighted_covariances$weighted_cross_cov_lag0_lag1[[a]] -
                      weighted_covariances$weighted_cross_cov_lag0_lag1[[b]])^2))
    d2 <- sqrt(sum((weighted_covariances$weighted_cross_cov_lag0_lag2[[a]] -
                      weighted_covariances$weighted_cross_cov_lag0_lag2[[b]])^2))
    d1 + d2
  }
  min_dist <- Inf
  if (k > 1L) {
    for (a in seq_len(k - 1L)) for (b in (a + 1L):k) min_dist <- min(min_dist, cov_dist(a, b))
  }
  reconstruction_error / (n * max(min_dist, .Machine$double.eps))
}

fcpca_fast <- function(ts, k, m = 1.5, startU = NULL, prep = NULL,
                       conver = 1e-3, rel_conver = 1e-7, maxit = 50,
                       verbose = TRUE, replicates = 1, variance_threshold = 0.95,
                       seed = NULL, damping = 0.8, fix_components = TRUE) {
  if (m <= 1) stop("m must be > 1.")
  if (damping <= 0 || damping > 1) stop("damping must be in (0, 1].")
  if (is.null(prep)) prep <- prepare_fcpca_data(ts)
  n <- length(prep$ts)
  best_result <- NULL; best_error <- Inf
  
  for (replicate_i in seq_len(replicates)) {
    if (verbose) cat(sprintf("\nReplicate %d:\n", replicate_i))
    U <- if (is.null(startU)) random_membership(n, k, seed = if (is.null(seed)) NULL else seed + replicate_i) else startU
    iteration <- 0L; current_error <- Inf; prev_error <- Inf; rel_obj_change <- Inf
    fixed_k1 <- NULL; fixed_k2 <- NULL; axes <- NULL; rec_error <- NULL
    
    repeat {
      iteration <- iteration + 1L
      weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
      if (fix_components && iteration == 1L) {
        axes <- projection_axes_fast(weighted_cov, variance_threshold = variance_threshold)
        fixed_k1 <- axes[[1]]$k1; fixed_k2 <- axes[[1]]$k2
      } else if (fix_components) {
        axes <- projection_axes_fast(weighted_cov, variance_threshold = variance_threshold, fixed_k1 = fixed_k1, fixed_k2 = fixed_k2)
      } else {
        axes <- projection_axes_fast(weighted_cov, variance_threshold = variance_threshold)
      }
      rec_error <- reconstruction_error_fast(prep, axes)
      U_raw <- update_membership_fast(rec_error, m)
      U_new <- damp_and_normalize(U, U_raw, damping = damping)
      diffU <- max(abs(U - U_new))
      prev_error <- current_error
      U <- U_new
      current_error <- sum((U^m) * rec_error)
      rel_obj_change <- relative_change(prev_error, current_error)
      if (verbose) cat(sprintf("Iteration %d: Obj = %.4f, dU = %.4g, relObj = %.4g, k1 = %d, k2 = %d\n",
                               iteration, current_error, diffU, rel_obj_change,
                               axes[[1]]$k1, axes[[1]]$k2))
      if ((diffU <= conver && rel_obj_change <= rel_conver) || iteration >= maxit) break
    }
    
    if (current_error < best_error) {
      best_error <- current_error
      final_weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
      best_result <- list(
        method = "FCPCA", membership_matrix = U, projection_axes = axes,
        iterations = iteration, converged = iteration < maxit,
        reconstruction_error = current_error, reconstruction_error_matrix = rec_error,
        hard_cluster = max.col(U), weighted_covariances = final_weighted_cov,
        S_value = compute_S_cov_fast(current_error, final_weighted_cov, n, k),
        n_components_lag01 = axes[[1]]$k1, n_components_lag02 = axes[[1]]$k2,
        sigma = prep$sigma1, sigma2 = prep$sigma2,
        combined_list = prep$combined1, combined_list2 = prep$combined2,
        norm1 = prep$norm1, norm2 = prep$norm2, prep = prep,
        k = k, m = m, damping = damping, conver = conver, rel_conver = rel_conver
      )
    }
  }
  best_result
}

rfcpca_e_fast_once <- function(ts, k, m, startU = NULL, prep = NULL,
                               conver = 1e-3, rel_conver = 1e-7, maxit = 50,
                               verbose = TRUE, variance_threshold = 0.95, seed = NULL,
                               damping = 0.8, fix_components = TRUE) {
  if (is.null(prep)) prep <- prepare_fcpca_data(ts)
  initial <- fcpca_fast(ts, k, m, startU, prep, conver, rel_conver, maxit,
                        verbose = FALSE, replicates = 1, variance_threshold, seed, damping, fix_components)
  U <- initial$membership_matrix
  rec_error <- initial$reconstruction_error_matrix
  beta <- 1 / (mean(apply(rec_error, 1, min)) + 1e-8)
  fixed_k1 <- initial$n_components_lag01; fixed_k2 <- initial$n_components_lag02
  iteration <- 0L; objective <- Inf; prev_obj <- Inf; rel_obj_change <- Inf; axes <- initial$projection_axes
  repeat {
    iteration <- iteration + 1L
    D <- pmax(1 - exp(-beta * rec_error), .Machine$double.eps)
    U_raw <- update_membership_fast(D, m)
    U_new <- damp_and_normalize(U, U_raw, damping)
    diffU <- max(abs(U - U_new)); U <- U_new
    weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
    axes <- projection_axes_fast(weighted_cov, variance_threshold, fixed_k1, fixed_k2)
    rec_error_new <- reconstruction_error_fast(prep, axes)
    D_new <- pmax(1 - exp(-beta * rec_error_new), .Machine$double.eps)
    prev_obj <- objective; objective <- sum((U^m) * D_new); rel_obj_change <- relative_change(prev_obj, objective)
    if (verbose) cat(sprintf("Iteration %d (RFCPCA-E): J_exp = %.6f, dU = %.6g, relObj = %.6g\n", iteration, objective, diffU, rel_obj_change))
    rec_error <- rec_error_new
    if ((diffU <= conver && rel_obj_change <= rel_conver) || iteration >= maxit) break
  }
  final_weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
  list(method = "RFCPCA-E", membership_matrix = U, projection_axes = axes,
       reconstruction_error_matrix = rec_error, weighted_covariances = final_weighted_cov,
       beta = beta, S_value = compute_S_cov_fast(sum((U^m) * rec_error), final_weighted_cov, nrow(U), ncol(U)),
       hard_cluster = max.col(U), iterations = iteration, converged = iteration < maxit,
       k = k, m = m, optimal_k = k, optimal_m = m, damping = damping,
       conver = conver, rel_conver = rel_conver)
}

find_optimal_delta_fast <- function(rec_err_mat, jump_threshold = 0.3, verbose = FALSE, lambda_min = 1e-4) {
  mean_vec <- colMeans(rec_err_mat)
  centered <- sweep(rec_err_mat, 2, mean_vec)
  delta_M <- sqrt(mean(rowSums(centered^2)))
  if (!is.finite(delta_M) || delta_M <= 0) delta_M <- sqrt(mean(apply(rec_err_mat, 1, min))) + 1e-8
  lambda <- 1; lambda_vec <- numeric(0); p_vec <- numeric(0); min_d2 <- apply(rec_err_mat, 1, min)
  repeat {
    delta <- lambda * delta_M
    p <- mean(min_d2 > delta^2)
    lambda_vec <- c(lambda_vec, lambda); p_vec <- c(p_vec, p)
    if (verbose) cat(sprintf("lambda = %.5f | p = %.3f\n", lambda, p))
    if (p > 0.5 || lambda < lambda_min) break
    lambda <- lambda / 2
  }
  trace <- data.frame(lambda = lambda_vec, p = p_vec)
  diffs <- diff(trace$p); jump_idx <- which(diffs > jump_threshold)[1]
  lambda_opt <- if (!is.na(jump_idx) && jump_idx > 1) trace$lambda[jump_idx] else min(trace$lambda)
  list(lambda_opt = lambda_opt, delta_opt = lambda_opt * delta_M, delta_M = delta_M, trace = trace)
}

noise_membership_update_fast <- function(rec_err_sq, delta_sq, m, eps = 1e-12) {
  D <- pmax(rec_err_sq, eps)
  n <- nrow(D); k <- ncol(D); pow <- 1 / (m - 1)
  U_valid <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    d_i <- D[i, ]
    for (c in seq_len(k)) {
      denom <- sum((d_i[c] / d_i)^pow) + (d_i[c] / delta_sq)^pow
      U_valid[i, c] <- 1 / max(denom, .Machine$double.xmin)
    }
  }
  u_noise <- pmax(0, 1 - rowSums(U_valid))
  cbind(U_valid, u_noise)
}

rfcpca_n_fast_once <- function(ts, k, m, startU = NULL, prep = NULL,
                               conver = 1e-3, rel_conver = 1e-7, maxit = 50,
                               verbose = TRUE, variance_threshold = 0.95, seed = NULL,
                               damping = 0.8, fix_components = TRUE, jump_threshold = 0.3) {
  if (is.null(prep)) prep <- prepare_fcpca_data(ts)
  initial <- fcpca_fast(ts, k, m, startU, prep, conver, rel_conver, maxit,
                        verbose = FALSE, replicates = 1, variance_threshold, seed, damping, fix_components)
  rec_error <- initial$reconstruction_error_matrix
  delta_res <- find_optimal_delta_fast(rec_error, jump_threshold = jump_threshold, verbose = FALSE)
  delta_sq <- max(delta_res$delta_opt^2, 1e-12)
  U <- noise_membership_update_fast(rec_error, delta_sq, m)
  U_valid <- U[, seq_len(k), drop = FALSE]
  fixed_k1 <- initial$n_components_lag01; fixed_k2 <- initial$n_components_lag02; axes <- initial$projection_axes
  iteration <- 0L; obj <- Inf; prev_obj <- Inf; rel_obj_change <- Inf
  repeat {
    iteration <- iteration + 1L
    weighted_cov <- weighted_covariances_fast(U_valid, m, prep$sigma1, prep$sigma2)
    axes <- projection_axes_fast(weighted_cov, variance_threshold, fixed_k1, fixed_k2)
    rec_error_new <- reconstruction_error_fast(prep, axes)
    U_raw <- noise_membership_update_fast(rec_error_new, delta_sq, m)
    U_raw_valid <- U_raw[, seq_len(k), drop = FALSE]
    U_valid_new <- damping * U_raw_valid + (1 - damping) * U_valid
    U_noise_new <- pmax(0, 1 - rowSums(U_valid_new))
    U_new <- cbind(U_valid_new, U_noise_new)
    diffU <- max(abs(U - U_new)); U <- U_new; U_valid <- U[, seq_len(k), drop = FALSE]
    prev_obj <- obj
    obj <- sum((U_valid^m) * rec_error_new) + sum((U[, k + 1]^m) * delta_sq)
    rel_obj_change <- relative_change(prev_obj, obj)
    if (verbose) cat(sprintf("Iteration %d (RFCPCA-N): Obj = %.6f, dU = %.6g, relObj = %.6g, delta = %.6g\n", iteration, obj, diffU, rel_obj_change, sqrt(delta_sq)))
    rec_error <- rec_error_new
    if ((diffU <= conver && rel_obj_change <= rel_conver) || iteration >= maxit) break
  }
  final_weighted_cov <- weighted_covariances_fast(U_valid, m, prep$sigma1, prep$sigma2)
  list(method = "RFCPCA-N", membership_matrix = U, valid_membership_matrix = U_valid,
       noise_membership = U[, k + 1], projection_axes = axes, reconstruction_error_matrix = rec_error,
       weighted_covariances = final_weighted_cov, delta_opt = sqrt(delta_sq), lambda_opt = delta_res$lambda_opt,
       trace = delta_res$trace, S_value = compute_S_cov_fast(sum((U_valid^m) * rec_error), final_weighted_cov, nrow(U), k),
       hard_cluster = max.col(U_valid), iterations = iteration, converged = iteration < maxit,
       k = k, m = m, optimal_k = k, optimal_m = m, damping = damping, conver = conver, rel_conver = rel_conver)
}

identify_trim_set_fast <- function(err_mat, U, m, alpha, rule = c("h", "e"), eps = 1e-12) {
  rule <- match.arg(rule); n <- nrow(err_mat); H <- max(1L, floor((1 - alpha) * n))
  if (rule == "h") {
    pow <- 1 / (1 - m)
    h <- rowSums((pmax(err_mat, eps))^pow)^(1 - m)
    ord <- order(h)
  } else {
    e <- rowSums((U^m) * pmax(err_mat, eps))
    ord <- order(e)
  }
  keep_id <- ord[seq_len(H)]
  trim <- rep(TRUE, n); trim[keep_id] <- FALSE
  trim
}

renormalize_trimmed_U <- function(U, trim_idx) {
  U_new <- U; U_new[trim_idx, ] <- 0
  rs <- rowSums(U_new); rs[!is.finite(rs) | rs <= 0] <- 1
  U_new / rs
}

trimmed_membership_update_fast <- function(rec_error, m, trim_idx, eps = 1e-12) {
  U_raw <- update_membership_fast(rec_error, m, eps = eps)
  renormalize_trimmed_U(U_raw, trim_idx)
}

rfcpca_t_fast_once <- function(ts, k, m, alpha, startU = NULL, prep = NULL,
                               conver = 1e-3, rel_conver = 1e-7, maxit = 50,
                               verbose = TRUE, variance_threshold = 0.95, seed = NULL,
                               damping = 0.8, fix_components = TRUE, trim_rule = c("h", "e")) {
  trim_rule <- match.arg(trim_rule)
  if (alpha < 0 || alpha >= 1) stop("alpha must be in [0, 1).")
  if (is.null(prep)) prep <- prepare_fcpca_data(ts)
  initial <- fcpca_fast(ts, k, m, startU, prep, conver, rel_conver, maxit,
                        verbose = FALSE, replicates = 1, variance_threshold, seed, damping, fix_components)
  U <- initial$membership_matrix; rec_error <- initial$reconstruction_error_matrix
  trim_idx <- identify_trim_set_fast(rec_error, U, m, alpha, rule = trim_rule)
  U <- renormalize_trimmed_U(U, trim_idx)
  fixed_k1 <- initial$n_components_lag01; fixed_k2 <- initial$n_components_lag02; axes <- initial$projection_axes
  iteration <- 0L; obj <- Inf; prev_obj <- Inf; rel_obj_change <- Inf
  repeat {
    iteration <- iteration + 1L
    weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
    axes <- projection_axes_fast(weighted_cov, variance_threshold, fixed_k1, fixed_k2)
    rec_error_new <- reconstruction_error_fast(prep, axes)
    trim_idx_new <- identify_trim_set_fast(rec_error_new, U, m, alpha, rule = trim_rule)
    U_raw <- trimmed_membership_update_fast(rec_error_new, m, trim_idx_new)
    U_damped <- damping * U_raw + (1 - damping) * U
    U_damped[trim_idx_new, ] <- 0
    U_new <- renormalize_trimmed_U(U_damped, trim_idx_new)
    diffU <- max(abs(U - U_new)); U <- U_new; trim_idx <- trim_idx_new
    prev_obj <- obj; obj <- sum((U^m) * rec_error_new); rel_obj_change <- relative_change(prev_obj, obj)
    if (verbose) cat(sprintf("Iteration %d (RFCPCA-T): Obj = %.6f, dU = %.6g, relObj = %.6g, keep = %d, alpha = %.2f\n", iteration, obj, diffU, rel_obj_change, sum(!trim_idx), alpha))
    rec_error <- rec_error_new
    if ((diffU <= conver && rel_obj_change <= rel_conver) || iteration >= maxit) break
  }
  final_weighted_cov <- weighted_covariances_fast(U, m, prep$sigma1, prep$sigma2)
  n_keep <- sum(!trim_idx)
  list(method = "RFCPCA-T", membership_matrix = U, projection_axes = axes,
       reconstruction_error_matrix = rec_error, weighted_covariances = final_weighted_cov,
       trim_set = trim_idx, kept_set = !trim_idx, alpha_used = alpha,
       S_value = compute_S_cov_fast(sum((U^m) * rec_error), final_weighted_cov, max(n_keep, 1L), k),
       hard_cluster = max.col(U), iterations = iteration, converged = iteration < maxit,
       k = k, m = m, alpha = alpha, optimal_k = k, optimal_m = m, optimal_alpha = alpha,
       damping = damping, conver = conver, rel_conver = rel_conver)
}

# Auto-selection wrappers ------------------------------------------------------

RFCPCA_E <- function(ts, k = NULL, m = NULL, startU = NULL, conver = 1e-3,
                     rel_conver = 1e-7, maxit = 50, verbose = TRUE, parallel = TRUE,
                     variance_threshold = 0.95, seed = TRUE, damping = 0.8,
                     fix_components = TRUE) {
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1L || length(m) > 1L
  k_range <- if (is.null(k)) 2:6 else k
  m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0) else m
  grid <- expand.grid(k = k_range, m = m_range)
  prep <- prepare_fcpca_data(ts); pmap <- .get_parallel_map(parallel, seed)
  run_one <- function(idx) {
    if (verbose) cat(sprintf("Running RFCPCA-E: k = %d, m = %.2f\n", grid$k[idx], grid$m[idx]))
    rfcpca_e_fast_once(ts, grid$k[idx], grid$m[idx], startU, prep, conver, rel_conver, maxit,
                       FALSE, variance_threshold, if (isTRUE(seed)) idx else if (is.numeric(seed)) seed + idx else NULL,
                       damping, fix_components)
  }
  results <- pmap(seq_len(nrow(grid)), run_one)
  if (!auto_mode) { final <- results[[1]]; final$all_results <- results; return(final) }
  S <- vapply(results, `[[`, numeric(1), "S_value"); best <- which.min(S)
  final <- results[[best]]; final$optimal_k <- grid$k[best]; final$optimal_m <- grid$m[best]
  final$all_results <- results; final$search_grid <- grid; final
}

RFCPCA_N <- function(ts, k = NULL, m = NULL, startU = NULL, conver = 1e-3,
                     rel_conver = 1e-7, maxit = 50, verbose = TRUE, parallel = TRUE,
                     variance_threshold = 0.95, seed = TRUE, damping = 0.8,
                     fix_components = TRUE, jump_threshold = 0.3) {
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1L || length(m) > 1L
  k_range <- if (is.null(k)) 2:6 else k
  m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0) else m
  grid <- expand.grid(k = k_range, m = m_range)
  prep <- prepare_fcpca_data(ts); pmap <- .get_parallel_map(parallel, seed)
  run_one <- function(idx) {
    if (verbose) cat(sprintf("Running RFCPCA-N: k = %d, m = %.2f\n", grid$k[idx], grid$m[idx]))
    rfcpca_n_fast_once(ts, grid$k[idx], grid$m[idx], startU, prep, conver, rel_conver, maxit,
                       FALSE, variance_threshold, if (isTRUE(seed)) idx else if (is.numeric(seed)) seed + idx else NULL,
                       damping, fix_components, jump_threshold)
  }
  results <- pmap(seq_len(nrow(grid)), run_one)
  if (!auto_mode) { final <- results[[1]]; final$all_results <- results; return(final) }
  S <- vapply(results, `[[`, numeric(1), "S_value"); best <- which.min(S)
  final <- results[[best]]; final$optimal_k <- grid$k[best]; final$optimal_m <- grid$m[best]
  final$all_results <- results; final$search_grid <- grid; final
}

RFCPCA_T <- function(ts, k = NULL, m = NULL, alpha = NULL, startU = NULL,
                     conver = 1e-3, rel_conver = 1e-7, maxit = 50, verbose = TRUE,
                     parallel = TRUE, variance_threshold = 0.95, seed = TRUE,
                     damping = 0.8, fix_components = TRUE, trim_rule = c("h", "e")) {
  trim_rule <- match.arg(trim_rule)
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1L || length(m) > 1L || is.null(alpha) || length(alpha) > 1L
  k_range <- if (is.null(k)) 2:6 else k
  m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0) else m
  alpha_range <- if (is.null(alpha)) c(0.1, 0.2, 0.3, 0.4, 0.5) else alpha
  grid <- expand.grid(k = k_range, m = m_range, alpha = alpha_range)
  prep <- prepare_fcpca_data(ts); pmap <- .get_parallel_map(parallel, seed)
  run_one <- function(idx) {
    if (verbose) cat(sprintf("Running RFCPCA-T: k = %d, m = %.2f, alpha = %.2f\n", grid$k[idx], grid$m[idx], grid$alpha[idx]))
    rfcpca_t_fast_once(ts, grid$k[idx], grid$m[idx], grid$alpha[idx], startU, prep, conver, rel_conver, maxit,
                       FALSE, variance_threshold, if (isTRUE(seed)) idx else if (is.numeric(seed)) seed + idx else NULL,
                       damping, fix_components, trim_rule)
  }
  results <- pmap(seq_len(nrow(grid)), run_one)
  if (!auto_mode) { final <- results[[1]]; final$all_results <- results; return(final) }
  S <- vapply(results, `[[`, numeric(1), "S_value"); best <- which.min(S)
  final <- results[[best]]; final$optimal_k <- grid$k[best]; final$optimal_m <- grid$m[best]; final$optimal_alpha <- grid$alpha[best]
  final$all_results <- results; final$search_grid <- grid; final
}

RFCPCA_ALL <- function(ts, k = 2, m = 1.6, alpha = 0.2, conver = 1e-3,
                       rel_conver = 1e-7, maxit = 50, damping = 0.8,
                       parallel = FALSE, seed = 123, variance_threshold = 0.95,
                       trim_rule = "h", run_E = TRUE, run_N = TRUE, run_T = TRUE) {
  out <- list()
  if (run_E) out$RFCPCA_E <- RFCPCA_E(ts, k, m, conver = conver, rel_conver = rel_conver, maxit = maxit,
                                      damping = damping, parallel = parallel, seed = seed, variance_threshold = variance_threshold)
  if (run_N) out$RFCPCA_N <- RFCPCA_N(ts, k, m, conver = conver, rel_conver = rel_conver, maxit = maxit,
                                      damping = damping, parallel = parallel, seed = seed, variance_threshold = variance_threshold)
  if (run_T) out$RFCPCA_T <- RFCPCA_T(ts, k, m, alpha, conver = conver, rel_conver = rel_conver, maxit = maxit,
                                      damping = damping, parallel = parallel, seed = seed, variance_threshold = variance_threshold,
                                      trim_rule = trim_rule)
  out
}

time_method <- function(expr) {
  t <- system.time({ value <- eval.parent(substitute(expr)) })
  list(result = value, elapsed_sec = unname(t["elapsed"]))
}

###############################################################################
# Example usage
###############################################################################
# source("RFCPCA_all_fast_stable.R")
#
# out_e <- time_method(RFCPCA_E(EEGsample_subject, k = 2, m = 1.6,
#                              parallel = FALSE, seed = 123))
# out_e$elapsed_sec
# res_e <- out_e$result
#
# out_n <- time_method(RFCPCA_N(EEGsample_subject, k = 2, m = 1.6,
#                              parallel = FALSE, seed = 1))
# res_n <- out_n$result
#
# out_t <- time_method(RFCPCA_T(EEGsample_subject, k = 2, m = 1.6, alpha = 0.2,
#                              parallel = FALSE, seed = 1))
# res_t <- out_t$result
#
# all_res <- RFCPCA_ALL(EEGsample_subject, k = 2, m = 1.6, alpha = 0.2)
