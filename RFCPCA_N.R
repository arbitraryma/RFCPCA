find_optimal_delta <- function(rec_err_mat, jump_threshold = 0.4, verbose = TRUE) {
  n <- nrow(rec_err_mat)
  
  # Step 1: Compute δ_M
  mean_vec <- colMeans(rec_err_mat)
  centered <- sweep(rec_err_mat, 2, mean_vec)
  delta_M_sq <- mean(rowSums(centered^2))
  delta_M <- sqrt(delta_M_sq)
  
  # Step 2: Iterate λ and compute proportion p
  lambda <- 1
  lambda_vec <- c()
  p_vec <- c()
  
  repeat {
    delta <- lambda * delta_M
    min_d2 <- apply(rec_err_mat, 1, min)
    p <- mean(min_d2 > delta^2)
    
    lambda_vec <- c(lambda_vec, lambda)
    p_vec <- c(p_vec, p)
    
    if (verbose) cat(sprintf("   λ = %.5f | p = %.3f\n", lambda, p))
    if (p > 0.5 || lambda < 1e-4) break
    lambda <- lambda / 2
  }
  
  trace <- data.frame(lambda = lambda_vec, p = p_vec)
  
  # Step 3: Jump detection
  diffs <- diff(trace$p)
  jump_idx <- which(diffs > jump_threshold)[1]
  if (!is.na(jump_idx) && jump_idx > 1) {
    lambda_opt <- trace$lambda[jump_idx]
  } else {
    lambda_opt <- min(trace$lambda)
    warning("No clear jump detected; using smallest lambda instead.")
  }
  
  delta_opt <- lambda_opt * delta_M
  
  list(
    lambda_opt = lambda_opt,
    delta_opt = delta_opt,
    delta_M = delta_M,
    trace = trace
  )
}








# ======================================================================
#  RFCPCA_N  --  Robust Fuzzy CPCA with Noise‑Cluster (single self‑contained function)
# ----------------------------------------------------------------------
#  * One exported function only:  RFCPCA_N()
#  * Internally contains helper closures for
#      - membership update (Eq. 18)
#      - single run with fixed δ
#      - λ‑halving + Pareto search for optimal δ
#  * Optional grid search over k ∈ 2:6 and m ∈ {1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.5}
#    with parallel execution via {furrr}.
#  * Requires existing helpers in your environment: fcpca(),
#    compute_weighted_cross_cov(), compute_S_cov(), etc.
# ======================================================================
# Finalized RFCPCA_N (Robust Fuzzy CPCA with Noise Cluster) with iterative projection updates

# Finalized RFCPCA_N (Robust Fuzzy CPCA with Noise Cluster) with iterative projection updates

# Finalized RFCPCA_N (Robust Fuzzy CPCA with Noise Cluster) with iterative projection updates

RFCPCA_N <- function(ts, 
                     k = NULL, 
                     m = NULL, 
                     startU = NULL, 
                     conver = 1e-3, 
                     maxit = 1000, 
                     verbose = TRUE, 
                     parallel = TRUE, 
                     seed = TRUE) {
  
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE)) stop("Install 'furrr' for parallel execution.")
    library(furrr)
    plan(multisession)
  }
  
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1 || length(m) > 1
  k_range <- if (is.null(k)) 2:6 else k
  m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5) else m
  param_grid <- expand.grid(k = k_range, m = m_range)
  
  membership_update <- function(rec_err_sq, delta_sq, m) {
    N <- nrow(rec_err_sq)
    k <- ncol(rec_err_sq)
    pow <- 1 / (m - 1)
    U_valid <- matrix(0, N, k)
    u_noise <- numeric(N)
    for (i in seq_len(N)) {
      d_i <- pmax(rec_err_sq[i, ], 1e-12)
      inv_delta_term <- (d_i / delta_sq)^pow
      for (c in seq_len(k)) {
        denom <- sum((d_i[c] / d_i)^pow) + inv_delta_term[c]
        denom <- ifelse(!is.finite(denom) || denom <= 0, .Machine$double.xmin, denom)
        U_valid[i, c] <- 1 / denom
      }
      row_sum <- sum(U_valid[i, ])
      u_noise[i] <- 1 - pmin(row_sum, 1)
    }
    cbind(U_valid, u_noise)
  }
  
  run_RFCPCA_N_once <- function(ts, k, m, ...) {
    fcpca_res <- fcpca(ts, k, m = m, startU = startU, conver = conver, maxit = maxit, verbose = FALSE)
    rec_err_sq <- fcpca_res$reconstruction_error_matrix
    delta_res <- find_optimal_delta(rec_err_sq, verbose = FALSE)
    delta_sq <- delta_res$delta_opt^2
    
    n <- nrow(rec_err_sq)
    sigma <- fcpca_res$sigma
    sigma2 <- fcpca_res$sigma2
    comb1 <- fcpca_res$combined_list
    comb2 <- fcpca_res$combined_list2
    k1 <- fcpca_res$n_components_lag01
    k2 <- fcpca_res$n_components_lag02
    
    U <- membership_update(rec_err_sq, delta_sq, m)
    projAx <- fcpca_res$projection_axes
    prev_obj <- Inf
    iteration <- 0
    diffU <- Inf
    
    repeat {
      iteration <- iteration + 1
      
      weighted_covariances <- compute_weighted_cross_cov(U[, 1:k], m, sigma, sigma2)
      
      projAx_new <- vector("list", k)
      for (s in 1:k) {
        svd1 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[s]])
        svd2 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[s]])
        U1 <- svd1$u[, 1:k1, drop = FALSE]
        U2 <- svd2$u[, 1:k2, drop = FALSE]
        projAx_new[[s]] <- list(lag01 = U1, lag02 = U2)
      }
      
      rec_error_new <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (s in 1:k) {
          rec_tt1 <- comb1[[i]] %*% projAx_new[[s]]$lag01 %*% t(projAx_new[[s]]$lag01)
          rec_tt2 <- comb2[[i]] %*% projAx_new[[s]]$lag02 %*% t(projAx_new[[s]]$lag02)
          err_tt1 <- norm(comb1[[i]] - rec_tt1, "F")^2
          err_tt2 <- norm(comb2[[i]] - rec_tt2, "F")^2
          rec_error_new[i, s] <- err_tt1 + err_tt2
        }
      }
      
      U_new <- membership_update(rec_error_new, delta_sq, m)
      diffU <- max(abs(U - U_new))
      U <- U_new
      projAx <- projAx_new
      rec_err_sq <- rec_error_new
      
      obj <- sum(U[, 1:k]^m * rec_err_sq)
      obj_change <- abs(prev_obj - obj)
      prev_obj <- obj
      
      if (verbose) cat(sprintf("Iter %d: Obj = %.4f | dU = %.5f | dObj = %.5f\n", iteration, obj, diffU, obj_change))
      if ((diffU < conver && obj_change < conver) || iteration >= maxit) break
    }
    
    S_val <- compute_S_cov(sum(U[, 1:k]^m * rec_err_sq), weighted_covariances, n, k)
    
    list(
      membership_matrix = U,
      reconstruction_error_matrix = rec_err_sq,
      delta_opt = sqrt(delta_sq),
      lambda_opt = delta_res$lambda_opt,
      k = k,
      m = m,
      S_value = S_val,
      projection_axes = projAx,
      weighted_covariances = weighted_covariances,
      hard_cluster = max.col(U[, 1:k, drop = FALSE]),
      trace = delta_res$trace
    )
  }
  
  if (auto_mode) {
    run_func <- if (parallel) furrr::future_map else purrr::map
    results_list <- run_func(1:nrow(param_grid), function(idx) {
      run_RFCPCA_N_once(ts, k = param_grid$k[idx], m = param_grid$m[idx])
    }, .options = furrr::furrr_options(seed = seed))
    
    S_vals <- sapply(results_list, function(x) x$S_value)
    best_idx <- which.min(S_vals)
    final_result <- results_list[[best_idx]]
    final_result$optimal_k <- param_grid$k[best_idx]
    final_result$optimal_m <- param_grid$m[best_idx]
    final_result$all_results <- results_list
  } else {
    final_result <- run_RFCPCA_N_once(ts, k, m)
    final_result$optimal_k <- k
    final_result$optimal_m <- m
    final_result$all_results <- list(final_result)
  }
  
  return(final_result)
}
