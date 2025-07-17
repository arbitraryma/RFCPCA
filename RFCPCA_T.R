###############################################################################
## Trimmed (Robust) FCPCA Function: RFCPCA_T
###############################################################################
RFCPCA_T <- function(ts,
                     k               = NULL,
                     m               = NULL,
                     alpha           = 0.10,   # trimming proportion
                     startU          = NULL,
                     conver          = 1e-3,
                     maxit           = 1000,
                     verbose         = TRUE,
                     parallel        = TRUE) {
  
  ## ----------------------------- 0)  PARALLEL SET-UP ---------------------- ##
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE))
      stop("Package ‘furrr’ is required for parallel execution.")
    library(furrr); plan(multisession)
  }
  
  ## ----------------------------- 1)  HYPER-GRID --------------------------- ##
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1 || length(m) > 1
  if (auto_mode) {
    k_range <- if (is.null(k)) 2:6 else k
    m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5) else m
  } else {
    k_range <- k; m_range <- m
  }
  param_grid <- expand.grid(k = k_range, m = m_range)
  
  ## ----------------------------- 2)  HELPERS ------------------------------ ##
  identify_trim_set <- function(err_mat, alpha) {
    min_err <- apply(err_mat, 1, min)
    h       <- ceiling((1 - alpha) * length(min_err))
    keep_id <- order(min_err)[1:h]
    trim    <- rep(TRUE, length(min_err))
    trim[keep_id] <- FALSE
    trim
  }
  
  renormalise_U <- function(U, trim_idx) {
    U_trim <- U
    U_trim[trim_idx, ] <- 0
    row_sums <- rowSums(U_trim)
    # avoid division by zero for fully-trimmed rows
    row_sums[row_sums == 0] <- 1
    U_trim <- U_trim / row_sums
    U_trim
  }
  
  run_once <- function(ts, k, m, alpha, startU, conver, maxit, verbose) {
    
    init <- fcpca(ts, k, m = m, startU = startU,
                  conver = conver, maxit = maxit,
                  verbose = FALSE, replicates = 1)
    
    U         <- init$membership_matrix
    projAx    <- init$projection_axes
    sigma     <- init$sigma
    sigma2    <- init$sigma2
    comb1     <- init$combined_list
    comb2     <- init$combined_list2
    k1        <- init$n_components_lag01
    k2        <- init$n_components_lag02
    k_clust   <- ncol(U)
    n         <- nrow(U)
    
    rec_err   <- init$reconstruction_error_matrix
    trim_idx  <- identify_trim_set(rec_err, alpha)
    U         <- renormalise_U(U, trim_idx)
    
    iter <- 0; diffU <- Inf; prev_obj <- Inf
    
    repeat {
      iter <- iter + 1
      
      ## 2a) Weighted covariances on kept objects
      w_cov <- compute_weighted_cross_cov(U, m, sigma, sigma2)
      
      ## 2b) Update axes (k1, k2 fixed from init)
      projAx <- lapply(1:k_clust, function(s) {
        svd1 <- svd(w_cov$weighted_cross_cov_lag0_lag1[[s]])
        svd2 <- svd(w_cov$weighted_cross_cov_lag0_lag2[[s]])
        list(lag01 = svd1$u[, 1:k1, drop = FALSE],
             lag02 = svd2$u[, 1:k2, drop = FALSE])
      })
      
      ## 2c) Re-compute reconstruction error matrix
      rec_err <- matrix(0, nrow = n, ncol = k_clust)
      for (i in 1:n) for (s in 1:k_clust) {
        r1 <- comb1[[i]] %*% projAx[[s]]$lag01 %*% t(projAx[[s]]$lag01)
        r2 <- comb2[[i]] %*% projAx[[s]]$lag02 %*% t(projAx[[s]]$lag02)
        rec_err[i, s] <- norm(comb1[[i]] - r1, "F")^2 +
          norm(comb2[[i]] - r2, "F")^2
      }
      
      ## 2d) Trim step (update trim set every iteration)
      trim_idx_new <- identify_trim_set(rec_err, alpha)
      keep_mask    <- !trim_idx_new
      
      ## 2e) Update memberships only for kept objects
      U_new <- U
      eps    <- 1e-10
      
      for (i in which(keep_mask)) {
        for (s in 1:k_clust) {
          denom <- sum((rec_err[i, s] / (rec_err[i, ] + eps))^(1 / (m - 1)))
          U_new[i, s] <- 1 / denom
        }
      }
      # fully trim rows
      U_new[trim_idx_new, ] <- 0
      # renormalise rows of kept objects
      U_new <- renormalise_U(U_new, trim_idx_new)
      
      
      diffU   <- max(abs(U - U_new))
      U       <- U_new
      trim_idx <- trim_idx_new
      
      ## 2f) objective on kept objects
      obj <- sum((U^m) * rec_err)
      if (verbose)
        cat(sprintf("Iter %3d:  Obj = %.6f   ΔU = %.6f   |keep| = %d\n",
                    iter, obj, diffU, sum(!trim_idx)))
      
      if ((diffU < conver && abs(prev_obj - obj) < conver) || iter >= maxit)
        break
      prev_obj <- obj
    }
    
    ## 2g) final diagnostics
    final_cov <- compute_weighted_cross_cov(U, m, sigma, sigma2)
    # S_val     <- compute_S_cov(sum(U^m * rec_err), final_cov, n, k_clust)
    S_val     <- compute_S_cov(sum(U^m * rec_err), final_cov, sum(!trim_idx), k_clust)
    
    list(
      membership_matrix = U,
      projection_axes   = projAx,
      reconstruction_error_matrix = rec_err,
      weighted_covariances = final_cov,
      trim_set          = trim_idx,
      alpha_used        = alpha,
      S_value           = S_val,
      k = k_clust, m = m
    )
  }
  
  ## ----------------------------- 3)  RUN GRID ----------------------------- ##
  exec <- function(idx) {
    run_once(ts,
             k       = param_grid$k[idx],
             m       = param_grid$m[idx],
             alpha   = alpha,
             startU  = startU,
             conver  = conver,
             maxit   = maxit,
             verbose = FALSE)
  }
  
  results <- if (auto_mode) {
    if (parallel) furrr::future_map(1:nrow(param_grid), exec,
                                    .options = furrr_options(seed = TRUE))
    else           lapply(1:nrow(param_grid), exec)
  } else {
    list(exec(1))
  }
  
  ## ----------------------------- 4)  SELECT BEST -------------------------- ##
  S_vec   <- sapply(results, `[[`, "S_value")
  best_ix <- which.min(S_vec)
  
  out                <- results[[best_ix]]
  out$optimal_k      <- param_grid$k[best_ix]
  out$optimal_m      <- param_grid$m[best_ix]
  out$all_results    <- results
  out
}
