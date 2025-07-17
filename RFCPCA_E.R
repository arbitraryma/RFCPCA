###############################################################################
## 1) Required Libraries
###############################################################################
library(base)
library(fungible)  

###############################################################################
## 2) Helper Functions  
###############################################################################

# 2.1) Compute the weighted cross-covariance matrices for each cluster
compute_weighted_cross_cov <- function(matrixU, m, cross_cov_lag0_lag1, cross_cov_lag0_lag2) {
  n <- nrow(matrixU)  # Number of MTS objects
  k <- ncol(matrixU)  # Number of clusters
  
  weighted_cross_cov_lag0_lag1 <- vector("list", k)
  weighted_cross_cov_lag0_lag2 <- vector("list", k)
  
  for (cluster in 1:k) {
    weights <- matrixU[, cluster]^m
    total_weight <- sum(weights)
    if (total_weight == 0) {
      stop("Total weight for cluster ", cluster, " is zero. Check your input matrixU or fuzziness parameter.")
    }
    
    weighted_cross_cov_lag0_lag1[[cluster]] <- Reduce(`+`, 
                                                      Map(function(cov_matrix, weight) weight * cov_matrix, cross_cov_lag0_lag1, weights)) / total_weight
    
    weighted_cross_cov_lag0_lag2[[cluster]] <- Reduce(`+`, 
                                                      Map(function(cov_matrix, weight) weight * cov_matrix, cross_cov_lag0_lag2, weights)) / total_weight
  }
  
  return(list(
    weighted_cross_cov_lag0_lag1 = weighted_cross_cov_lag0_lag1, 
    weighted_cross_cov_lag0_lag2 = weighted_cross_cov_lag0_lag2
  ))
}

# 2.2) Standardize each MTS object (center & scale columns)
standardize_mts <- function(mts_data) {
  lapply(mts_data, function(ts) {
    scale(ts, center = TRUE, scale = TRUE)
  })
}

# 2.3) Cross-covariance at lag 1 (returns block matrix + X_t + X_t+1)
cross_covariance_lag1 <- function(mts_data) {
  n <- nrow(mts_data)
  p <- ncol(mts_data)
  
  if (n < 2) {
    stop("The multivariate time series should have at least 2 time points.")
  }
  
  Gamma_0 <- cov(mts_data)
  X_t <- mts_data[1:(n-1), ]
  X_t_plus_1 <- mts_data[2:n, ]
  Gamma_1 <- cov(X_t, X_t_plus_1)
  
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_1),
    cbind(t(Gamma_1), Gamma_0)
  )
  output <- list(block_matrix, X_t, X_t_plus_1)
  return(output)
}

# 2.4) Cross-covariance at lag 2 (returns block matrix + X_t + X_t+2)
cross_covariance_lag2 <- function(mts_data) {
  n <- nrow(mts_data)
  p <- ncol(mts_data)
  
  if (n < 3) {
    stop("The multivariate time series should have at least 3 time points for lag=2.")
  }
  
  Gamma_0 <- cov(mts_data)
  X_t <- mts_data[1:(n-2), ]
  X_t_plus_2 <- mts_data[3:n, ]
  Gamma_2 <- cov(X_t, X_t_plus_2)
  
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_2),
    cbind(t(Gamma_2), Gamma_0)
  )
  output <- list(block_matrix, X_t, X_t_plus_2)
  return(output)
}

# 2.5) Combine X_t and X_t_plus_1 into one matrix (dimension ~ (n-1) x 2p)
combine_xt_xt1 <- function(X_t, X_t_plus_1) {
  combined_matrix <- cbind(X_t, X_t_plus_1)
  return(combined_matrix)
}

# 2.6) S-value computations (unused or used for cluster validity checks)
compute_S <- function(reconstruction_error, projection_axes, n, k) {
  numerator <- sum(reconstruction_error)
  
  P_list <- vector("list", k)
  for (c in seq_len(k)) {
    Clag1 <- projection_axes[[c]]$lag01
    Clag2 <- projection_axes[[c]]$lag02
    P_lag1 <- Clag1 %*% t(Clag1)
    P_lag2 <- Clag2 %*% t(Clag2)
    P_list[[c]] <- list(P_lag1 = P_lag1, P_lag2 = P_lag2)
  }
  
  subspace_dist <- function(Pc, Pm) {
    d_lag1 <- norm(Pc$P_lag1 - Pm$P_lag1, type = "F")
    d_lag2 <- norm(Pc$P_lag2 - Pm$P_lag2, type = "F")
    return(d_lag1 + d_lag2)
  }
  
  min_dist <- Inf
  for (c in seq_len(k)) {
    for (m in seq_len(k)) {
      if (m > c) {
        dist_val <- subspace_dist(P_list[[c]], P_list[[m]])
        if (dist_val < min_dist) {
          min_dist <- dist_val
        }
      }
    }
  }
  
  denominator <- min_dist
  S_value <- numerator / (n * denominator)
  return(S_value)
}

compute_S_cov <- function(reconstruction_error, weighted_covariances, n, k) {
  numerator <- reconstruction_error
  cov_dist <- function(cluster_c, cluster_m) {
    cov_c_lag1 <- weighted_covariances$weighted_cross_cov_lag0_lag1[[cluster_c]]
    cov_c_lag2 <- weighted_covariances$weighted_cross_cov_lag0_lag2[[cluster_c]]
    cov_m_lag1 <- weighted_covariances$weighted_cross_cov_lag0_lag1[[cluster_m]]
    cov_m_lag2 <- weighted_covariances$weighted_cross_cov_lag0_lag2[[cluster_m]]
    
    d_lag1 <- norm(cov_c_lag1 - cov_m_lag1, type = "F")
    d_lag2 <- norm(cov_c_lag2 - cov_m_lag2, type = "F")
    return(d_lag1 + d_lag2)
  }
  
  min_dist <- Inf
  for (c_idx in seq_len(k)) {
    for (m_idx in seq_len(k)) {
      if (m_idx > c_idx) {
        dist_val <- cov_dist(c_idx, m_idx)
        if (dist_val < min_dist) {
          min_dist <- dist_val
        }
      }
    }
  }
  
  denominator <- min_dist
  S_value <- numerator / (n * denominator)
  return(S_value)
}

###############################################################################
## 3) Revised FCPCA function with Additional Returned Outcomes
###############################################################################

fcpca <- function(ts, 
                  k, 
                  m = 1.5, 
                  startU = NULL, 
                  conver = 1e-3, 
                  maxit = 1000, 
                  verbose = TRUE, 
                  replicates = 1) {
  # ts: List of N MTS objects
  # k:  Number of clusters
  # m:  Fuzziness parameter
  # startU: Initial membership matrix
  # conver: Convergence criterion
  # maxit: Maximum number of iterations
  # verbose: Print iteration details
  # replicates: # of times to repeat with random initial U (choose best)
  
  # 3.1) Standardize MTS
  ts <- standardize_mts(ts)
  
  n <- length(ts)   # Number of MTS objects
  p <- ncol(ts[[1]])  # # of variables in each MTS
  
  # 3.2) For each MTS, compute cross-covariance at lag1 and lag2
  results  <- lapply(ts, cross_covariance_lag1)
  results2 <- lapply(ts, cross_covariance_lag2)
  
  # 3.3) Extract block matrices, X_t, X_t+1, etc.
  sigma  <- lapply(results,  function(x) x[[1]]) # block matrix (lag1)
  sigma2 <- lapply(results2, function(x) x[[1]]) # block matrix (lag2)
  
  tsxt  <- lapply(results,  function(x) x[[2]]) 
  tsxt1 <- lapply(results,  function(x) x[[3]])
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE)
  
  tsxt2  <- lapply(results2, function(x) x[[2]])
  tsxt12 <- lapply(results2, function(x) x[[3]])
  combined_list2 <- mapply(combine_xt_xt1, tsxt2, tsxt12, SIMPLIFY = FALSE)
  
  # 3.4) Possibly random initial membership matrix
  if (is.null(startU)) {
    U <- matrix(runif(n * k), nrow = n)
    U <- U / rowSums(U) # normalize each row
  } else {
    U <- startU
  }
  
  # 3.5) We'll keep track of the best solution over multiple replicates
  best_result <- NULL
  best_error  <- Inf
  
  for (replicate_i in seq_len(replicates)) {
    if (verbose) cat(sprintf("\nReplicate %d:\n", replicate_i))
    
    # Re-init membership if no startU provided
    if (is.null(startU)) {
      U <- matrix(runif(n * k), nrow = n)
      U <- U / rowSums(U)
    }
    
    iteration   <- 0
    diff        <- Inf
    prev_error  <- Inf
    current_err <- Inf
    
    # 3.6) Iteration
    while (iteration < maxit && (abs(prev_error - current_err) > conver || diff > conver)) {
      iteration <- iteration + 1
      
      # Weighted cross-cov
      weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
      
      # For cluster 1..k, do SVD on the weighted cross-cov to get principal axes
      projection_axes <- vector("list", k)
      
      for (c_idx in 1:k) {
        # SVD lag1
        svd_lag01 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[c_idx]])
        # pick #PCs that capture >=95% variance
        cumvar1 <- cumsum(svd_lag01$d) / sum(svd_lag01$d)
        n_components_lag01 <- which(cumvar1 >= 0.95)[1]
        if (is.na(n_components_lag01)) n_components_lag01 <- length(svd_lag01$d)
        
        proj_axes_lag01 <- svd_lag01$u[, 1:n_components_lag01, drop=FALSE]
        
        # SVD lag2
        svd_lag02 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[c_idx]])
        cumvar2 <- cumsum(svd_lag02$d) / sum(svd_lag02$d)
        n_components_lag02 <- which(cumvar2 >= 0.95)[1]
        if (is.na(n_components_lag02)) n_components_lag02 <- length(svd_lag02$d)
        
        proj_axes_lag02 <- svd_lag02$u[, 1:n_components_lag02, drop=FALSE]
        
        projection_axes[[c_idx]] <- list(
          lag01 = proj_axes_lag01,
          lag02 = proj_axes_lag02,
          # store #PCs used for each cluster in case we want them
          k1 = n_components_lag01,
          k2 = n_components_lag02
        )
      }
      
      # Compute reconstruction error for each series i, cluster c
      total_reconstruction_error <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c_idx in 1:k) {
          rec_tt1 <- combined_list[[i]]  %*% projection_axes[[c_idx]]$lag01 %*% 
            t(projection_axes[[c_idx]]$lag01)
          err_tt1 <- norm(combined_list[[i]] - rec_tt1, "F")^2
          
          rec_tt2 <- combined_list2[[i]] %*% projection_axes[[c_idx]]$lag02 %*% 
            t(projection_axes[[c_idx]]$lag02)
          err_tt2 <- norm(combined_list2[[i]] - rec_tt2, "F")^2
          
          total_reconstruction_error[i, c_idx] <- err_tt1 + err_tt2
        }
      }
      
      # Update membership matrix
      U_new <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c_idx in 1:k) {
          denom_term <- sum((total_reconstruction_error[i, c_idx] / (total_reconstruction_error[i, ]+1e-15))^(1/(m-1)))
          U_new[i, c_idx] <- 1 / denom_term
        }
      }
      
      diff       <- max(abs(U - U_new))
      prev_error <- current_err
      U          <- U_new
      current_err <- sum(U^m * total_reconstruction_error)
      
      if (verbose) {
        cat(sprintf("Iteration %d: Reconstruction Error = %.6f, Max Change in U = %.6f\n",
                    iteration, current_err, diff))
      }
    } # end while
    
    # Compare with best replicate
    if (current_err < best_error) {
      best_error <- current_err
      best_result <- list(
        membership_matrix = U,
        projection_axes   = projection_axes,
        iterations        = iteration,
        converged         = (iteration < maxit),
        reconstruction_error       = current_err,
        reconstruction_error_matrix= total_reconstruction_error,
        hard_cluster      = max.col(U),
        weighted_covariances = compute_weighted_cross_cov(U, m, sigma, sigma2)  # final
      )
      # We'll also store the final #PCs used in the best solution for lag1, lag2
      # For simplicity, we use cluster 1's k1, k2 as representative 
      best_result$k1 <- projection_axes[[1]]$k1
      best_result$k2 <- projection_axes[[1]]$k2
    }
  } # end replicates
  
  # 3.7) Add the sigma, sigma2, combined_list, combined_list2 to the returned best_result
  best_result$sigma          = sigma
  best_result$sigma2         = sigma2
  best_result$combined_list  = combined_list
  best_result$combined_list2 = combined_list2
  
  return(best_result)
}





###############################################################################
## Revised FCPCA Function with Additional Outputs
###############################################################################

fcpca <- function(ts, 
                  k, 
                  m = 1.5, 
                  startU = NULL, 
                  conver = 1e-3, 
                  maxit = 1000, 
                  verbose = TRUE, 
                  replicates = 1) {
  # ts: List of N MTS objects
  # k: Number of clusters
  # m: Fuzziness parameter
  # startU: Initial membership matrix
  # conver: Convergence criterion
  # maxit: Maximum number of iterations
  # verbose: Whether to display iteration details
  # replicates: Number of times to repeat clustering
  
  # 1) Standardize MTS objects
  ts <- standardize_mts(ts)
  
  n <- length(ts)     # Number of MTS objects
  p <- ncol(ts[[1]])  # Dimensionality (variables) of each MTS
  
  # 2) Apply the cross_covariance function to each MTS object
  #    to get lag-1 and lag-2 block matrices, plus X_t, X_{t+1}, etc.
  results  <- lapply(ts, cross_covariance_lag1) 
  results2 <- lapply(ts, cross_covariance_lag2)
  
  # 3) Extract the block matrices (sigma, sigma2) and combined lists
  sigma   <- lapply(results,  function(x) x[[1]])  # dimension 2p x 2p for lag-1
  sigma2  <- lapply(results2, function(x) x[[1]])  # dimension 2p x 2p for lag-2
  
  tsxt    <- lapply(results,  function(x) x[[2]])  # X_t
  tsxt1   <- lapply(results,  function(x) x[[3]])  # X_{t+1}
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE)
  
  tsxt2   <- lapply(results2, function(x) x[[2]])  # X_t
  tsxt12  <- lapply(results2, function(x) x[[3]])  # X_{t+2}
  combined_list2 <- mapply(combine_xt_xt1, tsxt2, tsxt12, SIMPLIFY = FALSE)
  
  # 4) Possibly random initialization of membership matrix
  if (is.null(startU)) {
    U <- matrix(runif(n * k), nrow = n)
    U <- U / rowSums(U)  # normalize each row
  } else {
    U <- startU
  }
  
  # 5) Compute weighted covariance for the initial membership to determine #PCs
  weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
  
  # 5a) For lag-1, pick #PCs from cluster 1's weighted cross-cov as reference
  svd_result_lag01 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[1]])
  n_components_lag01 <- which(cumsum(svd_result_lag01$d) / sum(svd_result_lag01$d) >= 0.95)[1]
  if (is.na(n_components_lag01)) {
    # fallback if the cumsum never goes above 0.95
    n_components_lag01 <- length(svd_result_lag01$d)
  }
  
  # 5b) For lag-2, pick #PCs similarly
  svd_result_lag02 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[1]])
  n_components_lag02 <- which(cumsum(svd_result_lag02$d) / sum(svd_result_lag02$d) >= 0.95)[1]
  if (is.na(n_components_lag02)) {
    n_components_lag02 <- length(svd_result_lag02$d)
  }
  
  # 6) Multiple replicates to mitigate poor initialization
  best_result <- NULL
  best_error  <- Inf
  
  for (replicate_i in 1:replicates) {
    if (verbose) cat(sprintf("\nReplicate %d:\n", replicate_i))
    
    # (Re-)Initialize membership if startU not given
    if (is.null(startU)) {
      U <- matrix(runif(n * k), nrow = n)
      U <- U / rowSums(U)
    } else {
      U <- startU
    }
    
    iteration     <- 0
    diff          <- Inf
    prev_error    <- Inf
    current_error <- Inf
    
    # 7) Main iteration loop
    while (iteration < maxit && (abs(prev_error - current_error) > conver || diff > conver)) {
      iteration <- iteration + 1
      
      # 7a) Compute weighted covariances
      weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
      
      # 7b) Update projection axes using the fixed #PCs (n_components_lag01, n_components_lag02)
      projection_axes <- vector("list", k)
      for (c in 1:k) {
        # SVD for lag-1
        svd_result_lag01 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[c]])
        proj_axes_lag01  <- svd_result_lag01$u[, 1:n_components_lag01, drop = FALSE]
        
        # SVD for lag-2
        svd_result_lag02 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[c]])
        proj_axes_lag02  <- svd_result_lag02$u[, 1:n_components_lag02, drop = FALSE]
        
        projection_axes[[c]] <- list(lag01 = proj_axes_lag01, lag02 = proj_axes_lag02)
      }
      
      # 7c) Compute reconstruction errors for each series i, cluster c
      total_reconstruction_error <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c in 1:k) {
          reconstructed_tt1 <- combined_list[[i]]  %*% projection_axes[[c]]$lag01 %*% 
            t(projection_axes[[c]]$lag01)
          err_tt1 <- norm(combined_list[[i]] - reconstructed_tt1, "F")^2
          
          reconstructed_tt2 <- combined_list2[[i]] %*% projection_axes[[c]]$lag02 %*% 
            t(projection_axes[[c]]$lag02)
          err_tt2 <- norm(combined_list2[[i]] - reconstructed_tt2, "F")^2
          
          total_reconstruction_error[i, c] <- err_tt1 + err_tt2
        }
      }
      
      # 7d) Update membership matrix
      U_new <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c in 1:k) {
          # fuzzy membership denominator
          denom <- sum((total_reconstruction_error[i, c] / (total_reconstruction_error[i, ] + 1e-10))^(1 / (m - 1)))
          U_new[i, c] <- 1 / denom
        }
      }
      
      # 7e) Convergence checks
      diff       <- max(abs(U - U_new))
      prev_error <- current_error
      U          <- U_new
      current_error <- sum(U^m * total_reconstruction_error)
      
      if (verbose) {
        cat(sprintf("Iteration %d: Reconstruction Error = %.6f, Max Change in U = %.6f\n",
                    iteration, current_error, diff))
      }
    } # end while
    
    # 7f) Check if this replicate is the best
    if (current_error < best_error) {
      best_error <- current_error
      best_result <- list(
        membership_matrix           = U,
        projection_axes             = projection_axes,
        iterations                  = iteration,
        converged                   = (iteration < maxit),
        reconstruction_error        = current_error,
        reconstruction_error_matrix = total_reconstruction_error,
        hard_cluster                = max.col(U),
        weighted_covariances        = weighted_covariances,
        S                           = compute_S(current_error, projection_axes, n, k)
      )
    }
  } # end for replicates
  
  # 8) Add the extra outputs (sigma, sigma2, combined_list, combined_list2, and #PCs) to best_result
  best_result$sigma          <- sigma
  best_result$sigma2         <- sigma2
  best_result$combined_list  <- combined_list
  best_result$combined_list2 <- combined_list2
  best_result$n_components_lag01 <- n_components_lag01
  best_result$n_components_lag02 <- n_components_lag02
  
  # 9) Return final best result
  return(best_result)
}







###############################################################################
## 4) Exponential (Robust) FCPCA Function: RFCPCA_E
###############################################################################

RFCPCA_E <- function(ts, k, m = 1.5, startU = NULL, conver = 1e-3, maxit = 1000, verbose = TRUE) {
  # 4.1) First run the standard FCPCA to get an initial solution
  initial_res <- fcpca(ts, k, m = m, startU = startU, conver = conver, maxit = maxit, 
                       verbose = verbose, replicates = 1)
  
  # Extract necessary items
  U      <- initial_res$membership_matrix
  projAx <- initial_res$projection_axes
  sigma  <- initial_res$sigma
  sigma2 <- initial_res$sigma2
  comb1  <- initial_res$combined_list
  comb2  <- initial_res$combined_list2
  k1     <- initial_res$n_components_lag01
  k2     <- initial_res$n_components_lag02
  
  n <- nrow(U)   # number of MTS
  # We'll also store final cluster count from fcpca
  k_clust <- ncol(U)
  
  # 4.2) Compute reconstruction errors for each series i, cluster s

  rec_error <- initial_res$reconstruction_error_matrix

  
  # 4.3) Compute beta as inverse of average minimal error
   min_err_per_series <- apply(rec_error, 1, min)
   beta <- 1 / mean(min_err_per_series)

  
  if (verbose) {
    cat("\nExponential FCPCA: Computed beta =", beta, "\n")
  }
  
  # 4.4) Now iterate using the exponential objective
  iteration <- 0
  diffU     <- Inf
  total_obj <- Inf
  prev_obj  <- Inf
  
  repeat {
    iteration <- iteration + 1
    
    # 4.4a) Update membership matrix using exponential distance
    # D_{is} = 1 - exp(-beta * rec_error[i,s])
    D_mat <- 1 - exp(-beta * rec_error)
    # Avoid exact zero
    D_mat[D_mat < .Machine$double.eps] <- .Machine$double.eps
    
    U_new <- matrix(0, nrow = n, ncol = k_clust)
    for (i in 1:n) {
      denom_i <- sum(D_mat[i, ]^(-1/(m - 1)))
      for (s in 1:k_clust) {
        U_new[i, s] <- D_mat[i, s]^(-1/(m - 1)) / denom_i
      }
    }
    
    diffU <- max(abs(U - U_new))
    U <- U_new
    
    # 4.4b) Recompute weighted cross-covariances with updated memberships
    weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
    
    # 4.4c) Update projection axes (SVD) with same k1, k2 as in initial fcpca
    projAx_new <- vector("list", k_clust)
    for (s in 1:k_clust) {
      W_lag1 <- weighted_covariances$weighted_cross_cov_lag0_lag1[[s]]
      W_lag2 <- weighted_covariances$weighted_cross_cov_lag0_lag2[[s]]
      
      svd1 <- svd(W_lag1)
      svd2 <- svd(W_lag2)
      # pick top k1, k2 from the *initially determined* dims
      U1 <- svd1$u[, 1:k1, drop=FALSE]
      U2 <- svd2$u[, 1:k2, drop=FALSE]
      
      projAx_new[[s]] <- list(lag01 = U1, lag02 = U2)
    }
    
    # 4.4d) Recompute reconstruction errors rec_error with updated axes
    rec_error_new <- matrix(0, nrow = n, ncol = k_clust)
    for (i in 1:n) {
      for (s in 1:k_clust) {
        rec_tt1 <- comb1[[i]] %*% projAx_new[[s]]$lag01 %*% t(projAx_new[[s]]$lag01)
        err_tt1 <- norm(comb1[[i]] - rec_tt1, "F")^2
        
        rec_tt2 <- comb2[[i]] %*% projAx_new[[s]]$lag02 %*% t(projAx_new[[s]]$lag02)
        err_tt2 <- norm(comb2[[i]] - rec_tt2, "F")^2
        
        rec_error_new[i, s] <- err_tt1 + err_tt2
      }
    }
    
    # 4.4e) Compute new exponential objective: sum_{i,s} u_{is}^m * [1 - exp(-beta * E_{is})]
    D_mat_new <- 1 - exp(-beta * rec_error_new)
    J_exp <- sum((U^m) * D_mat_new)
    
    obj_change <- abs(prev_obj - J_exp)
    prev_obj <- J_exp
    
    if (verbose) {
      cat(sprintf("Iteration %d (ExpFCPCA): J_exp = %.6f, MembershipChange = %.6f, ObjChange = %.6f\n",
                  iteration, J_exp, diffU, obj_change))
    }
    
    # Check for convergence
    if ((diffU < conver && obj_change < conver) || iteration >= maxit) {
      # finalize
      break
    }
    
    # Prepare for next iteration
    U_old <- U
    rec_error <- rec_error_new
    projAx <- projAx_new
  }
  
  # 4.5) Return final results of robust exponential FCPCA
  result <- list(
    beta                        = beta,
    membership_matrix           = U,
    reconstruction_error_matrix = rec_error,
    total_exponential_objective = prev_obj,
    projection_axes             = projAx,
    # Weighted cross-cov for the final iteration if needed
    weighted_covariances        = compute_weighted_cross_cov(U, m, sigma, sigma2),
    hard_cluster = max.col(U),
    # Possibly compute S from compute_S_cov or compute_S
    # e.g., we can compute S with total reconstruction error or with weighted cov
    # For example, let's store the sum of rec_error:
    S_value = compute_S_cov(sum(U^m * rec_error), 
                            #compute_weighted_cross_cov(U, m, sigma, sigma2), 
                            weighted_covariances,
                            n, 
                            k_clust)
  )
  
  return(result)
}




# the parallel computing, you can have fixed k and m, if at least one of them is not specified, 
# the RFCPCA_E will do the automatically selection 
RFCPCA_E <- function(ts, 
                     k = NULL, 
                     m = NULL, 
                     startU = NULL, 
                     conver = 1e-3, 
                     maxit = 1000, 
                     verbose = TRUE, 
                     parallel = TRUE) {
  
  # Load necessary packages
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' needed for parallel execution. Please install it.")
    }
    library(furrr)
    plan(multisession)  # Use multiple R sessions (good default)
  }
  
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1 || length(m) > 1
  
  # Set default search ranges if needed
  if (auto_mode) {
    k_range <- if (is.null(k)) 2:6 else k
    m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5) else m
  } else {
    k_range <- k
    m_range <- m
  }
  
  # Expand all (k, m) combinations
  param_grid <- expand.grid(k = k_range, m = m_range)
  
  # Helper function: run one RFCPCA_E pass
  run_RFCPCA_E_once <- function(ts, k, m, startU, conver, maxit, verbose) {
    initial_res <- fcpca(ts, k, m = m, startU = startU, conver = conver, maxit = maxit, verbose = FALSE, replicates = 1)
    
    U      <- initial_res$membership_matrix
    projAx <- initial_res$projection_axes
    sigma  <- initial_res$sigma
    sigma2 <- initial_res$sigma2
    comb1  <- initial_res$combined_list
    comb2  <- initial_res$combined_list2
    k1     <- initial_res$n_components_lag01
    k2     <- initial_res$n_components_lag02
    
    n <- nrow(U)
    k_clust <- ncol(U)
    
    rec_error <- initial_res$reconstruction_error_matrix
    
    min_err_per_series <- apply(rec_error, 1, min)
    beta <- 1 / (mean(min_err_per_series) + 1e-8)
    
    iteration <- 0
    diffU <- Inf
    prev_obj <- Inf
    total_obj <- Inf
    
    repeat {
      iteration <- iteration + 1
      
      D_mat <- 1 - exp(-beta * rec_error)
      D_mat[D_mat < .Machine$double.eps] <- .Machine$double.eps
      
      U_new <- matrix(0, nrow = n, ncol = k_clust)
      for (i in 1:n) {
        denom_i <- sum(D_mat[i, ]^(-1/(m - 1)))
        for (s in 1:k_clust) {
          U_new[i, s] <- D_mat[i, s]^(-1/(m - 1)) / denom_i
        }
      }
      
      diffU <- max(abs(U - U_new))
      U <- U_new
      
      weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
      
      projAx_new <- vector("list", k_clust)
      for (s in 1:k_clust) {
        svd1 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[s]])
        svd2 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[s]])
        U1 <- svd1$u[, 1:k1, drop = FALSE]
        U2 <- svd2$u[, 1:k2, drop = FALSE]
        projAx_new[[s]] <- list(lag01 = U1, lag02 = U2)
      }
      
      rec_error_new <- matrix(0, nrow = n, ncol = k_clust)
      for (i in 1:n) {
        for (s in 1:k_clust) {
          rec_tt1 <- comb1[[i]] %*% projAx_new[[s]]$lag01 %*% t(projAx_new[[s]]$lag01)
          err_tt1 <- norm(comb1[[i]] - rec_tt1, "F")^2
          
          rec_tt2 <- comb2[[i]] %*% projAx_new[[s]]$lag02 %*% t(projAx_new[[s]]$lag02)
          err_tt2 <- norm(comb2[[i]] - rec_tt2, "F")^2
          
          rec_error_new[i, s] <- err_tt1 + err_tt2
        }
      }
      
      D_mat_new <- 1 - exp(-beta * rec_error_new)
      total_obj <- sum((U^m) * D_mat_new)
      
      obj_change <- abs(prev_obj - total_obj)
      prev_obj <- total_obj
      
      if ((diffU < conver && obj_change < conver) || iteration >= maxit) {
        break
      }
      
      rec_error <- rec_error_new
      projAx <- projAx_new
    }
    
    final_weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
    S_value <- compute_S_cov(sum(U^m * rec_error), final_weighted_covariances, n, k_clust)
    
    return(list(
      membership_matrix = U,
      projection_axes = projAx,
      reconstruction_error_matrix = rec_error,
      weighted_covariances = final_weighted_covariances,
      beta = beta,
      S_value = S_value,
      k = k_clust,
      m = m
    ))
  }
  
  # --- Parallel vs Sequential Execution ---
  if (auto_mode) {
    if (parallel) {
      # --- Parallel search using future_map ---
      results_list <- furrr::future_map(1:nrow(param_grid), function(idx) {
        current_k <- param_grid$k[idx]
        current_m <- param_grid$m[idx]
        if (verbose) cat(sprintf("Running parallel: k = %d, m = %.2f\n", current_k, current_m))
        run_RFCPCA_E_once(ts, current_k, current_m, startU, conver, maxit, verbose = FALSE)
      }, .options = furrr_options(seed = TRUE))
    } else {
      # --- Sequential search with early stop ---
      results_list <- list()
      previous_S <- Inf
      for (idx in 1:nrow(param_grid)) {
        current_k <- param_grid$k[idx]
        current_m <- param_grid$m[idx]
        if (verbose) cat(sprintf("Running sequential: k = %d, m = %.2f\n", current_k, current_m))
        
        res <- run_RFCPCA_E_once(ts, current_k, current_m, startU, conver, maxit, verbose = FALSE)
        results_list[[idx]] <- res
        
        current_S <- res$S_value
        if (current_S > previous_S * 100) {
          if (verbose) cat("Early stopping triggered!\n")
          break
        }
        previous_S <- min(previous_S, current_S)
      }
    }
    
    # --- Find the best result ---
    S_vals <- sapply(results_list, function(x) x$S_value)
    best_idx <- which.min(S_vals)
    
    final_result <- results_list[[best_idx]]
    final_result$optimal_k <- param_grid$k[best_idx]
    final_result$optimal_m <- param_grid$m[best_idx]
    final_result$all_results <- results_list
  } else {
    # --- Just run one if k and m specified ---
    final_result <- run_RFCPCA_E_once(ts, k, m, startU, conver, maxit, verbose)
    final_result$optimal_k <- k
    final_result$optimal_m <- m
    final_result$all_results <- list(final_result)
  }
  
  return(final_result)
}
