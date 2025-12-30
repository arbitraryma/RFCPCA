## =========================================================
## End-to-end feature pipeline (WWW / MODWT only)
## Parallel + checkpoint + resume; E / N / T variants
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(future)
  library(future.apply)
  library(progressr)
})

## ---------- Thread hygiene (1 BLAS thread per worker) ----------
Sys.setenv(
  OMP_NUM_THREADS        = 1,
  MKL_NUM_THREADS        = 1,
  OPENBLAS_NUM_THREADS   = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS    = 1
)

## ---------- Parallel plan ----------
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(workers) || workers < 1) workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
handlers(global = TRUE)

## =========================================================
## 0) Minimal simulator (EEG-like, with bursts)
## =========================================================
simulate_eeg_mts <- function(n, p, fs = 100,
                             len_sec = 5,              # controls T
                             contam_prop = 0.20,       # fraction with bursts
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tlen <- as.integer(len_sec * fs)
  n1 <- floor(n/2); n2 <- n - n1
  labels <- c(rep(1L, n1), rep(2L, n2))
  
  phi_by_group <- c(`1` = 0.70, `2` = 0.40)
  sigma_eps    <- 1.0
  
  gen_one <- function(phi) {
    X <- matrix(0, nrow = Tlen, ncol = p)
    phi_j <- rnorm(p, mean = phi, sd = 0.05)
    for (j in seq_len(p)) {
      eps <- rnorm(Tlen, 0, sigma_eps)
      for (t in 2:Tlen) X[t, j] <- phi_j[j] * X[t-1, j] + eps[t]
    }
    X
  }
  
  ts_list <- vector("list", n)
  for (i in seq_len(n)) {
    g <- labels[i]
    ts_list[[i]] <- gen_one(phi_by_group[as.character(g)])
  }
  
  n_bad <- ceiling(contam_prop * n)
  if (n_bad > 0) {
    bad_idx <- sample.int(n, n_bad)
    for (i in bad_idx) {
      k_bursts <- sample(1:3, 1)
      for (b in seq_len(k_bursts)) {
        center <- sample(30:(Tlen-30), 1)
        width  <- sample(10:25, 1)
        amp    <- runif(1, 6, 10)
        win    <- max(1, center - width):min(Tlen, center + width)
        pulse  <- amp * dnorm(seq_along(win), mean = length(win)/2, sd = length(win)/8)
        ch <- sample.int(p, max(1, floor(p * runif(1, 0.1, 0.3))))
        ts_list[[i]][win, ch] <- ts_list[[i]][win, ch] +
          matrix(pulse, nrow = length(win), ncol = length(ch))
      }
      if (runif(1) < 0.3) {
        ch <- sample.int(p, floor(p * 0.2))
        ts_list[[i]][, ch] <- ts_list[[i]][, ch] + runif(1, -5, 5)
      }
    }
  }
  
  list(ts = ts_list, labels = labels)
}

## =========================================================
## 1) Lightweight feature extractors (fallbacks)
##    Replace these with your real WWW/MODWT functions when ready
## =========================================================
dis_www <- function(ts_list, features = TRUE) {
  # per-channel mean & sd concatenated
  feat_one <- function(X) c(colMeans(X), apply(X, 2, stats::sd))
  out <- do.call(rbind, lapply(ts_list, feat_one))
  scale(out)  # standardize for clustering
}
dis_modwt <- function(ts_list, features = TRUE) {
  # per-channel median & MAD concatenated
  feat_one <- function(X) c(apply(X, 2, stats::median), apply(X, 2, stats::mad))
  out <- do.call(rbind, lapply(ts_list, feat_one))
  scale(out)
}

## =========================================================
## 2) Fuzzy algorithms (dependency-free, fast)
## =========================================================

## -- Fuzzy C-means (C clusters), Euclidean --
fuzzy_c_means <- function(X, C = 2, m = 1.5, tol = 1e-4, maxit = 200, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X)
  N <- nrow(X)
  
  # kmeans++-like init: pick C centers
  centers <- X[sample.int(N, 1), , drop = FALSE]
  for (i in 2:C) {
    d2 <- apply(X, 1, function(x) min(colSums((t(centers) - x)^2)))
    prob <- d2 / sum(d2)
    centers <- rbind(centers, X[sample.int(N, 1, prob = prob), , drop = FALSE])
  }
  
  U <- matrix(0, N, C)
  dist2 <- function(M, ctrs) {
    # returns N x C squared distances
    sapply(seq_len(nrow(ctrs)), function(k) rowSums((M - ctrs[rep(k, N), ])^2))
  }
  
  for (it in 1:maxit) {
    D2 <- pmax(dist2(X, centers), .Machine$double.eps)
    pow <- 1 / (m - 1)
    inv <- (1 / D2)^pow
    U_new <- inv / rowSums(inv)
    
    # update centers
    Um <- U_new^m
    centers_new <- t(Um) %*% X / colSums(Um)
    
    # convergence checks
    du <- max(abs(U_new - U))
    dc <- max(abs(centers_new - centers))
    U <- U_new
    centers <- centers_new
    if (max(du, dc) < tol) break
  }
  list(U = U, centers = centers, iter = it, converged = (it < maxit))
}

## -- Trimmed fuzzy C-means (trim α worst by distance to nearest center) --
trimmed_fuzzy_c_means <- function(X, C = 2, m = 1.5, alpha = 0.2, tol = 1e-4, maxit = 200, seed = 1) {
  # clamp alpha to [0, 0.9] to avoid trimming all rows
  alpha <- max(0, min(alpha, 0.9))
  base <- fuzzy_c_means(X, C = C, m = m, tol = tol, maxit = maxit, seed = seed)
  centers <- base$centers
  # distance to nearest center
  D2 <- sapply(seq_len(nrow(centers)), function(k) rowSums((X - centers[rep(k, nrow(X)), ])^2))
  dmin <- apply(D2, 1, min)
  n_trim <- min(floor(alpha * nrow(X)), nrow(X) - 1L)  # ensure at least 1 kept
  trim_idx <- if (n_trim > 0) order(dmin, decreasing = TRUE)[seq_len(n_trim)] else integer(0)
  trimmed <- logical(nrow(X)); if (n_trim > 0) trimmed[trim_idx] <- TRUE
  list(U = base$U, centers = centers, trim_set = trimmed)
}

## -- FKM with noise (simple: mark points with max membership < τ as noise) --
FKM.noise <- function(X, k = 3, m = 1.5, tol = 1e-4, maxit = 200, seed = 1) {
  fc <- fuzzy_c_means(X, C = k, m = m, tol = tol, maxit = maxit, seed = seed)
  list(U = fc$U, centers = fc$centers)
}
predict_fkm_noise <- function(model, tau = 0.5) {
  U <- model$U
  maxu <- apply(U, 1, max)
  labels3 <- max.col(U, ties.method = "first")
  list(labels3 = labels3, noise_idx = which(maxu < tau))
}

## =========================================================
## 3) Pipeline helpers
## =========================================================
acc_best2 <- function(y_true, pred) {
  ok <- !is.na(pred); if (!any(ok)) return(NA_real_)
  y <- as.integer(y_true[ok]); p <- as.integer(pred[ok])
  Cmat <- table(factor(p, levels = sort(unique(p))),
                factor(y, levels = sort(unique(y))))
  if (nrow(Cmat) < 2 || ncol(Cmat) < 2) return(sum(y == p) / length(y))
  max((Cmat[1,1] + Cmat[2,2]), (Cmat[1,2] + Cmat[2,1])) / length(y)
}

eval_variants_on_features <- function(feat_mat, y, m_grid, alpha_trim = 0.2) {
  N <- length(y)
  out_rows <- vector("list", length(m_grid) * 3L); j <- 0L
  for (m in m_grid) {
    ## E: fuzzy c-means (2)
    resE <- try(fuzzy_c_means(feat_mat, C = 2, m = m), silent = TRUE)
    if (!inherits(resE, "try-error")) {
      labE <- max.col(resE$U, ties.method = "first")
      accE <- acc_best2(y, labE)
      outE <- 0
      j <- j + 1L; out_rows[[j]] <- tibble(variant="E", m=m, accuracy=accE, outlier_pct=outE, alpha=NA_real_)
    } else { j <- j + 1L; out_rows[[j]] <- tibble(variant="E", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=NA_real_) }
    
    ## N: noise model (k=3)
    resN <- try(FKM.noise(feat_mat, k = 3, m = m), silent = TRUE)
    if (!inherits(resN, "try-error")) {
      pr  <- predict_fkm_noise(resN, tau = 0.5)
      is_noise <- rep(FALSE, N); if (!is.null(pr$noise_idx)) is_noise[pr$noise_idx] <- TRUE
      keep  <- which(!is_noise)
      accN  <- if (length(keep)) acc_best2(y[keep], pr$labels3[keep]) else NA_real_
      outN  <- mean(is_noise)
      j <- j + 1L; out_rows[[j]] <- tibble(variant="N", m=m, accuracy=accN, outlier_pct=outN, alpha=NA_real_)
    } else { j <- j + 1L; out_rows[[j]] <- tibble(variant="N", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=NA_real_) }
    
    ## T: trimmed fuzzy c-means (2) — uses alpha_trim (now configurable)
    resT <- try(trimmed_fuzzy_c_means(feat_mat, C = 2, m = m, alpha = alpha_trim), silent = TRUE)
    if (!inherits(resT, "try-error")) {
      U      <- resT$U
      trim_v <- resT$trim_set
      keep   <- which(!trim_v)
      labT   <- max.col(U, ties.method = "first")
      accT   <- if (length(keep)) acc_best2(y[keep], labT[keep]) else NA_real_
      outT   <- mean(trim_v)
      j <- j + 1L; out_rows[[j]] <- tibble(variant="T", m=m, accuracy=accT, outlier_pct=outT, alpha=alpha_trim)
    } else { j <- j + 1L; out_rows[[j]] <- tibble(variant="T", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=alpha_trim) }
  }
  bind_rows(out_rows)
}

get_features <- function(family, ts_list) {
  if (family == "WWW")   return(dis_www(ts_list,  features = TRUE))
  if (family == "MODWT") return(dis_modwt(ts_list, features = TRUE))
  stop("Unknown family: ", family)
}

## =========================================================
## 4) Checkpoint paths & resume helpers
## =========================================================
ROOT        <- Sys.getenv("SCRATCH", ".")
CHKPT_ROOT  <- file.path(ROOT, "feat_chkpts")
RESULTS_DIR <- file.path(ROOT, "feat_results")
dir.create(CHKPT_ROOT,  showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

sc_dir       <- function(family, n, p) file.path(CHKPT_ROOT, family, sprintf("n%04d_p%03d", n, p))
rep_feats    <- function(dir, r) file.path(dir, sprintf("rep_%03d_feats.rds",  r))
rep_results  <- function(dir, r) file.path(dir, sprintf("rep_%03d_done.rds",   r))

scenario_status <- function(family, n, p, reps) {
  d <- sc_dir(family, n, p); dir.create(d, FALSE, TRUE)
  done_files <- list.files(d, "^rep_\\d+_done\\.rds$", full.names = FALSE)
  done_idx   <- sort(as.integer(sub("^rep_(\\d+)_done\\.rds$", "\\1", done_files)))
  list(dir = d, done = done_idx, missing = setdiff(seq_len(reps), done_idx))
}

reset_scenario <- function(family, n, p) {
  d <- sc_dir(family, n, p)
  if (dir.exists(d)) unlink(d, recursive = TRUE, force = TRUE)
  dir.create(d, FALSE, TRUE)
}
reset_all <- function(families, n_vals, p_vals) {
  grid <- expand_grid(family = families, n = n_vals, p = p_vals)
  for (i in seq_len(nrow(grid))) reset_scenario(grid$family[i], grid$n[i], grid$p[i])
  if (dir.exists(RESULTS_DIR)) unlink(RESULTS_DIR, recursive = TRUE, force = TRUE)
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
}

## =========================================================
## 5) Run one (family, n, p) with resume & feature cache
## =========================================================
run_family_np_resume <- function(family, n, p, m_grid, reps,
                                 base_seed = 2025, alpha_trim = 0.2,
                                 reset = FALSE, quiet = TRUE) {
  if (reset) reset_scenario(family, n, p)
  
  st <- scenario_status(family, n, p, reps)
  d  <- st$dir
  
  if (length(st$missing) > 0) {
    run_one <- function(r) {
      feats_file <- rep_feats(d, r)
      dat <- NULL
      if (file.exists(feats_file)) {
        dat <- readRDS(feats_file)
      } else {
        set.seed(base_seed + r)
        sim <- simulate_eeg_mts(
          n = n, p = p, fs = 100, len_sec = 5,
          contam_prop = 0.20, seed = base_seed + r
        )
        feats <- get_features(family, sim$ts)
        dat <- list(feats = as.matrix(feats), labels = as.integer(sim$labels))
        saveRDS(dat, feats_file)
      }
      
      res <- try(eval_variants_on_features(dat$feats, dat$labels, m_grid, alpha_trim), silent = TRUE)
      saveRDS(
        cbind(n = n, p = p, family = family,
              if (!inherits(res, "try-error")) res else
                tibble(variant=c("E","N","T"), m=NA_real_, accuracy=NA_real_,
                       outlier_pct=NA_real_, alpha=c(NA_real_, NA_real_, alpha_trim))),
        rep_results(d, r)
      )
      TRUE
    }
    
    with_progress({
      pbar <- progressor(steps = length(st$missing))
      future_lapply(
        st$missing,
        function(r){ on.exit(pbar(sprintf("rep %d", r))); run_one(r) },
        future.seed = TRUE
      )
    })
  }
  
  files <- file.path(d, sprintf("rep_%03d_done.rds", 1:reps))
  files <- files[file.exists(files)]
  out_list <- lapply(files, readRDS)
  out_list <- Filter(Negate(is.null), out_list)
  if (!length(out_list)) return(NULL)
  bind_rows(out_list)
}

## =========================================================
## 6) Sweep all scenarios
## =========================================================
run_all_features_resume <- function(families, n_vals, p_vals, m_grid, reps,
                                    alpha_trim = 0.2,
                                    reset_everything = FALSE, quiet = TRUE) {
  if (reset_everything) reset_all(families, n_vals, p_vals)
  
  grid <- tidyr::expand_grid(family = families, n = n_vals, p = p_vals)
  parts <- purrr::pmap(
    grid,
    ~ run_family_np_resume(..1, ..2, ..3, m_grid, reps,
                           alpha_trim = alpha_trim,
                           reset = FALSE, quiet = quiet)
  )
  raw  <- bind_rows(Filter(Negate(is.null), parts))
  
  raw_csv <- file.path(RESULTS_DIR, "feature_replicates.csv")
  write.csv(raw, raw_csv, row.names = FALSE)
  
  means <- raw %>%
    group_by(family, n, p, variant, m) %>%
    summarise(
      accuracy    = mean(accuracy, na.rm = TRUE),
      outlier_pct = mean(outlier_pct, na.rm = TRUE),
      alpha       = mean(alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(family, n, p, variant, m)
  
  means_csv <- file.path(RESULTS_DIR, "feature_means.csv")
  write.csv(means, means_csv, row.names = FALSE)
  means
}

## =======================
## 7) Config + Run
## =======================
families      <- c("WWW", "MODWT")     # QCD intentionally excluded
n_vals        <- c(400, 1000)
p_vals        <- c(32, 64, 128)        # adjust to 32, 64 for speed if needed
m_grid        <- c(1.5, 1.8, 2.0, 2.2)
reps          <- 2
TRIM_FRACTION <- 0.2                  # <<<<<< trim 40% in the T variant

MEANS <- run_all_features_resume(
  families, n_vals, p_vals, m_grid, reps,
  alpha_trim = TRIM_FRACTION,
  reset_everything = FALSE, quiet = TRUE
)
print(MEANS, n = Inf)
