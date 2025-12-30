## =========================================================
## WWW / MODWT Feature Pipeline with Fixed-α Trimming
## - Parallel + checkpoint + resume
## - Metrics for trimming and inlier-only clustering accuracy
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
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

## =======================
## 1) Config
## =======================
families      <- c("WWW", "MODWT")      # QCD intentionally excluded
n_vals        <- c(400, 1000)
p_vals        <- c(32, 64, 128)
m_grid        <- c(1.5, 1.8, 2.0, 2.2)
reps          <- 2

## evaluation thresholds
TAU_DOM     <- 0.7   # E: flag outlier if max(u_i•) < 0.7 (dominance)
GAMMA_NOISE <- 0.5   # N: flag outlier if noise-membership >= 0.5
NOISE_COL   <- 3L    # N: assume the 3rd fuzzy cluster is the noise cluster (k=3)

## trimming amount for T variant (fixed α)
TRIM_FRACTION <- 0.40

## checkpointing
ROOT <- Sys.getenv("SCRATCH", ".")
FRESH_RUN_ID <- TRUE  # TRUE -> new run folder (no resume); FALSE -> reuse tag
RUN_TAG <- if (FRESH_RUN_ID) format(Sys.time(), "%Y%m%d_%H%M%S") else "stable"

alpha_tag   <- gsub("\\.", "p", sprintf("alpha_%0.2f", TRIM_FRACTION))
BASE_TAG    <- paste0(alpha_tag, "__", RUN_TAG)
CHKPT_ROOT  <- file.path(ROOT, "feat_chkpts",  BASE_TAG)
RESULTS_DIR <- file.path(ROOT, "feat_results", BASE_TAG)
dir.create(CHKPT_ROOT,  showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## =======================
## 2) Simulator (EEG-like, with bursts) + ground truth outliers
## =======================
simulate_eeg_mts <- function(n, p, fs = 100, len_sec = 5, contam_prop = 0.20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Tlen <- as.integer(len_sec * fs)
  n1 <- floor(n/2); n2 <- n - n1
  labels <- c(rep(1L, n1), rep(2L, n2))
  
  phi_by_group <- c(`1` = 0.70, `2` = 0.40); sigma_eps <- 1.0
  
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
  
  is_outlier <- rep(FALSE, n)
  n_bad <- ceiling(contam_prop * n)
  if (n_bad > 0) {
    bad_idx <- sample.int(n, n_bad)
    is_outlier[bad_idx] <- TRUE
    for (i in bad_idx) {
      k_bursts <- sample(1:3, 1)
      for (b in seq_len(k_bursts)) {
        center <- sample(30:(Tlen-30), 1)
        width  <- sample(10:25, 1)
        amp    <- runif(1, 6, 10)
        win    <- max(1, center - width):min(Tlen, center + width)
        pulse  <- amp * dnorm(seq_along(win), mean = length(win)/2, sd = length(win)/8)
        ch     <- sample.int(p, max(1, floor(p * runif(0.1, 0.3))))
        ts_list[[i]][win, ch] <- ts_list[[i]][win, ch] +
          matrix(pulse, nrow = length(win), ncol = length(ch))
      }
      if (runif(1) < 0.3) { # DC shift on 20% channels
        ch <- sample.int(p, floor(p * 0.2))
        ts_list[[i]][, ch] <- ts_list[[i]][, ch] + runif(1, -5, 5)
      }
    }
  }
  
  list(ts = ts_list, labels = labels, is_outlier = is_outlier)
}

## =======================
## 3) Lightweight feature extractors (placeholders)
## =======================
dis_www <- function(ts_list, features = TRUE) {
  feat_one <- function(X) c(colMeans(X), apply(X, 2, stats::sd))
  out <- do.call(rbind, lapply(ts_list, feat_one))
  scale(out)
}
dis_modwt <- function(ts_list, features = TRUE) {
  feat_one <- function(X) c(apply(X, 2, stats::median), apply(X, 2, stats::mad))
  out <- do.call(rbind, lapply(ts_list, feat_one))
  scale(out)
}
get_features <- function(family, ts_list) {
  if (family == "WWW")   return(dis_www(ts_list,  TRUE))
  if (family == "MODWT") return(dis_modwt(ts_list, TRUE))
  stop("Unknown family: ", family)
}

## =======================
## 4) Fuzzy algorithms
## =======================
fuzzy_c_means <- function(X, C = 2, m = 1.5, tol = 1e-4, maxit = 200, seed = 1) {
  set.seed(seed); X <- as.matrix(X); N <- nrow(X)
  centers <- X[sample.int(N, 1), , drop = FALSE]
  for (i in 2:C) {
    d2 <- apply(X, 1, function(x) min(colSums((t(centers) - x)^2)))
    prob <- d2 / sum(d2)
    centers <- rbind(centers, X[sample.int(N, 1, prob = prob), , drop = FALSE])
  }
  U <- matrix(0, N, C)
  dist2 <- function(M, ctrs)
    sapply(seq_len(nrow(ctrs)), function(k) rowSums((M - ctrs[rep(k, N), ])^2))
  
  for (it in 1:maxit) {
    D2 <- pmax(dist2(X, centers), .Machine$double.eps)
    pow <- 1 / (m - 1)
    inv <- (1 / D2)^pow
    U_new <- inv / rowSums(inv)
    Um <- U_new^m
    centers_new <- t(Um) %*% X / colSums(Um)
    if (max(abs(U_new - U), abs(centers_new - centers)) < tol) {
      U <- U_new; centers <- centers_new; break
    }
    U <- U_new; centers <- centers_new
  }
  list(U = U, centers = centers)
}

## --- Fixed-fraction trimmed FCM (α = TRIM_FRACTION) ---
trimmed_fuzzy_c_means <- function(X, C = 2, m = 1.5, alpha = 0.4, tol = 1e-4, maxit = 200, seed = 1) {
  alpha <- max(0, min(alpha, 0.9))
  base <- fuzzy_c_means(X, C = C, m = m, tol = tol, maxit = maxit, seed = seed)
  centers <- base$centers
  D2 <- sapply(seq_len(nrow(centers)), function(k) rowSums((X - centers[rep(k, nrow(X)), ])^2))
  dmin <- apply(D2, 1, min)
  n_trim <- min(floor(alpha * nrow(X)), nrow(X) - 1L)
  trim_idx <- if (n_trim > 0) order(dmin, decreasing = TRUE)[seq_len(n_trim)] else integer(0)
  trimmed <- logical(nrow(X)); if (n_trim > 0) trimmed[trim_idx] <- TRUE
  list(U = base$U, centers = centers, trim_set = trimmed, alpha_used = mean(trimmed))
}

FKM.noise <- function(X, k = 3, m = 1.5, tol = 1e-4, maxit = 200, seed = 1) {
  fc <- fuzzy_c_means(X, C = k, m = m, tol = tol, maxit = maxit, seed = seed)
  list(U = fc$U, centers = fc$centers)
}

predict_fkm_noise <- function(model, noise_col = NOISE_COL, gamma = GAMMA_NOISE) {
  U <- model$U
  stopifnot(noise_col >= 1, noise_col <= ncol(U))
  is_noise <- U[, noise_col] >= gamma
  U_main   <- U[, setdiff(seq_len(ncol(U)), noise_col), drop = FALSE]
  labels_main <- max.col(U_main, ties.method = "first")
  list(labels_main = labels_main, is_noise = is_noise)
}

## =======================
## 5) Metrics
## =======================
rand_index <- function(y_true, y_pred) {
  y_true <- as.integer(y_true); y_pred <- as.integer(y_pred)
  n <- length(y_true); if (n < 2) return(NA_real_)
  M <- table(y_pred, y_true)
  tp <- sum(choose(M, 2)); P <- choose(n, 2)
  rowm <- rowSums(M); colm <- colSums(M)
  t1 <- sum(choose(rowm, 2)); t2 <- sum(choose(colm, 2))
  tn <- P - t1 - t2 + tp
  (tp + tn) / P
}

## =======================
## 6) Evaluate E / N / T for one feature matrix
## =======================
eval_variants_on_features <- function(feat_mat, y, m_grid,
                                      alpha_trim = TRIM_FRACTION,
                                      tau_dom = TAU_DOM,
                                      gamma_noise = GAMMA_NOISE,
                                      noise_col = NOISE_COL,
                                      is_outlier = NULL) {
  N <- length(y)
  if (!is.null(is_outlier)) stopifnot(length(is_outlier) == N)
  
  out_rows <- vector("list", length(m_grid) * 3L); j <- 0L
  
  for (m in m_grid) {
    ## --- E variant ---
    resE <- try(fuzzy_c_means(feat_mat, C = 2, m = m), silent = TRUE)
    if (!inherits(resE, "try-error")) {
      U    <- resE$U
      badE <- apply(U, 1, max) < tau_dom
      keep <- which(!badE)
      labE <- if (length(keep) > 1L) max.col(U[keep, , drop = FALSE], ties.method = "first") else integer(0)
      accE <- if (length(keep) > 1L) rand_index(y[keep], labE) else NA_real_
      outE <- mean(badE)
      j <- j + 1L; out_rows[[j]] <- tibble(
        variant="E", m=m, accuracy=accE, outlier_pct=outE, alpha=NA_real_,
        trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_
      )
    } else { j <- j + 1L; out_rows[[j]] <- tibble(
      variant="E", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=NA_real_,
      trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_
    )}
    
    ## --- N variant ---
    resN <- try(FKM.noise(feat_mat, k = 3, m = m), silent = TRUE)
    if (!inherits(resN, "try-error")) {
      pr   <- predict_fkm_noise(resN, noise_col = noise_col, gamma = gamma_noise)
      keep <- which(!pr$is_noise)
      accN <- if (length(keep) > 1L) rand_index(y[keep], pr$labels_main[keep]) else NA_real_
      outN <- mean(pr$is_noise)
      j <- j + 1L; out_rows[[j]] <- tibble(
        variant="N", m=m, accuracy=accN, outlier_pct=outN, alpha=NA_real_,
        trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_
      )
    } else { j <- j + 1L; out_rows[[j]] <- tibble(
      variant="N", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=NA_real_,
      trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_
    )}
    
    ## --- T variant (fixed α) ---
    resT <- try(trimmed_fuzzy_c_means(feat_mat, C = 2, m = m, alpha = alpha_trim), silent = TRUE)
    if (!inherits(resT, "try-error")) {
      U      <- resT$U
      trim_v <- resT$trim_set
      keep   <- which(!trim_v)
      
      labT <- if (length(keep) > 1L) max.col(U[keep, , drop = FALSE], ties.method = "first") else integer(0)
      accT <- if (length(keep) > 1L) rand_index(y[keep], labT) else NA_real_
      
      outT <- mean(trim_v)  # equals alpha_trim by design
      
      # New metrics vs ground truth outliers
      hit_pct   <- if (!is.null(is_outlier)) mean(is_outlier[trim_v]) else NA_real_
      rec_trim  <- if (!is.null(is_outlier) && any(is_outlier)) mean(trim_v[is_outlier]) else NA_real_
      
      acc_inl <- NA_real_
      if (!is.null(is_outlier) && length(keep) > 1L) {
        kept_inl <- keep[!is_outlier[keep]]      # non-trimmed true inliers
        if (length(kept_inl) > 1L) {
          map_idx <- match(kept_inl, keep)       # align to labT (which is on 'keep')
          acc_inl <- rand_index(y[kept_inl], labT[map_idx])
        }
      }
      
      j <- j + 1L; out_rows[[j]] <- tibble(
        variant="T", m=m, accuracy=accT, outlier_pct=outT, alpha=alpha_trim,
        trim_hit_pct=hit_pct, trim_recall=rec_trim, acc_kept_inliers=acc_inl
      )
    } else { j <- j + 1L; out_rows[[j]] <- tibble(
      variant="T", m=m, accuracy=NA_real_, outlier_pct=NA_real_, alpha=alpha_trim,
      trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_
    )}
  }
  
  bind_rows(out_rows)
}

## =======================
## 7) Checkpoints & resume helpers
## =======================
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

## =======================
## 8) Run one (family, n, p) with resume & feature cache
## =======================
run_family_np_resume <- function(family, n, p, m_grid, reps,
                                 base_seed = 2025, alpha_trim = TRIM_FRACTION,
                                 reset = FALSE, quiet = TRUE) {
  if (reset) reset_scenario(family, n, p)
  
  st <- scenario_status(family, n, p, reps)
  d  <- st$dir
  if (!quiet) cat(sprintf("[%-5s n=%d p=%d] done: %s | missing: %s\n",
                          family, n, p, toString(st$done), toString(st$missing)))
  
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
        dat <- list(
          feats = as.matrix(feats),
          labels = as.integer(sim$labels),
          is_outlier = as.logical(sim$is_outlier)
        )
        saveRDS(dat, feats_file)
      }
      
      res <- try(
        eval_variants_on_features(
          dat$feats, dat$labels, m_grid,
          alpha_trim = alpha_trim,
          tau_dom = TAU_DOM,
          gamma_noise = GAMMA_NOISE,
          noise_col = NOISE_COL,
          is_outlier = dat$is_outlier
        ),
        silent = TRUE
      )
      
      saveRDS(
        cbind(n = n, p = p, family = family,
              if (!inherits(res, "try-error")) res else
                tibble(variant=c("E","N","T"), m=NA_real_, accuracy=NA_real_,
                       outlier_pct=NA_real_, alpha=c(NA_real_, NA_real_, alpha_trim),
                       trim_hit_pct=NA_real_, trim_recall=NA_real_, acc_kept_inliers=NA_real_)),
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

## =======================
## 9) Sweep all scenarios
## =======================
run_all_features_resume <- function(families, n_vals, p_vals, m_grid, reps,
                                    alpha_trim = TRIM_FRACTION,
                                    reset_everything = FALSE, quiet = TRUE) {
  if (reset_everything) reset_all(families, n_vals, p_vals)
  
  grid  <- tidyr::expand_grid(family = families, n = n_vals, p = p_vals)
  parts <- purrr::pmap(
    grid,
    ~ run_family_np_resume(..1, ..2, ..3, m_grid, reps,
                           alpha_trim = alpha_trim,
                           reset = reset_everything,
                           quiet = quiet)
  )
  raw  <- bind_rows(Filter(Negate(is.null), parts))
  
  raw_csv <- file.path(RESULTS_DIR, "feature_replicates.csv")
  write.csv(raw, raw_csv, row.names = FALSE)
  
  means <- raw %>%
    group_by(family, n, p, variant, m) %>%
    summarise(
      accuracy          = mean(accuracy, na.rm = TRUE),
      acc_kept_inliers  = mean(acc_kept_inliers, na.rm = TRUE),
      outlier_pct       = mean(outlier_pct, na.rm = TRUE),
      trim_hit_pct      = mean(trim_hit_pct, na.rm = TRUE),
      trim_recall       = mean(trim_recall, na.rm = TRUE),
      alpha             = if (all(is.na(alpha))) NA_real_ else mean(alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(family, n, p, variant, m)
  
  write.csv(means, file.path(RESULTS_DIR, "feature_means.csv"), row.names = FALSE)
  means
}

## =======================
## 10) Run
## =======================
MEANS <- run_all_features_resume(
  families, n_vals, p_vals, m_grid, reps,
  alpha_trim = TRIM_FRACTION,
  reset_everything = FALSE,   # set TRUE to wipe this RUN_TAG's caches
  quiet = TRUE
)

print(MEANS, n = Inf)
cat("\nReplicate CSV:", file.path(RESULTS_DIR, "feature_replicates.csv"),
    "\nMeans CSV:     ", file.path(RESULTS_DIR, "feature_means.csv"), "\n")
