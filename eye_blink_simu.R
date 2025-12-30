## =========================================================
## RFCPCA (Eye-Blink) — Variable T per MTS (400..2000)
## Parallel + Checkpoint + Resume  |  40% contamination
## Re-runnable via FORCE_RERUN or RUN_ID
## =========================================================

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(progressr)
  library(R.utils)
  library(stats)
})

## ---------- Controls ----------
FORCE_RERUN <- TRUE           # set TRUE to re-run even if checkpoints exist
RUN_ID      <- ""             # set like "rerun1" to write to new folders

## ---------- Threading hygiene ----------
Sys.setenv(
  OMP_NUM_THREADS      = 1,
  MKL_NUM_THREADS      = 1,
  OPENBLAS_NUM_THREADS = 1,
  NUMEXPR_NUM_THREADS  = 1
)

## ---------- Parallel plan ----------
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(workers) || workers < 1) workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
handlers(global = TRUE)

## ---------- Paths ----------
ROOT <- Sys.getenv("SCRATCH", ".")
tag  <- if (nzchar(RUN_ID)) paste0("_", RUN_ID) else ""
BASE_CHKPT   <- file.path(ROOT, paste0("rfcpca_chkpts_blink_varlen_noalpha", tag))
BASE_RESULTS <- file.path(ROOT, paste0("rfcpca_results_blink_varlen_noalpha", tag))
dir.create(BASE_CHKPT,  showWarnings = FALSE, recursive = TRUE)
dir.create(BASE_RESULTS, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1) AR(2) EEG-like generator (5 bands) + Eye-blink contamination
## =========================================================

fband_arcoefs <- function(fband, samp.rate = 100) {
  if (fband == "delta") { center <-  2; sharp <- 0.05 }
  else if (fband == "theta") { center <-  6; sharp <- 0.05 }
  else if (fband == "alpha") { center <- 10; sharp <- 0.05 }
  else if (fband == "beta")  { center <- 22; sharp <- 0.08 }
  else if (fband == "gamma") { center <- 37; sharp <- 0.10 }
  else stop("Unknown band: ", fband)
  M <- exp(sharp)                  # pole radius
  psi <- center / samp.rate
  c((2 / M) * cos(2 * pi * psi), -(1 / M^2))
}

gen_latent <- function(fband, Tlen, fs = 100) {
  phis <- fband_arcoefs(fband, fs)
  z <- as.numeric(arima.sim(n = Tlen, model = list(ar = phis)))
  as.numeric(scale(z, center = TRUE, scale = TRUE))
}

gen_coef <- function(idx_major, p_major_sum) {
  stopifnot(all(idx_major %in% 1:5), p_major_sum > 0, p_major_sum < 1)
  w <- numeric(5); rem <- p_major_sum
  for (k in seq_along(idx_major)) {
    i <- idx_major[k]
    if (k < length(idx_major)) {
      w[i] <- runif(1, (p_major_sum/length(idx_major)) * 0.9,
                    (p_major_sum/length(idx_major)) * 1.1)
      rem <- rem - w[i]
    } else w[i] <- max(rem, 0)
  }
  idx_minor <- setdiff(1:5, idx_major)
  if (length(idx_minor)) {
    rem2 <- 1 - sum(w)
    for (k in seq_along(idx_minor)) {
      i <- idx_minor[k]
      if (k < length(idx_minor)) {
        w[i] <- runif(1, (rem2/length(idx_minor)) * 0.9,
                      (rem2/length(idx_minor)) * 1.1)
        rem2 <- rem2 - w[i]
      } else w[i] <- max(rem2, 0)
    }
  }
  w / sum(w)
}

add_eye_blinks <- function(X,
                           fs            = 100,
                           blinks_range  = 1:2,
                           width_ms      = c(200, 400),
                           ch_frac       = c(0.15, 0.30),
                           amp_factor    = c(4, 8),
                           frontal_idx   = NULL,
                           sign_random   = TRUE) {
  Tlen <- nrow(X); p <- ncol(X)
  if (is.null(frontal_idx)) frontal_idx <- seq_len(max(1, floor(0.25 * p)))
  nb <- if (length(blinks_range) > 1) sample(blinks_range, 1) else blinks_range
  for (b in seq_len(nb)) {
    L  <- round(runif(1, width_ms[1], width_ms[2]) * fs / 1000)
    L  <- max(10L, min(L, Tlen - 2L))
    st <- sample.int(Tlen - L + 1L, 1)
    k  <- max(1L, round(runif(1, ch_frac[1], ch_frac[2]) * length(frontal_idx)))
    ch <- sample(frontal_idx, k)
    t  <- 0:(L - 1)
    bump <- sin(pi * t / (L - 1))   # half-sine (0 → π)
    sgn  <- if (sign_random) sample(c(-1,1), 1) else 1
    for (j in ch) {
      A <- runif(1, amp_factor[1], amp_factor[2]) * sd(X[, j])
      X[st:(st + L - 1), j] <- X[st:(st + L - 1), j] + sgn * A * bump
    }
  }
  X
}

make_trial <- function(W, Tlen, fs = 100, noise_sd = 0) {
  bands <- c("delta","theta","alpha","beta","gamma")
  Z5 <- t(vapply(bands, function(b) gen_latent(b, Tlen, fs), numeric(Tlen)))  # 5 x T
  L  <- scale(t(Z5), center = TRUE, scale = TRUE)                             # T x 5
  X  <- L %*% W                                                               # T x p
  if (noise_sd > 0) X <- X + matrix(rnorm(length(X), 0, noise_sd), nrow = nrow(X))
  X
}

simulate_two_groups_mts_varlen <- function(p, fs = 100,
                                           n_per_group = 10,
                                           contam_prop = 0.40,
                                           contam_type = "blink",
                                           T_min = 400, T_max = 2000,
                                           T_sampler = c("uniform","grid"),
                                           grid_vals = NULL,
                                           blink_args = list(),
                                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  T_sampler <- match.arg(T_sampler)
  draw_T <- switch(
    T_sampler,
    "uniform" = function() sample(T_min:T_max, 1),
    "grid"    = function() {
      if (is.null(grid_vals) || !length(grid_vals))
        stop("grid_vals must be provided for T_sampler='grid'")
      sample(grid_vals, 1)
    }
  )
  W1 <- sapply(1:p, \(.) gen_coef(idx_major = c(1,3,5), p_major_sum = 0.95))
  W2 <- sapply(1:p, \(.) gen_coef(idx_major = c(2,4),   p_major_sum = 0.95))
  G1 <- vector("list", n_per_group); G2 <- vector("list", n_per_group)
  T1 <- integer(n_per_group); T2 <- integer(n_per_group)
  for (i in seq_len(n_per_group)) {
    T1[i] <- draw_T(); T2[i] <- draw_T()
    G1[[i]] <- make_trial(W1, T1[i], fs)
    G2[[i]] <- make_trial(W2, T2[i], fs)
  }
  kpg  <- ceiling(contam_prop * n_per_group)
  idx1 <- if (kpg > 0) sample.int(n_per_group, kpg) else integer(0)
  idx2 <- if (kpg > 0) sample.int(n_per_group, kpg) else integer(0)
  if (tolower(contam_type) == "blink") {
    for (i in idx1) G1[[i]] <- do.call(add_eye_blinks, c(list(X = G1[[i]], fs = fs), blink_args))
    for (i in idx2) G2[[i]] <- do.call(add_eye_blinks, c(list(X = G2[[i]], fs = fs), blink_args))
  } else stop("Unknown contam_type: ", contam_type)
  ts_list    <- c(G1, G2)
  labels     <- c(rep(1L, n_per_group), rep(2L, n_per_group))
  is_outlier <- rep(FALSE, length(ts_list))
  if (length(idx1)) is_outlier[idx1] <- TRUE
  if (length(idx2)) is_outlier[n_per_group + idx2] <- TRUE
  list(ts = ts_list, labels = labels, is_outlier = is_outlier,
       T_lengths = c(T1, T2))
}

## =========================================================
## 2) Evaluation helpers
## =========================================================

THETA_CRISP <- 0.70
NOISE_THR   <- 0.50

hard_from_U <- function(U, theta = THETA_CRISP) {
  crisp <- apply(U, 1, max)
  lab   <- max.col(U, ties.method = "first")
  lab[crisp < theta] <- NA_integer_
  list(pred = lab, crisp = crisp)
}

acc_best2 <- function(y_true, pred) {
  ok <- !is.na(pred); if (!any(ok)) return(NA_real_)
  y <- as.integer(y_true[ok]); p <- as.integer(pred[ok])
  C <- table(factor(p, levels = sort(unique(p))),
             factor(y, levels = sort(unique(y))))
  if (nrow(C) < 2 || ncol(C) < 2) return(mean(y == p))
  max((C[1,1] + C[2,2]), (C[1,2] + C[2,1])) / length(y)
}

as_trim_idx <- function(trim_set, N) {
  if (is.null(trim_set)) return(integer(0))
  if (is.logical(trim_set)) return(which(trim_set))
  if (is.numeric(trim_set) && length(trim_set) == N && all(trim_set %in% c(0,1)))
    return(which(trim_set == 1))
  as.integer(trim_set)
}

## safe mean to avoid NaN in summaries
safe_mean <- function(v) { v <- v[!is.na(v)]; if (length(v)==0) NA_real_ else mean(v) }

## One replicate
one_rep <- function(sim) {
  ts <- sim$ts; y <- sim$labels; out_true <- which(sim$is_outlier)
  acc_FCPCA <- out_FCPCA <- NA_real_
  acc_E     <- out_E     <- NA_real_
  acc_N     <- det_N     <- NA_real_
  acc_T     <- det_T     <- alpha_T <- NA_real_
  
  ## --- FCPCA (using your current function name) ---
  fitA <- try(fcpca_auto(ts, 2), silent = TRUE)   # change to FCPCA(ts,2) if needed
  if (!inherits(fitA, "try-error")) {
    U <- fitA$membership_matrix
    hf <- hard_from_U(U)
    pred_out <- hf$crisp < THETA_CRISP
    keep <- which(!pred_out)
    acc_FCPCA <- if (length(keep)) acc_best2(y[keep], hf$pred[keep]) else NA_real_
    out_FCPCA <- if (length(out_true)) mean(pred_out[out_true]) else NA_real_
  }
  
  ## --- RFCPCA-E ---
  fitE <- try(RFCPCA_E(ts, 2), silent = TRUE)
  if (!inherits(fitE, "try-error")) {
    U <- fitE$membership_matrix
    hf <- hard_from_U(U)
    pred_out <- hf$crisp < THETA_CRISP
    keep <- which(!pred_out)
    acc_E <- if (length(keep)) acc_best2(y[keep], hf$pred[keep]) else NA_real_
    out_E <- if (length(out_true)) mean(pred_out[out_true]) else NA_real_
  }
  
  ## --- RFCPCA-N ---
  fitN <- try(RFCPCA_N(ts, 2), silent = TRUE)
  if (!inherits(fitN, "try-error")) {
    U  <- fitN$membership_matrix
    nz <- grep("noise|u_noise", colnames(U), ignore.case = TRUE)
    if (!length(nz)) nz <- 3L
    nz <- nz[1]
    pred_out <- U[, nz] >= NOISE_THR
    keep <- which(!pred_out)
    Uc  <- U[, setdiff(seq_len(ncol(U)), nz), drop = FALSE]
    hf  <- hard_from_U(Uc)
    acc_N <- if (length(keep)) acc_best2(y[keep], hf$pred[keep]) else NA_real_
    det_N <- if (length(out_true)) mean(pred_out[out_true]) else NA_real_
  }
  
  ## --- RFCPCA-T (no alpha passed) ---
  fitT <- try(RFCPCA_T(ts, 2), silent = TRUE)
  if (!inherits(fitT, "try-error")) {
    U <- fitT$membership_matrix
    trim_idx <- as_trim_idx(fitT$trim_set, nrow(U))
    keep <- setdiff(seq_len(nrow(U)), trim_idx)
    hf <- hard_from_U(U)
    acc_T   <- if (length(keep)) acc_best2(y[keep], hf$pred[keep]) else NA_real_
    det_T   <- if (length(out_true)) mean(out_true %in% trim_idx) else NA_real_
    alpha_T <- if (!is.null(fitT$alpha_used)) fitT$alpha_used else length(trim_idx)/nrow(U)
  }
  
  c(acc_FCPCA, out_FCPCA,
    acc_E,     out_E,
    acc_N,     det_N,
    acc_T,     det_T, alpha_T)
}

## =========================================================
## 3) Checkpoint + Resume helpers (no alpha dimension)
## =========================================================

range_tag <- function(T_min, T_max) sprintf("Trange%04d_%04d", as.integer(T_min), as.integer(T_max))
sc_dir    <- function(p, T_min, T_max)
  file.path(BASE_CHKPT, range_tag(T_min, T_max), sprintf("p%03d", p))
rep_done_path <- function(d, r) file.path(d, sprintf("rep_%03d_done.rds", r))

scenario_status <- function(p, reps, T_min, T_max) {
  d <- sc_dir(p, T_min, T_max); dir.create(d, FALSE, TRUE)
  files <- list.files(d, "^rep_\\d+_done\\.rds$", full.names = FALSE)
  done  <- sort(as.integer(sub("^rep_(\\d+)_done\\.rds$", "\\1", files)))
  list(dir = d, done = done, missing = setdiff(seq_len(reps), done))
}

## =========================================================
## 4) Runner for one (p, Trange) — interrupt-safe
## =========================================================

run_scenario_means_resume_varlen <- function(p,
                                             reps = 50,
                                             T_min = 400, T_max = 2000,
                                             T_sampler = c("uniform","grid"),
                                             grid_vals = NULL,
                                             base_seed = 1000,
                                             timeout_sec = 1800,
                                             contam_prop = 0.40,
                                             force_rerun = FORCE_RERUN) {
  T_sampler <- match.arg(T_sampler)
  st <- scenario_status(p, reps, T_min, T_max); d <- st$dir
  
  ## if forcing rerun, remove existing rep_*** files for this scenario
  if (isTRUE(force_rerun) && length(st$done)) {
    old_files <- file.path(d, sprintf("rep_%03d_done.rds", st$done))
    file.remove(old_files)
    st <- scenario_status(p, reps, T_min, T_max)  # recompute
  }
  
  blink_args <- list(
    blinks_range = 1:2,
    width_ms     = c(200, 400),
    ch_frac      = c(0.15, 0.30),
    amp_factor   = c(4, 8),
    sign_random  = TRUE
  )
  
  if (length(st$missing) > 0) {
    message(sprintf("[p=%d, %d..%d] running %d/%d missing reps ...",
                    p, T_min, T_max, length(st$missing), reps))
    run_one <- function(r) {
      res <- try({
        sim <- simulate_two_groups_mts_varlen(
          p = p, fs = 100,
          n_per_group = 10,
          contam_prop = contam_prop,
          contam_type = "blink",
          T_min = T_min, T_max = T_max,
          T_sampler = T_sampler,
          grid_vals = grid_vals,
          blink_args = blink_args,
          seed = base_seed + r
        )
        withTimeout(one_rep(sim), timeout = timeout_sec, onTimeout = "warning")
      }, silent = TRUE)
      saveRDS(res, rep_done_path(d, r))
      TRUE
    }
    with_progress({
      pbar <- progressor(steps = length(st$missing))
      future_lapply(st$missing,
                    function(r) { on.exit(pbar(sprintf("rep %d", r))); run_one(r) },
                    future.seed = TRUE)
    })
  } else {
    message(sprintf("[p=%d, %d..%d] all %d reps already complete.",
                    p, T_min, T_max, reps))
  }
  
  files <- file.path(d, sprintf("rep_%03d_done.rds", 1:reps))
  files <- files[file.exists(files)]
  stopifnot(length(files) > 0)
  
  Mlist <- lapply(files, readRDS)
  Mlist <- Filter(function(x) is.numeric(x) && length(x) == 9, Mlist)
  if (!length(Mlist)) stop("No valid replicate results found for this scenario.")
  M <- do.call(rbind, Mlist)
  
  MEANS <- data.frame(
    p = p, T_min = T_min, T_max = T_max,
    method      = c("FCPCA", "RFCPCA_E", "RFCPCA_N", "RFCPCA_T"),
    accuracy    = c(safe_mean(M[,1]), safe_mean(M[,3]), safe_mean(M[,5]), safe_mean(M[,7])),
    outlier_pct = c(safe_mean(M[,2]), safe_mean(M[,4]), safe_mean(M[,6]), safe_mean(M[,8])),
    alpha_used  = c(NA, NA, NA, safe_mean(M[,9]))
  )
  num_cols <- sapply(MEANS, is.numeric)
  MEANS[num_cols] <- lapply(MEANS[num_cols], function(v) round(v, 3))
  
  out_dir <- file.path(BASE_RESULTS, range_tag(T_min, T_max))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(out_dir, sprintf("means_p%d.csv", p))
  write.csv(MEANS, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
  MEANS
}

## =========================================================
## 5) Sweep p; aggregate; resume-friendly
## =========================================================

run_all_eye_blink_varlen <- function(p_vals = c(32, 64, 128),
                                     reps = 50,
                                     T_min = 400, T_max = 2000,
                                     T_sampler = c("uniform","grid"),
                                     grid_vals = NULL,
                                     contam_prop = 0.40,
                                     force_rerun = FORCE_RERUN) {
  T_sampler <- match.arg(T_sampler)
  parts <- lapply(p_vals, function(p_) {
    run_scenario_means_resume_varlen(
      p = p_, reps = reps,
      T_min = T_min, T_max = T_max,
      T_sampler = T_sampler,
      grid_vals = grid_vals,
      contam_prop = contam_prop,
      force_rerun = force_rerun
    )
  })
  OUT <- do.call(rbind, parts)
  OUT$contam_prop <- contam_prop
  out_csv <- file.path(BASE_RESULTS, sprintf("all_means_varlen_%s.csv", T_sampler))
  write.csv(OUT, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
  OUT
}
 
OUT <- run_all_eye_blink_varlen(
  p_vals = c(32, 64, 128),
  reps = 50,
  T_min = 400, T_max = 2000,
  T_sampler = "uniform",
  contam_prop = 0.40,
  force_rerun = FORCE_RERUN
)
print(OUT)
