#!/usr/bin/env Rscript
# Published comparator: group-penalized Cox regression (grpreg::grpsurv;
# Breheny & Huang, 2015, Statistics and Computing 25:173-187), applied as in
# its documentation: group MCP on a cubic B-spline expansion of every additive
# candidate plus singleton groups for the linear candidates, with the
# penalty level chosen by the package's own 10-fold cross-validation
# (cv.grpsurv, lambda.min). Failure margin only; no screening.
#
# Usage:
#   Rscript grpsurv_comparator.R <scenario_data_dir> <first_file> <last_file> <out_csv>
# Data files are data_seed<seed>.csv.gz written by the MATLAB engine with
# columns: time, status, Z_1..Z_p, X_1..X_q (p = q).

suppressPackageStartupMessages({
  library(grpreg)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: grpsurv_comparator.R data_dir first last out_csv")
data_dir <- args[1]
first <- as.integer(args[2])
last <- as.integer(args[3])
out_csv <- args[4]

files <- sort(list.files(data_dir, pattern = "^data_seed[0-9]+\\.csv\\.gz$", full.names = TRUE))
seeds <- as.numeric(sub("^data_seed([0-9]+)\\.csv\\.gz$", "\\1", basename(files)))
files <- files[order(seeds)]
seeds <- sort(seeds)
if (last > length(files)) last <- length(files)
if (first > last) stop("empty file range")

# Truth (failure margin) for the discriminating / shared_balanced_bounded DGP
beta_true <- c(1.5, -1.2, 0.9); beta_active <- c(1, 3, 6)
add_active <- 1:4
g_true <- list(
  function(x) 1.6 * x,
  function(x) sqrt(10) * (6 * x^2 - 6 * x + 1),
  function(x) 2 * sin(2 * pi * x),
  function(x) 2 * cos(2 * pi * x))
grid_for <- function(j) if (j == 1) seq(-2, 2, length.out = 101) else seq(0, 1, length.out = 101)
trap_w <- function(m) { w <- rep(1, m); w[c(1, m)] <- 0.5; w / sum(w) }

interp_extrap <- function(x, y, xout) {
  o <- order(x); x <- x[o]; y <- y[o]
  ux <- unique(x); uy <- tapply(y, x, mean)[as.character(ux)]
  if (length(ux) < 2) return(rep(NA_real_, length(xout)))
  res <- approx(ux, uy, xout = xout, rule = 2)$y
  lo <- xout < ux[1]; hi <- xout > ux[length(ux)]
  if (any(lo)) res[lo] <- uy[1] + (xout[lo] - ux[1]) * (uy[2] - uy[1]) / (ux[2] - ux[1])
  k <- length(ux)
  if (any(hi)) res[hi] <- uy[k] + (xout[hi] - ux[k]) * (uy[k] - uy[k - 1]) / (ux[k] - ux[k - 1])
  res
}

spline_block <- function(x) {
  kn <- quantile(x, probs = seq(0, 1, length.out = 6), type = 7)
  inner <- kn[2:5]
  B <- bs(x, knots = inner, degree = 3, Boundary.knots = range(x), intercept = FALSE)
  B <- unclass(B)[, , drop = FALSE]
  sweep(B, 2, colMeans(B))
}

rows <- vector("list", last - first + 1)
curve_rows <- vector("list", last - first + 1)
for (ii in first:last) {
  t0 <- proc.time()[["elapsed"]]
  D <- as.matrix(read.csv(gzfile(files[ii]), header = FALSE))
  time <- D[, 1]; status <- D[, 2]
  pq <- (ncol(D) - 2) / 2
  Z <- D[, 2 + seq_len(pq), drop = FALSE]
  X <- D[, 2 + pq + seq_len(pq), drop = FALSE]
  blocks <- lapply(seq_len(pq), function(j) spline_block(X[, j]))
  k <- ncol(blocks[[1]])
  design <- cbind(Z, do.call(cbind, blocks))
  group <- c(seq_len(pq), rep(pq + seq_len(pq), each = k))
  set.seed(seeds[ii])
  cvfit <- cv.grpsurv(design, cbind(time, status), group = group,
                      penalty = "grMCP", nfolds = 10)
  b <- as.numeric(coef(cvfit$fit, lambda = cvfit$lambda.min))  # no intercept in Cox models
  b_lin <- b[seq_len(pq)]
  b_add <- matrix(b[pq + seq_len(pq * k)], nrow = k)
  lin_sel <- which(b_lin != 0)
  add_sel <- which(colSums(b_add != 0) > 0)
  # Function error for the four active failure-margin components
  ise <- numeric(4)
  cr <- vector("list", 4)
  for (j in add_active) {
    fit_obs <- blocks[[j]] %*% b_add[, j]
    grid <- grid_for(j)
    w <- trap_w(length(grid))
    fit_grid <- interp_extrap(X[, j], as.numeric(fit_obs), grid)
    truth <- g_true[[j]](grid)
    fit_grid <- fit_grid - sum(w * fit_grid)
    truth <- truth - sum(w * truth)
    ise[j] <- sum(w * (fit_grid - truth)^2)
    cr[[j]] <- data.frame(seed = seeds[ii], component = j, x = grid, fit = fit_grid)
  }
  curve_rows[[ii - first + 1]] <- do.call(rbind, cr)
  rows[[ii - first + 1]] <- data.frame(
    seed = seeds[ii], n = nrow(D), p = pq, events = sum(status),
    lambda_min = cvfit$lambda.min,
    lin_selected = length(lin_sel),
    lin_fp = sum(!(lin_sel %in% beta_active)),
    lin_fn = sum(!(beta_active %in% lin_sel)),
    add_selected = length(add_sel),
    add_fp = sum(!(add_sel %in% add_active)),
    add_fn = sum(!(add_active %in% add_sel)),
    beta1 = b_lin[1], beta3 = b_lin[3], beta6 = b_lin[6],
    mse_beta_active = mean((b_lin[beta_active] - beta_true)^2),
    imse_add = mean(ise),
    ise_g1 = ise[1], ise_g2 = ise[2], ise_g3 = ise[3], ise_g4 = ise[4],
    seconds = proc.time()[["elapsed"]] - t0)
  cat(sprintf("seed %d done in %.1f s\n", seeds[ii], rows[[ii - first + 1]]$seconds))
}
out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
# Centered fitted curves on the common grid (for ISB^2 / IVar summaries)
write.csv(do.call(rbind, curve_rows), sub("\\.csv$", "_curves.csv", out_csv), row.names = FALSE)
