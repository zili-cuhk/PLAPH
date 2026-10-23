suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

# Folder holding luad_pilot_clinical.csv, Z_RNA_var_top2000_os.csv and X_log2CNA_var_top1000_os.csv
root <- Sys.getenv("PLAPH_LUAD_ROOT", "plaph_luad_noimage_pilot_screening")
out_dir <- file.path(root, "pilot_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_grid <- c(300, 500, 800, 1000)
q_grid <- c(50, 100, 200)
spline_basis_k <- 5

stage_group <- function(stage) {
  out <- rep(NA_character_, length(stage))
  out[grepl("STAGE I", stage, fixed = TRUE)] <- "I"
  out[grepl("STAGE II", stage, fixed = TRUE)] <- "II"
  out[grepl("STAGE III", stage, fixed = TRUE)] <- "III"
  out[grepl("STAGE IV", stage, fixed = TRUE)] <- "IV"
  factor(out, levels = c("I", "II", "III", "IV"))
}

neglog10p <- function(p) {
  p <- ifelse(is.finite(p), pmax(p, .Machine$double.xmin), NA_real_)
  -log10(p)
}

cox_lrt_p <- function(fit_small, fit_large) {
  ll0 <- as.numeric(logLik(fit_small))
  ll1 <- as.numeric(logLik(fit_large))
  df0 <- attr(logLik(fit_small), "df")
  df1 <- attr(logLik(fit_large), "df")
  chisq <- 2 * (ll1 - ll0)
  df <- df1 - df0
  if (!is.finite(chisq) || chisq < 0 || df <= 0) return(NA_real_)
  pchisq(chisq, df = df, lower.tail = FALSE)
}

adjusted_linear_screen <- function(mat, clinical, time_col, event_col, label) {
  features <- setdiff(names(mat), "patient_id")
  rows <- vector("list", length(features))
  for (i in seq_along(features)) {
    f <- features[i]
    d <- data.frame(
      time = clinical[[time_col]],
      event = clinical[[event_col]],
      age = clinical$age,
      sex = factor(clinical$sex),
      stage_group = clinical$stage_group,
      x = mat[[f]]
    )
    d <- d[is.finite(d$time) & d$time > 0 & is.finite(d$event) &
             is.finite(d$age) & !is.na(d$sex) & !is.na(d$stage_group) & is.finite(d$x), ]
    if (nrow(d) < 100 || sum(d$event) < 20 || length(unique(d$x)) < 3) {
      rows[[i]] <- data.frame(feature = f, beta = NA_real_, hr = NA_real_, p = NA_real_,
                              n = nrow(d), events = sum(d$event), outcome = label)
      next
    }
    fit <- try(coxph(Surv(time, event) ~ x + age + sex + stage_group, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) {
      rows[[i]] <- data.frame(feature = f, beta = NA_real_, hr = NA_real_, p = NA_real_,
                              n = nrow(d), events = sum(d$event), outcome = label)
    } else {
      co <- summary(fit)$coefficients
      rows[[i]] <- data.frame(feature = f, beta = unname(co["x", "coef"]),
                              hr = unname(co["x", "exp(coef)"]),
                              p = unname(co["x", "Pr(>|z|)"]),
                              n = nrow(d), events = sum(d$event), outcome = label)
    }
  }
  ans <- do.call(rbind, rows)
  ans$fdr <- p.adjust(ans$p, method = "BH")
  ans[order(ans$p, na.last = TRUE), ]
}

adjusted_spline_screen <- function(mat, clinical, time_col, event_col, label) {
  features <- setdiff(names(mat), "patient_id")
  rows <- vector("list", length(features))
  for (i in seq_along(features)) {
    f <- features[i]
    d <- data.frame(
      time = clinical[[time_col]],
      event = clinical[[event_col]],
      age = clinical$age,
      sex = factor(clinical$sex),
      stage_group = clinical$stage_group,
      x = mat[[f]]
    )
    d <- d[is.finite(d$time) & d$time > 0 & is.finite(d$event) &
             is.finite(d$age) & !is.na(d$sex) & !is.na(d$stage_group) & is.finite(d$x), ]
    if (nrow(d) < 100 || sum(d$event) < 20 || length(unique(d$x)) < 6) {
      rows[[i]] <- data.frame(feature = f, p_overall = NA_real_, p_nonlinear = NA_real_,
                              n = nrow(d), events = sum(d$event), outcome = label)
      next
    }
    fit0 <- try(coxph(Surv(time, event) ~ age + sex + stage_group, data = d), silent = TRUE)
    fit1 <- try(coxph(Surv(time, event) ~ x + age + sex + stage_group, data = d), silent = TRUE)
    fit3 <- try(coxph(Surv(time, event) ~ ns(x, df = 3) + age + sex + stage_group, data = d), silent = TRUE)
    if (inherits(fit0, "try-error") || inherits(fit1, "try-error") || inherits(fit3, "try-error")) {
      rows[[i]] <- data.frame(feature = f, p_overall = NA_real_, p_nonlinear = NA_real_,
                              n = nrow(d), events = sum(d$event), outcome = label)
    } else {
      rows[[i]] <- data.frame(feature = f,
                              p_overall = cox_lrt_p(fit0, fit3),
                              p_nonlinear = cox_lrt_p(fit1, fit3),
                              n = nrow(d), events = sum(d$event), outcome = label)
    }
  }
  ans <- do.call(rbind, rows)
  ans$fdr_overall <- p.adjust(ans$p_overall, method = "BH")
  ans$fdr_nonlinear <- p.adjust(ans$p_nonlinear, method = "BH")
  ans[order(ans$p_overall, na.last = TRUE), ]
}

combine_linear_scores <- function(event_tab, censor_tab) {
  e <- event_tab[, c("feature", "p", "fdr", "beta", "hr")]
  c <- censor_tab[, c("feature", "p", "fdr", "beta", "hr")]
  names(e) <- c("feature", "p_event", "fdr_event", "beta_event", "hr_event")
  names(c) <- c("feature", "p_censor", "fdr_censor", "beta_censor", "hr_censor")
  out <- merge(e, c, by = "feature", all = TRUE)
  out$score_event <- neglog10p(out$p_event)
  out$score_censor <- neglog10p(out$p_censor)
  out$combined_score <- pmax(out$score_event, out$score_censor, na.rm = TRUE)
  out$combined_score[!is.finite(out$combined_score)] <- NA_real_
  out <- out[order(-out$combined_score, out$p_event, out$p_censor, na.last = TRUE), ]
  out$combined_rank <- seq_len(nrow(out))
  out
}

combine_spline_scores <- function(event_tab, censor_tab) {
  e <- event_tab[, c("feature", "p_overall", "p_nonlinear", "fdr_overall", "fdr_nonlinear")]
  c <- censor_tab[, c("feature", "p_overall", "p_nonlinear", "fdr_overall", "fdr_nonlinear")]
  names(e) <- c("feature", "p_event_overall", "p_event_nonlinear", "fdr_event_overall", "fdr_event_nonlinear")
  names(c) <- c("feature", "p_censor_overall", "p_censor_nonlinear", "fdr_censor_overall", "fdr_censor_nonlinear")
  out <- merge(e, c, by = "feature", all = TRUE)
  out$score_event <- neglog10p(out$p_event_overall)
  out$score_censor <- neglog10p(out$p_censor_overall)
  out$combined_score <- pmax(out$score_event, out$score_censor, na.rm = TRUE)
  out$combined_score[!is.finite(out$combined_score)] <- NA_real_
  out <- out[order(-out$combined_score, out$p_event_overall, out$p_censor_overall, na.last = TRUE), ]
  out$combined_rank <- seq_len(nrow(out))
  out
}

clinical <- read.csv(file.path(root, "luad_pilot_clinical.csv"), check.names = FALSE)
z <- read.csv(file.path(root, "Z_RNA_var_top2000_os.csv"), check.names = FALSE)
x <- read.csv(file.path(root, "X_log2CNA_var_top1000_os.csv"), check.names = FALSE)
stopifnot(identical(clinical$patient_id, z$patient_id), identical(clinical$patient_id, x$patient_id))

clinical$stage_group <- stage_group(clinical$stage)
clinical$os_censoring_indicator <- 1 - clinical$os_event
usable <- is.finite(clinical$os_months) & clinical$os_months > 0 & !is.na(clinical$stage_group)

cat("Running adjusted RNA event screen...\n")
rna_event <- adjusted_linear_screen(z, clinical, "os_months", "os_event", "OS event")
cat("Running adjusted RNA censoring screen...\n")
rna_censor <- adjusted_linear_screen(z, clinical, "os_months", "os_censoring_indicator", "OS censoring")
cat("Running adjusted CNA event spline screen...\n")
cna_event <- adjusted_spline_screen(x, clinical, "os_months", "os_event", "OS event")
cat("Running adjusted CNA censoring spline screen...\n")
cna_censor <- adjusted_spline_screen(x, clinical, "os_months", "os_censoring_indicator", "OS censoring")

rna_rank <- combine_linear_scores(rna_event, rna_censor)
cna_rank <- combine_spline_scores(cna_event, cna_censor)

write.csv(rna_event, file.path(out_dir, "rna_adjusted_event_screen_top2000.csv"), row.names = FALSE)
write.csv(rna_censor, file.path(out_dir, "rna_adjusted_censor_screen_top2000.csv"), row.names = FALSE)
write.csv(rna_rank, file.path(out_dir, "rna_adjusted_combined_rank_top2000.csv"), row.names = FALSE)
write.csv(cna_event, file.path(out_dir, "cna_adjusted_event_spline_screen_top1000.csv"), row.names = FALSE)
write.csv(cna_censor, file.path(out_dir, "cna_adjusted_censor_spline_screen_top1000.csv"), row.names = FALSE)
write.csv(cna_rank, file.path(out_dir, "cna_adjusted_combined_rank_top1000.csv"), row.names = FALSE)

grid_rows <- list()
idx <- 1
for (p in p_grid) {
  rna_keep <- head(rna_rank$feature, p)
  rna_event_top50 <- head(rna_event$feature, 50)
  rna_censor_top50 <- head(rna_censor$feature, 50)
  for (q in q_grid) {
    cna_keep <- head(cna_rank$feature, q)
    cna_event_top50 <- head(cna_event$feature, 50)
    cna_censor_top50 <- head(cna_censor$feature, 50)
    grid_rows[[idx]] <- data.frame(
      p = p,
      q = q,
      n = sum(usable),
      total_model_dimension_approx = p + q * spline_basis_k,
      p_gt_n = p > sum(usable),
      qk_gt_n = q * spline_basis_k > sum(usable),
      rna_top50_event_retained = sum(rna_event_top50 %in% rna_keep),
      rna_top50_censor_retained = sum(rna_censor_top50 %in% rna_keep),
      rna_event_censor_top50_overlap = length(intersect(rna_event_top50, rna_censor_top50)),
      cna_top50_event_retained = sum(cna_event_top50 %in% cna_keep),
      cna_top50_censor_retained = sum(cna_censor_top50 %in% cna_keep),
      cna_event_censor_top50_overlap = length(intersect(cna_event_top50, cna_censor_top50)),
      rna_event_fdr05_retained = sum(rna_event$feature[rna_event$fdr < 0.05] %in% rna_keep, na.rm = TRUE),
      rna_censor_fdr05_retained = sum(rna_censor$feature[rna_censor$fdr < 0.05] %in% rna_keep, na.rm = TRUE),
      cna_event_fdr10_retained = sum(cna_event$feature[cna_event$fdr_overall < 0.10] %in% cna_keep, na.rm = TRUE),
      cna_censor_fdr10_retained = sum(cna_censor$feature[cna_censor$fdr_overall < 0.10] %in% cna_keep, na.rm = TRUE)
    )
    idx <- idx + 1
  }
}
grid <- do.call(rbind, grid_rows)

grid$pilot_recommendation <- "sensitivity"
grid$pilot_recommendation[grid$p == 500 & grid$q == 100] <- "stable_fallback_or_sensitivity"
grid$pilot_recommendation[grid$p == 500 & grid$q == 200] <- "primary_candidate_if_PLAPH_converges"
grid$pilot_recommendation[grid$p == 1000 & grid$q == 100] <- "linear_high_dimension_sensitivity"

write.csv(grid, file.path(out_dir, "pilot_dimension_grid_summary.csv"), row.names = FALSE)

summary_tab <- data.frame(
  quantity = c("n_total", "n_with_stage_I_to_IV_for_adjusted_screening",
               "OS_events_used", "OS_censored_used",
               "RNA_unsupervised_dimension", "CNA_unsupervised_dimension",
               "spline_basis_K_for_dimension_accounting",
               "recommended_primary_dimension",
               "recommended_sensitivity_dimensions"),
  value = c(nrow(clinical), sum(usable), sum(clinical$os_event[usable]),
            sum(clinical$os_censoring_indicator[usable]),
            length(setdiff(names(z), "patient_id")), length(setdiff(names(x), "patient_id")),
            spline_basis_k,
            "p=500, q=200 if proposed PLAPH is numerically stable",
            "p=500, q=100 as stable fallback; p=1000, q=100 as linear high-dimensional sensitivity")
)
write.csv(summary_tab, file.path(out_dir, "pilot_dimension_experiment_summary.csv"), row.names = FALSE)

cat("Saved pilot dimension experiment to", out_dir, "\n")
cat("\nDimension grid summary:\n")
print(grid, row.names = FALSE)
cat("\nTop adjusted RNA combined ranks:\n")
print(head(rna_rank[, c("feature", "p_event", "p_censor", "combined_score", "combined_rank")], 10), row.names = FALSE)
cat("\nTop adjusted CNA combined ranks:\n")
print(head(cna_rank[, c("feature", "p_event_overall", "p_censor_overall", "combined_score", "combined_rank")], 10), row.names = FALSE)
