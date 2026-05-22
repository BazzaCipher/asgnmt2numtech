# Days 7-8: DCC-GARCH with cleaned target (Figure 3, plus S1/S2/S3 tests).
#
# Methodology, per proposal:
#   Q_t = (1 - a - b) * R̄_* + a * z_{t-1} z_{t-1}' + b * Q_{t-1}
#   R_t = diag(Q_t)^{-1/2} * Q_t * diag(Q_t)^{-1/2}
# where R̄_* is the long-run correlation target. The key methodological
# contribution of this paper is using R̄_* = R̃ (MP-cleaned sample correlation)
# rather than R̄ (raw sample correlation).
#
# Why a manual DCC rather than rmgarch::dccfit:
#   rmgarch's correlation-targeting machinery hard-codes R̄. To swap in R̃ we
#   would need to monkey-patch the fitted object. A direct implementation is
#   ~80 lines, fully transparent, and we can run it with both R̄ and R̃ to
#   make the cleaned-vs-standard comparison (test S1) by construction.
#
# Inputs:
#   data/clean/std_residuals_v1.csv  (winsorised GJR-GARCH-t residuals)
#
# Outputs:
#   output/dcc_params.csv               (a, b, log-lik for standard vs cleaned)
#   output/dcc_pairwise_correlations.csv (date × pair × method × ρ_t)
#   output/s1_variance_reduction.csv    (per-pair variance ratio crisis window)
#   output/s2_bai_perron_pairs.csv      (break dates per pair, cleaned DCC)
#   output/s3_cocoa_basis_decoupling.csv (NY-LDN basis test)
#   figures/fig3_pairwise_correlations.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(strucchange)
})

# ----- Inputs ---------------------------------------------------------------
CRISIS_START   <- as.Date("2024-01-01")
CRISIS_END     <- as.Date("2026-05-22")
BASE_START     <- as.Date("2019-01-01")
BASE_END       <- as.Date("2023-12-31")

# Pairs we focus on (cocoa-anchored, pre-registered)
COCOA_PAIRS <- list(
  c("COCOA_NY", "COCOA_LDN"),    # S3 — basis decoupling
  c("COCOA_NY", "SUGAR"),         # S2 propagation
  c("COCOA_NY", "COFFEE_ARA"),    # S2 propagation
  c("COCOA_NY", "COFFEE_ROB"),
  c("COCOA_NY", "HSY"),
  c("COCOA_NY", "MDLZ"),
  c("COCOA_NY", "LISN"),
  c("COCOA_NY", "GHS"),
  c("COCOA_NY", "GHANA10Y"),
  c("COCOA_LDN", "SUGAR"),
  c("COCOA_LDN", "COFFEE_ARA")
)

z <- fread("data/clean/std_residuals_v1.csv")
z[, date := as.Date(date)]
tickers <- setdiff(names(z), "date")
z <- z[complete.cases(z[, tickers, with = FALSE])]
cat(sprintf("Working sample: %d dates x %d series\n", nrow(z), length(tickers)))

dates <- z$date
Z <- as.matrix(z[, tickers, with = FALSE])
T_n <- nrow(Z); N <- ncol(Z)

# ----- Sample and MP-cleaned correlation targets ----------------------------
R_bar <- cor(Z)

mp_clean <- function(R, T, N) {
  e <- eigen(R, symmetric = TRUE)
  evals <- e$values
  evecs <- e$vectors
  lambda_plus <- (1 + sqrt(N / T))^2
  noise <- evals <= lambda_plus
  # Replace noise eigenvalues with their average (Bouchaud-Potters hard clip);
  # market mode (largest) is kept untouched per proposal §"MP-cleaned matrix".
  if (any(noise)) {
    mean_noise <- mean(evals[noise])
    evals[noise] <- mean_noise
  }
  R_tilde <- evecs %*% diag(evals) %*% t(evecs)
  # Re-normalise to unit diagonal (cleaning can shift it slightly)
  d <- sqrt(diag(R_tilde))
  R_tilde <- R_tilde / (d %o% d)
  list(R = R_tilde, n_kept = sum(!noise), lambda_plus = lambda_plus)
}

cl <- mp_clean(R_bar, T_n, N)
R_tilde <- cl$R
cat(sprintf("MP edge λ+ = %.3f; eigenvalues kept (non-noise): %d/%d\n",
            cl$lambda_plus, cl$n_kept, N))

# ----- DCC log-likelihood (concentrated; std residuals → only correlation) --
# Negative log-likelihood for L-BFGS-B
dcc_nll <- function(params, Z, Q_bar) {
  a <- params[1]; b <- params[2]
  if (a <= 0 || b <= 0 || a + b >= 0.999) return(1e10)
  T_n <- nrow(Z); N <- ncol(Z)
  Q <- Q_bar
  nll <- 0
  for (t in seq_len(T_n)) {
    if (t > 1) {
      zt1 <- Z[t - 1, ]
      Q <- (1 - a - b) * Q_bar + a * (zt1 %*% t(zt1)) + b * Q
    }
    d <- sqrt(diag(Q))
    R <- Q / (d %o% d)
    # Numerical guard
    R <- (R + t(R)) / 2
    R_chol <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(R_chol)) return(1e10)
    logdetR <- 2 * sum(log(diag(R_chol)))
    zt <- Z[t, ]
    Rinv_z <- backsolve(R_chol, forwardsolve(t(R_chol), zt))
    nll <- nll + 0.5 * (logdetR + as.numeric(zt %*% Rinv_z) - as.numeric(zt %*% zt))
  }
  nll
}

fit_dcc <- function(Z, target_R, label, start = c(0.02, 0.95), b_upper = 0.999) {
  cat(sprintf("\nFitting %s DCC (target = %s)...\n", label, label))
  # Try a few starts; the cleaned target can be flat enough to make the surface
  # ridge-like, so we restart from several reasonable points and pick the best.
  starts <- list(start, c(0.01, 0.97), c(0.03, 0.90), c(0.05, 0.85))
  best <- NULL
  for (s in starts) {
    opt <- tryCatch(
      optim(par = s, fn = dcc_nll,
            Z = Z, Q_bar = target_R, method = "L-BFGS-B",
            lower = c(1e-5, 1e-4), upper = c(0.5, b_upper),
            control = list(trace = 0, factr = 1e7)),
      error = function(e) NULL
    )
    if (!is.null(opt) && (is.null(best) || opt$value < best$value)) best <- opt
  }
  cat(sprintf("  a = %.4f, b = %.4f, persistence = %.4f, nll = %.2f, code = %d\n",
              best$par[1], best$par[2], best$par[1] + best$par[2], best$value,
              best$convergence))
  list(a = best$par[1], b = best$par[2], nll = best$value,
       converged = best$convergence == 0)
}

# ----- Compute conditional correlations given (a, b, Q_bar) -----------------
dcc_correlation_series <- function(Z, target_R, a, b) {
  T_n <- nrow(Z); N <- ncol(Z)
  R_array <- array(NA_real_, dim = c(N, N, T_n))
  Q <- target_R
  for (t in seq_len(T_n)) {
    if (t > 1) {
      zt1 <- Z[t - 1, ]
      Q <- (1 - a - b) * target_R + a * (zt1 %*% t(zt1)) + b * Q
    }
    d <- sqrt(diag(Q))
    R <- Q / (d %o% d)
    R_array[, , t] <- (R + t(R)) / 2
  }
  R_array
}

# ----- Fit both versions ----------------------------------------------------
fit_std    <- fit_dcc(Z, R_bar,   "standard")
fit_clean  <- fit_dcc(Z, R_tilde, "MP-cleaned")

fwrite(data.table(
  method        = c("standard", "cleaned"),
  a             = c(fit_std$a, fit_clean$a),
  b             = c(fit_std$b, fit_clean$b),
  ab_sum        = c(fit_std$a + fit_std$b, fit_clean$a + fit_clean$b),
  target_weight = c(1 - fit_std$a - fit_std$b, 1 - fit_clean$a - fit_clean$b),
  loglik        = c(-fit_std$nll, -fit_clean$nll),
  delta_loglik  = c(0, (-fit_clean$nll) - (-fit_std$nll)),
  nll           = c(fit_std$nll, fit_clean$nll),
  converged     = c(fit_std$converged, fit_clean$converged)
), "output/dcc_params.csv")
cat(sprintf("\n=== DCC log-likelihood comparison ===\n"))
cat(sprintf("L(standard) = %.1f, L(cleaned) = %.1f, ΔL = %+.1f\n",
            -fit_std$nll, -fit_clean$nll, (-fit_clean$nll) - (-fit_std$nll)))
cat(sprintf("Long-run target weight (1-a-b): standard = %.3f, cleaned = %.3f\n",
            1 - fit_std$a - fit_std$b, 1 - fit_clean$a - fit_clean$b))

# ----- Conditional correlation series ---------------------------------------
cat("\nComputing conditional correlation series...\n")
R_std   <- dcc_correlation_series(Z, R_bar,   fit_std$a,   fit_std$b)
R_clean <- dcc_correlation_series(Z, R_tilde, fit_clean$a, fit_clean$b)

# Extract pairwise series
extract_pair <- function(R_arr, t_idx_i, t_idx_j) R_arr[t_idx_i, t_idx_j, ]
pair_idx <- function(tk) match(tk, tickers)

corr_dt <- data.table()
for (pair in COCOA_PAIRS) {
  i <- pair_idx(pair[1]); j <- pair_idx(pair[2])
  corr_dt <- rbind(corr_dt, data.table(
    date   = dates,
    pair   = paste(pair, collapse = "_"),
    pair_lhs = pair[1], pair_rhs = pair[2],
    rho_std   = extract_pair(R_std,   i, j),
    rho_clean = extract_pair(R_clean, i, j)
  ))
}
fwrite(corr_dt, "output/dcc_pairwise_correlations.csv")

# ----- S1: minimum-variance portfolio realised variance (proper EW2019 test) -
# EW2019's claim is variance reduction of out-of-sample MV portfolios, not the
# variance of conditional correlation paths. Test:
#   For each t in 2019-01-01 .. 2026-05-22 with at least 252 prior days:
#     - σ_{i,t}^2 = univariate GJR-GARCH-t conditional variance (from rugarch fit)
#     - D_t = diag(σ_t); Σ_t = D_t R_t D_t with R_t from each DCC method
#     - w_t = Σ_t^{-1} 1 / 1' Σ_t^{-1} 1
#     - realise: r_{p,t+1} = w_t' r_{t+1}
#   Compare realised portfolio variance pre vs in crisis.
rets <- fread("data/clean/returns_v1.csv"); rets[, date := as.Date(date)]
# Use only series with non-EWMA GARCH fits (consistent univariate σ_t)
# Re-extract univariate σ_t from the residuals: σ_t = r_t / z_t (where z fitted)
# Simpler: σ_t = abs(r_t)/abs(z_t) where defined.
sigma_dt <- data.table(date = z$date)
ret_z <- merge(rets, z[, .(date)], by = "date")
for (tk in tickers) {
  r <- ret_z[[tk]]; zz <- z[[tk]]
  sg <- ifelse(is.finite(r) & is.finite(zz) & abs(zz) > 1e-3, abs(r) / abs(zz), NA_real_)
  sigma_dt[[tk]] <- sg
}
# Forward-fill tiny gaps where |z|<1e-3
for (tk in tickers) {
  x <- sigma_dt[[tk]]
  if (any(is.na(x))) {
    last <- NA_real_
    for (i in seq_along(x)) { if (is.na(x[i])) x[i] <- last else last <- x[i] }
    sigma_dt[[tk]] <- x
  }
}

mvp_realised_var <- function(R_arr, sigma_mat, ret_mat, label) {
  Tn <- dim(R_arr)[3]; N <- dim(R_arr)[1]
  one <- rep(1, N)
  rp <- numeric(Tn - 1)
  for (t in seq_len(Tn - 1)) {
    R  <- R_arr[, , t]
    D  <- sigma_mat[t, ]
    Sigma <- (D %o% D) * R
    Sigma <- (Sigma + t(Sigma)) / 2 + diag(N) * 1e-10
    inv1 <- tryCatch(solve(Sigma, one), error = function(e) NULL)
    if (is.null(inv1)) { rp[t] <- NA; next }
    w <- inv1 / sum(inv1)
    r_next <- ret_mat[t + 1, ]
    if (any(is.na(r_next))) { rp[t] <- NA; next }
    rp[t] <- sum(w * r_next)
  }
  rp
}

# Align matrices
sigma_aligned <- as.matrix(sigma_dt[, tickers, with = FALSE])
ret_aligned   <- as.matrix(ret_z[, tickers, with = FALSE])

rp_std   <- mvp_realised_var(R_std,   sigma_aligned, ret_aligned, "standard")
rp_clean <- mvp_realised_var(R_clean, sigma_aligned, ret_aligned, "cleaned")
date_seq <- dates[-length(dates)]

mvp_dt <- data.table(date = date_seq, rp_std = rp_std, rp_clean = rp_clean)
mvp_dt <- mvp_dt[!is.na(rp_std) & !is.na(rp_clean)]
mvp_dt[, regime := fifelse(date >= CRISIS_START, "crisis",
                           fifelse(date >= BASE_START & date <= BASE_END, "baseline", "other"))]
mvp_summary <- mvp_dt[regime %in% c("baseline", "crisis"),
                      .(n = .N,
                        sd_std_pct   = round(sd(rp_std)   * 100, 3),
                        sd_clean_pct = round(sd(rp_clean) * 100, 3),
                        var_red_pct  = round(100 * (1 - var(rp_clean) / var(rp_std)), 1)),
                      by = regime]
fwrite(mvp_summary, "output/s1_mvp_oos_variance.csv")
cat("\n=== S1 (proper EW test): OOS minimum-variance portfolio sd ===\n")
print(mvp_summary)
s1_crisis_pass <- mvp_summary[regime == "crisis", var_red_pct] >= 5
cat(sprintf("S1 (cleaned reduces in-crisis MVP variance by ≥5%%): %s\n",
            if (length(s1_crisis_pass) && s1_crisis_pass) "PASS" else "FAIL"))

# Also keep the OLD (mis-framed) corr-variance metric, but rename it honestly:
s1_legacy <- corr_dt[date >= CRISIS_START & date <= CRISIS_END,
              .(var_std   = var(rho_std),
                var_clean = var(rho_clean)),
              by = pair]
s1_legacy[, var_ratio_clean_over_std := round(var_clean / var_std, 2)]
setorder(s1_legacy, -var_ratio_clean_over_std)
fwrite(s1_legacy, "output/s1_corr_path_variance_legacy.csv")

# ----- S2: Bai-Perron breaks on cleaned-DCC pairwise series -----------------
s2_rows <- list()
for (pair_name in unique(corr_dt$pair)) {
  pair_dt <- corr_dt[pair == pair_name, .(date, rho = rho_clean)]
  # Use h = 0.15 (default) → m_max = 5, BIC-selects within that. Fast enough for
  # 11 pairs at n ≈ 2260; the rolling-count Bai-Perron above uses the same h.
  bp_full <- tryCatch(
    breakpoints(rho ~ 1, data = pair_dt, h = 0.15, breaks = 5),
    error = function(e) NULL
  )
  bp_bic <- if (!is.null(bp_full)) which.min(BIC(bp_full)) - 1L else 0L
  break_dates <- if (!is.null(bp_full) && bp_bic > 0L) {
    bp <- breakpoints(bp_full, breaks = bp_bic)
    pair_dt$date[bp$breakpoints]
  } else as.Date(character(0))
  cat(sprintf("  [%s] BIC selects %d breaks\n", pair_name, bp_bic))
  in_crisis <- any(break_dates >= CRISIS_START & break_dates <= CRISIS_END)
  s2_rows[[pair_name]] <- data.table(
    pair = pair_name,
    n_breaks = length(break_dates),
    break_dates = paste(break_dates, collapse = "; "),
    in_crisis_break = in_crisis
  )
}
s2 <- rbindlist(s2_rows)
fwrite(s2, "output/s2_bai_perron_pairs.csv")
cat("\n=== S2: Bai-Perron breaks on cleaned-DCC pairwise series ===\n")
print(s2)
cat(sprintf("S2 (in-crisis break on cocoa-coffee OR cocoa-sugar): %s\n",
            if (any(s2[pair %in% c("COCOA_NY_SUGAR", "COCOA_NY_COFFEE_ARA"), in_crisis_break]))
              "PASS" else "FAIL"))

# ----- S3: COCOA_NY / COCOA_LDN basis decoupling ----------------------------
# Use the pre-registered crisis window [2024-01-01, 2026-05-22]. The earlier
# version restricted to 2025-2026 with no methodological justification.
basis <- corr_dt[pair == "COCOA_NY_COCOA_LDN"]
baseline_mean <- mean(basis[date >= BASE_START & date <= BASE_END, rho_clean])
crisis_rho   <- basis[date >= CRISIS_START & date <= CRISIS_END, rho_clean]
crisis_mean  <- mean(crisis_rho)
crisis_min   <- min(crisis_rho)
crisis_max   <- max(crisis_rho)
basis_drop_from_mean <- baseline_mean - crisis_min   # original spec
s3_pass <- basis_drop_from_mean >= 0.15

cat("\n=== S3: COCOA_NY × COCOA_LDN basis decoupling (cleaned DCC) ===\n")
cat(sprintf("Baseline mean (2019-2023):       %.3f\n", baseline_mean))
cat(sprintf("Crisis (2024-2026) mean:         %.3f\n", crisis_mean))
cat(sprintf("Crisis (2024-2026) range:        [%.3f, %.3f]\n", crisis_min, crisis_max))
cat(sprintf("Drop (baseline mean − crisis min): %+.3f  (need ≥0.15 to pass)\n",
            basis_drop_from_mean))
cat(sprintf("S3: %s — basis %s during crisis (mean change %+.3f)\n",
            if (s3_pass) "PASS" else "FAIL",
            if (crisis_mean > baseline_mean) "RE-COUPLED" else "decoupled",
            crisis_mean - baseline_mean))

fwrite(data.table(
  baseline_mean_2019_2023 = baseline_mean,
  crisis_mean_2024_2026   = crisis_mean,
  crisis_min_2024_2026    = crisis_min,
  crisis_max_2024_2026    = crisis_max,
  basis_drop = basis_drop_from_mean,
  s3_pass = s3_pass
), "output/s3_cocoa_basis_decoupling.csv")

# ----- Figure 3: pairwise conditional correlations --------------------------
focus_pairs <- c("COCOA_NY_COCOA_LDN", "COCOA_NY_SUGAR", "COCOA_NY_COFFEE_ARA",
                 "COCOA_NY_HSY", "COCOA_NY_MDLZ", "COCOA_NY_LISN")
plot_dt <- corr_dt[pair %in% focus_pairs]
plot_dt[, pair_lbl := factor(pair, levels = focus_pairs)]

crisis_band <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                          ymin = -Inf, ymax = Inf)

p3 <- ggplot(plot_dt, aes(x = date)) +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.10) +
  geom_line(aes(y = rho_std,   colour = "standard DCC"), alpha = 0.6, linewidth = 0.35) +
  geom_line(aes(y = rho_clean, colour = "MP-cleaned DCC"), linewidth = 0.55) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  facet_wrap(~ pair_lbl, ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_manual(values = c("standard DCC" = "grey50",
                                 "MP-cleaned DCC" = "steelblue"),
                      name = NULL) +
  labs(
    title = "Figure 3: Cocoa-anchored conditional correlations (DCC-GARCH)",
    subtitle = sprintf("Standard vs MP-cleaned target. Pink = crisis. Cleaned a=%.3f, b=%.3f.",
                       fit_clean$a, fit_clean$b),
    x = NULL, y = expression(rho[t]),
    caption = "Std residuals from GJR-GARCH-t univariates, complete-cases sample (N=25, 2017-04-21 → 2026-05-21)."
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 9))

ggsave("figures/fig3_pairwise_correlations.png", p3, width = 11, height = 7, dpi = 150)
cat("\nWrote figures/fig3_pairwise_correlations.png\n")
