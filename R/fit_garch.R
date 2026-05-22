# Days 3-4: Univariate GJR-GARCH(1,1) per asset.
#
# Tier 1a: GJR-GARCH(1,1) skewed-Student-t (sstd)   — primary (Hansen 1994)
# Tier 1b: GJR-GARCH(1,1) symmetric Student-t (std) — fallback if sstd fails
# Tier 2:  GJR-GARCH(1,1) Gaussian (norm)
# Tier 3:  EGARCH(1,1) sstd
# Tier 4:  EWMA standardisation (lambda = 0.94)
#
# Mean equation: AR(0) by default. OJ and VND ALWAYS use AR(1); any other series
# whose Tier-1 fit shows Ljung-Box p < 0.01 on standardised residuals is
# automatically retried with AR(1) on the same distribution tier. Reason:
# residual autocorrelation in z propagates to off-diagonal correlation entries
# and can manufacture spurious supra-MP eigenvalues.
#
# Outputs:
#   data/clean/std_residuals_v1.csv  — winsorised z (wide, primary downstream input)
#   data/clean/std_residuals_norm.csv — Gaussian-innovations z (R4 robustness)
#   output/garch_fits.csv            — per-asset diagnostics (legacy filename)
#   output/garch_fits_sstd.csv       — same content under friend's filename
#   output/return_skewness.csv       — empirical skew/excess-kurt motivating sstd
#   output/near_igarch.csv           — series with persistence >= 0.99

suppressPackageStartupMessages({
  library(data.table)
  library(rugarch)
  library(parallel)
})

set.seed(42)

NCORES <- max(1L, min(8L, detectCores() - 1L))

RETS_IN    <- "data/clean/returns_v1.csv"
RESID_OUT  <- "data/clean/std_residuals_v1.csv"
RESID_NORM <- "data/clean/std_residuals_norm.csv"
FITS_OUT   <- "output/garch_fits.csv"
FITS_SSTD  <- "output/garch_fits_sstd.csv"
SKEW_OUT   <- "output/return_skewness.csv"
NEAR_IG    <- "output/near_igarch.csv"

ALWAYS_AR1 <- c("OJ", "VND")
LBZ_RETRY  <- 0.01     # if LB-z p < this, retry AR(0) fit with AR(1)
NEAR_IGARCH_THRESHOLD <- 0.99

rets <- fread(RETS_IN)
rets[, date := as.Date(date)]
tickers <- setdiff(names(rets), "date")
cat(sprintf("Returns panel: %d dates x %d series\n", nrow(rets), length(tickers)))

# -----  Empirical skewness / kurtosis (motivation table) --------------------
skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x); m <- mean(x); s <- sd(x)
  if (s == 0 || n < 3) return(NA_real_)
  mean(((x - m) / s)^3)
}
ex_kurt <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x); m <- mean(x); s <- sd(x)
  if (s == 0 || n < 4) return(NA_real_)
  mean(((x - m) / s)^4) - 3
}
skew_dt <- data.table(
  ticker      = tickers,
  n_obs       = sapply(tickers, function(t) sum(is.finite(rets[[t]]))),
  skew        = sapply(tickers, function(t) round(skewness(rets[[t]]), 3)),
  excess_kurt = sapply(tickers, function(t) round(ex_kurt(rets[[t]]),  3))
)
fwrite(skew_dt, SKEW_OUT)

# -----  Model specifications  ------------------------------------------------
make_spec <- function(model, distr, ar) {
  variance.model <- if (model == "eGARCH") {
    list(model = "eGARCH", garchOrder = c(1, 1))
  } else {
    list(model = "gjrGARCH", garchOrder = c(1, 1))
  }
  ugarchspec(
    variance.model = variance.model,
    mean.model     = list(armaOrder = c(ar, 0), include.mean = TRUE),
    distribution.model = distr
  )
}

# -----  Helpers  -------------------------------------------------------------
try_fit <- function(spec, x) {
  fit <- tryCatch(
    ugarchfit(spec, x, solver = "hybrid"),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) return(NULL)
  if (convergence(fit) != 0) return(NULL)
  z <- as.numeric(residuals(fit, standardize = TRUE))
  if (!all(is.finite(z))) return(NULL)
  if (sd(z) > 2)          return(NULL)
  if (max(abs(z)) > 30)   return(NULL)
  fit
}

ewma_residuals <- function(x, lambda = 0.94) {
  n <- length(x)
  sigma2 <- numeric(n)
  sigma2[1] <- var(x, na.rm = TRUE)
  for (t in 2:n) sigma2[t] <- lambda * sigma2[t - 1] + (1 - lambda) * x[t - 1]^2
  x / sqrt(sigma2)
}

extract_coef <- function(fit, name) {
  cf <- coef(fit)
  if (name %in% names(cf)) unname(cf[name]) else NA_real_
}

# LR test sstd vs std: 2 * (L_sstd - L_std) ~ chi^2(1) under H0: skew = 1
lr_sstd_vs_std <- function(fit_sstd, fit_std) {
  if (is.null(fit_sstd) || is.null(fit_std)) return(c(NA_real_, NA_real_))
  L1 <- tryCatch(likelihood(fit_sstd), error = function(e) NA_real_)
  L0 <- tryCatch(likelihood(fit_std),  error = function(e) NA_real_)
  if (!is.finite(L1) || !is.finite(L0)) return(c(NA_real_, NA_real_))
  stat <- 2 * (L1 - L0)
  if (stat < 0) stat <- 0  # numerical
  c(stat, pchisq(stat, df = 1, lower.tail = FALSE))
}

# -----  Fit per ticker (returns metadata for primary + norm variants) -------
fit_one <- function(ticker, x_raw) {
  valid_idx <- which(!is.na(x_raw))
  x <- x_raw[valid_idx]
  n_obs <- length(x_raw)
  empty <- list(
    model = "skipped", ar = NA_integer_, z = rep(NA_real_, n_obs),
    z_norm = rep(NA_real_, n_obs), info = NULL
  )
  if (length(x) < 100) return(empty)

  ar <- if (ticker %in% ALWAYS_AR1) 1L else 0L

  fit_sstd <- try_fit(make_spec("gjrGARCH", "sstd", ar), x)
  fit_std  <- try_fit(make_spec("gjrGARCH", "std",  ar), x)
  fit_norm <- try_fit(make_spec("gjrGARCH", "norm", ar), x)

  # LR test (only meaningful if both converge)
  lr <- lr_sstd_vs_std(fit_sstd, fit_std)

  fit <- fit_sstd; model_name <- "gjr_sstd"
  if (is.null(fit)) { fit <- fit_std;  model_name <- "gjr_t" }
  if (is.null(fit)) { fit <- fit_norm; model_name <- "gjr_n" }
  if (is.null(fit)) {
    fit <- try_fit(make_spec("eGARCH", "sstd", ar), x); model_name <- "egarch_sstd"
  }
  if (is.null(fit)) {
    z_clean <- ewma_residuals(x)
    z_full <- rep(NA_real_, n_obs); z_full[valid_idx] <- z_clean
    z_norm_full <- if (!is.null(fit_norm)) {
      zn <- as.numeric(residuals(fit_norm, standardize = TRUE))
      out <- rep(NA_real_, n_obs); out[valid_idx] <- zn; out
    } else rep(NA_real_, n_obs)
    info <- list(model = "ewma", ar = ar, n = length(x),
                 omega = NA, alpha1 = NA, beta1 = NA, gamma1 = NA,
                 skew = NA, shape = NA, persistence = NA, loglik = NA,
                 LB10_z = NA, LB10_z2 = NA,
                 lr_stat = lr[1], lr_pval = lr[2],
                 ar1_retry = FALSE)
    return(list(model = "ewma", ar = ar, z = z_full, z_norm = z_norm_full,
                info = info))
  }

  z_clean <- as.numeric(residuals(fit, standardize = TRUE))
  lb_z  <- tryCatch(Box.test(z_clean, lag = 10, type = "Ljung-Box")$p.value,
                    error = function(e) NA)

  # AR(1) auto-retry rule: if AR(0) Tier-1 fit shows residual autocorrelation,
  # retry the same distribution tier with AR(1).
  ar1_retry <- FALSE
  if (ar == 0L && !is.na(lb_z) && lb_z < LBZ_RETRY &&
      model_name %in% c("gjr_sstd", "gjr_t", "gjr_n", "egarch_sstd")) {
    distr <- switch(model_name,
                    gjr_sstd = "sstd", gjr_t = "std",
                    gjr_n = "norm", egarch_sstd = "sstd")
    model <- if (model_name == "egarch_sstd") "eGARCH" else "gjrGARCH"
    fit_retry <- try_fit(make_spec(model, distr, 1L), x)
    if (!is.null(fit_retry)) {
      z_retry <- as.numeric(residuals(fit_retry, standardize = TRUE))
      lb_retry <- tryCatch(Box.test(z_retry, lag = 10, type = "Ljung-Box")$p.value,
                           error = function(e) NA)
      if (is.finite(lb_retry) && lb_retry > lb_z) {
        fit <- fit_retry; z_clean <- z_retry; lb_z <- lb_retry
        ar <- 1L; ar1_retry <- TRUE
      }
    }
  }

  z_full <- rep(NA_real_, n_obs); z_full[valid_idx] <- z_clean
  z_norm_full <- if (!is.null(fit_norm)) {
    zn <- as.numeric(residuals(fit_norm, standardize = TRUE))
    out <- rep(NA_real_, n_obs); out[valid_idx] <- zn; out
  } else rep(NA_real_, n_obs)

  lb_z2 <- tryCatch(Box.test(z_clean^2, lag = 10, type = "Ljung-Box")$p.value,
                    error = function(e) NA)

  info <- list(
    model = model_name, ar = ar, n = length(x),
    omega  = extract_coef(fit, "omega"),
    alpha1 = extract_coef(fit, "alpha1"),
    beta1  = extract_coef(fit, "beta1"),
    gamma1 = extract_coef(fit, "gamma1"),
    skew   = extract_coef(fit, "skew"),
    shape  = extract_coef(fit, "shape"),
    persistence = tryCatch(persistence(fit), error = function(e) NA),
    loglik = tryCatch(likelihood(fit),       error = function(e) NA),
    LB10_z = lb_z, LB10_z2 = lb_z2,
    lr_stat = lr[1], lr_pval = lr[2],
    ar1_retry = ar1_retry
  )
  list(model = model_name, ar = ar, z = z_full, z_norm = z_norm_full,
       info = info)
}

# -----  Run (parallel over tickers) -----------------------------------------
cat(sprintf("Fitting %d series in parallel on %d cores...\n", length(tickers), NCORES))
results <- mclapply(tickers, function(ticker) {
  res <- fit_one(ticker, rets[[ticker]])
  list(ticker = ticker, model = res$model, ar = res$ar,
       z = res$z, z_norm = res$z_norm, info = res$info,
       ar1_retry = isTRUE(res$info$ar1_retry))
}, mc.cores = NCORES, mc.preschedule = FALSE)

# mclapply on fork-based parallel can return try-error objects on child crash;
# fall back to serial for any failures.
failed <- vapply(results, inherits, logical(1L), what = "try-error")
if (any(failed)) {
  cat(sprintf("WARN: %d series failed in parallel; retrying serially.\n", sum(failed)))
  for (i in which(failed)) {
    res <- fit_one(tickers[i], rets[[tickers[i]]])
    results[[i]] <- list(ticker = tickers[i], model = res$model, ar = res$ar,
                         z = res$z, z_norm = res$z_norm, info = res$info,
                         ar1_retry = isTRUE(res$info$ar1_retry))
  }
}

std_resid <- data.table(date = rets$date)
std_norm  <- data.table(date = rets$date)
diag_list <- list()
for (r in results) {
  std_resid[[r$ticker]] <- r$z
  std_norm[[r$ticker]]  <- r$z_norm
  diag_list[[r$ticker]] <- c(list(ticker = r$ticker), r$info)
  cat(sprintf("  %s [%s, AR(%d)%s]\n", r$ticker, r$model, r$ar,
              if (r$ar1_retry) " RETRY" else ""))
}

diag_dt <- rbindlist(diag_list, fill = TRUE)
setcolorder(diag_dt, c("ticker", "model", "ar", "n", "omega", "alpha1", "beta1",
                       "gamma1", "skew", "shape", "persistence", "loglik",
                       "LB10_z", "LB10_z2", "lr_stat", "lr_pval", "ar1_retry"))

# ----- Winsorise std residuals at ±5 ----------------------------------------
WINSORIZE_AT <- 5
winsorise_inplace <- function(dt) {
  n_w <- integer(length(tickers)); names(n_w) <- tickers
  for (ticker in tickers) {
    x <- dt[[ticker]]
    hi <- !is.na(x) & x >  WINSORIZE_AT
    lo <- !is.na(x) & x < -WINSORIZE_AT
    n_w[ticker] <- sum(hi) + sum(lo)
    x[hi] <-  WINSORIZE_AT
    x[lo] <- -WINSORIZE_AT
    dt[[ticker]] <- x
  }
  n_w
}
n_winsorized <- winsorise_inplace(std_resid)
winsorise_inplace(std_norm)
diag_dt[, n_winsorized := n_winsorized[ticker]]

dir.create("output", showWarnings = FALSE)
fwrite(std_resid, RESID_OUT)
fwrite(std_norm,  RESID_NORM)
fwrite(diag_dt,   FITS_OUT)
fwrite(diag_dt,   FITS_SSTD)   # friend's filename = same data

# ----- Near-IGARCH diagnostic ------------------------------------------------
near_ig <- diag_dt[!is.na(persistence) & persistence >= NEAR_IGARCH_THRESHOLD,
                   .(ticker, model, persistence = round(persistence, 4),
                     alpha1 = round(alpha1, 3), beta1 = round(beta1, 3),
                     gamma1 = round(gamma1, 3))]
setorder(near_ig, -persistence)
fwrite(near_ig, NEAR_IG)

cat(sprintf("\nWrote %s\n", RESID_OUT))
cat(sprintf("Wrote %s\n", RESID_NORM))
cat(sprintf("Wrote %s and %s\n", FITS_OUT, FITS_SSTD))
cat(sprintf("Wrote %s\n", SKEW_OUT))
cat(sprintf("Wrote %s (%d series with persistence >= %.2f)\n",
            NEAR_IG, nrow(near_ig), NEAR_IGARCH_THRESHOLD))

# -----  Console summary  -----------------------------------------------------
cat("\n=== Fit summary ===\n")
print(diag_dt[, .(ticker, model, ar, alpha1 = round(alpha1, 3),
                  beta1 = round(beta1, 3), gamma1 = round(gamma1, 3),
                  skew = round(skew, 2), shape = round(shape, 2),
                  persist = round(persistence, 3),
                  LB_z = round(LB10_z, 3), LB_z2 = round(LB10_z2, 3),
                  retry = ar1_retry)])

cat("\n=== Model tier counts ===\n")
print(diag_dt[, .N, by = model])

cat("\n=== LR test (sstd vs std), p < 0.05 ===\n")
sig_skew <- diag_dt[!is.na(lr_pval) & lr_pval < 0.05,
                    .(ticker, skew = round(skew, 3),
                      lr_stat = round(lr_stat, 2), lr_pval = round(lr_pval, 4))]
if (nrow(sig_skew)) print(sig_skew) else cat("(none)\n")

cat("\n=== AR(1) retry / always-AR(1) ===\n")
ar1_series <- diag_dt[ar == 1L,
                      .(ticker, model, ar1_retry,
                        LB_z = round(LB10_z, 4))]
print(ar1_series)

cat("\n=== Near-IGARCH series (persistence >= 0.99) ===\n")
print(near_ig)

cat("\n=== Diagnostic flags ===\n")
bad_z  <- diag_dt[!is.na(LB10_z)  & LB10_z  < 0.05, ticker]
bad_z2 <- diag_dt[!is.na(LB10_z2) & LB10_z2 < 0.05, ticker]
cat(sprintf("Series with residual autocorrelation (LB on z, p<0.05): %s\n",
            if (length(bad_z))  paste(bad_z,  collapse = ", ") else "none"))
cat(sprintf("Series with residual ARCH effects (LB on z^2, p<0.05): %s\n",
            if (length(bad_z2)) paste(bad_z2, collapse = ", ") else "none"))
