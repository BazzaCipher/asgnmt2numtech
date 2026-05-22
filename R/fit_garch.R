# Days 3-4: Univariate GJR-GARCH(1,1) with Student-t innovations per asset.
#
# Inputs:
#   data/clean/returns_v1.csv       (frozen returns panel)
#
# Outputs:
#   data/clean/std_residuals_v1.csv (wide z_{i,t}, same shape as returns)
#   output/garch_fits.csv           (per-asset fitted params + diagnostics)
#
# Fallback chain (per scope_and_target.md):
#   Tier 1: GJR-GARCH(1,1) Student-t
#   Tier 2: GJR-GARCH(1,1) Gaussian
#   Tier 3: EGARCH(1,1) Student-t
#   Tier 4: EWMA standardisation (lambda = 0.94)
#
# Each tier is tried only if the previous failed convergence or threw an error.

suppressPackageStartupMessages({
  library(data.table)
  library(rugarch)
})

set.seed(42)

RETS_IN   <- "data/clean/returns_v1.csv"
RESID_OUT <- "data/clean/std_residuals_v1.csv"
FITS_OUT  <- "output/garch_fits.csv"

rets <- fread(RETS_IN)
rets[, date := as.Date(date)]
tickers <- setdiff(names(rets), "date")
cat(sprintf("Returns panel: %d dates x %d series\n", nrow(rets), length(tickers)))

# -----  Model specifications  ------------------------------------------------
spec_gjr_t <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
spec_gjr_n <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)
spec_egarch_t <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

# -----  Helpers  -------------------------------------------------------------
try_fit <- function(spec, x) {
  fit <- tryCatch(
    ugarchfit(spec, x, solver = "hybrid"),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) return(NULL)
  if (convergence(fit) != 0) return(NULL)
  # Sanity check: a "good" fit produces std residuals with sd ~ 1 and no
  # absurd outliers. A degenerate Student-t fit (df very close to 2) can pass
  # convergence but produce z values in the millions.
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

# -----  Fit per ticker  ------------------------------------------------------
fit_one <- function(ticker, x_raw) {
  valid_idx <- which(!is.na(x_raw))
  x <- x_raw[valid_idx]
  if (length(x) < 100) {
    return(list(model = "skipped", z = rep(NA_real_, length(x_raw)), info = NULL))
  }

  # Tier 1: GJR-GARCH-t
  fit <- try_fit(spec_gjr_t, x); model_name <- "gjr_t"
  # Tier 2: GJR-GARCH-N
  if (is.null(fit)) { fit <- try_fit(spec_gjr_n, x); model_name <- "gjr_n" }
  # Tier 3: EGARCH-t
  if (is.null(fit)) { fit <- try_fit(spec_egarch_t, x); model_name <- "egarch_t" }
  # Tier 4: EWMA
  if (is.null(fit)) {
    z_clean <- ewma_residuals(x)
    z_full <- rep(NA_real_, length(x_raw)); z_full[valid_idx] <- z_clean
    info <- list(model = "ewma", n = length(x), omega = NA, alpha1 = NA, beta1 = NA,
                 gamma1 = NA, shape = NA, persistence = NA, loglik = NA,
                 LB10_z = NA, LB10_z2 = NA)
    return(list(model = "ewma", z = z_full, info = info))
  }

  z_clean <- as.numeric(residuals(fit, standardize = TRUE))
  z_full  <- rep(NA_real_, length(x_raw)); z_full[valid_idx] <- z_clean

  lb_z  <- tryCatch(Box.test(z_clean,    lag = 10, type = "Ljung-Box")$p.value, error = function(e) NA)
  lb_z2 <- tryCatch(Box.test(z_clean^2,  lag = 10, type = "Ljung-Box")$p.value, error = function(e) NA)

  info <- list(
    model = model_name, n = length(x),
    omega  = extract_coef(fit, "omega"),
    alpha1 = extract_coef(fit, "alpha1"),
    beta1  = extract_coef(fit, "beta1"),
    gamma1 = extract_coef(fit, "gamma1"),
    shape  = extract_coef(fit, "shape"),
    persistence = tryCatch(persistence(fit), error = function(e) NA),
    loglik = tryCatch(likelihood(fit),       error = function(e) NA),
    LB10_z = lb_z, LB10_z2 = lb_z2
  )
  list(model = model_name, z = z_full, info = info)
}

# -----  Run  -----------------------------------------------------------------
std_resid <- data.table(date = rets$date)
diag_list <- list()

for (ticker in tickers) {
  cat(sprintf("  Fitting %s ... ", ticker))
  res <- fit_one(ticker, rets[[ticker]])
  std_resid[[ticker]] <- res$z
  diag_list[[ticker]] <- c(list(ticker = ticker), res$info)
  cat(sprintf("[%s]\n", res$model))
}

diag_dt <- rbindlist(diag_list, fill = TRUE)
setcolorder(diag_dt, c("ticker", "model", "n", "omega", "alpha1", "beta1", "gamma1",
                       "shape", "persistence", "loglik", "LB10_z", "LB10_z2"))

# ----- Winsorise std residuals at ±5 ----------------------------------------
# Standard practice for RMT financial analysis (Bouchaud-Potters): a single
# extreme z can dominate a 252-day correlation window and distort the spectrum.
# We cap at ±5 (a 1-in-3.5M event under Gaussian; clearly a model-failure tail
# under any reasonable distribution). Affected counts logged in diagnostics.
WINSORIZE_AT <- 5
n_winsorized <- integer(length(tickers)); names(n_winsorized) <- tickers
for (ticker in tickers) {
  x <- std_resid[[ticker]]
  hi <- !is.na(x) & x >  WINSORIZE_AT
  lo <- !is.na(x) & x < -WINSORIZE_AT
  n_winsorized[ticker] <- sum(hi) + sum(lo)
  x[hi] <-  WINSORIZE_AT
  x[lo] <- -WINSORIZE_AT
  std_resid[[ticker]] <- x
}
diag_dt[, n_winsorized := n_winsorized[ticker]]

dir.create("output", showWarnings = FALSE)
fwrite(std_resid, RESID_OUT)
fwrite(diag_dt,   FITS_OUT)

cat(sprintf("\nWrote %s\n", RESID_OUT))
cat(sprintf("Wrote %s\n", FITS_OUT))

# -----  Console summary  -----------------------------------------------------
cat("\n=== Fit summary ===\n")
print(diag_dt[, .(ticker, model, n, alpha1 = round(alpha1, 3),
                  beta1 = round(beta1, 3), gamma1 = round(gamma1, 3),
                  shape = round(shape, 2), persist = round(persistence, 3),
                  LB_z = round(LB10_z, 3), LB_z2 = round(LB10_z2, 3))])

cat("\n=== Model tier counts ===\n")
print(diag_dt[, .N, by = model])

cat("\n=== Diagnostic flags ===\n")
bad_z  <- diag_dt[!is.na(LB10_z)  & LB10_z  < 0.05, ticker]
bad_z2 <- diag_dt[!is.na(LB10_z2) & LB10_z2 < 0.05, ticker]
cat(sprintf("Series with residual autocorrelation (Ljung-Box on z, p<0.05):  %s\n",
            if (length(bad_z))  paste(bad_z,  collapse = ", ") else "none"))
cat(sprintf("Series with residual ARCH effects (Ljung-Box on z^2, p<0.05): %s\n",
            if (length(bad_z2)) paste(bad_z2, collapse = ", ") else "none"))
