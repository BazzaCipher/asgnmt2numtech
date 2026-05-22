# Strict end-to-end diagnostics. No claims are accepted without a check.
# Prints numbers at every stage so plausibility can be verified by hand.

suppressPackageStartupMessages({
  library(data.table)
})
# Light-weight moments
skewness <- function(x) { x <- x[is.finite(x)]; n <- length(x); m <- mean(x); s <- sd(x); sum((x-m)^3)/(n*s^3) }
kurtosis <- function(x) { x <- x[is.finite(x)]; n <- length(x); m <- mean(x); s <- sd(x); sum((x-m)^4)/(n*s^4) }

options(width = 140)
sep <- function(s) cat(sprintf("\n========== %s ==========\n", s))

# ============================================================================
# STAGE 1: RAW LEVELS PANEL
# ============================================================================
sep("STAGE 1: levels panel (data/clean/panel_v1.csv)")
panel <- fread("data/clean/panel_v1.csv"); panel[, date := as.Date(date)]
cat(sprintf("rows = %d (dates), cols = %d (incl date)\n", nrow(panel), ncol(panel)))
cat(sprintf("date range: %s to %s\n", min(panel$date), max(panel$date)))
tickers <- setdiff(names(panel), "date")

# Per-series level statistics with focus on cocoa
lev_stats <- rbindlist(lapply(tickers, function(t) {
  x <- panel[[t]]; nx <- sum(!is.na(x))
  if (nx == 0) return(NULL)
  list(ticker = t, n = nx,
       min = min(x, na.rm = TRUE), p50 = median(x, na.rm = TRUE),
       max = max(x, na.rm = TRUE), max_date = panel$date[which.max(x)],
       last = tail(na.omit(x), 1))
}))
print(lev_stats)

cat("\n--- COCOA_NY peak check (paper claims 12565 on 2024-12-18) ---\n")
peak_row <- panel[which.max(COCOA_NY), .(date, COCOA_NY)]
print(peak_row)
cat(sprintf("In window 2024-11-15..2024-12-31, max = %.2f on %s\n",
            panel[date >= "2024-11-15" & date <= "2024-12-31", max(COCOA_NY, na.rm = TRUE)],
            panel[date >= "2024-11-15" & date <= "2024-12-31"][which.max(COCOA_NY), date]))

# ============================================================================
# STAGE 2: RETURNS PANEL
# ============================================================================
sep("STAGE 2: returns panel (data/clean/returns_v1.csv)")
rets <- fread("data/clean/returns_v1.csv"); rets[, date := as.Date(date)]
cat(sprintf("rows = %d, cols = %d\n", nrow(rets), ncol(rets)))

ret_stats <- rbindlist(lapply(tickers, function(t) {
  x <- rets[[t]]; xc <- na.omit(x); nx <- length(xc)
  if (nx == 0) return(NULL)
  list(ticker = t, n = nx,
       mean_bp  = round(mean(xc) * 1e4, 2),     # bp
       sd_pct   = round(sd(xc) * 100, 3),       # %
       skew     = round(skewness(xc), 2),
       kurt     = round(kurtosis(xc), 1),
       min      = round(min(xc), 3),
       max      = round(max(xc), 3),
       n_inf    = sum(!is.finite(xc)),
       n_zero   = sum(xc == 0))
}))
print(ret_stats)

# WTI patch check
cat("\n--- WTI patch check (2020-04-20, 21 should be NA) ---\n")
print(rets[date %in% as.Date(c("2020-04-17","2020-04-20","2020-04-21","2020-04-22")),
           .(date, WTI)])

# ============================================================================
# STAGE 3: GARCH FIT DIAGNOSTICS
# ============================================================================
sep("STAGE 3: GARCH fits (output/garch_fits.csv)")
fits <- fread("output/garch_fits.csv")
fits[, persistence_recalc := round(alpha1 + beta1 + 0.5 * gamma1, 4)]
fits[, neg_leverage := gamma1 < 0]
fits[, near_igarch  := persistence > 0.995]
fits[, LB_z_fail    := !is.na(LB10_z) & LB10_z < 0.05]
fits[, LB_z2_fail   := !is.na(LB10_z2) & LB10_z2 < 0.05]

print(fits[, .(ticker, model, alpha1 = round(alpha1, 3), beta1 = round(beta1, 3),
               gamma1 = round(gamma1, 3), persistence = round(persistence, 4),
               LB10_z = round(LB10_z, 4), LB10_z2 = round(LB10_z2, 4),
               n_wins = n_winsorized,
               neg_lev = neg_leverage, near_iG = near_igarch,
               LB_z_fail, LB_z2_fail)])

cat("\n--- Series flagged for at least one defect ---\n")
flagged <- fits[neg_leverage | near_igarch | LB_z_fail | LB_z2_fail | model != "gjr_t"]
print(flagged[, .(ticker, model, gamma1 = round(gamma1, 3),
                  persistence = round(persistence, 3),
                  LB10_z = round(LB10_z, 4), LB10_z2 = round(LB10_z2, 4))])
cat(sprintf("\n%d / %d series flagged.\n", nrow(flagged), nrow(fits)))

# ============================================================================
# STAGE 4: STANDARDISED RESIDUALS
# ============================================================================
sep("STAGE 4: standardised residuals (data/clean/std_residuals_v1.csv)")
z <- fread("data/clean/std_residuals_v1.csv"); z[, date := as.Date(date)]
zstats <- rbindlist(lapply(tickers, function(t) {
  x <- na.omit(z[[t]])
  list(ticker = t, n = length(x),
       mean = round(mean(x), 3), sd = round(sd(x), 3),
       skew = round(skewness(x), 2), kurt = round(kurtosis(x), 1),
       absmax = round(max(abs(x)), 2),
       n_eq_5 = sum(abs(x) == 5))
}))
print(zstats)

cat("\n--- z stats that should hold: mean~=0, sd~=1, skew small. Flags: ---\n")
print(zstats[abs(mean) > 0.1 | abs(sd - 1) > 0.05 | abs(skew) > 0.5 | absmax >= 5])

# ============================================================================
# STAGE 5: FULL-SAMPLE MP SPECTRUM
# ============================================================================
sep("STAGE 5: full-sample MP spectrum")
z_cc <- z[complete.cases(z[, tickers, with = FALSE])]
cat(sprintf("complete-cases rows = %d (from %s to %s)\n",
            nrow(z_cc), min(z_cc$date), max(z_cc$date)))
N  <- length(tickers); Tf <- nrow(z_cc)
q  <- N / Tf
mp_hi <- (1 + sqrt(q))^2; mp_lo <- (1 - sqrt(q))^2
R_full <- cor(as.matrix(z_cc[, tickers, with = FALSE]))
ev_full <- sort(eigen(R_full, symmetric = TRUE, only.values = TRUE)$values,
                decreasing = TRUE)
cat(sprintf("N=%d  T=%d  Q=T/N=%.1f  MP edges: [%.4f, %.4f]\n",
            N, Tf, 1/q, mp_lo, mp_hi))
cat("Top 8 eigenvalues:\n"); print(round(ev_full[1:8], 3))
cat(sprintf("count above MP edge (incl market): %d\n", sum(ev_full > mp_hi)))
cat(sprintf("count above MP edge (excl market): %d\n", sum(ev_full[-1] > mp_hi)))
cat(sprintf("variance explained by top 1 mode: %.1f%%\n", 100 * ev_full[1] / N))
cat(sprintf("variance explained by top 6 modes: %.1f%%\n", 100 * sum(ev_full[1:6]) / N))

# ============================================================================
# STAGE 6: SUBSAMPLE SPECTRA (pre-crisis vs in-crisis)
# ============================================================================
sep("STAGE 6: subsample spectra")
COCOA_BLOC <- c("COCOA_NY","COCOA_LDN","HSY","MDLZ","LISN","GHS","GHANA10Y")
ALT_BLOC   <- c("COCOA_NY","COCOA_LDN","HSY","MDLZ","LISN","SJM")   # no Ghana, no GHS
NARROW_BLOC <- c("COCOA_NY","COCOA_LDN","HSY","MDLZ","LISN")

analyze <- function(zsub, label) {
  N <- length(tickers); Tn <- nrow(zsub)
  q <- N / Tn; edge <- (1 + sqrt(q))^2
  R <- cor(as.matrix(zsub[, tickers, with = FALSE]))
  e <- eigen(R, symmetric = TRUE)
  ev <- e$values
  vecs <- e$vectors
  ord <- order(ev, decreasing = TRUE)
  ev <- ev[ord]; vecs <- vecs[, ord]
  n_supra <- sum(ev > edge)
  cat(sprintf("\n--- %s | T=%d N=%d MP+ = %.3f ---\n", label, Tn, N, edge))
  cat(sprintf("supra-MP count: %d (incl market), %d (excl market)\n",
              n_supra, sum(ev[-1] > edge)))
  cat("Top 6 eigenvalues:\n"); print(round(ev[1:6], 3))

  bloc_share <- function(bloc, k) {
    idx <- match(bloc, tickers)
    sum(vecs[idx, k]^2)
  }
  cat("Cocoa-bloc squared-loading share, per rank:\n")
  shares <- data.table(
    rank = 1:6,
    eigenvalue = round(ev[1:6], 3),
    above_MP = ev[1:6] > edge,
    bloc7  = round(sapply(1:6, function(k) bloc_share(COCOA_BLOC, k)), 3),
    bloc6_noGhana = round(sapply(1:6, function(k) bloc_share(ALT_BLOC, k)), 3),
    bloc5_narrow  = round(sapply(1:6, function(k) bloc_share(NARROW_BLOC, k)), 3)
  )
  print(shares)
  list(ev = ev, vecs = vecs, edge = edge, n = Tn)
}

pre <- analyze(z_cc[date >= "2019-01-01" & date <= "2023-12-31"], "PRE 2019-2023")
crs <- analyze(z_cc[date >= "2024-01-01" & date <= "2026-05-22"], "CRISIS 2024-2026")

cat("\n--- T3 (paper version: rank-1 bloc share, ≥0.50 AND doubles) ---\n")
pre_b1 <- sum(pre$vecs[match(COCOA_BLOC, tickers), 1]^2)
crs_b1 <- sum(crs$vecs[match(COCOA_BLOC, tickers), 1]^2)
cat(sprintf("rank-1 bloc share: pre=%.3f, crisis=%.3f, ratio=%.2fx, doubled=%s, >=0.50=%s\n",
            pre_b1, crs_b1, crs_b1/pre_b1, crs_b1 >= 2*pre_b1, crs_b1 >= 0.50))

cat("\n--- T3 (original spec: emerging eigenvalue, bloc share ≥0.55) ---\n")
pre_supra <- which(pre$ev > pre$edge)
crs_supra <- which(crs$ev > crs$edge)
cat(sprintf("pre supra ranks:    %s\n", paste(pre_supra, collapse = ", ")))
cat(sprintf("crisis supra ranks: %s\n", paste(crs_supra, collapse = ", ")))
emerging <- setdiff(crs_supra, pre_supra)
cat(sprintf("emerging (in crisis, not in pre): %s\n",
            if (length(emerging)) paste(emerging, collapse = ", ") else "(none)"))
if (length(emerging)) {
  k <- min(emerging)
  bs <- sum(crs$vecs[match(COCOA_BLOC, tickers), k]^2)
  cat(sprintf("Bloc share at emerging rank %d: %.3f (T3 threshold 0.55, pass=%s)\n",
              k, bs, bs >= 0.55))
}

# ============================================================================
# STAGE 7: ROLLING COUNT SANITY
# ============================================================================
sep("STAGE 7: rolling count series")
roll <- fread("output/rolling_count_v1.csv")
roll[, date := as.Date(date)]
cat(sprintf("rows = %d, n_assets values: %s\n", nrow(roll),
            paste(sort(unique(roll$n_assets)), collapse = ", ")))

cat("\nDistribution of count_above_mp_ex_market by n_assets:\n")
print(roll[, .(min = min(count_above_mp_ex_market, na.rm=TRUE),
               p25 = quantile(count_above_mp_ex_market, .25, na.rm=TRUE),
               p50 = median(count_above_mp_ex_market, na.rm=TRUE),
               p75 = quantile(count_above_mp_ex_market, .75, na.rm=TRUE),
               max = max(count_above_mp_ex_market, na.rm=TRUE),
               n = .N), by = n_assets])

cat("\nIn-crisis (W ending in [2024-01-01, 2026-05-22]) distribution:\n")
print(roll[date >= "2024-01-01" & date <= "2026-05-22",
           .(n=.N, mean=round(mean(count_above_mp_ex_market),2),
             median=median(count_above_mp_ex_market),
             min=min(count_above_mp_ex_market), max=max(count_above_mp_ex_market))])

cat("\nPre-crisis (W ending in [2019-01-01, 2023-12-31] & n_assets==25):\n")
print(roll[date >= "2019-01-01" & date <= "2023-12-31" & n_assets == 25,
           .(n=.N, mean=round(mean(count_above_mp_ex_market),2),
             median=median(count_above_mp_ex_market),
             min=min(count_above_mp_ex_market), max=max(count_above_mp_ex_market),
             share_eq_2=round(mean(count_above_mp_ex_market==2),3))])

# How many "fully in-crisis" windows are there? (W=252 ⇒ window start ≥ 2024-01-01)
# That means window-end ≥ 2024-01-01 + 252 trading days ≈ 2024-12-31
fully_in <- roll[date >= "2024-12-31" & date <= "2026-05-22"]
cat(sprintf("\nFully-in-crisis windows (end >= 2024-12-31): n=%d, median=%g, range=[%d,%d]\n",
            nrow(fully_in), median(fully_in$count_above_mp_ex_market),
            min(fully_in$count_above_mp_ex_market), max(fully_in$count_above_mp_ex_market)))

# ============================================================================
# STAGE 8: BOOTSTRAP CI verification
# ============================================================================
sep("STAGE 8: bootstrap CI degeneracy check")
roll_ci <- fread("output/bootstrap_rolling_ci.csv")
cat(sprintf("rows = %d\n", nrow(roll_ci)))
cat("Distribution of per-date q025, q500, q975 across windows:\n")
print(roll_ci[, .(q025_unique = length(unique(q025)), q025_min = min(q025), q025_max = max(q025),
                  q500_unique = length(unique(q500)), q500_min = min(q500), q500_max = max(q500),
                  q975_unique = length(unique(q975)), q975_min = min(q975), q975_max = max(q975))])
t2b <- fread("output/bootstrap_t2b.csv"); print(t2b)

# ============================================================================
# STAGE 9: BAI-PERRON
# ============================================================================
sep("STAGE 9: Bai-Perron break dates")
bp <- fread("output/bai_perron_breaks.csv"); print(bp)
cat("(All five breaks were FORCED — script uses `breaks=5` not BIC selection.)\n")

# ============================================================================
# STAGE 10: DCC fits
# ============================================================================
sep("STAGE 10: DCC fits and S1/S2/S3")
dcc <- fread("output/dcc_params.csv"); print(dcc)
cat("\nNote: cleaned DCC b=0.98 is the user-imposed CAP, not free optimum.\n")
s1 <- fread("output/s1_variance_reduction.csv"); print(s1)
cat("\nNote: 'variance' here is variance of the conditional correlation series,\n")
cat("NOT out-of-sample portfolio variance. EW2019's claim is about the latter.\n")
s3 <- fread("output/s3_cocoa_basis_decoupling.csv"); print(s3)
cat("\n========================================\n")
cat("DIAGNOSTICS COMPLETE\n")
cat("========================================\n")
