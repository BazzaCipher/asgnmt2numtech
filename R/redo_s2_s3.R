# Redo S2 (BIC Bai-Perron per pair) + S3 (full crisis window) from already-saved
# dcc_pairwise_correlations.csv. Much faster than re-running fit_dcc.R.

suppressPackageStartupMessages({
  library(data.table)
  library(strucchange)
})

CRISIS_START <- as.Date("2024-01-01")
CRISIS_END   <- as.Date("2026-05-22")
BASE_START   <- as.Date("2019-01-01")
BASE_END     <- as.Date("2023-12-31")

corr_dt <- fread("output/dcc_pairwise_correlations.csv")
corr_dt[, date := as.Date(date)]
cat(sprintf("Loaded %d rows, %d pairs\n", nrow(corr_dt), length(unique(corr_dt$pair))))

# ===== S2 ====================================================================
s2_rows <- list()
for (pair_name in unique(corr_dt$pair)) {
  pair_dt <- corr_dt[pair == pair_name, .(date, rho = rho_clean)]
  # h=0.15 ⇒ m_max=5, BIC selects ≤ 5; fast.
  bp_full <- tryCatch(
    breakpoints(rho ~ 1, data = pair_dt, h = 0.15, breaks = 5),
    error = function(e) NULL
  )
  bp_bic <- if (!is.null(bp_full)) which.min(BIC(bp_full)) - 1L else 0L
  break_dates <- if (!is.null(bp_full) && bp_bic > 0L) {
    bp <- breakpoints(bp_full, breaks = bp_bic)
    pair_dt$date[bp$breakpoints]
  } else as.Date(character(0))
  in_crisis <- any(break_dates >= CRISIS_START & break_dates <= CRISIS_END)
  s2_rows[[pair_name]] <- data.table(
    pair = pair_name,
    n_breaks = length(break_dates),
    break_dates = paste(break_dates, collapse = "; "),
    in_crisis_break = in_crisis
  )
  cat(sprintf("  [%s] BIC selects %d breaks\n", pair_name, bp_bic))
}
s2 <- rbindlist(s2_rows)
fwrite(s2, "output/s2_bai_perron_pairs.csv")
cat("\n=== S2 (BIC-selected breaks per pair) ===\n")
print(s2)
cat(sprintf("S2 (in-crisis break on cocoa-coffee OR cocoa-sugar): %s\n",
            if (any(s2[pair %in% c("COCOA_NY_SUGAR", "COCOA_NY_COFFEE_ARA"),
                        in_crisis_break])) "PASS" else "FAIL"))

# ===== S3 ====================================================================
basis <- corr_dt[pair == "COCOA_NY_COCOA_LDN"]
baseline_mean <- mean(basis[date >= BASE_START & date <= BASE_END, rho_clean])
crisis_rho   <- basis[date >= CRISIS_START & date <= CRISIS_END, rho_clean]
crisis_mean  <- mean(crisis_rho)
crisis_min   <- min(crisis_rho)
crisis_max   <- max(crisis_rho)
basis_drop_from_mean <- baseline_mean - crisis_min
s3_pass <- basis_drop_from_mean >= 0.15

cat("\n=== S3: COCOA_NY × COCOA_LDN basis decoupling (cleaned DCC) ===\n")
cat(sprintf("Baseline mean (2019-2023):       %.3f\n", baseline_mean))
cat(sprintf("Crisis (2024-2026) mean:         %.3f  (change %+.3f)\n",
            crisis_mean, crisis_mean - baseline_mean))
cat(sprintf("Crisis (2024-2026) range:        [%.3f, %.3f]\n", crisis_min, crisis_max))
cat(sprintf("Drop (baseline mean − crisis min): %+.3f  (need ≥0.15 to pass)\n",
            basis_drop_from_mean))
cat(sprintf("S3: %s — basis %s during crisis\n",
            if (s3_pass) "PASS" else "FAIL",
            if (crisis_mean > baseline_mean) "RE-COUPLED" else "decoupled"))

fwrite(data.table(
  baseline_mean_2019_2023 = baseline_mean,
  crisis_mean_2024_2026   = crisis_mean,
  crisis_min_2024_2026    = crisis_min,
  crisis_max_2024_2026    = crisis_max,
  basis_drop = basis_drop_from_mean,
  s3_pass = s3_pass
), "output/s3_cocoa_basis_decoupling.csv")
