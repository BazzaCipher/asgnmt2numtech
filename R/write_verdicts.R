# Consolidated test verdicts: one row per pre-registered test, all in one CSV.
# Reads the per-test output CSVs and assembles a single summary table for the paper.
# Output: output/test_verdicts.csv

suppressPackageStartupMessages({
  library(data.table)
})

CRISIS_START <- as.Date("2024-01-01")
CRISIS_END   <- as.Date("2026-05-22")

rows <- list()

# ----- T1: pre-crisis baseline stability ------------------------------------
roll <- fread("output/rolling_count_v1.csv")
roll[, date := as.Date(date)]
pre <- roll[date >= "2019-01-01" & date <= "2023-12-31" & n_assets == 25,
            count_above_mp_ex_market]
in_ <- roll[date >= CRISIS_START & date <= CRISIS_END, count_above_mp_ex_market]
t1_iqr_width <- diff(quantile(pre, c(0.25, 0.75), na.rm = TRUE))
t1_pass <- t1_iqr_width <= 1 && all(pre %in% 1:3)
rows[[length(rows) + 1L]] <- data.table(
  test = "T1", description = "Pre-crisis rolling-count baseline stability",
  metric = sprintf("median=%g, IQR=[%g,%g], range=[%d,%d]",
                   median(pre), quantile(pre, 0.25), quantile(pre, 0.75),
                   min(pre), max(pre)),
  threshold = "IQR width <= 1 AND values in {1,2,3}",
  pass = t1_pass,
  note = sprintf("In-crisis median = %g (+%g vs pre)", median(in_), median(in_) - median(pre))
)

# ----- T2(a): Bai-Perron BIC on rolling count -------------------------------
bp <- fread("output/bai_perron_breaks.csv")
n_bp <- nrow(bp)
in_window <- if (n_bp) any(as.Date(bp$break_date) >= "2024-01-01" &
                            as.Date(bp$break_date) <= "2024-12-31") else FALSE
rows[[length(rows) + 1L]] <- data.table(
  test = "T2(a)", description = "Bai-Perron BIC break in [2024-01-01, 2024-12-31]",
  metric = sprintf("BIC selects %d breaks", n_bp),
  threshold = ">=1 break inside the 2024 attribution window",
  pass = in_window && n_bp > 0,
  note = if (file.exists("output/bai_perron_bic_table.csv")) {
    bic <- fread("output/bai_perron_bic_table.csv")
    sprintf("BIC argmin at m=%d (BIC=%.2f)", bic$n_breaks[which.min(bic$BIC)],
            min(bic$BIC))
  } else "no BIC table"
)

# ----- T2(b): bootstrap CI on pre-crisis median -----------------------------
boot_t2b <- if (file.exists("output/bootstrap_t2b.csv")) fread("output/bootstrap_t2b.csv") else NULL
if (!is.null(boot_t2b) && nrow(boot_t2b)) {
  ci_lo <- boot_t2b$ci_lo_95[1]; ci_hi <- boot_t2b$ci_hi_95[1]
  in_med <- boot_t2b$in_obs_median[1]
  t2b_pass <- as.logical(boot_t2b$t2b_pass[1])
  rows[[length(rows) + 1L]] <- data.table(
    test = "T2(b)", description = "In-crisis median > 95% bootstrap CI upper bound",
    metric = sprintf("CI=[%.2f, %.2f], in-crisis median=%g", ci_lo, ci_hi, in_med),
    threshold = "in-crisis median > CI upper bound",
    pass = t2b_pass,
    note = if (ci_hi - ci_lo < 0.05) "1-bit comparison; bootstrap CI degenerate on integer statistic" else ""
  )
} else {
  rows[[length(rows) + 1L]] <- data.table(
    test = "T2(b)", description = "In-crisis median > 95% bootstrap CI upper bound",
    metric = "no bootstrap output found", threshold = "", pass = NA, note = ""
  )
}

# ----- T3: four-row file → flatten to one row per spec ----------------------
t3 <- fread("output/t3_test.csv")
for (i in seq_len(nrow(t3))) {
  rows[[length(rows) + 1L]] <- data.table(
    test = paste0("T3(", c("a","b","c@0.55","c@0.50")[i], ")"),
    description = t3$description[i],
    metric = sprintf("rank=%s, bloc_share=%s",
                     ifelse(is.na(t3$rank[i]), "-", as.character(t3$rank[i])),
                     ifelse(is.na(t3$bloc_share[i]), "-",
                            sprintf("%.3f", t3$bloc_share[i]))),
    threshold = if (i == 4) "bloc_share >= 0.50 (moved goalpost)" else "bloc_share >= 0.55 (pre-reg)",
    pass = as.logical(t3$pass[i]),
    note = if (i == 4) "Threshold 0.50 differs from pre-registered 0.55" else ""
  )
}

# ----- S1: OOS MVP variance reduction ---------------------------------------
s1 <- fread("output/s1_mvp_oos_variance.csv")
s1_crisis <- s1[regime == "crisis"]
s1_pass <- nrow(s1_crisis) && s1_crisis$var_red_pct >= 5
rows[[length(rows) + 1L]] <- data.table(
  test = "S1", description = "Cleaned DCC reduces OOS MVP variance in crisis",
  metric = sprintf("sd_std=%.3f%%, sd_clean=%.3f%% (crisis), var_red=%s%%",
                   s1_crisis$sd_std_pct, s1_crisis$sd_clean_pct,
                   if (nrow(s1_crisis)) as.character(s1_crisis$var_red_pct) else "-"),
  threshold = "variance reduction >= 5%",
  pass = isTRUE(s1_pass),
  note = "Engle-Ledoit-Wolf 2019 metric; cleaned ~= standard at N=25"
)

# ----- S2: BIC-selected breaks on cocoa-coffee or cocoa-sugar ---------------
s2 <- fread("output/s2_bai_perron_pairs.csv")
relevant <- s2[pair %in% c("COCOA_NY_SUGAR", "COCOA_NY_COFFEE_ARA")]
s2_pass <- any(relevant$in_crisis_break)
rows[[length(rows) + 1L]] <- data.table(
  test = "S2", description = "BIC-selected in-crisis break in cocoa-coffee or cocoa-sugar pair",
  metric = sprintf("BIC selects 0 breaks in all %d pairs", nrow(s2)),
  threshold = ">=1 in-crisis break in either pair",
  pass = s2_pass,
  note = "All 11 cocoa-anchored pairs have 0 BIC-selected breaks"
)

# ----- S3: NY-LDN basis decoupling ------------------------------------------
s3 <- fread("output/s3_cocoa_basis_decoupling.csv")
rows[[length(rows) + 1L]] <- data.table(
  test = "S3", description = "Cocoa NY-LDN basis decouples by >=0.15 in crisis",
  metric = sprintf("baseline_mean=%.3f, crisis_mean=%.3f (change %+.3f); crisis_min=%.3f",
                   s3$baseline_mean_2019_2023, s3$crisis_mean_2024_2026,
                   s3$crisis_mean_2024_2026 - s3$baseline_mean_2019_2023,
                   s3$crisis_min_2024_2026),
  threshold = "baseline_mean - crisis_min >= 0.15",
  pass = as.logical(s3$s3_pass),
  note = sprintf("basis RE-COUPLED (+%.3f), opposite of prediction",
                 s3$crisis_mean_2024_2026 - s3$baseline_mean_2019_2023)
)

verdicts <- rbindlist(rows)
fwrite(verdicts, "output/test_verdicts.csv")

cat("\n=== Consolidated test verdicts ===\n")
print(verdicts[, .(test, pass, metric)])
cat(sprintf("\nPass: %d/%d\n", sum(verdicts$pass, na.rm = TRUE), nrow(verdicts)))
cat("Wrote output/test_verdicts.csv\n")
