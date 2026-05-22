# Extra checks: items I noticed only after reading fit_dcc.R + loadings CSV.

suppressPackageStartupMessages({
  library(data.table)
  library(strucchange)
})
options(width = 140)
sep <- function(s) cat(sprintf("\n========== %s ==========\n", s))

# ============================================================================
# E1. DCC log-likelihood comparison: cleaned vs standard
# ============================================================================
sep("E1. DCC log-likelihood: cleaned target is WORSE than standard")
dcc <- fread("output/dcc_params.csv")
print(dcc)
cat(sprintf("\nL(standard) = %.1f,  L(cleaned) = %.1f\n", -dcc$nll[1], -dcc$nll[2]))
cat(sprintf("ΔL (clean − std) = %.1f  (negative ⇒ cleaning hurts fit)\n",
            (-dcc$nll[2]) - (-dcc$nll[1])))
cat("Same # params (a,b). With ΔL ≈ -780 the cleaned model is dominated\n")
cat("on AIC/BIC/likelihood-ratio. The paper does not report this.\n")

# ============================================================================
# E2. Cocoa basis correlation FULL 2024-2026 path (S3 used 2025-2026 only)
# ============================================================================
sep("E2. Cocoa basis (NY × LDN) DCC correlation across windows")
corr <- fread("output/dcc_pairwise_correlations.csv")
corr[, date := as.Date(date)]
basis <- corr[pair == "COCOA_NY_COCOA_LDN"]

# Baseline mean
bsl <- mean(basis[date >= "2019-01-01" & date <= "2023-12-31", rho_clean])

# Different crisis windows
windows <- list(
  full_2024_2026 = c("2024-01-01", "2026-05-22"),
  s3_paper_2025_2026 = c("2025-01-01", "2026-05-22"),
  spike_2024H1 = c("2024-01-01", "2024-06-30"),
  decline_2025H2 = c("2025-07-01", "2026-05-22")
)
out <- rbindlist(lapply(names(windows), function(nm) {
  w <- windows[[nm]]
  b <- basis[date >= as.Date(w[1]) & date <= as.Date(w[2]), rho_clean]
  data.table(window = nm, n = length(b),
             mean = round(mean(b), 4), min = round(min(b), 4),
             max = round(max(b), 4),
             drop_from_bsl_mean = round(bsl - mean(b), 4),
             drop_from_bsl_min  = round(bsl - min(b),  4))
}))
cat(sprintf("Baseline mean (2019-2023): %.4f\n\n", bsl))
print(out)
cat("\nNote: S3 used 2025-2026 (excluded 2024) — no methodological reason given.\n")
cat("Including 2024 changes the numbers materially.\n")

# ============================================================================
# E3. BIC-selected Bai-Perron on the rolling count series
# ============================================================================
sep("E3. Bai-Perron with BIC selection (vs forced breaks=5)")
roll <- fread("output/rolling_count_v1.csv")
roll[, date := as.Date(date)]
plot_dt <- roll[!is.na(count_above_mp_ex_market)]
bp_input <- plot_dt[, .(date, count = count_above_mp_ex_market)]

# Without forcing m, use BIC selection
bp <- breakpoints(count ~ 1, data = bp_input, h = 0.05)
cat("BIC table (number of breaks → BIC):\n")
bic_vec <- BIC(bp)
bic_table <- data.table(n_breaks = seq_along(bic_vec) - 1L,
                        BIC = round(as.numeric(bic_vec), 1))
print(bic_table)
opt_m <- which.min(bic_vec) - 1
cat(sprintf("\nBIC-optimal number of breaks: %d\n", opt_m))
if (opt_m > 0) {
  bp_opt <- breakpoints(bp, breaks = opt_m)
  break_dates <- bp_input$date[bp_opt$breakpoints]
  cat(sprintf("BIC-optimal break dates: %s\n", paste(break_dates, collapse = ", ")))
  in_crisis <- break_dates[break_dates >= "2024-01-01" & break_dates <= "2024-12-31"]
  cat(sprintf("Breaks inside pre-registered T2(a) window [2024-01-01, 2024-12-31]: %s\n",
              if (length(in_crisis)) paste(in_crisis, collapse = ", ") else "NONE"))
} else {
  cat("BIC selects ZERO breaks — T2(a) cannot be assessed against a null model.\n")
}

# ============================================================================
# E4. T3 re-checked with proper rolling-windows test
# ============================================================================
sep("E4. T3 rolling-windows test (≥50% / <50%)")
# Per docs/scope_and_target.md, T3 was supposed to identify eigenvalues that
# cross the MP edge in ≥50% of in-crisis windows but <50% of pre-crisis.
# Need to recompute per-window eigendecomposition and track per-rank supra-MP indicator.
WINDOW <- 252L
z <- fread("data/clean/std_residuals_v1.csv"); z[, date := as.Date(date)]
tickers <- setdiff(names(z), "date")
z_cc <- z[complete.cases(z[, tickers, with = FALSE])]
N <- length(tickers)

n_total <- nrow(z_cc)
starts <- seq_len(n_total - WINDOW + 1L)
window_end_dates <- z_cc$date[starts + WINDOW - 1L]
edge <- (1 + sqrt(N / WINDOW))^2

# For each window: compute all N eigenvalues, supra-MP indicator per rank
supra_indicator <- matrix(FALSE, nrow = length(starts), ncol = N)
zmat <- as.matrix(z_cc[, tickers, with = FALSE])
for (i in seq_along(starts)) {
  w <- zmat[starts[i]:(starts[i] + WINDOW - 1L), , drop = FALSE]
  R <- cor(w)
  ev <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  supra_indicator[i, ] <- ev > edge
}
pre_mask <- window_end_dates >= "2019-01-01" & window_end_dates <= "2023-12-31"
crs_mask <- window_end_dates >= "2024-01-01" & window_end_dates <= "2026-05-22"
pre_share <- colMeans(supra_indicator[pre_mask, ])
crs_share <- colMeans(supra_indicator[crs_mask, ])
t3_rolling <- data.table(rank = 1:N,
                         pre_supra_share = round(pre_share, 3),
                         crs_supra_share = round(crs_share, 3),
                         emerging = crs_share >= 0.5 & pre_share < 0.5)
print(t3_rolling)
emerging_ranks <- which(t3_rolling$emerging)
cat(sprintf("\nEmerging ranks (≥50%% crisis, <50%% pre): %s\n",
            if (length(emerging_ranks)) paste(emerging_ranks, collapse = ", ") else "NONE"))

# ============================================================================
# E5. Cocoa-bloc rank 3 was a cocoa factor PRE-CRISIS
# ============================================================================
sep("E5. Cocoa factor already existed pre-crisis at rank 3")
eig <- fread("output/eigenvector_loadings_subsamples.csv")
cat("Pre-crisis rank 3:  eigenvalue=%.2f, bloc share=%.3f, top-3=%s\n" |>
    sprintf(eig[sample == "pre_2019_2023" & rank == 3, eigenvalue],
            eig[sample == "pre_2019_2023" & rank == 3, cocoa_bloc_l2_frac],
            eig[sample == "pre_2019_2023" & rank == 3, top3_loadings]))
cat("Crisis rank 3:      eigenvalue=%.2f, bloc share=%.3f, top-3=%s\n" |>
    sprintf(eig[sample == "crisis_2024_2026" & rank == 3, eigenvalue],
            eig[sample == "crisis_2024_2026" & rank == 3, cocoa_bloc_l2_frac],
            eig[sample == "crisis_2024_2026" & rank == 3, top3_loadings]))
cat("Crisis rank 5:      eigenvalue=%.2f, bloc share=%.3f, top-3=%s\n" |>
    sprintf(eig[sample == "crisis_2024_2026" & rank == 5, eigenvalue],
            eig[sample == "crisis_2024_2026" & rank == 5, cocoa_bloc_l2_frac],
            eig[sample == "crisis_2024_2026" & rank == 5, top3_loadings]))
cat("⇒ The 'cocoa factor' is not new in crisis — rank 3 pre-crisis already loads\n")
cat("  COCOA_NY, COCOA_LDN, CORN with 0.70 bloc share. What changes is RANK 1.\n")

# ============================================================================
# E6. DCC params: actual cap vs. paper claim
# ============================================================================
sep("E6. DCC constraint reality check")
cat("Paper §3.2 (paper.typ:64): 'a + b < 0.98'\n")
cat("Script (fit_dcc.R:96):     'a + b >= 0.999' returns +Inf (so a+b<0.999)\n")
cat("Script (fit_dcc.R:130):    upper = c(0.5, 0.98)  (b cap, not a+b)\n")
cat("Realised cleaned fit:      a=0.006, b=0.98, a+b=0.986  (b at cap)\n")
cat("⇒ Paper's stated constraint is violated by the actual fit.\n")

cat("\n========================================\n")
cat("EXTRA DIAGNOSTICS COMPLETE\n")
cat("========================================\n")
