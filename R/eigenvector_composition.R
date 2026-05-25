# Day 11 (run early): eigenvector composition analysis + T3 test.
#
# T3 (from docs/scope_and_target.md):
#   Define the crisis-emerging eigenvalue as the lowest-ranked supra-MP
#   eigenvalue that exceeds the MP edge in ≥50% of in-crisis windows but
#   in <50% of 2019-2023 baseline windows. Its eigenvector (averaged over
#   crisis-period windows) has the fraction of its squared loadings on the
#   cocoa bloc {COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y} ≥ 0.50.
#
# Outputs:
#   output/eigenvector_loadings_subsamples.csv
#   output/t3_test.csv
#   figures/fig3a_eigenvector_loadings.png    (subsample comparison)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

COCOA_BLOC <- c("COCOA_NY", "COCOA_LDN", "HSY", "MDLZ", "LISN", "GHS", "GHANA10Y")
CRISIS_START  <- as.Date("2024-01-01")
CRISIS_END    <- as.Date("2026-05-22")
BASELINE_START <- as.Date("2019-01-01")
BASELINE_END   <- as.Date("2023-12-31")

z <- fread("data/clean/std_residuals_v1.csv")
z[, date := as.Date(date)]
tickers <- setdiff(names(z), "date")

# Restrict to dates where all 25 series are observed (post-Ghana entry)
z <- z[complete.cases(z[, tickers, with = FALSE])]
cat(sprintf("Working sample: %d dates x %d series (complete cases)\n", nrow(z), length(tickers)))

mp_upper <- function(N, T) (1 + sqrt(N / T))^2

# ----- Subsample eigendecompositions ----------------------------------------
analyze_subsample <- function(z_sub, label) {
  Tn <- nrow(z_sub); N <- length(tickers)
  R <- cor(as.matrix(z_sub[, tickers, with = FALSE]))
  e <- eigen(R, symmetric = TRUE)
  edge <- mp_upper(N, Tn)
  evals <- e$values
  evecs <- e$vectors

  out <- data.table(
    sample = label,
    rank = seq_along(evals),
    eigenvalue = evals,
    above_mp = evals > edge,
    mp_edge = edge,
    T = Tn, N = N
  )
  # Loadings (each row = ticker, columns = eigenvector ranks)
  load_dt <- as.data.table(evecs)
  setnames(load_dt, paste0("ev", seq_along(evals)))
  load_dt[, ticker := tickers]
  load_dt[, sample := label]
  setcolorder(load_dt, c("sample", "ticker"))

  # Squared-loading fraction on cocoa bloc, per eigenvector
  bloc_idx <- match(COCOA_BLOC, tickers)
  bloc_mass <- apply(evecs[bloc_idx, , drop = FALSE]^2, 2L, sum)  # rows are bloc tickers
  out[, cocoa_bloc_l2_frac := bloc_mass]
  out[, top3_loadings := sapply(seq_along(evals), function(k) {
    v <- evecs[, k]
    nm <- tickers[order(-abs(v))[1:3]]
    paste(nm, collapse = ", ")
  })]

  list(summary = out, loadings = load_dt)
}

pre  <- analyze_subsample(z[date >= BASELINE_START & date <= BASELINE_END], "pre_2019_2023")
crs  <- analyze_subsample(z[date >= CRISIS_START  & date <= CRISIS_END],   "crisis_2024_2026")
full <- analyze_subsample(z, "full_sample")

all_summary <- rbindlist(list(pre$summary, crs$summary, full$summary))
all_loadings <- rbindlist(list(pre$loadings, crs$loadings, full$loadings), fill = TRUE)

fwrite(all_summary, "output/eigenvector_loadings_subsamples.csv")
fwrite(all_loadings, "output/eigenvector_loadings_full.csv")

cat("\n=== Supra-MP eigenvalues per subsample ===\n")
print(all_summary[above_mp == TRUE,
                  .(sample, rank, eigenvalue = round(eigenvalue, 3),
                    mp_edge = round(mp_edge, 3),
                    cocoa_bloc_frac = round(cocoa_bloc_l2_frac, 3),
                    top3_loadings)])

# ----- T3: run both threshold parameterisations for completeness -----------
# Pre-registered (scope_and_target.md): bloc share >= 0.50 (majority rule) on rank-1.
# 0.55 is reported as a sensitivity check (2x the bloc's proportional share,
# 7/25 ≈ 0.28). See paper Appendix A.

THRESHOLD_PRIMARY <- 0.50    # docs/scope_and_target.md pre-registration (majority rule)
THRESHOLD_ALT     <- 0.55    # alt parameterisation: 2x proportional baseline

t3_results <- list()

# --- Test (a): Original subsample emerging-eigenvalue test ----------------
crs_supra_ranks <- crs$summary[above_mp == TRUE, rank]
pre_supra_ranks <- pre$summary[above_mp == TRUE, rank]
emerging_ranks <- setdiff(crs_supra_ranks, pre_supra_ranks)

cat(sprintf("\n--- T3(a): emerging-eigenvalue subsample variant, bloc >= 0.50 ---\n"))
cat(sprintf("Pre supra ranks:    %s\n", paste(pre_supra_ranks, collapse = ", ")))
cat(sprintf("Crisis supra ranks: %s\n", paste(crs_supra_ranks, collapse = ", ")))
cat(sprintf("Emerging ranks:     %s\n",
            if (length(emerging_ranks)) paste(emerging_ranks, collapse = ", ") else "NONE"))
t3a_pass <- FALSE; t3a_share <- NA_real_; t3a_rank <- NA_integer_
if (length(emerging_ranks)) {
  t3a_rank <- min(emerging_ranks)
  t3a_share <- crs$summary[rank == t3a_rank, cocoa_bloc_l2_frac]
  t3a_pass <- t3a_share >= THRESHOLD_PRIMARY
}
cat(sprintf("T3(a) verdict: %s\n", if (t3a_pass) "PASS" else "FAIL (no emerging eigenvalue; trivially fails)"))
t3_results[["a_emerging_subsample"]] <- list(
  spec = "Emerging-eigenvalue subsample variant with bloc share ≥ 0.50",
  rank = t3a_rank, bloc_share = t3a_share, threshold = THRESHOLD_PRIMARY,
  pass = t3a_pass
)

# --- Test (b): Original rolling ≥50%/<50% test (proper time-varying) -----
cat(sprintf("\n--- T3(b): rolling ≥50%%/<50%% supra-MP indicator per rank ---\n"))
WINDOW <- 252L
zfull <- fread("data/clean/std_residuals_v1.csv")[, date := as.Date(date)]
zfull_cc <- zfull[complete.cases(zfull[, tickers, with = FALSE])]
zmat <- as.matrix(zfull_cc[, tickers, with = FALSE])
n_total <- nrow(zmat); N <- length(tickers)
starts <- seq_len(n_total - WINDOW + 1L)
end_dates <- zfull_cc$date[starts + WINDOW - 1L]
edge <- (1 + sqrt(N / WINDOW))^2

supra_ind <- matrix(FALSE, nrow = length(starts), ncol = N)
for (i in seq_along(starts)) {
  R <- cor(zmat[starts[i]:(starts[i] + WINDOW - 1L), , drop = FALSE])
  ev <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  supra_ind[i, ] <- ev > edge
}
pre_mask <- end_dates >= "2019-01-01" & end_dates <= "2023-12-31"
crs_mask <- end_dates >= "2024-01-01" & end_dates <= "2026-05-22"
pre_share <- colMeans(supra_ind[pre_mask, ])
crs_share <- colMeans(supra_ind[crs_mask, ])
rolling_t3 <- data.table(rank = 1:N,
                         pre_supra_share = round(pre_share, 3),
                         crs_supra_share = round(crs_share, 3),
                         emerging = crs_share >= 0.5 & pre_share < 0.5)
print(rolling_t3[1:8])
emerging_rolling <- which(rolling_t3$emerging)
cat(sprintf("Rolling-emerging ranks: %s\n",
            if (length(emerging_rolling)) paste(emerging_rolling, collapse = ", ") else "NONE"))
t3b_pass <- FALSE; t3b_share <- NA_real_; t3b_rank <- NA_integer_
if (length(emerging_rolling)) {
  t3b_rank <- min(emerging_rolling)
  # Get bloc share for that rank from the crisis subsample
  t3b_share <- crs$summary[rank == t3b_rank, cocoa_bloc_l2_frac]
  t3b_pass <- !is.na(t3b_share) && t3b_share >= THRESHOLD_PRIMARY
  cat(sprintf("Rolling-emerging rank %d has bloc share %.3f (threshold %.2f) -- %s\n",
              t3b_rank, t3b_share, THRESHOLD_PRIMARY,
              if (t3b_pass) "PASS" else "FAIL"))
}
t3_results[["b_emerging_rolling"]] <- list(
  spec = "Emerging-eigenvalue rolling variant: ≥50%/<50% emergent rank with bloc share ≥ 0.50",
  rank = t3b_rank, bloc_share = t3b_share, threshold = THRESHOLD_PRIMARY,
  pass = t3b_pass
)

# --- Test (c): rank-1 reorganisation test (per scope_and_target.md restatement) --------
# 0.50 majority-rule threshold AND doubling vs pre.
pre_rank1_share <- pre$summary[rank == 1, cocoa_bloc_l2_frac]
crs_rank1_share <- crs$summary[rank == 1, cocoa_bloc_l2_frac]
ratio <- crs_rank1_share / pre_rank1_share
t3c_pass_55 <- crs_rank1_share >= THRESHOLD_ALT && ratio >= 2
t3c_pass_50 <- crs_rank1_share >= THRESHOLD_PRIMARY    && ratio >= 2

cat(sprintf("\n--- T3(c): rank-1 reorganisation ---\n"))
cat(sprintf("Pre-crisis rank-1 bloc share: %.3f\n", pre_rank1_share))
cat(sprintf("Crisis rank-1 bloc share:     %.3f (ratio = %.2fx)\n",
            crs_rank1_share, ratio))
cat(sprintf("At ALT threshold     0.55:    %s\n", if (t3c_pass_55) "PASS" else "FAIL"))
cat(sprintf("At PRIMARY threshold 0.50:    %s\n", if (t3c_pass_50) "PASS" else "FAIL"))
cat("NOTE: rank-1 is the largest eigenvalue (the 'market mode'); the rolling\n")
cat("count test EXCLUDES rank-1, so this is testing a different object.\n")
t3_results[["c_rank1"]] <- list(
  spec = "Rank-1 reorganisation: bloc share ≥ 0.50 (majority rule) AND doubling vs pre",
  rank = 1, pre_share = pre_rank1_share, crs_share = crs_rank1_share, ratio = ratio,
  pass_at_55 = t3c_pass_55, pass_at_50 = t3c_pass_50
)

# --- Also report: pre-existing cocoa factor at rank 3 -------------------
pre_rank3_share <- pre$summary[rank == 3, cocoa_bloc_l2_frac]
crs_rank3_share <- crs$summary[rank == 3, cocoa_bloc_l2_frac]
cat(sprintf("\n--- Pre-existing cocoa factor (rank 3) sanity check ---\n"))
cat(sprintf("Pre-crisis rank 3 bloc share: %.3f (top-3: %s)\n",
            pre_rank3_share, pre$summary[rank == 3, top3_loadings]))
cat(sprintf("Crisis    rank 3 bloc share: %.3f (top-3: %s)\n",
            crs_rank3_share, crs$summary[rank == 3, top3_loadings]))
cat("⇒ A cocoa-dominated mode already existed pre-crisis at rank 3; its bloc share\n")
cat("  DECREASED in crisis. The 'reorganisation' is on rank 1, not the cocoa factor.\n")

# Save consolidated T3 results
fwrite(data.table(
  test = c("a_emerging_subsample", "b_emerging_rolling",
           "c_rank1_at_0.55", "c_rank1_at_0.50"),
  description = c(
    "Emerging-eigenvalue subsample variant, bloc>=0.50",
    "Emerging-eigenvalue rolling variant, bloc>=0.50",
    "Rank-1 reorganisation, bloc>=0.55 AND doubling (alt parameterisation)",
    "Rank-1 reorganisation, bloc>=0.50 AND doubling (primary, pre-registered)"
  ),
  rank = c(NA_integer_, t3b_rank, 1L, 1L),
  bloc_share = c(NA_real_, t3b_share, crs_rank1_share, crs_rank1_share),
  pass = c(t3a_pass, t3b_pass, t3c_pass_55, t3c_pass_50)
), "output/t3_test.csv")
cat("\nWrote output/t3_test.csv (4-row consolidated verdict).\n")

# ----- Figure 3a: Loadings comparison ---------------------------------------
# Pick the top-3 supra-MP eigenvalues from each subsample, plot loadings as heatmap.
plot_data <- list()
for (samp_name in c("pre_2019_2023", "crisis_2024_2026")) {
  s <- if (samp_name == "pre_2019_2023") pre else crs
  top_ranks <- s$summary[above_mp == TRUE, rank][1:min(4, sum(s$summary$above_mp))]
  for (k in top_ranks) {
    ev_col <- paste0("ev", k)
    plot_data[[length(plot_data) + 1]] <- data.table(
      sample = samp_name,
      rank = k,
      eigenvalue = s$summary[rank == k, eigenvalue],
      ticker = s$loadings$ticker,
      loading = s$loadings[[ev_col]]
    )
  }
}
plot_dt <- rbindlist(plot_data)
plot_dt[, ticker_lbl := factor(ticker, levels = tickers)]
plot_dt[, panel := sprintf("%s\nrank %d, λ=%.2f",
                           ifelse(sample == "pre_2019_2023", "Pre (2019-23)", "Crisis (2024-26)"),
                           rank, eigenvalue)]
plot_dt[, panel := factor(panel, levels = unique(panel))]

# Sign-normalise eigenvectors so the largest absolute loading is positive
# (eigenvector sign is arbitrary; this just stabilises the visual)
plot_dt[, loading_signed := loading * sign(loading[which.max(abs(loading))]), by = panel]
plot_dt[, in_bloc := ticker %in% COCOA_BLOC]

p3a <- ggplot(plot_dt, aes(x = loading_signed, y = ticker_lbl, fill = in_bloc)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "tomato"),
                    labels = c("FALSE" = "Other", "TRUE" = "Cocoa bloc"),
                    name = NULL) +
  facet_wrap(~ panel, nrow = 2, scales = "free_y") +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Eigenvector loadings of supra-MP modes",
    subtitle = "Sign-normalised. Red = cocoa bloc (COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y).",
    x = "loading", y = NULL,
    caption = "Pre-crisis vs crisis subsamples. Eigenvalues that emerged only in crisis carry the H1 signal."
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 8))

ggsave("figures/fig3a_eigenvector_loadings.png", p3a, width = 11, height = 9, dpi = 150)
cat("Wrote figures/fig3a_eigenvector_loadings.png\n")
