# Days 9-10 (early): stationary block bootstrap of the rolling MP-edge count.
#
# T2(b): the in-crisis (2024-2026) median count must exceed the upper bound of
# the 95% bootstrap CI on the pre-crisis (2019-2023, N=25 windows only) median
# count.
#
# Method:
#   - Stationary bootstrap (Politis-Romano, geometric block length, mean=25) of
#     the pre-crisis standardised residuals. Block length 25 ≈ one trading month.
#   - For each replication, recompute the rolling 252-day correlation matrix,
#     count supra-MP eigenvalues (ex-market), take the median over the rep's
#     count series. → distribution of pre-crisis median counts → 95% CI.
#   - Compare to observed in-crisis median.
#
# Also produces a per-date 95% bootstrap CI band on the full pre-crisis count
# trajectory (visual sanity check).
#
# Outputs:
#   output/bootstrap_t2b.csv       (point estimates + CI + verdict)
#   output/bootstrap_rolling_ci.csv (per-date bootstrap quantiles)
#   figures/fig2b_count_with_ci.png (augmented Figure 2)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(parallel)
})

set.seed(20260522)

B           <- 1000L    # bootstrap replications (proposal target)
BLOCK_MEAN  <- 25L      # mean geometric block length
NCORES      <- max(1L, min(20L, detectCores() - 4L))
WINDOW      <- 252L
BASE_START  <- as.Date("2019-01-01")
BASE_END    <- as.Date("2023-12-31")
CRISIS_START <- as.Date("2024-01-01")
CRISIS_END   <- as.Date("2026-05-22")

# ----- Load data ------------------------------------------------------------
z <- fread("data/clean/std_residuals_v1.csv")
z[, date := as.Date(date)]
tickers <- setdiff(names(z), "date")

# Restrict to complete-cases (post-Ghana entry, N=25)
z_cc <- z[complete.cases(z[, tickers, with = FALSE])]
cat(sprintf("Complete-cases panel: %d dates (%s to %s)\n",
            nrow(z_cc), min(z_cc$date), max(z_cc$date)))

# Observed rolling counts (already computed in mp_analysis.R)
roll <- fread("output/rolling_count_v1.csv")
roll[, date := as.Date(date)]

pre_obs_counts <- roll[date >= BASE_START & date <= BASE_END & n_assets == 25,
                       count_above_mp_ex_market]
in_obs_counts  <- roll[date >= CRISIS_START & date <= CRISIS_END,
                       count_above_mp_ex_market]
pre_obs_median <- median(pre_obs_counts)
in_obs_median  <- median(in_obs_counts)
cat(sprintf("Observed pre-crisis median: %.1f (n=%d windows)\n",
            pre_obs_median, length(pre_obs_counts)))
cat(sprintf("Observed in-crisis median:  %.1f (n=%d windows)\n",
            in_obs_median, length(in_obs_counts)))

# ----- Stationary bootstrap helper ------------------------------------------
stationary_indices <- function(n, p) {
  idx <- integer(n)
  idx[1] <- sample.int(n, 1L)
  for (t in 2:n) {
    if (runif(1L) < p) idx[t] <- sample.int(n, 1L)
    else               idx[t] <- if (idx[t - 1L] == n) 1L else idx[t - 1L] + 1L
  }
  idx
}

# ----- Bootstrap loop -------------------------------------------------------
pre_z <- as.matrix(z_cc[date >= BASE_START & date <= BASE_END, tickers, with = FALSE])
n_pre <- nrow(pre_z)
N <- length(tickers)
p_block <- 1 / BLOCK_MEAN
mp_edge_fixed <- (1 + sqrt(N / WINDOW))^2   # N=25 throughout bootstrap (complete cases)

n_windows <- n_pre - WINDOW + 1L

one_rep <- function(b_seed) {
  set.seed(b_seed)
  idx <- stationary_indices(n_pre, p_block)
  z_b <- pre_z[idx, , drop = FALSE]
  counts <- integer(n_windows)
  for (i in seq_len(n_windows)) {
    w <- z_b[i:(i + WINDOW - 1L), , drop = FALSE]
    R <- cor(w)
    evals <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
    counts[i] <- sum(evals[-1L] > mp_edge_fixed)
  }
  counts
}

cat(sprintf("\nRunning %d bootstrap reps in parallel on %d cores (block mean=%d, N=%d, T=%d, windows=%d)...\n",
            B, NCORES, BLOCK_MEAN, N, n_pre, n_windows))
t0 <- proc.time()
seeds <- 20260522L + seq_len(B)
boot_list <- mclapply(seeds, one_rep, mc.cores = NCORES, mc.preschedule = TRUE)
cat(sprintf("Done in %.0fs (parallel speedup ~%.1fx vs single-core est)\n",
            (proc.time() - t0)[3], NCORES * 0.8))

boot_count_matrix <- do.call(rbind, boot_list)
boot_medians <- apply(boot_count_matrix, 1, median)

# ----- T2(b) verdict --------------------------------------------------------
ci <- quantile(boot_medians, c(0.025, 0.5, 0.975))
t2b_pass <- in_obs_median > ci["97.5%"]

cat("\n========== T2(b) verdict ==========\n")
cat(sprintf("Pre-crisis median count: observed %.1f\n", pre_obs_median))
cat(sprintf("Bootstrap distribution of pre-crisis median: median %.2f, 95%% CI [%.2f, %.2f]\n",
            ci["50%"], ci["2.5%"], ci["97.5%"]))
cat(sprintf("In-crisis median count: observed %.1f\n", in_obs_median))
cat(sprintf("T2(b): %s  (need in-crisis median > %.2f, observed %.1f)\n",
            if (t2b_pass) "PASS" else "FAIL", ci["97.5%"], in_obs_median))
cat("====================================\n")

fwrite(data.table(
  pre_obs_median = pre_obs_median,
  in_obs_median = in_obs_median,
  bootstrap_median = ci["50%"],
  ci_lo_95 = ci["2.5%"],
  ci_hi_95 = ci["97.5%"],
  t2b_pass = t2b_pass,
  B = B, block_mean = BLOCK_MEAN
), "output/bootstrap_t2b.csv")

# ----- Per-date bootstrap quantiles for visualisation ------------------------
boot_dates <- z_cc[date >= BASE_START & date <= BASE_END,
                   date][WINDOW:n_pre]   # length n_windows
roll_ci <- data.table(
  date = boot_dates,
  q025 = apply(boot_count_matrix, 2, quantile, 0.025),
  q500 = apply(boot_count_matrix, 2, quantile, 0.50),
  q975 = apply(boot_count_matrix, 2, quantile, 0.975)
)
fwrite(roll_ci, "output/bootstrap_rolling_ci.csv")

# ----- Figure 2b: distribution comparison (pre-crisis vs in-crisis) ---------
# The bootstrap CI on the pre-crisis median is degenerate ([2, 2]) because the
# statistic is an integer count, so a ribbon overlay carries no visual signal.
# A distribution comparison conveys the +1 shift directly.

dist_dt <- rbind(
  data.table(regime = "Pre-crisis (2019–2023)",  count = pre_obs_counts),
  data.table(regime = "In-crisis (2024–2026)",   count = in_obs_counts)
)
dist_dt[, regime := factor(regime,
                           levels = c("Pre-crisis (2019–2023)",
                                      "In-crisis (2024–2026)"))]
prop_dt <- dist_dt[, .(n = .N), by = .(regime, count)]
prop_dt[, prop := n / sum(n), by = regime]
# Ensure every regime × count combination is represented (so the dodge has paired bars)
all_counts <- sort(unique(prop_dt$count))
all_regimes <- levels(prop_dt$regime)
grid <- CJ(regime = all_regimes, count = all_counts)
grid[, regime := factor(regime, levels = all_regimes)]
prop_dt <- merge(grid, prop_dt, by = c("regime", "count"), all.x = TRUE)
prop_dt[is.na(prop), `:=`(n = 0L, prop = 0)]
prop_dt[, count := factor(count, levels = all_counts)]

medians_dt <- dist_dt[, .(median_count = as.numeric(median(count)),
                          n = .N,
                          mean_count = round(mean(count), 2)),
                      by = regime]

p2b <- ggplot(prop_dt, aes(x = count, y = prop, fill = regime)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("Pre-crisis (2019–2023)" = "grey60",
                               "In-crisis (2024–2026)"  = "tomato"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Figure 2b: Distribution of the rolling supra-MP count, pre-crisis vs in-crisis",
    subtitle = sprintf(
      "Pre-crisis (n = %d, median = %g, mean = %.2f).  In-crisis (n = %d, median = %g, mean = %.2f).  Bootstrap 95%% CI on the pre-crisis median: [%g, %g] (B = %d, block mean = %d).",
      medians_dt$n[1], medians_dt$median_count[1], medians_dt$mean_count[1],
      medians_dt$n[2], medians_dt$median_count[2], medians_dt$mean_count[2],
      ci["2.5%"], ci["97.5%"], B, BLOCK_MEAN),
    x = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }   (market mode excluded)")),
    y = "fraction of rolling windows",
    caption = sprintf(
      "Mass shifts from count = 2 (modal pre-crisis, %.0f%% of windows) to count = 3 (modal in-crisis, %.0f%% of windows).",
      100 * prop_dt[regime == "Pre-crisis (2019–2023)" & count == "2", prop],
      100 * prop_dt[regime == "In-crisis (2024–2026)"  & count == "3", prop]
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        plot.subtitle = element_text(size = 9))

ggsave("figures/fig2b_count_with_ci.png", p2b, width = 10, height = 5.2, dpi = 200)
cat("\nWrote figures/fig2b_count_with_ci.png\n")
