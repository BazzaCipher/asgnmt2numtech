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

# ----- Augmented Figure 2 with bootstrap band -------------------------------
plot_dt <- roll[!is.na(count_above_mp_ex_market) & n_assets == 25]
crisis_band <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                          ymin = -Inf, ymax = Inf)

SMOOTH_K <- 63L  # ~one trading quarter; smooths integer {1,2,3} step jaggedness
plot_dt[, count_smooth := frollmean(count_above_mp_ex_market, SMOOTH_K,
                                    align = "right", fill = NA)]
roll_ci[, `:=`(
  q025_smooth = frollmean(q025, SMOOTH_K, align = "right", fill = NA),
  q975_smooth = frollmean(q975, SMOOTH_K, align = "right", fill = NA)
)]

p2b <- ggplot() +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.10) +
  geom_ribbon(data = roll_ci[!is.na(q025_smooth)],
              aes(x = date, ymin = q025_smooth, ymax = q975_smooth),
              fill = "steelblue", alpha = 0.25) +
  geom_step(data = plot_dt,
            aes(x = date, y = count_above_mp_ex_market),
            colour = "grey70", linewidth = 0.3) +
  geom_line(data = plot_dt[!is.na(count_smooth)],
            aes(x = date, y = count_smooth),
            colour = "black", linewidth = 0.7) +
  geom_hline(yintercept = ci["97.5%"], colour = "steelblue",
             linetype = 2, alpha = 0.7) +
  annotate("text",
           x = as.Date("2018-06-01"),
           y = ci["97.5%"] + 0.05,
           label = sprintf("pre-crisis bootstrap 95%% upper bound = %.2f", ci["97.5%"]),
           hjust = 0, colour = "steelblue", size = 3.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Figure 2b: Rolling count with stationary-bootstrap 95% CI",
    subtitle = sprintf("Bootstrap CI on 2019-2023 baseline (N=25, block mean=%d, B=%d reps). Grey = raw integer count; black = %dd trailing mean.",
                       BLOCK_MEAN, B, SMOOTH_K),
    x = NULL, y = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }")),
    caption = sprintf("Blue band = per-date 95%% bootstrap CI; dashed line = bootstrap CI upper bound of the pre-crisis MEDIAN (%.2f). In-crisis median = %.1f.",
                      ci["97.5%"], in_obs_median)
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/fig2b_count_with_ci.png", p2b, width = 10, height = 4.8, dpi = 200)
cat("\nWrote figures/fig2b_count_with_ci.png\n")
