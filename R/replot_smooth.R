# Regenerate Figures 2, 2b, 4 with a 63-day trailing-mean overlay to smooth
# the integer-valued rolling MP-edge count ({1,2,3}). Reads cached CSVs only;
# no GARCH refit or bootstrap. Matches the plot blocks in
# mp_analysis.R / bootstrap_count.R / robustness.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

SMOOTH_K     <- 63L
CRISIS_START <- as.Date("2024-01-01")
CRISIS_END   <- as.Date("2026-05-22")
crisis_band  <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                           ymin = -Inf, ymax = Inf)

# ===== Figure 2 ==============================================================
roll <- fread("output/rolling_count_v1.csv")
roll[, date := as.Date(date)]
plot_dt <- roll[!is.na(count_above_mp_ex_market)]
plot_dt[, count_smooth := frollmean(count_above_mp_ex_market, SMOOTH_K,
                                    align = "right", fill = NA)]

bp_dt <- tryCatch(fread("output/bai_perron_breaks.csv"), error = function(e) NULL)
break_dates <- if (!is.null(bp_dt) && nrow(bp_dt)) as.Date(bp_dt$break_date) else as.Date(character(0))

n_tickers <- length(setdiff(names(fread("data/clean/std_residuals_v1.csv", nrows = 1)), "date"))

p2 <- ggplot(plot_dt, aes(x = date)) +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.12) +
  geom_step(aes(y = count_above_mp_ex_market), colour = "grey70", linewidth = 0.3) +
  geom_line(aes(y = count_smooth), colour = "black", linewidth = 0.7, na.rm = TRUE)
if (length(break_dates) > 0) {
  p2 <- p2 + geom_vline(xintercept = as.numeric(break_dates),
                        colour = "purple", linetype = 2, alpha = 0.7)
}
p2 <- p2 +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Figure 2: Rolling 252-day count of correlation eigenvalues above MP edge",
    subtitle = sprintf("Standardised GJR-GARCH-t residuals; market mode excluded. Grey = raw integer count; black = %dd trailing mean. Pink = crisis; purple = Bai-Perron breaks.", SMOOTH_K),
    x = NULL, y = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }")),
    caption = sprintf("Universe = %d series. MP edge λ+ = (1+√(N/T))² recomputed per window. Pre-Ghana-entry windows have N=24.",
                      n_tickers)
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/fig2_rolling_count.png", p2, width = 10, height = 4.8, dpi = 200)
cat("Wrote figures/fig2_rolling_count.png\n")

# ===== Figure 2b =============================================================
roll_ci <- fread("output/bootstrap_rolling_ci.csv")
roll_ci[, date := as.Date(date)]
roll_ci[, `:=`(
  q025_smooth = frollmean(q025, SMOOTH_K, align = "right", fill = NA),
  q975_smooth = frollmean(q975, SMOOTH_K, align = "right", fill = NA)
)]

t2b <- fread("output/bootstrap_t2b.csv")
ci_hi <- t2b$ci_hi_95[1L]
in_obs_median <- t2b$in_obs_median[1L]
B          <- t2b$B[1L]
BLOCK_MEAN <- t2b$block_mean[1L]

plot_dt_25 <- roll[!is.na(count_above_mp_ex_market) & n_assets == 25]
plot_dt_25[, count_smooth := frollmean(count_above_mp_ex_market, SMOOTH_K,
                                       align = "right", fill = NA)]

p2b <- ggplot() +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.10) +
  geom_ribbon(data = roll_ci[!is.na(q025_smooth)],
              aes(x = date, ymin = q025_smooth, ymax = q975_smooth),
              fill = "steelblue", alpha = 0.25) +
  geom_step(data = plot_dt_25,
            aes(x = date, y = count_above_mp_ex_market),
            colour = "grey70", linewidth = 0.3) +
  geom_line(data = plot_dt_25[!is.na(count_smooth)],
            aes(x = date, y = count_smooth),
            colour = "black", linewidth = 0.7) +
  geom_hline(yintercept = ci_hi, colour = "steelblue",
             linetype = 2, alpha = 0.7) +
  annotate("text",
           x = as.Date("2018-06-01"),
           y = ci_hi + 0.05,
           label = sprintf("pre-crisis bootstrap 95%% upper bound = %.2f", ci_hi),
           hjust = 0, colour = "steelblue", size = 3.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Figure 2b: Rolling count with stationary-bootstrap 95% CI",
    subtitle = sprintf("Bootstrap CI on 2019-2023 baseline (N=25, block mean=%d, B=%d reps). Grey = raw integer count; black = %dd trailing mean.",
                       BLOCK_MEAN, B, SMOOTH_K),
    x = NULL, y = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }")),
    caption = sprintf("Blue band = per-date 95%% bootstrap CI; dashed line = bootstrap CI upper bound of the pre-crisis MEDIAN (%.2f). In-crisis median = %.1f.",
                      ci_hi, in_obs_median)
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/fig2b_count_with_ci.png", p2b, width = 10, height = 4.8, dpi = 200)
cat("Wrote figures/fig2b_count_with_ci.png\n")

# ===== Figure 4 (recompute per-variant rolling counts) =======================
WINDOW_BASE <- 252L
WINDOW_LONG <- 504L

rolling_count_local <- function(z_dt, window, tickers, min_cov = 0.9) {
  z_mat <- as.matrix(z_dt[, tickers, with = FALSE])
  n_total <- nrow(z_mat)
  if (n_total < window + 10) return(NULL)
  starts <- seq_len(n_total - window + 1L)
  out <- data.table(
    date = z_dt$date[starts + window - 1L],
    n_assets = NA_integer_,
    count = NA_integer_
  )
  for (i in seq_along(starts)) {
    w <- z_mat[starts[i]:(starts[i] + window - 1L), , drop = FALSE]
    keep <- apply(w, 2, function(x)
      sum(!is.na(x)) >= min_cov * window && sd(x, na.rm = TRUE) > 0)
    w <- w[, keep, drop = FALSE]
    for (j in seq_len(ncol(w))) {
      na_idx <- is.na(w[, j])
      if (any(na_idx)) w[na_idx, j] <- mean(w[, j], na.rm = TRUE)
    }
    Nw <- ncol(w); if (Nw < 5L) next
    R <- cor(w)
    evals <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
    edge <- (1 + sqrt(Nw / window))^2
    out[i, `:=`(n_assets = Nw, count = sum(evals[-1L] > edge))]
  }
  out
}

z_full <- fread("data/clean/std_residuals_v1.csv"); z_full[, date := as.Date(date)]
z_norm <- fread("data/clean/std_residuals_norm.csv"); z_norm[, date := as.Date(date)]
all_tickers <- setdiff(names(z_full), "date")
tk_R1       <- setdiff(all_tickers, "GHANA10Y")

variants <- list(
  list(label = "baseline (W=252)",        z = z_full, tk = all_tickers, w = WINDOW_BASE),
  list(label = "excl GHANA10Y (W=252)",   z = z_full, tk = tk_R1,       w = WINDOW_BASE),
  list(label = "baseline (W=504)",        z = z_full, tk = all_tickers, w = WINDOW_LONG),
  list(label = "excl GHANA10Y (W=504)",   z = z_full, tk = tk_R1,       w = WINDOW_LONG),
  list(label = "Gaussian GJR (W=252)",    z = z_norm, tk = all_tickers, w = WINDOW_BASE)
)
rolls <- list()
for (v in variants) {
  cat(sprintf("  recomputing rolling count: %s\n", v$label))
  r <- rolling_count_local(v$z, v$w, v$tk)
  if (!is.null(r)) rolls[[length(rolls)+1]] <- r[, .(date, count, variant = v$label)]
}
roll_all <- rbindlist(rolls)
roll_all[, variant := factor(variant, levels = sapply(variants, `[[`, "label"))]
roll_all[, count_smooth := frollmean(count, SMOOTH_K, align = "right", fill = NA),
         by = variant]

p4 <- ggplot(roll_all[!is.na(count)], aes(x = date, y = count)) +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.10) +
  geom_step(colour = "grey70", linewidth = 0.3) +
  geom_line(aes(y = count_smooth), colour = "black",
            linewidth = 0.6, na.rm = TRUE) +
  facet_wrap(~ variant, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Figure 4: Rolling MP-edge count under robustness variants",
    subtitle = sprintf("Grey = raw integer count; black = %dd trailing mean. If the +1 crisis signal survives every variant, the result is robust.", SMOOTH_K),
    x = NULL, y = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }"))
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9))
ggsave("figures/fig4_robustness.png", p4, width = 11, height = 11, dpi = 200)
cat("Wrote figures/fig4_robustness.png\n")
