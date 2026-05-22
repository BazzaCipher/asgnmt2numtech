# Days 5-6: MP eigenvalue analysis + rolling supra-MP count (Figures 1 & 2)
#
# Inputs:
#   data/clean/std_residuals_v1.csv   (winsorised GJR-GARCH-t std residuals)
#
# Outputs:
#   figures/fig1_mp_spectrum.png            (full-sample spectrum vs MP density)
#   figures/fig2_rolling_count.png          (rolling count + Bai-Perron breaks)
#   output/mp_spectrum_fullsample.csv       (eigenvalues, eigenvectors of R̄)
#   output/rolling_count_v1.csv             (date, count_above_mp_ex_market, eigvals)
#   output/bai_perron_breaks.csv            (break dates from strucchange)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(strucchange)
})

RESID_IN <- "data/clean/std_residuals_v1.csv"
WINDOW   <- 252L
CRISIS_START <- as.Date("2024-01-01")
CRISIS_END   <- as.Date("2026-05-22")
MIN_WINDOW_COVERAGE <- 0.9

z <- fread(RESID_IN)
z[, date := as.Date(date)]
setorder(z, date)
tickers <- setdiff(names(z), "date")
cat(sprintf("Loaded std residuals: %d dates x %d series\n", nrow(z), length(tickers)))

# ----- Figure 1: Full-sample spectrum vs MP density -------------------------
# Restrict to the window where ALL series have data (post-Ghana entry).
z_complete <- z[complete.cases(z[, tickers, with = FALSE])]
T_full <- nrow(z_complete)
N_full <- length(tickers)
q_full <- N_full / T_full
cat(sprintf("Full-sample (complete-cases) spectrum: T=%d, N=%d, Q=T/N=%.1f\n",
            T_full, N_full, 1/q_full))

R_full <- cor(as.matrix(z_complete[, tickers, with = FALSE]))
eig <- eigen(R_full, symmetric = TRUE)
eigvals <- sort(eig$values, decreasing = TRUE)
mp_lo <- (1 - sqrt(q_full))^2
mp_hi <- (1 + sqrt(q_full))^2

# MP density
mp_density <- function(x, q) {
  out <- numeric(length(x))
  lo <- (1 - sqrt(q))^2; hi <- (1 + sqrt(q))^2
  in_range <- x > lo & x < hi
  out[in_range] <- sqrt((hi - x[in_range]) * (x[in_range] - lo)) / (2 * pi * q * x[in_range])
  out
}
x_grid <- seq(max(mp_lo * 0.95, 0.01), mp_hi * 1.05, length.out = 400)
mp_curve <- data.table(x = x_grid, density = mp_density(x_grid, q_full))

n_above <- sum(eigvals > mp_hi)
cat(sprintf("MP upper edge λ+ = %.3f; eigenvalues above edge: %d (incl. market mode)\n",
            mp_hi, n_above))
cat("Top 8 eigenvalues:\n"); print(round(eigvals[1:8], 3))

eig_dt <- data.table(rank = seq_along(eigvals), eigenvalue = eigvals,
                     above_mp = eigvals > mp_hi)
fwrite(eig_dt, "output/mp_spectrum_fullsample.csv")
fwrite(data.table(ticker = tickers, eig$vectors),
       "output/mp_eigenvectors_fullsample.csv")

p1 <- ggplot() +
  geom_histogram(data = data.table(eigenvalue = eigvals),
                 aes(x = eigenvalue, y = after_stat(density)),
                 binwidth = 0.08, fill = "grey80", colour = "grey40") +
  geom_line(data = mp_curve, aes(x = x, y = density),
            colour = "steelblue", linewidth = 0.9) +
  geom_vline(xintercept = c(mp_lo, mp_hi), linetype = 2, colour = "red", alpha = 0.7) +
  geom_vline(xintercept = eigvals[eigvals > mp_hi],
             linetype = 3, colour = "darkgreen", alpha = 0.8) +
  annotate("text", x = mp_hi, y = max(mp_curve$density) * 0.9,
           label = sprintf(" λ+ = %.3f", mp_hi), hjust = 0, colour = "red", size = 3.5) +
  annotate("text", x = max(eigvals) * 0.92, y = 0.05,
           label = sprintf("%d eigenvalues > λ+", n_above),
           hjust = 1, colour = "darkgreen", size = 3.5) +
  coord_cartesian(xlim = c(0, max(eigvals) * 1.05)) +
  labs(
    title = "Figure 1: Eigenvalue spectrum of full-sample correlation matrix",
    subtitle = sprintf("Standardised GJR-GARCH-t residuals, %s to %s (N=%d, T=%d, Q=T/N=%.1f)",
                       min(z_complete$date), max(z_complete$date), N_full, T_full, 1/q_full),
    x = "eigenvalue", y = "density",
    caption = "Blue curve = Marchenko-Pastur density; red dashes = MP edges; green dashes = supra-edge eigenvalues"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/fig1_mp_spectrum.png", p1, width = 9, height = 4.8, dpi = 150)
cat("Wrote figures/fig1_mp_spectrum.png\n")

# ----- Rolling 252-day count above MP edge ----------------------------------
mp_upper <- function(N, T) (1 + sqrt(N / T))^2

z_mat <- as.matrix(z[, tickers, with = FALSE])
n_total <- nrow(z_mat)
window_starts <- seq_len(n_total - WINDOW + 1L)

out <- data.table(
  date = z$date[window_starts + WINDOW - 1L],
  n_assets = NA_integer_,
  count_above_mp_raw = NA_integer_,
  count_above_mp_ex_market = NA_integer_,
  top_eig = NA_real_, second_eig = NA_real_, third_eig = NA_real_,
  fourth_eig = NA_real_, mp_edge = NA_real_
)

for (i in seq_along(window_starts)) {
  w <- z_mat[window_starts[i]:(window_starts[i] + WINDOW - 1L), , drop = FALSE]
  keep <- apply(w, 2L, function(x)
    sum(!is.na(x)) >= MIN_WINDOW_COVERAGE * WINDOW && sd(x, na.rm = TRUE) > 0)
  w <- w[, keep, drop = FALSE]
  for (j in seq_len(ncol(w))) {
    na_idx <- is.na(w[, j])
    if (any(na_idx)) w[na_idx, j] <- mean(w[, j], na.rm = TRUE)
  }
  Nw <- ncol(w); if (Nw < 5L) next
  R <- cor(w)
  evals <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  edge <- mp_upper(Nw, WINDOW)
  out[i, `:=`(
    n_assets = Nw,
    count_above_mp_raw       = sum(evals > edge),
    count_above_mp_ex_market = sum(evals[-1L] > edge),
    top_eig = evals[1L], second_eig = evals[2L],
    third_eig = evals[3L], fourth_eig = evals[4L],
    mp_edge = edge
  )]
}

fwrite(out, "output/rolling_count_v1.csv")
cat(sprintf("Wrote output/rolling_count_v1.csv (%d windows)\n", nrow(out)))

# ----- Bai-Perron breaks on the count series --------------------------------
plot_dt <- out[!is.na(count_above_mp_ex_market)]
bp_input <- plot_dt[, .(date, count = count_above_mp_ex_market)]
bp_full <- tryCatch(
  breakpoints(count ~ 1, data = bp_input, h = 0.05),
  error = function(e) NULL
)
bp_bic <- if (!is.null(bp_full)) which.min(BIC(bp_full)) - 1L else 0L
if (!is.null(bp_full) && bp_bic > 0L) {
  bp <- breakpoints(bp_full, breaks = bp_bic)
  break_dates <- bp_input$date[bp$breakpoints]
} else {
  break_dates <- as.Date(character(0))
}
cat(sprintf("Bai-Perron BIC selects %d breaks: %s\n",
            bp_bic,
            if (length(break_dates)) paste(break_dates, collapse = ", ") else "none"))
# Also save BIC table for diagnostics
if (!is.null(bp_full)) {
  bic_vec <- BIC(bp_full)
  fwrite(data.table(n_breaks = seq_along(bic_vec) - 1L,
                    BIC = as.numeric(bic_vec)),
         "output/bai_perron_bic_table.csv")
}
fwrite(data.table(break_date = break_dates), "output/bai_perron_breaks.csv")

# ----- Figure 2: rolling count with crisis band + breaks --------------------
crisis_band <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                          ymin = -Inf, ymax = Inf)

SMOOTH_K <- 63L  # ~one trading quarter; smooths integer {1,2,3} step jaggedness
plot_dt[, count_smooth := frollmean(count_above_mp_ex_market, SMOOTH_K,
                                    align = "right", fill = NA)]

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
                      length(tickers))
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/fig2_rolling_count.png", p2, width = 10, height = 4.8, dpi = 200)
cat("Wrote figures/fig2_rolling_count.png\n")

# ----- T1 check (pre-crisis baseline stability, 2019-2023, N=25) ------------
pre <- plot_dt[date >= as.Date("2019-01-01") & date <= as.Date("2023-12-31") &
               n_assets == 25, count_above_mp_ex_market]
in_  <- plot_dt[date >= CRISIS_START & date <= CRISIS_END,
                count_above_mp_ex_market]

q <- function(x, p) quantile(x, p, na.rm = TRUE)

cat("\n========== T1 / T2(eye-test) ON CLEANED PIPELINE ==========\n")
cat(sprintf("Pre-crisis (2019-2023, N=25):  median %.1f, IQR [%.1f, %.1f], range [%d, %d], n=%d\n",
            median(pre), q(pre,.25), q(pre,.75), min(pre), max(pre), length(pre)))
cat(sprintf("In-crisis  (2024-2026):        median %.1f, IQR [%.1f, %.1f], range [%d, %d], n=%d\n",
            median(in_), q(in_,.25), q(in_,.75), min(in_), max(in_), length(in_)))
cat(sprintf("T1 (IQR width <= 1 and contained in {1,2,3}):  %s\n",
            if (max(pre) - min(pre) <= 2 && all(pre %in% 1:3)) "PASS" else "FAIL"))
cat(sprintf("Diff (in − pre): median %+.1f\n", median(in_) - median(pre)))
cat("===========================================================\n")
