# Day-1 smoke test: rolling MP-edge eigenvalue count, raw returns (no GARCH layer).
#
# Goal: does the 2024-26 cocoa crisis produce a visible increase in the count of
# correlation eigenvalues above the Marchenko-Pastur upper edge?
#
# Method, per proposal §"Rolling eigenvalue count":
#   1. 252-day rolling window of returns.
#       - log returns for "level" series (futures, equities, FX, indices)
#       - first-differences for "yield" series (US10Y, Ghana10Y, both in %)
#   2. Sample correlation matrix R_t in each window.
#   3. Eigendecompose, remove the largest eigenvalue (market mode).
#   4. Count remaining eigenvalues above the MP upper edge λ+ = (1 + sqrt(N/T))^2.
#   5. Plot the count over time, shade the crisis window.

suppressPackageStartupMessages({
  library(data.table)
  library(zoo)
  library(ggplot2)
  library(scales)
})

WINDOW         <- 252L
CRISIS_START   <- as.Date("2024-01-01")
CRISIS_END     <- as.Date("2026-05-22")
PRE_CRISIS_END <- as.Date("2023-12-31")
MIN_WINDOW_COVERAGE <- 0.9   # drop a series from a window if <90% non-NA

# ----- Load -----------------------------------------------------------------
panel <- fread("data/clean/panel.csv")
panel[, date := as.Date(date)]
setorder(panel, date)

diag <- fread("data/clean/panel_diagnostics.csv")
yield_tickers <- diag[value_kind == "yield", ticker]
level_tickers <- diag[value_kind == "level", ticker]
all_tickers   <- c(level_tickers, yield_tickers)
all_tickers   <- intersect(all_tickers, names(panel))

cat(sprintf("Loaded panel: %d dates x %d series (%d level, %d yield)\n",
            nrow(panel), length(all_tickers),
            length(level_tickers), length(yield_tickers)))

# ----- Returns --------------------------------------------------------------
# Forward-fill small gaps (max 3 days) inside each series, then take returns.
filled <- copy(panel)
for (col in all_tickers) {
  filled[, (col) := na.locf(get(col), na.rm = FALSE, maxgap = 3)]
}

# Drop rows where the cocoa NY anchor is missing (collapses to its trading calendar)
filled <- filled[!is.na(COCOA_NY)]

rets <- copy(filled)
for (col in all_tickers) {
  x <- filled[[col]]
  if (col %in% yield_tickers) {
    r <- c(NA_real_, diff(x))          # first-difference of % yield
  } else {
    r <- c(NA_real_, diff(log(x)))     # log return
  }
  rets[, (col) := r]
}
rets <- rets[-1L]   # drop initial NA row

cat(sprintf("Returns panel: %d rows from %s to %s\n",
            nrow(rets), min(rets$date), max(rets$date)))

# ----- MP edge --------------------------------------------------------------
mp_upper <- function(N, T) (1 + sqrt(N / T))^2

# ----- Rolling eigenvalue count ---------------------------------------------
ret_mat <- as.matrix(rets[, all_tickers, with = FALSE])
n_total <- nrow(ret_mat)
stopifnot(n_total > WINDOW + 10)

window_starts <- seq_len(n_total - WINDOW + 1L)
out <- data.table(
  date = rets$date[window_starts + WINDOW - 1L],
  n_assets = NA_integer_,
  count_above_mp_raw = NA_integer_,         # including market mode
  count_above_mp_ex_market = NA_integer_,   # market mode removed
  top_eig = NA_real_,
  second_eig = NA_real_,
  mp_edge = NA_real_
)

for (i in seq_along(window_starts)) {
  w <- ret_mat[window_starts[i]:(window_starts[i] + WINDOW - 1L), , drop = FALSE]
  keep <- apply(w, 2L, function(x) {
    sum(!is.na(x)) >= MIN_WINDOW_COVERAGE * WINDOW && sd(x, na.rm = TRUE) > 0
  })
  w <- w[, keep, drop = FALSE]
  # Mean-impute the few remaining NAs in window
  for (j in seq_len(ncol(w))) {
    na_idx <- is.na(w[, j])
    if (any(na_idx)) w[na_idx, j] <- mean(w[, j], na.rm = TRUE)
  }
  Nw <- ncol(w)
  if (Nw < 5L) next
  R <- cor(w)
  evals <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  edge <- mp_upper(Nw, WINDOW)
  out[i, `:=`(
    n_assets = Nw,
    count_above_mp_raw = sum(evals > edge),
    count_above_mp_ex_market = sum(evals[-1L] > edge),
    top_eig = evals[1L],
    second_eig = evals[2L],
    mp_edge = edge
  )]
}

dir.create("output", showWarnings = FALSE)
fwrite(out, "output/day1_rolling_eig_count.csv")
cat(sprintf("Wrote output/day1_rolling_eig_count.csv (%d windows)\n", nrow(out)))

# ----- Plot -----------------------------------------------------------------
plot_dt <- out[!is.na(count_above_mp_ex_market)]
crisis_band <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                          ymin = -Inf, ymax = Inf)

p <- ggplot(plot_dt, aes(x = date)) +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.12) +
  geom_step(aes(y = count_above_mp_ex_market), colour = "black", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Day-1 smoke test: rolling 252-day eigenvalue count above MP edge",
    subtitle = "Largest eigenvalue (market mode) excluded. Pink band = 2024-26 cocoa crisis.",
    x = NULL,
    y = expression(paste("#{ ", lambda[k], " > ", lambda["+"], " }")),
    caption = sprintf("Universe = %d series (incl. Ghana 10Y from Apr-2017); raw returns; MP edge λ+ = (1+√(N/T))²",
                      length(all_tickers))
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

dir.create("figures", showWarnings = FALSE)
ggsave("figures/day1_rolling_eig_count.png", p, width = 10, height = 4.8, dpi = 150)
cat("Wrote figures/day1_rolling_eig_count.png\n")

# Also plot the top-2 eigenvalues themselves
p2 <- ggplot(plot_dt, aes(x = date)) +
  geom_rect(data = crisis_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "tomato", alpha = 0.12) +
  geom_line(aes(y = top_eig,    colour = "1st (market mode)"), linewidth = 0.6) +
  geom_line(aes(y = second_eig, colour = "2nd"), linewidth = 0.6) +
  geom_line(aes(y = mp_edge,    colour = "MP edge λ+"), linewidth = 0.4, linetype = 2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_colour_manual(values = c("1st (market mode)" = "steelblue",
                                 "2nd" = "darkorange",
                                 "MP edge λ+" = "grey40")) +
  labs(title = "Top eigenvalues of the rolling correlation matrix",
       y = "eigenvalue", x = NULL, colour = NULL) +
  theme_minimal(base_size = 11)
ggsave("figures/day1_top_eigenvalues.png", p2, width = 10, height = 4.5, dpi = 150)
cat("Wrote figures/day1_top_eigenvalues.png\n")

# ----- Verdict --------------------------------------------------------------
pre_crisis <- plot_dt[date <= PRE_CRISIS_END, count_above_mp_ex_market]
in_crisis  <- plot_dt[date >= CRISIS_START & date <= CRISIS_END, count_above_mp_ex_market]

q <- function(x, p) quantile(x, p, na.rm = TRUE)

cat("\n========== VERDICT ==========\n")
cat(sprintf("Window:        252 days\n"))
cat(sprintf("Universe:      %d series\n", length(all_tickers)))
cat(sprintf("Pre-crisis count_ex_market (2015-2023):  median %.1f, IQR [%.1f, %.1f], n=%d\n",
            median(pre_crisis), q(pre_crisis, .25), q(pre_crisis, .75), length(pre_crisis)))
cat(sprintf("In-crisis  count_ex_market (2024-2026):  median %.1f, IQR [%.1f, %.1f], n=%d\n",
            median(in_crisis), q(in_crisis, .25), q(in_crisis, .75), length(in_crisis)))
peak_idx <- which.max(plot_dt$count_above_mp_ex_market)
cat(sprintf("Peak count:    %d on %s\n",
            plot_dt$count_above_mp_ex_market[peak_idx],
            plot_dt$date[peak_idx]))
cat(sprintf("Diff (in-pre): median %+.1f\n",
            median(in_crisis) - median(pre_crisis)))
cat("=============================\n")
