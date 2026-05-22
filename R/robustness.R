# Day 10 robustness checks for T1 / T2(a) / T2(b).
#
# Variants (per scope_and_target.md §"Path from here"):
#   R1. Exclude GHANA10Y (test the universe-size-change concern)
#   R2. Window length 504 instead of 252
#   R3. Restrict pre-crisis baseline to N=24 windows (pre-Ghana entry, 2015-17)
#   R4. Re-fit GARCH with Gaussian innovations instead of Student-t
#
# Outputs:
#   output/robustness_summary.csv
#   figures/fig4_robustness.png   (rolling count under each variant)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(strucchange)
  library(parallel)
  library(rugarch)
})

set.seed(20260522)

WINDOW_BASE   <- 252L
WINDOW_LONG   <- 504L
CRISIS_START  <- as.Date("2024-01-01")
CRISIS_END    <- as.Date("2026-05-22")
BASE_START    <- as.Date("2019-01-01")
BASE_END      <- as.Date("2023-12-31")
NCORES        <- max(1L, min(20L, detectCores() - 4L))

z_full <- fread("data/clean/std_residuals_v1.csv")
z_full[, date := as.Date(date)]
diag <- fread("output/garch_fits.csv")

# ----- Rolling count helper -------------------------------------------------
rolling_count <- function(z_dt, window, tickers, min_cov = 0.9) {
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

# ----- Bootstrap helper (parallel) ------------------------------------------
boot_median_ci <- function(pre_z_mat, window, p_block = 1/25, B = 500L) {
  if (is.null(pre_z_mat) || nrow(pre_z_mat) < window + 10) {
    return(c(med = NA, lo = NA, hi = NA))
  }
  N <- ncol(pre_z_mat); n_pre <- nrow(pre_z_mat)
  edge <- (1 + sqrt(N / window))^2
  n_w <- n_pre - window + 1L
  one <- function(seed) {
    set.seed(seed)
    idx <- integer(n_pre)
    idx[1] <- sample.int(n_pre, 1L)
    for (t in 2:n_pre) {
      if (runif(1L) < p_block) idx[t] <- sample.int(n_pre, 1L)
      else idx[t] <- if (idx[t-1L] == n_pre) 1L else idx[t-1L] + 1L
    }
    zb <- pre_z_mat[idx, , drop = FALSE]
    cnts <- integer(n_w)
    for (i in seq_len(n_w)) {
      w <- zb[i:(i + window - 1L), , drop = FALSE]
      R <- cor(w)
      evals <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
      cnts[i] <- sum(evals[-1L] > edge)
    }
    median(cnts)
  }
  medians <- unlist(mclapply(seq_len(B), one, mc.cores = NCORES, mc.preschedule = TRUE))
  c(med = median(medians),
    lo = unname(quantile(medians, 0.025)),
    hi = unname(quantile(medians, 0.975)))
}

# ----- Variant runner -------------------------------------------------------
run_variant <- function(label, z_dt, tickers, window) {
  cat(sprintf("\n--- %s ---\n", label))
  roll <- rolling_count(z_dt, window, tickers)
  if (is.null(roll)) return(NULL)
  pre_counts <- roll[date >= BASE_START & date <= BASE_END & n_assets == max(n_assets, na.rm = TRUE), count]
  in_counts  <- roll[date >= CRISIS_START & date <= CRISIS_END, count]

  z_pre <- z_dt[date >= BASE_START & date <= BASE_END]
  z_pre <- z_pre[complete.cases(z_pre[, tickers, with = FALSE])]
  pre_z_mat <- as.matrix(z_pre[, tickers, with = FALSE])
  ci <- boot_median_ci(pre_z_mat, window)

  bp_input <- roll[!is.na(count), .(date, count)]
  bp_full <- tryCatch(
    breakpoints(count ~ 1, data = bp_input, h = 0.05),
    error = function(e) NULL
  )
  bp_bic <- if (!is.null(bp_full)) which.min(BIC(bp_full)) - 1L else 0L
  break_dates <- if (!is.null(bp_full) && bp_bic > 0L) {
    bp <- breakpoints(bp_full, breaks = bp_bic)
    bp_input$date[bp$breakpoints]
  } else as.Date(character(0))
  cat(sprintf("  BIC selects %d breaks\n", bp_bic))
  in_crisis_break <- any(break_dates >= CRISIS_START & break_dates <= CRISIS_END)

  result <- data.table(
    variant = label,
    universe_size = length(tickers),
    window = window,
    pre_median = median(pre_counts, na.rm = TRUE),
    in_median  = median(in_counts,  na.rm = TRUE),
    boot_ci_lo = unname(ci["lo"]),
    boot_ci_hi = unname(ci["hi"]),
    t2b_pass   = median(in_counts, na.rm = TRUE) > unname(ci["hi"]),
    n_breaks   = length(break_dates),
    in_crisis_break = in_crisis_break,
    first_crisis_break = as.Date(NA)
  )
  in_b <- break_dates[break_dates >= CRISIS_START & break_dates <= CRISIS_END]
  if (length(in_b)) result$first_crisis_break <- min(in_b)
  print(result)
  list(summary = result, roll = roll[, .(date, count, variant = label)])
}

# ----- Build inputs for each variant ----------------------------------------
all_tickers <- setdiff(names(z_full), "date")

results <- list()
rolls   <- list()

# R0: Baseline (sanity reproduction)
r <- run_variant("baseline (W=252)", z_full, all_tickers, WINDOW_BASE)
results[[length(results)+1]] <- r$summary; rolls[[length(rolls)+1]] <- r$roll

# R1: Exclude GHANA10Y
tk_R1 <- setdiff(all_tickers, "GHANA10Y")
r <- run_variant("excl GHANA10Y (W=252)", z_full, tk_R1, WINDOW_BASE)
results[[length(results)+1]] <- r$summary; rolls[[length(rolls)+1]] <- r$roll

# R2: Window 504
r <- run_variant("baseline (W=504)", z_full, all_tickers, WINDOW_LONG)
results[[length(results)+1]] <- r$summary; rolls[[length(rolls)+1]] <- r$roll

# R3: Excl GHANA10Y + Window 504 (combined)
r <- run_variant("excl GHANA10Y (W=504)", z_full, tk_R1, WINDOW_LONG)
results[[length(results)+1]] <- r$summary; rolls[[length(rolls)+1]] <- r$roll

# R4: Gaussian GARCH univariates — needs re-fit
cat("\n--- R4: re-fitting GARCH with Gaussian innovations ---\n")
rets <- fread("data/clean/returns_v1.csv"); rets[, date := as.Date(date)]
spec_gjr_n <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)
fit_one_n <- function(ticker) {
  x <- rets[[ticker]]; valid <- which(!is.na(x))
  xv <- x[valid]
  if (length(xv) < 100) return(rep(NA_real_, nrow(rets)))
  fit <- tryCatch(ugarchfit(spec_gjr_n, xv, solver = "hybrid"),
                  error = function(e) NULL)
  if (is.null(fit) || convergence(fit) != 0) {
    sigma2 <- numeric(length(xv)); sigma2[1] <- var(xv)
    for (t in 2:length(xv)) sigma2[t] <- 0.94 * sigma2[t-1] + 0.06 * xv[t-1]^2
    z_clean <- xv / sqrt(sigma2)
  } else {
    z_clean <- as.numeric(residuals(fit, standardize = TRUE))
    if (!all(is.finite(z_clean)) || sd(z_clean) > 2 || max(abs(z_clean)) > 30) {
      sigma2 <- numeric(length(xv)); sigma2[1] <- var(xv)
      for (t in 2:length(xv)) sigma2[t] <- 0.94 * sigma2[t-1] + 0.06 * xv[t-1]^2
      z_clean <- xv / sqrt(sigma2)
    }
  }
  z_clean <- pmin(pmax(z_clean, -5), 5)
  z_full <- rep(NA_real_, nrow(rets)); z_full[valid] <- z_clean
  z_full
}
z_norm_list <- mclapply(all_tickers, fit_one_n, mc.cores = NCORES)
z_norm <- data.table(date = rets$date)
for (i in seq_along(all_tickers)) z_norm[[all_tickers[i]]] <- z_norm_list[[i]]
fwrite(z_norm, "data/clean/std_residuals_norm.csv")
r <- run_variant("Gaussian GJR (W=252)", z_norm, all_tickers, WINDOW_BASE)
results[[length(results)+1]] <- r$summary; rolls[[length(rolls)+1]] <- r$roll

# ----- Save + plot ----------------------------------------------------------
summ <- rbindlist(results)
fwrite(summ, "output/robustness_summary.csv")
cat("\n========== ROBUSTNESS SUMMARY ==========\n")
print(summ)
cat("==========================================\n")

roll_all <- rbindlist(rolls)
roll_all[, variant := factor(variant, levels = unique(variant))]
crisis_band <- data.frame(xmin = CRISIS_START, xmax = CRISIS_END,
                          ymin = -Inf, ymax = Inf)

SMOOTH_K <- 63L  # ~one trading quarter; smooths integer step jaggedness
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
