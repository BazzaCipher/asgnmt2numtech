# Apply locked cleaning decisions from docs/scope_and_target.md to the parsed
# panel and write the frozen v1 outputs:
#
#   data/clean/panel_v1.csv      — wide levels, restricted to cocoa-NY calendar
#   data/clean/returns_v1.csv    — wide returns (log for levels, first-diff for yields)
#   data/clean/panel_v1_summary.csv — per-series coverage report on the frozen panel
#
# Locked cleaning rules:
#   • WTI returns on 2020-04-20 and 2020-04-21 = NA (front-month settled at -$37.63
#     on 2020-04-20; log return undefined).
#   • Trading-day calendar = dates where COCOA_NY is observed.
#   • Forward-fill within each series, max 3 days.
#   • Returns: log diff for level series; first-difference (in %-points) for yield series.

suppressPackageStartupMessages({
  library(data.table)
  library(zoo)
})

PANEL_IN   <- "data/clean/panel.csv"
DIAG_IN    <- "data/clean/panel_diagnostics.csv"
PANEL_OUT  <- "data/clean/panel_v1.csv"
RETS_OUT   <- "data/clean/returns_v1.csv"
SUMMARY_OUT <- "data/clean/panel_v1_summary.csv"

panel <- fread(PANEL_IN)
panel[, date := as.Date(date)]
setorder(panel, date)

diag <- fread(DIAG_IN)
yield_tickers <- diag[value_kind == "yield", ticker]
all_tickers   <- diag$ticker

# 1. Anchor on COCOA_NY trading days
panel <- panel[!is.na(COCOA_NY)]

# 2. Forward-fill small gaps inside each series (max 3 days)
for (col in all_tickers) {
  panel[, (col) := na.locf(get(col), na.rm = FALSE, maxgap = 3)]
}

# 3. WTI negative-price patch BEFORE computing returns
#    (we want log(WTI_t) - log(WTI_{t-1}) to be NA on the two affected days)
wti_bad_dates <- as.Date(c("2020-04-20", "2020-04-21"))

# 4. Compute returns
rets <- data.table(date = panel$date)
for (col in all_tickers) {
  x <- panel[[col]]
  if (col %in% yield_tickers) {
    r <- c(NA_real_, diff(x))           # first-difference of % yield
  } else {
    r <- c(NA_real_, diff(suppressWarnings(log(x))))  # log return (WTI may have negatives)
    if (col == "WTI") {
      r[panel$date %in% wti_bad_dates] <- NA_real_
    }
  }
  rets[[col]] <- r
}
rets <- rets[-1L]   # drop the initial NA row

# 5. Write frozen outputs
fwrite(panel, PANEL_OUT)
fwrite(rets,  RETS_OUT)

# 6. Per-series coverage summary on the frozen panel
summary_dt <- data.table(
  ticker = all_tickers,
  value_kind = sapply(all_tickers, function(t) diag[ticker == t, value_kind]),
  n_levels  = sapply(all_tickers, function(t) sum(!is.na(panel[[t]]))),
  n_returns = sapply(all_tickers, function(t) sum(!is.na(rets[[t]]))),
  start_date = sapply(all_tickers, function(t) format(min(panel$date[!is.na(panel[[t]])]))),
  end_date   = sapply(all_tickers, function(t) format(max(panel$date[!is.na(panel[[t]])])))
)
summary_dt[, coverage_pct := round(100 * n_returns / nrow(rets), 1)]
setorder(summary_dt, -coverage_pct)
fwrite(summary_dt, SUMMARY_OUT)

cat(sprintf("\nFrozen v1 outputs written:\n"))
cat(sprintf("  %s   (%d dates x %d series)\n", PANEL_OUT, nrow(panel), length(all_tickers)))
cat(sprintf("  %s   (%d dates x %d series)\n", RETS_OUT, nrow(rets), length(all_tickers)))
cat(sprintf("  %s\n", SUMMARY_OUT))
cat(sprintf("\nPanel date range: %s to %s (cocoa NY trading days)\n",
            min(panel$date), max(panel$date)))
cat("\nPer-series coverage on frozen returns panel:\n")
print(summary_dt)
