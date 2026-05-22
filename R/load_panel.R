suppressPackageStartupMessages({
  library(data.table)
  library(zoo)
})

load_price_panel <- function(path = NULL, roll_flags_path = NULL) {
  if (is.null(path)) {
    candidates <- c("data/raw/panel.parquet", "data/raw/panel.csv")
    path <- candidates[file.exists(candidates)][1]
    if (is.na(path)) stop("No panel.parquet or panel.csv found in data/raw/")
  }

  if (grepl("\\.parquet$", path)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Install 'arrow' to read parquet, or convert to CSV.")
    }
    dt <- as.data.table(arrow::read_parquet(path))
  } else {
    dt <- fread(path)
  }

  date_col <- intersect(c("date", "Date", "DATE"), names(dt))[1]
  if (is.na(date_col)) stop("Panel needs a 'date' column.")
  setnames(dt, date_col, "date")
  dt[, date := as.Date(date)]
  setorder(dt, date)

  ticker_cols <- setdiff(names(dt), "date")
  dt[, (ticker_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))),
     .SDcols = ticker_cols]

  # Excise futures roll days if a flag table is supplied.
  if (!is.null(roll_flags_path) && file.exists(roll_flags_path)) {
    flags <- fread(roll_flags_path)
    flags[, date := as.Date(date)]
    for (col in setdiff(names(flags), "date")) {
      if (col %in% ticker_cols) {
        roll_dates <- flags[get(col) == TRUE, date]
        dt[date %in% roll_dates, (col) := NA_real_]
      }
    }
  }

  list(prices = dt, tickers = ticker_cols)
}

compute_log_returns <- function(prices_dt, min_obs_frac = 0.6) {
  tickers <- setdiff(names(prices_dt), "date")

  # Forward-fill at most 5 days inside each series, then take log returns.
  filled <- copy(prices_dt)
  for (col in tickers) {
    filled[, (col) := na.locf(get(col), na.rm = FALSE, maxgap = 5)]
  }

  # Restrict to dates where the cocoa anchor (first column matching /COCOA/) is present.
  cocoa_col <- grep("^COCOA", tickers, value = TRUE)[1]
  if (!is.na(cocoa_col)) {
    filled <- filled[!is.na(get(cocoa_col))]
  }

  rets <- copy(filled)
  for (col in tickers) {
    p <- filled[[col]]
    r <- c(NA_real_, diff(log(p)))
    rets[, (col) := r]
  }
  rets <- rets[-1]  # drop initial NA row

  # Drop series with too few non-missing returns.
  keep <- vapply(tickers, function(col) {
    mean(!is.na(rets[[col]])) >= min_obs_frac
  }, logical(1))
  dropped <- tickers[!keep]
  if (length(dropped) > 0) {
    message(sprintf("Dropping low-coverage series: %s", paste(dropped, collapse = ", ")))
  }
  rets <- rets[, c("date", tickers[keep]), with = FALSE]

  rets
}
