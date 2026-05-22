# Parse the 24 Refinitiv xlsx files + 1 FRED VIX CSV in data/raw/ into one wide
# panel of daily levels written to data/clean/panel.csv.
#
# Each Refinitiv export has the same shape:
#   - lines 1-7: title / RIC / interval / period / (currency)
#   - then a price-distribution or statistics table
#   - then a header row whose first cell is "Exchange Date" or "Date"
#   - then the time-series in descending-date order
#
# We find the header row by name, grab the requested column, and sort ascending.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# Manifest: (filename, ticker, value_column, value_kind)
# value_kind = "level" (futures/equity/FX/index) or "yield" (in %)
MANIFEST <- data.table(
  filename = c(
    "ICE_US_Cocoa_Futures_Price History_20260522_1656.xlsx",
    "ICE_Europe_London_Cocoa_Price History_20260522_1655.xlsx",
    "ICE_US_Coffee_Price History_20260522_1657.xlsx",
    "LIFFE_Robusta_Coffee_Price History_20260522_1659.xlsx",
    "ICE_US_Sugar_Price History_20260522_1700.xlsx",
    "ICE_US_Cotton_Price History_20260522_1701.xlsx",
    "ICE_US_FCOJ_A_Price History_20260522_1702.xlsx",
    "CBoT_Corn_Price History_20260522_1703.xlsx",
    "CBoT_Soybeans_Price History_20260522_1704.xlsx",
    "CBoT_Wheat_Price History_20260522_1706.xlsx",
    "Palm_Oil_Price History_20260522_1653.xlsx",
    "CLc1 Price History_20260522_1659.xlsx",
    "DXY sPrice History_20260522_1700.xlsx",
    "US10YT Price History_20260522_1701.xlsx",
    "GHANA Price History_20260522_1711.xlsx",
    "BRL Price History_20260522_1702.xlsx",
    "VND Price History_20260522_1707.xlsx",
    "GHANA Price History_20260522_1702.xlsx",
    "HRSHY Price History_20260522_1703.xlsx",
    "MDLZ Price History_20260522_1704.xlsx",
    "NESN Price History_20260522_1707.xlsx",
    "SJM Price History_20260522_1708.xlsx",
    "SBUX Price History_20260522_1709.xlsx",
    "LISN Price History_20260522_1710.xlsx"
  ),
  ticker = c(
    "COCOA_NY", "COCOA_LDN", "COFFEE_ARA", "COFFEE_ROB",
    "SUGAR", "COTTON", "OJ",
    "CORN", "SOYB", "WHEAT", "PALMOIL",
    "WTI", "DXY", "US10Y",
    "GHANA10Y", "BRL", "VND", "GHS",
    "HSY", "MDLZ", "NESN", "SJM", "SBUX", "LISN"
  ),
  value_col = c(
    "Close", "Close", "Close", "Close",
    "Close", "Close", "Close",
    "Close", "Close", "Close", "Close",
    "Close", "Trade Price", "BidYld",
    "BidYld", "Bid", "Bid", "Bid",
    "Close", "Close", "Close", "Close", "Close", "Close"
  ),
  value_kind = c(
    rep("level", 12),       # commodity futures + WTI
    "level", "yield",       # DXY, US 10Y yield
    "yield",                # Ghana 10Y yield
    "level", "level", "level",  # FX
    rep("level", 6)         # equities
  )
)

# ---------------------------------------------------------------------------
# Read one Refinitiv xlsx and return data.table(date, value)
# ---------------------------------------------------------------------------
parse_refinitiv_file <- function(path, value_col) {
  # Read with no skipping, no col names — we'll find the header row ourselves.
  raw <- suppressMessages(read_excel(path, col_names = FALSE, col_types = "text"))
  setDT(raw)

  # Find the header row: first column == "Exchange Date" or "Date"
  first_col <- raw[[1]]
  header_idx <- which(first_col %in% c("Exchange Date", "Date"))[1]
  if (is.na(header_idx)) stop(sprintf("No 'Exchange Date'/'Date' header found in %s", basename(path)))

  # Re-read from header_idx with header row as col names
  df <- suppressMessages(read_excel(path, skip = header_idx - 1, col_names = TRUE,
                                    col_types = "text", .name_repair = "minimal"))
  setDT(df)
  nm <- names(df)

  # Find date column and the requested value column
  date_idx  <- which(nm %in% c("Exchange Date", "Date"))[1]
  value_idx <- which(nm == value_col)
  if (length(value_idx) == 0) {
    stop(sprintf("Column '%s' not found in %s. Available: %s",
                 value_col, basename(path), paste(nm, collapse = ", ")))
  }

  out <- data.table(
    date  = suppressWarnings(as.Date(as.numeric(df[[date_idx]]), origin = "1899-12-30")),
    value = suppressWarnings(as.numeric(df[[value_idx]]))
  )
  # Sometimes the date column already comes through as a date string instead of Excel serial
  if (all(is.na(out$date))) {
    out$date <- suppressWarnings(as.Date(df[[date_idx]]))
  }
  out <- out[!is.na(date) & !is.na(value)]
  setorder(out, date)
  # Deduplicate any double-counted dates (keep last)
  out <- out[, .(value = last(value)), by = date]
  out
}

# ---------------------------------------------------------------------------
# VIX from FRED — DATE,VIXCLS format, missing values are "."
# ---------------------------------------------------------------------------
parse_vix <- function(path) {
  v <- fread(path, na.strings = c(".", "NA", ""))
  date_col <- intersect(c("DATE", "observation_date", "Date", "date"), names(v))[1]
  val_col  <- setdiff(names(v), date_col)[1]
  out <- data.table(
    date  = as.Date(v[[date_col]]),
    value = as.numeric(v[[val_col]])
  )
  out <- out[!is.na(date) & !is.na(value)]
  setorder(out, date)
  out
}

# ---------------------------------------------------------------------------
# Build the wide panel
# ---------------------------------------------------------------------------
build_panel <- function(raw_dir = "data/raw") {
  series_list <- list()
  diag <- data.table()

  for (i in seq_len(nrow(MANIFEST))) {
    row <- MANIFEST[i]
    path <- file.path(raw_dir, row$filename)
    s <- parse_refinitiv_file(path, row$value_col)
    setnames(s, "value", row$ticker)
    series_list[[row$ticker]] <- s
    diag <- rbind(diag, data.table(
      ticker = row$ticker, value_col = row$value_col, value_kind = row$value_kind,
      n_obs = nrow(s), start = min(s$date), end = max(s$date)
    ))
  }

  # VIX
  vix <- parse_vix(file.path(raw_dir, "VIXCLS.csv"))
  setnames(vix, "value", "VIX")
  series_list[["VIX"]] <- vix
  diag <- rbind(diag, data.table(
    ticker = "VIX", value_col = "VIXCLS", value_kind = "level",
    n_obs = nrow(vix), start = min(vix$date), end = max(vix$date)
  ))

  # Outer join all series on date
  panel <- Reduce(function(a, b) merge(a, b, by = "date", all = TRUE), series_list)
  setorder(panel, date)

  list(panel = panel, diagnostics = diag)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  result <- build_panel()
  panel <- result$panel
  diag  <- result$diagnostics

  dir.create("data/clean", showWarnings = FALSE, recursive = TRUE)
  fwrite(panel, "data/clean/panel.csv")
  fwrite(diag,  "data/clean/panel_diagnostics.csv")

  cat(sprintf("\nWrote data/clean/panel.csv: %d dates x %d series\n",
              nrow(panel), ncol(panel) - 1L))
  cat(sprintf("Date range: %s to %s\n", min(panel$date), max(panel$date)))

  cat("\nSeries coverage:\n")
  diag[, coverage_pct := round(100 * n_obs / nrow(panel), 1)]
  print(diag[order(-n_obs)])

  # Head & tail of the panel
  cat("\nFirst 3 rows:\n"); print(head(panel, 3))
  cat("\nLast 3 rows:\n");  print(tail(panel, 3))
}
