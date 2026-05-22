# Data dictionary — Cocoa RMT-DCC project

*Frozen 2026-05-22. Source files in `data/raw/`; cleaned outputs in `data/clean/`.*

## Files

| File | Shape | Source |
|---|---|---|
| `data/clean/panel_v1.csv` | 2864 × 26 (date + 25 series) | Frozen levels, cocoa-NY calendar |
| `data/clean/returns_v1.csv` | 2863 × 26 | Frozen returns (log / first-diff) |
| `data/clean/std_residuals_v1.csv` | 2863 × 26 | GJR-GARCH-t standardised residuals, winsorised at ±5 |
| `data/clean/panel_v1_summary.csv` | per-series coverage report | |
| `data/clean/panel.csv` | 2974 × 26 | Raw merged panel before calendar restriction |
| `data/clean/panel_diagnostics.csv` | series manifest | Source filename, value column, value kind |
| `output/garch_fits.csv` | per-asset GARCH parameters + diagnostics | |

## Pipeline

1. `R/parse_refinitiv.R` reads the 24 `.xlsx` files + `VIXCLS.csv` → `panel.csv`.
2. `R/clean_panel.R` applies locked cleaning rules → `panel_v1.csv` and `returns_v1.csv`.
3. `R/fit_garch.R` fits a univariate GJR-GARCH(1,1) with Student-t innovations per asset (fallback chain: Student-t → Gaussian → EGARCH → EWMA). Standardised residuals are winsorised at ±5σ to prevent any single tail event from dominating the rolling correlation matrices. Outputs `std_residuals_v1.csv` and `garch_fits.csv`.

### GARCH fit summary

| Tier used | Count | Notes |
|---|---|---|
| GJR-GARCH(1,1) Student-t | 23/25 | Default; all parameters in stationary region (α+β+½γ < 1). |
| GJR-GARCH(1,1) Gaussian | 1/25 | VND only — Student-t shape collapsed (managed-currency jumps). |
| EWMA (λ=0.94) | 1/25 | GHANA10Y — all GARCH variants failed the residual quality check (sd > 2 or |z| > 30) due to the 2022-23 sovereign-debt restructuring period. |

Total winsorised observations: **150** out of ~71,500 (0.21%). Per-series winsorisation counts in `output/garch_fits.csv:n_winsorized`.

## Locked cleaning rules

| Rule | Detail |
|---|---|
| Trading-day calendar | Restrict to dates where `COCOA_NY` is observed (cocoa is the analysis anchor) |
| Forward fill | Within each series, max 3 consecutive missing days |
| Level returns | log diff: r_t = ln(P_t) − ln(P_{t−1}) |
| Yield returns | first-difference of % yield: r_t = y_t − y_{t−1}, in percentage points |
| WTI patch | Returns for 2020-04-20 and 2020-04-21 set to NA (front-month settled at -$37.63 on 2020-04-20; log undefined) |
| Date range | 2015-01-02 to 2026-05-21 (cocoa NY trading days only) |

## Series

| Ticker | Full name | RIC | Source file | Native unit | Value column | Value kind |
|---|---|---|---|---|---|---|
| `COCOA_NY` | ICE-US Cocoa Futures front-month | `CCc1` | `ICE_US_Cocoa_Futures_Price History_…xlsx` | USD/ton | Close | level |
| `COCOA_LDN` | ICE London Cocoa front-month | `LCCc1` | `ICE_Europe_London_Cocoa_Price History_…xlsx` | GBP/ton | Close | level |
| `COFFEE_ARA` | ICE-US Arabica Coffee front-month | `KCc1` | `ICE_US_Coffee_Price History_…xlsx` | USc/lb | Close | level |
| `COFFEE_ROB` | LIFFE Robusta Coffee front-month | `LRCc1` | `LIFFE_Robusta_Coffee_Price History_…xlsx` | USD/ton | Close | level |
| `SUGAR` | ICE-US Sugar #11 front-month | `SBc1` | `ICE_US_Sugar_Price History_…xlsx` | USc/lb | Close | level |
| `COTTON` | ICE-US Cotton #2 front-month | `CTc1` | `ICE_US_Cotton_Price History_…xlsx` | USc/lb | Close | level |
| `OJ` | ICE-US FCOJ-A front-month | `OJc1` | `ICE_US_FCOJ_A_Price History_…xlsx` | USc/lb | Close | level |
| `CORN` | CBOT Corn front-month | `Cc1` | `CBoT_Corn_Price History_…xlsx` | USc/bu | Close | level |
| `SOYB` | CBOT Soybeans front-month | `Sc1` | `CBoT_Soybeans_Price History_…xlsx` | USc/bu | Close | level |
| `WHEAT` | CBOT SRW Wheat front-month | `Wc1` | `CBoT_Wheat_Price History_…xlsx` | USc/bu | Close | level |
| `PALMOIL` | Bursa Malaysia Crude Palm Oil front-month | `FCPOc1` | `Palm_Oil_Price History_…xlsx` | MYR/ton | Close | level |
| `WTI` | NYMEX Light-Sweet Crude Oil front-month | `CLc1` | `CLc1 Price History_…xlsx` | USD/bbl | Close | level |
| `DXY` | ICE US Dollar Currency Index | `.DXY` | `DXY sPrice History_…xlsx` | index level | Trade Price | level |
| `VIX` | CBOE Volatility Index | `VIXCLS` | `VIXCLS.csv` (FRED) | index level (annualised % vol) | VIXCLS | level |
| `US10Y` | US Treasury 10-Year Benchmark Bid Yield | `US10YT=RR` | `US10YT Price History_…xlsx` | % yield | BidYld | **yield** |
| `BRL` | USD/Brazilian Real Spot | `BRL=` | `BRL Price History_…xlsx` | BRL per 1 USD | Bid | level |
| `VND` | USD/Vietnamese Dong Spot | `VND=` | `VND Price History_…xlsx` | VND per 1 USD | Bid | level |
| `GHS` | USD/Ghanaian Cedi Spot | `GHS=` | `GHANA Price History_…xlsx` (17:02) | GHS per 1 USD | Bid | level |
| `GHANA10Y` | Ghana 10-Year Benchmark Bid Yield | `GH10YT=RR` | `GHANA Price History_…xlsx` (17:11) | % yield | BidYld | **yield** |
| `HSY` | Hershey Co. | `HSY.N` | `HRSHY Price History_…xlsx` | USD | Close | level |
| `MDLZ` | Mondelez International | `MDLZ.OQ` | `MDLZ Price History_…xlsx` | USD | Close | level |
| `NESN` | Nestle SA (SIX Swiss) | `NESN.S` | `NESN Price History_…xlsx` | CHF | Close | level |
| `SJM` | J.M. Smucker Co. | `SJM.N` | `SJM Price History_…xlsx` | USD | Close | level |
| `SBUX` | Starbucks Corp. | `SBUX.OQ` | `SBUX Price History_…xlsx` | USD | Close | level |
| `LISN` | Chocoladefabriken Lindt & Sprüngli AG (Participation Cert.) | `LISN.S` | `LISN Price History_…xlsx` | CHF | Close | level |

## Coverage on the frozen returns panel

(see `data/clean/panel_v1_summary.csv` for the machine-readable version)

| Tier | Tickers | Coverage |
|---|---|---|
| Full | All except WTI and GHANA10Y | 100% (2863 / 2863 days) |
| Patched | WTI | 99.9% (2 days NA: 2020-04-20, 2020-04-21) |
| Partial | GHANA10Y | 79.0% (2262 days, starts 2017-04-20) |

`NESN` and `LISN` levels start one trading day later (2015-01-05) because the Swiss exchange was closed 2015-01-02; this affects 1 row of returns on each.

## Known quirks (referenced by the paper)

1. **WTI 2020-04-20:** Front-month settle = −$37.63/bbl due to the COVID-era Cushing storage glut. Real data point. Returns NA'd on 2020-04-20 and 2020-04-21 because the log return formula requires positive prices on both ends of the difference.
2. **Ghana 10Y short history:** GHANA10Y starts April 2017, ~74% of the panel period. Rolling 252-day windows entirely before April 2018 use N=24 series; from April 2018 onwards, N=25. The MP edge λ₊ = (1+√(N/T))² is recomputed per window, but the count is N-dependent. Day-10 robustness re-runs analysis excluding GHANA10Y to isolate the effect.
3. **Lindt PS price level:** `LISN.S` trades at ~CHF 95,000-105,000 per share (the participation certificate is one of Switzerland's highest-priced equities). Not a data error.
4. **FX convention:** `BRL`, `VND`, `GHS` are quoted as units of local currency per 1 USD (Refinitiv default). Higher values = USD strength. A shared USD-strength factor will show up as a common eigenvector component.
5. **Equity time-zone non-synchroneity:** US equities close after the European/Swiss market; correlations between (NESN, LISN) and (HSY, MDLZ, SJM, SBUX) at daily horizon may be biased downward by ~3-4 hours of stale information. Acknowledged limitation.
6. **Palm oil currency:** `FCPOc1` quoted in MYR (Malaysian Ringgit). Returns therefore include MYR/USD movement.
