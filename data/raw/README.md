# Raw data drop format

Drop the Refinitiv pull here as **one** wide-format file. Either is fine:

- `panel.csv` — header row of tickers, first column `date` (ISO `YYYY-MM-DD`).
- `panel.parquet` — same shape; one `date` column (Date or Timestamp), one column per series.

## Conventions

- Daily frequency. Calendar days where any market is closed are allowed; loader will forward-fill within each series and drop dates where the cocoa front-month is missing.
- Values are **levels** (settle prices, FX rates, yields, equity closes), not returns. Log returns are computed downstream.
- Missing values as empty string / `NA` / `NaN` — anything `readr` / `arrow` reads as missing.
- Front-month futures should be the continuous **back-adjusted** roll (or just the front contract with roll-day flagged separately). If the latter, also drop a `roll_flags.csv` with a `date` column and one boolean column per futures ticker (`TRUE` on roll days, which we will exclude).

## Series I expect (proposal §"Asset Universe")

| Category | Suggested column name |
|---|---|
| ICE London cocoa front | `COCOA_LDN` |
| ICE NY cocoa front | `COCOA_NY` |
| ICE Arabica coffee front | `COFFEE_ARA` |
| ICE Robusta coffee front | `COFFEE_ROB` |
| Sugar #11 | `SUGAR11` |
| Cotton | `COTTON` |
| Orange Juice | `OJ` |
| Class III Milk | `MILK3` |
| Corn | `CORN` |
| Soybeans | `SOYB` |
| Wheat | `WHEAT` |
| WTI Crude | `WTI` |
| USD Index | `DXY` |
| VIX | `VIX` |
| 10Y US Treasury yield | `US10Y` |
| BRL/USD | `BRL` |
| GHS/USD | `GHS` |
| XOF/USD (CFA, Ivory Coast) | `XOF` |
| VND/USD | `VND` |
| Hershey | `HSY` |
| Mondelez | `MDLZ` |
| Nestle | `NESN` |
| JM Smucker | `SJM` |
| Starbucks | `SBUX` |
| Lindt | `LISN` |

Don't worry about matching these column names exactly — the loader is column-agnostic and just uses whatever is in the file. The list is here so you can see what I'll expect when interpreting results.
