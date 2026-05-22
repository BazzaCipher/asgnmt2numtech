# Refinitiv pull guide

## Pull settings (apply to every series)

| Setting | Value |
|---|---|
| Frequency | Daily |
| Date range | 2015-01-01 → 2026-05-22 |
| Field — futures | Settlement Price (`TR.SETTLEMENTPRICE` in Workspace; `SETTLE` in old Eikon) |
| Field — equities | Close (`TR.PriceClose` adjusted; we want **adjusted close** so splits and dividends don't leak into returns) |
| Field — FX, yields, indices | Close / level |
| Currency | Native (no FX conversion for equities — Nestle in CHF, Lindt in CHF, the rest in USD; we work in log returns so currency unit doesn't matter as long as it's consistent within series) |
| Calendar | Trading days of the source exchange; we'll align on cocoa's trading days downstream |
| Futures continuation | **Back-adjusted continuous** front-month (see note below) |

### Futures continuation: which convention to pick

Two options in Workspace, and they matter:

1. **`LCCc1` (front-month, no back-adjust)** — easy to pull, but log returns on roll days are contaminated by the price gap between the expiring contract and the new front month. Use this **only** if you also export a roll-date flag so we can null those days.
2. **Back-adjusted continuous (`LCC1` in Datastream-style, or use Workspace's "Generic continuation, back-adjusted, ratio" setting)** — clean returns across rolls. **Preferred.** No flag file needed.

If your pull supports it, go with back-adjusted continuous for every futures series below. If not, give me the `c1` series plus a `roll_flags.csv` with one boolean column per futures ticker (TRUE on the roll day for that contract).

---

## Tickers

### Softs (ICE)
| Series | RIC (front-month) | Notes |
|---|---|---|
| ICE London cocoa | `LCCc1` | GBP/ton |
| ICE NY cocoa | `CCc1` | USD/ton |
| ICE Arabica coffee | `KCc1` | USc/lb |
| ICE Robusta coffee | `LRCc1` | USD/ton |
| Sugar #11 (raw) | `SBc1` | USc/lb |
| Cotton #2 | `CTc1` | USc/lb |
| FCOJ | `OJc1` | USc/lb |

### Dairy (CME)
| Series | RIC | Notes |
|---|---|---|
| Class III Milk | `DAc1` | USD/cwt |

### Grains (CBOT)
| Series | RIC | Notes |
|---|---|---|
| Corn | `Cc1` | USc/bu |
| Soybeans | `Sc1` | USc/bu |
| Wheat (SRW) | `Wc1` | USc/bu |

### Energy + macro
| Series | RIC | Notes |
|---|---|---|
| WTI Crude (NYMEX front) | `CLc1` | USD/bbl |
| USD Index (ICE) | `.DXY` | Spot index, not the futures |
| VIX | `.VIX` | CBOE |
| 10Y US Treasury yield | `US10YT=RR` | Yield in %, we'll convert to first-differences |

### FX (USD-base in Refinitiv convention; we'll invert downstream if needed — keep as-pulled)
| Series | RIC | Notes |
|---|---|---|
| BRL | `BRL=` | USD/BRL |
| GHS | `GHS=` | USD/GHS — *see note below* |
| XOF (CFA) | `XOF=` | USD/XOF — *see note below* |
| VND | `VND=` | USD/VND |

> **GHS/XOF caveat (re-flagged from chat).** XOF is hard-pegged to EUR via the CFA franc arrangement, so `XOF=` is basically a scaled `EUR=` inverse and won't carry an independent cocoa signal. GHS is managed-float but historically has long stretches of admin-driven flatness. I'd suggest replacing one or both with:
> - **Ghana 10Y USD sovereign bond yield** (RIC chain `GH10YT=RR` or via the bond list `0#GHANA=`) — more responsive to cocoa export revenue shocks
> - **Côte d'Ivoire EUR sovereign** (e.g. `XS1245045136` or via `0#CIV=`) — similar logic
>
> Pulling these as **mid-yield** daily would give the same shape of "West Africa stress" signal but with actual variation. If unsure, **pull both the FX and the bond yields** — costs nothing extra, and we can pick during the day-1 analysis.

### Equities
| Series | RIC | Exchange |
|---|---|---|
| Hershey | `HSY.N` | NYSE |
| Mondelez | `MDLZ.OQ` | NASDAQ |
| Nestle | `NESN.S` | SIX (CHF) |
| JM Smucker | `SJM.N` | NYSE |
| Starbucks | `SBUX.OQ` | NASDAQ |
| Lindt (PS) | `LISN.S` | SIX (CHF) |

> Lindt has two listings on SIX — the registered share (`LISP.S`) and the participation certificate (`LISN.S`). `LISN.S` is the more liquid; use that.

---

## Output format

One wide-format file, `panel.csv` or `panel.parquet`, into `data/raw/`:

- First column: `date` (ISO `YYYY-MM-DD`)
- Other columns: one per series, named however you like (the loader is column-agnostic — but the suggested names in `data/raw/README.md` will read more cleanly in figures)
- Values: levels, not returns

If you're pulling via the Excel add-in, a simple way:

```
=@TR(A2:A26, "TR.SettlementPrice; TR.PriceClose; ...",
     "SDate=2015-01-01 EDate=2026-05-22 Frq=D CH=Fd")
```

where `A2:A26` is a column of RICs. Workspace's "Data Item Browser → Pricing → Daily Price History" with multiple instruments also exports straight to a wide CSV.

---

## Sanity checks before you hand it over

Spot-check three things and you can save a round-trip:

1. Cocoa NY (`CCc1`) on **2024-12-18** should settle around **USD 12,500-12,900/ton** (the all-time high).
2. Cocoa NY on **2024-04-19** should be around **USD 11,000/ton** (April spike).
3. `.DXY` on **2023-01-03** should be around **103.8**.

If those are off by more than a few percent, the field/adjustment is wrong.
