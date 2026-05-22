# Cocoa RMT-DCC project

FNCE40003 major assignment. See [`cocoa_rmt_dcc_proposal.md`](cocoa_rmt_dcc_proposal.md) for the full brief.

## Layout

```
R/                 R scripts (load_panel, day1_smoke_test, …)
data/raw/          Refinitiv panel goes here (see data/raw/README.md)
data/clean/        Processed return panel, standardised residuals, etc.
figures/           Figure 1–3 PNGs
output/            CSVs of rolling eigenvalue count, conditional correlations, etc.
docs/              Draft writeup
```

## Day-1 smoke test (gating)

Run from project root:

```
Rscript R/day1_smoke_test.R
```

Inputs: `data/raw/panel.csv` or `data/raw/panel.parquet` (see `data/raw/README.md`).
Outputs: `output/day1_rolling_eig_count.csv`, `figures/day1_rolling_eig_count.png`, and a verdict to stdout comparing pre-crisis vs in-crisis eigenvalue counts.

If the in-crisis median count is not clearly above the pre-crisis median, the proposal's premortem says: stop and escalate before committing to the full GARCH pipeline.
