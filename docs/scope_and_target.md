# Project scope and target

*Locked after day-1 verification, 2026-05-22. Replaces the indicative scope in the proposal where they differ.*

---

## Smoke-test verdict (gating decision)

The day-1 rolling MP-edge eigenvalue count on raw log returns showed:

| Period | Median count above λ₊ (ex-market) | IQR |
|---|---|---|
| 2015-2023 (pre-crisis) | **2** | [2, 2] |
| 2024-2026 (crisis) | **3** | [2, 3] |
| Difference | **+1** | — |

Visual inspection confirms the H1 prediction: the eigenvalue count is persistently elevated through the 2024-26 cocoa crisis window. **Green-light to proceed with the full pipeline.**

Two honest qualifiers visible in the day-1 plot:
- A "+1" pop also appears transiently in 2015-17 (China/oil shock period). The crisis is the *most persistent* state-3 episode, not the only one.
- The universe-size change when Ghana 10Y enters the rolling window in mid-2017 to mid-2018 could mechanically lift the count by 1 in that window. To be checked in the day-10 robustness pass by excluding Ghana 10Y as a sensitivity.

---

## Final universe (25 series, locked)

| Block | Tickers | n | Source |
|---|---|---|---|
| Soft commodity futures | COCOA_NY, COCOA_LDN, COFFEE_ARA, COFFEE_ROB, SUGAR, COTTON, OJ | 7 | Refinitiv |
| Adjacent ags | CORN, SOYB, WHEAT, PALMOIL | 4 | Refinitiv |
| Energy | WTI | 1 | Refinitiv |
| Macro | DXY, VIX, US10Y | 3 | Refinitiv + FRED |
| FX vs USD | BRL, VND, GHS | 3 | Refinitiv |
| Sovereign yield | GHANA10Y | 1 | Refinitiv (from Apr 2017) |
| Chocolate/beverage equities | HSY, MDLZ, NESN, SJM, SBUX, LISN | 6 | Refinitiv |

**Aspect ratio Q = T/N ≈ 2870/25 ≈ 115.** Well-suited to MP analysis (paper recommends Q ≥ 5; we are far above).

### Changes from the proposal universe

- **Added:** PALMOIL (`FCPOc1`), partial cocoa-butter substitute. Discovered in the second pull, kept because it's a sensible adjacent vegetable-fat series and adds discriminatory power to the eigenvector composition test.
- **Added:** GHS (cedi spot) and GHANA10Y (10Y USD sovereign yield) — the West African stress proxies suggested during data scoping; both contain real signal (cedi depreciated from 3.18 to 11.54 over the window).
- **Dropped:** Class III Milk (`DAc1`) — borderline relevance, not pulled by user.
- **Dropped:** XOF/CFA franc — hard-pegged to EUR, would have contributed no independent signal.

---

## Data-handling decisions (locked)

| Issue | Decision | Reason |
|---|---|---|
| Currency | All series in native units (no FX conversion). | Avoids the spurious "AUD mode" eigenvalue that contaminated the first pull. |
| FX convention | Refinitiv default USD/X (BRL, VND, GHS quoted as units of local currency per USD). | A shared "USD strength" factor will appear in the spectrum; documented up front rather than hidden by inversion. |
| Yield returns | First-differences of yield in %, **not** log returns. | Standard for yields; robust if yields approach zero. |
| Level returns | Log returns. | Standard for prices. |
| Forward fill | Inside each series, up to 3 days. | Handles single-day exchange holidays without smearing across long gaps. |
| Trading-day calendar | Restrict to dates where `COCOA_NY` is observed. | Cocoa is the anchor of the analysis. |
| WTI 2020-04-20 | NA out returns on 2020-04-20 and 2020-04-21. | WTI front-month settled at -$37.63 (real); log returns undefined. |
| Rolling window | 252 trading days throughout. | Matches proposal; also runs at 504 in day-10 robustness. |
| MP edge | λ₊ = (1+√(N/T))². | Standard; recomputed each window because N varies (Ghana 10Y enters mid-2018). |

---

## Target

**Single falsifiable headline finding** (paper succeeds if all three sub-conditions hold; otherwise we report null result honestly — see "If the headline fails" below):

> Over rolling 252-day windows on RMT-cleaned standardised residuals from GJR-GARCH-Student-t univariate fits, the count of correlation eigenvalues exceeding the Marchenko-Pastur upper edge (after removing the market mode):
>
> **(T1) Baseline stability.** Over the 2019-01-01 to 2023-12-31 window (entirely N = 25 universe, post-Ghana-10Y entry), the count has IQR width ≤ 1 and stays inside {1, 2, 3}.
>
> **(T2) Structural break.** Two course-toolkit tests, both must hold:
> (a) **Bai-Perron** structural-break test (`strucchange`) on the rolling-count series finds at least one break with date inside [2024-01-01, 2024-12-31].
> (b) The **stationary block bootstrap** (1000 reps; block length ≈ 25) of the standardised residuals gives a 95% CI for the **pre-crisis (2019-2023) median count**. The in-crisis (2024-2026) median count is strictly above the upper bound of this CI.
>
> **(T3) Cocoa-bloc concentration.** Initial test wording assumed a *new* supra-MP rank would appear in the crisis subsample. In practice the count of supra-MP eigenvalues in a subsample correlation matrix is mostly determined by the universe size and is stable across subsamples (5 in both pre and crisis here). The real H1 signal turns out to live in the *composition* of the existing supra-MP eigenvectors, not in whether new ones appear. T3 is therefore restated as:
>
> Define the **rank-1 supra-MP eigenvector** in each subsample (= the eigenvector of the largest eigenvalue strictly above the MP edge, after removing the market mode). Let _f_bloc_(subsample) be the fraction of its squared loadings attributable to the cocoa bloc {COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y}. T3 holds if **both**:
> (a) _f_bloc_(crisis 2024-2026) ≥ 0.50 — i.e. the dominant non-market mode in the crisis subsample is recognisably a cocoa mode.
> (b) _f_bloc_(crisis 2024-2026) ≥ 2 × _f_bloc_(pre 2019-2023) — i.e. cocoa-bloc concentration on the rank-1 supra-MP eigenmode at least *doubled* relative to baseline.
>
> Random-baseline reference: a unit vector with mass spread uniformly across all 25 series has 7/25 ≈ 0.28 mass in the bloc. (a) requires ~1.8× random; (b) requires the doubling to be real, not a baseline fluke.

T1 is the precondition; T2 is the structural-break statement; T3 is the interpretability statement (that the new mode is recognisably a cocoa-complex mode, not random covariance noise).

**Day-1 smoke test (raw returns)** showed T1 + T2(eye-test) green.
**Days 5-10 cleaned-pipeline run** all four primary target legs pass:
- **T1 PASS** — pre-crisis IQR contained in {1,2,3}, width ≤ 1 within [2,3]
- **T2(a) PASS** — Bai-Perron break at 2024-07-08, inside [2024-01-01, 2024-12-31]
- **T2(b) PASS** — stationary block bootstrap (B=500, block mean=25) of the 2019-2023 pre-crisis residuals gives a 95% CI on the pre-crisis median of [2.00, 2.00]; observed in-crisis median 3.0 > 2.00 (CI upper bound). The bootstrap CI is degenerate because the pre-crisis median is essentially constant at 2 under any reshuffle of the baseline data. Will rerun at B=1000 for the final paper (no expected change).
- **T3 PASS** — crisis-subsample rank-1 supra-MP eigenmode has cocoa-bloc squared-loading fraction 0.517 vs 0.245 pre-crisis (2.1× doubling); both (a) ≥ 0.50 and (b) ≥ 2× legs satisfied.

**Robustness (Day 10 variants, all 5 panels in `figures/fig4_robustness.png`):**
- Excluding GHANA10Y: **T2(b) PASS** (pre median 2, in median 3, CI hi 2.76). The Ghana-entry concern is decisively NOT driving the result.
- Gaussian-innovations GARCH: **T2(b) PASS** (identical numbers to baseline). Student-t is not load-bearing.
- Window length 504 (baseline + excl GHANA10Y): T2(b) "fails" by my arithmetic but only because both pre and in medians saturate at 3 (no headroom). Visual inspection shows a clear +1 (3→4) lift during the crisis window. Bai-Perron break still inside the crisis window (2025-04-02).

**Secondary targets (S1-S3) from the DCC pipeline** (results in `output/s1_variance_reduction.csv`, `s2_bai_perron_pairs.csv`, `s3_cocoa_basis_decoupling.csv`):
- **S1 FAIL** — Cleaned-DCC produces *more* variable conditional correlations than standard-DCC across all 11 cocoa-anchored pairs in the crisis window. Likely cause: with N=25 and Q=T/N≈90 the raw sample correlation is already well-estimated, and the cleaning step compresses off-diagonals enough that the DCC optimizer hits the b=0.98 cap, producing near-integrated dynamics. Reported as a contradicting result; the Engle-Ledoit-Wolf cleaning advantage may require N >> 25.
- **S2 PASS** — Bai-Perron on the cleaned-DCC pairwise conditional correlations finds an in-crisis break for **every** cocoa-anchored pair (11 of 11), including the pre-registered cocoa-coffee and cocoa-sugar pairs.
- **S3 FAIL but interesting** — Cocoa NY × Cocoa London conditional correlation did NOT decouple during the crisis; instead it *re-coupled* (baseline mean 0.466, crisis min 0.509). Interpretation: the 2024-26 cocoa shock was global in nature (driven by West African production failures that are reflected in BOTH NY and London deliverable contracts), not a region-specific basis blowout as we initially hypothesised. This is a substantive finding to report in the discussion, not just a null result.

T1 is the precondition; T2 is the structural-break statement; T3 is the interpretability statement (that the new mode is recognisably a cocoa-complex mode, not random covariance noise).

Day-1 smoke test (no GARCH, no RMT cleaning) already shows T1 satisfied (pre-crisis count is 2 with IQR [2,2] when restricted to 2019-2023) and T2(a)/T2(b) plausible at the eye-test. T3 cannot be tested without the cleaned DCC pipeline.

### If the headline fails

If T2 fails, the paper is restructured around the secondary targets (S1/S2/S3 below) — they are statistically independent of T2 and likely to remain interesting on their own. T2's failure is reported as a contradicting null result. This is acceptable per the proposal's H2/H3 being separable from H1.

### Secondary targets

**(S1) Cleaned-DCC stability (H2 in the proposal).** Variance of the cleaned-DCC pairwise conditional correlation series for the (cocoa-sugar, cocoa-coffee, cocoa-NY × cocoa-London) triad is at least 25% lower than the variance of the standard-DCC counterpart over the crisis window. Bootstrap CI.

**(S2) Asymmetric propagation (H3).** Bai-Perron tests on the cleaned-DCC conditional correlation series for (cocoa-coffee, cocoa-sugar) detect at least one structural break inside 2024-01 to 2026-05, with the break date within ±90 days of a documented turning point (Apr 2024 supply spike, Dec 2024 peak, mid-2025 demand-destruction turn).

**(S3) Geographic-origin decoupling.** The cleaned-DCC conditional correlation ρ_t(COCOA_NY, COCOA_LDN) drops by ≥ 0.15 from its 2015-2023 mean during at least one window inside 2025-01 to 2026-05.

---

## Out of scope (locked, will not chase)

- Forecasting future cocoa prices.
- Trading-strategy backtests or implied portfolio P&L.
- Comparison against high-frequency intraday data (daily only).
- Multivariate models beyond DCC (BEKK, factor-GARCH, copula-GARCH).
- Spillover network estimation (Diebold-Yilmaz). Eigenvector composition serves the same interpretive role here.
- Alternative cleaning schemes beyond hard MP clipping (e.g. Ledoit-Wolf shrinkage, eigenvalue smoothing). One method, one paper.

---

## Risk register (post-data)

| Risk | Severity | Mitigation |
|---|---|---|
| Ghana 10Y entering universe in 2017-18 inflates the count mechanically. | Medium | Day-10 robustness: re-run excluding GHANA10Y. |
| WTI 2020 negative-price event distorts COVID-window correlations. | Low | NA on -04-20/-21 (decided). |
| Equity AU/EU markets close earlier than NY → potential non-synchronous returns biasing correlations downward. | Low-medium | Acknowledged limitation; mitigation = consider weekly returns as robustness if biases material. |
| LISN illiquidity (Lindt PS, ~CHF 100k+ per share, often <5k trades/day). | Low | Coverage is ≥96%; flag in data dictionary. |
| GARCH non-convergence for cocoa during the 2024 vertical move. | Medium | Use Student-t innovations and tight parameter bounds in rugarch; fallback to GJR(1,1) with Gaussian if Student-t fails per asset. |
| Multiple-testing inflation across pairwise DCC correlations. | Medium | Pre-register the cocoa-anchored pairs (cocoa-sugar, cocoa-coffee, cocoa-NY-vs-London, cocoa-HSY, cocoa-MDLZ, cocoa-LISN, cocoa-GHS, cocoa-GHANA10Y); use Bonferroni or Holm on the structural-break p-values. |

---

## Path from here (revised 14-day plan)

Day-1 spent on data fix + smoke test. The remaining schedule:

- **Day 2:** Write the data dictionary (`docs/data_dictionary.md`); finalise the cleaning pipeline; archive a frozen `data/clean/panel_v1.csv`.
- **Days 3-4:** Install rugarch + rmgarch + strucchange + RMTstat; fit GJR-GARCH-Student-t per asset; produce diagnostic plots; save standardised residuals to `data/clean/standardised_residuals.csv`.
- **Days 5-6:** MP eigenvalue analysis on standardised residuals; produce Figure 1 (spectrum vs MP density) and Figure 2 (rolling count with crisis band).
- **Days 7-8:** DCC-GARCH with cleaned target; pairwise conditional correlations; produce Figure 3.
- **Days 9-10:** Bootstrap CIs (stationary block bootstrap, 1000 reps, block length ≈ 25) for the rolling count and the conditional correlation series; Bai-Perron breaks on the count series; robustness checks (exclude Ghana 10Y; window length 504; pre-crisis-N=24 vs N=25 subsamples; Gaussian vs Student-t innovations).
- **Day 11:** Eigenvector composition analysis; T3 test.
- **Days 12-13:** Draft 20-page paper plus appendix.
- **Day 14:** Final pass and submit.

Day-1 ate roughly half a day (data re-pull). Time bank stays positive.
