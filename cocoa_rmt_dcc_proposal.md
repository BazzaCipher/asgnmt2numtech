# Identifying the 2024 to 2026 Cocoa Crisis as a Structural Break in Soft-Commodity Correlation Dynamics

*A Random-Matrix-Cleaned DCC-GARCH Approach*

**FNCE40003 Numerical Techniques in Finance. Major Assignment Proposal.**

---

## Research Question

Did the 2024 to 2026 cocoa crisis (peak-to-trough decline of approximately 66% in cocoa futures prices) produce an identifiable structural break in the dependence structure of the global soft-commodity complex? Specifically, does the crisis manifest as an emerging eigenvalue above the Marchenko-Pastur upper edge in the correlation matrix of soft-commodity futures and related macro variables?

## Why This Setup Works

The cocoa crisis is the largest amplitude soft-commodity shock in over fifty years. Cocoa futures rose from approximately USD 3,700 per ton in late 2023 to USD 12,906 per ton in late 2024, then collapsed to near USD 3,000 by early 2026. The shock has multiple identifiable drivers (El Niño, swollen shoot virus, Ghana-Ivory Coast policy divergence, demand destruction) and clean structural break dates. The cocoa-coffee-sugar correlation cluster has been documented in prior work (Arfaoui 2013, "Volatility in a Mug Cup" 2024) but never under the random-matrix-cleaned DCC framework introduced by Engle, Ledoit, and Wolf (2019), and never with the 2024 to 2026 data.

## Hypotheses

**H1 (eigenvalue emergence).** The empirical correlation matrix of the soft-commodity universe exhibits a stable count of supra-threshold eigenvalues (above the Marchenko-Pastur upper edge) during normal periods (2015 to 2023). During the 2024 to 2026 cocoa crisis, this count temporarily increases by one, with the additional eigenvalue's eigenvector loading on cocoa contracts, West African currencies, and chocolate-maker equities. The emergence reverses as the crisis resolves.

**H2 (cleaned DCC stability).** RMT-cleaned DCC-GARCH produces more stable conditional correlation estimates than standard DCC-GARCH during the high-volatility crisis window, measured by lower variance in conditional correlation series, lower implied minimum-variance portfolio turnover, and tighter bootstrap confidence bands.

**H3 (asymmetric propagation).** Cocoa-sugar and cocoa-coffee conditional correlations strengthened during the 2024 supply-driven rally (positive comovement under shared agricultural risk) but decoupled during the 2025 to 2026 demand-destruction phase, identifiable as a structural break in the time-varying conditional correlation series.

## Asset Universe

Approximately 25 series, daily, January 2015 through May 2026 (about 2,850 observations).

**Soft commodity futures:** ICE London cocoa front-month, ICE New York cocoa front-month, ICE Arabica coffee front-month, ICE Robusta coffee front-month, Sugar #11, Cotton, Orange Juice, Class III Milk. **Adjacent agricultural commodities:** Corn, Soybeans, Wheat. **Macro variables:** WTI Crude, USD Index, VIX, 10-year US Treasury yield. **Currencies:** BRL/USD, GHS/USD (Ghanaian cedi), XOF/USD (CFA franc, Ivory Coast), VND/USD. **Chocolate and beverage equities:** Hershey (HSY), Mondelez (MDLZ), Nestle (NESN), JM Smucker (SJM), Starbucks (SBUX), Lindt (LISN).

Aspect ratio Q = T/N is approximately 110, well-suited to MP analysis.

## Methodology

### 1. Univariate volatility (GJR-GARCH per asset)

For each asset, fit a GJR-GARCH(1,1) with Student-t innovations:

$$h_{i,t} = \omega_i + \alpha_i \epsilon_{i,t-1}^2 + \gamma_i \epsilon_{i,t-1}^2 \mathbb{1}_{\epsilon_{i,t-1} < 0} + \beta_i h_{i,t-1}$$

Asymmetric specification captures leverage effects expected in cocoa (sharper vol response to downward moves during the 2025 crash). Student-t handles fat tails. Standardize residuals: $z_{i,t} = \epsilon_{i,t} / \sqrt{h_{i,t}}$.

### 2. MP-cleaned correlation matrix

Compute the sample correlation matrix $\bar{R}$ of standardized residuals over a rolling 252-day window. Eigendecompose: $\bar{R} = V \Lambda V'$.

The Marchenko-Pastur upper edge for the correlation matrix is $\lambda_+ = (1 + \sqrt{N/T})^2$.

Hard MP clipping (Bouchaud-Potters): keep eigenvalues with $\lambda_k > \lambda_+$; replace eigenvalues with $\lambda_k \leq \lambda_+$ by their average. Reconstruct the cleaned matrix $\tilde{R} = V \tilde{\Lambda} V'$. Largest eigenvalue (market mode) is removed before MP fitting to avoid spectrum compression, following standard practice.

### 3. DCC-GARCH with cleaned target

DCC dynamics:

$$Q_t = (1 - a - b) \tilde{R} + a z_{t-1} z_{t-1}' + b Q_{t-1}$$

$$R_t = \text{diag}(Q_t)^{-1/2} Q_t \text{diag}(Q_t)^{-1/2}$$

The key change from standard DCC is the use of the cleaned $\tilde{R}$ rather than the raw sample $\bar{R}$ as the long-run correlation target. Parameters $a$ and $b$ estimated by Gaussian quasi-maximum likelihood.

### 4. Rolling eigenvalue count (H1 test)

On each day in the sample, compute the rolling 252-day correlation matrix of standardized residuals, count eigenvalues above the MP upper edge, and plot the time series of this count. Bai-Perron structural break tests on the count series to formally date regime changes. Block bootstrap of the standardized residuals (1000 replications) to obtain confidence bands.

### 5. Eigenvector composition (H1 interpretation)

For the eigenvalue emerging during the crisis window, examine its associated eigenvector. Report asset loadings. Test whether the loadings match the hypothesised cocoa-FX-equity pattern.

### 6. Conditional correlation regime tests (H3)

Extract pairwise conditional correlation series for cocoa-sugar, cocoa-coffee, cocoa-HSY, and similar pairs. Bai-Perron tests for structural breaks in each series. Compare break dates to crisis chronology.

## Course Toolkit Mapping

- **GJR-GARCH univariate volatility:** Assignments 1, 4, 6, 8
- **DCC-GARCH multivariate dynamics:** Assignments 1, 4, 6
- **PCA and eigenvalue decomposition:** Assignment 11
- **Bai-Perron structural break tests:** Assignment 7
- **Bootstrap inference:** Assignment 9
- **Random matrix theory cleaning:** new contribution layered on top of Assignment 4's DCC infrastructure

## Software

R, with `rugarch` for univariate GJR-GARCH (Student-t), `rmgarch` for DCC with custom targeting, base linear algebra for MP eigenvalue analysis, `strucchange` for Bai-Perron, `boot` for block bootstrap.

## Two-Week Timeline

- **Day 1:** Exploratory analysis. Pull 20-asset universe, compute rolling eigenvalue count without GARCH layer. Validate that the cocoa crisis produces a visible change in the eigenvalue count before committing further effort.
- **Days 2 to 3:** Data finalisation. Pull full universe from Refinitiv (backup: Yahoo, FRED, Investing.com). Build cleaned daily-return panel. Write data dictionary.
- **Days 4 to 5:** Univariate GJR-GARCH per asset. Standardised residuals. Diagnostic plots.
- **Days 6 to 7:** Sample correlation, MP analysis, rolling eigenvalue count. Generate Figure 1 (eigenvalue spectrum vs MP density) and Figure 2 (time series of eigenvalue count with Bai-Perron breaks).
- **Days 8 to 9:** DCC estimation with cleaned target. Pairwise conditional correlation extraction. Generate Figure 3 (cocoa-sugar and cocoa-coffee conditional correlations through the crisis).
- **Day 10:** Bootstrap confidence bands. Robustness checks (window length 252, 504; pre-crisis vs crisis subsamples; Gaussian vs Student-t innovations).
- **Day 11:** Eigenvector composition analysis. Identify asset loadings on the emerging eigenvalue.
- **Days 12 to 13:** Write. 20 pages plus appendix.
- **Day 14:** Final read-through, format, submit.

## Headline Figures

**Figure 1:** Eigenvalue spectrum of the full-sample correlation matrix, with theoretical Marchenko-Pastur density overlaid. Eigenvalues above the upper edge labelled and interpreted (market mode, ag-softs mode, currency mode, additional crisis-emerging mode).

**Figure 2:** Rolling 252-day count of supra-threshold eigenvalues, January 2015 to May 2026. Bai-Perron structural breaks marked. Bootstrap 95% confidence band. The emergence of an additional eigenvalue during the 2024 to 2026 window is the visual centrepiece of the paper.

**Figure 3:** Conditional correlation between cocoa and sugar, and cocoa and coffee, from RMT-cleaned DCC. Three panels: pre-crisis baseline (2015 to 2023), supply-rally phase (2024), demand-destruction phase (2025 to 2026). Visual support for H3.

## Premortem: Most Likely Failure Modes

1. **MP analysis shows no visible change during the crisis.** Highest-priority risk. Mitigated by day-1 exploratory analysis before committing.
2. **Largest eigenvalue dominates and washes out the spectrum.** Standard practice (remove market mode before MP fitting) handles this.
3. **GARCH instability in cocoa due to extreme returns.** Mitigated by Student-t innovations and parameter constraints.
4. **Front-month cocoa futures dominated by exchange-disruption noise during the peak crisis.** Mitigated by using rolled continuous futures and excluding flagged disruption days.
5. **Group dynamics or workload imbalance.** Mitigated by day-1 writing assignment.
6. **Data access delays.** Mitigated by maintaining a Yahoo-FRED backup throughout.

## Distinguishing Contribution

Three things together. First paper to apply RMT-cleaned DCC-GARCH (in the Engle-Ledoit-Wolf 2019 tradition) to the soft-commodity complex. First paper to identify the 2024 to 2026 cocoa crisis as a structural break in the eigenvalue spectrum of the correlation matrix. First paper to document the asymmetric supply-rally versus demand-destruction propagation pattern in the cocoa-coffee-sugar-equity system. The methodological novelty (RMT cleaning of the DCC target) is genuine but disciplined: one method, applied carefully, with robustness depth rather than estimator breadth.

---

*Proposal prepared May 2026. Submission deadline: two weeks from approval.*
