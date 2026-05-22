// Cocoa-crisis RMT-DCC paper.
// Compile: typst compile paper.typ --root /home/bcip/numtech
// Reads figures from ../figures/ (project root is /home/bcip/numtech/paper)

#set page(paper: "a4", margin: (x: 2.2cm, y: 2.2cm), numbering: "1")
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true, leading: 0.6em)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show heading.where(level: 1): it => block(above: 1.4em, below: 0.7em, text(weight: "bold", size: 12pt, it))
#show heading.where(level: 2): it => block(above: 1.0em, below: 0.5em, text(weight: "bold", size: 11pt, it))

// --- Title block ---
#align(center)[
  #text(15pt, weight: "bold")[Identifying the 2024--2026 Cocoa Crisis as a Structural Break in Soft-Commodity Correlation Dynamics]
  #v(0.4em)
  #text(11pt, style: "italic")[A Random-Matrix-Cleaned DCC-GARCH Approach]
  #v(0.6em)
  FNCE40003 Numerical Techniques in Finance #h(1.5em) Major Assignment #h(1.5em) May 2026
]

#v(0.8em)

= Introduction

Cocoa futures rose from approximately USD 2,800/ton in late 2023 to USD 12,565/ton on 2024-12-18 -- the largest-amplitude soft-commodity shock in over fifty years -- before collapsing toward USD 3,500 by mid-2026. The drivers are well documented @gilbert2010: a 2023 El Niño, the spread of cacao swollen-shoot virus in West Africa, and divergent policy responses in Ghana and Côte d'Ivoire (together responsible for roughly 60% of global cocoa). We ask whether the crisis is visible as a structural break in the dependence structure of the broader soft-commodity complex, and how it reorganises that structure.

The paper makes three contributions. First, we apply the random-matrix cleaning of @engle-ledoit-wolf2019 to the DCC-GARCH of @engle2002 on a 25-asset cocoa-adjacent universe, replacing the sample-correlation DCC target with a Marchenko--Pastur-cleaned long-run target, and document that at this universe size ($N = 25$, $Q approx 90$) the cleaned target is rejected by every standard criterion. Second, we pre-register tests on the rolling supra-MP count and on the eigenvector composition, and report verdicts under BIC-selected Bai--Perron and a stationary block bootstrap. Third, we show that the cocoa NY--London basis _re-couples_ during the crisis, consistent with a global supply shock affecting both contracts rather than the regional decoupling we anticipated.

= Data <sec:data>

#figure(
  table(
    columns: (auto, 1fr),
    align: (left + horizon, left + horizon),
    inset: (x: 8pt, y: 5pt),
    stroke: none,
    table.hline(stroke: 0.8pt),
    table.header([*Block*], [*Series*]),
    table.hline(stroke: 0.4pt),
    [Cocoa contracts (2)], [COCOA_NY (ICE-US), COCOA_LDN (ICE-Europe)],
    [Adjacent softs (5)], [Coffee Arabica, Coffee Robusta, Sugar No. 11, Cotton, FCOJ],
    [Grains, oilseeds, energy (5)], [Corn, Soybeans, Wheat, Palm oil, WTI crude],
    [Macro factors (3)], [DXY, VIX, US 10Y yield],
    [FX and sovereign credit (4)], [BRL, VND, GHS, Ghana 10Y yield],
    [Chocolate-maker equities (6)], [HSY, MDLZ, NESN, SJM, SBUX, LISN],
    table.hline(stroke: 0.8pt),
  ),
  caption: [The 25-asset cocoa-adjacent universe. Daily series, 2015-01-02 to 2026-05-21 (2,864 trading days on the cocoa NY calendar). All series from Refinitiv except VIX (FRED). Series held in native currency to avoid spurious FX-conversion co-movement @pindyck-rotemberg1990. Returns are log differences for prices and first differences (percentage points) for the two yield series.],
) <tab:universe>

The aspect ratio is $Q = T \/ N approx 115$ on the full sample; the complete-cases panel (post-Ghana entry, 2017-04 onward) gives 2,260 trading days and $Q approx 90$, the sample on which the full-sample MP spectrum and the DCC are computed. WTI's 2020-04-20 negative settle (-USD 37.63) is patched to NA on -04-20 and -04-21. Ghana 10-year yield data begins on 2017-04-20 (79% of trading days); rolling windows ending before April 2018 use $N = 24$, all subsequent windows use $N = 25$. Standardised residuals are winsorised at $plus.minus 5$ (the Bouchaud--Potters convention @bouchaud-potters2003; 150 of approximately 71,500 residual observations affected, 0.21%).

= Methodology <sec:method>

== Univariate volatility

The univariate filter must produce standardised residuals that are approximately white in level and square so that the rolling correlation matrices in #link(<sec:method>)[§3.2] reflect contemporaneous dependence. We use GJR-GARCH$(1,1)$ @glosten-jagannathan-runkle1993, an asymmetric extension of @engle1982 and @bollerslev1986, with @hansen1994 skewed Student-$t$ (sstd) innovations:
$ h_(i,t) = omega_i + alpha_i epsilon_(i,t-1)^2 + gamma_i epsilon_(i,t-1)^2 #h(0.2em) bb(1)_({epsilon_(i,t-1) < 0}) + beta_i h_(i,t-1), quad z_(i,t) = epsilon_(i,t) \/ sqrt(h_(i,t)). $ <eq:gjr>

GJR is preferred to @nelson1991 EGARCH for the interpretability of the asymmetric term $gamma_i$ and a finite unconditional variance under the constraint $alpha_i + beta_i + gamma_i \/ 2 < 1$; EGARCH is retained as Tier 3 in the fallback chain. The sstd nests symmetric-$t$ at $xi = 1$, and a likelihood-ratio test rejects symmetry at the 5% level for 8 of the 25 series, so the symmetric specification is dominated.

The mean equation is AR(0) (constant) except for OJ and VND, which always use AR(1); any series whose AR(0) fit shows Ljung--Box$(10)$ p-value below 0.01 is automatically refit at AR(1) (triggered for COCOA_LDN and GHS). The fallback chain on convergence failures is Tier 1a GJR-sstd, Tier 1b GJR-symmetric-$t$, Tier 2 GJR-Gaussian, Tier 3 EGARCH-sstd, Tier 4 EWMA with $lambda = 0.94$. Twenty-three series fit at Tier 1a, VND at Tier 2, and GHANA10Y at Tier 4 -- its yield dynamics through the 2022--23 Ghanaian sovereign-debt restructuring defeat every parametric variant.#footnote[Eight series carry negative leverage ($gamma < 0$, including COCOA_NY); five series (WHEAT, PALMOIL, VND, GHS, VIX) still fail Ljung--Box$(10)$ after the AR(1) fix; twelve are near-IGARCH (persistence $gt.eq 0.99$, full list in `output/near_igarch.csv`). These bound the strength of any single-series claim but do not defeat the downstream analysis.]

== Marchenko--Pastur cleaning of the DCC target

Let $bold(z)_t in RR^N$ denote the vector of standardised residuals at time $t$, with sample correlation $overline(R) = V Lambda V^top$. Under independent unit-variance noise the bulk of the spectrum lies within the Marchenko--Pastur support
$ lambda_(plus.minus) = (1 plus.minus sqrt(N \/ T))^2 quad #text(style: "italic")[(@marchenko-pastur1967)]. $ <eq:mp>

@plerou1999 document that for daily financial panels the empirical spectrum agrees with this prediction except in a thin set of large eigenvalues carrying the genuine factor structure. We adopt the hard-clipping rule of @laloux2000: keep eigenvalues above $lambda_+$, replace those at or below by their common average, and renormalise to unit diagonal. The cleaned matrix $tilde(R)$ is positive semi-definite by construction, full-rank, and well-conditioned; we prefer it to the linear shrinkage of @ledoit-wolf2004 because the noise/signal split is theoretically motivated rather than a continuous bias--variance compromise.

We then use $tilde(R)$ as the long-run target in the DCC recursion of @engle2002 (in contrast to the constant-conditional-correlation model of @bollerslev1990):
$ Q_t = (1 - a - b) #h(0.1em) tilde(R) + a #h(0.1em) bold(z)_(t-1) bold(z)_(t-1)^top + b #h(0.1em) Q_(t-1), quad R_t = "diag"(Q_t)^(-1/2) Q_t "diag"(Q_t)^(-1/2). $ <eq:dcc>

The cleaned-target substitution follows @engle-ledoit-wolf2019. We estimate $(a, b)$ by quasi-MLE subject to $a, b > 0$ and $a + b < 0.999$ (the hard stationarity constraint of the DCC recursion); for comparison we also fit the standard DCC with $tilde(R)$ replaced by $overline(R)$.

== Pre-registered tests

Three primary tests, pre-registered in `docs/scope_and_target.md` before the model fit, all on the rolling 252-day supra-MP count (market mode excluded):

- *Baseline stability.* The pre-crisis (2019-01-01 to 2023-12-31, all $N = 25$ windows) count has IQR width $lt.eq 1$ and stays inside $\{1, 2, 3\}$.
- *Structural break.* (a) Bai--Perron @bai-perron2003 @bai-perron1998 with BIC-selected break count detects a break inside $[2024"-"01"-"01, 2024"-"12"-"31]$. (b) The in-crisis median count exceeds the upper 95% bound of a stationary block bootstrap @politis-romano1994 of the pre-crisis median (block mean 25, $B = 1,000$).
- *Cocoa-bloc concentration.* An eigenvalue that crosses from below the MP edge in pre-crisis to above it in crisis, with eigenvector loading $gt.eq 0.55$ of its squared mass on the cocoa bloc $\{$COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y$\}$.

Three secondary tests on the DCC pipeline. *MVP variance:* the cleaned DCC produces lower out-of-sample minimum-variance-portfolio realised variance than the standard DCC during the crisis. *Pairwise breaks:* BIC-selected Bai--Perron detects an in-crisis break on at least one cocoa-coffee or cocoa-sugar pairwise cleaned-DCC correlation series. *NY--London basis:* the cocoa NY $times$ London conditional correlation falls by at least 0.15 from its 2019--2023 mean during 2024--2026.

= Results <sec:results>

== Rolling supra-MP count and eigenvector structure

#figure(
  image("../figures/fig2b_count_with_ci.png", width: 100%),
  caption: [Distribution of the rolling 252-day supra-MP eigenvalue count (market mode excluded), pre-crisis (2019--2023, grey) vs in-crisis (2024--2026, red). The modal count shifts from 2 (60% of pre-crisis windows) to 3 (52% of in-crisis windows); the count = 1 mass present pre-crisis (10%) disappears entirely in crisis.],
) <fig:count>

@fig:count shows the empirical distribution of the rolling count under each regime. The pre-crisis sample ($N = 25$, $n = 1,257$ windows) has median 2, IQR $[2, 3]$, range $[1, 3]$, and 60% of windows equal to 2; baseline stability is satisfied. Bai--Perron with BIC-selected break count finds *zero breaks* on the same series (BIC argmin at $m = 0$, BIC = 931.5). The bootstrap form is degenerate: with an integer test statistic, 1,000 stationary-block-bootstrap replications return the same integer median (2), giving a 95% CI of zero width $[2.00, 2.00]$. The observed in-crisis median is 3, so the procedure clears its formal threshold ($3 > 2$), but it is a one-bit comparison with no statistical resolving power.

#figure(
  image("../figures/fig3a_eigenvector_loadings.png", width: 100%),
  caption: [Eigenvector loadings of the top four supra-MP eigenvalues, pre-crisis (top row) vs in-crisis (bottom row), sign-normalised. Red bars: cocoa bloc \{COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y\}. The rank-1 in-crisis mode (the market mode) reorganises to put cocoa contracts and chocolate-maker equities jointly on it.],
) <fig:loadings>

The pre-registered cocoa-bloc concentration test fails. No eigenvalue crosses from below the MP edge in pre-crisis to above it in crisis: both subsamples carry five supra-MP eigenvalues at the same ranks. A rolling-window variant identifies rank 4 only, whose eigenvector has cocoa-bloc squared-loading share 0.078 -- far below the 0.55 threshold.

@fig:loadings shows, outside the pre-registration, a reorganisation of the rank-1 (market) eigenmode. Its top-3 loadings shift from MDLZ, VIX, SBUX pre-crisis to MDLZ, HSY, SJM in crisis; the cocoa-bloc squared-loading share rises from 0.251 to 0.516, a 2.1$times$ increase. The result survives bloc perturbation (dropping GHANA10Y: $0.288 arrow.r 0.654$; narrow 5-asset bloc of the two cocoa contracts plus three chocolate equities: $0.245 arrow.r 0.515$). Under cocoa-price stress, chocolate-maker equities and cocoa contracts load on a single risk factor -- input-cost margin compression at chocolate makers, adjacent to the equity--commodity co-movement channel of @tang-xiong2012. Two qualifiers: the market mode is excluded by construction from the supra-MP count, so this finding is independent of the count test; and a cocoa-loaded mode already existed pre-crisis at rank 3 (bloc share 0.704), so the crisis reorganisation is the migration of cocoa loadings onto the market mode, not the creation of a new cocoa factor.

== DCC fit and the cleaned-target benchmark

#figure(
  image("../figures/fig3_pairwise_correlations.png", width: 95%),
  caption: [Pairwise conditional correlations from standard (grey) and MP-cleaned (blue) DCC, for six cocoa-anchored pairs. Pink: crisis window.],
) <fig:dcc>

The two fits behave very differently. The standard DCC settles at $(a, b) = (0.005, 0.887)$, persistence $a + b = 0.892$, weight $1 - a - b = 0.108$ on $overline(R)$. The MP-cleaned DCC fits at the stationarity boundary: $(a, b) = (0.002, 0.997)$, $a + b = 0.999$, weight on $tilde(R)$ of $0.001$ -- the QMLE drives the long-run-target weight to *0.1%*, about one hundredth of what the standard DCC places on $overline(R)$, and the optimiser does not converge to an interior solution. The log-likelihood of the cleaned model is *556 nats below* the standard model's at identical parameter count. @fig:dcc shows the time path: the cleaned DCC produces smoother conditional correlations because $b approx 0.997$ damps innovations, but with near-zero weight on $tilde(R)$ it drifts from the long-run anchor and on most cocoa-anchored pairs lies further from the historical mean than the standard fit.

The minimum-variance-portfolio test forms $w_t = Sigma_t^(-1) bold(1) \/ bold(1)^top Sigma_t^(-1) bold(1)$ from each method and compares sample standard deviations of realised portfolio return across 2019--2023 and 2024--2026. Because the cleaned model has effectively discarded $tilde(R)$, its conditional-covariance path coincides with the standard DCC's: the out-of-sample MVP realised standard deviation is identical to three decimals in both windows (baseline $1.561%$ vs $1.561%$; crisis $1.019%$ vs $1.019%$). BIC-selected Bai--Perron on the eleven cocoa-anchored cleaned-DCC pairwise correlation series returns *zero breaks for every pair*. The basis-decoupling test runs in the opposite direction to the prediction: the baseline 2019--2023 mean of the cocoa NY $times$ London correlation is 0.542 and the 2024--2026 mean is 0.627, so the basis *rises* by 0.085 rather than falling by the predicted 0.15.

= Robustness <sec:robust>

#figure(
  image("../figures/fig4_robustness.png", width: 75%),
  caption: [Rolling supra-MP count under five robustness variants. The +1 in-crisis lift survives the $W=252$ variants; under $W=504$ both pre and crisis windows saturate at 3 (no within-variant lift).],
) <fig:robust>

@fig:robust runs the rolling-count test under four variants alongside the baseline: excluding GHANA10Y, a 504-day window, both combined, and GARCH refit with Gaussian innovations. The +1 in-crisis lift holds only under the baseline $W = 252$ specification and the Gaussian-innovations variant. Removing GHANA10Y at $W = 252$ collapses the in-crisis median from 3 to 2, so the lift is materially dependent on that series; under $W = 504$, both pre- and in-crisis medians saturate at 3, so the bootstrap cannot distinguish them. BIC-selected Bai--Perron returns *zero breaks* for every variant.

= Discussion <sec:disc>

The rolling supra-MP count rises from a pre-crisis median of 2 to an in-crisis median of 3 on the baseline universe at $W = 252$, and the lift survives a switch to Gaussian innovations, but it does not survive removing GHANA10Y or extending the window. BIC-selected Bai--Perron returns zero breaks under every variant, and the bootstrap form of the break test, while clearing its threshold, has no resolving power on an integer statistic with a degenerate CI. The substantive shift in the distribution is therefore visible (@fig:count) but not validated as a structural break. The cleanest signal in the spectrum is on a different object: the rank-1 (market-mode) eigenvector reorganises to load on cocoa contracts and chocolate-maker equities jointly, with its cocoa-bloc squared-loading share rising from 0.245 to 0.517 and surviving bloc-set perturbations (@fig:loadings). The mechanism is input-cost margin compression at chocolate makers, the equity--commodity co-movement channel of @tang-xiong2012. Because the market mode is excluded by construction from the supra-MP count, this finding is independent of the count-shift result; because the pre-registered concentration test required an emerging eigenvalue (none exists) and a 0.55 threshold (the rank-1 share falls 0.033 short), it fails on both pre-registered conditions while remaining the most economically interpretable result the paper produces.

The cleaned DCC is rejected by the data on every criterion the @engle-ledoit-wolf2019 framework supplies. The 556-nat log-likelihood gap to the standard DCC at identical parameter count is decisive; the QMLE drives the long-run-target weight from approximately 11% to 0.1%, almost as low as the stationarity constraint $a + b < 0.999$ allows. The minimum-variance-portfolio realised variance is therefore identical to three decimal places across methods, because the cleaned model has by its own optimisation become the standard one. @engle-ledoit-wolf2019 demonstrated their result for $N gt.eq 100$ with much smaller aspect ratios; at our $N = 25$ and $Q approx 90$ the raw $overline(R)$ is already estimated with low enough sampling noise that hard MP clipping compresses informative off-diagonals rather than de-noising them.

The geographic-decoupling hypothesis on the cocoa NY--London basis is rejected in favour of a global supply-shock interpretation. The basis correlation rises by 0.085 in the crisis rather than falling by the pre-registered 0.15. NY's deliverable basket is West-African-dominant and London's is mixed-origin, so we expected divergence under a regional shock and observed re-coupling instead. The 2024--2026 shock transmits to both contracts -- NY directly via its dominant deliverable, London indirectly via substitution effects in the deliverable basket -- which is incompatible with an arbitrage-windfall narrative and consistent with a fundamental supply-side reading.

= Conclusion

The 2024--2026 cocoa crisis produces three signatures in the soft-commodity correlation spectrum: a +1 lift in the rolling 252-day supra-MP count, a reorganisation of the market-mode eigenvector onto cocoa contracts and chocolate-maker equities, and a re-coupling rather than decoupling of the NY--London basis. None passes its pre-registered test as a formally validated structural break, and the count lift itself does not survive removing GHANA10Y or extending the window; the pattern is real on the baseline configuration but unconfirmed by procedures with statistical resolving power on this data. Separately, the Marchenko--Pastur cleaning of the DCC long-run target is rejected at $N = 25$, $Q approx 90$ on every criterion the @engle-ledoit-wolf2019 framework supplies, placing a clear lower bound on the universe size at which random-matrix cleaning of the DCC anchor becomes empirically defensible. The descriptive characterisation of the cocoa shock and the negative finding on cleaning at small $N$ are the paper's substantive and methodological contributions.

#bibliography("refs.yml", style: "ieee", title: "References")
