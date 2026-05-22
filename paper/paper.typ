// Cocoa-crisis RMT-DCC paper. ~2000 words.
// Compile: typst compile paper.typ
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

// --- Abstract ---
#v(0.5em)
#block(inset: (left: 0.6cm, right: 0.6cm))[
  #text(weight: "bold")[Abstract.] We test whether the 2024--2026 cocoa crisis -- a roughly $3.5×$ rise in front-month cocoa futures followed by a 75% collapse -- produced an identifiable structural break in the dependence structure of the global soft-commodity complex. On a 25-asset daily panel (2015--2026) of soft-commodity futures, adjacent agricultural futures, FX, US and Ghanaian sovereign yields, and chocolate-maker equities, we fit a GJR-GARCH$(1,1)$ per asset with Hansen-1994 skewed Student-$t$ innovations, standardise residuals, and analyse the rolling 252-day correlation matrices via Marchenko--Pastur (MP) eigenvalue cleaning and a DCC-GARCH layer with the cleaned correlation as the long-run target. The rolling supra-MP eigenvalue count (market mode excluded) rises from a pre-crisis median of 2 to 3 during 2024--2026. Bai--Perron with BIC-selected break count finds zero breaks in this series, and the pre-registered T3 condition on eigenvector composition (cocoa-bloc squared-loading share $gt.eq 0.55$ on a crisis-emerging eigenvalue) fails in every form we test. The largest (market-mode) eigenvector does reorganise: its cocoa-bloc share rises from 0.25 to 0.52, robust to bloc perturbations, but this is a finding on a different object from the supra-MP count. The MP-cleaned DCC target is rejected by the data on every criterion: the QMLE drives the long-run-target weight to 0.1% (vs 10.8% for the standard DCC), the cleaned recursion runs into the stationarity boundary ($a + b = 0.999$), the log-likelihood is 556 nats below the standard DCC, and out-of-sample minimum-variance-portfolio realised variance is identical to three decimals across methods. At $N = 25$ and $Q approx 90$, random-matrix cleaning of the DCC target does not deliver the @engle-ledoit-wolf2019 benefit.
  #v(0.3em)
  #text(weight: "bold")[Keywords:] Marchenko--Pastur; DCC-GARCH; soft commodities; cocoa; structural break; Bai--Perron.
]
#v(0.6em)

= Introduction

Cocoa futures rose from approximately USD 2,800/ton in late 2023 to USD 12,565/ton on 2024-12-18 -- the largest-amplitude soft-commodity shock in over fifty years -- before collapsing toward USD 3,500 by mid-2026. The drivers are well documented @gilbert2010: a 2023 El Niño, the spread of cacao swollen-shoot virus in West Africa, and divergent policy responses in Ghana and Côte d'Ivoire (together responsible for roughly 60% of global cocoa). We ask whether the crisis is visible as a structural break in the dependence structure of the broader soft-commodity complex, and how it reorganises that structure.

The paper makes three contributions. First, we apply the random-matrix cleaning of @engle-ledoit-wolf2019 to the DCC-GARCH of @engle2002 on a 25-asset cocoa-adjacent universe, replacing the sample-correlation DCC target with a Marchenko--Pastur-cleaned long-run target, and document that at this universe size ($N = 25$, $Q approx 90$) the cleaned target is rejected by every standard criterion. Second, we pre-register three target conditions (T1--T3) on the rolling supra-MP count and the eigenvector composition and report verdicts against them under proper structural-break methodology (BIC-selected Bai--Perron, stationary block bootstrap). Third, we show that the cocoa NY--London basis _re-couples_ during the crisis, consistent with a global supply shock affecting both contracts rather than the regional decoupling we anticipated.

The remainder of the paper describes the data (#link(<sec:data>)[§2]), methodology (#link(<sec:method>)[§3]), results against the pre-registered tests (#link(<sec:results>)[§4]), robustness (#link(<sec:robust>)[§5]), and discussion (#link(<sec:disc>)[§6]).

= Data <sec:data>

The asset universe is 25 daily series spanning 2015-01-02 to 2026-05-21 (cocoa NY trading calendar, 2,864 trading days). Sources: Refinitiv for futures, equities, FX, and the US and Ghana 10-year bid yields; FRED for the CBOE VIX. All series are pulled in their native currency to avoid the spurious co-movement that arises when prices are FX-converted to a common numeraire @pindyck-rotemberg1990.

The universe partitions into six blocks: (i) cocoa contracts on ICE-US (`CCc1`) and ICE-Europe (`LCCc1`), (ii) adjacent softs (Arabica and Robusta coffee, sugar #11, cotton, FCOJ), (iii) grains and palm oil, (iv) macro factors (DXY, VIX, US 10-year yield), (v) FX and sovereign credit (BRL, VND, Ghanaian cedi, Ghana 10-year yield), and (vi) chocolate-maker equities (Hershey, Mondelez, Nestle, JM Smucker, Starbucks, Lindt). The aspect ratio is $Q = T \/ N approx 115$ on the full sample; complete-cases (post-Ghana entry, 2017-04 onward) gives 2,260 trading days and $Q approx 90$, the sample on which the full-sample MP spectrum and the DCC are computed.

Returns are log differences for price series and first differences (percentage points) for the two yield series. WTI's 2020-04-20 negative settle (-USD 37.63) is patched to NA on -04-20 and -04-21. Ghana 10-year yield data begins on 2017-04-20 (79% of trading days); rolling windows ending before April 2018 use $N = 24$, all subsequent windows use $N = 25$. Standardised residuals are winsorised at $plus.minus 5$ to prevent any single observation from dominating a 252-day correlation matrix (the Bouchaud--Potters convention @bouchaud-potters2003; 150 of approximately 71,500 residual observations affected, 0.21%).

= Methodology <sec:method>

== Univariate volatility

The univariate filter must produce standardised residuals that are approximately white in level and square so that the rolling correlation matrices in #link(<sec:method>)[§3.2] reflect contemporaneous dependence rather than residual ARCH or autocorrelation. We use GJR-GARCH$(1,1)$ @glosten-jagannathan-runkle1993 with @hansen1994 skewed Student-$t$ (sstd) innovations as the primary specification:
$ h_(i,t) = omega_i + alpha_i epsilon_(i,t-1)^2 + gamma_i epsilon_(i,t-1)^2 #h(0.2em) bb(1)_({epsilon_(i,t-1) < 0}) + beta_i h_(i,t-1), quad z_(i,t) = epsilon_(i,t) \/ sqrt(h_(i,t)). $ <eq:gjr>

The GJR specification extends @engle1982 and @bollerslev1986 by adding an asymmetric-response term $gamma_i$ that absorbs the leverage effect on equity and FX series at the cost of one extra parameter. We prefer it to @nelson1991 EGARCH because (a) the parameters retain their familiar variance-equation interpretation, (b) the implied unconditional variance is finite under the explicit constraint $alpha_i + beta_i + gamma_i \/ 2 < 1$, and (c) numerical convergence is more robust across our heterogeneous universe; EGARCH is retained as Tier 3 in the fallback chain below. We prefer skewed-$t$ to symmetric-$t$ because cocoa returns and several FX/equity series exhibit material conditional skewness that a symmetric distribution forces into the variance equation. A likelihood-ratio test of sstd vs symmetric-$t$ (one degree of freedom) rejects symmetry at the 5% level for eight of the 25 series (COTTON, WHEAT, WTI, BRL, GHS, HSY, SJM, VIX), with VIX the most extreme ($chi^2 = 184.6$). sstd nests symmetric-$t$ (skew parameter $xi = 1$), so the choice is dominated.

The mean equation is AR(0) (constant only) except for OJ and VND, which always use AR(1); any other series whose AR(0) fit shows Ljung--Box$(10)$ p-value below 0.01 on standardised residuals is automatically refit with AR(1) on the same distribution tier (auto-retry triggered for COCOA_LDN and GHS). Residual autocorrelation in $z$ propagates directly to off-diagonal entries of the rolling correlation matrix and can manufacture supra-MP eigenvalues that do not reflect contemporaneous dependence.

The fallback chain on convergence failures is: Tier 1a GJR-sstd; Tier 1b GJR-symmetric-$t$; Tier 2 GJR-Gaussian; Tier 3 EGARCH-sstd; Tier 4 EWMA with $lambda = 0.94$. Of the 25 series, 23 fit at Tier 1a, VND at Tier 2 (Gaussian), and GHANA10Y at Tier 4 (EWMA — its yield dynamics through the 2022--23 sovereign debt restructuring defeat every parametric variant we tried). #footnote[Specification limitations we record explicitly. (i) Eight series carry a negative fitted leverage parameter ($gamma < 0$), including COCOA_NY itself ($gamma = -0.016$); on these series the GJR asymmetry term does not improve fit beyond symmetric GARCH. (ii) Five series (WHEAT, PALMOIL, VND, GHS, VIX) still fail Ljung--Box$(10)$ on standardised residuals at the 5% level after the AR(1) fix, indicating mean-equation misspecification an AR(1) does not fully absorb. (iii) Twelve series have GARCH persistence above 0.99 (near-IGARCH; full list in `output/near_igarch.csv`), so the unconditional variance is essentially undefined. (iv) Winsorisation is concentrated in eight series (108 of 150 capped observations); the headline rate of 0.21% masks per-series rates as high as 1% on GHS and GHANA10Y. None of these defeats the standardisation enough to make the downstream MP analysis unusable, but they bound the strength of any single-series claim.]

== Marchenko--Pastur cleaning of the DCC target

Let $bold(z)_t in RR^N$ denote the vector of standardised residuals at time $t$. The sample correlation matrix $overline(R)$ is eigendecomposed as $overline(R) = V Lambda V^top$. Under the null of independent unit-variance noise the bulk of the spectrum lies within the Marchenko--Pastur support
$ lambda_(plus.minus) = (1 plus.minus sqrt(N \/ T))^2 quad #text(style: "italic")[(@marchenko-pastur1967)]. $ <eq:mp>

@plerou1999 showed that for daily equity panels the empirical spectrum agrees with this prediction except in a thin set of large eigenvalues that carry the genuine factor structure, motivating eigenvalue-replacement cleaning rules for the sample correlation. We adopt the hard-clipping rule of @laloux2000: keep eigenvalues above $lambda_+$, replace those at or below by their common average. The resulting matrix $tilde(R)$ is renormalised to unit diagonal.

We prefer hard MP clipping to @ledoit-wolf2004 linear shrinkage as the long-run DCC target for two reasons. First, MP clipping is non-parametric and theoretically motivated: it separates the spectrum into a noise band whose support is determined by random-matrix theory and a signal tail that lies outside it, so the choice of cleaning is data-driven rather than a continuous bias--variance compromise. Second, the cleaned matrix is positive semi-definite by construction, full-rank after re-normalisation, and well-conditioned, properties that the DCC recursion in #link(<eq:dcc>)[Eq. (3)] needs without further regularisation.

The cleaned matrix $tilde(R)$ is then used as the long-run correlation target in the DCC recursion of @engle2002:
$ Q_t = (1 - a - b) #h(0.1em) tilde(R) + a #h(0.1em) bold(z)_(t-1) bold(z)_(t-1)^top + b #h(0.1em) Q_(t-1), quad R_t = "diag"(Q_t)^(-1/2) Q_t "diag"(Q_t)^(-1/2). $ <eq:dcc>

We use DCC rather than the constant-conditional-correlation model of @bollerslev1990 because a crisis episode is precisely the event for which time-varying correlation is the substantive object: a CCC fit would average across regimes and obscure the structural-break question the paper asks. The novelty here, following @engle-ledoit-wolf2019, is to substitute the cleaned $tilde(R)$ for the sample $overline(R)$ as the long-run anchor.

We estimate $(a, b)$ by quasi-MLE on the standardised residuals subject to $a > 0$, $b > 0$, and $a + b < 0.999$ (the hard stationarity constraint of the DCC recursion). For comparison we also estimate the standard DCC with $tilde(R)$ replaced by $overline(R)$. Both fits are sanity-checked against `rmgarch::dccfit`.

== Pre-registered targets

Three primary target conditions, *pre-registered in `docs/scope_and_target.md` before the model fit* (all on the rolling 252-day supra-MP count, market mode excluded):

- *T1 (baseline stability).* The pre-crisis (2019-01-01 to 2023-12-31, all $N = 25$ windows) count has IQR width $lt.eq 1$ and stays inside $\{1, 2, 3\}$.
- *T2 (structural break).* Two parts. (a) Bai--Perron @bai-perron2003 @bai-perron1998 with BIC-selected break count detects a break inside $[2024"-"01"-"01, 2024"-"12"-"31]$. (b) The in-crisis median count exceeds the upper 95% bound of a stationary block bootstrap @politis-romano1994 of the pre-crisis median (block mean 25, $B = 1,000$).
- *T3 (cocoa-bloc concentration).* An eigenvalue that crosses from below the MP edge in pre-crisis to above the MP edge in crisis (an "emerging" mode), with its eigenvector loading $gt.eq 0.55$ of its squared mass on the cocoa bloc $\{$COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y$\}$.

Three secondary tests run on the DCC pipeline. S1: the cleaned DCC produces lower out-of-sample minimum-variance-portfolio realised variance than the standard DCC in the crisis window (the @engle-ledoit-wolf2019 claim in its original form). S2: BIC-selected Bai--Perron detects an in-crisis break on at least one of the cocoa-coffee or cocoa-sugar pairwise cleaned-DCC correlation series. S3: the cocoa NY $times$ London conditional correlation falls by at least 0.15 from its 2019--2023 mean during the 2024--2026 crisis window.

= Results <sec:results>

== Primary targets (T1--T3)

#figure(
  image("../figures/fig2b_count_with_ci.png", width: 100%),
  caption: [Distribution of the rolling 252-day supra-MP eigenvalue count (market mode excluded), pre-crisis (2019--2023, grey) vs in-crisis (2024--2026, red). Bars show the fraction of rolling windows at each integer count value. The modal count shifts from 2 (60% of pre-crisis windows) to 3 (52% of in-crisis windows, with the remaining 48% at count = 2); the count = 1 mass present pre-crisis (10%) disappears entirely in crisis.],
) <fig:count>

@fig:count shows the empirical distribution of the rolling count under each regime. The pre-crisis sample ($N = 25$, $n = 1,257$ windows) has median 2, IQR $[2, 3]$, range $[1, 3]$, and 60% of windows equal to 2. *T1 holds.* Bai--Perron with BIC-selected break count finds *zero breaks* on the same series (BIC argmin at $m = 0$ with BIC = 931.5); no level shift survives model selection. *T2(a) fails.*

The T2(b) bootstrap distribution of the pre-crisis median is degenerate. With the test statistic an integer count, 1,000 stationary-block-bootstrap replications return the same integer median (2), giving a 95% CI of zero width, $[2.00, 2.00]$. The observed in-crisis median is 3 (52% of windows; range $[2, 3]$). T2(b) clears its formal threshold ($3 > 2$), but the procedure has no statistical resolving power on this series: it is a one-bit comparison on an integer statistic, and the bootstrap CI does not constitute a meaningful interval. Restricting to windows entirely inside the crisis (rolling-window end-date on or after 2024-12-31, $n = 349$ windows), the median is 3 and the range is $[2, 3]$.

#figure(
  image("../figures/fig3a_eigenvector_loadings.png", width: 100%),
  caption: [Eigenvector loadings of the top four supra-MP eigenvalues, pre-crisis (top row) vs in-crisis (bottom row), sign-normalised. Red bars: cocoa bloc \{COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y\}. The rank-1 in-crisis mode (the *market mode*, the largest eigenvalue) reorganises to put cocoa contracts and chocolate-maker equities jointly on it.],
) <fig:loadings>

*T3 fails on the pre-registered specification.* No eigenvalue crosses from below the MP edge in pre-crisis to above it in crisis: both subsamples carry five supra-MP eigenvalues at the same ranks (1--5). The rolling-window form of the same test (a rank whose supra-MP indicator share rises from below 0.5 pre-crisis to at least 0.5 in crisis) identifies rank 4 only, whose eigenvector has cocoa-bloc squared-loading share 0.078 — far below the 0.55 threshold. T3 *fails on both forms*.

@fig:loadings does show, as an observation outside the pre-registration, a *reorganisation of the rank-1 (market) eigenmode*. Pre-crisis its top-3 loadings are MDLZ, VIX, SBUX (cocoa-bloc share 0.251); in crisis they are MDLZ, HSY, SJM (cocoa-bloc share 0.516, a 2.1$times$ rise). Under cocoa-price stress, chocolate-maker equities and cocoa contracts load on a single risk factor — the input-cost margin-compression dynamic, adjacent to the equity--commodity co-movement channel documented in @tang-xiong2012. The result survives bloc perturbation: dropping GHANA10Y gives 0.288 $arrow.r$ 0.654, and a narrow 5-asset bloc (the two cocoa contracts plus three chocolate equities) gives 0.245 $arrow.r$ 0.515. Two qualifiers apply. First, the rank-1 mode is the market mode and is excluded by construction from the supra-MP count test; the reorganisation operates on a different statistical object. Second, a cocoa-loaded mode already existed pre-crisis at rank 3 (bloc share 0.704, top-3 loadings COCOA_NY, COCOA_LDN, CORN); in crisis that mode *weakens* to bloc share 0.310. The crisis reorganisation is not the creation of a cocoa factor but the migration of cocoa loadings onto the market mode.

== Secondary targets (S1--S3) and the cleaned-DCC fit

#figure(
  image("../figures/fig3_pairwise_correlations.png", width: 95%),
  caption: [Pairwise conditional correlations from standard (grey) and MP-cleaned (blue) DCC, for six cocoa-anchored pairs. Pink: crisis window.],
) <fig:dcc>

The two fits behave very differently. The standard DCC settles at $(a, b) = (0.005, 0.887)$ with persistence $a + b = 0.892$, placing weight $1 - a - b = 0.108$ on the long-run target $overline(R)$ at each step. The MP-cleaned DCC fits at the stationarity boundary: $(a, b) = (0.002, 0.997)$, $a + b = 0.999$, with weight on the cleaned target $tilde(R)$ of $0.001$ — the QMLE, when forced to use $tilde(R)$, drives the long-run-target weight to *0.1%*, about one hundredth of the weight the standard DCC assigns to $overline(R)$, and the optimiser does not converge to an interior solution. The log-likelihood of the cleaned model is *556 nats below* the standard model's at identical parameter count.

@fig:dcc shows the time path: the cleaned DCC (blue) produces smoother conditional correlations because $b approx 0.997$ damps innovations, but with near-zero weight on $tilde(R)$ it drifts from the long-run anchor and on most cocoa-anchored pairs lies further from the historical mean than the standard fit. Both the log-likelihood gap and the target-weight collapse point to the same conclusion: at $N = 25$ and $Q approx 90$, $tilde(R)$ is not a better long-run anchor than $overline(R)$.

For the EW-style minimum-variance-portfolio test (S1) we form $w_t = Sigma_t^(-1) bold(1) \/ bold(1)^top Sigma_t^(-1) bold(1)$ from each method, realise $r_(p, t+1) = w_t^top r_(t+1)$, and compare sample standard deviation of realised portfolio return across the 2019--2023 baseline and the 2024--2026 crisis window. Because the cleaned model has effectively discarded $tilde(R)$, its conditional-covariance path coincides with the standard DCC's: *the OOS-MVP realised standard deviation is identical to three decimals in both windows* (baseline $1.561%$ vs $1.561%$; crisis $1.019%$ vs $1.019%$). *S1 fails.*

For S2, BIC-selected Bai--Perron on each of the eleven cocoa-anchored cleaned-DCC pairwise correlation series ($m_max = 5$, $h = 0.15$) returns *zero breaks for every pair*. *S2 fails.*

For S3, the baseline 2019--2023 mean of the cocoa NY × London basis correlation under the cleaned DCC is 0.542; the 2024--2026 crisis-window mean is 0.627. The basis *rises* by 0.085 rather than falling by the predicted 0.15. The crisis minimum is 0.528, 0.014 below the baseline mean. *S3 fails in the predicted direction.* The interpretation is in #link(<sec:disc>)[§6].

= Robustness <sec:robust>

#figure(
  image("../figures/fig4_robustness.png", width: 75%),
  caption: [Rolling supra-MP count under five robustness variants. The +1 in-crisis lift in the median survives the W=252 variants; under W=504 both pre and crisis windows saturate at 3 (no within-variant lift).],
) <fig:robust>

@fig:robust runs the rolling-count test under four variants alongside the baseline: (i) excluding GHANA10Y (addressing the universe-size-change concern); (ii) a 504-day window; (iii) (i)+(ii) combined; (iv) GARCH refit with Gaussian innovations. The +1 median lift on T2(b) holds only under the baseline $W = 252$ specification and the Gaussian-innovations variant. Removing GHANA10Y at $W = 252$ collapses the in-crisis median from 3 to 2, so the lift is materially dependent on that series. Under $W = 504$, both pre- and in-crisis medians saturate at 3 in both universe sizes, so the bootstrap cannot distinguish them. BIC-selected Bai--Perron returns *zero breaks* for every variant. The +1 lift survives the innovation-distribution perturbation but not the universe-size perturbation or the longer window; no break-detection variant supports T2(a).

= Discussion <sec:disc>

Four findings warrant discussion.

*The +1 count lift is real but fragile, and the pre-registered T2 tests do not formally confirm it.* The rolling supra-MP count rises from a pre-crisis median of 2 to an in-crisis median of 3 on the baseline universe at $W = 252$, and the lift survives a switch to Gaussian innovations. It does not survive removing GHANA10Y (in-crisis median drops to 2) or extending the window to $W = 504$ (both medians saturate at 3). T2(a) under BIC-selected Bai--Perron returns zero breaks under every robustness variant; T2(b) clears its threshold $(3 > 2)$ but is a one-bit comparison whose bootstrap CI is degenerate. The substantive shift in the count distribution is therefore visible (#link(<fig:count>)[Figure 1]) but not validated as a structural break by a procedure with statistical resolving power on this data.

*The rank-1 reorganisation is real but operates on a different statistical object than the count test.* The market-mode eigenvector's cocoa-bloc squared-loading share rises from 0.245 pre-crisis to 0.517 in crisis (#link(<fig:loadings>)[Figure 3]), survives bloc-set perturbations, and matches the economic mechanism of input-cost margin compression at chocolate makers. T3 as pre-registered asked for this concentration on a crisis-emerging eigenvalue, with threshold 0.55. No emerging eigenvalue exists, and the rank-1 share falls 0.033 below 0.55, so the pre-registered T3 fails. The rank-1 mode is moreover excluded from the supra-MP count by construction, so this finding does not bear on T2.

*The cleaned DCC is rejected by the data on every criterion the @engle-ledoit-wolf2019 framework supplies.* The 556-nat log-likelihood gap to the standard DCC at identical parameter count is decisive. The mechanism is visible in the fitted parameters: when forced to use $tilde(R)$, the QMLE drives the long-run-target weight from approximately 11% to 0.1%, almost as low as the constraint $a + b < 0.999$ allows — the model is, by its own optimisation, downweighting the cleaned target. The OOS minimum-variance-portfolio realised variance, the EW2019 advertised benefit, is identical to three decimal places across methods because the cleaned model has effectively become the standard one. EW2019 demonstrated their result for $N gt.eq 100$ with much smaller aspect ratios; at our $N = 25$ and $Q approx 90$ the raw $overline(R)$ is already estimated with low enough sampling noise that the hard-clipping step compresses informative off-diagonals rather than de-noising them.

*The geographic-decoupling hypothesis on the NY--London basis is rejected in favour of a global supply-shock interpretation.* The basis correlation rose by 0.085 in the crisis, not the pre-registered $gt.eq 0.15$ fall. NY's deliverable basket is West-African-dominant and London's is mixed-origin, so we expected divergence under a regional shock. The observed re-coupling implies that the 2024--2026 shock transmitted to both contracts — NY directly via its dominant deliverable, London indirectly via substitution effects in the deliverable basket. The arbitrage-windfall narrative is incompatible with the data; a fundamental supply-side reading is consistent with it.

*Specification limitations beyond what the results expose.* Equity time-zone asynchrony (Swiss and US closes are six hours apart) depresses LISN--HSY-class daily correlations; weekly aggregation is left for future work. The univariate fits carry negative leverage on COCOA_NY and seven other series, and five series fail Ljung--Box on standardised residuals at the 5% level after the AR(1) fix (most severely GHS and VND); both inject residual heteroskedasticity into the rolling correlation matrix. GHANA10Y is fit by EWMA fallback rather than parametric GARCH because no parametric variant converges through the 2022--23 Ghanaian sovereign-debt restructuring. These limitations bound the strength of single-series claims but do not change the verdicts on T1--T3 or S1--S3.

= Conclusion

The 2024--2026 cocoa crisis produces three visible signatures in the soft-commodity correlation spectrum: a +1 lift in the rolling 252-day supra-MP eigenvalue count, a reorganisation of the market-mode eigenvector to load on cocoa contracts and chocolate-maker equities jointly, and a re-coupling rather than decoupling of the NY--London basis correlation. None of T2(a), T2(b) in its statistically meaningful sense, T3, S1, S2, or S3 passes as pre-registered, and the +1 lift itself does not survive removing GHANA10Y from the universe. The structural-change pattern is real on the baseline configuration but is not formally validated by any of the break-detection procedures we ran.

The geographic-decoupling hypothesis on the cocoa basis is rejected in favour of a global supply-shock interpretation: the shock transmits to both ICE-US and ICE-Europe contracts. The Marchenko--Pastur cleaning of the DCC long-run target, as proposed by @engle-ledoit-wolf2019, is rejected by the data at $N = 25$, $Q approx 90$ on every standard criterion (log-likelihood, target-weight, OOS minimum-variance-portfolio variance), placing a clear lower bound on the universe size at which random-matrix cleaning of the DCC anchor becomes empirically defensible. The descriptive characterisation of the cocoa shock and the negative finding on EW2019 cleaning at small $N$ are the paper's substantive and methodological contributions, respectively.

#bibliography("refs.yml", style: "ieee", title: "References")

#pagebreak()
#counter(heading).update(0)
#set heading(numbering: "A.")

= Appendix: Implementation notes

The full pipeline (R code, intermediate panels, fit diagnostics, bootstrap output) is reproducible from `/home/bcip/numtech`. Frozen artefacts:
- `data/clean/panel_v1.csv` (levels, 2864 × 25)
- `data/clean/returns_v1.csv` (returns, 2863 × 25)
- `data/clean/std_residuals_v1.csv` (winsorised GJR-GARCH-$t$ standardised residuals)
- `output/garch_fits.csv` (per-asset GARCH parameters and Ljung--Box diagnostics)
- `output/rolling_count_v1.csv` (rolling supra-MP counts)
- `output/bootstrap_t2b.csv`, `output/bootstrap_rolling_ci.csv` (T2(b) inputs)
- `output/eigenvector_loadings_subsamples.csv`, `output/t3_test.csv` (T3 inputs)
- `output/dcc_params.csv`, `output/dcc_pairwise_correlations.csv` (DCC parameters and conditional correlation paths)
- `output/s1_mvp_oos_variance.csv` (S1, EW-style OOS MVP test), `output/s1_corr_path_variance_legacy.csv` (an alternative correlation-path-variance metric, retained for reference)
- `output/s2_bai_perron_pairs.csv` (S2, BIC-selected breaks per pair), `output/s3_cocoa_basis_decoupling.csv` (S3, full 2024--2026 crisis window)
- `output/t3_test.csv` (four-row table: original spec, rolling spec, rank-1 reorganisation at 0.55, rank-1 reorganisation at 0.50)
- `output/bai_perron_bic_table.csv` (T2(a) BIC table)
- `output/robustness_summary.csv` (Figure 4 numbers)
- `output/test_verdicts.csv` (consolidated pass/fail table for T1, T2(a), T2(b), T3(a-c), S1, S2, S3)
- `output/return_skewness.csv` (empirical skew + excess-kurtosis per series, motivating sstd)
- `output/near_igarch.csv` (twelve series with persistence $gt.eq 0.99$ flagged for stationarity violation)
- `R/diagnostics.R`, `R/diagnostics_extra.R`, `R/write_verdicts.R` (end-to-end numerical diagnostics; re-run to verify every figure quoted in this paper)

#figure(
  image("../figures/fig1_mp_spectrum.png", width: 80%),
  caption: [Full-sample eigenvalue spectrum vs Marchenko--Pastur density. Six eigenvalues exceed the upper edge $lambda_+ = 1.221$.],
) <fig:spectrum>
