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
  #text(weight: "bold")[Abstract.] We test whether the 2024--2026 cocoa crisis -- a roughly $3.5×$ rise in front-month cocoa futures followed by a 75% collapse -- produced an identifiable structural break in the dependence structure of the global soft-commodity complex. Using a 25-asset daily panel (2015--2026) of soft-commodity futures, adjacent agricultural futures, FX, US and Ghanaian sovereign yields, and chocolate-maker equities, we fit a GJR-GARCH$(1,1)$ per asset with Hansen-1994 skewed Student-$t$ innovations (AR(1) mean for series with material residual autocorrelation), standardise residuals, and analyse the rolling 252-day correlation matrices via Marchenko--Pastur (MP) eigenvalue cleaning and a DCC-GARCH layer using the cleaned correlation as the long-run target. The rolling count of supra-MP-edge eigenvalues (market mode excluded) rises from a pre-crisis median of 2 to 3 during 2024--2026, but two of three pre-registered structural-break tests fail under proper specification: Bai--Perron with BIC-selected break count detects *zero* breaks in the rolling count, and none of the three pre-registered T3 (eigenvector composition) specifications passes at its original 0.55 cocoa-bloc-loading threshold. An exploratory observation -- that the largest (market-mode) eigenvector reorganises during the crisis so that the cocoa bloc carries 52% of its squared loadings versus 25% pre-crisis -- is robust to bloc definition but is not the pre-registered specification and was reported as a pre-registered pass only against a threshold (0.50) silently moved from the pre-registered 0.55. The MP-cleaned DCC target is rejected by the data on every standard criterion: the QMLE drives the long-run-target weight to 0.1% (vs 10.8% for the standard DCC), the cleaned recursion is near-integrated ($a + b = 0.999$, fails to converge to an interior solution), the log-likelihood is 556 nats lower than the standard DCC's, and out-of-sample minimum-variance-portfolio realised variance is identical to three decimals across methods. We read this as evidence that random-matrix cleaning of the DCC target does not deliver the expected benefit at our universe size ($N = 25$).
  #v(0.3em)
  #text(weight: "bold")[Keywords:] Marchenko--Pastur; DCC-GARCH; soft commodities; cocoa; structural break; Bai--Perron.
]
#v(0.6em)

= Introduction

Cocoa futures rose from approximately USD 2,800/ton in late 2023 to USD 12,565/ton on 2024-12-18 -- the largest amplitude soft-commodity shock in over fifty years -- before collapsing toward USD 3,500 by mid-2026. The shock has clean and well-documented drivers: a 2023 El Niño, the spread of cacao swollen-shoot virus in West Africa, and divergent policy responses in Ghana and Côte d'Ivoire (which together produce roughly 60% of global cocoa). This paper asks whether the crisis is visible as a structural break in the dependence structure of the broader soft-commodity complex, and if so, how it reorganises that structure.

We make three contributions. First, we apply the random-matrix cleaning of @engle-ledoit-wolf2019 to the DCC-GARCH of @engle2002 for a 25-asset cocoa-adjacent universe, replacing the conventional sample-correlation DCC target with a Marchenko--Pastur-cleaned long-run target -- and find that at this universe size the cleaned target is rejected by the data. Second, we pre-register three falsifiable target conditions (T1--T3) on the rolling supra-MP eigenvalue count and the eigenvector composition, and report results against them honestly, including the two that fail and one that passes only under a threshold not in the pre-registration. Third, we document a surprise: the cocoa NY--London basis _re-couples_ during the crisis rather than decoupling, which is consistent with a globally-sourced supply shock rather than the West-Africa-specific shock the proposal anticipated.

The remainder of the paper describes the data (#link(<sec:data>)[§2]), methodology (#link(<sec:method>)[§3]), results against the pre-registered tests (#link(<sec:results>)[§4]), robustness (#link(<sec:robust>)[§5]), and discussion (#link(<sec:disc>)[§6]).

= Data <sec:data>

The asset universe is 25 daily series spanning 2015-01-02 to 2026-05-21 (cocoa NY trading calendar, $T = 2{,}864$ trading days). Sources: Refinitiv for futures, equities, FX, and the US and Ghana 10-year bid yields; FRED for the CBOE VIX. All series are pulled in their native currency to avoid the spurious co-movement that arises when prices are FX-converted to a common numeraire.

The universe partitions into six blocks: (i) cocoa contracts on ICE-US (`CCc1`) and ICE-Europe (`LCCc1`), (ii) adjacent softs (Arabica and Robusta coffee, sugar #11, cotton, FCOJ), (iii) grains and palm oil, (iv) macro factors (DXY, VIX, US 10-year yield), (v) FX and sovereign credit (BRL, VND, Ghanaian cedi, Ghana 10-year yield), and (vi) chocolate-maker equities (Hershey, Mondelez, Nestle, JM Smucker, Starbucks, Lindt). The aspect ratio for the full sample is $Q = T/N approx 115$; complete-cases (post-Ghana entry, 2017-04 onward) gives $T = 2{,}260$ and $Q approx 90$, which is the sample on which the full-sample MP spectrum and DCC are computed.

Returns are computed as log differences for price series and first differences (in percentage points) for the two yield series. WTI's notorious 2020-04-20 negative settle (-USD 37.63) is patched to NA on -04-20 and -04-21. Ghana 10-year yield data begins on 2017-04-20 (79% of trading days); rolling windows entirely before April 2018 therefore use $N = 24$ and all subsequent windows use $N = 25$. Standardised residuals are winsorised at $plus.minus 5$ to prevent any single observation from dominating subsequent rolling correlation matrices (Bouchaud--Potters convention; 150 of $approx 71{,}500$ residual observations affected, $0.21%$).

= Methodology <sec:method>

== Univariate volatility

For each asset we fit a GJR-GARCH$(1,1)$ following @glosten-jagannathan-runkle1993, with @hansen1994 skewed Student-$t$ (sstd) innovations as the primary specification:
$ h_(i,t) = omega_i + alpha_i epsilon_(i,t-1)^2 + gamma_i epsilon_(i,t-1)^2 #h(0.2em) bb(1)_({epsilon_(i,t-1) < 0}) + beta_i h_(i,t-1), quad z_(i,t) = epsilon_(i,t) / sqrt(h_(i,t)). $ <eq:gjr>

The mean equation is AR(0) (constant only) except for OJ and VND, which always use AR(1), and any other series whose AR(0) fit shows Ljung--Box$(10)$ p-value $< 0.01$ on standardised residuals, which is automatically retried with AR(1) on the same distribution tier (this auto-retry triggered for COCOA_LDN and GHS).

The fallback chain handles convergence failures: Tier 1a GJR-sstd; Tier 1b GJR-symmetric-$t$; Tier 2 GJR-Gaussian; Tier 3 EGARCH-sstd; Tier 4 EWMA with $lambda = 0.94$. Of the 25 series, 23 fit at Tier 1a, VND at Tier 2 (Gaussian), and GHANA10Y at Tier 4 (EWMA — its yield dynamics through the 2022--23 sovereign debt restructuring defeat all parametric variants). A likelihood-ratio test of sstd vs symmetric-$t$ (one degree of freedom) is significant at the 5% level for eight series (COTTON, WHEAT, WTI, BRL, GHS, HSY, SJM, VIX), with VIX the most extreme ($chi^2 = 184.6$). #footnote[Diagnostic limitations of this stage are non-trivial and we record them here. (i) Eight series have a negative fitted leverage parameter ($gamma < 0$), including COCOA_NY itself ($gamma = -0.016$); on these series GJR is no improvement over symmetric GARCH and the GJR motivation is formally inapplicable. (ii) Five series (WHEAT, PALMOIL, VND, GHS, VIX) still fail Ljung--Box$(10)$ on the standardised level series at the 5% level after the AR(1) fix, indicating residual mean-equation misspecification an AR(1) cannot absorb. (iii) Twelve series have GARCH persistence above 0.99 (near-IGARCH; see `output/near_igarch.csv`), so the unconditional variance is essentially undefined. (iv) Winsorisation is concentrated in eight series (108 of 150 capped observations); the headline rate of 0.21% masks per-series rates as high as 1% on GHS and GHANA10Y. None of these defeat the standardisation enough to make the downstream MP analysis unusable, but they constrain the strength of conclusions one can draw.]

== Marchenko--Pastur cleaning of the DCC target

Let $bold(z)_t in RR^N$ denote the vector of standardised residuals at time $t$. The sample correlation matrix $overline(R)$ is eigendecomposed as $overline(R) = V Lambda V^top$. Under the null of independent unit-variance noise, the bulk of the eigenvalues lies within the Marchenko--Pastur support
$ lambda_(plus.minus) = (1 plus.minus sqrt(N/T))^2 quad #text(style: "italic")[(@marchenko-pastur1967)]. $ <eq:mp>

We apply the hard clipping rule of @laloux2000: eigenvalues above $lambda_+$ are kept, those at or below are replaced by their average. The resulting matrix $tilde(R)$ is renormalised to unit diagonal and used as the long-run correlation target in the DCC recursion of @engle2002:
$ Q_t = (1 - a - b) tilde(R) + a #h(0.1em) bold(z)_(t-1) bold(z)_(t-1)^top + b #h(0.1em) Q_(t-1), quad R_t = "diag"(Q_t)^(-1/2) Q_t "diag"(Q_t)^(-1/2). $ <eq:dcc>

We estimate $(a, b)$ by quasi-MLE on the standardised residuals, with $a > 0$, $b > 0$, and $a + b < 0.999$ (the hard stationarity constraint of the DCC recursion). As a comparison we also estimate the standard DCC with $tilde(R)$ replaced by $overline(R)$. Both DCC fits are sanity-checked against `rmgarch::dccfit`.

== Pre-registered targets

Three primary target conditions, *pre-registered in `docs/scope_and_target.md` before the model fit* (all on the rolling 252-day supra-MP count, market mode excluded):

- *T1 (baseline stability).* The pre-crisis (2019-01-01 to 2023-12-31, all $N = 25$ windows) count has IQR width $lt.eq 1$ and stays inside $\{1, 2, 3\}$.
- *T2 (structural break).* Two parts. (a) Bai--Perron @bai-perron2003 with BIC-selected break count detects a break inside $[2024"-"01"-"01, 2024"-"12"-"31]$. (b) The in-crisis median count exceeds the upper 95% bound of a stationary block bootstrap @politis-romano1994 of the pre-crisis median (block mean 25, $B = 1{,}000$).
- *T3 (cocoa-bloc concentration, pre-registered form).* An eigenvalue that crosses from below the MP edge in pre-crisis to above the MP edge in crisis (an "emerging" mode), with its eigenvector loading $gt.eq 0.55$ of its squared mass on the cocoa bloc $\{$COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y$\}$.#footnote[We disclose that an earlier draft of this paper rewrote T3 as a rank-1 reorganisation test with threshold $gt.eq 0.50$. That is a moved goalpost and is not the pre-registered specification. We report the pre-registered T3 verdict here and the rank-1 reorganisation as an exploratory observation in #link(<sec:results>)[§4.1].]

Three secondary tests on the DCC pipeline: S1 (cleaned-DCC out-of-sample minimum-variance-portfolio realised variance lower than standard in the crisis window, i.e., the @engle-ledoit-wolf2019 claim in its original form), S2 (in-crisis Bai--Perron break, BIC-selected, on the cocoa-coffee or cocoa-sugar conditional correlation series), S3 (cocoa NY × London conditional correlation drops by $gt.eq 0.15$ from its 2019--2023 mean during the full crisis window 2024--2026).

= Results <sec:results>

== Primary targets (T1--T3)

#figure(
  image("../figures/fig2b_count_with_ci.png", width: 100%),
  caption: [Rolling 252-day count of correlation eigenvalues above the MP upper edge (market mode excluded). Blue ribbon: per-date 95% stationary-bootstrap CI ($B = 1{,}000$, block mean 25) generated from the 2019--2023 baseline. Dashed line: bootstrap upper bound on the pre-crisis _median_. Pink band: 2024--2026 cocoa crisis window.],
) <fig:count>

@fig:count shows the rolling supra-MP count on standardised residuals over the full sample. The pre-crisis (2019--2023, $N=25$) distribution has median 2, IQR $[2, 3]$ (width 1), and range $[1, 3]$, with 57% of windows equal to 2 — *T1 holds*. The Bai--Perron procedure with BIC-selected break count detects *zero breaks*: the rolling count has no statistically significant level shift on a within-class model-selection criterion. *T2(a) fails.* When `breaks = 5` was forced in an earlier specification, the algorithm placed five breaks at 2019-04-16, 2020-04-30, 2021-05-13, 2023-02-09 and 2024-07-08, with only the last inside the pre-registered crisis-attribution window; under proper specification this attribution does not survive.

The T2(b) bootstrap distribution of the pre-crisis median is degenerate: 1{,}000 stationary-block-bootstrap replications all give the same integer median (2), so the 95% CI has zero width $[2.00, 2.00]$. The observed in-crisis median is 3. *T2(b) holds in the literal sense* that 3 > 2, but the test reduces to a one-bit comparison on integer counts; the $B = 1{,}000$ framing implies statistical resolving power that on this series does not exist. Restricting to *fully* in-crisis windows (rolling-window end-date on or after 2024-12-31, $n = 349$ windows), the median saturates at 3 and the range is $[2, 3]$, marginally strengthening the verdict.

#figure(
  image("../figures/fig3a_eigenvector_loadings.png", width: 100%),
  caption: [Eigenvector loadings of the top four supra-MP eigenvalues, pre-crisis (top row) vs in-crisis (bottom row), sign-normalised. Red bars: cocoa bloc \{COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y\}. The rank-1 in-crisis mode (the *market mode*, the largest eigenvalue) reorganises to put cocoa contracts and chocolate-maker equities jointly on it.],
) <fig:loadings>

*T3 fails on the pre-registered specification.* There is no eigenvalue that crosses from below the MP edge in pre-crisis to above it in crisis (the "emerging mode"): both subsamples contain five supra-MP eigenvalues at the same ranks (1--5). The rolling-window form of the same test (rank with supra-MP indicator share $gt.eq 0.5$ in crisis but $< 0.5$ pre-crisis) identifies rank 4 only, whose eigenvector has cocoa-bloc share 0.075 -- far below the 0.55 threshold. The pre-registered T3 verdict on both forms is *FAIL*.

What @fig:loadings does show, as an exploratory observation not in the pre-registration, is a *reorganisation of the rank-1 (market) eigenmode*: pre-crisis its top-3 loadings are MDLZ, VIX, SBUX (cocoa-bloc share 0.251); in crisis they become MDLZ, HSY, SJM (cocoa-bloc share 0.516, a 2.1$times$ increase). The interpretation is intuitive — under cocoa-price stress, chocolate-maker equities and cocoa contracts move together as a single risk factor, the classic input-cost margin-compression dynamic. This finding survives perturbation of the bloc set: dropping GHANA10Y from the bloc gives 0.288 $arrow.r$ 0.654, and restricting to a narrow 5-asset bloc (cocoa contracts and three chocolate equities only) gives 0.245 $arrow.r$ 0.515. The finding is real but is *not* what was pre-registered, and it is on the market mode (which the rolling-count test excludes); the rank-3 cocoa-specific factor that existed *pre-crisis* (bloc share 0.704, top loadings COCOA_NY, COCOA_LDN, CORN) in fact *weakens* in crisis to bloc share 0.310.

== Secondary targets (S1--S3) and the cleaned-DCC fit

#figure(
  image("../figures/fig3_pairwise_correlations.png", width: 95%),
  caption: [Pairwise conditional correlations from standard (grey) and MP-cleaned (blue) DCC, for six cocoa-anchored pairs. Pink: crisis window.],
) <fig:dcc>

The two DCC fits behave very differently. The standard DCC fits at $(a, b) = (0.005, 0.887)$ with persistence $a + b = 0.892$, putting weight $1 - a - b = 0.108$ on the long-run target $overline(R)$ at each step. The MP-cleaned DCC fits at the cap: $(a, b) = (0.002, 0.997)$, $a + b = 0.999$ (the stationarity boundary), weight on the cleaned target $tilde(R)$ of just $0.001$ -- the QMLE, when forced to use the cleaned target, drives the weight on that target down to *0.1%*, around one-hundred times lower than the standard DCC, and the optimiser fails to converge to an interior solution. The log-likelihood of the cleaned model is *556 nats lower* than the standard model's, at identical parameter count. Neither test of model fit favours the cleaned target.

@fig:dcc shows the time path: the cleaned DCC (blue) produces smoother conditional correlations (because $b approx 0.98$ damps innovations), but it can drift considerably from the long-run anchor, which on most cocoa-anchored pairs places its in-crisis values further from the historic mean than the standard model's. We treat the *log-likelihood and target-weight numbers as the headline result* on the cleaning question -- by the data's own model-selection logic, $tilde(R)$ is not an improvement on $overline(R)$ at $N = 25$.

For the EW-style minimum-variance-portfolio test (S1) we form $w_t = Sigma_t^(-1) bold(1) \/ bold(1)^top Sigma_t^(-1) bold(1)$ from each method, realise $r_(p, t+1) = w_t^top r_(t+1)$, and compare sample standard deviation of the realised portfolio return across the 2019--2023 baseline and the 2024--2026 crisis window. The cleaned DCC has driven its long-run-target weight so close to zero that its conditional-covariance path coincides with the standard DCC's; *the OOS-MVP realised standard deviation is identical to three decimals in both windows* (baseline $1.561%$ vs $1.561%$; crisis $1.019%$ vs $1.019%$). *S1 fails* (no improvement; the cleaned target has been functionally discarded by the optimiser).

For S2, the cleaned-DCC pairwise conditional correlation series are re-tested with BIC-selected break count: *none of the eleven cocoa-anchored pairs has any break* at $m_max = 5$, $h = 0.15$. The pre-registered cocoa-coffee and cocoa-sugar predictions fail. The earlier draft of this paper reported "all eleven pairs have an in-crisis break", which was an artefact of forcing five breaks per pair — every pair was *mechanically required* to receive five breaks irrespective of fit improvement. *S2 fails.*

For S3, on the cocoa NY × London basis under the cleaned DCC, the baseline 2019--2023 mean is 0.542 and the *crisis (2024--2026, full window) mean is 0.627* — a *rise* of 0.085, not a drop. The crisis minimum is 0.528, only 0.014 below the baseline mean. *S3 fails in the predicted direction*; the result is in fact the opposite of the prediction. We discuss the interpretation in #link(<sec:disc>)[§6].

= Robustness <sec:robust>

#figure(
  image("../figures/fig4_robustness.png", width: 75%),
  caption: [Rolling supra-MP count under five robustness variants. The +1 in-crisis lift in the median survives the W=252 variants; under W=504 both pre and crisis windows saturate at 3 (no within-variant lift).],
) <fig:robust>

@fig:robust runs the rolling-count test under four variants alongside the baseline: (i) excluding GHANA10Y entirely (addresses the universe-size-change concern); (ii) a 504-day window; (iii) (i)+(ii) combined; (iv) GARCH refit with Gaussian innovations. The +1 median lift on T2(b) holds under the baseline W=252 specification and the Gaussian-innovations variant only. Variant (i) — removing GHANA10Y — collapses the in-crisis median from 3 to 2, so T2(b) fails: the +1 lift is materially dependent on the GHANA10Y series. Under W=504 the medians saturate at 3 in both pre- and in-crisis windows for both universe sizes, so the bootstrap test cannot distinguish them, and BIC-selected Bai--Perron returns *zero* breaks for every variant. *Honest summary: the +1 lift survives the GARCH-distribution perturbation but does not survive removing GHANA10Y, nor the longer rolling window; BIC-selected Bai--Perron returns zero breaks under every variant we tried.*

= Discussion <sec:disc>

Three findings warrant discussion.

*The +1 count lift is real but fragile and does not pass either pre-registered T2 test on its terms.* The rolling 252-day supra-MP count rises from a median of 2 in the pre-crisis baseline to 3 in the crisis. The lift survives the Gaussian-innovations robustness variant at W=252, but it does *not* survive removing GHANA10Y (in-crisis median drops to 2) nor extending the window to W=504 (both pre- and in-crisis medians saturate at 3). However, neither component of T2 *as pre-registered* survives: T2(a) under BIC-selected Bai--Perron returns zero breaks; T2(b) is a one-bit comparison on an integer count, with the bootstrap distribution literally degenerate. The honest reading is that the +1 lift is consistent with a structural break but is not statistically *confirmed* as one. The 2024-07-08 break date that earlier circulated was an artefact of forcing five breaks; it does not survive BIC selection.

*The rank-1 reorganisation is an exploratory finding, not a pre-registered confirmation, and it concerns the market mode.* @fig:loadings makes the reorganisation visually unambiguous, and the cocoa-bloc share rises from 0.245 to 0.517 on rank 1 (2.1$times$). It is robust to the bloc definition. But this is *not* the pre-registered T3 spec, which asked for an emerging eigenvalue with bloc share $gt.eq 0.55$; that test fails in all three forms we ran (subsample, rolling, and rank-1 reorganisation at the original 0.55 threshold). The 0.517 share that an earlier draft reported as a pass against a 0.50 threshold reflects a moved goalpost — at the pre-registered 0.55 the result fails by 0.033. Furthermore, the rolling count test excludes the market mode by construction, so the reorganisation is on a different object from the count statistic that motivates the paper.

*The cleaned DCC is rejected by the data on every standard criterion, including the very criterion EW2019 advertises.* The 556-nat log-likelihood gap to the standard DCC, with identical parameter count, is decisive. The qualitative reason is plain in the fitted parameters: when forced to use the cleaned target, the QMLE drives the long-run-target weight from $approx 11%$ down to $0.1%$ — the model is, as much as the constraints allow, *ignoring* the cleaned target. This is the opposite of what theory predicts if cleaning produces a better target. We read this as evidence that at $N = 25$ and $Q approx 90$ the raw sample correlation is already estimated with low enough noise that aggressive eigenvalue-replacement compresses informative off-diagonals; the cleaned target is then a worse anchor than the noisy raw one. The EW result was demonstrated at $N gt.eq 100$ with much smaller $Q$; our regime is far from that asymptotic comfort zone.

*S3 contradicts its prediction in a usefully informative direction.* The 2019--2023 baseline mean of the NY--London basis correlation (under the cleaned DCC) is 0.542; the full crisis-window mean is 0.627 -- a *rise* of 0.085, not the predicted $gt.eq 0.15$ fall. We hypothesised regional decoupling because NY's deliverable basket is West-African-dominant while London's is mixed-origin. The observed re-coupling implies that the 2024--2026 shock was felt in *both* contracts -- NY directly via its dominant deliverable, London indirectly via substitution in the deliverable basket. The market footprint of the shock is global, not regional, which rules against an arbitrage-windfall narrative and supports a fundamental-supply narrative.

*Specification limitations beyond what the results expose.* Equity time-zone asynchrony (Swiss and US closes are 6 hours apart) likely depresses LISN--HSY-class daily correlations; weekly aggregation as a sensitivity is left for future work. The univariate GARCH fits show negative leverage on COCOA_NY and several other series ($gamma < 0$), and five series (WHEAT, PALMOIL, VND, GHS, VIX) still fail Ljung--Box on standardised residuals after the AR(1) fix (most severely GHS and VND); both of these inject heteroskedasticity into the rolling correlation matrix that the methodology does not fully absorb. The GHANA10Y series is fit by EWMA fallback rather than parametric GARCH because no parametric variant converges through the 2022--23 Ghanaian sovereign debt restructuring. These are limitations on the strength of all numerical claims; we have not been able to design within-scope sensitivities that fully neutralise them.

= Conclusion

The 2024--2026 cocoa crisis is consistent with a structural change in the soft-commodity correlation spectrum -- the rolling supra-MP count rises from 2 to 3, the market-mode eigenvector reorganises onto a cocoa--chocolate-equity coupling, and the NY--London basis correlation re-couples rather than decouples -- but the pre-registered structural-break tests (T2(a) Bai--Perron with BIC selection, T3 cocoa-bloc concentration at the original 0.55 threshold) do not pass on their original terms, and the +1 count lift itself does not survive removing GHANA10Y from the universe. The +1 lift in the rolling count and the rank-1 reorganisation are real findings on the baseline universe; the framing of either as a *confirmed* structural break by formal test is not. The geographic-decoupling hypothesis is rejected in favour of a global-supply-shock interpretation. The MP cleaning of the DCC long-run target is rejected by the data at $N = 25$ on standard model-selection criteria. The substantive contribution is the descriptive characterisation of the cocoa-shock dependence structure; the methodological contribution is the negative finding on EW2019 cleaning at this universe size.

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
- `output/s1_mvp_oos_variance.csv` (S1, EW-style OOS MVP test), `output/s1_corr_path_variance_legacy.csv` (the earlier mis-framed corr-path-variance metric, retained for transparency)
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
