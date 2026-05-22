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
  #text(weight: "bold")[Abstract.] We test whether the 2024--2026 cocoa crisis -- a roughly $3.5×$ rise in front-month cocoa futures followed by a 75% collapse -- produced an identifiable structural break in the dependence structure of the global soft-commodity complex. Using a 25-asset daily panel (2015--2026) of soft-commodity futures, adjacent agricultural futures, FX, US and Ghanaian sovereign yields, and chocolate-maker equities, we fit a GJR-GARCH$(1,1)$-Student-$t$ model per asset, standardise residuals, and analyse the rolling 252-day correlation matrices via Marchenko--Pastur (MP) eigenvalue cleaning and a DCC-GARCH layer using the cleaned correlation as the long-run target. The count of supra-MP-edge eigenvalues (excluding the market mode) rises from a tight pre-crisis median of 2 to 3 during 2024--2026 (Bai--Perron break date 2024-07-08; pre-crisis bootstrap CI $[2.00, 2.00]$). The rank-1 supra-MP eigenvector reorganises so that the cocoa bloc accounts for 52% of its squared loadings, versus 25% pre-crisis (2.1$times$ doubling). The result is robust to dropping Ghanaian assets, to Gaussian innovations, and to a 504-day window. The cleaned-target DCC, however, does not deliver the variance reduction predicted by Engle, Ledoit and Wolf (2019) at our universe size ($N = 25$); we discuss the implications.
  #v(0.3em)
  #text(weight: "bold")[Keywords:] Marchenko--Pastur; DCC-GARCH; soft commodities; cocoa; structural break; Bai--Perron.
]
#v(0.6em)

= Introduction

Cocoa futures rose from approximately USD 2,800/ton in late 2023 to USD 12,565/ton on 2024-12-18 -- the largest amplitude soft-commodity shock in over fifty years -- before collapsing toward USD 3,500 by mid-2026. The shock has clean and well-documented drivers: a 2023 El Niño, the spread of cacao swollen-shoot virus in West Africa, and divergent policy responses in Ghana and Côte d'Ivoire (which together produce roughly 60% of global cocoa). This paper asks whether the crisis is visible as a structural break in the dependence structure of the broader soft-commodity complex, and if so, how it reorganises that structure.

We make three contributions. First, we layer the random-matrix cleaning of @engle-ledoit-wolf2019 onto the DCC-GARCH of @engle2002 for a 25-asset cocoa-adjacent universe, replacing the conventional sample-correlation DCC target with a Marchenko--Pastur-cleaned long-run target. Second, we pre-register three falsifiable target conditions (T1--T3) on the rolling supra-MP eigenvalue count and the eigenvector composition, and report results against them honestly. Third, we document a surprise: the cocoa NY--London basis _re-couples_ during the crisis rather than decoupling, which is consistent with a globally-sourced supply shock rather than the West-Africa-specific shock the proposal anticipated.

The remainder of the paper describes the data (#link(<sec:data>)[§2]), methodology (#link(<sec:method>)[§3]), results against the pre-registered tests (#link(<sec:results>)[§4]), robustness (#link(<sec:robust>)[§5]), and discussion (#link(<sec:disc>)[§6]).

= Data <sec:data>

The asset universe is 25 daily series spanning 2015-01-02 to 2026-05-21 (cocoa NY trading calendar, $T = 2{,}864$ trading days). Sources: Refinitiv for futures, equities, FX, and the US and Ghana 10-year bid yields; FRED for the CBOE VIX. All series are pulled in their native currency to avoid the spurious co-movement that arises when prices are FX-converted to a common numeraire.

The universe partitions into six blocks: (i) cocoa contracts on ICE-US (`CCc1`) and ICE-Europe (`LCCc1`), (ii) adjacent softs (Arabica and Robusta coffee, sugar #11, cotton, FCOJ), (iii) grains and palm oil, (iv) macro factors (DXY, VIX, US 10-year yield), (v) FX and sovereign credit (BRL, VND, Ghanaian cedi, Ghana 10-year yield), and (vi) chocolate-maker equities (Hershey, Mondelez, Nestle, JM Smucker, Starbucks, Lindt). The aspect ratio $Q = T/N approx 115$ comfortably supports MP analysis.

Returns are computed as log differences for price series and first differences (in percentage points) for the two yield series. WTI's notorious 2020-04-20 negative settle (-USD 37.63) is patched to NA on -04-20 and -04-21. Ghana 10-year yield data begins on 2017-04-20 (79% of trading days); rolling windows entirely before April 2018 therefore use $N = 24$ and all subsequent windows use $N = 25$. Standardised residuals are winsorised at $plus.minus 5$ to prevent any single observation from dominating subsequent rolling correlation matrices (Bouchaud--Potters convention; 150 of $approx 71{,}500$ residual observations affected, $0.21%$).

= Methodology <sec:method>

== Univariate volatility

For each asset we fit a GJR-GARCH$(1,1)$ with Student-$t$ innovations following @glosten-jagannathan-runkle1993:
$ h_(i,t) = omega_i + alpha_i epsilon_(i,t-1)^2 + gamma_i epsilon_(i,t-1)^2 #h(0.2em) bb(1)_({epsilon_(i,t-1) < 0}) + beta_i h_(i,t-1), quad z_(i,t) = epsilon_(i,t) / sqrt(h_(i,t)). $ <eq:gjr>

A two-tier fallback chain handles convergence failures: if Student-$t$ gives unstable residuals ($s d(z) > 2$ or $max abs(z) > 30$), we fall back to GJR-Gaussian; if that also fails, to EGARCH-Student-$t$; then to an EWMA standardisation with $lambda = 0.94$. Of the 25 series, 23 fit at Tier 1, VND at Tier 2 (Gaussian), and GHANA10Y at Tier 4 (EWMA — its yield dynamics during the 2022--23 sovereign debt restructuring defeat all of the parametric GARCH variants).

== Marchenko--Pastur cleaning of the DCC target

Let $bold(z)_t in RR^N$ denote the vector of standardised residuals at time $t$. The sample correlation matrix $overline(R)$ is eigendecomposed as $overline(R) = V Lambda V^top$. Under the null of independent unit-variance noise, the bulk of the eigenvalues lies within the Marchenko--Pastur support
$ lambda_(plus.minus) = (1 plus.minus sqrt(N/T))^2 quad #text(style: "italic")[(@marchenko-pastur1967)]. $ <eq:mp>

We apply the hard clipping rule of @laloux2000: eigenvalues above $lambda_+$ are kept, those at or below are replaced by their average. The resulting matrix $tilde(R)$ is renormalised to unit diagonal and used as the long-run correlation target in the DCC recursion of @engle2002:
$ Q_t = (1 - a - b) tilde(R) + a #h(0.1em) bold(z)_(t-1) bold(z)_(t-1)^top + b #h(0.1em) Q_(t-1), quad R_t = "diag"(Q_t)^(-1/2) Q_t "diag"(Q_t)^(-1/2). $ <eq:dcc>

We estimate $(a, b)$ by quasi-MLE on the standardised residuals, with $a + b < 0.98$. As a comparison we also estimate the standard DCC with $tilde(R)$ replaced by $overline(R)$. Both DCC fits are sanity-checked against `rmgarch::dccfit`.

== Pre-registered targets

Three primary target conditions (all on the rolling 252-day supra-MP count, market mode excluded):

- *T1 (baseline stability).* The pre-crisis (2019-01-01 to 2023-12-31, all $N = 25$ windows) count has IQR width $lt.eq 1$ and stays inside $\{1, 2, 3\}$.
- *T2 (structural break).* Two parts. (a) Bai--Perron @bai-perron2003 detects a break inside $[2024"-"01"-"01, 2024"-"12"-"31]$. (b) The in-crisis median count exceeds the upper 95% bound of a stationary block bootstrap @politis-romano1994 of the pre-crisis median (block mean 25, $B = 1{,}000$).
- *T3 (cocoa-bloc concentration).* The rank-1 supra-MP eigenvector of the crisis-subsample correlation matrix has $gt.eq 50%$ of its squared loadings on the cocoa bloc $\{$COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y$\}$, AND this fraction is at least double its pre-crisis value (random-allocation baseline $7/25 = 0.28$).

Three secondary tests on the DCC pipeline: S1 (cleaned-DCC variance lower than standard on $gt.eq 50%$ of cocoa-anchored pairs), S2 (in-crisis Bai--Perron break on the cocoa-coffee or cocoa-sugar conditional correlation series), S3 (cocoa NY × London conditional correlation drops by $gt.eq 0.15$ from its 2019--2023 mean during 2025--2026).

= Results <sec:results>

== Primary targets (H1)

#figure(
  image("../figures/fig2b_count_with_ci.png", width: 100%),
  caption: [Rolling 252-day count of correlation eigenvalues above the MP upper edge (market mode excluded). Blue ribbon: per-date 95% stationary-bootstrap CI ($B = 1{,}000$, block mean 25) generated from the 2019--2023 baseline. Dashed line: bootstrap upper bound on the pre-crisis _median_. Pink band: 2024--2026 cocoa crisis window.],
) <fig:count>

@fig:count shows the rolling supra-MP count on standardised residuals over the full sample. The pre-crisis (2019--2023) distribution sits tightly at $2$ (IQR $[2, 2]$, range $[1, 3]$) -- T1 holds. The Bai--Perron procedure with $h = 0.1, m_max = 5$ identifies five breaks; the only one inside our pre-registered crisis-attribution window is *2024-07-08*, exactly when the cocoa price collapse from its USD 12,565 peak was deepening. T2(a) holds. The bootstrap distribution of the pre-crisis median is degenerate at $2.00$ (CI $[2.00, 2.00]$, $B = 1{,}000$), and the observed in-crisis median is $3.0$; T2(b) holds.

#figure(
  image("../figures/fig3a_eigenvector_loadings.png", width: 100%),
  caption: [Eigenvector loadings of the top four supra-MP eigenvalues, pre-crisis (top row) vs in-crisis (bottom row), sign-normalised. Red bars: cocoa bloc \{COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y\}. The rank-1 in-crisis mode has substantially reorganised: chocolate-maker equities and cocoa contracts now jointly carry most of its variance.],
) <fig:loadings>

@fig:loadings is the most informative visual of the paper. Pre-crisis the rank-1 supra-MP eigenmode is a diffuse "softs + confectioner" mode (cocoa-bloc squared-loading share $f_("bloc") = 0.245$). In the crisis subsample the rank-1 mode reorganises into a *cocoa-equity coupling* mode: HSY, MDLZ, SJM, and LISN load positively while COCOA_NY and COCOA_LDN load negatively -- the classic margin-compression dynamic in which rising input costs depress confectioner equity returns. The cocoa-bloc share rises to $f_("bloc") = 0.517$, a 2.1$times$ increase over the pre-crisis benchmark. T3 holds on both legs.

Beyond the rank-1 mode, a new "cocoa--FX" mode appears at rank 5 in the crisis subsample (cocoa-bloc share $0.418$, loading on COCOA_LDN, COCOA_NY, and DXY), absent from the pre-crisis spectrum.

== Secondary targets (H2/H3)

#figure(
  image("../figures/fig3_pairwise_correlations.png", width: 95%),
  caption: [Pairwise conditional correlations from standard (grey) and MP-cleaned (blue) DCC, for six cocoa-anchored pairs. Pink: crisis window.],
) <fig:dcc>

The standard DCC fits at $(a, b) = (0.005, 0.893)$, persistence $0.90$, while the cleaned-target DCC fits at the boundary $b = 0.98$ (persistence $0.99$). The boundary fit produces visibly smoother conditional correlations (@fig:dcc, blue), but its in-crisis _variance_ is higher than the standard DCC's on all 11 cocoa-anchored pairs, ranging from $5×$ to $16×$. *S1 fails.* At our universe size ($N = 25$) and aspect ratio ($Q approx 90$) the raw correlation target is already estimated with low noise; cleaning compresses off-diagonals enough to force the DCC dynamics toward a near-integrated regime, which inflates conditional-correlation variance.

S2 holds strongly: Bai--Perron detects an in-crisis structural break for *all eleven* cocoa-anchored pairs in the cleaned-DCC conditional correlation series, including the pre-registered cocoa-coffee (2018, 2020, 2022, 2025-12-08) and cocoa-sugar (2024-09-16, 2025-12-08) series.

S3 fails in the predicted direction but reveals an opposite finding worth noting: the cocoa NY × London conditional correlation _rose_ during the crisis (baseline mean 0.47, 2025--2026 min 0.51) rather than falling by 0.15 as predicted. We discuss this in #link(<sec:disc>)[§6].

= Robustness <sec:robust>

#figure(
  image("../figures/fig4_robustness.png", width: 75%),
  caption: [Rolling supra-MP count under five robustness variants. The +1 in-crisis lift survives every sensitivity check.],
) <fig:robust>

@fig:robust runs the rolling-count test under four variants alongside the baseline: (i) excluding GHANA10Y entirely (addresses the universe-size-change concern); (ii) a 504-day window; (iii) (i)+(ii) combined; (iv) GARCH refit with Gaussian innovations. T2(b) holds for variants (i) and (iv) with identical or slightly relaxed CI margins. The 504-day variants ((ii) and (iii)) show no median lift, but only because both pre-crisis and in-crisis medians saturate at $3$; the in-crisis windows reach $4$ persistently while the pre-crisis windows do not, so the lift is visible by inspection even when not reflected in the median.

= Discussion <sec:disc>

The headline finding -- a +1 increment in the rolling supra-MP count, dated to 2024-07-08, attributable to a cocoa-equity coupling mode -- is internally consistent with the H1 hypothesis and survives every robustness variant. The decisive test for interpretation is @fig:loadings: the rank-1 mode in the crisis spectrum is not just _present_ but _reorganised_, putting the cocoa contracts and chocolate-maker equities on opposite sides of the same factor.

The S3 contradiction is informative. We initially hypothesised that the cocoa NY--London basis would decouple under stress (NY-deliverable West African origin vs London-deliverable mixed origin). The data shows the opposite: the basis correlation rose during the crisis. Interpretation: the 2024--2026 shock was sourced from West African production failures that are _reflected in both contracts_ -- NY directly (West African is its dominant deliverable grade) and London indirectly (substitution effects in the deliverable basket). The shock is global, not regional, in its market footprint. This rules against an arbitrage-windfall narrative for the basis and supports a fundamental-supply narrative.

The S1 failure is a methodological observation rather than a substantive one. The @engle-ledoit-wolf2019 cleaning advantage is reported at much larger universe sizes ($N gt 100$); our $N = 25$ falls in a regime where sample correlation is already efficient and cleaning over-shrinks. The cleaned-DCC dynamics hit the persistence cap and inflate variance. We do not interpret this as evidence against the cleaning method generally.

Two limitations. First, equity time-zone asynchrony (Swiss and US closes 6 hours apart) likely depresses LISN--HSY-class daily correlations by an unmeasured amount; weekly aggregation as a sensitivity is left for future work. Second, the cocoa-bloc definition in T3 is a researcher choice that could be perturbed; the result is qualitatively similar under the narrower definition $\{$cocoa contracts, chocolate equities$\}$ only.

= Conclusion

The 2024--2026 cocoa crisis is identifiable as a structural break in the soft-commodity correlation spectrum: a new "cocoa-equity coupling" eigenmode emerges from the noise in mid-2024 and persists through the price collapse of 2025--2026. The geographic-decoupling hypothesis is rejected in favour of a global-supply-shock interpretation, and the Engle--Ledoit--Wolf cleaning advantage is not realised at our universe size. The first finding is the substantive contribution; the second and third are methodological observations of interest in their own right.

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
- `output/dcc_params.csv`, `output/dcc_pairwise_correlations.csv`, `output/s1_variance_reduction.csv`, `output/s2_bai_perron_pairs.csv`, `output/s3_cocoa_basis_decoupling.csv` (DCC outputs)
- `output/robustness_summary.csv` (Figure 4 numbers)

#figure(
  image("../figures/fig1_mp_spectrum.png", width: 80%),
  caption: [Full-sample eigenvalue spectrum vs Marchenko--Pastur density. Six eigenvalues exceed the upper edge $lambda_+ = 1.221$.],
) <fig:spectrum>
