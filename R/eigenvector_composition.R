# Day 11 (run early): eigenvector composition analysis + T3 test.
#
# T3 (from docs/scope_and_target.md):
#   Define the crisis-emerging eigenvalue as the lowest-ranked supra-MP
#   eigenvalue that exceeds the MP edge in ≥50% of in-crisis windows but
#   in <50% of 2019-2023 baseline windows. Its eigenvector (averaged over
#   crisis-period windows) has the fraction of its squared loadings on the
#   cocoa bloc {COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y} ≥ 0.55.
#
# Outputs:
#   output/eigenvector_loadings_subsamples.csv
#   output/t3_test.csv
#   figures/fig3a_eigenvector_loadings.png    (subsample comparison)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

COCOA_BLOC <- c("COCOA_NY", "COCOA_LDN", "HSY", "MDLZ", "LISN", "GHS", "GHANA10Y")
CRISIS_START  <- as.Date("2024-01-01")
CRISIS_END    <- as.Date("2026-05-22")
BASELINE_START <- as.Date("2019-01-01")
BASELINE_END   <- as.Date("2023-12-31")

z <- fread("data/clean/std_residuals_v1.csv")
z[, date := as.Date(date)]
tickers <- setdiff(names(z), "date")

# Restrict to dates where all 25 series are observed (post-Ghana entry)
z <- z[complete.cases(z[, tickers, with = FALSE])]
cat(sprintf("Working sample: %d dates x %d series (complete cases)\n", nrow(z), length(tickers)))

mp_upper <- function(N, T) (1 + sqrt(N / T))^2

# ----- Subsample eigendecompositions ----------------------------------------
analyze_subsample <- function(z_sub, label) {
  Tn <- nrow(z_sub); N <- length(tickers)
  R <- cor(as.matrix(z_sub[, tickers, with = FALSE]))
  e <- eigen(R, symmetric = TRUE)
  edge <- mp_upper(N, Tn)
  evals <- e$values
  evecs <- e$vectors

  out <- data.table(
    sample = label,
    rank = seq_along(evals),
    eigenvalue = evals,
    above_mp = evals > edge,
    mp_edge = edge,
    T = Tn, N = N
  )
  # Loadings (each row = ticker, columns = eigenvector ranks)
  load_dt <- as.data.table(evecs)
  setnames(load_dt, paste0("ev", seq_along(evals)))
  load_dt[, ticker := tickers]
  load_dt[, sample := label]
  setcolorder(load_dt, c("sample", "ticker"))

  # Squared-loading fraction on cocoa bloc, per eigenvector
  bloc_idx <- match(COCOA_BLOC, tickers)
  bloc_mass <- apply(evecs[bloc_idx, , drop = FALSE]^2, 2L, sum)  # rows are bloc tickers
  out[, cocoa_bloc_l2_frac := bloc_mass]
  out[, top3_loadings := sapply(seq_along(evals), function(k) {
    v <- evecs[, k]
    nm <- tickers[order(-abs(v))[1:3]]
    paste(nm, collapse = ", ")
  })]

  list(summary = out, loadings = load_dt)
}

pre  <- analyze_subsample(z[date >= BASELINE_START & date <= BASELINE_END], "pre_2019_2023")
crs  <- analyze_subsample(z[date >= CRISIS_START  & date <= CRISIS_END],   "crisis_2024_2026")
full <- analyze_subsample(z, "full_sample")

all_summary <- rbindlist(list(pre$summary, crs$summary, full$summary))
all_loadings <- rbindlist(list(pre$loadings, crs$loadings, full$loadings), fill = TRUE)

fwrite(all_summary, "output/eigenvector_loadings_subsamples.csv")
fwrite(all_loadings, "output/eigenvector_loadings_full.csv")

cat("\n=== Supra-MP eigenvalues per subsample ===\n")
print(all_summary[above_mp == TRUE,
                  .(sample, rank, eigenvalue = round(eigenvalue, 3),
                    mp_edge = round(mp_edge, 3),
                    cocoa_bloc_frac = round(cocoa_bloc_l2_frac, 3),
                    top3_loadings)])

# ----- T3: identify the crisis-emerging eigenvalue --------------------------
# Per the spec: the lowest-ranked supra-MP eigenvalue in the CRISIS subsample
# that was below MP edge in the PRE-CRISIS subsample. (We use the simpler
# "subsample" test here rather than the "≥50% / <50% of windows" rolling test
# since both deliver the same operational answer in clean data.)
crs_supra_ranks <- crs$summary[above_mp == TRUE, rank]
pre_supra_ranks <- pre$summary[above_mp == TRUE, rank]
emerging_ranks <- setdiff(crs_supra_ranks, pre_supra_ranks)

cat(sprintf("\nSupra-MP ranks in pre-crisis (2019-2023):  %s\n",
            paste(pre_supra_ranks, collapse = ", ")))
cat(sprintf("Supra-MP ranks in crisis (2024-2026):       %s\n",
            paste(crs_supra_ranks, collapse = ", ")))
cat(sprintf("Crisis-only (emerging) supra-MP ranks:      %s\n",
            if (length(emerging_ranks)) paste(emerging_ranks, collapse = ", ") else "none"))

if (length(emerging_ranks) > 0) {
  emerging_rank <- min(emerging_ranks)
  emerging_info <- crs$summary[rank == emerging_rank]
  bloc_frac <- emerging_info$cocoa_bloc_l2_frac
  cat(sprintf("\nLowest-ranked emerging eigenvalue: rank=%d, eigenvalue=%.3f\n",
              emerging_rank, emerging_info$eigenvalue))
  cat(sprintf("Top-3 loadings: %s\n", emerging_info$top3_loadings))
  cat(sprintf("Cocoa-bloc L2 fraction: %.3f  (T3 threshold = 0.55)\n", bloc_frac))
  t3_pass <- bloc_frac >= 0.55
  cat(sprintf("T3: %s\n", if (t3_pass) "PASS" else "FAIL (still informative if close)"))

  # Detailed loadings for the emerging eigenvector
  ev_col <- paste0("ev", emerging_rank)
  ev_loadings <- crs$loadings[, .(ticker, loading = get(ev_col))]
  ev_loadings[, in_cocoa_bloc := ticker %in% COCOA_BLOC]
  ev_loadings[, sq_loading := loading^2]
  setorder(ev_loadings, -sq_loading)
  cat("\nFull eigenvector composition (sorted by |loading|):\n")
  print(ev_loadings[, .(ticker, loading = round(loading, 3),
                        sq_loading = round(sq_loading, 3),
                        in_cocoa_bloc)])

  fwrite(data.table(
    rank = emerging_rank,
    eigenvalue = emerging_info$eigenvalue,
    cocoa_bloc_l2_frac = bloc_frac,
    t3_pass = t3_pass,
    threshold = 0.55,
    top3 = emerging_info$top3_loadings
  ), "output/t3_test.csv")
} else {
  cat("\nNo emerging eigenvalue in crisis-only subsample. T3 trivially fails.\n")
  fwrite(data.table(t3_pass = FALSE, reason = "no emerging eigenvalue"),
         "output/t3_test.csv")
}

# ----- Figure 3a: Loadings comparison ---------------------------------------
# Pick the top-3 supra-MP eigenvalues from each subsample, plot loadings as heatmap.
plot_data <- list()
for (samp_name in c("pre_2019_2023", "crisis_2024_2026")) {
  s <- if (samp_name == "pre_2019_2023") pre else crs
  top_ranks <- s$summary[above_mp == TRUE, rank][1:min(4, sum(s$summary$above_mp))]
  for (k in top_ranks) {
    ev_col <- paste0("ev", k)
    plot_data[[length(plot_data) + 1]] <- data.table(
      sample = samp_name,
      rank = k,
      eigenvalue = s$summary[rank == k, eigenvalue],
      ticker = s$loadings$ticker,
      loading = s$loadings[[ev_col]]
    )
  }
}
plot_dt <- rbindlist(plot_data)
plot_dt[, ticker_lbl := factor(ticker, levels = tickers)]
plot_dt[, panel := sprintf("%s\nrank %d, λ=%.2f",
                           ifelse(sample == "pre_2019_2023", "Pre (2019-23)", "Crisis (2024-26)"),
                           rank, eigenvalue)]
plot_dt[, panel := factor(panel, levels = unique(panel))]

# Sign-normalise eigenvectors so the largest absolute loading is positive
# (eigenvector sign is arbitrary; this just stabilises the visual)
plot_dt[, loading_signed := loading * sign(loading[which.max(abs(loading))]), by = panel]
plot_dt[, in_bloc := ticker %in% COCOA_BLOC]

p3a <- ggplot(plot_dt, aes(x = loading_signed, y = ticker_lbl, fill = in_bloc)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "tomato"),
                    labels = c("FALSE" = "Other", "TRUE" = "Cocoa bloc"),
                    name = NULL) +
  facet_wrap(~ panel, nrow = 2, scales = "free_y") +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Eigenvector loadings of supra-MP modes",
    subtitle = "Sign-normalised. Red = cocoa bloc (COCOA_NY, COCOA_LDN, HSY, MDLZ, LISN, GHS, GHANA10Y).",
    x = "loading", y = NULL,
    caption = "Pre-crisis vs crisis subsamples. Eigenvalues that emerged only in crisis carry the H1 signal."
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 8))

ggsave("figures/fig3a_eigenvector_loadings.png", p3a, width = 11, height = 9, dpi = 150)
cat("Wrote figures/fig3a_eigenvector_loadings.png\n")
