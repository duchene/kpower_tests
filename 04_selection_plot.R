# kpower simulation tests -- Step 4: Family-selection plot
#
# For each (true_type, true_K, tag, seq_length, kpower_picked), counts the
# fraction of reps in which kpower's empirical fit selected each family.
# This is built from `best_mix_type` in summary.csv -- the cross-family
# BIC winner on the actual empirical alignment per rep.
#
# We do NOT reconstruct selection fractions from survey_result.rds bootstrap
# tables anymore: those mix BIC scores across independently-generated
# bootstrap alignments and are systematically biased. See CLAUDE.md
# "+H bootstrap absolute-BIC bias" for the diagnosis.
#
# Outputs:
#   results/selection_K2.pdf, results/selection_K3.pdf (one per K present)
# Each PDF is a 6-panel grid:
#   rows    = true family (+R, +H, +T)
#   columns = class (sep, close)
#   per panel: x = sequence length (log), y = fraction of reps where
#              kpower picked each candidate family, with Wilson 95% CI.
# Ideal: the line matching the panel's true family rises to 1 with length
# while the others fall to 0.

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
# Results directory: defaults to results/. Set KPOWER_RESULTS_DIR (absolute, or
# relative to this script) to plot an archived run, e.g. results_B50.
OUT_BASE    <- local({
  o <- Sys.getenv("KPOWER_RESULTS_DIR", "results")
  if (startsWith(o, "/")) o else file.path(SCRIPT_DIR, o)
})
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

# Shared publication theme, palettes and save helpers.
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv. Run 02_run_kpower_tests.R first.")

wilson_ci <- function(k, n, z = 1.959964) {
  p <- ifelse(n > 0, k / n, NA_real_)
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}

raw <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
if (!"rep" %in% names(raw)) raw$rep <- 1L

# Build counts: for each (scenario, candidate family), how many reps did
# kpower pick that family?
all_families <- c("+R", "+H", "+T")
rows <- lapply(all_families, function(fa) {
  r <- raw
  r$picked <- fa
  r$picked_this <- as.integer(raw$best_mix_type == fa)
  r[, c("scenario", "true_type", "true_K", "tag", "seq_length",
        "rep", "picked", "picked_this")]
})
long <- do.call(rbind, rows)

agg <- aggregate(
  picked_this ~ scenario + true_type + true_K + tag + seq_length + picked,
  data = long, FUN = function(x) c(
    k = sum(x, na.rm = TRUE),
    n = sum(!is.na(x))
  )
)
d <- data.frame(
  agg[, !names(agg) %in% "picked_this"],
  k = as.integer(agg$picked_this[, "k"]),
  n_reps = as.integer(agg$picked_this[, "n"])
)
d$frac  <- ifelse(d$n_reps > 0, d$k / d$n_reps, NA_real_)
ci <- wilson_ci(d$k, d$n_reps)
d$ci_lo <- ci$lo
d$ci_hi <- ci$hi
d$true_type <- factor(d$true_type, levels = c("+R", "+H", "+T"))
d$tag       <- factor(d$tag,       levels = c("sep", "close"))
d$picked    <- factor(d$picked,    levels = c("+R", "+H", "+T"))

n_reps_used <- max(d$n_reps, na.rm = TRUE)
message(sprintf("Selection plots: aggregating across up to %d reps", n_reps_used))

# Plot one panel per K present in data (e.g. 2 and 3)
true_K_values <- sort(unique(d$true_K))
for (K in true_K_values) {
  sub <- d[d$true_K == K, ]
  # A faint rule at 1.0 gives the eye the target each panel's true-family line
  # should approach, without adding a second scale.
  # CIs as dodged bars rather than ribbons: three overlapping translucent
  # bands mixed into hues that matched none of the three families.
  dodge <- position_dodge(width = 0.09)   # log10 units: ~5% of x
  p <- ggplot(sub, aes(x = seq_length, y = frac,
                       colour = picked, group = picked)) +
    geom_hline(yintercept = 1, linewidth = 0.3, linetype = "22",
               colour = INK[["grid"]]) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), position = dodge,
                   linewidth = 0.4, alpha = 0.55) +
    geom_line(linewidth = 0.5, position = dodge) +
    geom_point(size = 1.1, position = dodge) +
    facet_grid(true_type ~ tag,
               labeller = labeller(
                 true_type = function(x) paste("true", x),
                 tag       = c(sep = "Separated classes",
                               close = "Close classes"))) +
    scale_x_length() + scale_y_prop() +
    scale_colour_manual(values = PALETTE_FAMILY, name = "Family selected") +
    labs(x = "Sequence length (sites)",
         y = "Replicates selecting this family",
         caption = sprintf("Mean across %d replicates; bars are Wilson 95%% CIs, dodged for legibility.",
                           n_reps_used)) +
    fig_labs(title = sprintf("Empirical family selection (K = %d)", K)) +
    theme_kpower()

  save_fig(p, file.path(OUT_BASE, sprintf("selection_K%d.pdf", K)),
           width = "onehalf", height = 4.4)
}

message("\nDone.")
