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
OUT_BASE    <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

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

family_palette <- c("+R" = "#d6604d", "+H" = "#4393c3", "+T" = "#2ca02c")

x_log <- scale_x_log10(breaks = c(100, 300, 1000, 3000, 6000),
                       labels = c("100", "300", "1k", "3k", "6k"))
y_pwr <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom")

# Plot one panel per K present in data (e.g. 2 and 3)
true_K_values <- sort(unique(d$true_K))
for (K in true_K_values) {
  sub <- d[d$true_K == K, ]
  p <- ggplot(sub, aes(x = seq_length, y = frac,
                       colour = picked, fill = picked,
                       group  = picked)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_grid(true_type ~ tag,
               labeller = labeller(true_type = function(x) paste("true:", x),
                                   tag       = function(x) paste("class:", x))) +
    x_log + y_pwr +
    scale_colour_manual(values = family_palette, name = "kpower picked") +
    scale_fill_manual(values = family_palette, guide = "none") +
    labs(x = "Sequence length (sites, log scale)",
         y = "Fraction of reps kpower picks this family",
         title = sprintf("Empirical family selection (K = %d, mean across %d reps)",
                         K, n_reps_used)) +
    base_theme

  out <- file.path(OUT_BASE, sprintf("selection_K%d.pdf", K))
  ggsave(out, p, width = 8, height = 7)
  message("Wrote ", basename(out))
}

message("\nDone.")
