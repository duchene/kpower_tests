# kpower simulation tests -- Step 3: Summary plots from results/summary.csv
#
# Aggregates across reps. The plotted "power" is the fraction of reps in
# which kpower's empirical fit correctly identifies BOTH the true family
# AND the true K, computed from the `correct_family` and `correct_K`
# columns in summary.csv. These reflect the cross-family BIC comparison
# on the actual empirical alignment for each rep -- the apples-to-apples
# call kpower makes when presented with one alignment.
#
# We do NOT use `power_BIC` from summary.csv: that metric naively merges
# bootstrap IC tables across families, but each family generates its own
# bootstrap alignments, so comparing min(BIC) across the merged table is
# comparing scores on different data and is systematically biased toward
# the more flexible family. See CLAUDE.md "+H bootstrap absolute-BIC bias"
# for the diagnosis.
#
# Per-scenario aggregates:
#   mean_correct  = mean of (correct_family AND correct_K) across reps
#   se_correct    = sqrt(p*(1-p)/n_reps)   (binomial standard error)
#   ci_lo / hi    = Wilson 95% CI on (correct count, n_reps)
#
# Outputs:
#   results/power_heatmap.pdf
#     rows = (family, class), columns = length, fill = mean success rate.
#   results/power_K{2,3}_by_{family,class,combined}.pdf  (6 PDFs)
#     line plots: x = length (log), y = mean success rate, ribbon = Wilson 95% CI.

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE    <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv. Run 02_run_kpower_tests.R first.")

# Wilson 95% CI for a binomial proportion
wilson_ci <- function(k, n, z = 1.959964) {
  p <- ifelse(n > 0, k / n, NA_real_)
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}

raw <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
if (!"rep" %in% names(raw)) raw$rep <- 1L   # back-compat
# Coerce in case CSV stored TRUE/FALSE as strings
raw$correct_family <- as.logical(raw$correct_family)
raw$correct_K      <- as.logical(raw$correct_K)
raw$success        <- raw$correct_family & raw$correct_K

# --- Aggregate across reps --------------------------------------------------
agg <- aggregate(
  success ~ scenario + true_type + true_K + tag + seq_length,
  data = raw, FUN = function(x) c(
    k = sum(x, na.rm = TRUE),
    n = sum(!is.na(x))
  )
)
d <- data.frame(
  agg[, !names(agg) %in% "success"],
  k = as.integer(agg$success[, "k"]),
  n_reps = as.integer(agg$success[, "n"])
)
d$mean_success <- ifelse(d$n_reps > 0, d$k / d$n_reps, NA_real_)
ci <- wilson_ci(d$k, d$n_reps)
d$ci_lo <- ci$lo
d$ci_hi <- ci$hi
d$true_type    <- factor(d$true_type, levels = c("+R", "+H", "+T"))
d$tag          <- factor(d$tag,       levels = c("sep", "close"))
d$true_K       <- as.integer(d$true_K)
d$family_class <- interaction(d$true_type, d$tag, sep = " ")
d$row_label    <- factor(
  paste0(d$true_type, " ", d$tag),
  levels = c("+R sep", "+R close", "+H sep", "+H close", "+T sep", "+T close")
)

n_reps_used <- max(d$n_reps, na.rm = TRUE)
message(sprintf("Aggregating across up to %d reps per scenario", n_reps_used))

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


# 1. Heatmap ----------------------------------------------------------------
heatmap_plot <- ggplot(d, aes(x = factor(seq_length),
                              y = row_label,
                              fill = mean_success)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mean_success)),
            size = 3, colour = "black") +
  facet_wrap(~ paste("K =", true_K), nrow = 1) +
  scale_fill_gradient2(low = "#fde0dd", mid = "#fa9fb5", high = "#7a0177",
                       midpoint = 0.5, limits = c(0, 1),
                       breaks = seq(0, 1, 0.25),
                       name = "P(correct\nfamily & K)") +
  scale_y_discrete(limits = rev) +
  labs(x = "Sequence length (sites)", y = NULL,
       title = sprintf("Correct family & K recovery (mean across %d reps, BIC)",
                       n_reps_used)) +
  base_theme +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_BASE, "power_heatmap.pdf"),
       heatmap_plot, width = 9, height = 4)
message("Heatmap -> power_heatmap.pdf")


# 2. Line plots -------------------------------------------------------------
x_log <- scale_x_log10(breaks = c(100, 300, 1000, 3000, 6000),
                       labels = c("100", "300", "1k", "3k", "6k"))
y_pwr <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

plot_one <- function(df, K, color_var, color_lab, title) {
  sub <- df[df$true_K == K, ]
  ggplot(sub, aes(x = seq_length, y = mean_success,
                  colour = .data[[color_var]],
                  fill   = .data[[color_var]],
                  group  = interaction(true_type, tag))) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    x_log + y_pwr +
    labs(x = "Sequence length (sites, log scale)",
         y = "P(correct family & K)",
         colour = color_lab, fill = color_lab, title = title) +
    base_theme + theme(legend.position = "right")
}

save_pdf <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_BASE, file), p, width = w, height = h)
  message("  ", file)
}

# Loop over whatever K values are present in the data (e.g. 2, 3)
true_K_values <- sort(unique(d$true_K))
for (K in true_K_values) {
  message(sprintf("Line plots for K = %d", K))
  save_pdf(
    plot_one(d, K, "true_type", "Family",
             sprintf("Recovery vs length (K = %d) -- by family", K)) +
      aes(linetype = tag) + labs(linetype = "Class"),
    sprintf("power_K%d_by_family.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "tag", "Class",
             sprintf("Recovery vs length (K = %d) -- by class", K)) +
      aes(linetype = true_type) + labs(linetype = "Family"),
    sprintf("power_K%d_by_class.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "family_class", "Family x Class",
             sprintf("Recovery vs length (K = %d) -- combined", K)),
    sprintf("power_K%d_combined.pdf", K)
  )
}

# 3. Per-rep line plots -----------------------------------------------------
# One line per rep (binary 0/1 success across lengths), faceted by family x
# class so each panel shows ~n_reps lines. Mean line overlaid in bold.

# Per-rep success: binary 1/0 per (scenario, rep)
per_rep <- raw[, c("scenario", "rep", "true_type", "true_K", "tag",
                   "seq_length", "success")]
per_rep$success_num <- as.integer(per_rep$success)
per_rep$true_type   <- factor(per_rep$true_type, levels = c("+R", "+H", "+T"))
per_rep$tag         <- factor(per_rep$tag,       levels = c("sep", "close"))
per_rep$row_label   <- factor(
  paste0(per_rep$true_type, " ", per_rep$tag),
  levels = c("+R sep", "+R close", "+H sep", "+H close", "+T sep", "+T close")
)
# matching mean rows (built earlier as `d`)
d$row_label <- factor(d$row_label, levels = levels(per_rep$row_label))

family_palette <- c("+R" = "#d6604d", "+H" = "#4393c3", "+T" = "#2ca02c")

per_rep_panel <- function(K) {
  pr <- per_rep[per_rep$true_K == K, ]
  mn <- d[d$true_K == K, ]
  ggplot() +
    geom_line(data = pr,
              aes(x = seq_length, y = success_num,
                  group = rep, colour = true_type),
              alpha = 0.30, linewidth = 0.4) +
    geom_line(data = mn,
              aes(x = seq_length, y = mean_success, colour = true_type),
              linewidth = 1.2) +
    geom_point(data = mn,
               aes(x = seq_length, y = mean_success, colour = true_type),
               size = 2) +
    facet_wrap(~ row_label, ncol = 2,
               labeller = labeller(row_label = function(x) x)) +
    scale_colour_manual(values = family_palette, guide = "none") +
    x_log +
    scale_y_continuous(limits = c(-0.05, 1.05),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "Sequence length (sites, log scale)",
         y = "Success (per rep: 0 or 1; bold = mean)",
         title = sprintf("Per-rep recovery across length (K = %d, n=%d reps)",
                         K, n_reps_used)) +
    base_theme +
    theme(strip.text = element_text(face = "bold"))
}

for (K in true_K_values) {
  save_pdf(per_rep_panel(K),
           sprintf("power_K%d_per_rep.pdf", K),
           w = 8, h = 7)
}

message("\nDone.")
