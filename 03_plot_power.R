# kpower simulation tests -- Step 3: Summary plots from results/summary.csv
#
# Aggregates across reps (Stage 2: multiple realisations per scenario).
# For each scenario:
#   mean_power  = mean of power_BIC across reps
#   se_power    = standard error of the mean across reps (sd / sqrt(n_reps))
#   ci_lo / hi  = mean +/- 1.96 * se (between-rep 95% CI; Wald)
#
# This is the *outer* CI -- it captures alignment-realisation variation,
# which is the variance component the within-rep bootstrap CI cannot see.
# The within-rep Wilson CI on B is still informative but is no longer
# the headline -- it's now nested inside the across-rep aggregate.
#
# Outputs:
#   results/power_heatmap.pdf
#     rows = (family, class), columns = length, fill = mean BIC power.
#   results/power_K{2,4}_by_{family,class,combined}.pdf  (6 PDFs)
#     line plots: x = length (log), y = mean BIC power, ribbon = +/- 1.96 SE.

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE    <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv. Run 02_run_kpower_tests.R first.")

raw <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
if (!"rep" %in% names(raw)) raw$rep <- 1L   # back-compat with pre-Stage-2 csv

# --- Aggregate across reps --------------------------------------------------
agg <- aggregate(
  power_BIC ~ scenario + true_type + true_K + tag + seq_length,
  data = raw, FUN = function(x) c(
    mean   = mean(x, na.rm = TRUE),
    sd     = if (sum(!is.na(x)) > 1) sd(x, na.rm = TRUE) else 0,
    n      = sum(!is.na(x))
  )
)
# aggregate returns a matrix in $power_BIC; flatten
d <- data.frame(
  agg[, !names(agg) %in% "power_BIC"],
  mean_power = agg$power_BIC[, "mean"],
  sd_power   = agg$power_BIC[, "sd"],
  n_reps     = as.integer(agg$power_BIC[, "n"])
)
d$se_power <- ifelse(d$n_reps > 1, d$sd_power / sqrt(d$n_reps), NA_real_)
d$ci_lo    <- pmax(0, d$mean_power - 1.959964 * d$se_power)
d$ci_hi    <- pmin(1, d$mean_power + 1.959964 * d$se_power)
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
                              fill = mean_power)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mean_power)),
            size = 3, colour = "black") +
  facet_wrap(~ paste("K =", true_K), nrow = 1) +
  scale_fill_gradient2(low = "#fde0dd", mid = "#fa9fb5", high = "#7a0177",
                       midpoint = 0.5, limits = c(0, 1),
                       breaks = seq(0, 1, 0.25), name = "BIC power\n(mean)") +
  scale_y_discrete(limits = rev) +
  labs(x = "Sequence length (sites)", y = NULL,
       title = sprintf("Power to recover true (family, K) -- mean BIC across %d reps",
                       n_reps_used)) +
  base_theme +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_BASE, "power_heatmap.pdf"),
       heatmap_plot, width = 9, height = 4)
message("Heatmap -> power_heatmap.pdf")


# 2. Line plots -------------------------------------------------------------
x_log <- scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000),
                       labels = c("100", "300", "1k", "3k", "10k"))
y_pwr <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

plot_one <- function(df, K, color_var, color_lab, title) {
  sub <- df[df$true_K == K, ]
  ggplot(sub, aes(x = seq_length, y = mean_power,
                  colour = .data[[color_var]],
                  fill   = .data[[color_var]],
                  group  = interaction(true_type, tag))) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    x_log + y_pwr +
    labs(x = "Sequence length (sites, log scale)",
         y = "BIC power (mean across reps)",
         colour = color_lab, fill = color_lab, title = title) +
    base_theme + theme(legend.position = "right")
}

save_pdf <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_BASE, file), p, width = w, height = h)
  message("  ", file)
}

for (K in c(2, 4)) {
  message(sprintf("Line plots for K = %d", K))
  save_pdf(
    plot_one(d, K, "true_type", "Family",
             sprintf("BIC power vs length (K = %d) -- by family", K)) +
      aes(linetype = tag) + labs(linetype = "Class"),
    sprintf("power_K%d_by_family.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "tag", "Class",
             sprintf("BIC power vs length (K = %d) -- by class", K)) +
      aes(linetype = true_type) + labs(linetype = "Family"),
    sprintf("power_K%d_by_class.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "family_class", "Family x Class",
             sprintf("BIC power vs length (K = %d) -- combined", K)),
    sprintf("power_K%d_combined.pdf", K)
  )
}

message("\nDone.")
