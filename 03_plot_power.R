# kpower simulation tests -- Step 3: Summary plots from results/summary.csv
#
# Two outputs:
#   1. results/power_heatmap.pdf
#      Grid view: rows = (family, class), columns = length, fill = BIC
#      power. Faceted by K=2 / K=4.
#   2. results/power_K{2,4}_by_{family,class,combined}.pdf  (6 PDFs)
#      Line plots: x=length (log), y=BIC power, 6 lines per panel
#      (3 families x 2 classes), coloured by family / class / both.

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE    <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv. Run 02_run_kpower_tests.R first.")

B <- 20   # number of bootstrap replicates used in 02

# Wilson 95% CI for a binomial proportion -- robust for small n / p near 0 or 1
wilson_ci <- function(p, n, z = 1.959964) {
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}

d <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
ci <- wilson_ci(d$power_BIC, B)
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

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


# 1. Heatmap ----------------------------------------------------------------
heatmap_plot <- ggplot(d, aes(x = factor(seq_length),
                              y = row_label,
                              fill = power_BIC)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * power_BIC)),
            size = 3, colour = "black") +
  facet_wrap(~ paste("K =", true_K), nrow = 1) +
  scale_fill_gradient2(low = "#fde0dd", mid = "#fa9fb5", high = "#7a0177",
                       midpoint = 0.5, limits = c(0, 1),
                       breaks = seq(0, 1, 0.25), name = "BIC power") +
  scale_y_discrete(limits = rev) +
  labs(x = "Sequence length (sites)", y = NULL,
       title = "Power to recover true (family, K) -- by BIC") +
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
  ggplot(sub, aes(x = seq_length, y = power_BIC,
                  colour = .data[[color_var]],
                  fill   = .data[[color_var]],
                  group  = interaction(true_type, tag))) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    x_log + y_pwr +
    labs(x = "Sequence length (sites, log scale)",
         y = "BIC power",
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
