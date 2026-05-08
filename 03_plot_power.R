# kpower simulation tests -- Step 3: Plot BIC ground-truth power
#
# Produces 6 PDFs from results/summary.csv:
#   power_K2_by_family.pdf      power_K4_by_family.pdf
#   power_K2_by_class.pdf       power_K4_by_class.pdf
#   power_K2_combined.pdf       power_K4_combined.pdf
#
# y = BIC ground-truth power (fraction of B bootstraps where the best
#     (family, K) across all families == true (family, K))
# x = sequence length (log10 scale)

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv found. Run 02_run_kpower_tests.R first.")

d <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)
d$true_type <- factor(d$true_type, levels = c("+R", "+H", "+T"))
d$tag       <- factor(d$tag,       levels = c("sep", "close"))
d$true_K    <- as.integer(d$true_K)
d$family_class <- interaction(d$true_type, d$tag, sep = " ")

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

x_log <- scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000),
                       labels = c("100", "300", "1k", "3k", "10k"))
y_pwr <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

plot_one <- function(df, K, color_var, color_lab, title) {
  sub <- df[df$true_K == K, ]
  ggplot(sub, aes(x = seq_length, y = power_BIC,
                  colour = .data[[color_var]],
                  group  = interaction(true_type, tag))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    x_log + y_pwr +
    labs(x = "Sequence length (sites, log scale)",
         y = "BIC power (fraction of bootstraps recovering true family + K)",
         colour = color_lab,
         title  = title) +
    base_theme
}

save_pdf <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_BASE, file), p, width = w, height = h)
  message("  ", file)
}

for (K in c(2, 4)) {
  message(sprintf("Plots for K = %d", K))
  save_pdf(
    plot_one(d, K, "true_type", "Family",
             sprintf("BIC power vs length (K = %d) -- by family", K)) +
      aes(linetype = tag) +
      labs(linetype = "Class"),
    sprintf("power_K%d_by_family.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "tag", "Class",
             sprintf("BIC power vs length (K = %d) -- by class", K)) +
      aes(linetype = true_type) +
      labs(linetype = "Family"),
    sprintf("power_K%d_by_class.pdf", K)
  )
  save_pdf(
    plot_one(d, K, "family_class", "Family x Class",
             sprintf("BIC power vs length (K = %d) -- combined", K)),
    sprintf("power_K%d_combined.pdf", K)
  )
}

message("\nDone. Plots in: ", OUT_BASE)
