# kpower simulation tests -- Step 4: Family selection plots
#
# For each scenario, computes the fraction of B=20 bootstrap replicates
# that selected each mixture family (regardless of K, BIC-best across
# all surveyed families and Ks). Produces two PDFs:
#
#   results/selection_K2.pdf
#   results/selection_K4.pdf
#
# Each PDF is a 6-panel grid:
#   rows    = true family (+R, +H, +T)
#   columns = class (sep, close)
#   per panel: x = sequence length (log), y = selection rate (0..1)
#   three lines, one per family kpower could pick (+R/+H/+T).
#
# The panel's own row colour is the "correct" family -- ideally rises to 1
# with length while the others fall to 0.

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE <- file.path(SCRIPT_DIR, "results")

rds_files <- list.files(OUT_BASE, pattern = "^survey_result\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if (length(rds_files) == 0)
  stop("No survey_result.rds files yet. Run 02 first.")


#' For one survey result, return data frame with one row per family with
#' columns: mix_type, frac (fraction of bootstrap reps where this family
#' won the BIC across all surveyed families and K values).
selection_fractions <- function(surv) {
  fams <- surv$families
  ok   <- names(fams)[!vapply(fams, is.null, logical(1))]
  if (length(ok) == 0)
    return(data.frame(mix_type = character(), frac = numeric()))

  long <- do.call(rbind, lapply(ok, function(mt) {
    df <- fams[[mt]]$sim_ic
    if (is.null(df)) return(NULL)
    df$mix_type <- mt
    df[, c("replicate", "mix_type", "K", "BIC")]
  }))
  if (is.null(long) || nrow(long) == 0)
    return(data.frame(mix_type = character(), frac = numeric()))

  reps <- unique(long$replicate)
  winners <- vapply(reps, function(r) {
    sub <- long[long$replicate == r & is.finite(long$BIC), ]
    if (nrow(sub) == 0) return(NA_character_)
    sub$mix_type[which.min(sub$BIC)]
  }, character(1))
  winners <- winners[!is.na(winners)]

  all_families <- c("+R", "+H", "+T")
  tab <- table(factor(winners, levels = all_families))
  data.frame(mix_type = names(tab),
             frac     = as.numeric(tab) / length(winners),
             stringsAsFactors = FALSE)
}


# Build long data frame across all scenarios
rows <- list()
for (rds in rds_files) {
  scenario <- basename(dirname(rds))
  surv     <- readRDS(rds)
  params   <- readLines(file.path(SCRIPT_DIR, "alignments", scenario,
                                  "sim_params.txt"))
  plist <- setNames(sub("^[^:]+:\\s*", "", params), sub(":.*", "", params))

  sf <- selection_fractions(surv)
  if (nrow(sf) == 0) next
  sf$scenario   <- scenario
  sf$true_type  <- plist[["type"]]
  sf$true_K     <- as.integer(plist[["K"]])
  sf$tag        <- plist[["tag"]]
  sf$seq_length <- as.integer(plist[["length"]])
  rows[[length(rows) + 1]] <- sf
}

d <- do.call(rbind, rows)
d$true_type <- factor(d$true_type, levels = c("+R", "+H", "+T"))
d$tag       <- factor(d$tag,       levels = c("sep", "close"))
d$mix_type  <- factor(d$mix_type,  levels = c("+R", "+H", "+T"))

family_palette <- c("+R" = "#d6604d", "+H" = "#4393c3", "+T" = "#2ca02c")

x_log <- scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000),
                       labels = c("100", "300", "1k", "3k", "10k"))
y_pwr <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom")

for (K in c(2, 4)) {
  sub <- d[d$true_K == K, ]
  p <- ggplot(sub, aes(x = seq_length, y = frac,
                       colour = mix_type, group = mix_type)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_grid(true_type ~ tag,
               labeller = labeller(true_type = function(x) paste("true:", x),
                                   tag       = function(x) paste("class:", x))) +
    x_log + y_pwr +
    scale_colour_manual(values = family_palette, name = "kpower picked") +
    labs(x = "Sequence length (sites, log scale)",
         y = "Fraction of bootstraps selecting this family",
         title = sprintf("Family selection across bootstraps (K = %d)", K)) +
    base_theme

  out <- file.path(OUT_BASE, sprintf("selection_K%d.pdf", K))
  ggsave(out, p, width = 8, height = 7)
  message("Wrote ", basename(out))
}

message("\nDone.")
