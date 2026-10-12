# kpower simulation tests -- Step 5: Calibration of kpower's power estimate
#
# QUESTION: when kpower reports "power = p" for a selection, does an
# independent dataset from the same process actually reproduce that selection
# ~p of the time? A power tool is only useful if its self-reported power is
# honest.
#
# Two ingredients, matched like-for-like WITHIN the true family (so this is
# free of the cross-family absolute-BIC bias described in CLAUDE.md):
#
#   PREDICTED (x): kpower's own parametric-bootstrap estimate. For each rep we
#     read the true family's per-bootstrap IC table (families[[true]]$sim_ic in
#     survey_result.rds) and take the fraction of the B bootstrap replicates
#     whose within-family min-BIC K equals the true K. Averaged over reps ->
#     kpower's predicted power to recover true K, per scenario.
#
#   OBSERVED (y): the empirical truth. Because the suite runs N independent
#     replicate alignments per scenario (fresh AliSim draw, same known truth),
#     the fraction of them whose empirical within-true-family K_best equals the
#     true K IS the actual power. Read from summary.csv's R_K / H_K / T_K.
#
# If kpower is well calibrated the points fall on the y = x diagonal. Points
# below the line => kpower OVER-states its power; above => it UNDER-states.
#
# NOTE: predicted power needs the per-family bootstrap IC tables, which live in
# the survey_result.rds files. Scenarios whose RDS have been cleaned up are
# skipped (with a message). Keep the RDS from the current rerun to get the full
# grid; this script consumes whatever RDS are present.
#
# Output: results/calibration.pdf

library(ggplot2)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE    <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV <- file.path(OUT_BASE, "summary.csv")
CRIT        <- "BIC"

if (!file.exists(SUMMARY_CSV))
  stop("No summary.csv. Run 02_run_kpower_tests.R first.")

# Column in summary.csv holding each family's own empirical K_best
fam_K_col <- c("+R" = "R_K", "+H" = "H_K", "+T" = "T_K")

#' Within-family predicted power from a survey RDS: fraction of the true
#' family's bootstrap replicates whose min-IC K (within that family) == K_target.
within_family_power <- function(rds_path, mix_type, K_target, crit = "BIC") {
  s <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(s) || is.null(s$families[[mix_type]]$sim_ic)) return(NA_real_)
  df <- s$families[[mix_type]]$sim_ic
  reps <- unique(df$replicate)
  hit <- vapply(reps, function(r) {
    sub <- df[df$replicate == r & is.finite(df[[crit]]), ]
    if (nrow(sub) == 0) return(NA)
    isTRUE(sub$K[which.min(sub[[crit]])] == K_target)
  }, logical(1))
  mean(hit, na.rm = TRUE)
}

raw <- read.csv(SUMMARY_CSV, stringsAsFactors = FALSE)

# --- OBSERVED: empirical within-true-family K recovery, per rep -------------
raw$emp_K <- vapply(seq_len(nrow(raw)), function(i)
  suppressWarnings(as.integer(raw[i, fam_K_col[raw$true_type[i]]])), integer(1))
raw$observed_hit <- as.integer(raw$emp_K == raw$true_K)

# --- PREDICTED: within-family bootstrap power from each rep's RDS -----------
raw$rds <- file.path(OUT_BASE, raw$scenario,
                     sprintf("rep_%02d", raw$rep), "survey_result.rds")
have <- file.exists(raw$rds)
message(sprintf("RDS present for %d / %d reps.", sum(have), nrow(raw)))
raw$predicted <- NA_real_
raw$predicted[have] <- vapply(which(have), function(i)
  within_family_power(raw$rds[i], raw$true_type[i], raw$true_K[i], CRIT),
  numeric(1))

usable <- raw[!is.na(raw$predicted), ]
if (nrow(usable) == 0)
  stop("No RDS files found, so predicted power cannot be computed. ",
       "Re-run 02 (which saves survey_result.rds per rep) and keep the RDS.")

# --- Aggregate to one point per scenario -----------------------------------
scn <- do.call(rbind, by(usable, usable$scenario, function(g) data.frame(
  scenario   = g$scenario[1],
  true_type  = g$true_type[1],
  true_K     = g$true_K[1],
  tag        = g$tag[1],
  seq_length = g$seq_length[1],
  n          = nrow(g),
  predicted  = mean(g$predicted,    na.rm = TRUE),
  observed   = mean(g$observed_hit, na.rm = TRUE)
)))
scn$true_type <- factor(scn$true_type, levels = c("+R", "+H", "+T"))

family_palette <- c("+R" = "#440154", "+H" = "#21918C", "+T" = "#7AD151")

message(sprintf("Calibration over %d scenarios (%d reps with RDS).",
                nrow(scn), sum(have)))
print(scn[order(scn$predicted),
          c("scenario", "n", "predicted", "observed")], row.names = FALSE)

p <- ggplot(scn, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = true_type, shape = factor(true_K), size = n),
             alpha = 0.85) +
  scale_colour_manual(values = family_palette, name = "True family", drop = FALSE) +
  scale_shape_discrete(name = "True K") +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(x = "Predicted power (kpower bootstrap, within true family)",
       y = "Observed recovery (independent replicates)",
       title = "Calibration of kpower's power estimate",
       subtitle = "On the y = x line = honest; below = over-confident, above = conservative") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

ggsave(file.path(OUT_BASE, "calibration.pdf"), p, width = 6.5, height = 6)
message("Wrote calibration.pdf")
