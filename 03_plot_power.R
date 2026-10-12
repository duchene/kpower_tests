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
#     rows = (family, class), columns = length, fill = mean success rate;
#     one panel per K, stacked vertically, with a vertical colourbar at right.
#   results/power_K{2,3}_by_{family,class,combined}
#     line plots: x = length (log), y = mean success rate, shaded band =
#     Wilson 95% CI (drawn on the combined and per-rep panels only).
#
# Every figure is written as .pdf, .svg and .png -- see save_fig() in
# 00_plot_theme.R.

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

# --- Cross-family recovery under each IC (AIC / AICc / BIC) -----------------
# summary.csv only records the honest empirical "correct family & K" call under
# BIC (best_mix_type / best_K). To draw the same honest call under AIC and AICc
# we recompute it from each rep's survey_result.rds: per family take the K that
# minimises the IC on the *empirical* alignment, then pick the family with the
# lowest such IC (all families scored on the SAME empirical data -- apples to
# apples, unlike the biased merged-bootstrap power_AIC/power_AICc columns).
# Result is cached to results/correct_by_ic.csv; set REBUILD_CORRECT_BY_IC=1 to
# force a rescan.
CBI_CSV   <- file.path(OUT_BASE, "correct_by_ic.csv")
IC_LEVELS <- c("AIC", "AICc", "BIC")

correct_from_survey <- function(rds_path, true_type, true_K) {
  s <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(s) || is.null(s$families)) return(NULL)
  do.call(rbind, lapply(IC_LEVELS, function(ic) {
    # best (K, IC value) per family on the empirical alignment
    per_fam <- lapply(s$families, function(f) {
      emp <- f$empirical
      if (is.null(emp) || !ic %in% names(emp) || !nrow(emp)) return(NULL)
      i <- which.min(emp[[ic]])
      data.frame(fam = f$mix_type, bestK = emp$K[i], val = emp[[ic]][i])
    })
    per_fam <- do.call(rbind, per_fam[!vapply(per_fam, is.null, logical(1))])
    if (is.null(per_fam) || !nrow(per_fam)) return(NULL)
    win <- per_fam[which.min(per_fam$val), ]
    data.frame(ic = ic,
               correct_family = win$fam == true_type,
               correct_K      = win$bestK == true_K)
  }))
}

if (file.exists(CBI_CSV) && Sys.getenv("REBUILD_CORRECT_BY_IC") != "1") {
  cbi <- read.csv(CBI_CSV, stringsAsFactors = FALSE)
  message("Loaded cached correct_by_ic.csv (", nrow(cbi), " rows)")
} else {
  message("Scanning survey_result.rds files for AIC/AICc/BIC recovery ...")
  meta <- unique(raw[, c("scenario", "rep", "true_type", "true_K",
                         "tag", "seq_length")])
  rows <- vector("list", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    m   <- meta[i, ]
    rds <- file.path(OUT_BASE, m$scenario,
                     sprintf("rep_%02d", m$rep), "survey_result.rds")
    if (!file.exists(rds)) next
    r <- correct_from_survey(rds, m$true_type, m$true_K)
    if (is.null(r)) next
    rows[[i]] <- data.frame(
      scenario = m$scenario, rep = m$rep, true_type = m$true_type,
      true_K = m$true_K, tag = m$tag, seq_length = m$seq_length, r,
      row.names = NULL)
  }
  cbi <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(cbi) || !nrow(cbi))
    stop("No survey_result.rds files found under ", OUT_BASE)
  cbi$success <- cbi$correct_family & cbi$correct_K
  write.csv(cbi, CBI_CSV, row.names = FALSE)
  message("Wrote ", basename(CBI_CSV), " (", nrow(cbi), " rows)")

  # Sanity check: BIC recovery from RDS should match summary.csv's columns.
  chk <- merge(
    cbi[cbi$ic == "BIC", c("scenario", "rep", "success")],
    data.frame(scenario = raw$scenario, rep = raw$rep,
               success_csv = raw$success),
    by = c("scenario", "rep"))
  mism <- sum(chk$success != chk$success_csv, na.rm = TRUE)
  if (mism > 0)
    message(sprintf("  NOTE: %d/%d BIC calls differ from summary.csv",
                    mism, nrow(chk)))
}
cbi$success <- cbi$correct_family & cbi$correct_K

# Reps run elsewhere may arrive as summary.csv rows only (no survey_result.rds,
# e.g. the L6000 server run). Their BIC call is in summary.csv, so fall back to
# it; AIC/AICc can't be recovered without the RDS and stay missing.
no_rds <- !(paste(raw$scenario, raw$rep) %in%
            paste(cbi$scenario[cbi$ic == "BIC"], cbi$rep[cbi$ic == "BIC"]))
if (any(no_rds)) {
  message(sprintf("  %d rep(s) without survey_result.rds: BIC call taken ",
                  sum(no_rds)), "from summary.csv, AIC/AICc left missing")
  cbi <- rbind(cbi, data.frame(
    raw[no_rds, c("scenario", "rep", "true_type", "true_K", "tag",
                  "seq_length")],
    ic = "BIC", correct_family = raw$correct_family[no_rds],
    correct_K = raw$correct_K[no_rds], success = raw$success[no_rds],
    row.names = NULL))
}

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

# 1. Heatmaps (one per IC) --------------------------------------------------
# Aggregate a per-(scenario, rep) success frame into per-scenario means with
# the same row_label / true_K structure used above.
agg_success <- function(df) {
  a <- aggregate(
    success ~ scenario + true_type + true_K + tag + seq_length,
    data = df, FUN = function(x) c(k = sum(x, na.rm = TRUE),
                                   n = sum(!is.na(x))))
  out <- data.frame(a[, !names(a) %in% "success"],
                    k = as.integer(a$success[, "k"]),
                    n_reps = as.integer(a$success[, "n"]))
  out$mean_success <- ifelse(out$n_reps > 0, out$k / out$n_reps, NA_real_)
  out$true_type <- factor(out$true_type, levels = c("+R", "+H", "+T"))
  out$tag       <- factor(out$tag,       levels = c("sep", "close"))
  out$true_K    <- as.integer(out$true_K)
  out$row_label <- factor(
    paste0(out$true_type, " ", out$tag),
    levels = c("+R sep", "+R close", "+H sep", "+H close", "+T sep", "+T close"))
  out
}

heatmap_for_ic <- function(ic) {
  hd <- agg_success(cbi[cbi$ic == ic, ])
  n_used <- max(hd$n_reps, na.rm = TRUE)
  # Recovery rate is a magnitude on [0, 1] with no meaningful midpoint, so it
  # takes a single-hue sequential ramp (SEQ_COLOURS) -- not a diverging one,
  # which would imply a polarity at 0.5 that the quantity does not have. Tile
  # labels take their ink from the same vector, so they flip to white once the
  # fill is dark enough that near-black text would fail contrast.
  hd$label_ink <- label_ink_for_fill(hd$mean_success)
  p <- ggplot(hd, aes(x = factor(seq_length), y = row_label,
                      fill = mean_success)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", 100 * mean_success),
                  colour = I(label_ink)), size = 2.4) +
    facet_wrap(~ paste("K =", true_K), ncol = 1) +
    scale_fill_gradientn(colours = SEQ_COLOURS,
                         limits = c(0, 1), breaks = seq(0, 1, 0.25),
                         labels = label_percent(accuracy = 1),
                         name = "Correct family & K") +
    scale_y_discrete(limits = rev) +
    guides(fill = guide_colourbar(title.position = "top", barwidth = 0.4,
                                  barheight = 8, ticks.colour = NA)) +
    fig_labs(x = "Sequence length (sites)", y = NULL,
             title = sprintf("Correct family & K recovery (%s, %d reps)",
                             ic, n_used)) +
    theme_kpower(grid = "none") +
    theme(axis.ticks = element_blank(),
          legend.position  = "right",
          legend.direction = "vertical",
          legend.key.width = unit(8, "pt"))
  save_fig(p, file.path(OUT_BASE, sprintf("power_heatmap_%s.pdf", ic)),
           width = "single", height = 6.0)
  if (ic == "BIC")                         # keep back-compat filename
    save_fig(p, file.path(OUT_BASE, "power_heatmap.pdf"),
             width = "single", height = 6.0)
}

for (ic in IC_LEVELS) heatmap_for_ic(ic)


# 2. Line plots -------------------------------------------------------------
# Every series carries two encodings (hue + dash, or hue + panel) so identity
# survives greyscale printing and colour-vision deficiency.
Y_LAB <- "Correct family & K"

# Six series in one panel: confidence bands are omitted here because six
# translucent ribbons overlap into colours that map to no series. The CIs are
# shown in the faceted "combined" and per-replicate figures, where each panel
# carries few enough series to draw them honestly.
power_lines <- function(sub, colour_var, linetype_var, palette,
                        colour_lab, linetype_lab,
                        linetypes = c("solid", "22")) {
  ggplot(sub, aes(x = seq_length, y = mean_success,
                  colour   = .data[[colour_var]],
                  linetype = .data[[linetype_var]],
                  group    = interaction(true_type, tag))) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = palette, name = colour_lab) +
    scale_linetype_manual(values = linetypes, name = linetype_lab) +
    scale_x_length() + scale_y_prop() +
    labs(x = "Sequence length (sites)", y = Y_LAB,
         caption = sprintf("Mean across %d replicates.", n_reps_used)) +
    theme_kpower() +
    theme(legend.box = "horizontal")
}

# Loop over whatever K values are present in the data (e.g. 2, 3)
true_K_values <- sort(unique(d$true_K))
for (K in true_K_values) {
  message(sprintf("Line plots for K = %d", K))
  sub <- d[d$true_K == K, ]

  # by family: hue = generating family, dash = class separation
  p_fam <- power_lines(sub, "true_type", "tag", PALETTE_FAMILY,
                       "Family", "Class") +
    fig_labs(title = sprintf("Recovery vs length (K = %d), by family", K))
  save_fig(p_fam, file.path(OUT_BASE, sprintf("power_K%d_by_family.pdf", K)),
           width = "single", height = 2.9)

  # by class: hue = class separation, dash = family
  p_cls <- power_lines(sub, "tag", "true_type", PALETTE_CLASS,
                       "Class", "Family",
                       linetypes = c("solid", "22", "11")) +
    fig_labs(title = sprintf("Recovery vs length (K = %d), by class", K))
  save_fig(p_cls, file.path(OUT_BASE, sprintf("power_K%d_by_class.pdf", K)),
           width = "single", height = 2.9)

  # combined: the six family x class series as small multiples rather than six
  # cycled hues -- three validated hues, one panel per class.
  # The CI is drawn as a band behind each line rather than as per-point bars.
  # Dodged linerange bars needed the x positions nudged off their true values
  # to stop three families' bars colliding, and read as six separate glyphs per
  # length; a band sits at the true x, needs no dodge, and stays visually
  # subordinate to the line it belongs to.
  p_comb <- ggplot(sub, aes(x = seq_length, y = mean_success,
                            colour = true_type, group = true_type)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = true_type),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.1) +
    facet_wrap(~ tag, nrow = 1,
               labeller = labeller(tag = c(sep = "Separated classes",
                                           close = "Close classes"))) +
    scale_colour_manual(values = PALETTE_FAMILY, name = "Family") +
    scale_fill_manual(values = PALETTE_FAMILY, guide = "none") +
    scale_x_length() + scale_y_prop() +
    labs(x = "Sequence length (sites)", y = Y_LAB,
         caption = sprintf("Mean across %d replicates; bands are Wilson 95%% CIs.",
                           n_reps_used)) +
    fig_labs(title = sprintf("Recovery vs length (K = %d)", K)) +
    theme_kpower()
  save_fig(p_comb, file.path(OUT_BASE, sprintf("power_K%d_combined.pdf", K)),
           width = "onehalf", height = 2.9)
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

per_rep_panel <- function(K) {
  pr <- per_rep[per_rep$true_K == K, ]
  mn <- d[d$true_K == K, ]
  # Per-rep outcomes are binary, so joining them with lines drew a stack of
  # identical segments that hid how many reps sat on each value. Draw them as
  # points jittered along x only (never along y, which would invent
  # intermediate success values), with the mean and its Wilson CI on top.
  ggplot() +
    geom_point(data = pr,
               aes(x = seq_length, y = success_num, colour = true_type),
               position = position_jitter(width = 0.03, height = 0, seed = 1),
               alpha = 0.35, size = 0.6) +
    geom_ribbon(data = mn,
                aes(x = seq_length, ymin = ci_lo, ymax = ci_hi,
                    fill = true_type),
                alpha = 0.15) +
    geom_line(data = mn,
              aes(x = seq_length, y = mean_success, colour = true_type),
              linewidth = 0.7) +
    geom_point(data = mn,
               aes(x = seq_length, y = mean_success, colour = true_type),
               size = 1.3) +
    facet_wrap(~ row_label, ncol = 3,
               labeller = labeller(row_label = function(x) x)) +
    scale_colour_manual(values = PALETTE_FAMILY, guide = "none") +
    scale_fill_manual(values = PALETTE_FAMILY, guide = "none") +
    scale_x_length() +
    scale_y_prop(limits = c(-0.06, 1.06), breaks = c(0, 0.5, 1)) +
    labs(x = "Sequence length (sites)",
         y = "Correct family & K",
         caption = sprintf("Points: individual replicates (n = %d per scenario, jittered horizontally). Line: mean with Wilson 95%% CI.",
                           n_reps_used)) +
    fig_labs(title = sprintf("Per-replicate recovery across length (K = %d)", K)) +
    theme_kpower()
}

for (K in true_K_values) {
  save_fig(per_rep_panel(K),
           file.path(OUT_BASE, sprintf("power_K%d_per_rep.pdf", K)),
           width = "double", height = 3.8)
}

message("\nDone.")
