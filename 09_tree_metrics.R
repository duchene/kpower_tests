# kpower simulation tests -- Step 5: Tree-shape metrics (stemminess, RRT-CoV)
#
# For every fitted tree produced during the survey -- both the empirical fit and
# each parametric-bootstrap replicate, at every K -- we compute two tree-shape
# statistics and plot them "through K" in the same empirical-vs-reps overlay
# style as the IC-profile / BIC-vs-K plots:
#
#   stemminess (= treeness) = sum(internal branch lengths) / sum(all branch lengths)
#   rrt_cov                 = CoV (sd/mean) of root-to-tip distances across tips,
#                             after midpoint-rooting the (unrooted) fitted tree.
#
# Trees are read from each rep's survey_runs:
#   results/<scenario>/rep_<NN>/survey_runs/survey_<FAM>_linked/
#     empirical/empirical_K<k>/empirical_K<k>.treefile      (source = "empirical")
#     sim_fits/sim<NNNN>_K<k>/sim<NNNN>_K<k>.treefile        (source = "sim")
#
# The (slow) extraction is cached to results/tree_metrics.csv; set
# REBUILD_TREE_METRICS=1 to force a rescan. Plots are always redrawn.
#
# Outputs (multi-page PDFs, one page per scenario):
#   results/stemminess.pdf
#   results/rrt_cov.pdf

suppressMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
})

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
# Results directory: defaults to results/. Set KPOWER_RESULTS_DIR (absolute, or
# relative to this script) to plot an archived run, e.g. results_B50.
OUT_BASE   <- local({
  o <- Sys.getenv("KPOWER_RESULTS_DIR", "results")
  if (startsWith(o, "/")) o else file.path(SCRIPT_DIR, o)
})
TM_CSV     <- file.path(OUT_BASE, "tree_metrics.csv")

if (!dir.exists(OUT_BASE))
  stop("No results/ directory. Run 02_run_kpower_tests.R first.")

# Shared publication theme, palettes and save helpers.
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))

FAM_MAP <- c(R = "+R", H = "+H", T = "+T")

# --- Metrics for a single tree (averaged over class trees if multiPhylo) ----
tree_metrics_one <- function(tr) {
  el <- tr$edge.length
  nt <- Ntip(tr)
  if (is.null(el) || !length(el) || sum(el) <= 0)
    return(c(stemminess = NA_real_, rrt_cov = NA_real_))
  internal   <- tr$edge[, 2] > nt              # child node is internal
  stemminess <- sum(el[internal]) / sum(el)
  rt <- tryCatch(phangorn::midpoint(tr), error = function(e) NULL)
  rrt_cov <- NA_real_
  if (!is.null(rt)) {
    d <- node.depth.edgelength(rt)[seq_len(Ntip(rt))]   # root-to-tip distances
    m <- mean(d)
    if (is.finite(m) && m > 0) rrt_cov <- stats::sd(d) / m
  }
  c(stemminess = stemminess, rrt_cov = rrt_cov)
}

tree_metrics <- function(path) {
  tr <- tryCatch(read.tree(path), error = function(e) NULL)
  if (is.null(tr)) return(c(stemminess = NA_real_, rrt_cov = NA_real_))
  if (inherits(tr, "multiPhylo")) {           # +T class trees -> average
    m <- vapply(tr, tree_metrics_one, numeric(2))
    return(c(stemminess = mean(m["stemminess", ], na.rm = TRUE),
             rrt_cov    = mean(m["rrt_cov", ],    na.rm = TRUE)))
  }
  tree_metrics_one(tr)
}

# --- Extraction (cached) ----------------------------------------------------
if (file.exists(TM_CSV) && Sys.getenv("REBUILD_TREE_METRICS") != "1") {
  tm <- read.csv(TM_CSV, stringsAsFactors = FALSE)
  message("Loaded cached tree_metrics.csv (", nrow(tm), " rows)")
} else {
  scen_dirs <- list.dirs(OUT_BASE, recursive = FALSE)
  scen_dirs <- scen_dirs[grepl("_K[0-9]+_", basename(scen_dirs))]
  message("Extracting tree metrics from ", length(scen_dirs), " scenarios ...")
  per_scen <- vector("list", length(scen_dirs))
  for (si in seq_along(scen_dirs)) {
    sc  <- basename(scen_dirs[si])
    fs  <- list.files(scen_dirs[si], pattern = "\\.treefile$",
                      recursive = TRUE, full.names = TRUE)
    # Keep only per-K model fits; exclude +T window_trees (candidate trees).
    fs  <- fs[grepl("survey_[RHT]_linked", fs) &
              grepl("/(empirical|sim_fits)/", fs)]
    if (!length(fs)) next
    node   <- sub("\\.treefile$", "", basename(fs))
    famL   <- sub(".*survey_([RHT])_linked.*", "\\1", fs)
    rep    <- as.integer(sub(".*/rep_([0-9]+)/.*", "\\1", fs))
    K      <- as.integer(sub(".*_K([0-9]+)$", "\\1", node))
    source <- ifelse(grepl("^empirical", node), "empirical", "sim")
    sim_id <- ifelse(source == "sim", sub("_K[0-9]+$", "", node), NA_character_)
    vals   <- vapply(fs, tree_metrics, numeric(2))
    per_scen[[si]] <- data.frame(
      scenario = sc, rep = rep, family = FAM_MAP[famL], K = K,
      source = source, sim_id = sim_id,
      stemminess = vals["stemminess", ], rrt_cov = vals["rrt_cov", ],
      row.names = NULL, stringsAsFactors = FALSE)
    message(sprintf("  [%2d/%2d] %-18s %d trees",
                    si, length(scen_dirs), sc, length(fs)))
  }
  tm <- do.call(rbind, per_scen[!vapply(per_scen, is.null, logical(1))])
  if (is.null(tm) || !nrow(tm))
    stop("No treefiles found under ", OUT_BASE)
  write.csv(tm, TM_CSV, row.names = FALSE)
  message("Wrote ", basename(TM_CSV), " (", nrow(tm), " rows)")
}

# Drop any non-fit trees that slipped in (e.g. +T window_trees -> NA K).
n_bad <- sum(is.na(tm$K))
if (n_bad) message("Dropping ", n_bad, " non-fit trees (NA K)")
tm <- tm[!is.na(tm$K), ]
tm$family <- factor(tm$family, levels = c("+R", "+H", "+T"))

# --- Scenario metadata (for titles) from the scenario name ------------------
# e.g. "H_K3_sep_L1000" -> true family +H, K 3, class sep, length 1000
parse_scen <- function(sc) {
  p <- strsplit(sc, "_")[[1]]
  list(true_type = FAM_MAP[p[1]],
       true_K    = as.integer(sub("^K", "", p[2])),
       tag       = p[3],
       seq_length = as.integer(sub("^L", "", p[4])))
}

# --- Plot one scenario page for a given metric ------------------------------
# These are diagnostic atlases -- one page per scenario -- so each page keeps
# its title: it is how you find your way through 48 pages. Manuscript figures
# (03, 04) suppress titles in favour of captions.
plot_scenario <- function(sc, metric, metric_lab) {
  d  <- tm[tm$scenario == sc & !is.na(tm[[metric]]), ]
  if (!nrow(d)) return(NULL)
  meta <- parse_scen(sc)
  emp  <- d[d$source == "empirical", ]
  sim  <- d[d$source == "sim", ]
  ggplot(mapping = aes(x = K, y = .data[[metric]], colour = family)) +
    geom_line(data = sim,
              aes(group = interaction(rep, sim_id, family)),
              alpha = 0.10, linewidth = 0.25) +
    geom_line(data = emp,
              aes(group = interaction(rep, family)),
              linewidth = 0.7) +
    scale_colour_manual(values = PALETTE_FAMILY, name = "Fitted family",
                        drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(d$K))) +
    labs(x = "K (number of classes)", y = metric_lab,
         title = sprintf("%s  (true %s, K = %d, %s, L = %d)",
                         sc, meta$true_type, meta$true_K,
                         meta$tag, meta$seq_length),
         subtitle = "Heavy lines: empirical fits. Faint lines: bootstrap replicates.") +
    theme_kpower(base_size = 9, grid = "both")
}

scenarios <- sort(unique(tm$scenario))

save_metric_pdf <- function(metric, metric_lab, file) {
  path <- file.path(OUT_BASE, file)
  open_atlas(path, width = 6.5, height = 4.5)   # cairo_pdf: embeds fonts
  on.exit(grDevices::dev.off())
  n <- 0L
  for (sc in scenarios) {
    p <- plot_scenario(sc, metric, metric_lab)
    if (!is.null(p)) { print(p); n <- n + 1L }
  }
  message(sprintf("%s -> %s (%d scenario pages)", metric, file, n))
}

save_metric_pdf("stemminess", "Stemminess (treeness)", "stemminess.pdf")
save_metric_pdf("rrt_cov",    "Root-to-tip CoV",        "rrt_cov.pdf")

message("\nDone.")
