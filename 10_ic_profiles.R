# kpower simulation tests -- Step 6: IC-through-K overlays (AIC / AICc / BIC)
#
# The per-rep ic_profile_{R,H,T}.pdf files produced during the survey show BIC
# vs K (empirical line + thin bootstrap-replicate lines) for a single rep and
# family. This script draws the same "IC through K" overlay for ALL three ICs,
# summarised per scenario: every family on one panel, all 10 outer reps'
# empirical profiles overlaid in bold, and their bootstrap replicates pooled as
# thin lines.
#
# Data come from each rep's survey_result.rds:
#   $families[[fam]]$empirical  (K, AIC, AICc, BIC)
#   $families[[fam]]$sim_ic     (K, AIC, AICc, BIC, replicate)
#
# The pooled IC table is cached to results/ic_tables.csv; set REBUILD_IC_TABLES=1
# to force a rescan. Plots are always redrawn.
#
# Outputs (multi-page PDFs, one page per scenario):
#   results/ic_profiles_AIC.pdf, ic_profiles_AICc.pdf, ic_profiles_BIC.pdf

suppressMessages(library(ggplot2))

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
# Results directory: defaults to results/. Set KPOWER_RESULTS_DIR (absolute, or
# relative to this script) to plot an archived run, e.g. results_B50.
OUT_BASE  <- local({
  o <- Sys.getenv("KPOWER_RESULTS_DIR", "results")
  if (startsWith(o, "/")) o else file.path(SCRIPT_DIR, o)
})
IC_CSV    <- file.path(OUT_BASE, "ic_tables.csv")
IC_LEVELS <- c("AIC", "AICc", "BIC")

if (!dir.exists(OUT_BASE))
  stop("No results/ directory. Run 02_run_kpower_tests.R first.")

# Shared publication theme, palettes and save helpers.
source(file.path(SCRIPT_DIR, "00_plot_theme.R"))

FAM_MAP <- c(R = "+R", H = "+H", T = "+T")

keep <- c("K", "AIC", "AICc", "BIC")

# --- Pooled IC tables (cached) ----------------------------------------------
if (file.exists(IC_CSV) && Sys.getenv("REBUILD_IC_TABLES") != "1") {
  icd <- read.csv(IC_CSV, stringsAsFactors = FALSE)
  message("Loaded cached ic_tables.csv (", nrow(icd), " rows)")
} else {
  scen_dirs <- list.dirs(OUT_BASE, recursive = FALSE)
  scen_dirs <- scen_dirs[grepl("_K[0-9]+_", basename(scen_dirs))]
  message("Reading survey_result.rds from ", length(scen_dirs), " scenarios ...")
  per_scen <- vector("list", length(scen_dirs))
  for (si in seq_along(scen_dirs)) {
    sc   <- basename(scen_dirs[si])
    reps <- list.dirs(scen_dirs[si], recursive = FALSE)
    reps <- reps[grepl("rep_[0-9]+$", reps)]
    rows <- list()
    for (rp in reps) {
      rep_n <- as.integer(sub(".*rep_([0-9]+)$", "\\1", rp))
      rds   <- file.path(rp, "survey_result.rds")
      if (!file.exists(rds)) next
      s <- tryCatch(readRDS(rds), error = function(e) NULL)
      if (is.null(s) || is.null(s$families)) next
      for (f in s$families) {
        emp <- f$empirical
        sim <- f$sim_ic
        if (!is.null(emp) && nrow(emp)) {
          e <- emp[, keep]; e$source <- "empirical"; e$line <- "emp"
          e$scenario <- sc; e$rep <- rep_n; e$family <- f$mix_type
          rows[[length(rows) + 1]] <- e
        }
        if (!is.null(sim) && nrow(sim)) {
          m <- sim[, keep]; m$source <- "sim"
          m$line <- paste0("sim", sim$replicate)
          m$scenario <- sc; m$rep <- rep_n; m$family <- f$mix_type
          rows[[length(rows) + 1]] <- m
        }
      }
    }
    if (length(rows))
      per_scen[[si]] <- do.call(rbind, rows)
    message(sprintf("  [%2d/%2d] %-18s", si, length(scen_dirs), sc))
  }
  icd <- do.call(rbind, per_scen[!vapply(per_scen, is.null, logical(1))])
  if (is.null(icd) || !nrow(icd))
    stop("No survey_result.rds files found under ", OUT_BASE)
  write.csv(icd, IC_CSV, row.names = FALSE)
  message("Wrote ", basename(IC_CSV), " (", nrow(icd), " rows)")
}

# Failed bootstrap fits store an IQ-TREE error string in place of numbers, which
# makes the numeric columns load as character. Coerce and drop unusable rows.
for (col in keep) icd[[col]] <- suppressWarnings(as.numeric(icd[[col]]))
n_bad <- sum(is.na(icd$K))
if (n_bad) message("Dropping ", n_bad, " rows from failed fits")
icd <- icd[!is.na(icd$K), ]
icd$family <- factor(icd$family, levels = c("+R", "+H", "+T"))

parse_scen <- function(sc) {
  p <- strsplit(sc, "_")[[1]]
  list(true_type = FAM_MAP[p[1]],
       true_K    = as.integer(sub("^K", "", p[2])),
       tag       = p[3],
       seq_length = as.integer(sub("^L", "", p[4])))
}

# Diagnostic atlas: one page per scenario, so every page keeps its title.
plot_scenario <- function(sc, ic) {
  d <- icd[icd$scenario == sc & !is.na(icd[[ic]]), ]
  if (!nrow(d)) return(NULL)
  meta <- parse_scen(sc)
  emp  <- d[d$source == "empirical", ]
  sim  <- d[d$source == "sim", ]
  ggplot(mapping = aes(x = K, y = .data[[ic]], colour = family)) +
    geom_line(data = sim,
              aes(group = interaction(rep, line, family)),
              alpha = 0.08, linewidth = 0.25) +
    geom_line(data = emp,
              aes(group = interaction(rep, family)),
              linewidth = 0.7) +
    scale_colour_manual(values = PALETTE_FAMILY, name = "Fitted family",
                        drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(d$K))) +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(x = "K (number of classes)", y = ic,
         title = sprintf("%s  (true %s, K = %d, %s, L = %d) -- %s",
                         sc, meta$true_type, meta$true_K,
                         meta$tag, meta$seq_length, ic),
         subtitle = "Heavy lines: empirical fits. Faint lines: bootstrap replicates.") +
    theme_kpower(base_size = 9, grid = "both")
}

scenarios <- sort(unique(icd$scenario))

for (ic in IC_LEVELS) {
  path <- file.path(OUT_BASE, sprintf("ic_profiles_%s.pdf", ic))
  open_atlas(path, width = 6.5, height = 4.5)   # cairo_pdf: embeds fonts
  n <- 0L
  for (sc in scenarios) {
    p <- plot_scenario(sc, ic)
    if (!is.null(p)) { print(p); n <- n + 1L }
  }
  grDevices::dev.off()
  message(sprintf("%s -> ic_profiles_%s.pdf (%d scenario pages)", ic, ic, n))
}

message("\nDone.")
