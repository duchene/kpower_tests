#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# kpower filovirus test — Step 2: BIC support for K across model families
#
# Runs kpower_survey() separately on each clade-restricted alignment built by
# 01_prep_alignment.R, across the three requested families:
#   +R  FreeRate  — rate heterogeneity across sites (RHAS)
#   +H  GHOST     — rate + branch-length classes (heterotachy)
#   +T  MAST      — mixture across sites and trees
#
# For each family it fits K = 1..K_MAX to the empirical data (BIC profile),
# then runs a parametric bootstrap (B replicates from the K_best model,
# refitting all K on each) to estimate how reliably BIC recovers K_best.
#
# Everything is written under test_data/ — nothing outside it is touched.
#
#   out: <clade>_kpower_run/        IQ-TREE working files
#        <clade>_survey.rds         full kpower_survey object
#        <clade>_ic_profile_*.pdf   kpower's per-family bootstrap figures
#        bic_profiles.csv           empirical lnL/df/AIC/AICc/BIC, both clades
#        family_comparison.csv      K_best + power under all three ICs
#        bic_support.pdf            empirical dBIC vs K, faceted by clade
#        survey_summary.txt         printed summary
#
# Usage:  Rscript 02_run_kpower_survey.R 2>&1 | tee survey_run.log
# ---------------------------------------------------------------------------

# The system-library kpower is an older build without kpower_survey(), so use
# the GitHub HEAD installed into test_data/Rlib.
HERE <- dirname(getwd())
.libPaths(c(file.path(dirname(HERE), "Rlib"), .libPaths()))  # Look in parent for Rlib

suppressPackageStartupMessages({
  library(kpower)
  library(ggplot2)
})
stopifnot("kpower_survey" %in% getNamespaceExports("kpower"))

# --- Configuration ---------------------------------------------------------
DATASETS <- list(
  list(prefix = "ebola",   label = "Orthoebolavirus"),
  list(prefix = "marburg", label = "Orthomarburgvirus + Mengla/Dehong")
)

IQTREE <- Sys.which("iqtree3")
if (!nzchar(IQTREE))
  IQTREE <- path.expand("~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3")

K_MAX      <- 5           # evaluate K = 1..5 (per-family override below)
IC         <- "BIC"       # primary criterion (AIC/AICc also reported)
B          <- 100L        # bootstrap replicates

# MAST is far more expensive than the other families and its K=5 fits
# dominate total runtime. At B=10 both clades put K_best well below 5
# (ebola K=1, marburg K=3, with K=5 worse by 1148 and 332 BIC units), so
# capping +T at 4 keeps K_best off the boundary while dropping the most
# expensive fits.
K_MAX_BY_FAMILY <- c("+R" = 5L, "+H" = 5L, "+T" = 4L)
MIX_TYPES  <- c("+R", "+H", "+T")
BASE_MODEL <- "GTR"
FIXED_TREE <- NULL        # NULL = --fast heuristic search, as in kpower_tests
FAST_TREES <- TRUE        # --fast for MAST window trees (skip ModelFinder)
N_CORES    <- 6L          # parallel R workers
THREADS    <- 2L          # IQ-TREE threads per run (12 of 16 cores)
TIMEOUT    <- 14400L      # 4 h per IQ-TREE run
SEED       <- 1L

FAMILY_LABELS <- c("+R" = "FreeRate (RHAS)", "+H" = "GHOST", "+T" = "MAST")

check_iqtree(IQTREE)

run_one <- function(ds) {
  aln <- file.path(HERE, paste0(ds$prefix, "_gb_strict.fasta"))
  if (!file.exists(aln))
    stop("Alignment not found: ", aln,
         "\nRun 01_prep_alignment.R then 01b_gblocks_filter.R first.")

  lines <- readLines(aln, warn = FALSE)
  ntax  <- sum(startsWith(lines, ">"))
  nsite <- nchar(lines[2])

  message(sprintf("\n\n#########################################################"))
  message(sprintf("### %s: %d taxa x %d sites", ds$label, ntax, nsite))
  message(sprintf("#########################################################"))

  outdir <- file.path(HERE, paste0(ds$prefix, "_kpower_run"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # One kpower_survey() call per family, checkpointed to its own RDS. A long
  # run killed part-way (MAST can take hours) then resumes at the family it
  # died on instead of recomputing everything.
  t0 <- Sys.time()
  per_family <- list()
  for (mt in MIX_TYPES) {
    ck <- file.path(HERE, sprintf("%s_%s_B%d_K%d_survey.rds", ds$prefix,
                                  gsub("[+*]", "", mt), B, K_MAX_BY_FAMILY[[mt]]))
    if (file.exists(ck)) {
      message("\n[resume] ", ds$label, " ", mt, " - using ", basename(ck))
      per_family[[mt]] <- readRDS(ck)
      next
    }
    one <- kpower_survey(
      alignment  = aln,
      K_max      = K_MAX_BY_FAMILY[[mt]],
      mix_types  = mt,
      base_model = BASE_MODEL,
      ic         = IC,
      fixed_tree = FIXED_TREE,
      fast_trees = FAST_TREES,
      B          = B,
      seed       = SEED,
      outdir     = outdir,
      iqtree_bin = IQTREE,
      n_cores    = N_CORES,
      threads    = THREADS,
      timeout    = TIMEOUT
    )
    saveRDS(one, ck)
    per_family[[mt]] <- one
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  message(sprintf("\n%s finished in %.2f h", ds$label, elapsed))

  # Stitch the per-family surveys back into one kpower_survey object.
  # Build the merged lists by hand: c() on a named list of named lists makes
  # compound names ("+R.+R"), which then miss on lookup by mix_type.
  surv <- per_family[[1]]
  surv$families <- stats::setNames(
    lapply(per_family, function(x) x$families[[1]]), names(per_family))
  surv$plots <- local({
    out <- list()
    for (mt in names(per_family)) {
      p <- per_family[[mt]]$plots
      if (length(p)) out[[mt]] <- p[[1]]
    }
    out
  })
  surv$comparison <- do.call(rbind, lapply(per_family, `[[`, "comparison"))
  surv$best <- local({
    ok <- surv$comparison[!is.na(surv$comparison[[IC]]), , drop = FALSE]
    if (!nrow(ok)) list(mix_type = NA, K = NA, ic_value = NA)
    else {
      i <- which.min(ok[[IC]])
      list(mix_type = ok$mix_type[i], K = ok$K_best[i], ic_value = ok[[IC]][i])
    }
  })

  saveRDS(surv, file.path(HERE, paste0(ds$prefix, "_survey.rds")))

  for (mt in names(surv$plots))
    ggsave(file.path(HERE, sprintf("%s_ic_profile_%s.pdf", ds$prefix,
                                   gsub("[+*]", "", mt))),
           surv$plots[[mt]], width = 7, height = 5)

  prof <- do.call(rbind, lapply(names(surv$families), function(mt) {
    fam <- surv$families[[mt]]
    if (is.null(fam)) return(NULL)
    e <- fam$empirical
    data.frame(clade = ds$prefix, clade_label = ds$label,
               ntax = ntax, nsite = nsite,
               family = mt, family_label = FAMILY_LABELS[[mt]],
               K = e$K, df = e$df, lnL = e$lnL,
               AIC = e$AIC, AICc = e$AICc, BIC = e$BIC,
               K_best = fam$K_best, stringsAsFactors = FALSE)
  }))

  comp <- surv$comparison
  comp$clade <- ds$prefix

  list(surv = surv, profiles = prof, comparison = comp,
       elapsed = elapsed, ntax = ntax, nsite = nsite, ds = ds)
}

results <- lapply(DATASETS, run_one)

# --- Combined tables --------------------------------------------------------
profiles <- do.call(rbind, lapply(results, `[[`, "profiles"))
if (is.null(profiles)) stop("All families failed in all clades.")

key <- paste(profiles$clade, profiles$family)
profiles$dBIC_within <- ave(profiles$BIC, key, FUN = function(x) x - min(x, na.rm = TRUE))
profiles$dBIC_clade  <- ave(profiles$BIC, profiles$clade,
                            FUN = function(x) x - min(x, na.rm = TRUE))

write.csv(profiles, file.path(HERE, "bic_profiles.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(results, `[[`, "comparison")),
          file.path(HERE, "family_comparison.csv"), row.names = FALSE)

# --- Combined figure --------------------------------------------------------
best_pts <- profiles[profiles$K == profiles$K_best, ]

p <- ggplot(profiles, aes(K, dBIC_clade, colour = family_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_point(data = best_pts, shape = 21, size = 5, stroke = 1.2,
             fill = NA, show.legend = FALSE) +
  facet_wrap(~ clade_label, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(profiles$K))) +
  labs(x = "Number of mixture categories (K)",
       y = expression(Delta * "BIC from best model in that clade"),
       colour = "Model family",
       title = "BIC support for K across mixture model families",
       subtitle = "Circles mark K_best per family") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(HERE, "bic_support.pdf"), p, width = 10, height = 5)

# --- Summary ----------------------------------------------------------------
sink(file.path(HERE, "survey_summary.txt"), split = TRUE)

cat("kpower survey - Filoviridae, by genus\n")
cat("  K range : 1 ..", K_MAX, "   B =", B, "   primary IC =", IC, "\n\n")

for (r in results) {
  cat("=====", r$ds$label, "-", r$ntax, "taxa x", r$nsite, "sites",
      sprintf("(%.2f h)", r$elapsed), "=====\n")

  pr <- r$profiles
  show <- pr[, c("family", "K", "df", "lnL", "BIC")]
  show$dBIC <- ave(pr$BIC, pr$family, FUN = function(x) x - min(x, na.rm = TRUE))
  show$lnL  <- round(show$lnL, 1)
  show$BIC  <- round(show$BIC, 1)
  show$dBIC <- round(show$dBIC, 1)
  print(show, row.names = FALSE)

  cat("\n")
  print(r$surv)

  cat("\n  Power under all three criteria:\n")
  for (mt in names(r$surv$families)) {
    fam <- r$surv$families[[mt]]
    if (is.null(fam)) { cat(sprintf("    %-3s failed\n", mt)); next }
    pa <- fam$power_all
    cat(sprintf("    %-3s  AIC: K=%d %5.1f%%   AICc: K=%d %5.1f%%   BIC: K=%d %5.1f%%\n",
                mt, pa$AIC$K_best, 100 * pa$AIC$power,
                pa$AICc$K_best, 100 * pa$AICc$power,
                pa$BIC$K_best, 100 * pa$BIC$power))
  }
  cat("\n\n")
}

sink()

message("\nWrote to ", HERE, ":")
message("  <clade>_survey.rds, <clade>_ic_profile_*.pdf")
message("  bic_profiles.csv, family_comparison.csv, bic_support.pdf, survey_summary.txt")
