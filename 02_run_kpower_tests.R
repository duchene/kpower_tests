# kpower simulation tests -- Step 2: Survey families on each alignment
#
# For each (scenario, rep) realisation, run kpower_survey() across +R, +H, +T
# and compute "ground-truth power" per IC: the proportion of bootstrap
# replicates in which the best (family, K) across all families matches the
# true (family, K) used to generate the data.
#
# Writes results/summary.csv with one row per (scenario, rep), with per-IC
# power and a `rep` column. Resumable: any (scenario, rep) already in
# summary.csv is skipped.
#
# Run 01_simulate_test_data.R first.

library(kpower)
library(ggplot2)

# --- Configuration ----------------------------------------------------------
# 48 scenarios x N_REPS realisations. With B=50, N_REPS=10, K_MAX=5,
# expect ~7-8 days of runtime. Resumable; safe to stop.
K_MAX      <- 5        # evaluate K = 1..K_MAX
B          <- 50       # bootstrap replicates per family
N_CORES    <- 10       # parallel R workers (bootstrap refits)
THREADS    <- 2        # IQ-TREE threads per run
FIXED_TREE <- NULL     # NULL = --fast heuristic (IQ-TREE 3.0.1 ARM BioNJ bug)
FAST_TREES <- TRUE     # GTR+R --fast for MAST candidate trees (skip MFP)
MIX_TYPES  <- c("+R", "+H", "+T")

IQTREE <- Sys.which("iqtree2")
#if (!nzchar(IQTREE))
#  IQTREE <- path.expand("~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3")

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
ALIGN_BASE <- file.path(SCRIPT_DIR, "alignments")
OUT_BASE   <- file.path(SCRIPT_DIR, "results")
SUMMARY_CSV  <- file.path(OUT_BASE, "summary.csv")
FAILURES_LOG <- file.path(OUT_BASE, "failures.log")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

#' Append a one-line entry to the failures log
log_failure <- function(label, msg) {
  cat(sprintf("[%s]  %s  --  %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label, msg),
      file = FAILURES_LOG, append = TRUE)
}

check_iqtree(IQTREE)


# --- Helpers ----------------------------------------------------------------

#' Compute fraction of bootstrap replicates whose minimum IC across all
#' (family, K) combinations matches (true_type, true_K).
ground_truth_power <- function(families, true_type, true_K) {
  parts <- lapply(names(families), function(mt) {
    fam <- families[[mt]]
    if (is.null(fam) || is.null(fam$sim_ic)) return(NULL)
    df <- fam$sim_ic
    df$mix_type <- mt
    df[, c("replicate", "mix_type", "K", "AIC", "AICc", "BIC")]
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0)
    return(c(AIC = NA_real_, AICc = NA_real_, BIC = NA_real_))
  long <- do.call(rbind, parts)

  reps <- unique(long$replicate)
  vapply(c("AIC", "AICc", "BIC"), function(crit) {
    correct <- vapply(reps, function(r) {
      sub <- long[long$replicate == r & is.finite(long[[crit]]), ]
      if (nrow(sub) == 0) return(NA)
      i <- which.min(sub[[crit]])
      isTRUE(sub$mix_type[i] == true_type) && isTRUE(sub$K[i] == true_K)
    }, logical(1))
    mean(correct, na.rm = TRUE)
  }, numeric(1))
}

#' Parse sim_params.txt into a named list
read_params <- function(path) {
  lines <- readLines(path)
  keys  <- sub(":.*", "", lines)
  vals  <- sub("^[^:]+:\\s*", "", lines)
  setNames(as.list(vals), keys)
}


# --- Discover (scenario, rep) tuples ----------------------------------------
# alignments/<scenario>/rep_NN/sim_params.txt
rep_dirs <- Sys.glob(file.path(ALIGN_BASE, "*", "rep_*"))
rep_dirs <- rep_dirs[file.exists(file.path(rep_dirs, "sim_params.txt"))]
if (length(rep_dirs) == 0)
  stop("No rep dirs found under ", ALIGN_BASE,
       ". Run 01_simulate_test_data.R first.")

# Build a tidy index: scenario, rep, dir. Sort by rep first then scenario
# so an interrupted run yields complete "passes" (all scenarios at rep 1
# before any at rep 2).
scenario  <- basename(dirname(rep_dirs))
rep       <- as.integer(sub("^rep_", "", basename(rep_dirs)))
idx       <- order(rep, scenario)
rep_dirs  <- rep_dirs[idx]
scenario  <- scenario[idx]
rep       <- rep[idx]


# --- Resume support ---------------------------------------------------------
existing <- if (file.exists(SUMMARY_CSV))
  read.csv(SUMMARY_CSV, stringsAsFactors = FALSE) else NULL
done_keys <- if (!is.null(existing) && "rep" %in% names(existing))
  paste(existing$scenario, existing$rep, sep = "::") else character()

todo_keys <- paste(scenario, rep, sep = "::")
todo_mask <- !(todo_keys %in% done_keys)
message(sprintf("Found %d (scenario, rep) tuples; %d already done; %d to run.",
                length(rep_dirs), sum(!todo_mask), sum(todo_mask)))

todo_dirs <- rep_dirs[todo_mask]
todo_scen <- scenario[todo_mask]
todo_rep  <- rep[todo_mask]


# --- Main loop --------------------------------------------------------------
for (i in seq_along(todo_dirs)) {
  sc_dir <- todo_dirs[i]
  label  <- todo_scen[i]
  r      <- todo_rep[i]
  params <- read_params(file.path(sc_dir, "sim_params.txt"))
  true_K    <- as.integer(params$K)
  true_type <- params$type
  tag       <- params$tag
  len       <- as.integer(params$length)

  align_file <- file.path(sc_dir, "sim.phy")
  if (!file.exists(align_file)) {
    warning("No sim.phy for ", label, " rep=", r, " -- skipping."); next
  }

  message(sprintf("\n[%d/%d] %s rep=%02d  (true: %s K=%d, %s, L=%d)",
                  i, length(todo_dirs), label, r,
                  true_type, true_K, tag, len))

  res_dir <- file.path(OUT_BASE, label, sprintf("rep_%02d", r))
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  surv <- tryCatch(
    kpower_survey(
      alignment  = align_file,
      K_max      = K_MAX,
      mix_types  = MIX_TYPES,
      ic         = "BIC",
      fixed_tree = FIXED_TREE,
      fast_trees = FAST_TREES,
      B          = B,
      outdir     = file.path(res_dir, "survey_runs"),
      iqtree_bin = IQTREE,
      n_cores    = N_CORES,
      threads    = THREADS
    ),
    error = function(e) {
      msg <- paste("survey crashed:", conditionMessage(e))
      message("!!! FAILED ", label, " rep=", r, " -- ", msg)
      log_failure(paste0(label, "/rep_", sprintf("%02d", r)), msg)
      NULL
    }
  )
  if (is.null(surv)) next

  null_fams <- names(surv$families)[vapply(surv$families, is.null,
                                           logical(1))]
  if (length(null_fams) > 0) {
    msg <- paste("family(ies) returned NULL:",
                 paste(null_fams, collapse = ", "))
    message("!!! FAMILY FAIL ", label, " rep=", r, " -- ", msg)
    log_failure(paste0(label, "/rep_", sprintf("%02d", r)), msg)
  }

  saveRDS(surv, file.path(res_dir, "survey_result.rds"))

  for (mt in names(surv$plots)) {
    mt_safe <- gsub("[+*]", "", mt)
    tryCatch(
      ggsave(file.path(res_dir, paste0("ic_profile_", mt_safe, ".pdf")),
             plot = surv$plots[[mt]], width = 7, height = 5),
      error = function(e) warning("Plot save failed for ", mt, ": ",
                                  e$message)
    )
  }

  comp <- surv$comparison
  best <- surv$best
  gtp  <- ground_truth_power(surv$families, true_type, true_K)

  row <- data.frame(
    scenario       = label,
    rep            = r,
    true_type      = true_type,
    true_K         = true_K,
    tag            = tag,
    seq_length     = len,
    best_mix_type  = best$mix_type,
    best_K         = best$K,
    best_BIC       = round(best$ic_value, 1),
    correct_family = (best$mix_type == true_type),
    correct_K      = (best$K == true_K),
    power_AIC      = round(unname(gtp["AIC"]),  4),
    power_AICc     = round(unname(gtp["AICc"]), 4),
    power_BIC      = round(unname(gtp["BIC"]),  4),
    elapsed_sec    = if (!is.null(surv$elapsed)) round(surv$elapsed, 1)
                     else NA_real_,
    stringsAsFactors = FALSE
  )

  for (mt in MIX_TYPES) {
    mt_safe <- gsub("[+*]", "", mt)
    cr <- comp[comp$mix_type == mt, ]
    row[[paste0(mt_safe, "_K")]]   <- if (nrow(cr) == 1) cr$K_best else NA
    row[[paste0(mt_safe, "_BIC")]] <- if (nrow(cr) == 1) round(cr$BIC, 1)
                                      else NA
  }

  message(sprintf("  best=%s K=%d | gt-power AIC=%.2f AICc=%.2f BIC=%.2f%s",
                  best$mix_type, best$K,
                  gtp["AIC"], gtp["AICc"], gtp["BIC"],
                  if (!is.null(surv$elapsed))
                    sprintf("  | elapsed %.0fs", surv$elapsed) else ""))

  write.table(
    row, SUMMARY_CSV,
    sep = ",", row.names = FALSE,
    col.names = !file.exists(SUMMARY_CSV),
    append    = file.exists(SUMMARY_CSV)
  )
}

message("\nResults: ", SUMMARY_CSV)
