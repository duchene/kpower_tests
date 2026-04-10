# kpower simulation tests — Step 2: Survey all model families on each alignment
#
# For each simulated alignment, runs kpower_survey() across all model families
# (+R, *R, +H, *H, +T, *T) to test whether:
#   1. The true generating model family is identified as best
#   2. The true K is recovered within that family
#   3. Simple processes aren't mistaken for complex ones and vice versa
#
# Run 01_simulate_test_data.R first.

library(kpower)
library(ggplot2)

# --- Configuration -----------------------------------------------------------
IQTREE  <- Sys.which("iqtree3")
if (!nzchar(IQTREE))
  IQTREE <- path.expand("~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3")

K_MAX      <- 4       # evaluate K = 1..K_MAX
B          <- 20      # bootstrap replicates per family (use 100+ for production)
IC         <- "BIC"
N_CORES    <- 4       # parallel R workers for bootstrap refits
THREADS    <- 1       # IQ-TREE threads per run
FIXED_TREE <- NULL    # NULL = --fast heuristic (IQ-TREE 3.0.1 ARM BioNJ bug)
FAST_TREES <- TRUE    # Use --fast for MAST candidate trees (skip MFP)

# Which model families to survey
MIX_TYPES <- c("+R", "+H", "+T")

ALIGN_BASE <- path.expand("~/Dropbox/Research/hiv_mast/kpower_tests/alignments")
OUT_BASE   <- path.expand("~/Dropbox/Research/hiv_mast/kpower_tests/results")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

check_iqtree(IQTREE)

# --- Discover scenarios ------------------------------------------------------
scenario_dirs <- sort(list.dirs(ALIGN_BASE, recursive = FALSE, full.names = TRUE))
scenario_dirs <- scenario_dirs[file.exists(file.path(scenario_dirs, "sim_params.txt"))]

if (length(scenario_dirs) == 0)
  stop("No scenarios found. Run 01_simulate_test_data.R first.")

message(sprintf("Found %d scenarios. Surveying %s on each (B=%d).\n",
                length(scenario_dirs), paste(MIX_TYPES, collapse = ", "), B))

summary_rows <- vector("list", length(scenario_dirs))

for (i in seq_along(scenario_dirs)) {
  sc_dir <- scenario_dirs[i]
  label  <- basename(sc_dir)

  # Parse sim_params.txt
  params <- readLines(file.path(sc_dir, "sim_params.txt"))
  plist  <- setNames(sub("^[^:]+:\\s*", "", params), sub(":.*", "", params))
  true_K    <- as.integer(plist[["K"]])
  true_type <- plist[["type"]]
  tag       <- plist[["tag"]]
  len       <- as.integer(plist[["length"]])

  align_file <- file.path(sc_dir, "sim.phy")
  if (!file.exists(align_file)) {
    warning("No sim.phy for ", label, " -- skipping.")
    next
  }

  message(sprintf("\n[%d/%d] %s  (true: %s K=%d, %s, L=%d)",
                  i, length(scenario_dirs), label,
                  true_type, true_K, tag, len))

  res_dir <- file.path(OUT_BASE, label)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  surv <- tryCatch(
    kpower_survey(
      alignment  = align_file,
      K_max      = K_MAX,
      mix_types  = MIX_TYPES,
      ic         = IC,
      fixed_tree = FIXED_TREE,
      fast_trees = FAST_TREES,
      B          = B,
      outdir     = file.path(res_dir, "survey_runs"),
      iqtree_bin = IQTREE,
      n_cores    = N_CORES,
      threads    = THREADS
    ),
    error = function(e) {
      warning("Survey failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(surv)) next

  saveRDS(surv, file.path(res_dir, "survey_result.rds"))

  # Save per-family plots (tryCatch: some may fail if IQ-TREE errors
  # produced non-numeric K values in the plot data)
  for (mt in names(surv$plots)) {
    mt_safe <- gsub("[+*]", "", mt)
    tryCatch(
      ggsave(
        filename = file.path(res_dir, paste0("ic_profile_", mt_safe, ".pdf")),
        plot     = surv$plots[[mt]],
        width    = 7, height = 5
      ),
      error = function(e) warning("Plot save failed for ", mt, ": ", e$message)
    )
  }

  # Extract comparison row
  comp <- surv$comparison
  best <- surv$best

  summary_rows[[i]] <- data.frame(
    scenario       = label,
    true_type      = true_type,
    true_K         = true_K,
    tag            = tag,
    seq_length     = len,
    best_mix_type  = best$mix_type,
    best_K         = best$K,
    best_ic        = round(best$ic_value, 1),
    correct_family = (best$mix_type == true_type),
    correct_K      = (best$K == true_K),
    stringsAsFactors = FALSE
  )

  # Append per-family detail columns
  for (mt in MIX_TYPES) {
    mt_safe <- gsub("[+*]", "", mt)
    row <- comp[comp$mix_type == mt, ]
    if (nrow(row) == 1) {
      summary_rows[[i]][[paste0(mt_safe, "_K")]]     <- row$K_best
      summary_rows[[i]][[paste0(mt_safe, "_BIC")]]    <- round(row$BIC, 1)
      summary_rows[[i]][[paste0(mt_safe, "_power")]]  <- round(row$power_BIC * 100, 1)
    } else {
      summary_rows[[i]][[paste0(mt_safe, "_K")]]     <- NA
      summary_rows[[i]][[paste0(mt_safe, "_BIC")]]    <- NA
      summary_rows[[i]][[paste0(mt_safe, "_power")]]  <- NA
    }
  }

  msg_parts <- vapply(MIX_TYPES, function(mt) {
    mt_safe <- gsub("[+*]", "", mt)
    sprintf("%s:K=%s", mt,
            ifelse(is.na(summary_rows[[i]][[paste0(mt_safe, "_K")]]),
                   "fail",
                   as.character(summary_rows[[i]][[paste0(mt_safe, "_K")]])))
  }, character(1))

  message(sprintf("  Best: %s K=%d | %s | %s",
                  best$mix_type, best$K,
                  if (best$mix_type == true_type) "FAMILY OK" else "FAMILY WRONG",
                  paste(msg_parts, collapse = "  ")))
}

# --- Summary table -----------------------------------------------------------
summary_tbl <- do.call(rbind, Filter(Negate(is.null), summary_rows))
message("\n\n===== SUMMARY =====")
print(summary_tbl[, c("scenario", "true_type", "true_K", "tag", "seq_length",
                       "best_mix_type", "best_K", "correct_family", "correct_K")])
write.csv(summary_tbl, file.path(OUT_BASE, "summary.csv"), row.names = FALSE)

message(sprintf("\nFamily correct: %d/%d (%.0f%%)",
                sum(summary_tbl$correct_family),
                nrow(summary_tbl),
                mean(summary_tbl$correct_family) * 100))
message(sprintf("K correct:      %d/%d (%.0f%%)",
                sum(summary_tbl$correct_K),
                nrow(summary_tbl),
                mean(summary_tbl$correct_K) * 100))

message("\nResults saved to: ", OUT_BASE)
