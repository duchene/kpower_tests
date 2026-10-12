# kpower simulation tests -- Step 2d: merge the rerun summaries into results/
#
# The +T reruns write to their own directories so results/summary.csv is never
# edited in place by a running job:
#
#   results_T_multinom/  02b -- T_* rows at L <= 1000  (done, 120 rows)
#   results_T_stale/     02c -- the 97 remaining stale rows
#
# This folds them into results/summary.csv, newest row per (scenario, rep)
# wins, and keeps a timestamped backup. Column sets differ between jobs (the
# rerun rows carry <fam>_Kmaxeff, <fam>_ceiling, T_rep_dropped,
# T_rep_below_Kb), so the union is taken and missing cells are NA rather than
# dropping columns to the intersection.
#
#   Rscript 02d_merge_summaries.R           # report what would change
#   Rscript 02d_merge_summaries.R --write   # apply it
#
# Note the +T reruns are not comparable to the +T rows they replace: those
# predate both 65590ef (per-replicate candidate trees) and 5e41937
# (per-replicate class counts). Rows that survive the merge unreplaced are
# ones where +T chose K = 1, whose bootstrap never used the MAST simulation
# path at all.

write_out <- "--write" %in% commandArgs(TRUE)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
TARGET <- file.path(SCRIPT_DIR, "results", "summary.csv")
SOURCES <- file.path(SCRIPT_DIR,
                     c("results_T_multinom", "results_T_stale"),
                     "summary.csv")

if (!file.exists(TARGET)) stop("No ", TARGET)

read_summary <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df$source <- basename(dirname(path))
  df
}

target <- read_summary(TARGET)
srcs   <- lapply(SOURCES[file.exists(SOURCES)], read_summary)
if (length(srcs) == 0) stop("No rerun summaries found alongside ", TARGET)

# Union of columns, in target-first order
all_cols <- unique(c(names(target), unlist(lapply(srcs, names))))
pad <- function(df) {
  for (cl in setdiff(all_cols, names(df))) df[[cl]] <- NA
  df[, all_cols, drop = FALSE]
}
target <- pad(target)
srcs   <- lapply(srcs, pad)

new    <- do.call(rbind, srcs)
key    <- function(df) paste(df$scenario, df$rep)
dupes  <- duplicated(key(new))
if (any(dupes)) {
  warning(sum(dupes), " (scenario, rep) row(s) appear in more than one ",
          "rerun directory; keeping the first.")
  new <- new[!dupes, ]
}

replaced <- key(target) %in% key(new)
added    <- !(key(new) %in% key(target))

message(sprintf("target      : %d rows, %d columns", nrow(target),
                ncol(target) - 1L))
message(sprintf("rerun rows  : %d (%s)", nrow(new),
                paste(unique(new$source), collapse = ", ")))
message(sprintf("replacing   : %d rows", sum(replaced)))
message(sprintf("adding      : %d rows", sum(added)))

merged <- rbind(target[!replaced, ], new)
merged <- merged[order(merged$seq_length, merged$scenario, merged$rep), ]

still_stale <- with(merged, !is.na(T_K) & T_K >= 2 &
                    (is.na(source) | source == "results"))
message(sprintf("still stale after merge (T_K >= 2, original +T bootstrap): %d",
                sum(still_stale)))
if (any(still_stale))
  print(unique(merged[still_stale, c("scenario")]))

if (!write_out) {
  message("\nDry run. Re-run with --write to apply.")
  quit(save = "no")
}

backup <- sub("\\.csv$", format(Sys.time(), "_%Y%m%d_%H%M%S.csv"), TARGET)
file.copy(TARGET, backup)
message("backup      : ", backup)
merged$source <- NULL
write.csv(merged, TARGET, row.names = FALSE)
message("wrote       : ", TARGET, " (", nrow(merged), " rows)")
message("\nCaches built from summary.csv are now stale; rebuild with ",
        "REBUILD_CORRECT_BY_IC=1 Rscript 03_plot_power.R (likewise ",
        "REBUILD_TREE_METRICS / REBUILD_IC_TABLES for 05 and 06).")
