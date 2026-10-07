#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# kpower filovirus test — Step 1b: Gblocks block selection
#
# Replaces the per-column missing-data filter from 01_prep_alignment.R with
# Gblocks. The missing-data filter judges each column independently, so it
# happily keeps an isolated well-occupied column sitting inside a region that
# is not reliably aligned. Gblocks instead selects contiguous BLOCKS that are
# conserved and flanked by conserved positions, which is the right notion of
# "sensible region" here.
#
# Runs two parameter sets on each clade so they can be compared:
#
#   strict   Gblocks defaults: b2 = 85% of taxa, b3 = 8, b4 = 10, b5 = none
#            (no gaps tolerated anywhere in a retained block)
#   relaxed  Talavera & Castresana (2007) Syst Biol 56:564: b2 = b1 (50%),
#            b3 = 10, b4 = 5, b5 = half (gaps in up to half the taxa)
#
# Input is <clade>_align2.fasta — the realigned, taxon-pruned alignment,
# BEFORE the old missing-data column filter, so Gblocks sees the real thing.
#
#   out: <clade>_gb_<setting>.fasta / .phy   filtered alignments
#        <clade>_gb_<setting>.html           Gblocks' own block visualisation
#        gblocks_summary.txt                 comparison table
#
# Usage:  Rscript 01b_gblocks_filter.R
# ---------------------------------------------------------------------------

HERE    <- path.expand("~/Dropbox/Research/kpower/test_data")
GBLOCKS <- file.path(HERE, "tools", "bin", "Gblocks")
OUT_LOG <- file.path(HERE, "gblocks_summary.txt")

CLADES <- c("ebola", "marburg")
VALID  <- c("A", "C", "G", "T", "U")
MIN_COMP_SITES <- 500

log_lines <- character(0)
say <- function(...) { m <- paste0(...); message(m); log_lines <<- c(log_lines, m) }

read_fasta_full <- function(path) {
  lines   <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  ends    <- c(headers[-1] - 1L, length(lines))
  seqs <- vapply(seq_along(headers), function(i) {
    body <- lines[(headers[i] + 1L):ends[i]]
    paste(body[nzchar(trimws(body))], collapse = "")
  }, character(1))
  data.frame(acc = sub("^>\\s*", "", sub("\\s+.*$", "", lines[headers])),
             seq = toupper(gsub("[^A-Za-z-]", "", seqs)),
             stringsAsFactors = FALSE)
}

write_fasta <- function(acc, seq, path) {
  out <- character(2L * length(acc))
  out[c(TRUE, FALSE)] <- paste0(">", acc)
  out[c(FALSE, TRUE)] <- seq
  writeLines(out, path)
}

p_summary <- function(aln) {
  m <- do.call(rbind, strsplit(aln$seq, "", fixed = TRUE))
  cm <- matrix(NA_integer_, nrow(m), ncol(m))
  cm[m == "A"] <- 1L; cm[m == "C"] <- 2L
  cm[m == "G"] <- 3L; cm[m == "T" | m == "U"] <- 4L
  n <- nrow(cm); v <- c()
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    ok <- !is.na(cm[i, ]) & !is.na(cm[j, ])
    if (sum(ok) >= MIN_COMP_SITES) v <- c(v, sum(cm[i, ok] != cm[j, ok]) / sum(ok))
  }
  list(missing = mean(is.na(cm)), median_p = median(v), max_p = max(v))
}

# Gblocks writes its output next to the input, so stage a copy per setting.
run_gblocks <- function(clade, setting, params) {
  src <- file.path(HERE, paste0(clade, "_align2.fasta"))
  if (!file.exists(src)) stop("Missing ", src, " - run 01_prep_alignment.R first.")

  work <- file.path(HERE, paste0(clade, "_gb_", setting, ".work.fasta"))
  file.copy(src, work, overwrite = TRUE)

  out <- system2(GBLOCKS, c(shQuote(work), params), stdout = TRUE, stderr = TRUE)
  nblock <- grep("Gblocks alignment", out, value = TRUE)

  gb_fa <- paste0(work, "-gb.fa")
  if (!file.exists(gb_fa)) {
    say("  GBLOCKS FAILED for ", clade, " / ", setting)
    say(paste("   ", out[seq_len(min(12, length(out)))], collapse = "\n"))
    return(NULL)
  }

  aln <- read_fasta_full(gb_fa)
  dest_fa  <- file.path(HERE, paste0(clade, "_gb_", setting, ".fasta"))
  dest_phy <- file.path(HERE, paste0(clade, "_gb_", setting, ".phy"))
  write_fasta(aln$acc, aln$seq, dest_fa)
  writeLines(c(paste(nrow(aln), nchar(aln$seq[1])), paste(aln$acc, aln$seq)), dest_phy)
  file.rename(paste0(work, "-gb.html"),
              file.path(HERE, paste0(clade, "_gb_", setting, ".html")))
  file.remove(work, gb_fa)

  st <- p_summary(aln)
  say(sprintf("  %-7s  %5d sites  %s", setting, nchar(aln$seq[1]),
              sub(".*positions", "", nblock)))
  say(sprintf("           missing %.2f%%   p-dist median %.3f  max %.3f",
              100 * st$missing, st$median_p, st$max_p))

  data.frame(clade = clade, setting = setting, taxa = nrow(aln),
             sites = nchar(aln$seq[1]), missing = st$missing,
             median_p = st$median_p, max_p = st$max_p,
             stringsAsFactors = FALSE)
}

say("=== Gblocks block selection: ", format(Sys.time()), " ===")
if (!file.exists(GBLOCKS)) stop("Gblocks not found at ", GBLOCKS)

rows <- list()
for (clade in CLADES) {
  src <- file.path(HERE, paste0(clade, "_align2.fasta"))
  a0  <- read_fasta_full(src)
  n   <- nrow(a0)
  b1  <- floor(n / 2) + 1          # min seqs for a conserved position
  b2_strict  <- ceiling(0.85 * n)  # Gblocks default flank stringency
  b2_relaxed <- b1                 # Talavera & Castresana relaxed

  say("\n", clade, ": ", n, " taxa x ", nchar(a0$seq[1]), " positions (input)")
  say("  b1 = ", b1, "  b2 strict = ", b2_strict, "  b2 relaxed = ", b2_relaxed)

  rows[[length(rows) + 1]] <- run_gblocks(clade, "strict", c(
    "-t=d", paste0("-b1=", b1), paste0("-b2=", b2_strict),
    "-b3=8", "-b4=10", "-b5=n", "-e=-gb", "-p=y"))

  rows[[length(rows) + 1]] <- run_gblocks(clade, "relaxed", c(
    "-t=d", paste0("-b1=", b1), paste0("-b2=", b2_relaxed),
    "-b3=10", "-b4=5", "-b5=h", "-e=-gb", "-p=y"))
}

res <- do.call(rbind, rows)
say("\n\n=== Comparison ===")
say(sprintf("  %-8s %-8s %5s %7s %9s %9s %8s", "clade", "setting", "taxa",
            "sites", "missing", "median_p", "max_p"))
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  say(sprintf("  %-8s %-8s %5d %7d %8.2f%% %9.3f %8.3f",
              r$clade, r$setting, r$taxa, r$sites, 100 * r$missing,
              r$median_p, r$max_p))
}

say("\nOpen the .html files to see which blocks were selected.")
say("Then pick a setting and point 02_run_kpower_survey.R at")
say("  <clade>_gb_<setting>.fasta")

writeLines(log_lines, OUT_LOG)
