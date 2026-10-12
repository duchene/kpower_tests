#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# kpower filovirus test — Step 1: build two clade-restricted alignments
#
# v3. The pass-1 distance matrix (v2) showed Filoviridae cannot be a single
# nucleotide alignment at realistic distances:
#
#   Orthoebolavirus lineages      0.30-0.40 between each other   OK
#   Marburg + Mengla + Dehong     0.38-0.43 between each other   OK
#   ebolavirus <-> marburgvirus   0.51-0.52                      saturated
#   Lloviu <-> everything         0.46-0.53                      no home
#
# So each genus is prepared and aligned separately. Lloviu is excluded from
# both: it is an outgroup to each at a distance nucleotide data cannot carry.
# Fish/reptile filoviruses are excluded for the same reason, more severely.
#
# Taxon counts are kept deliberately small. Sequences are retained only when
# they differ from everything already kept by at least MIN_DIFFS substitutions,
# so a lineage with little real variation contributes few taxa rather than many
# near-identical ones. Padding the set with duplicates adds branches of near-
# zero length, which is exactly what makes GHOST and MAST fits degenerate.
#
#   out: <clade>_subsample.fasta   initial subsample (unaligned)
#        <clade>_align1.fasta      first MAFFT pass
#        <clade>_align2.fasta      realignment of retained taxa
#        <clade>_filtered.fasta    site-filtered final alignment -> kpower
#        <clade>_filtered.phy
#        <clade>_distances.csv     final pairwise p-distance matrix
#        prep_summary.txt          full log of every decision, both clades
#
# Usage:  Rscript 01_prep_alignment.R
# ---------------------------------------------------------------------------

HERE    <- path.expand("~/Dropbox/Research/kpower/test_data")
RAW     <- file.path(HERE, "ncbi_virus_filoviridae.fasta")
OUT_LOG <- file.path(HERE, "prep_summary.txt")

MAFFT      <- Sys.which("mafft")
THREADS    <- 6
MAFFT_ARGS <- c("--retree", "2", "--maxiterate", "2", "--reorder")

MIN_LEN_FRAC  <- 0.95
DROP_PATTERNS <- "mutant|unverified|synthetic construct|vaccine|recombinant"

MAX_PDIST      <- 0.45   # saturation gate
MIN_DIFFS      <- 50     # a retained taxon must differ by >= this many subs
MIN_COMP_SITES <- 500

MAX_COL_MISSING   <- 0.50
MAX_TAXON_MISSING <- 0.50
VALID <- c("A", "C", "G", "T", "U")

# Two clades, each internally alignable. Caps are set so that after
# deduplication each lands well under 50 taxa.
CLADES <- list(
  list(
    prefix = "ebola",
    label  = "Orthoebolavirus",
    cap    = 15,
    groups = c(
      Sudan      = "sudan",
      Bundibugyo = "bundibugyo",
      Reston     = "reston",
      Bombali    = "bombali",
      TaiForest  = "tai forest|ta.. forest|cote d|ivory coast",
      Zaire      = "zaire|ebola virus|ebov"
    )
  ),
  list(
    prefix = "marburg",
    label  = "Orthomarburgvirus + Mengla/Dehong",
    cap    = 20,
    groups = c(
      Mengla  = "mengla|dianlovirus",
      Dehong  = "dehong|wanding",
      Marburg = "marburg|lake victoria|ravn"
    )
  )
)

# Excluded from both clades (see header).
EXCLUDE <- paste0("wenling|tapajos|oberland|fiwi|kander|huangjiao|xilang|",
                  "frogfish|thamnaconus|lloviu")

log_lines <- character(0)
say <- function(...) { msg <- paste0(...); message(msg); log_lines <<- c(log_lines, msg) }

# --- I/O --------------------------------------------------------------------
read_fasta_full <- function(path) {
  lines   <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  ends    <- c(headers[-1] - 1L, length(lines))
  seqs <- vapply(seq_along(headers), function(i) {
    body <- lines[(headers[i] + 1L):ends[i]]
    paste(body[nzchar(trimws(body))], collapse = "")
  }, character(1))
  desc <- sub("^>", "", lines[headers])
  data.frame(acc = sub("\\s+.*$", "", desc), desc = desc,
             seq = toupper(gsub("\\s", "", seqs)), stringsAsFactors = FALSE)
}

write_fasta <- function(acc, seq, path) {
  out <- character(2L * length(acc))
  out[c(TRUE, FALSE)] <- paste0(">", acc)
  out[c(FALSE, TRUE)] <- seq
  writeLines(out, path)
}

run_mafft <- function(infile, outfile, tag) {
  say("\nMAFFT (", paste(MAFFT_ARGS, collapse = " "), ", ", THREADS,
      " threads) - ", tag, " ...")
  t0 <- Sys.time()
  st <- system2(MAFFT, c(MAFFT_ARGS, "--thread", THREADS, shQuote(infile)),
                stdout = outfile, stderr = file.path(HERE, "mafft.log"))
  if (st != 0) stop("MAFFT failed (status ", st, ") - see mafft.log")
  say("  done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  read_fasta_full(outfile)
}

# --- Distances --------------------------------------------------------------
code_matrix <- function(acc, seqs) {
  m <- do.call(rbind, strsplit(seqs, "", fixed = TRUE))
  out <- matrix(NA_integer_, nrow(m), ncol(m), dimnames = list(acc, NULL))
  out[m == "A"] <- 1L; out[m == "C"] <- 2L
  out[m == "G"] <- 3L; out[m == "T" | m == "U"] <- 4L
  out
}

# Returns both the proportion and the raw substitution count: the saturation
# gate is a proportion, the duplicate gate is an absolute number of changes.
distances <- function(cm) {
  n <- nrow(cm)
  P <- D <- matrix(NA_real_, n, n, dimnames = list(rownames(cm), rownames(cm)))
  diag(P) <- 0; diag(D) <- 0
  for (i in seq_len(n - 1L)) {
    a <- cm[i, ]
    for (j in (i + 1L):n) {
      b  <- cm[j, ]
      ok <- !is.na(a) & !is.na(b)
      ns <- sum(ok)
      if (ns >= MIN_COMP_SITES) {
        nd <- sum(a[ok] != b[ok])
        P[i, j] <- P[j, i] <- nd / ns
        D[i, j] <- D[j, i] <- nd
      }
    }
  }
  list(p = P, diffs = D)
}

summarise_distances <- function(P, groups, header) {
  off <- P[upper.tri(P)]
  say("\n", header)
  say(sprintf("  pairwise p-distance: median %.3f  mean %.3f  max %.3f  (%d pairs)",
              median(off, na.rm = TRUE), mean(off, na.rm = TRUE),
              max(off, na.rm = TRUE), sum(!is.na(off))))
  say(sprintf("  pairs above %.2f: %d  (%.1f%%)", MAX_PDIST,
              sum(off > MAX_PDIST, na.rm = TRUE),
              100 * mean(off > MAX_PDIST, na.rm = TRUE)))
  gl <- sort(unique(groups))
  say("  median p-distance between lineages:")
  say(paste0("    ", sprintf("%-11s", ""), paste(sprintf("%8s", gl), collapse = "")))
  for (g1 in gl) {
    cells <- vapply(gl, function(g2) {
      v <- P[groups == g1, groups == g2, drop = FALSE]
      v <- if (g1 == g2) v[upper.tri(v)] else as.vector(v)
      if (!length(v) || all(is.na(v))) "       -" else sprintf("%8.3f", median(v, na.rm = TRUE))
    }, character(1))
    say(paste0("    ", sprintf("%-11s", g1), paste(cells, collapse = "")))
  }
}

# --- One clade --------------------------------------------------------------
prep_clade <- function(clade, pool) {
  px <- clade$prefix
  f <- function(suffix) file.path(HERE, paste0(px, suffix))

  say("\n\n###########################################################")
  say("### ", clade$label, "  (", px, ")")
  say("###########################################################")

  dat <- pool
  dat$group <- NA_character_
  for (g in names(clade$groups))
    dat$group[is.na(dat$group) &
                grepl(clade$groups[[g]], dat$desc, ignore.case = TRUE)] <- g
  dat <- dat[!is.na(dat$group), , drop = FALSE]
  say("\nPool for this clade: ", nrow(dat), " sequences")

  sel <- do.call(rbind, lapply(split(dat, dat$group), function(sub) {
    full <- sub[sub$len >= MIN_LEN_FRAC * max(sub$len), , drop = FALSE]
    if (!nrow(full)) full <- sub
    if (nrow(full) <= clade$cap) return(full)
    full <- full[order(full$acc), , drop = FALSE]
    full[unique(round(seq(1, nrow(full), length.out = clade$cap))), , drop = FALSE]
  }))
  sel <- sel[order(sel$group, sel$acc), , drop = FALSE]
  med_raw <- median(sel$len)

  say("Subsampled ", nrow(sel), " near-complete genomes:")
  for (g in names(sort(table(sel$group), decreasing = TRUE)))
    say(sprintf("  %-11s %3d   (%d-%d bp)", g, sum(sel$group == g),
                min(sel$len[sel$group == g]), max(sel$len[sel$group == g])))
  write_fasta(sel$acc, sel$seq, f("_subsample.fasta"))

  # Pass 1
  a1 <- run_mafft(f("_subsample.fasta"), f("_align1.fasta"), "pass 1")
  say("  aligned length ", nchar(a1$seq[1]), " vs median input ", med_raw,
      sprintf("  (expansion %.2fx)", nchar(a1$seq[1]) / med_raw))
  g1 <- sel$group[match(a1$acc, sel$acc)]
  d1 <- distances(code_matrix(a1$acc, a1$seq))
  summarise_distances(d1$p, g1, "Distances after pass 1:")

  keep <- rep(TRUE, nrow(a1)); names(keep) <- a1$acc

  # Saturation prune
  repeat {
    idx  <- which(keep)
    sub  <- d1$p[idx, idx, drop = FALSE]
    viol <- sub > MAX_PDIST; viol[is.na(viol)] <- FALSE
    if (!any(viol)) break
    nv    <- rowSums(viol)
    worst <- which(nv == max(nv))
    if (length(worst) > 1)
      worst <- worst[which.max(rowMeans(sub[worst, , drop = FALSE], na.rm = TRUE))]
    acc <- rownames(sub)[worst]
    say(sprintf("  prune saturated: %-12s (%s, %d pairs > %.2f)",
                acc, g1[match(acc, a1$acc)], max(nv), MAX_PDIST))
    keep[acc] <- FALSE
  }
  say("Saturation prune removed ", sum(!keep), " taxa")

  # Duplicate prune: keep only sequences that differ by >= MIN_DIFFS subs
  miss <- vapply(a1$seq, function(s) {
    ch <- strsplit(s, "", fixed = TRUE)[[1]]; mean(!(ch %in% VALID))
  }, numeric(1))
  names(miss) <- a1$acc

  n_dup <- 0
  repeat {
    idx <- which(keep)
    sub <- d1$diffs[idx, idx, drop = FALSE]
    sub[lower.tri(sub, diag = TRUE)] <- NA
    if (!any(sub < MIN_DIFFS, na.rm = TRUE)) break
    w  <- which(sub == min(sub, na.rm = TRUE), arr.ind = TRUE)[1, ]
    pr <- rownames(sub)[c(w[1], w[2])]
    keep[pr[which.max(miss[pr])]] <- FALSE
    n_dup <- n_dup + 1
  }
  say("Duplicate prune removed ", n_dup, " taxa (< ", MIN_DIFFS, " substitutions apart)")

  retained <- a1$acc[keep]
  gk <- g1[keep]
  say("\nRetained ", length(retained), " taxa:")
  for (g in names(sort(table(gk), decreasing = TRUE)))
    say(sprintf("  %-11s %3d", g, sum(gk == g)))
  if (any(!(unique(g1) %in% gk)))
    say("  lineages dropped entirely: ", paste(setdiff(unique(g1), gk), collapse = ", "))

  # Pass 2
  write_fasta(retained, sel$seq[match(retained, sel$acc)], f("_retained.fasta"))
  a2 <- run_mafft(f("_retained.fasta"), f("_align2.fasta"), "pass 2, retained taxa only")
  say("  aligned length ", nchar(a2$seq[1]), " vs median input ", med_raw,
      sprintf("  (expansion %.2fx)", nchar(a2$seq[1]) / med_raw))

  # Missing-data filtering
  m <- do.call(rbind, strsplit(a2$seq, "", fixed = TRUE))
  rownames(m) <- a2$acc
  is_miss <- matrix(!(m %in% VALID), nrow = nrow(m), dimnames = dimnames(m))
  say("\nMissing data before site filtering: ", sprintf("%.1f%%", 100 * mean(is_miss)))

  keep_col <- colMeans(is_miss) <= MAX_COL_MISSING
  say("Columns kept at <=", 100 * MAX_COL_MISSING, "% missing: ",
      sum(keep_col), " / ", ncol(is_miss), sprintf(" (%.1f%%)", 100 * mean(keep_col)))

  keep_tax <- rowMeans(is_miss[, keep_col, drop = FALSE]) <= MAX_TAXON_MISSING
  if (any(!keep_tax)) {
    say("Taxa dropped at >", 100 * MAX_TAXON_MISSING, "% missing: ", sum(!keep_tax))
    keep_col <- colMeans(is_miss[keep_tax, , drop = FALSE]) <= MAX_COL_MISSING
    say("Columns kept after taxon removal: ", sum(keep_col))
  }

  final <- m[keep_tax, keep_col, drop = FALSE]
  fseq  <- apply(final, 1, paste, collapse = "")
  gf    <- g1[match(names(fseq), a1$acc)]

  d2 <- distances(code_matrix(names(fseq), unname(fseq)))
  summarise_distances(d2$p, gf, "Distances in the FINAL filtered alignment:")

  fmiss <- matrix(!(final %in% VALID), nrow = nrow(final))
  say(sprintf("\nColumn occupancy: %.1f%% of kept columns are >=90%% occupied",
              100 * mean(1 - colMeans(fmiss) >= 0.9)))
  say(sprintf("Mean pairwise identity: %.1f%%",
              100 * (1 - mean(d2$p[upper.tri(d2$p)], na.rm = TRUE))))
  say("\nFINAL (", px, "): ", length(fseq), " taxa x ", ncol(final), " sites")
  say("  missing data: ", sprintf("%.1f%%", 100 * mean(fmiss)))
  say(sprintf("  expansion vs median genome: %.2fx", ncol(final) / med_raw))

  write_fasta(names(fseq), unname(fseq), f("_filtered.fasta"))
  writeLines(c(paste(length(fseq), ncol(final)), paste(names(fseq), unname(fseq))),
             f("_filtered.phy"))
  write.csv(round(d2$p, 5), f("_distances.csv"))
  say("Wrote ", px, "_filtered.fasta / .phy / _distances.csv")

  data.frame(clade = px, taxa = length(fseq), sites = ncol(final),
             median_p = median(d2$p[upper.tri(d2$p)], na.rm = TRUE),
             max_p = max(d2$p[upper.tri(d2$p)], na.rm = TRUE),
             missing = mean(fmiss), expansion = ncol(final) / med_raw,
             stringsAsFactors = FALSE)
}

# --- Run both clades --------------------------------------------------------
say("=== kpower filovirus prep v3: ", format(Sys.time()), " ===\n")

raw <- read_fasta_full(RAW)
raw$len <- nchar(raw$seq)
say("Read ", nrow(raw), " sequences (", min(raw$len), "-", max(raw$len), " bp)")

excl <- grepl(EXCLUDE, raw$desc, ignore.case = TRUE)
qc   <- grepl(DROP_PATTERNS, raw$desc, ignore.case = TRUE)
say("  excluded ", sum(excl), " unalignable records (fish/reptile filoviruses, Lloviu)")
say("  excluded ", sum(qc & !excl), " engineered/unverified records")
pool <- raw[!excl & !qc, , drop = FALSE]

summaries <- do.call(rbind, lapply(CLADES, prep_clade, pool = pool))

say("\n\n=== Summary ===")
for (i in seq_len(nrow(summaries))) {
  s <- summaries[i, ]
  say(sprintf("  %-8s %3d taxa x %5d sites   p-dist median %.3f max %.3f   missing %.1f%%   expansion %.2fx",
              s$clade, s$taxa, s$sites, s$median_p, s$max_p,
              100 * s$missing, s$expansion))
}
say("\nGates: expansion < ~1.2x, max p-dist <= ", MAX_PDIST, ", missing well under 10%")

writeLines(log_lines, OUT_LOG)
