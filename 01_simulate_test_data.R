# kpower simulation tests -- Step 1: Simulate all test alignments
#
# Design (48 scenarios = 3 x 2 x 2 x 4):
#   - 3 families:       +R, +H, +T
#   - 2 class types:    "sep" (well-separated), "close" (subtle)
#   - 2 K values:       2, 4
#   - 4 seq lengths:    100, 300, 1000, 3000 (log-spaced)
#   - 20 taxa, shared random tree (mean RTT = 0.45 subs/site)
#
# A single "spread" parameter per family per class type is applied
# identically to K=2 and K=4:
#   +R: rates are a geometric series from 1/s to s, equal proportions,
#       rescaled so weighted mean rate = 1 (preserves overall tree length)
#       sep   -> s = 5     (25x fast/slow ratio)
#       close -> s = 1.86  (~3.5x ratio)
#   +H: each class has independently-perturbed branch lengths via lognormal
#       deviates with sd = log_sd
#       sep   -> log_sd = 1.5
#       close -> log_sd = 0.3
#   +T: K class trees, each separated from a base tree by n_nni NNI moves
#       sep   -> n_nni = 10  (RF~16, ~47% of max for 20 taxa)
#       close -> n_nni = 1   (RF~2)
#
# Output: kpower_tests/alignments/<scenario>/sim.phy + sim_params.txt

library(ape)
library(phangorn)

# --- Configuration ----------------------------------------------------------
set.seed(42)

IQTREE  <- Sys.which("iqtree2")
#if (!nzchar(IQTREE))
#  IQTREE <- path.expand("~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3")

N_TAXA          <- 20
LENGTHS         <- c(100, 300, 1000, 3000)
TARGET_MEAN_RTT <- 0.45            # mean root-to-tip subs/site (Duchene 2017)
KS              <- c(2, 4)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE <- file.path(SCRIPT_DIR, "alignments")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Spread parameters per class type (identical across K within a family)
SPREAD <- list(
  R = list(sep = 5.0,  close = 1.86),   # geometric spread for rates
  H = list(sep = 1.5,  close = 0.3),    # log_sd for branch perturbation
  T = list(sep = 10L,  close = 1L)      # NNI moves between class trees
)

# Shared GTR -- strong ti/tv bias, realistic for nuclear DNA
GTR_RATES <- "5.0,1.0,1.0,1.0,5.0,1.0"   # AC,AG,AT,CG,CT,GT
BASE_FREQ <- "0.30,0.20,0.20,0.30"
MODEL_STR <- paste0("GTR{", GTR_RATES, "}+FU{", BASE_FREQ, "}")


# --- Shared tree ------------------------------------------------------------
tree <- rtree(N_TAXA, rooted = FALSE)
rtt  <- dist.nodes(tree)[Ntip(tree) + 1, seq_len(Ntip(tree))]
tree$edge.length <- tree$edge.length * (TARGET_MEAN_RTT / mean(rtt))
TREE_FILE <- file.path(OUT_BASE, "test_tree.nwk")
write.tree(tree, TREE_FILE)
message(sprintf("Tree: mean RTT = %.3f, total length = %.3f",
                mean(dist.nodes(tree)[Ntip(tree) + 1, seq_len(Ntip(tree))]),
                sum(tree$edge.length)))


# =============================================================================
# Helpers
# =============================================================================

#' Run AliSim for a single alignment
alisim_one <- function(model_str, tree_file, n_sites, prefix, seed) {
  args <- c(
    "--alisim", prefix,
    "-m",       model_str,
    "-t",       tree_file,
    "--length", as.character(n_sites),
    "--seed",   as.character(seed),
    "--redo"
  )
  result <- processx::run(IQTREE, args, error_on_status = FALSE)
  if (result$status != 0)
    stop("AliSim failed:\n", result$stderr)
  for (ext in c(".phy", ".fa", ".fasta")) {
    f <- paste0(prefix, ext)
    if (file.exists(f)) return(f)
  }
  stop("No AliSim output found at: ", prefix)
}

#' Read a sequential PHYLIP alignment as named character vector
read_phy <- function(path) {
  lines    <- readLines(path)
  nonempty <- lines[nzchar(trimws(lines))]
  ntaxa    <- as.integer(strsplit(trimws(nonempty[1]), "\\s+")[[1]][1])
  taxa     <- character(ntaxa)
  seqs     <- character(ntaxa)
  for (i in seq_len(ntaxa)) {
    parts   <- strsplit(trimws(nonempty[i + 1]), "\\s+")[[1]]
    taxa[i] <- parts[1]
    seqs[i] <- paste(parts[-1], collapse = "")
  }
  names(seqs) <- taxa
  seqs
}

#' Concatenate multiple PHYLIP alignments site-wise
concat_phy <- function(phy_files, outfile) {
  alns     <- lapply(phy_files, read_phy)
  taxa     <- names(alns[[1]])
  combined <- rep("", length(taxa))
  names(combined) <- taxa
  for (aln in alns) combined <- paste0(combined, aln[taxa])
  out <- c(paste(length(taxa), nchar(combined[1])))
  for (i in seq_along(taxa)) out <- c(out, paste(taxa[i], combined[i]))
  writeLines(out, outfile)
  outfile
}

#' Heterotachous tree: independently perturb each branch by a lognormal deviate
heterotachy_tree <- function(tr, log_sd = 1.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  noise <- exp(rnorm(length(tr$edge.length), mean = 0, sd = log_sd))
  tr$edge.length <- tr$edge.length * noise
  tr
}

#' Topologically distinct tree via random NNI moves (proper bipartition change)
nni_tree <- function(tr, n_moves = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rNNI(tr, moves = n_moves)
}

#' Build +R{K} model string with geometric rates + equal proportions, mean=1
#'
#' Rates form a geometric series from 1/s to s (K points), then rescaled so
#' that the weighted mean rate equals 1 (equal proportions = 1/K).
build_freerate_model <- function(base_model, K, s) {
  if (K == 1) {
    raw <- 1
  } else {
    raw <- s ^ seq(-1, 1, length.out = K)
  }
  rates <- raw / mean(raw)                 # equal weights -> mean = 1
  props <- rep(1 / K, K)
  interleaved <- as.vector(rbind(props, rates))
  paste0(base_model, "+R", K, "{",
         paste(format(interleaved, trim = TRUE), collapse = ","), "}")
}

#' Simulate a multi-class mixture alignment by per-class AliSim + concat
simulate_mixture <- function(class_trees, weights, model_str, len,
                             sim_dir, seed_base) {
  K <- length(weights)
  n_per_class <- round(weights * len)
  residual <- len - sum(n_per_class)
  if (residual != 0)
    n_per_class[which.max(n_per_class)] <-
      n_per_class[which.max(n_per_class)] + residual

  phy_files <- character(K)
  for (k in seq_len(K)) {
    tree_file <- file.path(sim_dir, paste0("class_", k, "_tree.nwk"))
    write.tree(class_trees[[k]], tree_file)
    phy_files[k] <- alisim_one(
      model_str = model_str,
      tree_file = tree_file,
      n_sites   = n_per_class[k],
      prefix    = file.path(sim_dir, paste0("class_", k)),
      seed      = seed_base + k
    )
  }
  concat_phy(phy_files, file.path(sim_dir, "sim.phy"))
}

write_params <- function(sim_dir, fields) {
  writeLines(
    vapply(names(fields), function(k) paste0(k, ": ", fields[[k]]),
           character(1)),
    file.path(sim_dir, "sim_params.txt")
  )
}


# =============================================================================
# +R (FreeRate) -- 20 scenarios
# =============================================================================
message("\n========== FreeRate (+R) scenarios ==========")
seed_counter <- 1000
for (K in KS) for (tag in c("sep", "close")) {
  s          <- SPREAD$R[[tag]]
  model_full <- build_freerate_model(MODEL_STR, K, s)
  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("R_K%d_%s_L%d", K, tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
    message(sprintf("[%s]  %s", label, model_full))
    alisim_one(model_full, TREE_FILE, len,
               file.path(sim_dir, "sim"), seed_counter)
    write_params(sim_dir, list(
      type = "+R", K = K, tag = tag, length = len, seed = seed_counter,
      spread_s = s, model = model_full, tree = TREE_FILE
    ))
  }
}


# =============================================================================
# +H (GHOST heterotachy) -- 20 scenarios
# =============================================================================
message("\n========== GHOST (+H) scenarios ==========")
seed_counter <- 2000
for (K in KS) for (tag in c("sep", "close")) {
  log_sd <- SPREAD$H[[tag]]
  # K independent class trees -- same topology, perturbed branch lengths
  tree_seeds  <- 5000 + 100 * K + which(c("sep", "close") == tag) * 10 +
                 seq_len(K)
  class_trees <- lapply(tree_seeds, function(s)
    heterotachy_tree(tree, log_sd = log_sd, seed = s))
  weights <- rep(1 / K, K)

  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("H_K%d_%s_L%d", K, tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
    message(sprintf("[%s]  log_sd=%.2f  K=%d", label, log_sd, K))
    simulate_mixture(class_trees, weights, MODEL_STR, len, sim_dir,
                     seed_counter)
    write_params(sim_dir, list(
      type = "+H", K = K, tag = tag, length = len, seed = seed_counter,
      log_sd = log_sd, weights = paste(weights, collapse = ","),
      tree_seeds = paste(tree_seeds, collapse = ","),
      model = MODEL_STR, tree = TREE_FILE
    ))
  }
}


# =============================================================================
# +T (MAST tree-mixture) -- 20 scenarios
# =============================================================================
message("\n========== MAST (+T) scenarios ==========")
seed_counter <- 3000
for (K in KS) for (tag in c("sep", "close")) {
  n_nni <- SPREAD$T[[tag]]
  # K class trees: tree_1 = original; tree_k (k>=2) = original + n_nni NNI
  # moves with distinct seeds, so each pair is roughly n_nni-divergent.
  class_trees <- vector("list", K)
  class_trees[[1]] <- tree
  tree_seeds <- 6000 + 100 * K + which(c("sep", "close") == tag) * 10 +
                seq_len(K - 1)
  for (k in seq_len(K - 1))
    class_trees[[k + 1]] <- nni_tree(tree, n_moves = n_nni,
                                     seed = tree_seeds[k])
  weights <- rep(1 / K, K)

  rfs <- vapply(class_trees[-1], function(t) RF.dist(tree, t), integer(1))
  message(sprintf("[+T K=%d %s]  n_nni=%d  RF(orig, class_k)=%s  (max=%d)",
                  K, tag, n_nni, paste(rfs, collapse = ","),
                  2 * (Ntip(tree) - 3)))

  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("T_K%d_%s_L%d", K, tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)
    simulate_mixture(class_trees, weights, MODEL_STR, len, sim_dir,
                     seed_counter)
    write_params(sim_dir, list(
      type = "+T", K = K, tag = tag, length = len, seed = seed_counter,
      n_nni = n_nni, weights = paste(weights, collapse = ","),
      RF_to_class1 = paste(c(0, rfs), collapse = ","),
      tree_seeds = paste(c(NA, tree_seeds), collapse = ","),
      model = MODEL_STR, tree = TREE_FILE
    ))
  }
}

message(sprintf("\nDone. %d scenarios written to: %s",
                3 * length(KS) * 2 * length(LENGTHS), OUT_BASE))
