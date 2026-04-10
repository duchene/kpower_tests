# kpower simulation tests — Step 1: Simulate all test alignments
#
# Generates test alignments under known FreeRate (+R), GHOST (+H), and
# MAST (+T) models with explicit parameters.
#
# Design:
#   - "sep" : well-separated classes (power should be high)
#   - "close": closely-spaced classes (power should be low)
#   - Two sequence lengths (100, 1000) per scenario
#
# +R: rate categories with known proportions and rates
# +H: same topology, per-class branch-length scaling (heterotachy)
# +T: different topologies per class (tree mixture)
#
# Output: kpower_tests/alignments/<scenario>/sim.phy + sim_params.txt

library(ape)
library(phangorn)

# --- Configuration -----------------------------------------------------------
set.seed(42)
IQTREE  <- Sys.which("iqtree3")
if (!nzchar(IQTREE))
  IQTREE <- path.expand("~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3")

N_TAXA   <- 20
LENGTHS  <- c(100, 1000)
TARGET_MEAN_RTT <- 0.45   # mean root-to-tip subs/site (Duchene et al. 2017)
OUT_BASE <- path.expand("~/Dropbox/Research/hiv_mast/kpower_tests/alignments")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Shared GTR parameters — strong ti/tv bias, realistic for nuclear DNA
GTR_RATES <- "5.0,1.0,1.0,1.0,5.0,1.0"   # AC,AG,AT,CG,CT,GT
BASE_FREQ <- "0.30,0.20,0.20,0.30"
MODEL_STR <- paste0("GTR{", GTR_RATES, "}+FU{", BASE_FREQ, "}")

# --- Generate a random tree --------------------------------------------------
# Scale so that the mean root-to-tip distance equals TARGET_MEAN_RTT
# (following Duchene et al. 2017, Syst Biol 66:769–785).
tree <- rtree(N_TAXA, rooted = FALSE)
rtt  <- dist.nodes(tree)[Ntip(tree) + 1, seq_len(Ntip(tree))]
tree$edge.length <- tree$edge.length * (TARGET_MEAN_RTT / mean(rtt))
message(sprintf("Tree scaled: mean RTT = %.3f, total length = %.3f",
                mean(dist.nodes(tree)[Ntip(tree) + 1, seq_len(Ntip(tree))]),
                sum(tree$edge.length)))
TREE_FILE <- file.path(OUT_BASE, "test_tree.nwk")
write.tree(tree, TREE_FILE)
message("Tree written to: ", TREE_FILE)

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
  # Return path to output alignment
  for (ext in c(".phy", ".fa", ".fasta")) {
    f <- paste0(prefix, ext)
    if (file.exists(f)) return(f)
  }
  stop("No AliSim output found at: ", prefix)
}

#' Read a sequential PHYLIP alignment as named character vector
read_phy <- function(path) {
  lines <- readLines(path)
  nonempty <- lines[nzchar(trimws(lines))]
  ntaxa <- as.integer(strsplit(trimws(nonempty[1]), "\\s+")[[1]][1])
  taxa <- character(ntaxa)
  seqs <- character(ntaxa)
  for (i in seq_len(ntaxa)) {
    parts <- strsplit(trimws(nonempty[i + 1]), "\\s+")[[1]]
    taxa[i] <- parts[1]
    seqs[i] <- paste(parts[-1], collapse = "")
  }
  names(seqs) <- taxa
  seqs
}

#' Concatenate multiple PHYLIP alignments site-wise
concat_phy <- function(phy_files, outfile) {
  alns <- lapply(phy_files, read_phy)
  taxa <- names(alns[[1]])
  combined <- rep("", length(taxa))
  names(combined) <- taxa
  for (aln in alns) combined <- paste0(combined, aln[taxa])
  out <- c(paste(length(taxa), nchar(combined[1])))
  for (i in seq_along(taxa)) out <- c(out, paste(taxa[i], combined[i]))
  writeLines(out, outfile)
  outfile
}

#' Scale all branch lengths by a factor
scale_tree <- function(tr, factor) {
  tr$edge.length <- tr$edge.length * factor
  tr
}

#' Generate a heterotachous tree by independently perturbing branch lengths
#'
#' Multiplies each branch length by an independent lognormal deviate with
#' the given sd on the log scale, then scales the total tree length to
#' `target_length`. This produces a tree with the same topology but
#' different relative branch-length proportions — true heterotachy.
#'
#' @param tr An ape phylo object.
#' @param log_sd SD of log-normal noise per branch. Larger = more heterotachy.
#'   0.5 = moderate, 1.5 = extreme.
#' @param target_length Desired total tree length (sum of branch lengths).
#' @param seed Random seed for reproducibility.
heterotachy_tree <- function(tr, log_sd = 1.0, target_length = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_edges <- length(tr$edge.length)
  # Multiply each branch by an independent lognormal deviate
  noise <- exp(rnorm(n_edges, mean = 0, sd = log_sd))
  tr$edge.length <- tr$edge.length * noise
  # Rescale to target total length if requested
  if (!is.null(target_length)) {
    tr$edge.length <- tr$edge.length / sum(tr$edge.length) * target_length
  }
  tr
}

#' Create a topologically distinct tree via random NNI moves
#'
#' Uses phangorn::rNNI which performs proper nearest-neighbour interchange
#' on the unrooted topology (changes bipartitions, nonzero RF distance).
nni_tree <- function(tr, n_moves = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rNNI(tr, moves = n_moves)
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
    message(sprintf("  Class %d: %d sites (weight %.3f)",
                    k, n_per_class[k], weights[k]))
  }
  concat_phy(phy_files, file.path(sim_dir, "sim.phy"))
}


# =============================================================================
# +R (FreeRate) scenarios
# =============================================================================
message("\n========== FreeRate (+R) scenarios ==========")

freerate_scenarios <- list(
  list(K = 2, tag = "sep",
       rates = "0.1,2.5", props = "0.5,0.5"),
  list(K = 2, tag = "close",
       rates = "0.7,1.3", props = "0.5,0.5"),
  list(K = 4, tag = "sep",
       rates = "0.01,0.2,1.0,5.0", props = "0.25,0.25,0.25,0.25"),
  list(K = 4, tag = "close",
       rates = "0.5,0.8,1.2,1.6", props = "0.25,0.25,0.25,0.25")
)

seed_counter <- 100
for (sc in freerate_scenarios) {
  K <- sc$K
  props_vec <- as.numeric(strsplit(sc$props, ",")[[1]])
  rates_vec <- as.numeric(strsplit(sc$rates, ",")[[1]])
  interleaved <- as.vector(rbind(props_vec, rates_vec))
  rk_str <- paste0("+R", K, "{", paste(interleaved, collapse = ","), "}")
  model_full <- paste0(MODEL_STR, rk_str)

  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("R_K%d_%s_L%d", K, sc$tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

    message(sprintf("Simulating %s: %s", label, model_full))
    alisim_one(model_full, TREE_FILE, len,
               file.path(sim_dir, "sim"), seed_counter)

    writeLines(
      c(paste("model:", model_full),
        paste("type:", "+R"),
        paste("tree:", TREE_FILE),
        paste("length:", len),
        paste("seed:", seed_counter),
        paste("K:", K),
        paste("rates:", sc$rates),
        paste("proportions:", sc$props),
        paste("tag:", sc$tag)),
      file.path(sim_dir, "sim_params.txt")
    )
    message("  -> OK")
  }
}


# =============================================================================
# +H (GHOST heterotachy) scenarios
# =============================================================================
message("\n========== GHOST (+H) scenarios ==========")

# Each class gets independently perturbed branch lengths (same topology).
# "sep"  : high log_sd = very different relative branch proportions across classes
# "close": low log_sd  = similar relative proportions (subtle heterotachy)
# Classes also naturally differ in overall tree length — no rescaling.

ghost_scenarios <- list(
  list(K = 2, tag = "sep",   log_sd = 1.5,
       weights = c(0.5, 0.5), seeds = c(51, 52)),
  list(K = 2, tag = "close", log_sd = 0.3,
       weights = c(0.5, 0.5), seeds = c(61, 62)),
  list(K = 3, tag = "sep",   log_sd = 1.5,
       weights = c(0.3, 0.4, 0.3), seeds = c(71, 72, 73))
)

seed_counter <- 200
for (sc in ghost_scenarios) {
  class_trees <- lapply(sc$seeds, function(s) {
    heterotachy_tree(tree, log_sd = sc$log_sd, seed = s)
  })

  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("H_K%d_%s_L%d", sc$K, sc$tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

    message(sprintf("Simulating %s: K=%d, log_sd=%.1f, weights=%s",
                    label, sc$K, sc$log_sd,
                    paste(sc$weights, collapse = ",")))

    simulate_mixture(class_trees, sc$weights, MODEL_STR, len,
                     sim_dir, seed_counter)

    writeLines(
      c(paste("model:", MODEL_STR),
        paste("type:", "+H"),
        paste("tree:", TREE_FILE),
        paste("length:", len),
        paste("seed:", seed_counter),
        paste("K:", sc$K),
        paste("log_sd:", sc$log_sd),
        paste("weights:", paste(sc$weights, collapse = ",")),
        paste("tag:", sc$tag)),
      file.path(sim_dir, "sim_params.txt")
    )
    message("  -> OK")
  }
}


# =============================================================================
# +T (MAST tree-mixture) scenarios
# =============================================================================
message("\n========== MAST (+T) scenarios ==========")

tree_close <- nni_tree(tree, n_moves = 1, seed = 10)
tree_sep   <- nni_tree(tree, n_moves = 10, seed = 20)
message(sprintf("  RF(orig, close): %d   RF(orig, sep): %d   (max %d)",
                RF.dist(tree, tree_close), RF.dist(tree, tree_sep),
                2 * (Ntip(tree) - 3)))

mast_scenarios <- list(
  list(K = 2, tag = "sep",
       trees = list(tree, tree_sep), weights = c(0.5, 0.5)),
  list(K = 2, tag = "close",
       trees = list(tree, tree_close), weights = c(0.5, 0.5))
)

seed_counter <- 300
for (sc in mast_scenarios) {
  for (len in LENGTHS) {
    seed_counter <- seed_counter + 1
    label   <- sprintf("T_K%d_%s_L%d", sc$K, sc$tag, len)
    sim_dir <- file.path(OUT_BASE, label)
    dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

    message(sprintf("Simulating %s: K=%d, weights=%s",
                    label, sc$K, paste(sc$weights, collapse = ",")))

    simulate_mixture(sc$trees, sc$weights, MODEL_STR, len,
                     sim_dir, seed_counter)

    writeLines(
      c(paste("model:", MODEL_STR),
        paste("type:", "+T"),
        paste("tree:", TREE_FILE),
        paste("length:", len),
        paste("seed:", seed_counter),
        paste("K:", sc$K),
        paste("weights:", paste(sc$weights, collapse = ",")),
        paste("tag:", sc$tag)),
      file.path(sim_dir, "sim_params.txt")
    )
    message("  -> OK")
  }
}

message("\nAll test alignments written to: ", OUT_BASE)
message(sprintf("Total scenarios: %d +R, %d +H, %d +T",
                length(freerate_scenarios) * length(LENGTHS),
                length(ghost_scenarios) * length(LENGTHS),
                length(mast_scenarios) * length(LENGTHS)))
