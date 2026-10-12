# kpower simulation tests -- Step 6: Tree-estimation accuracy
#
# QUESTION: how well does IQ-TREE recover the TRUE simulation tree from the
# simulated alignment, and does it improve with alignment length?
#
#   TRUE tree      -- the shared base tree we simulated on (seed 42, 20 taxa).
#   ESTIMATED tree -- IQ-TREE's single-tree ML fit (empirical K=1) on the
#                     alignment. One per (scenario, rep). At K=1 there is no
#                     mixture, so this single ML tree is the natural estimate
#                     regardless of which family kpower was fitting.
#
# Points are labelled by the TRUE (generating) family, taken from the scenario
# name prefix (R_/H_/T_), NOT by the fitting model.
#
# Two comparisons, each valid only where the generating classes share the
# measured property with the base tree:
#   TOPOLOGY (+R, +H): both have ONE true topology (the base). Metric = RF.
#   TREE LENGTH (+R, +T): both preserve the base tree's TOTAL length. Metric =
#     (est - true)/true total tree length. (+H excluded: its per-class branch-
#     length perturbation gives no single true length.)
#
# Data may live under results/ and/or old_R/ (extra scenarios pulled from the
# server). Both are scanned.
#
# Output: results/tree_accuracy.pdf, results/tree_length_accuracy.pdf

suppressMessages({library(ape); library(phangorn); library(ggplot2)})

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
OUT_BASE <- file.path(SCRIPT_DIR, "results")

# --- TRUE tree: regenerate exactly as 01_simulate_test_data.R does ----------
set.seed(42)
true_tree <- rtree(20, rooted = FALSE)
rtt <- dist.nodes(true_tree)[Ntip(true_tree) + 1, seq_len(Ntip(true_tree))]
true_tree$edge.length <- true_tree$edge.length * (0.45 / mean(rtt))
TL_true <- sum(true_tree$edge.length)

# --- ESTIMATED trees: the single ML tree (empirical K=1) for every dataset ---
# survey_R_linked/empirical/empirical_K1 exists for every dataset and, being a
# K=1 fit, is the plain single ML tree (family-independent).
search_dirs <- c("results", "old_R")
est_files <- unlist(lapply(search_dirs, function(d)
  Sys.glob(file.path(SCRIPT_DIR, d, "*", "rep_*", "survey_runs",
    "survey_R_linked", "empirical", "empirical_K1", "empirical_K1.treefile"))))
est_files <- unique(est_files)

if (length(est_files) == 0)
  stop("No empirical_K1 treefiles found. Pull results from server.")

parse_row <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  scn   <- parts[grep("^[RHT]_K[0-9]", parts)][1]       # scenario folder name
  rep   <- as.integer(sub(".*/rep_([0-9]+)/.*", "\\1", path))
  len   <- as.integer(sub(".*_L([0-9]+)$", "\\1", scn))
  true  <- paste0("+", substr(scn, 1, 1))               # R/H/T from prefix
  est   <- tryCatch(read.tree(path), error = function(e) NULL)
  if (is.null(est) || is.null(scn)) return(NULL)
  data.frame(
    scenario  = scn, true_family = true, seq_length = len, rep = rep,
    K         = as.integer(sub(".*_K([0-9]).*", "\\1", scn)),
    spread    = ifelse(grepl("close", scn), "close", "sep"),
    nrf       = RF.dist(true_tree, est, normalize = TRUE),
    tl_relerr = (sum(est$edge.length) - TL_true) / TL_true
  )
}

df <- do.call(rbind, lapply(est_files, parse_row))
df$len_lab <- factor(paste0("L", df$seq_length),
                     levels = paste0("L", sort(unique(df$seq_length))))
message(sprintf("Loaded %d single-ML estimates.", nrow(df)))
cat("Datasets by true family:\n"); print(table(df$true_family))

fam_levels  <- c("+R", "+H", "+T")
spread_pal  <- c("close" = "#5AAE61", "sep" = "#762A83")  # weak vs strong signal
dodge <- position_dodge(width = 0.75)

# Violins split by class-signal spread (close = weak, sep = strong); K2 and K3
# pooled within each spread (10 reps each -> 20 dots per violin).
# Panels: rows = true K (mixture classes in the DATA), cols = true family.
# Within each panel: x = length, violins split by class spread (close/sep).
# One violin = 10 reps of that (family, K, spread, length).
violin_plot <- function(d, yvar, ylab, title, sub, pctfmt = FALSE) {
  d$true_family <- droplevels(factor(d$true_family, levels = fam_levels))
  d$spread      <- factor(d$spread, levels = c("close", "sep"))
  g <- ggplot(d, aes(len_lab, .data[[yvar]], fill = spread)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.35,
                colour = NA, position = dodge) +
    geom_point(aes(colour = spread), position = position_jitterdodge(
      jitter.width = 0.15, dodge.width = 0.75),
      size = 1.2, alpha = 0.7, show.legend = FALSE) +
    facet_grid(K ~ true_family, labeller = labeller(K = function(x) paste0("K=", x))) +
    scale_fill_manual(values = spread_pal, name = "Class spread") +
    scale_colour_manual(values = spread_pal) +
    labs(x = "Alignment length (sites)", y = ylab, title = title, subtitle = sub) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
  if (pctfmt) g <- g +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    scale_y_continuous(labels = function(x) sprintf("%+.0f%%", 100 * x))
  g
}

# --- TOPOLOGY: single-true-topology families (+R, +H) -----------------------
dfTopo <- df[df$true_family %in% c("+R", "+H"), ]
pT <- violin_plot(dfTopo, "nrf",
  "Normalised RF distance to true tree",
  "IQ-TREE topology recovery vs alignment length",
  "By TRUE family x K; each dot = one replicate; 0 = perfect recovery")
ggsave(file.path(OUT_BASE, "tree_accuracy.pdf"), pT, width = 7, height = 5.5)
message("Wrote tree_accuracy.pdf")

# --- TREE LENGTH: total-length-preserving families (+R, +H, +T) -------------
# All three rescale every class to the base tree's TOTAL length (01: +H via the
# base_TL rescale in heterotachy_tree; +T via NNI; +R uses the base directly),
# so total length has one true value (= base_TL) for every family.
dfLen <- df[df$true_family %in% c("+R", "+H", "+T"), ]
pL <- violin_plot(dfLen, "tl_relerr",
  "Total tree-length error  (est - true) / true",
  "IQ-TREE branch-length recovery vs alignment length",
  "By TRUE family x K; each dot = one replicate; 0 = unbiased, +/- = over / under",
  pctfmt = TRUE)
ggsave(file.path(OUT_BASE, "tree_length_accuracy.pdf"), pL, width = 9, height = 5.5)
message("Wrote tree_length_accuracy.pdf")
