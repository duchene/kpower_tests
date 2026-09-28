# kpower simulation tests -- Step 7: does kpower's CHOSEN model recover the
# true total tree length?
#
# For each dataset kpower picks a family (best_mix_type) and K (best_K) by BIC.
# We read THAT chosen model's estimated tree(s) and compare total tree length to
# the base tree (every true class preserves base_TL, so it is the shared truth).
# Mixture picks (+H, +T) give several trees; we take the mean component length.
#
# Estimate = kpower-chosen tree only (no K=1 baseline). Metric = (est-true)/true.
# Output: results/tree_length_kpower_chosen.pdf

suppressMessages({library(ape); library(phangorn); library(ggplot2)})

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd())
OUT_BASE <- file.path(SCRIPT_DIR, "results")

# TRUE tree total length (seed 42, as in 01)
set.seed(42)
tt <- rtree(20, rooted = FALSE)
rtt <- dist.nodes(tt)[Ntip(tt) + 1, seq_len(Ntip(tt))]
tt$edge.length <- tt$edge.length * (0.45 / mean(rtt))
TL_true <- sum(tt$edge.length)

meta <- read.csv(file.path(OUT_BASE, "summary.csv.bak"), stringsAsFactors = FALSE)

# survey subdir + within-family layout for each family's empirical K fit
chosen_treefile <- function(base_dir, scn, rep, fam, K) {
  rd <- file.path(base_dir, scn, sprintf("rep_%02d", rep), "survey_runs")
  lab <- sprintf("empirical_K%d", K)
  switch(fam,
    "+R" = file.path(rd, "survey_R_linked", "empirical", lab, paste0(lab, ".treefile")),
    "+H" = file.path(rd, "survey_H_linked", "empirical", lab, paste0(lab, ".treefile")),
    "+T" = file.path(rd, "survey_T_linked", "mast_empirical", lab, paste0(lab, ".treefile")),
    NA_character_)
}

# Total tree length of the chosen model's estimate.
#  +R : single tree -> its length.
#  +H : GHOST writes [summary tree, then K per-class trees]; the FIRST tree is
#       the class-weight-weighted-average tree (== IQ-TREE's reported "Total tree
#       length"), so use it. Per-class trees can explode at ~0 weight, so a naive
#       mean is wrong.
#  +T : MAST writes K component trees, no summary; each true class preserves
#       base_TL, so mean of component lengths is the summary.
summary_TL <- function(path, fam) {
  if (is.na(path) || !file.exists(path)) return(NA_real_)
  tr <- tryCatch(read.tree(path), error = function(e) NULL)
  if (is.null(tr)) return(NA_real_)
  if (!inherits(tr, "multiPhylo")) return(sum(tr$edge.length))
  tl <- vapply(tr, function(x) sum(x$edge.length), numeric(1))
  if (fam == "+T") mean(tl) else tl[1]      # +T: mean components; +H: summary (first)
}

meta$base_dir <- ifelse(startsWith(meta$scenario, "R_"),
                        file.path(SCRIPT_DIR, "old_R"), OUT_BASE)
meta$path <- mapply(chosen_treefile, meta$base_dir, meta$scenario, meta$rep,
                    meta$best_mix_type, meta$best_K)
meta$tl   <- mapply(summary_TL, meta$path, meta$best_mix_type)
meta$tl_relerr <- (meta$tl - TL_true) / TL_true

n_miss <- sum(is.na(meta$tl))
if (n_miss > 0) message(sprintf("%d/%d datasets missing chosen treefile (skipped).",
                                n_miss, nrow(meta)))
d <- meta[!is.na(meta$tl_relerr), ]
d$true_family <- factor(d$true_type, levels = c("+R", "+H", "+T"))
d$spread   <- factor(d$tag, levels = c("close", "sep"))
d$len_lab  <- factor(paste0("L", d$seq_length),
                     levels = paste0("L", sort(unique(d$seq_length))))

cat("Mean total-length error of kpower's CHOSEN model, by true family:\n")
print(aggregate(tl_relerr ~ true_family, d, mean), row.names = FALSE, digits = 3)

spread_pal <- c("close" = "#5AAE61", "sep" = "#762A83")
p <- ggplot(d, aes(len_lab, tl_relerr, fill = spread)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.35, colour = NA,
              position = position_dodge(width = 0.75)) +
  geom_point(aes(colour = spread), position = position_jitterdodge(
    jitter.width = 0.15, dodge.width = 0.75), size = 1.2, alpha = 0.7,
    show.legend = FALSE) +
  facet_grid(true_K ~ true_family,
             labeller = labeller(true_K = function(x) paste0("K=", x))) +
  scale_fill_manual(values = spread_pal, name = "Class spread") +
  scale_colour_manual(values = spread_pal) +
  scale_y_continuous(labels = function(x) sprintf("%+.0f%%", 100 * x)) +
  labs(x = "Alignment length (sites)",
       y = "Total tree-length error  (est - true) / true",
       title = "Total-length recovery of kpower's CHOSEN model",
       subtitle = "Estimate = tree(s) from best_mix_type + best_K; vs base tree; 0 = unbiased") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_BASE, "tree_length_kpower_chosen.pdf"), p, width = 9, height = 5.5)
message("Wrote tree_length_kpower_chosen.pdf")
