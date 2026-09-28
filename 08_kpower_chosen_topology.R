# kpower simulation tests -- Step 8: does kpower's CHOSEN model recover the
# true TOPOLOGY? (companion to 07, which did total length.)
#
# Estimate = tree(s) from kpower's selected family+K (best_mix_type, best_K).
# Metric   = normalised RF to the base tree.
# Scope    = +R and +H TRUTHS only: these have ONE true topology (the base), so
#            base is the valid reference. +T-truth has K real topologies and
#            needs each rep's class_*_tree.nwk (not pulled) -- excluded here.
# When kpower's PICK is a tree mixture (+T), its treefile has several topologies;
# we summarise with the MEAN RF over components.
#
# Output: results/tree_topology_kpower_chosen.pdf

suppressMessages({library(ape); library(phangorn); library(ggplot2)})

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd())
OUT_BASE <- file.path(SCRIPT_DIR, "results")

set.seed(42)
tt <- rtree(20, rooted = FALSE)
rtt <- dist.nodes(tt)[Ntip(tt) + 1, seq_len(Ntip(tt))]
tt$edge.length <- tt$edge.length * (0.45 / mean(rtt))

meta <- read.csv(file.path(OUT_BASE, "summary.csv.bak"), stringsAsFactors = FALSE)

chosen_treefile <- function(base_dir, scn, rep, fam, K) {
  rd  <- file.path(base_dir, scn, sprintf("rep_%02d", rep), "survey_runs")
  lab <- sprintf("empirical_K%d", K)
  switch(fam,
    "+R" = file.path(rd, "survey_R_linked", "empirical", lab, paste0(lab, ".treefile")),
    "+H" = file.path(rd, "survey_H_linked", "empirical", lab, paste0(lab, ".treefile")),
    "+T" = file.path(rd, "survey_T_linked", "mast_empirical", lab, paste0(lab, ".treefile")),
    NA_character_)
}

# mean normalised RF to the base tree over the tree(s) in a treefile
mean_nrf <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_real_)
  tr <- tryCatch(read.tree(path), error = function(e) NULL)
  if (is.null(tr)) return(NA_real_)
  if (inherits(tr, "multiPhylo"))
    mean(vapply(tr, function(x) RF.dist(tt, x, normalize = TRUE), numeric(1)))
  else RF.dist(tt, tr, normalize = TRUE)
}

meta <- meta[meta$true_type %in% c("+R", "+H"), ]        # single-true-topology only
meta$base_dir <- ifelse(startsWith(meta$scenario, "R_"),
                        file.path(SCRIPT_DIR, "old_R"), OUT_BASE)
meta$path <- mapply(chosen_treefile, meta$base_dir, meta$scenario, meta$rep,
                    meta$best_mix_type, meta$best_K)
meta$nrf  <- vapply(meta$path, mean_nrf, numeric(1))

n_miss <- sum(is.na(meta$nrf))
if (n_miss > 0) message(sprintf("%d/%d datasets missing chosen treefile (skipped).",
                                n_miss, nrow(meta)))
d <- meta[!is.na(meta$nrf), ]
d$true_family <- factor(d$true_type, levels = c("+R", "+H"))
d$spread  <- factor(d$tag, levels = c("close", "sep"))
d$len_lab <- factor(paste0("L", d$seq_length),
                    levels = paste0("L", sort(unique(d$seq_length))))

cat("Mean normalised RF of kpower's CHOSEN model to base, by true family:\n")
print(aggregate(nrf ~ true_family, d, mean), row.names = FALSE, digits = 3)

spread_pal <- c("close" = "#5AAE61", "sep" = "#762A83")
p <- ggplot(d, aes(len_lab, nrf, fill = spread)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.35, colour = NA,
              position = position_dodge(width = 0.75)) +
  geom_point(aes(colour = spread), position = position_jitterdodge(
    jitter.width = 0.15, dodge.width = 0.75), size = 1.2, alpha = 0.7,
    show.legend = FALSE) +
  facet_grid(true_K ~ true_family,
             labeller = labeller(true_K = function(x) paste0("K=", x))) +
  scale_fill_manual(values = spread_pal, name = "Class spread") +
  scale_colour_manual(values = spread_pal) +
  labs(x = "Alignment length (sites)",
       y = "Normalised RF distance to true tree",
       title = "Topology recovery of kpower's CHOSEN model",
       subtitle = "Estimate = tree(s) from best_mix_type + best_K; vs base tree; 0 = perfect") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_BASE, "tree_topology_kpower_chosen.pdf"), p, width = 7, height = 5.5)
message("Wrote tree_topology_kpower_chosen.pdf")
