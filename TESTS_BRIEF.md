# kpower validation tests — Brief

## Purpose

Validate that `kpower_survey()` can correctly identify which mixture model
family (+R, +H, +T) generated an alignment and recover the true number of
mixture categories K, across controlled simulation scenarios with known
ground truth.

---

## Test design

### Overview

- 18 scenarios: 8 FreeRate (+R), 6 GHOST (+H), 4 MAST (+T)
- Each scenario produces one simulated alignment (`sim.phy`)
- Each alignment is surveyed with `kpower_survey()` across +R, +H, +T
  (K = 1..4, B = 20 bootstrap replicates)
- Primary IC: BIC
- 20 taxa, two alignment lengths (100 and 1000 sites) per condition
- "sep" scenarios: well-separated classes (high signal)
- "close" scenarios: subtle class differences (low signal)

### Scripts

- `01_simulate_test_data.R` — generates all test alignments via IQ-TREE
  AliSim. Depends on `ape` and `phangorn` for tree manipulation, `processx`
  for IQ-TREE.
- `02_run_kpower_tests.R` — runs `kpower_survey()` on each alignment and
  writes `results/summary.csv` plus per-scenario RDS and PDF files.

### Shared parameters

- GTR exchangeability: A-C=5, A-G=1, A-T=1, C-G=1, C-T=5, G-T=1
  (strong transition bias, realistic for nuclear DNA)
- Base frequencies: A=0.30, C=0.20, G=0.20, T=0.30
- Tree: 20-taxon random tree (`ape::rtree`), scaled so that the mean
  root-to-tip distance = 0.45 expected substitutions/site (following
  Duchêné et al. 2017, Syst Biol 66:769–785). This yields a total
  tree length of ~3.9.

### Test runner settings

```r
K_MAX      <- 4       # evaluate K = 1..4
B          <- 20      # bootstrap replicates (low, for speed)
N_CORES    <- 4       # parallel R workers
THREADS    <- 1       # IQ-TREE threads per run
FIXED_TREE <- NULL    # --fast heuristic (IQ-TREE 3.0.1 ARM BioNJ bug)
FAST_TREES <- TRUE    # GTR+R --fast for MAST candidate trees
MIX_TYPES  <- c("+R", "+H", "+T")
```

---

## Scenario details

### FreeRate (+R) — 8 scenarios

Simulated directly with AliSim using explicit `+R{K}{prop,rate,...}` syntax.

| Scenario | K | Tag | Sites | Rates | Proportions |
|---|---|---|---|---|---|
| R_K2_sep_L100 | 2 | sep | 100 | 0.1, 2.5 | 0.5, 0.5 |
| R_K2_sep_L1000 | 2 | sep | 1000 | 0.1, 2.5 | 0.5, 0.5 |
| R_K2_close_L100 | 2 | close | 100 | 0.7, 1.3 | 0.5, 0.5 |
| R_K2_close_L1000 | 2 | close | 1000 | 0.7, 1.3 | 0.5, 0.5 |
| R_K4_sep_L100 | 4 | sep | 100 | 0.01, 0.2, 1.0, 5.0 | equal |
| R_K4_sep_L1000 | 4 | sep | 1000 | 0.01, 0.2, 1.0, 5.0 | equal |
| R_K4_close_L100 | 4 | close | 100 | 0.5, 0.8, 1.2, 1.6 | equal |
| R_K4_close_L1000 | 4 | close | 1000 | 0.5, 0.8, 1.2, 1.6 | equal |

### GHOST (+H) — 6 scenarios

Simulated by per-class AliSim + concatenation. Each class uses the same
tree topology but with **independently perturbed branch lengths** — each
branch is multiplied by an independent lognormal deviate
(`exp(N(0, log_sd))`). This creates true heterotachy where relative
branch-length proportions differ across classes, not just overall rate.

**Important**: an earlier version incorrectly used uniform tree scaling
(all branches multiplied by the same factor), which is indistinguishable
from +R rate variation. This was fixed by using `heterotachy_tree()` with
independent per-branch perturbations.

Classes are NOT rescaled to the same total tree length — they naturally
differ in both relative branch proportions and overall rate.

| Scenario | K | Tag | Sites | log_sd | Weights | Seeds |
|---|---|---|---|---|---|---|
| H_K2_sep_L100 | 2 | sep | 100 | 1.5 | 0.5, 0.5 | 51, 52 |
| H_K2_sep_L1000 | 2 | sep | 1000 | 1.5 | 0.5, 0.5 | 51, 52 |
| H_K2_close_L100 | 2 | close | 100 | 0.3 | 0.5, 0.5 | 61, 62 |
| H_K2_close_L1000 | 2 | close | 1000 | 0.3 | 0.5, 0.5 | 61, 62 |
| H_K3_sep_L100 | 3 | sep | 100 | 1.5 | 0.3, 0.4, 0.3 | 71, 72, 73 |
| H_K3_sep_L1000 | 3 | sep | 1000 | 1.5 | 0.3, 0.4, 0.3 | 71, 72, 73 |

### MAST (+T) — 4 scenarios

Simulated by per-class AliSim + concatenation with distinct tree
topologies. Topological variation created via `phangorn::rNNI()` which
performs proper nearest-neighbour interchange on the unrooted topology,
producing trees with genuinely different bipartition sets (nonzero
Robinson-Foulds distance). Max possible RF for 20 taxa = 34.

- "sep": 10 NNI moves (RF = 16, 47% of max — topologically distant)
- "close": 1 NNI move (RF = 2, minimal perturbation)

| Scenario | K | Tag | Sites | NNI moves | RF dist | Weights |
|---|---|---|---|---|---|---|
| T_K2_sep_L100 | 2 | sep | 100 | 10 | 16 | 0.5, 0.5 |
| T_K2_sep_L1000 | 2 | sep | 1000 | 10 | 16 | 0.5, 0.5 |
| T_K2_close_L100 | 2 | close | 100 | 1 | 2 | 0.5, 0.5 |
| T_K2_close_L1000 | 2 | close | 1000 | 1 | 2 | 0.5, 0.5 |

---

## Results (run 2026-04-09)

Settings: mean RTT = 0.45 (TL ~3.9), `phangorn::rNNI()` for +T trees
(RF=2 close, RF=16 sep), alignment lengths 100 and 1000, B=20, K_MAX=4.
Post-bugfix run: MAST bootstrap concatenation and error handling fixed.

### Full table

| # | Scenario | True | BIC Best | +R K (pwr) | +H K (pwr) | +T K (pwr) | Family | K |
|---|---|---|---|---|---|---|---|---|
| 1 | H_K2_close_L100 | +H K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 2 | H_K2_close_L1000 | +H K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 3 | H_K2_sep_L100 | +H K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 4 | H_K2_sep_L1000 | +H K=2 | +H K=2 | 1 (100%) | 2 (100%) | fail | OK | OK |
| 5 | H_K3_sep_L100 | +H K=3 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 6 | H_K3_sep_L1000 | +H K=3 | +H K=3 | 1 (100%) | 3 (100%) | fail | OK | OK |
| 7 | R_K2_close_L100 | +R K=2 | +R K=2 | 2 (0%) | 1 (100%) | 1 (100%) | OK | OK |
| 8 | R_K2_close_L1000 | +R K=2 | +R K=3 | 3 (0%) | 1 (100%) | 1 (100%) | OK | WRONG |
| 9 | R_K2_sep_L100 | +R K=2 | +R K=2 | 2 (100%) | 2 (100%) | 1 (100%) | OK | OK |
| 10 | R_K2_sep_L1000 | +R K=2 | +R K=2 | 2 (95%) | 2 (100%) | 1 (100%) | OK | OK |
| 11 | R_K4_close_L100 | +R K=4 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | OK | WRONG |
| 12 | R_K4_close_L1000 | +R K=4 | +R K=3 | 3 (5%) | 1 (100%) | 1 (100%) | OK | WRONG |
| 13 | R_K4_sep_L100 | +R K=4 | +R K=3 | 3 (10%) | 2 (75%) | 1 (100%) | OK | WRONG |
| 14 | R_K4_sep_L1000 | +R K=4 | +R K=3 | 3 (100%) | 2 (100%) | 1 (100%) | OK | WRONG |
| 15 | T_K2_close_L100 | +T K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 16 | T_K2_close_L1000 | +T K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 17 | T_K2_sep_L100 | +T K=2 | +R K=1 | 1 (100%) | 1 (100%) | 1 (100%) | WRONG | WRONG |
| 18 | T_K2_sep_L1000 | +T K=2 | +T K=2 | 1 (100%) | 1 (100%) | 2 (100%) | OK | OK |

**Family correct: 11/18 (61%)  |  K correct: 6/18 (33%)**

Note: "fail" for +T on +H sep L=1000 scenarios means IQ-TREE MAST
fits crashed (numerical underflow), so +T was excluded from comparison.
Power is the proportion of B=20 bootstrap replicates that recover the
same K_best; high power with a wrong K_best means the wrong answer is
reproducible.

### Key findings

**FreeRate (+R) scenarios:**
- Family correctly identified in all 8 cases (+R always wins on BIC).
- K=2 "sep" perfectly recovered at both 100 and 1000 sites (power
  100% at L=100, 95% at L=1000).
- K=2 "close" (rates 0.7 vs 1.3): correctly recovered at L=100
  (power 0% — right K but fragile), but BIC overshoots to K=3 at
  L=1000 — more data exposes the subtle rate difference and overfits.
- K=4 scenarios: L=100 can only detect K=1 or K=3 (too few sites for
  4 rate classes). L=1000 consistently picks K=3 — BIC under-selects,
  likely a K_MAX=4 ceiling effect. Worth revisiting with K_MAX=6+.
- The L=100 vs L=1000 contrast works as intended: short alignments
  clearly lack power, producing K=1 for complex scenarios.

**GHOST (+H) scenarios:**
- +H correctly identified for "sep" at L=1000 (K=2 and K=3 both
  recovered perfectly, power 100%). At L=100, the signal is
  insufficient — all families return K=1 (tied BIC with +R).
- "close" (log_sd=0.3): too subtle at both lengths — all families
  pick K=1. Expected.
- +T consistently fails on +H data (crashes or picks K=1).
- The L=100/L=1000 contrast is clear: heterotachy requires substantial
  data to detect.

**MAST (+T) scenarios:**
- With 10 NNI moves (RF=16, 47% of max), **MAST correctly identifies
  +T K=2 at L=1000** (BIC=28369 vs +R K=1 BIC=28415, power=100%).
  This is the first successful end-to-end MAST detection in the test
  suite.
- At L=100 (50 sites per class), no signal is detected — expected.
- "close" (1 NNI move, RF=2): not detected at either length — the
  topological difference is too small relative to tree estimation error.
- The previous run with 5 NNI moves (RF=6) found no MAST signal even
  at L=1000. This confirms that substantial topological divergence
  (~50% of max RF) is needed for MAST detection with 20 taxa and
  1000 sites.

---

## Known issues and workarounds

### IQ-TREE 3.0.1 ARM BioNJ crash

The macOS ARM build of IQ-TREE 3.0.1 crashes during Jukes-Cantor distance
computation when `-t BIONJ` is used (with or without `--tree-fix`, with
any rate model). This affects:
- `fixed_tree = "NJ"` in `kpower()` (the default)
- The MAST K=1 single-tree fit

**Workaround**: use `fixed_tree = NULL` which triggers `--fast` heuristic
search instead. This is slower but avoids the crash.

The binary also required macOS quarantine removal:
```
xattr -d com.apple.quarantine ~/Desktop/Software/iqtree-3.0.1-macOS/bin/iqtree3
```

### MAST candidate tree speed

MFP (ModelFinder + full heuristic search) on each window is very slow.
The `fast_trees = TRUE` option uses `GTR+R --fast` instead, which is
much faster at the cost of candidate tree quality. Acceptable for power
assessment.

### +T tree perturbation validity (FIXED)

The original `perturb_tree()` function swapped child pointers in the
`ape` edge table, which only reordered the traversal without changing
any bipartitions (RF = 0). Replaced with `phangorn::rNNI()` which
performs proper NNI on the unrooted topology. Verified: 1 NNI move
gives RF = 2, 10 moves gives RF = 16 (max possible = 34 for 20 taxa).

### MAST bootstrap crash (FIXED)

Two bugs in kpower caused the +T bootstrap to crash:

1. **`concatenate_alignments()` dropped taxon names** (`alignment_io.R`):
   `paste0()` on named vectors drops names. Fixed by using `result[] <-`
   (subassignment) instead of `result <-` to preserve the names attribute.
   This caused IQ-TREE to reject the simulated alignments ("Please rename
   sequences").

2. **`assess_mast_power()` didn't detect `mclapply` errors**
   (`power_assessment.R`): `mclapply` returns `try-error` objects (not
   `error`), so the `inherits(x, "error")` check missed them. Failed
   replicates (character strings) were passed to `do.call(rbind, ...)`,
   producing an atomic matrix, and then `sim_ic$replicate` triggered
   "$ operator is invalid for atomic vectors". Fixed by checking for
   both `try-error` and `error`, warning and dropping failed replicates
   instead of crashing.

---

## File structure

```
kpower_tests/
  01_simulate_test_data.R     # Simulate all +R, +H, +T alignments
  02_run_kpower_tests.R       # Run kpower_survey on each, write summary
  TESTS_BRIEF.md              # This file
  run_log.txt                 # Console output from last test run
  alignments/
    test_tree.nwk             # Shared 20-taxon tree
    R_K2_sep_L100/            # One dir per scenario
      sim.phy                 # Simulated alignment
      sim_params.txt          # Generating parameters
    H_K2_sep_L100/
      sim.phy
      sim_params.txt
      class_1_tree.nwk        # Per-class trees (for +H and +T)
      class_2_tree.nwk
      ...
  results/
    summary.csv               # Full results table
    R_K2_sep_L100/
      survey_result.rds       # kpower_survey object
      ic_profile_R.pdf        # Per-family IC profile plots
      ic_profile_H.pdf
      ic_profile_T.pdf
      survey_runs/            # All IQ-TREE output
```

---

## Planned improvements

### Coverage gaps

- **+T K=3 scenarios**: add a 3-topology MAST scenario (e.g., 3 trees
  each separated by ~10 NNI moves) to test whether kpower recovers
  K=3 tree mixtures.
- **Intermediate +T signal**: the current suite jumps from RF=2
  (undetectable) to RF=16 (strong). Add a "medium" scenario (e.g.,
  5 NNI moves, RF~6) at L=2000–5000 to map the detection threshold.
- **Longer +T close alignments**: RF=2 fails at L=1000; test at
  L=5000 to determine if it is recoverable with more data.
- **K_MAX ceiling for +R K=4**: BIC always picks K=3. Rerun with
  K_MAX=6 to check whether this is a ceiling effect or a genuine BIC
  limitation.

### Cross-family confusion

- **+H absorbed by +R?** H_K2_close is selected as +R K=1 — is weak
  heterotachy genuinely indistinguishable from homogeneous, or does
  +R at higher K absorb the signal? Fit +R K=1..6 on H_close data.
- **+T misidentified as +H?** Topological heterogeneity could mimic
  branch-length heterogeneity. Systematically compare +H vs +T BIC
  on tree-mixture data at varying RF distances.

### Regression tests

- **Minimal MAST round-trip test**: simulate a K=2 MAST alignment,
  run `kpower_survey()` with `mix_types = "+T"` and B=5, assert it
  completes without error and returns power > 0. This would catch
  regressions in the concatenation and error-handling code fixed in
  `alignment_io.R` and `power_assessment.R`.
- **Alignment concatenation unit test**: verify that
  `concatenate_alignments()` preserves taxon names when inputs have
  different taxon orderings.

### Extended model families

- **Unlinked models (*R, *H, *T)**: the current suite only tests
  linked families (+). Unlinked models have separate substitution
  parameters per class and may behave differently for power and
  family identification.

### Broader simulation conditions

- **Tree shape**: all scenarios use one random 20-taxon tree. Add
  balanced and caterpillar trees to test robustness to topology.
- **More taxa**: test with 50–100 taxa to check scaling and to see
  if MAST detection improves with more phylogenetic information.
- **Production bootstrap**: increase B to 100–1000 for publishable
  power estimates with confidence intervals.
