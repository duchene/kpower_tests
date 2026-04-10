# kpower validation tests

Simulation-based validation suite for
[kpower](https://github.com/duchene/kpower), an R package for assessing
the statistical power of mixture-model selection in phylogenetics.

## Overview

These tests simulate sequence alignments under known mixture models
(FreeRate +R, GHOST +H, MAST +T) with controlled parameters, then run
`kpower_survey()` to check whether the correct model family and number
of mixture classes (K) are recovered.

The test suite is designed to:

1. **Verify family identification** -- does `kpower_survey()` select the
   model family that generated the data?
2. **Verify K recovery** -- does it recover the true number of mixture
   categories?
3. **Characterise power boundaries** -- how do alignment length, class
   separation, and topological divergence affect detection?
4. **Catch regressions** -- confirm that bug fixes in kpower's MAST
   bootstrap pipeline remain functional.

## Quick start

### Prerequisites

- R (>= 4.0) with packages: `kpower`, `ape`, `phangorn`, `processx`,
  `ggplot2`
- [IQ-TREE 3](http://www.iqtree.org/) (tested with 3.0.1)

### Running

```bash
# 1. Simulate all test alignments (~1 min)
Rscript 01_simulate_test_data.R

# 2. Run kpower_survey on each alignment (~2-3 hours)
Rscript 02_run_kpower_tests.R
```

Results are written to `results/summary.csv` and per-scenario RDS/PDF
files in `results/<scenario>/`.

## Test design

**18 scenarios** across three model families, two signal strengths
("sep" and "close"), and two alignment lengths (100 and 1000 sites):

| Family | Scenarios | Signal |
|---|---|---|
| **+R** (FreeRate) | 8 (K=2 and K=4) | Rate categories with known proportions |
| **+H** (GHOST) | 6 (K=2 and K=3) | Per-branch lognormal heterotachy |
| **+T** (MAST) | 4 (K=2) | Distinct topologies via NNI perturbation |

**Shared parameters:**
- 20 taxa, random tree scaled to mean root-to-tip = 0.45 subs/site
  (Duchene et al. 2017, Syst Biol 66:769-785)
- GTR model with strong transition bias (A-C=5, A-G=1, ...)
- B = 20 bootstrap replicates per family (low, for speed)
- K_MAX = 4

See [TESTS_BRIEF.md](TESTS_BRIEF.md) for full scenario details,
parameter tables, results, and discussion.

## Current results

Family correct: **11/18 (61%)** | K correct: **6/18 (33%)**

Detection succeeds for well-separated scenarios at L=1000:

| Condition | L=100 | L=1000 |
|---|---|---|
| +R K=2 sep (rates 0.1 vs 2.5) | OK | OK |
| +H K=2 sep (log_sd=1.5) | -- | OK |
| +H K=3 sep (log_sd=1.5) | -- | OK |
| +T K=2 sep (RF=16/34) | -- | OK |

Short alignments (L=100) generally lack power. "Close" scenarios
(subtle class differences) are not detected -- as expected for a
power assessment tool.

## Bugs found and fixed

These tests identified two bugs in kpower that broke the MAST (+T)
bootstrap pipeline:

1. **`concatenate_alignments()` dropped taxon names** -- `paste0()` on
   named vectors drops names; fixed with subassignment (`result[] <-`).
2. **`assess_mast_power()` missed `mclapply` errors** -- `mclapply`
   returns `try-error` objects, not `error`; fixed to check both classes
   and degrade gracefully.

Both bugs only affected the MAST bootstrap phase (Phase 2). MAST
empirical fits (Phase 1) and all +R/+H results were unaffected.

## Repository structure

```
kpower_tests/
  README.md                     # This file
  TESTS_BRIEF.md                # Detailed design, results, and discussion
  01_simulate_test_data.R       # Simulate all +R, +H, +T alignments
  02_run_kpower_tests.R         # Run kpower_survey on each, write summary
  alignments/                   # Simulated data (not tracked in git)
    test_tree.nwk
    <scenario>/sim.phy
  results/                      # Survey output (not tracked in git)
    summary.csv
    <scenario>/survey_result.rds
```

## References

- Duchene, D.A., Duchene, S., & Ho, S.Y.W. (2017). New Statistical
  Criteria Detect Phylogenetic Bias Caused by Compositional
  Heterogeneity. *Systematic Biology*, 66(5), 769-785.
- Crotty, S.M., et al. (2020). GHOST: Recovering Historical Signal
  from Heterotachously Evolved Sequence Alignments. *Systematic
  Biology*, 69(2), 249-264.
- Woodhams, M.D., et al. (2024). MAST: Mixture Across Sites and Trees.
  *Systematic Biology*, 73(2), 375-391.
