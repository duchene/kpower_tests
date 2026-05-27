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

**48 scenarios** = 3 model families x 2 values of K (2 and 4) x 2 signal
strengths ("sep" and "close") x 4 alignment lengths (100, 300, 1000,
3000 sites):

| Family | Scenarios | Signal |
|---|---|---|
| **+R** (FreeRate) | 20 (K=2 and K=4) | Rate categories with known proportions |
| **+H** (GHOST) | 20 (K=2 and K=4) | Per-branch lognormal heterotachy |
| **+T** (MAST) | 20 (K=2 and K=4) | Distinct topologies via NNI perturbation |

**Shared parameters:**
- 20 taxa, random tree scaled to mean root-to-tip = 0.45 subs/site
  (Duchene et al. 2017, Syst Biol 66:769-785)
- GTR model with strong transition bias (A-C=5, A-G=1, ...)
- B = 20 bootstrap replicates per family (low, for speed)
- K_MAX = 5

See [TESTS_BRIEF.md](TESTS_BRIEF.md) for full scenario details,
parameter tables, results, and discussion.

## Current results

> Note: the numbers below are from the earlier 5-length grid (max
> L=6000, 60 scenarios). The suite now runs 4 lengths (max L=3000, 48
> scenarios); rerun `02`-`04` to regenerate. The L=6000 column will
> disappear and the per-family tallies drop from /20 to /16.

Family correct: **32/60 (53%)** | K correct: **20/60 (33%)**

By family (family-correct / K-correct out of 20 each):

| Family | Family correct | K correct | Notes |
|---|---|---|---|
| **+R** | 19/20 (95%) | 10/20 | Family almost always right; K caps at 3 for true K=4 |
| **+H** | 9/20 (45%) | 7/20 | Detected only for "sep" at L>=300 |
| **+T** | 4/20 (20%) | 3/20 | Rarely detected, but rock-solid when it is (power ~0.95) |

Detection is driven by both signal strength and length. Family-correct
rate by length: 25% (L=100), 50% (300), 50% (1000), 75% (3000), 67%
(6000). By signal: 73% for "sep" vs 33% for "close".

Short alignments (L=100) lack power across the board. "Close" scenarios
(subtle class differences) are mostly undetected -- as expected for a
power tool -- with one striking exception: +T close (RF=2) reaches 90%
power at L=6000, showing even a tiny topological split is recoverable
given enough data.

See `results/power_heatmap.pdf` for the full power grid and
`results/selection_K{2,4}.pdf` for per-bootstrap family selection.

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
