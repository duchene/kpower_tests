# kpower validation tests

Simulation-based validation suite for
[kpower](https://github.com/duchene/kpower), an R package for assessing
the statistical power of mixture-model selection in phylogenetics.

## Overview

This repository holds two complementary analyses:

### 1. Simulation-based validation (TESTS_BRIEF.md, scripts 01–02)

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

### 2. Empirical analysis on real filovirus genomes (FILOVIRUS_ANALYSIS.md, scripts 03–04)

Real-world application of kpower to NCBI filovirus genomic data. Demonstrates
the power assessment framework on two genera-restricted alignments, addressing
the key finding: **BIC selected a model family, but the data lack the signal to
reliably recover that choice** (Susko et al. 2023 applies here). Results include:

- **Ebola (27 taxa):** FreeRate K=3 wins; BIC recovers it 95% of the time; GHOST
  K=3 only 53%. Power asymmetry for heterotachy-like signals.
- **Marburg (16 taxa):** MAST K=3 wins; BIC's margin over FreeRate K=3 is only 5.0
  units; FreeRate power is 1% (data lack the signal). Canonical example of weak
  signal despite model selection.

See [FILOVIRUS_ANALYSIS.md](FILOVIRUS_ANALYSIS.md) for full results, caveats, and
how to reproduce the analysis.

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
  README.md                                 # This file
  TESTS_BRIEF.md                            # Simulation-based validation details
  FILOVIRUS_ANALYSIS.md                     # Empirical analysis details
  
  ## Simulation-based validation (01–02)
  01_simulate_test_data.R                   # Simulate all +R, +H, +T alignments
  02_run_kpower_tests.R                     # Run kpower_survey on each, write summary
  
  ## Empirical filovirus analysis (03–04)
  03_filovirus_prep_alignment.R             # Download, filter, and align filoviruses
  03b_filovirus_gblocks_filter.R            # Remove gappy regions with Gblocks
  04_filovirus_kpower_survey.R              # Run kpower_survey (B=100) on each genus
  
  alignments/                               # Alignment data (not tracked in git)
    test_tree.nwk
    <scenario>/sim.phy                      # Test suite alignments
    filovirus/
      ncbi_virus_filoviridae.fasta          # Input NCBI data (16 MB)
      ebola_*.fasta                         # Intermediate steps
      marburg_*.fasta
  
  results/                                  # Survey output (not tracked in git)
    summary.csv                             # Test suite summary
    <scenario>/survey_result.rds            # Test suite results
    filovirus/
      survey_summary.txt                    # Human-readable summary
      bic_profiles.csv                      # All empirical IC scores
      family_comparison.csv                 # K_best + power table
      bic_support.pdf                       # Cross-family BIC profile plot
      *_ic_profile_*.pdf                    # Per-family bootstrap figures
      *_survey.rds                          # Full kpower_survey results
      *_B100_K*.rds                         # Per-family checkpoints
```

## References

- Crotty, S.M., et al. (2020). GHOST: Recovering Historical Signal
  from Heterotachously Evolved Sequence Alignments. *Systematic
  Biology*, 69(2), 249-264.
- Duchene, D.A., Duchene, S., & Ho, S.Y.W. (2017). New Statistical
  Criteria Detect Phylogenetic Bias Caused by Compositional
  Heterogeneity. *Systematic Biology*, 66(5), 769-785.
- Susko, E., Leigh, J.W., & Doron-Faigenboim, A. (2023). Comparing
  the fit of complex models of sequence evolution: Cross-validation
  and the discrete Kolmogorov-Smirnov test. *Molecular Biology and
  Evolution*, 40(1), msac239.
- Woodhams, M.D., et al. (2024). MAST: Mixture Across Sites and Trees.
  *Systematic Biology*, 73(2), 375-391.
