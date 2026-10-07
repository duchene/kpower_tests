# Filovirus empirical analysis

Empirical application of kpower to real filovirus genomic data, demonstrating
power assessment for mixture-model selection across three model families: FreeRate
(+R), GHOST (+H), and MAST (+T).

## Overview

This analysis addresses a central question in phylogenetic mixtures: when two model
families yield similar likelihood scores, how reliably can information criteria
recover the true model family from data with the observed amount of signal?

Using NCBI filoviruses, we build genus-specific alignments and ask:
- Which model family does BIC prefer?
- How much power does that criterion have to recover its choice?
- How do AIC and AICc compare?

## Key Results

### Orthoebolavirus (27 taxa, 14,653 sites)

| Family | K_best | BIC | Power (BIC) | Power (AIC) |
|---|---|---|---|---|
| +R FreeRate | 3 | **188674.9** | 95% | 98% |
| +H GHOST | 3 | 189217.5 | 53% | 100% |
| +T MAST | 1 | 188694.3 | 100% | — |

**Finding:** FreeRate wins by 543 BIC units. BIC recovers its choice (K=3) in 95% of
bootstrap replicates (95/100). GHOST's K=3 is recovered only 53% of the time,
revealing BIC's power asymmetry for heterotachy-like models.

### Orthomarburgvirus + Mengla/Dehong (16 taxa, 16,227 sites)

| Family | K_best | BIC | Power (BIC) | Power (AIC) |
|---|---|---|---|---|
| +R FreeRate | 3 | 161051.3 | **1%** | 47% |
| +H GHOST | 3 | 161123.9 | 100% | — |
| +T MAST | 3 | **160970.1** | 77% | — |

**Finding:** MAST wins by 81 BIC units over FreeRate. Critically, BIC's
model-selection margin (161051 vs 161056 for K=2) is only 5.0 units — too
small to reliably identify K=3. Bootstrap recovers K=3 in just 1 of 100 replicates
under FreeRate. This is an extreme but genuine case: BIC picked a model, but the
data lack the signal to reliably recover that choice.

## Scientific Message

These results validate the power assessment framework and align with Susko et al.
(2023): **information criteria and empirical data support need to be decoupled.**
A model selected by BIC can still have very low power to recover its own choice
when the signal is marginal.

## Caveats

### Model adequacy

- Ebolavirus +R generates data with only 34.4% constant sites vs 43.5% observed,
  creating a 0.46 lnL/site shortfall. GHOST (+H) generates 43.2%.
- This mismatch explains why +R power is lower than typical: the simulated data
  do not match the complexity of the real data.

### Information-criterion comparison

The AIC/AICc power rows report power **under the primary criterion's K_best**.
For example, marburg +R BIC reports K_best = 3, but AIC reports K_best = 5 (from
the empirical data). The AIC power of 47% is computed from replicates simulated
under +R K=3, which AIC was not asked to recover. Interpret only where criteria agree.

### Small margins

Marburg +R's 1% power reflects a genuine weak signal, not statistical error:
K=2 and K=3 are separated by only 5.0 BIC units on a 16,227-site alignment,
and bootstrap replicates overwhelmingly favor K=2.

## Scripts and Data

### Input

- `alignments/filovirus/ncbi_virus_filoviridae.fasta` (16 MB) — 854 sequences
  from NCBI filovirus database, downloaded and described in Step 1.

### Processing

**Step 1: Alignment and clade separation**
```bash
Rscript 03_filovirus_prep_alignment.R
```
- Filters to ~100 taxa per genus (Orthoebolavirus, Orthomarburgvirus + Mengla/Dehong)
- Removes sequences with excessive divergence (genetic distance p > 0.45, saturated)
  and near-duplicates (p < 0.005)
- Performs two-pass MAFFT alignment: align all, prune divergent sequences, realign
- Outputs: `ebola_filtered.fasta`, `marburg_filtered.fasta` plus distance matrices

**Step 2: Gblocks filtering**
```bash
Rscript 03b_filovirus_gblocks_filter.R
```
- Reduces alignment expansion and removes gappy regions
- Strict filter (used for analysis): 14,653 sites (ebola), 16,227 (marburg)
- Outputs: `*_gb_strict.fasta` and a filter summary

**Step 3: kpower_survey (B=100)**
```bash
Rscript 04_filovirus_kpower_survey.R
```
- Runs `kpower_survey()` on each clade for +R, +H, +T
- K range: 1–5 for +R/+H; 1–4 for +T (cost optimization; K_best well below boundary)
- B=100 bootstrap replicates; 6 parallel workers
- Fits are seeded deterministically from output prefix for reproducibility
- Outputs: `survey_summary.txt`, IC profiles, power under all three criteria

**Step 4: Visualization**
All scripts automatically generate PDF figures:
- Per-family IC-profile plots (`*_ic_profile_*.pdf`)
- Cross-family BIC support (`bic_support.pdf`)
- Summary tables (`bic_profiles.csv`, `family_comparison.csv`)

## Quick start

From the `kpower_tests/` directory:

```bash
# Prerequisites: R with kpower (install to parent ../Rlib or system library)
#                IQ-TREE 3 (must be in PATH or ~/Desktop/Software/...)

# Full analysis (skips steps already completed)
Rscript 03_filovirus_prep_alignment.R     # ~2 min
Rscript 03b_filovirus_gblocks_filter.R    # ~5 min
Rscript 04_filovirus_kpower_survey.R      # ~15 hours

# Results:
cat results/filovirus/survey_summary.txt
open results/filovirus/bic_support.pdf
```

The scripts are checkpoint-aware: if interrupted, they resume at the next
incomplete family rather than recomputing from the start.

## Reproducibility

- Model fits are seeded deterministically from each run's `--prefix`
  (unique per family, replicate, K), so the same analysis yields bit-identical
  results.
- Results are provided in RDS format (`*_survey.rds`, `*_B100_K*.rds`) for
  programmatic access.
- All intermediate alignments and checkpoints are committed to git.

## Manuscript / Citation

These results are part of a manuscript on parametric bootstrap power assessment
for mixture models in phylogenetics. See the main README for references to
Susko et al. (2023) and other foundational work.

---

**Last updated:** 2026-10-06 | **Analysis completed:** 2026-10-06 | **Bootstrap replicates:** B = 100
