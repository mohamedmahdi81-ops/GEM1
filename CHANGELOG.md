# Changelog

All notable changes to GEM1 are documented in this file.

This project has not yet been formally released on GitHub or archived on
Zenodo. This changelog documents the state of the codebase as prepared for
that eventual release, not a published version history.

## [1.0.1] - 2026-08-11

- Corrected author name and title/abstract metadata in CITATION.cff;
  removed tool-specific references from code comments and documentation.

## [Unreleased] — Manuscript / research release preparation

This is the current manuscript/research release preparation state of GEM1,
prepared alongside the manuscript describing the framework and its
NAFLD/MASLD case study. It has **not yet been published as a GitHub release
or archived on Zenodo** — this entry documents what is implemented and
reproducible in the repository at this stage, not a formal software release
announcement.

### Paired-counterfactual perturbation replaces the unpaired estimator (2026-09)

The manuscript's reaction-dependence method has changed. The original
"structural perturbation robustness" item below (Step 8) and the
"empirical weight calibration" item (Step 9, the 0.15/0.15/0.70 composite
score) are **superseded and suspended** — see the top-level `README.md`'s
"Current method" section for the full rationale and methodology. Summary:

- The original (unpaired) perturbation estimator was found not to provide
  reproducible reaction-level rankings (`audit/item3_flux_sampling_perturbation.py`:
  Spearman ~0.005 between independent replicates at the production trial
  count; further characterized by `scripts/stage2_convergence_pilot.py`).
- Replaced by a **paired-counterfactual** estimator
  (`scripts/08b_paired_perturbation_production.py`): each reaction is
  evaluated within the *same* randomized perturbed background, with and
  without that reaction, across M=800 backgrounds per cohort/group
  (dropout probability q=0.30, TAU=1e-6, obj_frac=0.99). All 10 cohort/group
  analyses were accepted: 8,000/8,000 backgrounds completed, zero unexpected
  solver statuses, zero monotonicity violations, 20,734,137 total paired
  reaction evaluations.
- The 0.15/0.15/0.70 calibrated composite score is **suspended** — not
  recalibrated, not reapplied to the new scores, not part of the final
  analysis. Structural consensus, flux determinacy, and perturbation
  dependence are reported as separate, non-fused properties instead.
- A 44-reaction curated-subset pilot (`scripts/09_paired_perturbation_pilot.py`)
  preceded the full-network production run; a read-only post-production
  synthesis (`scripts/_paired_perturbation_synthesis.py`) and an advanced
  statistics/figure package (`scripts/_stats_core_analysis.py`,
  `scripts/_fig_A_architecture.py` through `_fig_F_pilot_mc.py`, written up
  in `data/paired_perturbation_statistics/STATISTICS_REPORT.md`) characterize
  the production results numerically. No external biological validation is
  claimed.

### Implemented pipeline (Version 1 — see above for the current method)

- Base metabolic model acquisition (Human-GEM, via SysBioChalmers upstream
  source — see `models/README.md` for provenance and licensing notes).
- Acquisition and processing of three independent, platform-diverse GEO
  cohorts (GSE89632, GSE126848, GSE135251) for the NAFLD/MASLD case study.
- Context-specific metabolic model extraction across four independent
  reconstruction algorithms (FASTCORE, GIMME, iMAT, tINIT) for each
  cohort × disease-group combination.
- Reconstruction-algorithm consensus scoring (structural agreement across
  FASTCORE/iMAT/tINIT; GIMME evaluated separately as a permissive-
  reconstruction robustness check).
- Hierarchical confidence engine combining reconstruction consensus, flux
  sampling uncertainty, and structural perturbation robustness.
- Probabilistic flux sampling (Monte Carlo, OptGP) and structural
  perturbation robustness testing (Monte Carlo) as independent uncertainty
  evidence layers.
- Empirical weight calibration via nested leave-one-cohort-out /
  leave-one-biomarker-out (LOCO/LOBO) cross-validation against a
  literature-curated, multi-source biomarker panel (5 axes: serine/glycine
  one-carbon metabolism, urea cycle, branched-chain amino acid catabolism,
  de novo lipogenesis, choline/phosphatidylcholine metabolism).
- Figure generation and biomarker ranking reporting.

### Reproducibility materials

- An independent arithmetic/provenance audit (`audit/`) that re-derives
  key pipeline quantities directly from available source data and
  repository outputs — without
  importing the pipeline's own scoring/calibration code — and documents
  discrepancy investigations, methodology-sensitivity checks, and open
  questions found in the course of that audit. See `audit/README.md` for
  the full trail.
- Repository preparation for publication: MIT license, `.gitignore`
  excluding raw/regenerable/private material, documentation of upstream
  Human-GEM provenance, and a `CITATION.cff` for the repository.

### Known limitations, as documented by the audit

- Several manuscript-reported statistics (Table 3 comparator IQR values,
  Section 3.8 enrichment p-values, a subsystem-level Kendall's W figure)
  could not be traced to a source script anywhere in this repository as of
  this release; the audit records first-source recomputations and
  discrepancy investigations for these rather than treating them as
  independently confirmed.
- The exact upstream Human-GEM release/commit used historically cannot be
  recovered from this repository (see `models/README.md`).
- This changelog entry does not itself assert which specific reported
  results are or are not correct — see `audit/README.md` and
  `audit/audit_report.csv` for the item-by-item findings.

### Not yet done

- No GitHub release, tag, or Zenodo DOI has been created for this project.
- `environment.yml` / `environment-troppo.yml` are not yet finalized in
  the repository.
