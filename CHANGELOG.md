# Changelog

All notable changes to GEM1 are documented in this file.

This project has not yet been formally released on GitHub or archived on
Zenodo. This changelog documents the state of the codebase as prepared for
that eventual release, not a published version history.

## [Unreleased] — Manuscript / research release preparation

This is the current manuscript/research release preparation state of GEM1,
prepared alongside the manuscript describing the framework and its
NAFLD/MASLD case study. It has **not yet been published as a GitHub release
or archived on Zenodo** — this entry documents what is implemented and
reproducible in the repository at this stage, not a formal software release
announcement.

### Implemented pipeline

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
