# `scripts/` — the GEM1 pipeline

This directory contains GEM1's full, executed pipeline. Every script listed
below has been run to completion, and its output is committed under `data/`
and `results/` elsewhere in this repository. The pipeline's results have
additionally been independently re-verified from raw data by the audit in
`audit/` (see `audit/README.md`) — that audit is a from-scratch
recomputation, not a re-run of these scripts, and documents exactly which
quantities were confirmed, which were revised, and which remain open
questions.

## Environments

Two separate conda environments are required — see the top-level
`README.md` and `environment.yml` / `environment-troppo.yml` for exact,
verified package versions. In short: Steps 1-2 and 6-13 run in the main
environment; Step 3 (context-specific model extraction) requires a
dedicated Python 3.10/3.11 environment because of `troppo`'s pinned
dependencies.

## Pipeline, in order

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_load_base_model.py` | Download Human-GEM, load via COBRApy, sanity-check reaction/gene/metabolite counts and default FBA feasibility. |
| 2 | `02_fetch_geo_cohorts.py` | Fetch the three GEO cohorts (GSE89632, GSE126848, GSE135251) via `GEOparse`; write per-cohort sample metadata. |
| 2b | `02b_build_rnaseq_expression_matrices.py` | Build CPM-normalized RNA-seq expression matrices for the two RNA-seq cohorts. |
| 3 | `03_build_context_specific_models.py` | Extract context-specific models via four independent reconstruction algorithms (FASTCORE, GIMME, iMAT, tINIT) per cohort × disease group, using `troppo`. **Run this one in the `gem1-troppo` environment.** |
| 4 | `04_consensus_scoring.py` | Reconstruction-algorithm consensus scoring (structural agreement across FASTCORE/iMAT/tINIT; GIMME evaluated separately). |
| 5 | `05_confidence_engine.py` | Hierarchical confidence engine, Layer 1 (reconstruction consensus). |
| 6 | `06_flux_analysis.py` | FBA/pFBA/FVA on each group's consensus-active reaction set. |
| 7 | `07_flux_sampling.py` | Monte Carlo flux sampling (OptGP) for flux-uncertainty evidence. |
| 8 | `08_perturbation_testing.py` | **Superseded.** Original (unpaired) Monte Carlo structural perturbation robustness testing — see "Current method: paired-counterfactual perturbation" below for its replacement. |
| 9 | `09_calibration.py` | **Superseded.** Empirical weight calibration (0.15/0.15/0.70) via nested LOCO/LOBO cross-validation against the literature-curated biomarker panel — suspended, see below. |
| 10 | `10_join_calibration_metadata.py` | Join calibrated scores with reaction metadata. |
| 11 | `11_biomarker_ranking_by_group.py` | Per-group biomarker ranking summary. |
| 12 | `12_figure*.py`, `12_supplementary_s1_escher_maps.py` | Manuscript figures and supplementary outputs. |

Long-running steps (3, 7, 8) were run as detached background processes so
they survive terminal closure — see the in-script documentation and
`docs/GEM1-roadmap-schedule.md` for the exact launch pattern used during
development, and `audit/run_sequential.ps1` for a working example of
chaining multiple long steps sequentially (required on memory-constrained
machines — running these concurrently can exhaust available RAM/page file).

## Current method: paired-counterfactual perturbation

Replaces Step 8's unpaired estimator (found not to provide reproducible
reaction-level rankings — see `audit/item3_flux_sampling_perturbation.py`,
Spearman ~0.005 between independent replicates) and does not feed into
Step 9's suspended composite score. Reuses Steps 1-6's outputs directly. See
the top-level `README.md`'s "Current method" section for the full
methodology, frozen-parameter list, and production summary statistics.

| Script | Purpose |
|---|---|
| `08b_paired_perturbation_production.py` | Frozen-specification (M=800, TAU=1e-6, magnitude=0.30, obj_frac=0.99) paired-counterfactual production run. One invocation per `cohort:group` (10 total); checkpointed/resumable. |
| `09_paired_perturbation_pilot.py` | 44-reaction curated-subset pilot preceding the full-network production run; configurable `--n-backgrounds`/`--seed` for the pilot's own convergence checks. Uses `proposed_pilot_reaction_subset_v2.csv`. |
| `stage2_convergence_pilot.py` | Follow-up investigation of how many Monte Carlo trials the *old* (Step 8) estimator needs for a stable per-reaction estimate, given the reproducibility failure the audit found. Motivates the paired-counterfactual redesign; does not itself produce the redesign. |
| `_paired_perturbation_synthesis.py` | Read-only, post-production (all 10 groups ACCEPTED) cross-group synthesis tables -> `data/paired_perturbation_synthesis/`. |
| `_stats_core_analysis.py` | Read-only advanced statistics: sparsity/occurrence, positive-effect magnitude, eligibility-aware recurrence permutation nulls, cross-cohort reproducibility, biomarker-reaction permutation tests, Monte Carlo precision, pilot-vs-production Bland-Altman -> `data/paired_perturbation_statistics/tables/`. Full methodology, exact results, and every test's assumptions/limitations are written out in `data/paired_perturbation_statistics/STATISTICS_REPORT.md`. |
| `_fig_A_architecture.py` ... `_fig_F_pilot_mc.py` | One script per manuscript figure (A-F), each reading only from `_stats_core_analysis.py`'s tables (or `reaction_scores.csv` directly) and writing PDF+SVG+PNG plus a source-data CSV to `data/paired_perturbation_statistics/figures/`. |

These are kept in the public release despite the `_`-prefix dev-script
naming convention used elsewhere in this directory (see `.gitignore`'s
explicit exceptions) — same rationale as `_floor_sensitivity_check.py`
below: real, cited manuscript outputs, not one-off diagnostics.

## `_floor_sensitivity_check.py`

A standalone sensitivity-analysis script (not part of the main numbered
sequence) that reruns Step 9's calibration weight grid search at alternative
floor values (0.05, 0.0) instead of the production floor (0.15). Its output,
`data/calibration/floor_sensitivity_results.json`, is a real, cited result
(see `docs/MANUSCRIPT_FIGURES.md`) — kept in the public release for that
reason, unlike the other one-off diagnostic/dev scripts used during
development (not included in this repository; see `.gitignore`).

## Known open items

See `audit/README.md` for a full, itemized account of what has and hasn't
been independently re-verified, including several manuscript-adjacent
figures that could not be traced to any script in this repository (treated
as open findings, not silently assumed correct) and a documented weight-
calibration instability finding (items 14/16/17 in the audit).
