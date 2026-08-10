# GEM1 — Step 10 Confidence Engine: Validation Report

**Date**: 2026-07-17
**Auditor**: Independent verification performed in the cloud coworking session, directly against raw files on Mohammed's machine — not self-reported by the Claude Code session that built the pipeline.
**Scope**: Steps 9 (consensus scoring) and 10 (hierarchical confidence engine v1), covering all 10 cohort/group combinations, all 12,931 reactions per group (129,310 rows total).

## Purpose

Before implementing Step 11 (empirical weight calibration), confirm that Step 10's confidence engine is internally correct — no formula errors, no hidden assumptions, no duplicated/drifted logic between output files, no mishandled edge cases — so that calibration and future uncertainty layers (Steps 7, 8) are built on a verified foundation rather than an untested one.

## Method

All checks were performed by independently reloading the raw Step 5 model JSONs (`data/context_specific_models/*.json`), the Human-GEM base model (`models/Human-GEM.json`), and the Step 9/10 outputs, then recomputing every derived value from first principles and comparing against what the pipeline produced. This audit did not trust any self-reported summary from the Claude Code session that built the pipeline — every number below was independently recomputed from source data.

## Results

### 1. Reaction universe integrity
- `Human-GEM.json` contains exactly 12,931 reactions, confirming the documented base model size.
- All 39 raw Step 5 active-reaction lists (across all cohorts/groups/algorithms) are subsets of this universe — no orphan reaction IDs anywhere in the pipeline.
- All 10 confidence CSVs contain exactly this 12,931-reaction universe: no duplicates, no missing reactions, no extras.

### 2. Algorithm activity flags
- `fastcore_active`, `imat_active`, `tinit_active`, `gimme_active` in every one of the 129,310 rows across all 10 groups were independently recomputed from the raw Step 5 JSONs and compared. **Zero mismatches.**

### 3. Score formula
- `n_core_algorithms_active` correctly counts active core algorithms per reaction (FASTCORE+iMAT+tINIT, or iMAT+tINIT for the one reduced group).
- `n_core_algorithms_available` correctly reflects 3 (nine groups) or 2 (`GSE126848/healthy`, which has no FASTCORE model).
- `primary_consensus_score` and `final_confidence_score` both equal `n_core_algorithms_active / n_core_algorithms_available` exactly, on every row, with no floating-point discrepancy exceeding 1e-9. Both values are identical to each other on every row — confirming Layer 2 (GIMME) is not blended into the score.
- All scores fall within [0, 1] with no exceptions.
- Spot-checked invariant: every `strict_consensus` row scores exactly 1.0; every `no_support` row scores exactly 0.0.

### 4. Consensus tier logic
Observed (n_available, n_active) → tier mapping across all rows, confirmed consistent with the documented design in every case:

| n_available | n_active | Tier | Row count |
|---|---|---|---|
| 3 | 3 | strict_consensus | 35,864 |
| 3 | 2 | majority_consensus | 16,886 |
| 3 | 1 | minority_signal | 33,458 |
| 3 | 0 | no_support | 30,171 |
| 2 | 2 | strict_consensus | 4,611 |
| 2 | 1 | minority_signal | 3,949 |
| 2 | 0 | no_support | 4,371 |

Edge case verified: `majority_consensus` (exactly 2-of-3) correctly never appears for the 2-algorithm reduced group, since that tier is only meaningful when 3 algorithms are available. No mislabeling observed anywhere.

### 5. GIMME corroboration status
The 4-way classification (`corroborated` / `gimme_only` / `unsupported_by_gimme` / `absent_from_all`) was independently rederived from `gimme_active` × `n_core_algorithms_active` and matched exactly on all 129,310 rows. GIMME's status is stored but structurally cannot affect `final_confidence_score` (confirmed directly, not just by code inspection).

### 6. Cross-file consistency
- `confidence_manifest.csv`: every tier count and GIMME-status count, in every one of the 10 group rows, matches direct recomputation from the corresponding per-group CSV.
- `consensus_manifest.csv` (Step 9): `n_strict_consensus` matches the confidence engine's own strict-tier count in every group.
- Step 9's per-group `*_consensus.json` files: `strict_consensus_reactions` lists match the confidence engine's derived strict-consensus reaction sets exactly (zero symmetric difference), in every group checked.
- `all_groups_confidence.csv` (129,310 rows): confirmed to be an exact row-signature match of the concatenation of the 10 per-group files — not an independently (and possibly divergent) regenerated file.

### 7. Provenance fields
- `reaction_name` and `subsystem` columns match `Human-GEM.json`'s source values exactly, including which reactions have legitimately blank names upstream (5,714 of 12,931 — a base-model annotation sparsity, not a pipeline mapping defect).

## Findings requiring attention (not defects)

- **168–828 `unsupported_by_gimme` reactions per group**: reactions in the FASTCORE+iMAT+tINIT consensus that GIMME — normally the more permissive method — does not include. The audit confirms these counts and their derivation are computed correctly; whether they represent a biologically meaningful pattern (e.g., concentrated in specific subsystems — transport, arachidonic acid, steroid, and folate metabolism were the top subsystems in one profiled group) is worth investigating before Step 11's calibration, but it is not a bug in Step 9 or 10.

## Not covered by this audit

This audit verified Step 9/10's internal logic and data integrity. It did not re-verify: the correctness of the underlying troppo/Gurobi solver behavior from Step 5 (separately verified when Step 5 completed), the biological validity of Human-GEM's reaction/subsystem annotations themselves, or anything from Steps 6–8/11–13, which are not yet built.

## Verdict

**Step 9 and Step 10 are internally correct.** Every formula, every column definition, every cross-file consistency check, and the one identified edge case (the 2-algorithm reduced group) were independently verified against source data with zero discrepancies across the full dataset (129,310 rows, 10 groups). The pipeline is ready to serve as the foundation for Step 11 calibration and future Layer 3 uncertainty/robustness scores (Steps 7–8).

**Caveat carried forward, not resolved by this audit**: the confidence engine's *outputs* are verified correct; whether `v1`'s primary_consensus_score (unweighted, uncalibrated 3-algorithm agreement fraction) is itself the *right* score to calibrate in Step 11 — versus some alternative weighting — remains a modeling choice for Step 11, not something this audit can settle.
