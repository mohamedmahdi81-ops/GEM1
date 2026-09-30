# GEM1

**A confidence-aware framework for genome-scale metabolic signature discovery.**

GEM1 integrates genome-scale metabolic modeling (GSMM) with probabilistic flux
uncertainty estimation, cross-reconstruction-algorithm consensus, and a
hierarchical, empirically-calibrated confidence score, to rank candidate
disease biomarkers by how reproducible and robust they are across
independent cohorts and reconstruction methods — not by statistical
significance alone.

GEM1 is demonstrated in this repository using NAFLD/MASLD (non-alcoholic
fatty liver disease / metabolic dysfunction-associated steatotic liver
disease) as a case study. The framework itself is disease-agnostic: the
disease-specific choices (GEO cohorts, calibration biomarkers) live in the
case-study data/config layer, not hard-coded into the core pipeline logic.

> **Two method generations exist in this repository.** The current,
> manuscript-reported reaction-dependence estimator is the **paired-
> counterfactual perturbation** method described immediately below
> ("Current method"). An earlier, structurally different **unpaired**
> perturbation estimator and its downstream **0.15/0.15/0.70 calibrated
> composite confidence score** (Steps 8-11 in the "Version 1 pipeline
> components" section further down) are retained in this repository only as
> historical/provenance material — **they are not part of the final
> analysis** and are explicitly marked as suspended below.

## Current method — paired-counterfactual perturbation

For each context-specific group model (cohort x disease state), GEM1's
current reaction-dependence estimator evaluates every dropout-eligible
reaction under **paired counterfactual perturbation**: a randomized
perturbed metabolic background (a Bernoulli(q=0.30) random dropout of other
reactions) is drawn, the group's biomass objective is solved once *with* the
reaction of interest still present in that background and once *with it
also removed*, and the paired difference is recorded. This is repeated for
**M = 800** independent randomized backgrounds per cohort/group, so each
reaction's estimate is an average over 800 paired evaluations of the *same*
randomized backgrounds — not 800 unrelated single-arm perturbation runs.

- **PR_objective** (primary score) — the mean paired change in the biomass
  objective across a reaction's eligible backgrounds (`with reaction` minus
  `without reaction`, i.e. the objective lost when that reaction is dropped
  out of an already-perturbed background).
- **Non-evaluable (NE) reactions are tracked as their own category**,
  distinct from reactions that *were* evaluated but showed no detectable
  effect (`delta_objective <= TAU`, TAU=1e-6). A reaction is NE for a given
  group if it is not a dropout candidate for that group at all (i.e. it is
  outside that group's consensus active set, or it is part of the group's
  growth-support set and therefore never a candidate for removal) — NE is
  never conflated with "evaluated, zero effect" anywhere in this repository's
  data or reporting.
- **Structural consensus** (cross-reconstruction-algorithm agreement, Step 4
  below), **flux determinacy** (whether a reaction's flux is pinned or
  variable at the FBA optimum, Step 6 below), and **perturbation
  dependence** (this paired-counterfactual estimator) are reported as three
  **distinct, non-fused properties** of each reaction — this method does not
  combine them into a single composite score (see "Suspended: the 0.15/
  0.15/0.70 composite score" below).

**Production run summary** (frozen specification: M=800, perturbation
dropout probability q=0.30, TAU=1e-6, feasibility threshold obj_frac=0.99):
all **10 of 10** cohort/group analyses were accepted, **8,000/8,000**
backgrounds completed, with **zero unexpected solver statuses** and **zero
monotonicity-invariant violations** across the entire run. In total the
production analysis comprised **20,734,137 paired reaction evaluations**
(the sum of each group's per-reaction eligible-background counts across all
10 groups' `reaction_scores.csv`) — this is a count of *paired in silico LP
evaluations*, not a count of distinct biological reactions (Human-GEM has
12,931 reactions; the 5,718-reaction full-network union of dropout
candidates across all 10 groups is likewise far smaller than 20.7 million).

No external biological validation of these results is claimed anywhere in
this repository. Cross-cohort concordance, biomarker-reaction enrichment,
and pilot-vs-production reproducibility are reported as internal
statistical/numerical characterizations only (see
`data/paired_perturbation_statistics/STATISTICS_REPORT.md`).

**Suspended: the 0.15/0.15/0.70 composite score.** The Version 1 pipeline's
calibrated confidence score (Step 9, `scripts/09_calibration.py`,
`data/calibration/`) fused reconstruction consensus, flux-sampling
uncertainty, and the *old, unpaired* perturbation estimator into a single
weighted score. That old perturbation estimator was independently found not
to provide reproducible reaction-level rankings — `audit/item3_flux_sampling_perturbation.py`
found its per-reaction score essentially uncorrelated (Spearman ~0.005)
between two independent replicates at the production trial count, and a
follow-up convergence investigation (`scripts/stage2_convergence_pilot.py`)
characterizes this instability directly. The 0.15/0.15/0.70 weighting is
**suspended** — it is not recalibrated, not reapplied to the new
paired-counterfactual scores, and not part of the final analysis. No fused
score currently exists for the paired-counterfactual method; consensus,
flux, and perturbation-dependence properties are reported separately (see
above).

Scripts: `scripts/08b_paired_perturbation_production.py` (frozen production
specification, checkpointed/resumable, one run per cohort/group),
`scripts/09_paired_perturbation_pilot.py` (44-reaction curated-subset pilot
preceding the full-network production run, `scripts/proposed_pilot_reaction_subset_v2.csv`).
Post-production, read-only cross-group synthesis and statistics/figure
packages: `scripts/_paired_perturbation_synthesis.py` ->
`data/paired_perturbation_synthesis/`, and
`scripts/_stats_core_analysis.py` + `scripts/_fig_A_architecture.py`
through `scripts/_fig_F_pilot_mc.py` -> `data/paired_perturbation_statistics/`
(tables, figures in PDF/SVG/PNG with source-data CSVs, and
`STATISTICS_REPORT.md`, the full numeric write-up including every
inferential test's assumptions/limitations).

## Version 1 pipeline components (historical)

The list below documents the originally-implemented pipeline, including the
two components superseded above (Steps 8-9, and the Step 10-11 outputs
downstream of Step 9's composite score). It is retained for provenance and
because Steps 1-7 (model acquisition through flux sampling) are still used
as-is by the current method.

1. **Base model acquisition** (`scripts/01_load_base_model.py`) — Human-GEM,
   a genome-scale metabolic reconstruction of a generic human cell (see
   [`models/README.md`](models/README.md) for upstream provenance).
2. **Cohort acquisition and expression matrix construction**
   (`scripts/02_fetch_geo_cohorts.py`, `scripts/02b_build_rnaseq_expression_matrices.py`)
   — three independent, platform-diverse GEO cohorts (see below).
3. **Context-specific model extraction** (`scripts/03_build_context_specific_models.py`)
   — four independent reconstruction algorithms (FASTCORE, GIMME, iMAT,
   tINIT, via [troppo](https://github.com/BioSystemsUM/troppo)) applied to
   each cohort × disease-group combination.
4. **Reconstruction-algorithm consensus scoring** (`scripts/04_consensus_scoring.py`)
   — structural agreement across FASTCORE/iMAT/tINIT, with GIMME evaluated
   separately as a permissive-reconstruction robustness check.
5. **Hierarchical confidence engine** (`scripts/05_confidence_engine.py`).
6. **Flux analysis** (`scripts/06_flux_analysis.py`) — FBA/pFBA/FVA on each
   group's consensus-active reaction set.
7. **Probabilistic flux sampling** (`scripts/07_flux_sampling.py`) — Monte
   Carlo flux-sampling uncertainty (OptGP).
8. **_Superseded_ — unpaired perturbation robustness testing**
   (`scripts/08_perturbation_testing.py`) — the original Monte Carlo
   structural perturbation analysis. Replaced by the paired-counterfactual
   method above; kept for provenance only (see "Suspended: the 0.15/0.15/
   0.70 composite score" above for why).
9. **_Superseded_ — empirical weight calibration** (`scripts/09_calibration.py`)
   — nested leave-one-cohort-out / leave-one-biomarker-out (LOCO/LOBO)
   cross-validated weight fitting of the 0.15/0.15/0.70 composite score
   against a literature-curated biomarker panel. Suspended; not part of the
   final analysis.
10. **_Superseded_ — metadata joining and biomarker ranking**
    (`scripts/10_join_calibration_metadata.py`, `scripts/11_biomarker_ranking_by_group.py`)
    — downstream of the suspended Step 9 composite score.
11. **Figures and reporting** (`scripts/12_*.py`) — Version 1 manuscript
    figures; see `docs/MANUSCRIPT_FIGURES.md` (frozen V1 spec, not updated
    for the paired-counterfactual method).
12. **Independent arithmetic/provenance audit** (`audit/`) — a from-scratch,
    read-only re-verification of the pipeline's raw-data-to-result
    computations, run without importing any of the pipeline's own
    scoring/calibration code. See [`audit/README.md`](audit/README.md) for
    the full audit trail and findings.

## Repository structure

```
GEM1/
  README.md                 <- this file
  models/
    README.md                <- Human-GEM provenance/licensing (model files not redistributed)
  data/
    geo_cohorts/              <- cohort sample metadata + gzip-compressed expression matrices (see below)
    context_specific_models/  <- Step 3 output: per-cohort x group x algorithm active-reaction sets
    consensus_scores/         <- Step 4 output
    confidence_scores/        <- Step 5 output
    flux_analysis/            <- Step 6 output (FBA/pFBA/FVA)
    flux_sampling/             <- Step 7 output (Monte Carlo flux sampling)
    perturbation_testing/     <- Step 8 output (superseded, unpaired estimator)
    calibration/               <- Step 9-11 output (superseded, suspended composite score)
    paired_perturbation_pilot/        <- current method: 44-reaction curated pilot
    paired_perturbation_production/   <- current method: full production (10 groups; see below)
    paired_perturbation_synthesis/    <- current method: read-only cross-group synthesis tables
    paired_perturbation_statistics/   <- current method: advanced statistics + figure package
                                          (tables/, figures/, STATISTICS_REPORT.md)
    engineering_scaling_check/        <- pre-production scale/runtime validation (current method)
  results/
    figures/                   <- Version 1 manuscript figures
  scripts/                    <- the numbered V1 pipeline (01-12) plus the current-method
                                  scripts (08b, 09_paired_perturbation_pilot,
                                  stage2_convergence_pilot, _paired_perturbation_synthesis,
                                  _stats_core_analysis, _fig_A..F_*)
  audit/                      <- independent arithmetic/provenance audit (code, findings, evidence)
  docs/                        <- project documentation and design records
```

Each `paired_perturbation_production/{cohort}__{group}/` directory contains
`manifest.json` (frozen config + provenance + ACCEPTED/REJECTED status),
`checkpoint.json` (completed-background indices), and `reaction_scores.csv`
(the final per-reaction PR_objective and related statistics — this
repository's actual current-method result). The raw per-background and
per-reaction-pair Monte Carlo shards (~335MB of gzip pair data plus ~2.7MB
of per-background CSVs across all 10 groups) are **not** tracked in this
repository — see `.gitignore` — they are fully regenerable by rerunning
`scripts/08b_paired_perturbation_production.py --group <cohort>:<group>`
against the tracked context-specific models and consensus scores.

## Installation

Two separate Python environments are required (see `docs/GEM1-roadmap-schedule.md`
for the full rationale): a main environment for Steps 1-2 and 6-13 (and all
`audit/` scripts), and a dedicated environment for Step 3 (context-specific
model extraction), which depends on `troppo` and its pinned dependencies
that require Python 3.10/3.11.

```bash
conda env create -f environment.yml           # main pipeline + audit
conda env create -f environment-troppo.yml     # Step 3 extraction only
```

Both files list exact, verified package versions (the actual output of
`pip list --format=freeze` against the working environments used to produce
this repository's results) — see the comments in each file for solver notes
and known required fixes (e.g. the `optlang` upgrade needed in the troppo
environment).

A solver capable of handling both LP and MILP problems is required; the
pipeline defaults to GLPK (bundled with `cobra`) but is documented and
tested against Gurobi (with an academic license) for practical runtime on
the MILP-based reconstruction algorithms (iMAT, tINIT).

## Reproducing the pipeline

Run the numbered scripts in `scripts/` in order (01 → 12), switching conda
environments as noted above for Step 3. Long-running steps (context-specific
model extraction, flux sampling, perturbation testing) are designed to be
launched as detached background processes so they survive terminal closure —
see the in-script documentation and `scripts/README.md` for the exact
invocation pattern used during development.

To independently re-verify the pipeline's results from raw data, see
[`audit/README.md`](audit/README.md) and `audit/independent_audit.py`.

### Reproducing the current (paired-counterfactual) method

Requires Steps 1-6's outputs (base model, context-specific models, consensus
scores, growth-support set, flux analysis) already present, as produced by
the Version 1 pipeline above — the paired-counterfactual scripts reuse those
directly rather than recomputing them.

```bash
# One run per cohort:group (10 total; resumes from checkpoint.json by
# default, so an interrupted run can simply be re-invoked identically):
python scripts/08b_paired_perturbation_production.py --group GSE126848:NASH
# ... repeat for GSE126848:healthy, GSE126848:obese_no_NAFLD, GSE126848:steatosis,
#     GSE135251:NASH, GSE135251:healthy, GSE135251:steatosis,
#     GSE89632:NASH, GSE89632:healthy, GSE89632:steatosis

# 44-reaction curated-subset pilot that preceded the full-network production run:
python scripts/09_paired_perturbation_pilot.py --group GSE126848:obese_no_NAFLD \
    --reaction-subset-csv scripts/proposed_pilot_reaction_subset_v2.csv --seed 101

# Read-only, post-production (all 10 groups must be ACCEPTED first):
python scripts/_paired_perturbation_synthesis.py   # -> data/paired_perturbation_synthesis/
python scripts/_stats_core_analysis.py             # -> data/paired_perturbation_statistics/tables/
python scripts/_fig_A_architecture.py              # -> data/paired_perturbation_statistics/figures/
python scripts/_fig_B_positive_effect.py
python scripts/_fig_C_recurrence_null.py
python scripts/_fig_D_cross_cohort.py
python scripts/_fig_E_biomarker_heatmap.py
python scripts/_fig_F_pilot_mc.py
```

`08b_paired_perturbation_production.py`'s M=800/TAU=1e-6/magnitude=0.30/
obj_frac=0.99 are frozen constants, not CLI-configurable, by design (see the
script's own docstring). Each production run takes on the order of 15-20
hours of solver time on the hardware this was run on; it is designed to
survive terminal closure when launched as a detached background process
(the same long-running-step pattern documented above for Steps 3/7/8).

## Source data

### GEO cohorts

GEM1's NAFLD/MASLD case study uses three independent, platform-diverse,
publicly available GEO cohorts:

| Accession | Platform | Samples | Groups |
|---|---|---|---|
| [GSE89632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89632) | Illumina microarray | 63 | healthy 24 / steatosis 20 / NASH 19 |
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) | RNA-seq | 57 | healthy 14 / obese (no NAFLD) 12 / steatosis 15 / NASH 16 |
| [GSE135251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251) | RNA-seq | 216 | healthy 10 / steatosis 51 / NASH 155 |

**The three GEO expression matrices ARE included in this repository**, as
gzip-compressed files (originals are large — up to ~123 MB per cohort
uncompressed — and GitHub enforces a 100 MB per-file limit, so the matrices
are distributed compressed rather than omitted):

```
data/geo_cohorts/GSE89632_expression_matrix.csv.gz
data/geo_cohorts/GSE126848_expression_matrix.csv.gz
data/geo_cohorts/GSE135251_expression_matrix.csv.gz
```

Each `.csv.gz` was verified byte-identical to its original uncompressed
`.csv` via SHA-256 hash comparison (hash of the decompressed content matches
the hash of the original file exactly) before the uncompressed originals
were removed from the working tree. To decompress:

```bash
gzip -dk data/geo_cohorts/GSE89632_expression_matrix.csv.gz      # -k keeps the .gz copy
# repeat for GSE126848 / GSE135251
```

Small, derived sample metadata (disease-group labels) is also included in
`data/geo_cohorts/` directly as plain CSV.

If you would rather regenerate the matrices from scratch instead of using
the tracked compressed copies, run `scripts/02_fetch_geo_cohorts.py` (which
downloads each cohort via `GEOparse` using the accessions above) and
`scripts/02b_build_rnaseq_expression_matrices.py` (which builds the
CPM-normalized RNA-seq matrices for GSE126848/GSE135251).

### Human-GEM (base metabolic model)

GEM1 uses [Human-GEM](https://github.com/SysBioChalmers/Human-GEM), a
third-party genome-scale metabolic model maintained by SysBioChalmers, as
its base reconstruction. **`models/Human-GEM.xml` and `models/Human-GEM.json`
are not redistributed in this repository** — see
[`models/README.md`](models/README.md) for the full explanation, including
an important caveat: the exact upstream Human-GEM release/commit used by
this project's historical pipeline run could not be recovered from this
repository (the acquisition script points at the `main` branch, not a
pinned tag or commit). Run `scripts/01_load_base_model.py` to fetch the
model yourself.

## License

GEM1's own code and documentation are distributed under the MIT License —
see [`LICENSE`](LICENSE).

Human-GEM is a separate, third-party model with its own license; see
[`models/README.md`](models/README.md).

## Citing this work

See [`CITATION.cff`](CITATION.cff). A Zenodo DOI will be added once this
repository has an archived release — **no DOI exists yet**; do not cite one
until this section is updated.

### Zenodo DOI

*Placeholder — to be added once this repository is published and archived
on Zenodo. No DOI currently exists for this project.*
