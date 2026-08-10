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

## Methodological components implemented in this repository

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
8. **Perturbation robustness testing** (`scripts/08_perturbation_testing.py`)
   — Monte Carlo structural perturbation analysis.
9. **Empirical weight calibration** (`scripts/09_calibration.py`) — nested
   leave-one-cohort-out / leave-one-biomarker-out (LOCO/LOBO) cross-validated
   weight fitting against a literature-curated biomarker panel.
10. **Metadata joining and biomarker ranking**
    (`scripts/10_join_calibration_metadata.py`, `scripts/11_biomarker_ranking_by_group.py`).
11. **Figures and reporting** (`scripts/12_*.py`).
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
    perturbation_testing/     <- Step 8 output (Monte Carlo perturbation)
    calibration/               <- Step 9-11 output: calibrated confidence scores, biomarker rankings
  results/
    figures/                   <- final manuscript figures
  scripts/                    <- the numbered pipeline, 01 through 12
  audit/                      <- independent arithmetic/provenance audit (code, findings, evidence)
  docs/                        <- project documentation and design records
```

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
