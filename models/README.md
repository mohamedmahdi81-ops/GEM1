# `models/` — Human-GEM (upstream third-party model)

This directory is where GEM1's base metabolic model, **Human-GEM**, is
expected to live locally (`Human-GEM.xml` and `Human-GEM.json`).

## Human-GEM is not redistributed in this repository

`Human-GEM.xml` and `Human-GEM.json` are **excluded from GEM1's Git
repository** (see `.gitignore`). Reasons:

1. **They are a third-party model**, not GEM1's own output. Human-GEM is
   developed and maintained independently by
   [SysBioChalmers](https://github.com/SysBioChalmers/Human-GEM) and should
   be obtained from its own upstream source, not vendored into a downstream
   project's repository.
2. **Their exact upstream release/commit cannot currently be established
   from this repository.** GEM1's acquisition script
   (`scripts/01_load_base_model.py`) downloads the model from:

   ```
   https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.xml
   ```

   This points at the `main` branch, not a pinned git tag, commit SHA, or
   release version — so the file that was actually downloaded for GEM1's
   historical pipeline run cannot be tied to a specific, citable Human-GEM
   release after the fact. **No version number, commit hash, or DOI is
   claimed for the historically-used copy** — stating otherwise would be a
   guess, not a verified fact.
3. **Size** — `Human-GEM.xml` is ~41 MB and `Human-GEM.json` ~10 MB,
   unnecessary bulk for a Git repository when the model is easily
   re-obtained.

## How to obtain Human-GEM

Run `scripts/01_load_base_model.py` from the repository root (with the main
GEM1 conda environment active — see the top-level `README.md`). This will:

1. Download `Human-GEM.xml` from the SysBioChalmers/Human-GEM `main` branch
   into this directory.
2. Load and sanity-check it via `cobra` (reaction/metabolite/gene counts,
   default FBA feasibility).
3. Save a JSON copy (`Human-GEM.json`) for faster reloading in later
   pipeline steps.

**Note**: because the download script points at a live, unpinned `main`
branch, re-running it today may produce a different file than the one
historically used to generate the results in `data/` and `results/` in
this repository, if Human-GEM has been updated upstream since. If your
downloaded model's reaction/metabolite/gene counts differ substantially
from those recorded in this project's documentation (12,931 reactions /
8,461 metabolites / 2,848 genes, as observed at the time of the historical
run), the upstream model has changed — that is expected behavior of an
unpinned reference, not a bug in this download script.

## Upstream source, attribution, and licensing

- **Project**: Human-GEM — a genome-scale metabolic model of a generic
  human cell.
- **Maintainer**: [SysBioChalmers](https://github.com/SysBioChalmers/Human-GEM)
  (Chalmers University of Technology), part of the broader
  [Metabolic Atlas](https://metabolicatlas.org/) initiative.
- **Source repository**: https://github.com/SysBioChalmers/Human-GEM
- **License**: Human-GEM is distributed by SysBioChalmers under its own
  license, stated in the upstream repository — consult
  https://github.com/SysBioChalmers/Human-GEM for the current, authoritative
  license terms before redistributing or reusing the model itself. GEM1's
  own MIT license (see the top-level `LICENSE` file) applies only to GEM1's
  own code and documentation, **not** to Human-GEM.
- **Citation**: if you use Human-GEM, cite the upstream Human-GEM
  publication(s) as directed by the SysBioChalmers repository — GEM1 does
  not restate that citation here to avoid transcribing bibliographic
  details that should come directly from the authoritative source.

If you need to pin an exact Human-GEM version for strict reproducibility
going forward (rather than relying on `main`), obtain a tagged release
directly from the SysBioChalmers/Human-GEM GitHub releases page and adjust
`scripts/01_load_base_model.py`'s `MODEL_URL` accordingly.
