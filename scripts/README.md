# GEM1 — Step 2/4/5 starter scripts

A note on where this came from: I tried to actually run this pipeline in this
cloud session, but this sandbox's network policy blocks PyPI (`pip install`)
and GitHub raw file downloads — only a narrow allowlist of hosts is reachable
here. So I could not install cobra/troppo/GEOparse or download Human-GEM in
this container. If you want future sessions to actually execute GSMM code
in the cloud sandbox rather than just producing scripts for you to run
locally, the network policy for this environment would need to allow
`pypi.org`, `files.pythonhosted.org`, and `raw.githubusercontent.com` — that's
an environment setting, see https://code.claude.com/docs/en/claude-code-on-the-web
for how egress policy is configured. Otherwise, run these three scripts on
your own machine, where you already have Python + cobra installed.

## What's here

1. **`01_load_base_model.py`** — downloads Human-GEM, loads it in COBRApy,
   verifies reaction/gene/metabolite counts, sanity-checks that it grows
   under default FBA.
2. **`02_fetch_geo_cohorts.py`** — pulls the three selected NAFLD/MASLD
   cohorts (GSE89632, GSE126848, GSE135251) via GEOparse and writes out a
   sample-metadata CSV per cohort for you to hand-label with disease groups.
3. **`03_build_context_specific_models.py`** — runs all four reconstruction
   algorithms (FASTCORE, GIMME, iMAT, tINIT/INIT) via `troppo`, per disease
   group, per cohort. This is Step 5 of the roadmap.

## Before running

```bash
pip install cobra troppo GEOparse pandas numpy
```

`troppo` was **not** in your original install list (you installed cobra,
escher, GEOparse, etc., but not troppo/corda) — the roadmap's Step 5 needs it
for GIMME/iMAT/FASTCORE/INIT.

## Run order

```bash
python 01_load_base_model.py
python 02_fetch_geo_cohorts.py
# --> open data/geo_cohorts/*_sample_metadata.csv, add a 'disease_group'
#     column by hand for each sample (healthy / steatosis / NASH / fibrosis)
python 03_build_context_specific_models.py
```

## Known gaps to close before trusting the output

- **Gene ID mapping**: Human-GEM genes are identified by Ensembl/Entrez IDs
  (check `[g.id for g in model.genes][:5]` after Step 1). GSE89632 is a
  probe-level Illumina microarray and GSE126848/GSE135251 are RNA-seq —
  each needs its own probe/gene-ID -> model-gene-ID mapping before the
  expression scores in Step 5 are meaningful. This script does not do that
  mapping for you.
- **troppo API drift**: the exact `ReconstructionWrapper` call signatures in
  `03_build_context_specific_models.py` are correct as of troppo's public
  example notebooks at the time this was written, but I could not run this
  script to confirm against your installed troppo version (network-blocked
  in this session). Test one cohort/group/algorithm combination first.
- **GSE135251 sample count**: I could not verify the exact current sample
  count for this series from this sandbox; the fetch script will print it —
  record it in the novelty dossier once you have it.
