Paste the message below as your first prompt to Claude Code, run from inside the GEM1 folder (`cd` there first, then start `claude`).

---

Read CLAUDE.md in this folder fully before doing anything else — it has the project context, locked-in decisions, and environment requirements. Then:

1. Confirm the folder structure from CLAUDE.md exists (docs/, models/, data/geo_cohorts/, data/context_specific_models/, scripts/, notebooks/) and create any missing pieces. docs/ and scripts/ should already be populated.

2. Check my existing conda environment: run `conda env list` and `python --version` inside whichever env has cobra/GEOparse/etc. already installed, and tell me its name and Python version. This becomes `gem1-main` in CLAUDE.md's terms — rename or just note the actual name, don't recreate it if it already has the right packages.

3. Create a second, separate conda environment called `gem1-troppo` pinned to Python 3.10 (not 3.11, not 3.12 — troppo's dependencies won't build on 3.12, see CLAUDE.md for why). Install troppo, GEOparse, pandas, numpy into it. Confirm it installs cleanly — if it doesn't, stop and show me the exact error rather than guessing around it.

4. In `gem1-main`, run `scripts/01_load_base_model.py` to download Human-GEM and verify it loads (reaction/gene/metabolite counts should roughly match what's noted in CLAUDE.md — flag me if they're substantially different, the upstream model may have been updated).

5. In `gem1-main`, run `scripts/02_fetch_geo_cohorts.py` to pull the three GEO cohorts (GSE89632, GSE126848, GSE135251). This produces a `*_sample_metadata.csv` per cohort with each sample's raw GEO characteristics text, but no `disease_group` column yet — do NOT guess and fill in disease-group labels yourself. Instead, show me the distinct `characteristics` values per cohort so I can tell you exactly how to map them to healthy/steatosis/NASH/fibrosis groups, then add the column once I've confirmed.

6. Once cohorts are labeled, switch to `gem1-troppo` and run `scripts/03_build_context_specific_models.py` — but first test it against just ONE cohort/group/algorithm combination (comment out the rest or add a quick filter) and confirm it produces a sensible non-empty model before running the full matrix across all 3 cohorts x groups x 4 algorithms. The troppo API calls in that script were written without being able to test-run them, so expect to need to debug the exact `ReconstructionWrapper` method signatures against your installed troppo version.

7. As each of these steps completes, update the "Roadmap status" section in CLAUDE.md so it reflects reality, and tell me a summary of what happened, what's left, and any decisions you need from me (like the disease-group labeling in step 5 above).

Do not proceed past step 3 if the `gem1-troppo` environment doesn't install cleanly — report the exact error back to me instead of trying workarounds I haven't approved.
