# GEM1 — Project Roadmap & Schedule

Mohammed Mahdi · Genome-scale metabolic modeling / confidence-scored biomarker discovery
Last updated: 2026-07-18 — iMAT/tINIT growth-capability gap found and fixed (re-extraction pending execution); Steps 9-10 flagged stale until re-run

## Why the split changed
Steps 2, 4, and 5 were originally tagged generically as "Sandbox" or "Coworker + Sandbox." In practice, the cloud sandbox's network policy blocks PyPI installs and GitHub raw / GEO downloads, so any step that needs to *pull external data or packages in* has to run on Mohammed's local machine instead, via Claude Code. Sandbox handles anything self-contained; Me (local PC, via Claude Code) handles anything needing a live download or full-scale compute; Coworker (this Claude session) writes/debugs scripts, does literature and design work, verifies everything independently against raw files on Mohammed's machine, and — per Mohammed's explicit instruction (2026-07-18) — now makes all file edits directly rather than routing them through Claude Code.

## Schedule

| # | Step | Software / Tools | Status |
|---|------|-------------------|--------|
| 1 | Environment setup | Python (conda/pip), COBRApy, Escher, troppo, GEOparse, Jupyter | Done |
| 2 | Base model acquisition | BiGG Database (Human-GEM), COBRApy model I/O | Done — Human-GEM downloaded and verified locally (12,931 reactions / 8,461 metabolites / 2,848 genes) |
| 3 | Toy-model pipeline validation | COBRApy built-in test models (E. coli core, Salmonella) | Done |
| 4 | Disease dataset acquisition | GEOparse, public GEO repositories | Done — GSE89632 (63), GSE126848 (57), GSE135251 (216) fetched and disease_group-labeled |
| 5 | Context-specific model extraction (GIMME / iMAT / FASTCORE / tINIT) | troppo, corda | **Done, but iMAT/tINIT outputs are being corrected (2026-07-18).** Full 3-cohort x 10-group x 4-algorithm matrix (40 combinations) run to completion via a detached `Start-Process`. 39/40 succeeded; 1 failed cleanly (`GSE126848/healthy/fastcore` — legitimate FASTCORE "impossible to build model" outcome). GIMME's initial ~5x cross-cohort anomaly was root-caused (a real `exp_vector == -1` code bug) and fixed; all 39 models re-verified post-fix. **New finding (2026-07-18): iMAT's and tINIT's extracted reaction sets cannot sustain biomass flux on their own — a real methodological gap, not a bug in later steps. Fix implemented in `03_build_context_specific_models.py`, re-extraction of the 20 iMAT/tINIT combinations pending execution. See "Step 5 growth-capability fix" section below.** |
| 6 | Flux analyses (FBA, pFBA, FVA) | COBRApy (cobra.flux_analysis) | Pending — blocked on the growth-capability fix above; do not run Step 6 on iMAT/tINIT reaction sets until the corrected extraction completes and is independently re-verified |
| 7 | Probabilistic flux sampling | COBRApy sampling (OptGPSampler / ACHR) | Pending |
| 8 | Monte Carlo perturbation testing | Custom Python (perturbation harness) | Pending |
| 9 | Reconstruction-algorithm consensus scoring | Custom Python comparison logic | **Implemented (2026-07-17), but stale — must be regenerated once the corrected iMAT/tINIT sets are in.** Consensus computed from FASTCORE + iMAT + tINIT only; GIMME evaluated separately as a permissive-reconstruction robustness check (see decision + quantified evidence below). |
| 10 | Hierarchical confidence engine | Custom Python (weighted scoring module) | **Implemented, v1, fully audited (2026-07-17), but stale — must be regenerated after Step 9 is regenerated.** The audit's verdict on internal logic/formula correctness still holds; the underlying numbers do not, since they were computed from the pre-fix iMAT/tINIT sets. |
| 11 | Empirical weight calibration | scikit-learn, statsmodels | Pending |
| 12 | Visualization & confidence reporting | Escher, matplotlib / plotly | Pending |
| 13 | Biological interpretation & manuscript synthesis | Domain expertise + literature | Pending (Mohammed's, not Claude Code's) |

## Step 5 growth-capability fix (2026-07-18)

**Problem discovered**: while preparing Step 6, constraining the full Human-GEM model to only iMAT's or only tINIT's extracted active-reaction set and maximizing biomass flux (`MAR13082`) yielded zero — these reaction sets, as produced by the original Step 5 run, cannot sustain growth. Real methodological gap, not a Step 6 bug: neither algorithm guarantees the objective reaction survives extraction unless explicitly protected. FASTCORE and GIMME are unaffected — FASTCORE's core-reaction definition already includes biomass-support structurally, and GIMME's high permissiveness (80-96% of the network) retains growth capability incidentally.

**Decision**: of three options surfaced (re-run Step 5 for iMAT/tINIT forcing biomass essential; run Step 6 on individual/full-model sets instead of consensus; report zero-growth as a Step 6 finding rather than fixing it), chose to fix at the source — a non-growth-capable "consensus" model would be a correctness problem for every downstream step (6, 7, 8, and the confidence engine), not just a reporting nuance.

**Fix, derived independently** (troppo/cobra can't be installed in the cloud sandbox — used `scipy.optimize.linprog` directly against `models/Human-GEM.json`):
1. Built the full stoichiometric matrix (8,461 x 12,931, 55,198 nonzeros) from `Human-GEM.json`.
2. Unconstrained FBA maximizing `MAR13082` → objective 124.8681483774457 — exactly matches the pre-existing Step 5 smoke-test value (124.87), confirming the independent LP setup.
3. Minimal growth-supporting set: biomass ≥ 90% of max (same `obj_frac=0.9` convention used elsewhere), minimize total flux (L1 norm, split-variable trick) → **519 reactions**, achieving biomass flux 112.38133353970113 (exactly 90% of max). Both LPs solved to HiGHS optimal.
4. Saved as `data/context_specific_models/minimal_growth_support.json` (519 MAR IDs including `MAR13082`).

**Methodology check performed before running the fix (2026-07-18)**, prompted by outside review questioning whether forcing biomass/essential reactions into iMAT/tINIT's reconstruction deviates from published methodology (vs. a post-reconstruction gap-fill instead). Checked against primary literature:
- Original INIT (Agren et al. 2012, *PLOS Comp Biol*) does not use a biomass equation or guarantee growth capability — not the relevant precedent.
- **tINIT**, the task-driven successor actually used for Human-GEM-family models (Agren et al. 2014, *Mol Syst Biol*; also the Human-GEM/ftINIT extraction guide), works in two stages: (1) identify reactions hard-required for a predefined set of essential metabolic tasks (including growth-related tasks) and force them in as mandatory *during* reconstruction — "cannot be removed during optimization" — then (2) gap-fill only whichever individual tasks still fail after that. Forcing an essential/growth-supporting set in at reconstruction time is not a deviation from tINIT's design — it is tINIT's design. Running tINIT without any essential-reaction protection (the original Step 5 pass) was the incomplete version.
- Confirmed against troppo's own source: `tINITProperties` has a first-class `essential_reactions` parameter (floors lower bound to nonzero, excludes from removable-candidate set, always included in result) and `IMATProperties` has an analogous `core` parameter (unioned directly into iMAT's high-expression set) — purpose-built mechanisms for exactly this, not an invented workaround.
- **Conclusion**: revised the implementation to use these first-class hard-constraint parameters instead of the score-inflation originally implemented — more scientifically defensible (matches tINIT's own published design) and mechanically stronger (a real constraint, not a bias).

**Implementation** in `scripts/03_build_context_specific_models.py` (edited directly by the cloud session, revised 2026-07-18): `FORCE_REGENERATE = {"imat", "tinit"}` limits re-extraction to just those two algorithms; `run_imat()` passes `IMATProperties(..., core=growth_support_idx)`; `run_tinit()` passes `tINITProperties(..., essential_reactions=growth_support_idx)`; `main()` loads the growth-support set, maps IDs to `model_reader.r_ids` indices, and the manifest-skip check reads `if os.path.exists(out_path) and algo not in FORCE_REGENERATE:` so FASTCORE/GIMME's existing outputs are left untouched while all 20 iMAT/tINIT combinations regenerate.

**Validation plan once the corrected re-extraction completes** (per outside review's checks, still worth doing regardless of the mechanism question above): (1) independently re-verify all 10 groups' new iMAT/tINIT sets sustain ≥90% max biomass flux; (2) confirm the fix doesn't add so many reactions that it swamps the expression-driven signal; (3) confirm disease-group differences are preserved (groups shouldn't converge toward near-identical models); (4) only then regenerate Step 9/10 and re-audit.

**Status**: script changes (including the 2026-07-18 revision to `core`/`essential_reactions`) written and committed to disk on Mohammed's machine. Re-extraction execution is pending — next action is a minimal, execution-only Claude Code prompt to run the script (no design decisions left for it to make).

## Environment — hard-won fixes (do not rediscover these)

`gem1-troppo` (Python 3.10, separate from the main env — troppo's pinned `cobamp==0.2.1`/`cobra==0.24.0` won't build on 3.12) needed two rounds of real debugging beyond the base install:

1. **cobamp/optlang GLPK `_remove_variables` bug**: optlang 1.5.2's GLPK backend reads `variable._index` without checking `variable.problem is not None` first, causing a `TypeError` in `intArray___setitem__` when removing an unattached variable. Root-caused by reading the installed source directly.
2. **optlang 1.5.2 GLPK `Objective._get_expression()` desync bug**: during FastCC's repeated reversible-reaction-flip iterations, optlang's Python-side variables container desyncs from GLPK's live column count, causing `IndexError: list index out of range`. Confirmed to be a genuine optlang 1.5.2 defect (swiglpk-version hypothesis was tested and ruled out).
3. **Fix for both**: `cobamp==0.2.1` pins `optlang==1.5.2` exactly, which is the buggy version with no viable older fallback. Upgraded past the pin to `optlang==1.8.3`, which resolved both bugs cleanly (verified: identical FASTCORE output, exit code 0, no monkeypatch needed). If rebuilding `gem1-troppo` from scratch: install per the base Step 1 instructions, then `pip install "optlang==1.8.3"` on top.
4. **iMAT/tINIT need explicit Gurobi MIP caps.** Neither algorithm's troppo code sets `MIPGap`/`TimeLimit` by default and both fail to converge in practical time on Human-GEM's scale (iMAT ran >20 min without converging in testing). Fix: `MIPGap=0.02` (2%) + `TimeLimit=180` converges both in 2-3 min with a stable, usable solution. Implemented as a monkeypatch (`patch_gurobi_solve_limits()`) on `cobamp.core.optimization.LinearSystemOptimizer.optimize`, the one choke point every algorithm's solve passes through — called once at the top of `main()` in `03_build_context_specific_models.py`.
5. **Background jobs die if tied to the launching Claude Code session.** A script run as a normal background shell inside a Claude Code conversation gets killed when that Claude Code process exits (e.g. relaunching with `--continue`, closing the terminal). Fix: launch as a fully detached OS process instead — `Start-Process -FilePath "<gem1-troppo>\python.exe" -ArgumentList "scripts\03_build_context_specific_models.py" -RedirectStandardOutput "...log" -RedirectStandardError "...err.log" -WindowStyle Hidden` — this process is not a child of the terminal/Claude Code and survives either closing. Use this pattern for the pending iMAT/tINIT re-extraction and for Steps 6-8 too; they'll likely be long-running.
6. **GIMME `exp_vector == -1` bug**: caused GIMME's expression-based reaction inclusion to behave inconsistently between microarray and RNA-seq inputs. Fixed 2026-07-17; all GIMME outputs regenerated post-fix.

## Verified Step 5 result (2026-07-17, post-GIMME-fix — iMAT/tINIT counts below are pre-growth-capability-fix and will change)
40 cohort x group x algorithm combinations attempted, 39 succeeded. One failure: `GSE126848/healthy/fastcore` — "Inconsistent irreversible core reactions, impossible to build model" (legitimate FASTCORE outcome, not a crash). Manifest: `data/context_specific_models/extraction_manifest.csv`; 39 model JSONs alongside it.

Active-reaction counts by algorithm (across all successful combinations, post-GIMME-fix, pre-growth-capability-fix):
- FASTCORE: 6000-6614 (mean 6214, n=9) — final, unaffected by the growth-capability fix
- iMAT: 5078-5555 (mean 5274, n=10) — **will change**, being re-extracted
- tINIT: 7666-8246 (mean 7933, n=10) — **will change**, being re-extracted
- GIMME: 10,397-12,419 (n=10) — final, unaffected by the growth-capability fix; now platform-consistent in mechanism, but structurally far more permissive than the other three (see below); this is documented, expected GIMME behavior in the literature, not a bug.

## GIMME resolution + Step 9 consensus design (locked in 2026-07-17)

**Root cause of the original ~5x cross-cohort gap**: a real code bug (`exp_vector == -1` comparison) that behaved differently on RNA-seq vs. microarray expression value distributions. Fixed; GIMME's output is now platform-consistent in mechanism. However, post-fix GIMME is still structurally far more permissive than FASTCORE/iMAT/tINIT — it follows a "keep unless evidence says remove" logic vs. the other three's "build up a minimal justified core," a known, literature-documented property of GIMME, not a pipeline artifact.

**Quantified overlap** (computed from all 39 model JSONs; full data in `data/context_specific_models/gimme_consensus_overlap.csv` — reuse this file for a manuscript supplementary table/figure rather than recomputing):

| Cohort | Group | 3-algo consensus | GIMME size | % consensus retained by GIMME | % of GIMME that's GIMME-only |
|---|---|---:|---:|---:|---:|
| GSE89632 | healthy | 3,918 | 12,419 | 99.6% | 21.8% |
| GSE89632 | NASH | 3,728 | 12,308 | 99.8% | 20.2% |
| GSE89632 | steatosis | 3,709 | 12,279 | 99.8% | 20.3% |
| GSE126848 | healthy | 4,611 | 10,397 | 97.2% | 20.9% |
| GSE126848 | NASH | 4,090 | 10,764 | 97.8% | 19.8% |
| GSE126848 | obese_no_NAFLD | 4,012 | 10,514 | 96.4% | 18.3% |
| GSE126848 | steatosis | 4,313 | 10,442 | 97.9% | 16.4% |
| GSE135251 | healthy | 4,178 | 10,559 | 97.1% | 19.0% |
| GSE135251 | NASH | 4,003 | 10,696 | 97.0% | 17.7% |
| GSE135251 | steatosis | 3,913 | 10,662 | 97.5% | 18.3% |

Average: **98.0% of the 3-algorithm consensus is retained by GIMME** (range 96.4-99.8%); **19.3% of GIMME's own reactions are GIMME-only** (range 16.4-21.8%). GIMME's total footprint is 80-96% of the model's full 12,931 reactions.

Note: this table's "3-algo consensus" column will shift once iMAT/tINIT are re-extracted under the growth-capability fix; the table will be regenerated alongside the Step 9 re-run.

**Methodology decision (locked in with Mohammed, 2026-07-17)**: GIMME is NOT folded into the Step 9 cross-algorithm consensus as an equal-weight vote — its near-total inclusiveness would systematically inflate agreement scores for reactions the other three are more selective about, diluting the discriminative signal the consensus score is meant to capture. Instead: **primary consensus is derived from FASTCORE, iMAT, and tINIT only. GIMME is evaluated independently as a permissive-reconstruction robustness check** (which reactions stay included even under a less restrictive strategy), feeding into Step 10's hierarchical, decomposable confidence engine as its own evidence axis rather than being blended into the 3-algorithm structural consensus.

**Manuscript-language consequence, flagged, not yet changed**: the novelty dossier's and CLAUDE.md's "Core novelty" framing both currently describe a flat "multi-algorithm reconstruction consensus (iMAT/GIMME/FASTCORE/INIT)" — this needs updating to reflect the 3-algorithm structural consensus + GIMME-as-independent-robustness-check design once Mohammed signs off on the specific wording.

## Decisions locked in
- **Application framing**: GEM1 is a disease-agnostic framework; NAFLD/MASLD is its first validation case study (not its scope). See GEM1-novelty-dossier.md.
- **Base model (Step 2)**: Human-GEM, downloaded and verified locally.
- **GEO cohorts (Step 4)** — all fetched, labeled, and confirmed:
  - **GSE89632** — microarray, 63 samples (healthy 24 / steatosis 20 / NASH 19). Gene-ID mapping done (2777/2848 genes mapped).
  - **GSE126848** — RNA-seq, 57 samples (healthy 14 / obese_no_NAFLD 12 / steatosis 15 / NASH 16). Gene-ID mapping done.
  - **GSE135251** — RNA-seq, 216 samples (healthy 10 / steatosis 51 / NASH 155). Gene-ID mapping done.
- **Step 9 consensus design**: FASTCORE + iMAT + tINIT form the primary consensus; GIMME evaluated separately as a robustness check (see above).
- **Growth-capability fix (2026-07-18)**: iMAT/tINIT re-extraction hard-constrains a 519-reaction minimal growth-supporting set (LP-derived) via troppo's `core`/`essential_reactions` parameters — confirmed to match tINIT's own published task-driven design, not a pipeline-specific workaround; FASTCORE/GIMME outputs are unaffected and untouched.
- **Editing ownership (2026-07-18)**: all file edits to the GEM1 project (scripts, CLAUDE.md, docs) are made directly by the cloud coworking session via the device bridge, not routed through or self-reported by Claude Code — this followed a recurring pattern of Claude Code's self-reported CLAUDE.md writes not actually landing on disk.

## Notes carried from Step 2 discussion
- PC specs on file: Intel i5-12400F (6C/12T, ~2.5GHz), 16GB RAM, RTX 2060 Super (GPU not used by COBRApy solvers), Windows 11 Pro.
- Gurobi academic license confirmed active and working (license 2844697) — used throughout the Step 5 full run.

## Future domain roadmap (proposed, not yet locked for programming)
Mohammed proposed organizing future case studies by metabolic domain rather than by disease popularity: Hepatic (MASLD → Fibrosis → HCC) → Systemic (T2D, Obesity, Metabolic Syndrome) → Cancer (HCC, CRC, Breast, Lung) → Immune (RA, Lupus, Psoriasis — needs its own gap-check and possibly a different base model). Not yet reflected in the Step 1-13 schedule above.
