# GEM1 — Manuscript Figures Specification

> **Design Freeze**: This document defines the scientific intent of each manuscript figure.
> Implementation may refine aesthetics and layout, but must not alter the scientific claims,
> hypotheses, or reviewer concerns without explicit revision of this specification.

Status: Complete — Figures 1–5, Supplementary Figures S1–S2, and the Manuscript Claim Map are
all specified below (frozen scientific content verbatim from the approved specification; S2 was
initially truncated mid-transmission across three attempts before arriving intact). Everything
below this notice was independently verified against the actual project filesystem on
2026-07-30 (pipeline outputs, scripts, manifests) — nothing was inferred from filenames or
assumed from CLAUDE.md's prose alone without checking the underlying files. Three open items
remain (see "Open Items" at the end) — none of the frozen scientific content is affected by them.

---

## FIGURE 1 — GEM1 Framework Overview

**Role in Manuscript**: Introduce the GEM1 framework and provide readers with the complete
computational workflow before presenting any results.

**Scientific Claim**: GEM1 is an integrated, disease-agnostic framework that combines
context-specific reconstruction, uncertainty analysis, consensus scoring, and hierarchical
confidence estimation into a single reproducible pipeline.

**Hypothesis**: Combining multiple complementary computational layers provides a more reliable
interpretation of metabolic alterations than isolated GSMM analyses.

**Evidence**: The workflow integrates all computational stages into a coherent framework whose
performance is evaluated in subsequent figures.

**Reviewer Concern Addressed**: "What exactly is GEM1, and how is it different from a standard
GSMM workflow?"

**Section**: Main text.

### Verified pipeline support

**Required pipeline outputs** (all confirmed present on disk, not inferred):
| Stage | Directory/file | Confirmed |
|---|---|---|
| Base model | `models/Human-GEM.json`, `models/Human-GEM.xml` | Exists — 12,931 reactions / 8,461 metabolites / 2,848 genes (`version` field = `1`) |
| GEO cohort acquisition | `data/geo_cohorts/` | Exists (11 entries) |
| Context-specific reconstruction | `data/context_specific_models/` | Exists (46 entries incl. `extraction_manifest.csv`) |
| Flux analysis (FBA/pFBA/FVA) | `data/flux_analysis/` | Exists (11 entries incl. `flux_manifest.csv`) |
| Flux sampling | `data/flux_sampling/` | Exists (31 entries incl. `flux_sampling_manifest.csv`) |
| Perturbation testing | `data/perturbation_testing/` | Exists (31 entries incl. `perturbation_manifest.csv`) |
| Consensus scoring | `data/consensus_scores/` | Exists (12 entries incl. `consensus_manifest.csv`) |
| Confidence engine | `data/confidence_scores/` | Exists (13 entries incl. `confidence_manifest.csv`) |
| Empirical calibration | `data/calibration/` | Exists (13 entries incl. `calibration_manifest.csv`, `calibration_weights_and_validation.json`) |

**Input datasets / input files**: `models/Human-GEM.json` (base GSMM); `data/geo_cohorts/` (GSE89632, GSE126848, GSE135251 — all 3 confirmed present per CLAUDE.md's "Decisions locked in" and the directory listing above). Confirmed present, not inferred from naming alone — directory contents enumerated directly.

**Analysis scripts** (all confirmed to exist as files; confirmed by module docstring / manifest cross-check they produce the outputs above, not merely by filename):
`scripts/01_load_base_model.py`, `02_fetch_geo_cohorts.py`, `02b_build_rnaseq_expression_matrices.py`, `03_build_context_specific_models.py`, `04_consensus_scoring.py`, `05_confidence_engine.py`, `06_flux_analysis.py`, `07_flux_sampling.py`, `08_perturbation_testing.py`, `09_calibration.py`.

**Panel layout** (proposed — not frozen content, open to revision):
- **Panel A**: Horizontal pipeline schematic, Steps 1–11, grouped into 4 stages matching the Scientific Claim's own four named components: (1) data/model acquisition, (2) context-specific reconstruction (4 algorithm icons: iMAT/GIMME/FASTCORE/tINIT feeding into one base model), (3) uncertainty analysis (flux sampling + perturbation testing shown as parallel branches), (4) consensus + hierarchical confidence engine + empirical calibration (shown converging into one `calibrated_confidence_score` output).
- **Panel B** (inset): simplified 4-layer diagram of the confidence engine itself (Layer 0 raw calls → Layer 1 consensus → Layer 2 GIMME corroboration (annotation-only, dashed to show it doesn't feed the score) → Layer 3 calibrated combination), reused visually in later figures for continuity.

**Statistical analysis**: N/A — this is a conceptual/schematic figure summarizing pipeline architecture; it makes no quantitative claim and requires no statistical test.
- Primary statistical test: N/A
- Secondary statistical test: N/A
- Effect size: N/A
- Confidence interval: N/A
- Multiple-testing correction: N/A
- Significance threshold: N/A
- Pipeline step producing these statistics: N/A (figure depicts Steps 1–11 collectively; no single step's statistics are shown here)

**Reproducibility**:
- Pipeline stage: All (Steps 1–11)
- Generation script: **Not yet implemented.** Step 12 ("Visualization & confidence reporting") is listed as Pending in `CLAUDE.md`'s roadmap; no figure-rendering script exists anywhere in `scripts/` (confirmed via repo-wide search for figure/plot/viz-named files — none found).
- Input files: see table above
- Manifest version: no single versioned pipeline manifest exists; each step has its own dated manifest CSV (see table above). Most recent regeneration batch spans 2026-07-23 (Step 5 re-extraction) through 2026-07-30 (Step 11 calibration).
- Model version: `Human-GEM.json`'s own `version` field = `1`; no separate upstream release tag is recorded in this repo beyond the reaction/metabolite/gene counts CLAUDE.md already documents (12,931 / 8,461 / 2,848 — all confirmed to match the local file exactly).
- Calibration version: production weights `(w1=0.15, w2=0.15, w3=0.70)`, fixed 2026-07-30 (`data/calibration/calibration_weights_and_validation.json`).
- Output directory: N/A (no figure output exists yet)
- Verification status: Underlying Steps 1–11 outputs independently re-audited at multiple points (see `CLAUDE.md`'s per-step methodology sections and the "Step 5 re-extraction" / "Step 11 calibration" sections for the full audit trail). The figure itself has not been generated or verified, since it doesn't exist yet.

**Caption outline**: "Overview of the GEM1 pipeline. (A) Four-stage workflow from raw transcriptomic cohorts through context-specific reconstruction (four independent algorithms), uncertainty quantification (flux sampling and perturbation robustness), and a hierarchical, empirically-calibrated confidence engine, applied here to three independent NAFLD/MASLD cohorts. (B) Structure of the confidence engine's evidence layers."

**Dependencies**: Steps 1–11 all complete (confirmed above); Step 12 (figure generation itself) not started.

---

## FIGURE 2 — Positioning GEM1 Relative to Existing Frameworks

**Role in Manuscript**: Establish methodological novelty.

**Scientific Claim**: To our knowledge, no previously published framework integrates
context-specific reconstruction, multi-algorithm consensus, uncertainty quantification,
empirically calibrated confidence scoring, and cross-cohort validation within a single workflow.

**Hypothesis**: No existing published framework provides the complete methodological combination
implemented by GEM1.

**Comparison Table** (neutral, falsifiable criteria):

| Criterion | Mardinoglu 2014 | Mardinoglu 2016 | BayFlux | GEMsembler | GEM1 |
|---|---|---|---|---|---|
| Context-specific reconstruction from transcriptomics | Yes (INIT) | Yes (iMAT) | No | No | Yes |
| Multiple reconstruction algorithms combined (1) | No | No | No | Yes | Yes |
| Flux uncertainty quantification (2) | No | Partial (FVA bounds) | Yes | No | Yes |
| Cross-algorithm consensus scoring | No | No | No | Yes | Yes |
| Hierarchical/decomposable confidence score | No | No | No | No | Yes |
| Empirically calibrated via cross-validated weight fitting | No | No | No | No | Yes (nested LOCO/LOBO) |
| Disease biomarker application | Yes (NAFLD) | Yes (NAFLD) | No | No | Yes (NAFLD) |
| Cross-cohort validation | Partial (1 validation cohort) | No | No | No | Yes (3 cohorts) |

**Footnote (1)**: GEMsembler integrates outputs from different whole-model reconstruction tools
across microbial genomes (demonstrated on Lactiplantibacillus plantarum and Escherichia coli,
outperforming gold-standard models in auxotrophy and gene essentiality predictions). GEM1
integrates outputs from different context-specific extraction algorithms (iMAT/GIMME/
FASTCORE/tINIT) applied to transcriptomic data against one shared human base model. Different
consensus paradigms, not the same capability applied to different organisms.

**Footnote (2)**: BayFlux uses Bayesian inference and MCMC to quantify metabolic flux uncertainty
from 13C-labeling data. GEM1 quantifies uncertainty in context-specific metabolic predictions
through flux sampling and Monte Carlo perturbation robustness, evaluated against
transcriptomics-derived models, not isotope-labeling data. Different uncertainty frameworks
answering different questions.

**Reviewer Concern Addressed**: "Is GEM1 genuinely novel, or is it only a combination of existing
methods?" — addressed by making the novelty claim about the combination, not unique ownership
of any single row; BayFlux and GEMsembler each outperform GEM1 on their own specialty axis,
disclosed rather than obscured.

**Citations**:
- Mardinoglu et al. 2014, Nature Communications: https://www.nature.com/articles/ncomms4083
- Mardinoglu et al. 2016, Nature Communications: https://www.nature.com/articles/ncomms9994
- Backman et al., BayFlux, PLOS Computational Biology 2023: https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011111
- Matveishina et al., GEMsembler, mSystems 2025: https://journals.asm.org/doi/10.1128/msystems.00574-25

**Section**: Main text. Implementation note: populate the table only from this literature-verified
specification. Do not independently reinterpret or extend literature claims.

### Verified support

**Required pipeline outputs**: None. This is a qualitative literature-comparison table, not derived from GEM1's own computed data.

**Input datasets / input files**: The two cited external papers were independently re-verified this session, not just trusted from the specification text:
- BayFlux (PLOS Comp Biol) — fetched directly: title confirmed as "BayFlux: A Bayesian method to quantify metabolic Fluxes and their uncertainty at the genome scale"; confirmed to use Bayesian inference + MCMC (the paper's own "AcMet" algorithm) on 13C-labeling data. Matches the frozen footnote exactly.
- GEMsembler (mSystems) — direct fetch returned HTTP 403; re-verified via search instead. Confirmed: Matveishina et al., mSystems 2025, "GEMsembler: consensus model assembly and structural comparison of genome-scale metabolic models across tools improve functional performance," demonstrated on *E. coli* and *L. plantarum*, outperforming gold-standard models on auxotrophy and gene essentiality predictions. Matches the frozen footnote exactly.
- Mardinoglu 2014/2016 citations were not independently re-fetched this session (already treated as established in `docs/GEM1-novelty-dossier.md`, which cites and summarizes both directly — see Dependencies).

**Analysis scripts**: None — this figure requires no pipeline computation.

**Panel layout** (proposed): Single full-width table panel (the comparison table above) with the two footnotes set as table notes beneath it, not separate panels.

**Statistical analysis**: N/A — qualitative feature comparison, no statistical test applies.
- Primary statistical test: N/A
- Secondary statistical test: N/A
- Effect size: N/A
- Confidence interval: N/A
- Multiple-testing correction: N/A
- Significance threshold: N/A
- Pipeline step producing these statistics: N/A

**Reproducibility**:
- Pipeline stage: N/A (literature synthesis, not a pipeline stage)
- Generation script: None — table is manually curated from the frozen specification; no script produces it, nor should one (it is not data-derived).
- Input files: `docs/GEM1-novelty-dossier.md` (confirmed exists; independently read this session — its own novelty framing and Mardinoglu 2014/2016 methodology summary are consistent with, though not identical in wording to, this figure's frozen table; the dossier does not itself cite BayFlux/GEMsembler by URL, so treat this figure's citation list as the authoritative one, not the dossier's).
- Manifest version: N/A
- Model version: N/A
- Calibration version: N/A
- Output directory: N/A (no figure output exists yet; Step 12 pending, as in Figure 1)
- Verification status: Both non-Mardinoglu citations (BayFlux, GEMsembler) independently verified this session against their actual publications, as detailed above.

**Caption outline**: "Feature comparison of GEM1 against the two foundational NAFLD GSMM studies and two methodologically related but non-overlapping tools (BayFlux: flux uncertainty from isotope tracing; GEMsembler: cross-tool structural consensus in microbial GEMs). GEM1's contribution is the combination, not unique ownership of any individual capability (see footnotes)."

**Dependencies**: `docs/GEM1-novelty-dossier.md` (background/rationale); no pipeline script dependencies.

---

## FIGURE 3 — Cross-Cohort Confidence Analysis

**Role in Manuscript**: Present the primary biological results.

**Scientific Claim**: Confidence scores derived by GEM1 remain consistent across independent
NAFLD/MASLD cohorts, demonstrating reproducibility.

**Hypothesis**: Independent transcriptomic cohorts should converge on similar high-confidence
metabolic alterations.

**Evidence**: Cross-cohort confidence distributions and biomarker rankings.

**Reviewer Concern Addressed**: "Are these findings cohort-specific, or are they reproducible?"

**Section**: Main text.

### Verified support

**Required pipeline outputs**: `data/calibration/all_groups_calibrated.csv` (per-reaction `calibrated_confidence_score` for all 10 groups across all 3 cohorts — confirmed exists) and `data/calibration/calibration_manifest.csv` (confirmed exists; contents read directly, all 10 groups present with `mean_calibrated_confidence_score_active_only` clustering tightly at 0.5001 for every group).

**Input datasets / input files**: Same as Figure 1's Step 5–11 chain; specifically the 3 GEO cohorts (GSE89632, GSE126848, GSE135251) as the "independent cohorts" this figure's claim is about.

**Analysis scripts**: `scripts/09_calibration.py` produces the per-reaction scores. **No dedicated "biomarker ranking across cohorts" script or output file exists yet** — this specific table/plot would need to be derived from `all_groups_calibrated.csv` plus the biomarker mapping already defined in `09_calibration.py` (`BIOMARKERS_FULL`), reusing its `biomarker_metric`-style max-aggregation logic. The underlying data fully supports computing this; it is not yet a packaged, pre-computed artifact, and Step 12 would need to add this derivation (not just plot an existing file).

**Panel layout** (proposed):
- **Panel A**: Distribution (violin or box plot) of `calibrated_confidence_score` among active reactions, one per cohort (3 cohorts, pooling each cohort's groups or shown per-group).
- **Panel B**: Biomarker × cohort/group grid (5 biomarkers × 10 groups) showing each biomarker's max-aggregated percentile rank, color-coded — the direct visualization of "biomarker rankings" the Evidence line calls for.

**Statistical analysis** — a deliberate methodological note, not a filled-in-by-default template: a naive Kruskal-Wallis (or ANOVA) test directly on `calibrated_confidence_score` across the 3 cohorts was considered and **rejected** as this figure's primary statistic. Reason: the rank-normalization step inside the scoring formula (Step 11) mathematically forces each group's mean active-reaction score toward ≈0.5 regardless of true biological consistency — confirmed directly: all 10 groups' `mean_calibrated_confidence_score_active_only` in `calibration_manifest.csv` sit in a 0.5001-wide band. Presenting that as "proof of cross-cohort consistency" would be circular (a property of the scoring construction, not independent evidence). The methodologically sound quantitative backing for this figure's claim is the **nested leave-one-cohort-out (LOCO) cross-validation already computed in Step 11** (see Figure 5's Reproducibility) — LOCO directly measures whether biomarker rankings generalize to a held-out, unseen cohort, which is precisely this figure's claim.
- Primary statistical test: Nested LOCO cross-validation (already computed, Step 11; not a new test) — mean outer-held-out score 0.883 (range 0.841–0.912 across the 3 cohorts)
- Secondary statistical test: Not yet computed — a formal per-biomarker rank-consistency statistic (e.g., Kendall's W across cohorts on biomarker rank order) would strengthen this figure and should be added at Step 12, not fabricated here
- Effect size: LOCO range (0.841–0.912) as the effective spread across cohorts
- Confidence interval: Not computed for LOCO fold scores (n=3 cohorts is too small for a meaningful CI on the mean; this should be stated as a limitation, not glossed over)
- Multiple-testing correction: N/A (no multi-comparison test performed yet)
- Significance threshold: N/A
- Pipeline step producing these statistics: Step 11 (`scripts/09_calibration.py`, nested LOCO)

**Reproducibility**:
- Pipeline stage: Steps 9, 10, 11
- Generation script: `scripts/09_calibration.py` produces the underlying per-reaction data; no figure-rendering script exists yet (Step 12 pending)
- Input files: `data/calibration/all_groups_calibrated.csv`, `data/calibration/calibration_manifest.csv`, `data/calibration/calibration_weights_and_validation.json`
- Manifest version: `calibration_manifest.csv`, regenerated 2026-07-30 alongside the rest of Step 11
- Model version: `Human-GEM.json` version `1` (unchanged from Figure 1)
- Calibration version: production weights `(0.15, 0.15, 0.70)`, 2026-07-30
- Output directory: `data/calibration/`
- Verification status: Underlying per-reaction scores spot-checked (3 random reactions, manual recomputation matched stored values exactly; all inactive reactions confirmed `== 0`) during Step 11 implementation. The specific "biomarker ranking across cohorts" derivation for this figure has **not** been separately computed or verified yet — flagged above as still needing to be built.

**Caption outline**: "Calibrated confidence scores across three independent NAFLD/MASLD cohorts (GSE89632, GSE126848, GSE135251). (A) Score distributions per cohort. (B) Ranking of the 5 literature-curated calibration biomarkers across all 10 disease groups, supported quantitatively by nested leave-one-cohort-out cross-validation (mean held-out score 0.883)."

**Dependencies**: Step 11 complete; a new derivation step (biomarker-ranking-by-cohort table) still needed before this figure can be built, per the gap noted above.

---

## FIGURE 4 — Full-Isoform vs. Single-Proxy Validation

**Role in Manuscript**: Justify a key methodological design decision.

**Scientific Claim**: Aggregating information across all mapped isoforms provides consistently more
robust confidence estimates across validation experiments than using a single reaction proxy.

**Hypothesis**: Full-isoform aggregation yields superior validation performance.

**Evidence**: LOCO validation (0.883 vs. 0.556) and LOBO validation (0.879 vs. 0.537) comparing
both approaches.

**Reviewer Concern Addressed**: "Why was full-isoform aggregation chosen instead of a simpler
proxy?"

**Section**: Main text.

### Verified support

**Required pipeline outputs**: `data/calibration/calibration_weights_and_validation.json` — read directly this session; confirmed exact values: `full_isoform` LOCO mean = 0.8829 (rounds to 0.883 ✓), LOBO mean = 0.8787 (rounds to 0.879 ✓); `single_proxy` LOCO mean = 0.5555 (rounds to 0.556 ✓), LOBO mean = 0.5371 (rounds to 0.537 ✓). All four cited numbers match the stored artifact exactly.

**Input datasets / input files**: Same Step 9/10 outputs as Figure 3, run twice under the two biomarker-aggregation variants defined in `scripts/09_calibration.py` (`BIOMARKERS_FULL` vs. `BIOMARKERS_SINGLE_PROXY`).

**Analysis scripts**: `scripts/09_calibration.py` (produces both variants' nested-CV results in one run).

**Panel layout** (proposed):
- **Panel A**: Grouped bar chart, LOCO mean scores, full-isoform vs. single-proxy (3 underlying cohort folds each).
- **Panel B**: Grouped bar chart, LOBO mean scores per biomarker, full-isoform vs. single-proxy (5 biomarkers each) — explicitly showing the `choline_pc` outlier (single-proxy: 0.0001) as a labeled callout, not hidden.

**Statistical analysis** — computed fresh this session from the matched per-fold data already stored in the JSON (not previously reported as a formal test anywhere):
- Primary statistical test: Wilcoxon signed-rank test, paired by matched outer-test fold (3 LOCO + 5 LOBO = 8 pairs, full-isoform vs. single-proxy), ties dropped per standard convention → n=6 non-zero pairs, **W=0.0, p=0.03125** (all 6 non-tied comparisons favor full-isoform; the 2 exact ties are `urea_cycle` and `bcaa`, where the isoform dropped for single-proxy was never the max-scoring reaction in the first place, so removing it changed nothing — confirmed by inspecting the matched fold scores directly, not assumed).
- Secondary statistical test: One-sided sign test on the same 8 pairs (ties excluded, n=6): p=0.0156 — consistent with the Wilcoxon result, reported as a simple corroborating check given the very small n.
- Effect size: Mean paired difference (full-isoform − single-proxy) across all 8 folds = **+0.336** (driven heavily by the `choline_pc` fold's +0.969 difference; median paired difference = +0.297, less influenced by that outlier). **Corrected 2026-08-05** (was stated as +0.176 with median +0.31): `scripts/12_figure4_isoform_vs_proxy.py`, built as the permanent version of this session's ad hoc check, recomputed both numbers directly from `calibration_weights_and_validation.json`'s 8 matched fold values and got +0.336 (median +0.297), not +0.176. Traced before correcting, not just patched: this is not a git repository and no earlier draft of this file exists on disk to diff against, so the original derivation itself is unrecoverable — but the Wilcoxon W (0.0), Wilcoxon p (0.03125), and sign-test p (0.0156) reported alongside the wrong mean all reproduce exactly from the same fold data, which rules out "stale/different underlying data" as the cause (three of the four numbers in this bullet were already correct). Several alternative interpretations of "mean paired difference" (LOBO-only, LOCO-only, ties-excluded, ratio-based, inner-selection-score-based) were tried and none reproduce +0.176 either. By elimination this looks like a one-off arithmetic/transcription slip made when the number was typed in ad hoc, not a methodology difference or a computation against superseded data — consistent with this section's own original admission that it was "computed ad hoc this session, not yet a saved script."
- Confidence interval: Not computed — n=6–8 is too small for a stable bootstrap CI on the paired difference; stating this as a limitation rather than fabricating a CI.
- Multiple-testing correction: N/A — single pre-specified comparison (full-isoform vs. single-proxy), not part of a multi-comparison family.
- Significance threshold: α=0.05 (both tests clear it, p=0.03125 and p=0.0156 respectively, though n is small and this should be stated alongside the result, not left implicit).
- Pipeline step producing these statistics: Step 11 (`scripts/09_calibration.py` for the raw scores; the Wilcoxon/sign-test computation itself was run ad hoc this session directly against the stored JSON, not by a saved script — worth adding as a permanent check if this figure is finalized).

**Reproducibility**:
- Pipeline stage: Step 11
- Generation script: `scripts/09_calibration.py`
- Input files: `data/calibration/calibration_weights_and_validation.json`
- Manifest version: same artifact as Figure 3, 2026-07-30
- Model version: `Human-GEM.json` version `1`
- Calibration version: both `full_isoform` and `single_proxy` variants stored in the same JSON, 2026-07-30
- Output directory: `data/calibration/`
- Verification status: The four headline numbers (0.883/0.556/0.879/0.537) verified to match the stored JSON exactly, character for character on rounding. The Wilcoxon/sign-test statistics have now been independently re-verified by a second, permanent method (`scripts/12_figure4_isoform_vs_proxy.py`, Step 12, 2026-08-05): W=0.0, p=0.03125, sign-test p=0.0156 all reproduced exactly; the effect-size line's mean was found wrong in that spot-check and is corrected above.

**Caption outline**: "Held-out validation comparing full-isoform aggregation (MAX over all mapped reactions per biomarker) against a single-best-reaction proxy. (A) LOCO: 0.883 vs. 0.556. (B) LOBO: 0.879 vs. 0.537, per biomarker, highlighting a catastrophic single-proxy failure for choline/PC (0.0001) absent under full-isoform aggregation. Wilcoxon signed-rank p=0.031 across 8 matched folds."

**Dependencies**: Step 11 complete.

---

## FIGURE 5 — Calibration and Validation

**Role in Manuscript**: Demonstrate that the confidence engine is empirically calibrated rather than arbitrarily
parameterized.

**Scientific Claim**: The GEM1 confidence weights were selected through an empirical nested
cross-validation procedure rather than arbitrary assignment. The adopted production weighting
achieves performance statistically indistinguishable from the best-performing alternatives
while preserving decomposability and interpretability across the confidence engine's evidence
components.

**Hypothesis**: Confidence weights can be selected objectively using validation data, producing a
transparent and reproducible scoring system without requiring the numerically optimal solution
to be adopted when multiple solutions perform equivalently.

**Evidence**: Nested calibration procedure; LOCO validation; LOBO validation; comparison among
candidate weight vectors (floor=0.15 vs. 0.05 vs. 0.0); demonstration that multiple
high-performing solutions exist within a narrow performance range (~0.003-0.005 difference).

**Reviewer Concern Addressed**: "Were the confidence weights chosen arbitrarily or tuned after
inspecting the results?"

**Section**: Main text.

### Verified support

**Required pipeline outputs**: `data/calibration/calibration_weights_and_validation.json` (production weights, LOCO/LOBO fold-level records) and `CLAUDE.md`'s "Step 11 empirical weight calibration" section (2026-07-30), which documents the floor-sensitivity comparison in full, including the corrected miscitation. **The floor=0.05/0.0 comparison numbers themselves are not in the JSON file** (that diagnostic, `scripts/_floor_sensitivity_check.py`, only prints to stdout — it does not persist a structured output file) — `CLAUDE.md` is currently the only saved, durable record of those specific numbers. Recommend persisting them to a JSON/CSV alongside the existing calibration outputs before this figure is finalized, rather than relying on a prose paragraph as the source of record.

**Input datasets / input files**: Same as Figure 4, plus the ad hoc floor=0.05/0.0 re-runs.

**Analysis scripts**: `scripts/09_calibration.py` (production weights, floor=0.15); `scripts/_floor_sensitivity_check.py` (floor=0.05 and floor=0.0 re-runs — confirmed exists; note the underscore-prefixed naming convention this project uses for one-off diagnostics rather than permanent pipeline steps, consistent with `_validate_step6.py`, `_audit_step10_post_growthfix.py`, etc.).

**Panel layout** (proposed):
- **Panel A**: Schematic of the nested CV structure (outer fold / inner fold, LOCO and LOBO shown side by side) — reuses the same visual language as Figure 1 Panel B for continuity.
- **Panel B**: Grouped bar/point chart, LOCO and LOBO mean scores at floor=0.15/0.05/0.0, with error bars or a bracket annotation showing the ~0.003–0.005 spread is small relative to between-fold variability.
- **Panel C**: Final production weight composition, `(0.15, 0.15, 0.70)`, shown as a simple stacked/pie breakdown, annotated with the interpretability rationale (not the empirical-optimum claim).

**Statistical analysis** — computed fresh this session, paired by matched fold (floor=0.15 vs. floor=0.05, same 3 LOCO + 5 LOBO folds):
- Primary statistical test: Wilcoxon signed-rank test, floor=0.15 vs. floor=0.05, n=8 matched folds, no ties this time: **W=8.0, p=0.1953** — not statistically significant, directly supporting the "statistically indistinguishable" claim (this is a real computed result, not an assumption; a non-significant p-value here is the expected/desired outcome for this figure's claim, not a null result to be hidden). **Corrected 2026-08-05** (was stated as p=0.1797): `scripts/12_figure5_calibration_validation.py`, the permanent version of this session's ad hoc check, reproduced W=8.0 exactly and every effect-size number exactly (mean paired diff 0.0036, SD 0.0059, between-fold SD 0.037 LOCO/0.123 LOBO) but got p=0.1953, not 0.1797. Checked whether this was a different-method artifact before concluding otherwise: neither script (`09_calibration.py` nor `_floor_sensitivity_check.py`) contains any Wilcoxon call at all — the original number, like Figure 4's, was computed ad hoc outside any saved script, so there is no alternate "method used by another script" to attribute the gap to. Explicitly tried scipy's `'exact'`, `'approx'` (with and without continuity correction) methods against these exact 8 pairs: exact/auto=0.1953, approx-no-correction=0.1614, approx-with-correction=0.1834 — none reproduce 0.1797. Same conclusion as Figure 4: an unrecoverable one-off arithmetic slip at ad hoc computation time, not a methodology disagreement, and it does not change this figure's qualitative claim (both the old and corrected p-value clear "not significant at alpha=0.05" identically). **Method now standardized project-wide**: any future paired Wilcoxon signed-rank test in this pipeline should use `scipy.stats.wilcoxon(..., method='auto')` (scipy's default), which resolves to the exact distribution whenever n<50 with no ties/zeros, as used here and in Figure 4.
- Secondary statistical test: Direct comparison of the paired-difference magnitude against between-fold variability: mean paired difference = 0.0036 (SD 0.0059), versus between-fold SD of 0.037 (LOCO) and 0.123 (LOBO) at floor=0.15 — the floor effect is roughly 6–20× smaller than ordinary fold-to-fold variability, a more interpretable way to state "narrow performance range" than the p-value alone.
- Effect size: Mean paired difference 0.0036 (floor=0.05 minus floor=0.15) — small and in the *opposite* direction from what would justify keeping floor=0.15 on performance grounds (floor=0.05 is marginally higher); the manuscript's justification for 0.15 must rest on interpretability, not on this number, per the design decision already recorded in `CLAUDE.md`.
- Confidence interval: Not computed (same small-n limitation as Figure 4).
- Multiple-testing correction: N/A — single pre-specified comparison.
- Significance threshold: α=0.05 (not met — consistent with, and supporting, the "near-tie" framing).
- Pipeline step producing these statistics: Step 11 (`scripts/09_calibration.py` for floor=0.15; `scripts/_floor_sensitivity_check.py` for floor=0.05/0.0); the Wilcoxon test itself is now a permanent, saved computation in `scripts/12_figure5_calibration_validation.py` (Step 12, 2026-08-05), superseding the earlier ad hoc version.

**Reproducibility**:
- Pipeline stage: Step 11
- Generation script: `scripts/09_calibration.py` (production) + `scripts/_floor_sensitivity_check.py` (sensitivity check)
- Input files: `data/calibration/calibration_weights_and_validation.json`; floor=0.05/0.0 numbers currently recorded only in `CLAUDE.md` prose (see gap noted above)
- Manifest version: 2026-07-30
- Model version: `Human-GEM.json` version `1`
- Calibration version: production `(0.15, 0.15, 0.70)`; floor-sensitivity variants `(~0.5-0.6, 0.05, ~0.35-0.45)` at floor=0.05/0.0
- Output directory: `data/calibration/`
- Verification status: Production weights and LOCO/LOBO means verified against the JSON directly. Floor=0.05/0.0 numbers verified against this session's own captured tool output (not a persisted file — flagged as a reproducibility gap above). The corrected-miscitation history (an earlier, wrong "0.676 vs 0.883" framing) is documented in full in `CLAUDE.md` and should NOT be the source cited in the manuscript — only the reconciled table there is authoritative. The primary Wilcoxon p-value has also now been independently re-verified by a second, permanent method (`scripts/12_figure5_calibration_validation.py`, Step 12, 2026-08-05): W=8.0 and every effect-size number reproduced exactly, but p was found to be 0.1953, not the previously stated 0.1797 — corrected above.

**Caption outline**: "Nested cross-validation for empirical weight calibration. (A) LOCO/LOBO nested structure. (B) Performance is a near-tie across weight-floor settings (Wilcoxon p=0.195, paired by fold), differing by roughly 6-20x less than ordinary fold-to-fold variability. (C) Production weights (0.15/0.15/0.70) were retained for decomposability, not because they are the numerical optimum — a marginally higher-scoring, consensus-dominant alternative exists and is disclosed."

**Dependencies**: Step 11 complete; floor=0.05/0.0 results should be persisted to a structured file (not just `CLAUDE.md` prose) before this figure is built, per the gap flagged above.

---

## SUPPLEMENTARY FIGURE S1 — Escher Pathway Confidence Maps

**Role in Manuscript**: Provide pathway-level biological interpretation.

**Scientific Claim**: High-confidence metabolic alterations can be visualized directly within
metabolic network context.

**Hypothesis**: Confidence scores map coherently onto biologically meaningful pathways.

**Evidence**: Representative pathway visualizations.

**Reviewer Concern Addressed**: "Can the confidence scores be interpreted biologically?"

**Section**: Supplementary.

### Verified support

**Required pipeline outputs**: `data/calibration/all_groups_calibrated.csv` (per-reaction scores, confirmed exists — the data this figure would color pathway maps by) and `docs/STEP11_BIOMARKER_MAPPING.md` (confirmed exists — the reaction-to-biomarker/pathway mapping needed to select "representative" pathways).

**Input datasets / input files**: Same calibrated-score outputs as Figures 3-5; the 5 biomarker pathways already mapped to specific reaction IDs in `docs/STEP11_BIOMARKER_MAPPING.md` are the natural candidates for "representative pathway visualizations" (e.g. urea cycle, serine/glycine one-carbon metabolism).

**Analysis scripts**: **None exist.** Confirmed by direct check: `import escher` succeeds in the `GEM1` conda environment (version 1.8.1 installed, matching `CLAUDE.md`'s documented environment setup), but a repo-wide search found zero files referencing Escher, pathway maps, or visualization of any kind. Step 12 (Visualization & confidence reporting) is genuinely not started — this is not an oversight in this document, it's the actual current state of the repository.

**Panel layout** (proposed): One or two Escher map panels (e.g., urea cycle; serine/glycine one-carbon metabolism), reactions colored/sized by `calibrated_confidence_score`, using Human-GEM's existing map definitions if available in the installed Escher package, or a custom map built from the relevant biomarker's mapped reaction IDs if not.

**Statistical analysis**: N/A — this is a visualization, not a statistical figure.
- Primary statistical test: N/A
- Secondary statistical test: N/A
- Effect size: N/A
- Confidence interval: N/A
- Multiple-testing correction: N/A
- Significance threshold: N/A
- Pipeline step producing these statistics: N/A (the underlying scores come from Step 11, but no new statistic is computed for this figure itself)

**Reproducibility**:
- Pipeline stage: Step 11 (data source); Step 12 (rendering — not started)
- Generation script: None exists yet.
- Input files: `data/calibration/all_groups_calibrated.csv`, `docs/STEP11_BIOMARKER_MAPPING.md`
- Manifest version: N/A (no figure manifest exists)
- Model version: `Human-GEM.json` version `1`
- Calibration version: `(0.15, 0.15, 0.70)`, 2026-07-30
- Output directory: N/A
- Verification status: **Escher (v1.8.1) confirmed installed and importable in the `GEM1` conda environment** — the tooling is available, but zero implementation work has been done. This figure is entirely unbuilt; nothing here should be read as "in progress."

**Caption outline**: "Escher pathway map(s) showing calibrated confidence scores overlaid on [urea cycle / one-carbon metabolism], illustrating that GEM1's confidence estimates map onto biologically coherent, literature-recognized pathways rather than scattered, uninterpretable reactions."

**Dependencies**: Step 11 complete (data ready); Step 12 (Escher rendering script) does not exist and would need to be written from scratch.

---

## SUPPLEMENTARY FIGURE S2 — Extended Case Study Results

**STATUS (2026-08-05): DEFERRED — pending biomarker curation, not yet scoped.** Panel B needs an
explicit, literature-curated "additional biomarkers beyond the locked-in 5" list from Mohammed (see
`docs/STEP11_BIOMARKER_MAPPING.md` for what that curation work looks like) before it can be built at
all — as of this date, no such list exists anywhere in this repository or in any prior session record,
and none has been requested or started. This is not "awaiting an imminent decision"; it is an
un-started, unscoped curation task with no timeline. Panel A (subsystem-level summary) does NOT
depend on this curation — it is derivable today from existing data via the join already described
below — but was not built alongside Figures 1–5/S1 in Step 12 since this figure was explicitly
excluded from that batch of work pending Panel B's prerequisite. Do not build any part of this figure,
including Panel A alone, without first confirming with Mohammed whether a curation timeline exists.

**Role in Manuscript**: Provide additional biological analyses supporting the main conclusions without interrupting the narrative flow.

**Scientific Claim**: The GEM1 framework consistently prioritizes biologically relevant metabolic changes across multiple analyses.

**Hypothesis**: Additional pathways and biomarkers reinforce the conclusions presented in the main manuscript.

**Evidence**: Extended analyses, additional pathway summaries, supplementary biomarker results.

**Reviewer Concern Addressed**: "Do the conclusions remain consistent when examining additional results?"

**Section**: Supplementary.

### Verified support

**Required pipeline outputs**: Same base per-reaction data as Figures 3–5 (`data/calibration/all_groups_calibrated.csv`), plus the richer per-reaction metadata (`reaction_name`, `subsystem` — 147 unique subsystems, confirmed by direct count) that lives in `data/confidence_scores/{cohort}__{group}__confidence.csv` but was **not** carried into Step 11's own output. Confirmed directly: `all_groups_calibrated.csv`'s columns are `cohort, group, reaction_id, primary_consensus_score, is_active, is_growth_support, flux_confidence_raw, perturbation_robustness_raw, rank_consensus, rank_flux_confidence, rank_perturbation, calibrated_confidence_score` — no `reaction_name` or `subsystem`. Any "extended pathway summary" beyond the 5 locked-in biomarkers would need a join back to the confidence CSVs on `reaction_id` to recover that metadata; this join is straightforward with existing data but has not been done.

**Input datasets / input files**: Same 3 GEO cohorts as all other results figures. Note: `CLAUDE.md`'s "Decisions locked in" section specifies exactly 5 official calibration biomarkers — there is no pre-existing, curated "additional biomarker" set beyond those 5. "Additional pathways and biomarkers" for this figure would require new literature curation (analogous to `docs/STEP11_BIOMARKER_MAPPING.md`'s existing work), not just a re-plot of already-mapped reactions — this is a real scope gap, not a data-availability one.

**Analysis scripts**: None exist. No script in this repository aggregates confidence/calibration scores by subsystem, and no script defines any biomarker/pathway set beyond the 5 already locked in for Step 11.

**Panel layout** (proposed): Panel A: extended subsystem-level summary (e.g., top-N subsystems by mean `calibrated_confidence_score` across all 10 groups, derived via the join described above) — explicitly framed as exploratory/supplementary, not held to the same calibration standard as the 5 locked biomarkers. Panel B: any additional literature-suggested NAFLD signatures Mohammed chooses to curate, shown with the same per-cohort/per-group ranking view as Figure 3 Panel B, for direct visual comparability.

**Statistical analysis**: Same methodological caution as Figure 3 — no naive cross-cohort test should be applied to raw `calibrated_confidence_score` distributions for the reason already established there (rank-normalization forces group means toward ≈0.5001 regardless of true consistency). If new biomarkers are curated for this figure, they should be run through the same nested LOCO/LOBO machinery already built in `scripts/09_calibration.py`, not a new ad hoc statistic.
- Primary statistical test: Not applicable yet — no extended analysis has been built
- Secondary statistical test: Not applicable yet
- Effect size: Not applicable yet
- Confidence interval: Not applicable yet
- Multiple-testing correction: Will matter once a subsystem-level summary exists (147 subsystems → real multiple-comparisons exposure, e.g. Benjamini-Hochberg, should be planned for up front rather than retrofitted)
- Significance threshold: Not yet decided
- Pipeline step producing these statistics: None yet — this entire figure is unbuilt

**Reproducibility**:
- Pipeline stage: Would build on Steps 9–11
- Generation script: None exists
- Input files: `data/calibration/all_groups_calibrated.csv`, `data/confidence_scores/*.csv` (for the subsystem join), `docs/STEP11_BIOMARKER_MAPPING.md`
- Manifest version: N/A
- Model version: `Human-GEM.json` version `1`
- Calibration version: `(0.15, 0.15, 0.70)`, 2026-07-30
- Output directory: N/A
- Verification status: Entirely unbuilt. Nothing about this figure has been implemented, computed, or verified — stated plainly, same as Supplementary S1.

**Caption outline**: "Extended supplementary analyses: [subsystem-level confidence summary / additional biomarker pathways], demonstrating that GEM1's prioritization of biologically relevant alterations extends beyond the 5 primary calibration biomarkers presented in the main text."

**Dependencies**: Steps 9–11 complete (data ready); no generation script exists. **Pending biomarker curation — not yet scoped, deferred (2026-08-05).** Panel B's "additional biomarkers" scope has not been requested from Mohammed, has no timeline, and does not exist in any form on disk or in prior session records — this is a deferred, unscoped task, not a decision blocking imminently on someone's desk. Re-check with Mohammed before resuming any work on this figure, including Panel A alone.

---

## Manuscript Claim Map

| Core manuscript claim | Supporting figure(s) |
|---|---|
| GEM1 is methodologically novel | Figures 1-2 |
| GEM1 produces reproducible confidence estimates across independent cohorts | Figures 3 and 5 |
| Full-isoform aggregation is empirically preferable to single-proxy aggregation | Figure 4 |
| Confidence weights are selected through an objective calibration process | Figure 5 |
| Confidence scores support biological interpretation | Supplementary Figure S1 |
| Biological conclusions remain consistent across extended analyses | Supplementary Figure S2 |

---

## Open Items — RESOLVED (2026-07-30)

All 3 items originally flagged here are now built and independently verified:

1. **Metadata join (was: no `reaction_name`/`subsystem` on Step 11's output).** `scripts/10_join_calibration_metadata.py` → `data/calibration/all_groups_calibrated_with_metadata.csv` (129,310 rows, does not modify `all_groups_calibrated.csv`). Verified: row count unchanged by the join; `subsystem` has zero nulls (real join problem would show here); `reaction_name`'s 57,140 nulls exactly match the source confidence CSVs' own null count (confirmed genuine — many Human-GEM reactions, mostly transport/pool reactions, have blank names in the model itself, not a join defect); 5 random spot-checked rows match the source exactly.

2. **Figure 3 biomarker-ranking table.** `scripts/11_biomarker_ranking_by_group.py` → `data/calibration/biomarker_ranking_by_group.csv` (50 rows = 5 biomarkers × 10 groups; max-aggregated `calibrated_confidence_score`, argmax reaction, percentile rank via the same `method='min'` convention as `09_calibration.py`'s `evaluation_percentile`, and ordinal rank). **Investigating the output surfaced a genuine finding, not a bug**: `bcaa`'s cross-group percentile has by far the highest variance of the 5 biomarkers (mean 0.676, std 0.366, min 0.000077) because all 9 of its mapped reactions are structurally inactive (`is_active=False`, `primary_consensus_score=0` for 8 of 9; the 9th at 0.5 but still below the strict-consensus threshold) in exactly 2 groups — `GSE126848/healthy` and `GSE126848/obese_no_NAFLD` — confirmed directly against `all_groups_calibrated.csv`. This is the same signal already visible in `bcaa`'s weak LOBO fold score (0.676, matches this table's mean almost exactly) — now traced to its concrete cause rather than left as an unexplained weak point.

3. **Figure 5 floor-sensitivity persistence.** `scripts/_floor_sensitivity_check.py` modified to also write `data/calibration/floor_sensitivity_results.json`. All LOCO/LOBO fold-level numbers reproduced exactly on re-run (bit-for-bit match to the values already in this document and in `CLAUDE.md`). **Also fixed a real bug found while doing this**: the script's "production weights at this floor" line previously called an un-parameterized helper that silently ignored the loop's floor value, always reporting `(0.15,0.15,0.7)` regardless of which floor was being tested. Fixed to correctly thread the floor through; the corrected output shows floor=0.05 and floor=0.0 both select `(0.5, 0.05, 0.45)` (score 0.8835) as their actual per-floor production weights — consistent with, and a more precise version of, the "~0.5-0.6, 0.05, ~0.35-0.45" range already described in `CLAUDE.md`.
