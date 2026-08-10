"""
GEM1 independent audit -- shared constants and helpers.

Everything in this module reads only *raw* source files (Human-GEM.json,
Step 5 troppo JSONs, GEO expression matrices, all_groups_calibrated_with_metadata.csv)
or performs generic, standard-library-grade computation (cobra/scipy solver calls,
pandas ranking). It never imports any of GEM1's own pipeline scripts
(scripts/04_consensus_scoring.py, 05_confidence_engine.py, 06_flux_analysis.py,
07_flux_sampling.py, 08_perturbation_testing.py, 09_calibration.py,
10_join_calibration_metadata.py, 11_biomarker_ranking_by_group.py) -- those are the
things being audited, not tools this audit is allowed to lean on.

All formulas reimplemented here (rank normalization, evaluation percentile,
flux_confidence_raw, perturbation_robustness_raw, combined_score) are restated from
the plain-English/algebraic definitions traced out of the pipeline scripts during
audit setup -- they are re-derived, not copy-pasted function bodies.
"""

import json
import os

import numpy as np
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUDIT_DIR = os.path.join(REPO_ROOT, "audit")
DETAILS_DIR = os.path.join(AUDIT_DIR, "details")
LOGS_DIR = os.path.join(AUDIT_DIR, "logs")
os.makedirs(DETAILS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

MODEL_PATH = os.path.join(REPO_ROOT, "models", "Human-GEM.json")
CONTEXT_MODELS_DIR = os.path.join(REPO_ROOT, "data", "context_specific_models")
GEO_DIR = os.path.join(REPO_ROOT, "data", "geo_cohorts")
CALIBRATION_DIR = os.path.join(REPO_ROOT, "data", "calibration")
GROWTH_SUPPORT_PATH = os.path.join(CONTEXT_MODELS_DIR, "minimal_growth_support.json")
CALIBRATED_METADATA_CSV = os.path.join(CALIBRATION_DIR, "all_groups_calibrated_with_metadata.csv")
CALIBRATION_WEIGHTS_JSON = os.path.join(CALIBRATION_DIR, "calibration_weights_and_validation.json")

BIOMASS_REACTION = "MAR13082"
CORE_ALGORITHMS = ["fastcore", "imat", "tinit"]

# GSE126848/healthy has no fastcore model (documented, legitimate FASTCORE failure);
# every other (cohort, group) has all 3 core algorithms.
ALL_GROUPS = [
    ("GSE89632", "healthy"), ("GSE89632", "steatosis"), ("GSE89632", "NASH"),
    ("GSE126848", "healthy"), ("GSE126848", "obese_no_NAFLD"),
    ("GSE126848", "steatosis"), ("GSE126848", "NASH"),
    ("GSE135251", "healthy"), ("GSE135251", "steatosis"), ("GSE135251", "NASH"),
]

# Transcribed verbatim from docs/STEP11_BIOMARKER_MAPPING.md (a specification
# document, not pipeline code) and cross-checked against the persisted
# calibration_weights_and_validation.json's biomarker_mapping block.
BIOMARKERS_FULL = {
    "serine_glycine_shmt": ["MAR03845", "MAR04792"],
    "urea_cycle": ["MAR03873", "MAR03809", "MAR03811", "MAR03813", "MAR03816", "MAR08426"],
    "bcaa": ["MAR03744", "MAR03747", "MAR03765", "MAR06923", "MAR03777", "MAR03778",
             "MAR06416", "MAR06419", "MAR06421"],
    "dnl_scd1_fasn": ["MAR02150", "MAR02182", "MAR00146"],
    "choline_pc": ["MAR00636", "MAR00638", "MAR00653", "MAR01603", "MAR01606"],
}


def load_human_gem_raw():
    with open(MODEL_PATH) as f:
        return json.load(f)


def load_growth_support_ids():
    with open(GROWTH_SUPPORT_PATH) as f:
        d = json.load(f)
    return d


def load_step5_active_set(cohort, group, algorithm):
    """Raw Step 5 troppo output -- explicitly an allowed audit input."""
    path = os.path.join(CONTEXT_MODELS_DIR, f"{cohort}__{group}__{algorithm}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        d = json.load(f)
    return set(d["active_reactions"])


def compute_consensus_from_raw(cohort, group):
    """
    Independently recompute strict/majority reconstruction consensus straight
    from the raw Step 5 JSONs (set intersection / vote-count over
    fastcore+imat+tinit), rather than trusting data/consensus_scores/*.json.
    Mirrors the *definition* documented in 04_consensus_scoring.py's docstring
    (strict = intersection of available core algorithms, majority = count>=2),
    not its code.
    """
    available = {}
    for algo in CORE_ALGORITHMS:
        s = load_step5_active_set(cohort, group, algo)
        if s is not None:
            available[algo] = s
    if len(available) < 2:
        raise ValueError(f"{cohort}/{group}: fewer than 2 core algorithms available.")
    sets = list(available.values())
    strict = set.intersection(*sets)
    from collections import Counter
    counter = Counter()
    for s in sets:
        counter.update(s)
    majority = {rid for rid, c in counter.items() if c >= 2}
    return {
        "algorithms_used": sorted(available.keys()),
        "strict_consensus": strict,
        "majority_consensus": majority,
    }


def load_calibrated_metadata():
    return pd.read_csv(CALIBRATED_METADATA_CSV)


def rank_pct_within_active(df, value_col, group_cols=("cohort", "group"), active_col="is_active"):
    """
    Reimplementation of the rank-normalization step: percentile rank
    (0-1, method='average'), computed only among is_active rows, within each
    (cohort, group). Returns a Series aligned to df.index, NaN outside is_active.
    """
    out = pd.Series(np.nan, index=df.index)
    for _, sub in df[df[active_col]].groupby(list(group_cols)):
        out.loc[sub.index] = sub[value_col].rank(pct=True, method="average")
    return out


def evaluation_percentile(df, score_col, group_cols=("cohort", "group")):
    """
    Reimplementation of the 'evaluation percentile': percentile rank
    (0-1, method='min') of score_col among ALL reactions in each
    (cohort, group) -- ties (e.g. all zero-scored inactive reactions) sink to
    the bottom via method='min'.
    """
    out = pd.Series(np.nan, index=df.index)
    for _, sub in df.groupby(list(group_cols)):
        out.loc[sub.index] = sub[score_col].rank(pct=True, method="min")
    return out


def biomarker_metric(df, score_col, biomarker_ids, groups, group_cols=("cohort", "group")):
    """
    Mean, over all (biomarker, group) cells, of the MAX evaluation-percentile
    among that biomarker's mapped reactions present in that group.
    Reimplements the documented aggregation rule (OR-like max-over-isoforms).
    """
    cells = []
    for cohort, group in groups:
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        for _bm_name, rxn_ids in biomarker_ids.items():
            vals = sub.loc[sub["reaction_id"].isin(rxn_ids), score_col]
            if len(vals) == 0:
                continue
            cells.append(vals.max())
    return float(np.mean(cells)) if cells else float("nan")
