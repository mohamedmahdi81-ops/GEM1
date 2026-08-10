"""
GEM1 - Step 11: Empirical weight calibration.

Combines Step 10's per-reaction consensus score, Step 7's flux-sampling
uncertainty, and Step 8's perturbation robustness into a single
`calibrated_confidence_score`, with weights fit empirically against the 5
literature-curated calibration biomarkers (see docs/STEP11_BIOMARKER_MAPPING.md).

Design approved by Mohammed (see conversation record, 2026-07-28) before any of
this was written:

1. COMBINATION FORMULA: within each group, rank-normalize each of the three
   raw signals to a percentile in (0,1] among that group's "active" reactions
   -- defined as strict_consensus_reactions UNION growth_support, i.e.
   exactly the set Step 6 built its constrained model from (NOT
   "primary_consensus_score > 0": majority/minority-consensus-tier reactions
   have nonzero consensus but were never part of the model Steps 6/7/8
   solved, so no flux/perturbation data exists for them) -- then take a
   weighted linear combination
       calibrated = w1*rank(consensus) + w2*rank(flux_confidence) + w3*rank(robustness)
   with w1+w2+w3=1 and each wi >= 0.15 (a floor constraint -- with only 5
   biomarkers to calibrate against, an unconstrained fit could collapse onto
   a single axis and quietly defeat the whole point of a 3-axis score).
   Reactions outside that active set get calibrated_confidence_score fixed
   at 0 -- mirrors Step 10's own "Layer 1 is foundational" rule, and
   sidesteps missing flux/perturbation data for reactions Steps 6-8 never
   built a model around in the first place.

2. WEIGHT SELECTION IS NESTED, NOT CIRCULAR: weights are never chosen by
   directly maximizing performance on the fold being reported as "held out".
   Leave-one-cohort-out (LOCO, 3 outer folds) and leave-one-biomarker-out
   (LOBO, 5 outer folds) each have their own INNER cross-validation (over the
   outer-training data only) to select weights; those fixed weights are then
   applied exactly once to the untouched outer-test fold. The reported
   LOCO/LOBO numbers are generalization ESTIMATES, not a source of the
   deployed weights. A separate, final grid search over ALL data (no holdout)
   produces the actual production weights used to populate
   calibrated_confidence_score for every reaction. Two different artifacts,
   not one -- see the "Two separate deliverables" point in the design
   conversation.

   STATED LIMITATION: LOCO's inner loop is thin (only 2 inner folds, forced
   by there being only 3 cohorts total) -- carried forward per Mohammed's
   explicit instruction, same treatment as the "5 biomarkers is small"
   caveat below.

3. ISOFORM AGGREGATION: docs/STEP11_BIOMARKER_MAPPING.md's biomarkers range
   from 1 mapped reaction (SCD1) to ~9 (BCAAs). Treating every mapped
   reaction as an independent vote in the calibration objective would let
   BCAAs/urea-cycle dominate purely by reaction count. Each biomarker is
   collapsed to ONE scalar per group via MAX over its mapped reactions (an
   OR-like "is there a confidently-active reaction in this pathway"
   semantics, consistent with how the GPRs themselves use OR logic for
   isoforms) before it touches any calibration objective. A secondary,
   single-reaction-proxy variant (one best-justified reaction per biomarker,
   collapsing isoform/compartment redundancy but NOT collapsing distinct
   pathway steps) is run alongside as a sensitivity check, not the primary
   result.

STATED LIMITATIONS (report honestly, do not hide):
 - Only 5 biomarkers and 10 groups to calibrate/validate against -- both
   LOCO and LOBO nested-CV estimates will have wide uncertainty.
 - LOCO's inner loop has only 2 folds (3 cohorts total, one held out as the
   outer test leaves only 2 for inner selection).
 - The single-reaction-proxy set includes two judgment calls not explicitly
   pinned down in docs/STEP11_BIOMARKER_MAPPING.md (which BCAT isoform, which
   single FASN step) -- flagged in BIOMARKERS_SINGLE_PROXY below, not
   silently decided.

Run from gem1-main (pandas/numpy/scipy only -- no solver calls, this step is
pure post-hoc analysis of Steps 6-10's already-computed outputs).
"""

import os
import json
import itertools
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CONFIDENCE_DIR = os.path.join(BASE_DIR, "data", "confidence_scores")
FLUX_SAMPLING_DIR = os.path.join(BASE_DIR, "data", "flux_sampling")
FLUX_ANALYSIS_DIR = os.path.join(BASE_DIR, "data", "flux_analysis")
PERTURBATION_DIR = os.path.join(BASE_DIR, "data", "perturbation_testing")
CONTEXT_MODELS_DIR = os.path.join(BASE_DIR, "data", "context_specific_models")
OUT_DIR = os.path.join(BASE_DIR, "data", "calibration")
os.makedirs(OUT_DIR, exist_ok=True)

GROUPS = {
    "GSE89632": ["healthy", "steatosis", "NASH"],
    "GSE126848": ["healthy", "obese_no_NAFLD", "steatosis", "NASH"],
    "GSE135251": ["healthy", "steatosis", "NASH"],
}
ALL_GROUPS = [(cohort, group) for cohort, groups in GROUPS.items() for group in groups]

WEIGHT_FLOOR = 0.15
GRID_STEP = 0.05

# ---------------------------------------------------------------------------
# Biomarker -> reaction-ID mapping, transcribed from
# docs/STEP11_BIOMARKER_MAPPING.md (final, resolved version, 2026-07-21).
# Full isoform/compartment set = primary. Single-reaction proxy = secondary
# sensitivity check.
# ---------------------------------------------------------------------------

BIOMARKERS_FULL = {
    "serine_glycine_shmt": ["MAR03845", "MAR04792"],
    "urea_cycle": ["MAR03873", "MAR03809", "MAR03811", "MAR03813", "MAR03816", "MAR08426"],
    "bcaa": ["MAR03744", "MAR03747", "MAR03765", "MAR06923", "MAR03777", "MAR03778",
             "MAR06416", "MAR06419", "MAR06421"],
    "dnl_scd1_fasn": ["MAR02150", "MAR02182", "MAR00146"],
    "choline_pc": ["MAR00636", "MAR00638", "MAR00653", "MAR01603", "MAR01606"],
}

# Single-reaction-proxy set. Doc explicitly names the proxy for SHMT (SHMT2),
# ARG1 vs ARG2, and choline (CCT alone). Two choices are NOT explicitly
# pinned down in the doc and are this script's own reasonable-default
# judgment calls, flagged here rather than decided silently:
#   - BCAT isoform: doc states no preference between BCAT1/BCAT2 for a
#     single-proxy case. Defaulting to BCAT2 (mitochondrial) since it sits in
#     the same compartment as the downstream BCKDH step it feeds, forming one
#     coherent mitochondrial sub-pathway -- not a doc-stated preference.
#   - FASN single step: doc doesn't say which of the two FASN reactions
#     (MAR02150 first-step vs MAR02182 terminal/palmitate-releasing step) to
#     keep alone. Defaulting to MAR02182 (terminal, product-defining step),
#     matching the same "rate-limiting/terminal step is the single proxy"
#     pattern the doc already uses for choline (CCT) and urea cycle (ARG1).
BIOMARKERS_SINGLE_PROXY = {
    "serine_glycine_shmt": ["MAR04792"],  # SHMT2 -- doc-stated preference
    "urea_cycle": ["MAR03873", "MAR03809", "MAR03811", "MAR03813", "MAR03816"],  # ARG1 -- doc-stated
    "bcaa": ["MAR03765", "MAR03744", "MAR03778", "MAR06416", "MAR06419", "MAR06421"],  # BCAT2 -- OUR default, not doc-stated
    "dnl_scd1_fasn": ["MAR02182", "MAR00146"],  # FASN terminal step only -- OUR default, not doc-stated
    "choline_pc": ["MAR00638"],  # CCT alone -- doc-stated preference
}

BIOMARKER_SETS = {
    "full_isoform": BIOMARKERS_FULL,
    "single_proxy": BIOMARKERS_SINGLE_PROXY,
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_growth_support():
    with open(os.path.join(CONTEXT_MODELS_DIR, "minimal_growth_support.json")) as f:
        return set(json.load(f)["support_reaction_ids"])


def load_strict_consensus(cohort, group):
    with open(os.path.join(BASE_DIR, "data", "consensus_scores", f"{cohort}__{group}__consensus.json")) as f:
        return set(json.load(f)["strict_consensus_reactions"])


def load_group_raw_scores(cohort, group, growth_support):
    """Assemble one row per reaction (12,931 total) with the three raw
    signals, for a single (cohort, group).

    "Active" here means "in Step 6's constrained model" -- i.e.
    strict_consensus_reactions UNION growth_support -- NOT
    "primary_consensus_score > 0". Majority/minority-consensus-tier
    reactions (1 or 2 of 3 algorithms active, but not all 3) have a nonzero
    primary_consensus_score but were never part of the model Steps 6/7/8
    built and solved, so no flux-sampling or perturbation data exists for
    them by design. Those reactions get calibrated_confidence_score fixed at
    0, same as true no_support reactions -- this pipeline's flux/robustness
    evidence is scoped to the strict-consensus model, not the looser
    consensus tiers, and Step 11 respects that same scope rather than
    inventing data Steps 6-8 never computed."""
    conf = pd.read_csv(os.path.join(CONFIDENCE_DIR, f"{cohort}__{group}__confidence.csv"))
    conf = conf[["reaction_id", "primary_consensus_score"]].copy()

    strict_consensus = load_strict_consensus(cohort, group)
    active_set = strict_consensus | growth_support

    flux_summary = pd.read_csv(os.path.join(FLUX_SAMPLING_DIR, f"{cohort}__{group}__summary.csv"))
    flux_summary = flux_summary[["reaction_id", "p5", "p95"]].set_index("reaction_id")

    with open(os.path.join(FLUX_ANALYSIS_DIR, f"{cohort}__{group}__flux.json")) as f:
        fva = json.load(f)["fva"]
    fva_df = pd.DataFrame.from_dict(fva, orient="index")  # columns: minimum, maximum
    fva_df.index.name = "reaction_id"

    perturbation = pd.read_csv(os.path.join(PERTURBATION_DIR, f"{cohort}__{group}__reactions.csv"))
    perturbation = perturbation[["reaction_id", "perturbation_robustness_score"]].set_index("reaction_id")

    df = conf.set_index("reaction_id")
    df["is_active"] = df.index.isin(active_set)
    df["is_growth_support"] = df.index.isin(growth_support)

    # Flux confidence: 1 - (p95-p5)/(fva_max-fva_min), clipped [0,1].
    # fva_max == fva_min (fully pinned reaction) -> confidence 1.0.
    joined = df.join(flux_summary, how="left").join(fva_df, how="left")
    span = joined["maximum"] - joined["minimum"]
    sample_spread = joined["p95"] - joined["p5"]
    flux_confidence = np.where(
        span > 0,
        1.0 - (sample_spread / span),
        1.0,
    )
    flux_confidence = np.clip(flux_confidence, 0.0, 1.0)
    joined["flux_confidence_raw"] = np.where(joined["is_active"], flux_confidence, np.nan)

    # Perturbation robustness, rescaled to [0,1]. Growth-support reactions
    # never appear in the Step 8 candidate set (never eligible for removal)
    # -> fixed at 1.0. Non-growth-support active reactions get the Step 8
    # score, rescaled. Inactive reactions -> NaN (never combined).
    joined = joined.join(perturbation, how="left")
    robustness_rescaled = (joined["perturbation_robustness_score"] + 1.0) / 2.0
    robustness_final = np.where(
        joined["is_growth_support"], 1.0,
        np.where(joined["is_active"], robustness_rescaled, np.nan),
    )
    # Sanity check: every active, non-growth-support reaction must have had
    # a real Step 8 row (candidates == strict_consensus - growth_support
    # should cover exactly this set). Surface a loud failure if not, rather
    # than silently defaulting.
    unexpected_missing = joined["is_active"] & ~joined["is_growth_support"] & joined["perturbation_robustness_score"].isna()
    if unexpected_missing.any():
        raise ValueError(
            f"{cohort}/{group}: {unexpected_missing.sum()} active, non-growth-support "
            f"reactions have no Step 8 perturbation score -- expected exact coverage "
            f"(candidates = strict_consensus - growth_support). Investigate before trusting "
            f"calibration output for this group."
        )
    joined["perturbation_robustness_raw"] = robustness_final

    joined["cohort"] = cohort
    joined["group"] = group
    return joined.reset_index()[[
        "cohort", "group", "reaction_id", "primary_consensus_score", "is_active",
        "is_growth_support", "flux_confidence_raw", "perturbation_robustness_raw",
    ]]


def add_rank_normalization(df):
    """Percentile-rank each of the 3 raw signals, computed ONLY among active
    (is_active -- strict_consensus UNION growth_support) reactions within
    each group -- this is the population the combination formula actually
    operates on. Standard
    average-rank tie-breaking (fine here; no huge zero-inflated tie block
    since we've already restricted to the active subset)."""
    df = df.copy()
    for col, rank_col in [
        ("primary_consensus_score", "rank_consensus"),
        ("flux_confidence_raw", "rank_flux_confidence"),
        ("perturbation_robustness_raw", "rank_perturbation"),
    ]:
        df[rank_col] = np.nan
        active_mask = df["is_active"]
        df.loc[active_mask, rank_col] = (
            df.loc[active_mask]
            .groupby(["cohort", "group"])[col]
            .rank(pct=True, method="average")
        )
    return df


def combined_score(df, weights):
    w1, w2, w3 = weights
    combined = (
        w1 * df["rank_consensus"].fillna(0)
        + w2 * df["rank_flux_confidence"].fillna(0)
        + w3 * df["rank_perturbation"].fillna(0)
    )
    # Reactions outside the active set (is_active == False) are fixed at 0, regardless
    # of the (NaN -> 0-filled) rank terms above -- this line is the actual
    # authority, the fillna(0) above is just to keep the arithmetic defined.
    combined = np.where(df["is_active"], combined, 0.0)
    return combined


def evaluation_percentile(df, weights):
    """Percentile rank of calibrated_confidence_score among ALL reactions in
    each group (not just the active subset) -- inactive reactions (score
    fixed at 0) must land at the bottom of the distribution, not in the
    middle of a huge tied-at-zero block, hence method='min' for ties (all
    tied reactions get the lowest rank in their tie group)."""
    df = df.copy()
    df["_combined"] = combined_score(df, weights)
    df["_pct"] = (
        df.groupby(["cohort", "group"])["_combined"]
        .rank(pct=True, method="min")
    )
    return df


def biomarker_metric(df_with_pct, biomarker_mapping, group_subset=None, biomarker_subset=None):
    """Mean, over (biomarker, group) cells in the requested subset, of the
    MAX percentile among that biomarker's mapped reactions present in that
    group. This is THE evaluation metric used at every level (inner
    selection criterion and outer generalization estimate alike) -- only the
    subset of (biomarker, group) cells averaged over differs between calls."""
    groups = group_subset if group_subset is not None else ALL_GROUPS
    biomarkers = biomarker_subset if biomarker_subset is not None else list(biomarker_mapping.keys())

    cells = []
    for cohort, group in groups:
        sub = df_with_pct[(df_with_pct["cohort"] == cohort) & (df_with_pct["group"] == group)]
        sub = sub.set_index("reaction_id")["_pct"]
        for name in biomarkers:
            rxn_ids = biomarker_mapping[name]
            present = sub.reindex(rxn_ids).dropna()
            if len(present) == 0:
                continue  # biomarker has zero mapped reactions in this network slice; skip, don't zero-fill silently
            cells.append(present.max())
    if not cells:
        return np.nan
    return float(np.mean(cells))


def generate_weight_grid(step=GRID_STEP, floor=WEIGHT_FLOOR):
    grid = []
    w1 = floor
    while w1 <= 1 - 2 * floor + 1e-9:
        w2 = floor
        while w2 <= 1 - w1 - floor + 1e-9:
            w3 = round(1 - w1 - w2, 10)
            if w3 >= floor - 1e-9:
                grid.append((round(w1, 4), round(w2, 4), round(w3, 4)))
            w2 = round(w2 + step, 10)
        w1 = round(w1 + step, 10)
    return grid


def select_best_weights(df_with_ranks, biomarker_mapping, group_subset, biomarker_subset=None):
    """Grid search: pick the weight vector maximizing biomarker_metric over
    the given (group_subset, biomarker_subset). This is called ONLY on
    training/inner data -- never on the fold being reported as held out."""
    grid = generate_weight_grid()
    best_weights, best_score = None, -np.inf
    for w in grid:
        df_pct = evaluation_percentile(df_with_ranks, w)
        score = biomarker_metric(df_pct, biomarker_mapping, group_subset, biomarker_subset)
        if score > best_score:
            best_score, best_weights = score, w
    return best_weights, best_score


# ---------------------------------------------------------------------------
# Nested cross-validation
# ---------------------------------------------------------------------------

def nested_loco(df_with_ranks, biomarker_mapping):
    """Leave-one-cohort-out, nested. Outer: 3 folds (one cohort held out).
    Inner: leave-one-cohort-out over the remaining 2 cohorts (2 inner folds
    -- thin, forced by only 3 cohorts existing; stated limitation, not
    hidden)."""
    cohorts = list(GROUPS.keys())
    outer_results = []
    for outer_test_cohort in cohorts:
        train_cohorts = [c for c in cohorts if c != outer_test_cohort]
        outer_test_groups = [(c, g) for c, g in ALL_GROUPS if c == outer_test_cohort]

        inner_scores_per_candidate = {}
        grid = generate_weight_grid()
        for w in grid:
            inner_vals = []
            for inner_val_cohort in train_cohorts:
                inner_train_cohorts = [c for c in train_cohorts if c != inner_val_cohort]
                inner_val_groups = [(c, g) for c, g in ALL_GROUPS if c == inner_val_cohort]
                df_pct = evaluation_percentile(df_with_ranks, w)
                score = biomarker_metric(df_pct, biomarker_mapping, inner_val_groups)
                inner_vals.append(score)
            inner_scores_per_candidate[w] = float(np.mean(inner_vals))
        best_inner_weights = max(inner_scores_per_candidate, key=inner_scores_per_candidate.get)
        best_inner_score = inner_scores_per_candidate[best_inner_weights]

        df_pct_outer = evaluation_percentile(df_with_ranks, best_inner_weights)
        outer_score = biomarker_metric(df_pct_outer, biomarker_mapping, outer_test_groups)

        outer_results.append({
            "outer_test_cohort": outer_test_cohort,
            "selected_weights": best_inner_weights,
            "inner_selection_score": best_inner_score,
            "outer_held_out_score": outer_score,
        })

    mean_outer = float(np.mean([r["outer_held_out_score"] for r in outer_results]))
    return {"folds": outer_results, "mean_outer_held_out_score": mean_outer,
            "note": "Inner loop has only 2 folds (3 cohorts total) -- thin, stated limitation."}


def nested_lobo(df_with_ranks, biomarker_mapping):
    """Leave-one-biomarker-out, nested. Outer: 5 folds (one biomarker held
    out). Inner: leave-one-biomarker-out over the remaining 4 biomarkers (4
    inner folds)."""
    biomarker_names = list(biomarker_mapping.keys())
    outer_results = []
    for outer_test_bio in biomarker_names:
        train_bios = [b for b in biomarker_names if b != outer_test_bio]

        inner_scores_per_candidate = {}
        grid = generate_weight_grid()
        for w in grid:
            inner_vals = []
            for inner_val_bio in train_bios:
                inner_train_bios = [b for b in train_bios if b != inner_val_bio]
                df_pct = evaluation_percentile(df_with_ranks, w)
                score = biomarker_metric(df_pct, biomarker_mapping, ALL_GROUPS, [inner_val_bio])
                inner_vals.append(score)
            inner_scores_per_candidate[w] = float(np.mean(inner_vals))
        best_inner_weights = max(inner_scores_per_candidate, key=inner_scores_per_candidate.get)
        best_inner_score = inner_scores_per_candidate[best_inner_weights]

        df_pct_outer = evaluation_percentile(df_with_ranks, best_inner_weights)
        outer_score = biomarker_metric(df_pct_outer, biomarker_mapping, ALL_GROUPS, [outer_test_bio])

        outer_results.append({
            "outer_test_biomarker": outer_test_bio,
            "selected_weights": best_inner_weights,
            "inner_selection_score": best_inner_score,
            "outer_held_out_score": outer_score,
        })

    mean_outer = float(np.mean([r["outer_held_out_score"] for r in outer_results]))
    return {"folds": outer_results, "mean_outer_held_out_score": mean_outer}


def fit_production_weights(df_with_ranks, biomarker_mapping):
    """Final weights for actual deployment: same grid-search selection
    procedure, but run once over ALL data (no outer holdout at all). A
    separate artifact from the nested-CV generalization estimates above --
    see module docstring."""
    weights, score = select_best_weights(df_with_ranks, biomarker_mapping, ALL_GROUPS)
    return {"weights": weights, "selection_score": score}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading growth-support set and assembling per-group raw score tables...")
    growth_support = load_growth_support()

    all_rows = []
    for cohort, group in ALL_GROUPS:
        rows = load_group_raw_scores(cohort, group, growth_support)
        all_rows.append(rows)
        n_active = int(rows["is_active"].sum())
        print(f"  {cohort}/{group}: {len(rows)} reactions, {n_active} active")
    df = pd.concat(all_rows, ignore_index=True)
    df = add_rank_normalization(df)

    results_by_biomarker_set = {}
    for set_name, biomarker_mapping in BIOMARKER_SETS.items():
        print(f"\n=== Biomarker set: {set_name} ({len(biomarker_mapping)} biomarkers) ===")

        print("Running nested LOCO (leave-one-cohort-out)...")
        loco = nested_loco(df, biomarker_mapping)
        print(f"  mean outer held-out score: {loco['mean_outer_held_out_score']:.4f}")
        for fold in loco["folds"]:
            print(f"    outer_test={fold['outer_test_cohort']:10} "
                  f"selected_weights={fold['selected_weights']} "
                  f"inner_score={fold['inner_selection_score']:.4f} "
                  f"outer_score={fold['outer_held_out_score']:.4f}")

        print("Running nested LOBO (leave-one-biomarker-out)...")
        lobo = nested_lobo(df, biomarker_mapping)
        print(f"  mean outer held-out score: {lobo['mean_outer_held_out_score']:.4f}")
        for fold in lobo["folds"]:
            print(f"    outer_test={fold['outer_test_biomarker']:24} "
                  f"selected_weights={fold['selected_weights']} "
                  f"inner_score={fold['inner_selection_score']:.4f} "
                  f"outer_score={fold['outer_held_out_score']:.4f}")

        print("Fitting final production weights (all data, no holdout)...")
        production = fit_production_weights(df, biomarker_mapping)
        print(f"  production weights: {production['weights']} "
              f"(selection score {production['selection_score']:.4f})")

        results_by_biomarker_set[set_name] = {
            "biomarker_mapping": biomarker_mapping,
            "loco": loco,
            "lobo": lobo,
            "production": production,
        }

    # Apply the PRIMARY (full_isoform) production weights to populate the
    # actual per-reaction calibrated_confidence_score used downstream.
    primary_weights = results_by_biomarker_set["full_isoform"]["production"]["weights"]
    df["calibrated_confidence_score"] = combined_score(df, primary_weights)

    for cohort, group in ALL_GROUPS:
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        out_path = os.path.join(OUT_DIR, f"{cohort}__{group}__calibrated.csv")
        sub.to_csv(out_path, index=False)

    combined_path = os.path.join(OUT_DIR, "all_groups_calibrated.csv")
    df.to_csv(combined_path, index=False)

    manifest_rows = []
    for cohort, group in ALL_GROUPS:
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        manifest_rows.append({
            "cohort": cohort,
            "group": group,
            "n_reactions_total": len(sub),
            "n_active": int(sub["is_active"].sum()),
            "mean_calibrated_confidence_score": round(sub["calibrated_confidence_score"].mean(), 4),
            "mean_calibrated_confidence_score_active_only": round(
                sub.loc[sub["is_active"], "calibrated_confidence_score"].mean(), 4),
        })
    manifest_path = os.path.join(OUT_DIR, "calibration_manifest.csv")
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    weights_path = os.path.join(OUT_DIR, "calibration_weights_and_validation.json")
    with open(weights_path, "w") as f:
        json.dump({
            "production_weights_used": {"w1_consensus": primary_weights[0],
                                         "w2_flux_confidence": primary_weights[1],
                                         "w3_perturbation": primary_weights[2]},
            "results_by_biomarker_set": results_by_biomarker_set,
            "weight_floor": WEIGHT_FLOOR,
            "grid_step": GRID_STEP,
            "limitations": [
                "Only 5 biomarkers and 10 groups to calibrate/validate against -- "
                "both LOCO and LOBO nested-CV estimates have wide uncertainty.",
                "LOCO's inner loop has only 2 folds (3 cohorts total, one held out "
                "as the outer test leaves only 2 for inner selection).",
                "single_proxy biomarker set includes two judgment calls not "
                "explicitly pinned down in docs/STEP11_BIOMARKER_MAPPING.md "
                "(BCAT isoform choice, single FASN step choice) -- see "
                "BIOMARKERS_SINGLE_PROXY comments in this script.",
            ],
        }, f, indent=2, default=str)

    print(f"\nPer-group calibrated CSVs + combined table: {combined_path}")
    print(f"Manifest: {manifest_path}")
    print(f"Weights + full nested-CV validation record: {weights_path}")
    print(f"\nProduction weights (full_isoform set, used to populate calibrated_confidence_score): {primary_weights}")


if __name__ == "__main__":
    main()
