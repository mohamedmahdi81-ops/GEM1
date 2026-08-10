"""
GEM1 - Step 10: Hierarchical Confidence Engine (v1 -- reconstruction-consensus
layer only).

This is GEM1's core methodological contribution: a decomposable, per-reaction
confidence score with a fully auditable decision trace, built in explicit
layers so additional evidence (flux-sampling uncertainty, perturbation
robustness) can be added later without changing the output schema. Full
formula definitions and design rationale are documented in CLAUDE.md under
"Step 10 methodology" -- this docstring gives the short version; CLAUDE.md is
the source of truth if the two ever disagree.

Layer 0 (raw evidence): per-algorithm active/inactive calls for FASTCORE,
iMAT, tINIT (the "core" algorithms) and GIMME, taken directly from Step 5's
JSONs in data/context_specific_models/ -- no re-derivation, no solver calls.

Layer 1 (primary reconstruction-consensus score): built ONLY from FASTCORE +
iMAT + tINIT (GSE126848/healthy uses iMAT + tINIT only, since FASTCORE has no
model for that group -- see Step 9). This is the entire primary confidence
score: primary_consensus_score = n_core_algorithms_active /
n_core_algorithms_available, in {0, 1/3, 2/3, 1} for 3-algorithm groups or
{0, 1/2, 1} for the one 2-algorithm group. final_confidence_score is set
identically equal to primary_consensus_score -- by construction, nothing in
Layer 2 can change it.

Layer 2 (complementary validation, GIMME): reported as a separate,
non-score-affecting annotation (gimme_active, gimme_corroboration_status).
GIMME is far more inclusive than the three core algorithms (80-96% of the
network active vs. 39-64%, see Step 9's overlap analysis), so treating its
"active" calls as equally strong evidence would be indefensible; it is used
here only to flag whether the primary consensus is corroborated by an
independent, more permissive method, never to raise a reaction's score.

Layer 3 (reserved, not yet computed): flux_sampling_uncertainty_score (Step
7, pending), perturbation_robustness_score (Step 8, pending),
calibrated_confidence_score (Step 11, pending -- the empirically-calibrated
combination of layers 1/3). Columns exist now, populated with None, so this
script's output schema does not need to change when those steps land.

Reproducibility: pure function of Step 5 (data/context_specific_models/*.json)
and Step 9 (data/consensus_scores/consensus_manifest.csv) outputs. No RNG, no
solver calls, deterministic given those inputs -- verified by cross-checking
every group's strict/majority consensus counts against Step 9's own manifest
before writing any output (see assert in main()).

Run from gem1-main (needs cobra for reaction name/subsystem metadata only;
Human-GEM is loaded, never re-solved).
"""

import os
import json
import pandas as pd
import cobra

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "models", "Human-GEM.json")
CONTEXT_MODELS_DIR = os.path.join(BASE_DIR, "data", "context_specific_models")
CONSENSUS_DIR = os.path.join(BASE_DIR, "data", "consensus_scores")
OUT_DIR = os.path.join(BASE_DIR, "data", "confidence_scores")
os.makedirs(OUT_DIR, exist_ok=True)

CORE_ALGORITHMS = ["fastcore", "imat", "tinit"]


def load_active_set(cohort, group, algorithm):
    path = os.path.join(CONTEXT_MODELS_DIR, f"{cohort}__{group}__{algorithm}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return set(json.load(f)["active_reactions"])


def build_reaction_metadata(model):
    return {r.id: (r.name, r.subsystem) for r in model.reactions}


def consensus_tier(n_active, n_available):
    if n_active == n_available:
        return "strict_consensus"
    elif n_active >= 2:
        return "majority_consensus"
    elif n_active == 1:
        return "minority_signal"
    else:
        return "no_support"


def gimme_corroboration(primary_score, gimme_active):
    if gimme_active is None:
        return "gimme_unavailable"
    if primary_score > 0 and gimme_active:
        return "corroborated"
    if primary_score == 0 and gimme_active:
        return "gimme_only"
    if primary_score > 0 and not gimme_active:
        return "unsupported_by_gimme"
    return "absent_from_all"


def score_group(cohort, group, reaction_meta):
    core_sets = {algo: load_active_set(cohort, group, algo) for algo in CORE_ALGORITHMS}
    available = {algo: s for algo, s in core_sets.items() if s is not None}
    missing = [algo for algo, s in core_sets.items() if s is None]
    n_available = len(available)
    if n_available < 2:
        raise ValueError(f"{cohort}/{group}: only {n_available} core algorithms "
                          f"available, need >=2 to score.")
    consensus_def = "strict" if n_available == 3 else "reduced"

    gimme_set = load_active_set(cohort, group, "gimme")

    rows = []
    for r_id, (name, subsystem) in reaction_meta.items():
        votes = {algo: (r_id in s) for algo, s in available.items()}
        n_active = sum(votes.values())
        primary_score = n_active / n_available
        tier = consensus_tier(n_active, n_available)

        gimme_active = (r_id in gimme_set) if gimme_set is not None else None
        corroboration = gimme_corroboration(primary_score, gimme_active)

        rows.append({
            "cohort": cohort,
            "group": group,
            "reaction_id": r_id,
            "reaction_name": name,
            "subsystem": subsystem,
            "fastcore_active": votes.get("fastcore"),
            "imat_active": votes.get("imat"),
            "tinit_active": votes.get("tinit"),
            "n_core_algorithms_available": n_available,
            "n_core_algorithms_active": n_active,
            "consensus_definition_used": consensus_def,
            "algorithms_missing": ",".join(missing) if missing else "",
            "consensus_tier": tier,
            "primary_consensus_score": primary_score,
            "gimme_active": gimme_active,
            "gimme_corroboration_status": corroboration,
            # By construction, identical to primary_consensus_score -- GIMME
            # (Layer 2) never modifies it. Kept as its own column so the
            # invariant is visible in the output itself, not just in code.
            "final_confidence_score": primary_score,
            "flux_sampling_uncertainty_score": None,  # reserved -- Step 7, pending
            "perturbation_robustness_score": None,     # reserved -- Step 8, pending
            "calibrated_confidence_score": None,        # reserved -- Step 11, pending
        })

    return pd.DataFrame(rows), consensus_def


def main():
    print("Loading Human-GEM for reaction metadata (name, subsystem) only "
          "-- not re-solving anything...")
    model = cobra.io.load_json_model(MODEL_PATH)
    reaction_meta = build_reaction_metadata(model)
    print(f"  {len(reaction_meta)} reactions loaded.")

    consensus_manifest = pd.read_csv(os.path.join(CONSENSUS_DIR, "consensus_manifest.csv"))

    all_group_frames = []
    summary_rows = []
    for _, row in consensus_manifest.iterrows():
        cohort, group = row["cohort"], row["group"]
        df, consensus_def = score_group(cohort, group, reaction_meta)

        # Correctness check: this script's own strict/majority tier counts
        # must exactly match Step 9's already-published consensus counts,
        # since both are derived from the same underlying JSONs. A mismatch
        # here means a bug in this script, not a data discrepancy.
        n_strict_here = int((df["consensus_tier"] == "strict_consensus").sum())
        n_majority_here = int(df["consensus_tier"].isin(
            ["strict_consensus", "majority_consensus"]).sum())
        assert n_strict_here == row["n_strict_consensus"], (
            f"{cohort}/{group}: strict consensus mismatch vs Step 9 manifest "
            f"({n_strict_here} != {row['n_strict_consensus']})"
        )
        assert n_majority_here == row["n_majority_consensus"], (
            f"{cohort}/{group}: majority consensus mismatch vs Step 9 manifest "
            f"({n_majority_here} != {row['n_majority_consensus']})"
        )
        assert consensus_def == row["consensus_definition_used"], (
            f"{cohort}/{group}: consensus_definition_used mismatch "
            f"({consensus_def} != {row['consensus_definition_used']})"
        )

        out_path = os.path.join(OUT_DIR, f"{cohort}__{group}__confidence.csv")
        df.to_csv(out_path, index=False)
        all_group_frames.append(df)

        tier_counts = df["consensus_tier"].value_counts().to_dict()
        corrob_counts = df["gimme_corroboration_status"].value_counts().to_dict()
        summary_rows.append({
            "cohort": cohort,
            "group": group,
            "consensus_definition_used": consensus_def,
            "n_reactions_total": len(df),
            "n_strict_consensus": tier_counts.get("strict_consensus", 0),
            "n_majority_consensus": tier_counts.get("majority_consensus", 0),
            "n_minority_signal": tier_counts.get("minority_signal", 0),
            "n_no_support": tier_counts.get("no_support", 0),
            "n_gimme_corroborated": corrob_counts.get("corroborated", 0),
            "n_gimme_only": corrob_counts.get("gimme_only", 0),
            "n_unsupported_by_gimme": corrob_counts.get("unsupported_by_gimme", 0),
            "n_absent_from_all": corrob_counts.get("absent_from_all", 0),
            "mean_final_confidence_score": round(df["final_confidence_score"].mean(), 4),
        })
        print(f"{cohort}/{group}: verified vs Step 9 (strict={n_strict_here}, "
              f"majority={n_majority_here}, def={consensus_def}) -> {out_path}")

    combined = pd.concat(all_group_frames, ignore_index=True)
    combined_path = os.path.join(OUT_DIR, "all_groups_confidence.csv")
    combined.to_csv(combined_path, index=False)

    summary_df = pd.DataFrame(summary_rows)
    manifest_path = os.path.join(OUT_DIR, "confidence_manifest.csv")
    summary_df.to_csv(manifest_path, index=False)

    print(f"\nCombined per-reaction table ({len(combined)} rows): {combined_path}")
    print(f"Group-level summary manifest: {manifest_path}")
    print("\nStep 10 (reconstruction-consensus layer) complete. Not proceeding "
          "further automatically -- flux_sampling_uncertainty_score, "
          "perturbation_robustness_score, and calibrated_confidence_score "
          "columns are reserved (None) pending Steps 7, 8, and 11.")


if __name__ == "__main__":
    main()
