"""
GEM1 - Step 9: Reconstruction-algorithm consensus scoring.

Design locked in 2026-07-17: the primary consensus is built
from FASTCORE + iMAT + tINIT only (strict all-3-agree). GIMME is excluded
from the vote itself -- its active-reaction sets are 80-96% of the whole
network (vs. 39-64% for the other three), so folding it into an equal-weight
vote would dilute rather than sharpen the consensus signal. Instead, GIMME is
reported per group as a corroboration check against the primary consensus,
reusing the overlap numbers already computed in
data/context_specific_models/gimme_consensus_overlap.csv (verified accurate
against the raw JSONs in that same session -- not recomputed here).

Pure set comparison over the already-extracted JSONs in
data/context_specific_models/ -- no solver calls, no troppo/cobra dependency,
runs in gem1-main (or any env with pandas).
"""

import os
import json
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
MODELS_DIR = os.path.join(BASE_DIR, "data", "context_specific_models")
OUT_DIR = os.path.join(BASE_DIR, "data", "consensus_scores")
os.makedirs(OUT_DIR, exist_ok=True)

CORE_ALGORITHMS = ["fastcore", "imat", "tinit"]
GIMME_OVERLAP_CSV = os.path.join(MODELS_DIR, "gimme_consensus_overlap.csv")


def load_active_set(cohort, group, algorithm):
    path = os.path.join(MODELS_DIR, f"{cohort}__{group}__{algorithm}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return set(json.load(f)["active_reactions"])


def load_manifest_combos():
    manifest_path = os.path.join(MODELS_DIR, "extraction_manifest.csv")
    manifest = pd.read_csv(manifest_path)
    return manifest[["cohort", "group"]].drop_duplicates().values.tolist()


def compute_consensus(cohort, group):
    algo_sets = {algo: load_active_set(cohort, group, algo) for algo in CORE_ALGORITHMS}
    available = {algo: s for algo, s in algo_sets.items() if s is not None}
    missing = [algo for algo, s in algo_sets.items() if s is None]

    n_available = len(available)
    if n_available < 2:
        raise ValueError(
            f"{cohort}/{group}: only {n_available} of {CORE_ALGORITHMS} available, "
            "need at least 2 to compute any consensus."
        )

    strict_consensus = set.intersection(*available.values())
    consensus_definition_used = "strict" if n_available == 3 else "reduced"

    from collections import Counter
    vote_counts = Counter()
    for s in available.values():
        vote_counts.update(s)
    majority_thresh = 2
    majority_consensus = set(r for r, c in vote_counts.items() if c >= majority_thresh)

    return {
        "strict_consensus": strict_consensus,
        "majority_consensus": majority_consensus,
        "consensus_definition_used": consensus_definition_used,
        "algorithms_used": sorted(available.keys()),
        "algorithms_missing": missing,
    }


def load_gimme_overlap():
    df = pd.read_csv(GIMME_OVERLAP_CSV)
    df = df.rename(columns={
        "pct_consensus_retained_by_gimme": "gimme_pct_retained",
        "pct_gimme_reactions_that_are_gimme_only": "gimme_pct_only",
    })
    return df.set_index(["cohort", "group"])[["gimme_pct_retained", "gimme_pct_only"]]


def main():
    combos = load_manifest_combos()
    gimme_overlap = load_gimme_overlap()

    manifest_rows = []
    for cohort, group in combos:
        result = compute_consensus(cohort, group)

        strict_ids = sorted(result["strict_consensus"])
        json_path = os.path.join(OUT_DIR, f"{cohort}__{group}__consensus.json")
        with open(json_path, "w") as f:
            json.dump({
                "cohort": cohort,
                "group": group,
                "consensus_definition_used": result["consensus_definition_used"],
                "algorithms_used": result["algorithms_used"],
                "algorithms_missing": result["algorithms_missing"],
                "strict_consensus_reactions": strict_ids,
                "n_strict_consensus": len(strict_ids),
                "majority_consensus_reactions": sorted(result["majority_consensus"]),
                "n_majority_consensus": len(result["majority_consensus"]),
            }, f)

        try:
            gimme_pct_retained, gimme_pct_only = gimme_overlap.loc[(cohort, group)]
        except KeyError:
            gimme_pct_retained, gimme_pct_only = None, None
            print(f"  WARNING: no GIMME overlap row found for {cohort}/{group} in "
                  f"{GIMME_OVERLAP_CSV} -- leaving blank.")

        flag = (
            f"REDUCED (missing: {', '.join(result['algorithms_missing'])})"
            if result["consensus_definition_used"] == "reduced"
            else "strict (all 3 algorithms available)"
        )
        print(f"{cohort}/{group}: n_strict={len(strict_ids)} "
              f"n_majority={len(result['majority_consensus'])} "
              f"[{flag}] gimme_retained={gimme_pct_retained} gimme_only={gimme_pct_only}")

        manifest_rows.append({
            "cohort": cohort,
            "group": group,
            "n_strict_consensus": len(strict_ids),
            "n_majority_consensus": len(result["majority_consensus"]),
            "consensus_definition_used": result["consensus_definition_used"],
            "gimme_pct_retained": gimme_pct_retained,
            "gimme_pct_only": gimme_pct_only,
        })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(OUT_DIR, "consensus_manifest.csv")
    manifest_df.to_csv(manifest_path, index=False)
    print(f"\nConsensus manifest saved to {manifest_path}")
    print(f"Per-group consensus reaction-ID JSONs saved to {OUT_DIR}")
    print("Next: Step 10 (hierarchical confidence engine) consumes "
          "strict_consensus_reactions from each group's JSON.")


if __name__ == "__main__":
    main()
