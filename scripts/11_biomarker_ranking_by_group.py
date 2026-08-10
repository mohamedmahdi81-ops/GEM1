"""
GEM1 - Step 11 follow-up: biomarker-ranking-across-cohorts table for
manuscript Figure 3.

Gap identified in docs/MANUSCRIPT_FIGURES.md (Figure 3 verification,
2026-07-30): the evidence for "confidence scores remain consistent across
independent cohorts" requires a per-biomarker, per-group RANK (not raw score
-- see the circularity already documented for Figure 3: rank-normalization
forces group-mean calibrated_confidence_score toward ~0.5001 regardless of
true consistency, so raw-score comparisons across cohorts would be
uninformative). No script previously computed or persisted this.

For each of the 5 locked-in biomarkers (docs/STEP11_BIOMARKER_MAPPING.md),
in each of the 10 groups: take the MAX calibrated_confidence_score among its
mapped reactions (same aggregation rule already used throughout Step 11 --
see scripts/09_calibration.py's biomarker_metric), then rank that value
against ALL 12,931 reactions in the group. Percentile uses method='min'
(matches scripts/09_calibration.py's evaluation_percentile exactly, so this
table is consistent with the same methodology already used to fit and
validate the calibration weights, not a new ad hoc ranking rule).

Reads scripts/09_calibration.py's BIOMARKERS_FULL mapping directly (single
source of truth) rather than re-transcribing reaction IDs here.
"""
import os
import sys
import json
import importlib
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
calib = importlib.import_module("09_calibration")

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CALIBRATION_DIR = os.path.join(BASE_DIR, "data", "calibration")
CALIBRATED_PATH = os.path.join(CALIBRATION_DIR, "all_groups_calibrated.csv")
OUT_PATH = os.path.join(CALIBRATION_DIR, "biomarker_ranking_by_group.csv")

BIOMARKERS = calib.BIOMARKERS_FULL
ALL_GROUPS = calib.ALL_GROUPS


def main():
    df = pd.read_csv(CALIBRATED_PATH)
    print(f"Loaded {CALIBRATED_PATH}: {len(df)} rows")

    n_reactions_per_group = df.groupby(["cohort", "group"]).size()
    assert (n_reactions_per_group == 12931).all(), (
        "Expected exactly 12,931 reactions per group -- found a different count, "
        "investigate before trusting rankings."
    )

    # Percentile rank of calibrated_confidence_score among ALL reactions in each
    # group, method='min' so tied (e.g. all the zero-scored, inactive) reactions
    # get the lowest rank in their tie group -- identical methodology to
    # scripts/09_calibration.py's evaluation_percentile.
    df["_pct"] = df.groupby(["cohort", "group"])["calibrated_confidence_score"].rank(
        pct=True, method="min"
    )
    df["_ordinal_rank_desc"] = df.groupby(["cohort", "group"])["calibrated_confidence_score"].rank(
        ascending=False, method="min"
    )

    rows = []
    for biomarker_name, reaction_ids in BIOMARKERS.items():
        for cohort, group in ALL_GROUPS:
            sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
            mapped = sub[sub["reaction_id"].isin(reaction_ids)]
            if mapped.empty:
                raise ValueError(
                    f"{biomarker_name} has ZERO of its {len(reaction_ids)} mapped reactions "
                    f"present in {cohort}/{group} -- unexpected, since all groups' calibrated "
                    f"tables cover the full 12,931-reaction network. Investigate."
                )
            best = mapped.loc[mapped["calibrated_confidence_score"].idxmax()]
            rows.append({
                "biomarker": biomarker_name,
                "cohort": cohort,
                "group": group,
                "n_mapped_reactions": len(reaction_ids),
                "n_mapped_reactions_present": len(mapped),
                "max_calibrated_confidence_score": best["calibrated_confidence_score"],
                "argmax_reaction_id": best["reaction_id"],
                "percentile_rank_in_group": round(best["_pct"], 6),
                "ordinal_rank_in_group": int(best["_ordinal_rank_desc"]),
                "n_reactions_in_group": int(len(sub)),
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {OUT_PATH}: {len(out_df)} rows "
          f"({len(BIOMARKERS)} biomarkers x {len(ALL_GROUPS)} groups)")

    print("\nPer-biomarker cross-group percentile-rank summary (consistency check):")
    summary = out_df.groupby("biomarker")["percentile_rank_in_group"].agg(["mean", "std", "min", "max"])
    print(summary.to_string())


if __name__ == "__main__":
    main()
