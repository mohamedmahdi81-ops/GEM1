"""
GEM1 - Step 11 follow-up: join reaction_name/subsystem metadata onto Step 11's
calibration output.

Gap identified while building docs/MANUSCRIPT_FIGURES.md (Supplementary Figure
S2 verification, 2026-07-30): `data/calibration/all_groups_calibrated.csv`
(scripts/09_calibration.py's output) carries only reaction_id, not the
descriptive reaction_name/subsystem metadata that Step 10's own confidence
CSVs (`data/confidence_scores/{cohort}__{group}__confidence.csv`) already
have -- any subsystem-level or "extended" analysis needs that metadata joined
back on. This script does exactly that join and nothing else.

Does NOT modify all_groups_calibrated.csv or any Step 10 file -- writes a new,
separate output file only.
"""
import os
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CALIBRATION_DIR = os.path.join(BASE_DIR, "data", "calibration")
CONFIDENCE_DIR = os.path.join(BASE_DIR, "data", "confidence_scores")

CALIBRATED_PATH = os.path.join(CALIBRATION_DIR, "all_groups_calibrated.csv")
OUT_PATH = os.path.join(CALIBRATION_DIR, "all_groups_calibrated_with_metadata.csv")


def main():
    calibrated = pd.read_csv(CALIBRATED_PATH)
    print(f"Loaded {CALIBRATED_PATH}: {len(calibrated)} rows, columns: {list(calibrated.columns)}")

    metadata_frames = []
    for (cohort, group), _ in calibrated.groupby(["cohort", "group"]):
        conf_path = os.path.join(CONFIDENCE_DIR, f"{cohort}__{group}__confidence.csv")
        conf = pd.read_csv(conf_path, usecols=["reaction_id", "reaction_name", "subsystem"])
        conf["cohort"] = cohort
        conf["group"] = group
        metadata_frames.append(conf)
    metadata = pd.concat(metadata_frames, ignore_index=True)
    print(f"Loaded metadata for {metadata[['cohort', 'group']].drop_duplicates().shape[0]} groups, "
          f"{len(metadata)} total (reaction_id, cohort, group) rows")

    before_rows = len(calibrated)
    joined = calibrated.merge(
        metadata, on=["cohort", "group", "reaction_id"], how="left", validate="one_to_one"
    )
    assert len(joined) == before_rows, (
        f"Row count changed by the join ({before_rows} -> {len(joined)}) -- "
        f"the merge should be exactly one_to_one, this would indicate a duplicate-key bug."
    )

    # subsystem is populated for every reaction in Human-GEM -- any missing value here
    # would indicate a real join problem. reaction_name, however, is genuinely blank for
    # some reactions in the source model itself (mostly transport/pool reactions) -- e.g.
    # 5,714 of 12,931 reactions have NaN reaction_name in the raw confidence CSVs, confirmed
    # directly against data/confidence_scores/GSE89632__healthy__confidence.csv before
    # relaxing this check. So only subsystem is asserted to have zero missing values; the
    # reaction_name null count is compared against the SAME null count already present in
    # the source metadata (not zero), to confirm the join didn't introduce any NEW nulls.
    n_missing_subsystem = joined["subsystem"].isna().sum()
    if n_missing_subsystem:
        raise ValueError(
            f"Join left {n_missing_subsystem} rows with no subsystem -- subsystem is "
            f"populated for every reaction in Human-GEM, so this indicates a real join "
            f"problem. Investigate before trusting output."
        )
    n_missing_name_joined = joined["reaction_name"].isna().sum()
    n_missing_name_source = metadata["reaction_name"].isna().sum()
    if n_missing_name_joined != n_missing_name_source:
        raise ValueError(
            f"Join changed the reaction_name null count ({n_missing_name_source} in source "
            f"metadata vs. {n_missing_name_joined} after join) -- these should be identical "
            f"since the join is one_to_one. Investigate before trusting output."
        )
    print(f"reaction_name nulls: {n_missing_name_joined} (matches source metadata's own "
          f"{n_missing_name_source} -- genuine blank names in Human-GEM, not a join defect)")

    joined.to_csv(OUT_PATH, index=False)
    print(f"\nWrote {OUT_PATH}: {len(joined)} rows, columns: {list(joined.columns)}")
    print(f"n_unique_subsystems: {joined['subsystem'].nunique()}")


if __name__ == "__main__":
    main()
