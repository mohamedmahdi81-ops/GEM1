"""
Item 9 (priority 3) -- Section 3.4 component-level SD sensitivity.

Manuscript claims:
  calibrated (0.15/0.15/0.70)          SD=0.2002
  perturbation-dominant (0.00/0.00/1.00) SD=0.2900
  equal weighting (0.33/0.33/0.33)       SD=0.1236

The calibrated (0.15/0.15/0.70) case is directly checkable from the existing
`calibrated_confidence_score` column of all_groups_calibrated_with_metadata.csv
(an allowed raw input -- no scoring/calibration module imported, this is the
persisted data, not the code that produced it). The other two architectures
are reconstructed the same rank-based way item2_iqr_reconciliation.py used for
Table 3: weighted sum of rank_consensus/rank_flux_confidence/rank_perturbation,
under the same grid of population (active-only vs all-rows-zero-filled, pooled
vs per-group-mean) and representation (raw weighted-rank-sum vs
evaluation-percentile) choices, since no script anywhere defines which one the
manuscript figure used (same "no traceable source" situation as Table 3).

Fast (~seconds); run synchronously, no detached process needed.
"""

import json
import os
import time

import numpy as np
import pandas as pd

from audit_common import DETAILS_DIR, evaluation_percentile, load_calibrated_metadata

REPORTED_SD = {
    "calibrated": 0.2002,
    "perturbation_dominant": 0.2900,
    "equal_weighting": 0.1236,
}
WEIGHTS = {
    "calibrated": (0.15, 0.15, 0.70),
    "perturbation_dominant": (0.0, 0.0, 1.0),
    "equal_weighting": (1.0 / 3, 1.0 / 3, 1.0 / 3),
}
TOL = 0.01


def weighted_rank_score(df, w1, w2, w3, zero_fill_inactive):
    score = (
        w1 * df["rank_consensus"]
        + w2 * df["rank_flux_confidence"]
        + w3 * df["rank_perturbation"]
    )
    if zero_fill_inactive:
        return np.where(df["is_active"], score.fillna(0), 0.0)
    return np.where(df["is_active"], score, np.nan)


def sd(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(np.std(x)) if len(x) else float("nan")


def run():
    t0 = time.time()
    df = load_calibrated_metadata()

    # --- Direct check: calibrated_confidence_score column as persisted ---
    direct_calibrated_sd_pooled = sd(df["calibrated_confidence_score"])
    per_group_sds = []
    for _, sub in df.groupby(["cohort", "group"]):
        per_group_sds.append(sd(sub["calibrated_confidence_score"]))
    direct_calibrated_sd_per_group_mean = float(np.nanmean(per_group_sds))

    grid = {}
    for arch_name, (w1, w2, w3) in WEIGHTS.items():
        for population in ("active_only", "all_rows_zero_filled"):
            zero_fill = population == "all_rows_zero_filled"
            raw_score = weighted_rank_score(df, w1, w2, w3, zero_fill)
            df["_raw"] = raw_score
            ep_score = evaluation_percentile(df, "_raw")

            grid[f"{arch_name}__{population}__raw__pooled"] = sd(raw_score)
            grid[f"{arch_name}__{population}__eval_pct__pooled"] = sd(ep_score)

            per_group_raw, per_group_ep = [], []
            for _, sub in df.assign(_raw=raw_score, _ep=ep_score).groupby(["cohort", "group"]):
                per_group_raw.append(sd(sub["_raw"]))
                per_group_ep.append(sd(sub["_ep"]))
            grid[f"{arch_name}__{population}__raw__per_group_mean"] = float(np.nanmean(per_group_raw))
            grid[f"{arch_name}__{population}__eval_pct__per_group_mean"] = float(np.nanmean(per_group_ep))

    df.drop(columns=["_raw"], inplace=True, errors="ignore")

    grid["calibrated__DIRECT_COLUMN__pooled"] = direct_calibrated_sd_pooled
    grid["calibrated__DIRECT_COLUMN__per_group_mean"] = direct_calibrated_sd_per_group_mean

    def closest_cell(target, prefix=None):
        candidates = {k: v for k, v in grid.items() if prefix is None or k.startswith(prefix)}
        k = min(candidates, key=lambda k: abs(candidates[k] - target))
        return {"cell": k, "value": candidates[k], "abs_diff": abs(candidates[k] - target)}

    closest = {arch: closest_cell(target, prefix=arch if arch != "calibrated" else None)
               for arch, target in REPORTED_SD.items()}

    def find_matches(target, tol=TOL):
        return {k: v for k, v in grid.items() if abs(v - target) < tol}

    matches = {arch: find_matches(target) for arch, target in REPORTED_SD.items()}

    direct_discrepancy = min(
        abs(direct_calibrated_sd_pooled - REPORTED_SD["calibrated"]),
        abs(direct_calibrated_sd_per_group_mean - REPORTED_SD["calibrated"]),
    )

    result = {
        "item": 9,
        "description": "Section 3.4 component-level SD sensitivity, independently recomputed from raw rank components",
        "reported_sd": REPORTED_SD,
        "direct_column_check": {
            "calibrated_confidence_score_sd_pooled": direct_calibrated_sd_pooled,
            "calibrated_confidence_score_sd_per_group_mean": direct_calibrated_sd_per_group_mean,
            "reported": REPORTED_SD["calibrated"],
            "discrepancy": direct_discrepancy,
        },
        "full_grid": grid,
        "closest_cell_per_architecture": closest,
        "cells_matching_reported_within_tol": matches,
        "max_abs_discrepancy": max(v["abs_diff"] for v in closest.values()),
        "pass_threshold": TOL,
        "status": "PASS" if all(v["abs_diff"] < TOL for v in closest.values()) else "PARTIAL",
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item9_sd_sensitivity.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(json.dumps(closest, indent=2))
    print(f"Status: {result['status']}")
    return result


if __name__ == "__main__":
    run()
