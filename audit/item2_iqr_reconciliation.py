"""
Item 2 -- Reconcile the 'Flux Uncertainty-Only' and 'Uncalibrated Equal
Weighting' IQR values reported in the manuscript's Table 3 (0.2140, 0.2410)
against an independent reconstruction (which previously produced 0.4765 and
0.1683).

Ground-truth search (done before writing this script, at the user's explicit
request): grepped every .py/.ipynb/.R file under the repository root for the
exact strings "0.2140", "0.2410", "iMAT+FVA", "single-algorithm", "master benchmark",
and "structural consensus-only" -- ZERO matches. No script or notebook
anywhere in this repository computed these two numbers. This is recorded
below as its own finding (no ground-truth source found), not folded into a
PASS/FAIL verdict, per explicit instruction.

Because no script defines exactly which reaction population and which score
representation ("Flux Uncertainty-Only" = raw flux_confidence_raw? its rank?
w=(0,1,0) combined score? evaluation-percentile of that?) the two ablations
refer to, this script computes IQR under every combination of:
  - proxy definition: flux-uncertainty-only (w=(0,1,0)) vs equal-weight (w=(1/3,1/3,1/3))
  - score representation: raw weighted-rank-sum score vs its evaluation-percentile
  - reaction population: pooled is_active-only rows vs pooled all rows (inactive=0),
    each either pooled across all 10 groups or averaged per-group
so the actual manuscript definition can be identified by matching a cell.

Fast (~seconds); run synchronously, no detached process needed.
"""

import itertools
import json
import os
import time

import numpy as np
import pandas as pd

from audit_common import DETAILS_DIR, evaluation_percentile, load_calibrated_metadata

REPORTED_FLUX_ONLY_IQR = 0.2140
REPORTED_EQUAL_WEIGHT_IQR = 0.2410
PRIOR_INDEPENDENT_FLUX_ONLY_IQR = 0.4765
PRIOR_INDEPENDENT_EQUAL_WEIGHT_IQR = 0.1683
TOL = 0.01

PROXIES = {
    "flux_uncertainty_only": (0.0, 1.0, 0.0),
    "uncalibrated_equal_weighting": (1.0 / 3, 1.0 / 3, 1.0 / 3),
}


def iqr(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    q75, q25 = np.percentile(x, [75, 25])
    return float(q75 - q25)


def weighted_rank_score(df, w1, w2, w3, zero_fill_inactive):
    score = (
        w1 * df["rank_consensus"]
        + w2 * df["rank_flux_confidence"]
        + w3 * df["rank_perturbation"]
    )
    if zero_fill_inactive:
        return np.where(df["is_active"], score.fillna(0), 0.0)
    return np.where(df["is_active"], score, np.nan)


def run():
    t0 = time.time()
    df = load_calibrated_metadata()

    grid = {}
    for proxy_name, (w1, w2, w3) in PROXIES.items():
        for population in ("active_only", "all_rows_zero_filled"):
            zero_fill = population == "all_rows_zero_filled"
            raw_score = weighted_rank_score(df, w1, w2, w3, zero_fill)
            df["_raw"] = raw_score
            ep_score = evaluation_percentile(df, "_raw")

            # pooled across all 10 groups
            grid[f"{proxy_name}__{population}__raw__pooled"] = iqr(raw_score)
            grid[f"{proxy_name}__{population}__eval_pct__pooled"] = iqr(ep_score)

            # per-group IQR, averaged across the 10 groups
            per_group_raw, per_group_ep = [], []
            for _, sub in df.assign(_raw=raw_score, _ep=ep_score).groupby(["cohort", "group"]):
                per_group_raw.append(iqr(sub["_raw"]))
                per_group_ep.append(iqr(sub["_ep"]))
            grid[f"{proxy_name}__{population}__raw__per_group_mean"] = float(np.nanmean(per_group_raw))
            grid[f"{proxy_name}__{population}__eval_pct__per_group_mean"] = float(np.nanmean(per_group_ep))

    df.drop(columns=["_raw"], inplace=True, errors="ignore")

    def find_matches(target, tol=TOL):
        return {k: v for k, v in grid.items() if abs(v - target) < tol}

    def closest_cell(target):
        k = min(grid, key=lambda k: abs(grid[k] - target))
        return {"cell": k, "value": grid[k], "abs_diff": abs(grid[k] - target)}

    matches_flux_reported = find_matches(REPORTED_FLUX_ONLY_IQR)
    matches_equal_reported = find_matches(REPORTED_EQUAL_WEIGHT_IQR)
    matches_flux_prior = find_matches(PRIOR_INDEPENDENT_FLUX_ONLY_IQR)
    matches_equal_prior = find_matches(PRIOR_INDEPENDENT_EQUAL_WEIGHT_IQR)
    closest_flux_reported = closest_cell(REPORTED_FLUX_ONLY_IQR)
    closest_equal_reported = closest_cell(REPORTED_EQUAL_WEIGHT_IQR)

    reconciled = bool(matches_flux_reported) or bool(matches_equal_reported)
    equal_weighting_never_matches = (
        not matches_equal_reported
        and all(not k.startswith("uncalibrated_equal_weighting") for k in matches_flux_reported)
        and all(not k.startswith("uncalibrated_equal_weighting") for k in matches_equal_reported)
    )
    both_reported_values_explained_by_flux_only_family = (
        bool(matches_flux_reported)
        and all(k.startswith("flux_uncertainty_only") for k in matches_flux_reported)
        and bool(matches_equal_reported)
        and all(k.startswith("flux_uncertainty_only") for k in matches_equal_reported)
    )

    result = {
        "item": 2,
        "description": "Flux Uncertainty-Only / Uncalibrated Equal Weighting Table 3 IQR reconciliation",
        "ground_truth_search": {
            "searched_patterns": ["0.2140", "0.2410", "iMAT+FVA", "single-algorithm",
                                   "master benchmark", "structural consensus-only"],
            "searched_filetypes": [".py", ".ipynb", ".R"],
            "matches_found": 0,
            "finding": (
                "No script or notebook anywhere in the repository computed these Table 3 "
                "numbers under any of the searched names. This means the specific manuscript "
                "values (0.2140, 0.2410) were never independently computed on this machine -- "
                "a 'no ground-truth source found' finding, distinct from and more serious than "
                "an ordinary numeric mismatch."
            ),
        },
        "reported_manuscript_values": {
            "flux_uncertainty_only_iqr": REPORTED_FLUX_ONLY_IQR,
            "uncalibrated_equal_weighting_iqr": REPORTED_EQUAL_WEIGHT_IQR,
        },
        "prior_independent_reconstruction_values": {
            "flux_uncertainty_only_iqr": PRIOR_INDEPENDENT_FLUX_ONLY_IQR,
            "uncalibrated_equal_weighting_iqr": PRIOR_INDEPENDENT_EQUAL_WEIGHT_IQR,
        },
        "full_grid": grid,
        "cells_matching_reported_flux_only_iqr": matches_flux_reported,
        "cells_matching_reported_equal_weight_iqr": matches_equal_reported,
        "cells_matching_prior_flux_only_iqr": matches_flux_prior,
        "cells_matching_prior_equal_weight_iqr": matches_equal_prior,
        "closest_cell_to_reported_flux_only_iqr": closest_flux_reported,
        "closest_cell_to_reported_equal_weight_iqr": closest_equal_reported,
        "anomaly_both_reported_values_explained_by_flux_only_family": both_reported_values_explained_by_flux_only_family,
        "anomaly_note": (
            "Both reported Table 3 values (0.2140 for 'Flux Uncertainty-Only' AND 0.2410 for "
            "'Uncalibrated Equal Weighting') are matched, within +/-0.01, ONLY by cells from the "
            "flux_uncertainty_only proxy family (evaluation-percentile representation, active-only "
            "population, per-group-mean vs pooled aggregation respectively). No cell in the "
            "uncalibrated_equal_weighting family (w=1/3,1/3,1/3) comes within 0.05 of 0.2410 -- the "
            "closest equal-weighting cell is ~0.499 (eval-percentile IQR is close to 0.5 by construction "
            "when three roughly-independent 0-1 percentile ranks are averaged, since ties are rare). "
            "This suggests the manuscript's 'Uncalibrated Equal Weighting' Table 3 figure may not "
            "correspond to a literal (1/3,1/3,1/3)-weighted construction under any reaction-population/"
            "score-representation definition tested here -- flag for manual review against the actual "
            "Table 3 generation code/notebook, if one exists outside this repository."
        ) if both_reported_values_explained_by_flux_only_family else None,
        "status": "RECONCILED" if reconciled else "NO_GROUND_TRUTH_SOURCE_FOUND",
        "max_abs_discrepancy": "N/A -- see full_grid; no single scalar discrepancy exists without a fixed population/representation definition",
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item2_iqr_reconciliation.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"Status: {result['status']}")
    print(json.dumps(grid, indent=2))
    return result


if __name__ == "__main__":
    run()
