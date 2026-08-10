"""
Item 17 -- Per-axis unconstrained weight preferences. Follow-up to items
14/16 (unconstrained weight instability): runs the same unconstrained
(floor=0, step=0.05) grid search separately for EACH of the 5 biomarker
axes individually, to characterize which axes pull the fit in which
direction.

This is a single-axis PRODUCTION fit (no CV/holdout) for each axis -- i.e.
select_best_weights(df, ALL_GROUPS, {axis: reactions}, grid, cache), the
unconstrained optimum when biomarker_metric is scored against that one axis
alone across all 10 groups. urea_cycle and bcaa's single-axis results were
already computed in item16 (as the "trained on X alone" legs of its simple
LOBO) and are carried forward here unchanged, not recomputed, to avoid
redundant work and to guarantee consistency with what was already reported.

DNL uses the SAME FASN-inclusive mapping item14/item16 used
(audit_common.BIOMARKERS_FULL["dnl_scd1_fasn"]), NOT the SCD1-only mapping
from items 12/13, per explicit instruction to stay consistent with the rest
of this comparison.

Does not modify any existing GEM1 output file. Fast (~seconds -- single-axis
fits are far cheaper than the LOCO/LOBO versions in items 14/16, no nested
CV needed here since there's nothing left to hold out with 1 axis).
"""

import json
import os
import time

from audit_common import ALL_GROUPS, BIOMARKERS_FULL, DETAILS_DIR, load_calibrated_metadata
from item5_loco_lobo_refit import generate_weight_grid, select_best_weights

AXES = ["serine_glycine_shmt", "urea_cycle", "bcaa", "dnl_scd1_fasn", "choline_pc"]

# Carried forward from item16 (trained on that single axis alone) -- not recomputed.
CARRIED_FORWARD = {
    "urea_cycle": {"weights": [0.45, 0.30, 0.25], "score": None, "source": "item16 simple-LOBO 'held-out=bcaa' leg (trained on urea_cycle alone)"},
    "bcaa": {"weights": [0.25, 0.70, 0.05], "score": None, "source": "item16 simple-LOBO 'held-out=urea_cycle' leg (trained on bcaa alone)"},
}


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    grid = generate_weight_grid(floor=0.0, step=0.05)
    cache = {}
    print(f"Grid size: {len(grid)}")

    results = {}
    for axis in AXES:
        if axis in CARRIED_FORWARD:
            w = tuple(CARRIED_FORWARD[axis]["weights"])
            results[axis] = {
                "weights": list(w), "score": None,
                "recomputed": False, "source": CARRIED_FORWARD[axis]["source"],
            }
            print(f"  {axis}: {list(w)}  (carried forward from item16, not recomputed)")
            continue
        axis_ids = {axis: BIOMARKERS_FULL[axis]}
        w, score = select_best_weights(df, ALL_GROUPS, axis_ids, grid, cache)
        results[axis] = {"weights": list(w), "score": score, "recomputed": True, "source": "this script"}
        print(f"  {axis}: {list(w)}  (score={score:.4f})")

    # Cross-check urea_cycle/bcaa by recomputing them too (cheap now that the cache is warm),
    # to confirm the carried-forward values are still exactly reproducible from this script's
    # own code path, not just trusted at face value.
    cross_check = {}
    for axis in ["urea_cycle", "bcaa"]:
        axis_ids = {axis: BIOMARKERS_FULL[axis]}
        w, score = select_best_weights(df, ALL_GROUPS, axis_ids, grid, cache)
        matches = list(w) == results[axis]["weights"]
        cross_check[axis] = {"recomputed_weights": list(w), "matches_carried_forward": matches}
        print(f"  cross-check {axis}: recomputed={list(w)} matches_carried_forward={matches}")

    # Simple pattern description: cluster axes by which component dominates (>=0.5).
    def dominant_component(w):
        labels = ["w_C", "w_F", "w_P"]
        i = max(range(3), key=lambda i: w[i])
        return labels[i] if w[i] >= 0.5 else "mixed/no-dominant"

    grouping = {}
    for axis in AXES:
        grouping[axis] = dominant_component(results[axis]["weights"])

    result = {
        "item": 17,
        "description": "Per-axis unconstrained weight preferences -- single-axis production fits for all 5 biomarker axes",
        "methodology": {
            "floor": 0.0, "grid_step": 0.05, "grid_size": len(grid),
            "fit_type": "single-axis production fit (no CV/holdout), select_best_weights(df, ALL_GROUPS, {axis: reactions}, grid, cache)",
            "dnl_mapping_used": "FASN-inclusive (audit_common.BIOMARKERS_FULL['dnl_scd1_fasn']) -- same as item14/item16, NOT the SCD1-only mapping from items 12/13",
        },
        "per_axis_weights": results,
        "cross_check_urea_bcaa": cross_check,
        "dominant_component_by_axis": grouping,
        "reference_multi_axis_fits": {
            "5-axis (item14)": {"weights": [0.528, 0.05, 0.422], "note": "mean across production+LOCO+LOBO folds, consensus-dominant"},
            "2-axis urea_cycle+bcaa (item16)": {"weights": [0.358, 0.558, 0.083], "note": "mean across production+LOCO+simple-LOBO folds, flux-confidence-dominant, BCAA-driven"},
        },
        "status": "RESOLVED",
        "max_abs_discrepancy": "N/A -- descriptive characterization, not a reconciliation against a prior figure",
        "elapsed_seconds": time.time() - t0,
    }

    out_path = os.path.join(DETAILS_DIR, "item17_per_axis_weight_preferences.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {out_path}")
    print("\nSummary table:")
    print(f"{'Axis':25s} {'w_C':>6s} {'w_F':>6s} {'w_P':>6s}  Dominant")
    for axis in AXES:
        w = results[axis]["weights"]
        print(f"{axis:25s} {w[0]:6.2f} {w[1]:6.2f} {w[2]:6.2f}  {grouping[axis]}")
    return result


if __name__ == "__main__":
    run()
