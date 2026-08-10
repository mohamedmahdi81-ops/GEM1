"""
Item 14 -- Unconstrained weight re-fit. Re-runs the exact same LOCO/LOBO
grid-search + nested-CV machinery as item5_loco_lobo_refit.py (imported
directly here -- this is this audit's own code, not GEM1's scoring modules,
so reuse is appropriate and keeps the comparison apples-to-apples), but with
the 0.15-per-weight floor removed: each weight can range from 0 to 1 in 0.05
steps, still constrained to sum to 1.

Purpose: determine whether the reported production weights (0.15, 0.15, 0.70)
reflect a genuine interior optimum that the 0.15 floor happens not to bind
on, or whether the floor was actively constraining the search away from a
true optimum that wants w_C and/or w_F below 0.15 (i.e. toward a boundary
solution with w_P closer to or at 1).

Does not modify any existing GEM1 output file. Fast (~seconds); synchronous.
"""

import json
import os
import time

from item5_loco_lobo_refit import (
    generate_weight_grid, select_best_weights, nested_loco, nested_lobo,
)
from audit_common import ALL_GROUPS, BIOMARKERS_FULL, CALIBRATION_WEIGHTS_JSON, DETAILS_DIR, load_calibrated_metadata

REPORTED_WEIGHTS = (0.15, 0.15, 0.70)
UNCONSTRAINED_FLOOR = 0.0
GRID_STEP = 0.05
BINDING_THRESHOLD = 0.15  # if unconstrained optimum < this for w_C or w_F, the floor was binding


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    grid = generate_weight_grid(floor=UNCONSTRAINED_FLOOR, step=GRID_STEP)
    print(f"Unconstrained weight grid size: {len(grid)} (vs. 78 under the 0.15 floor)")
    cache = {}

    production_w, production_score = select_best_weights(df, ALL_GROUPS, BIOMARKERS_FULL, grid, cache)
    print(f"Unconstrained production weights (no holdout): {production_w} (score={production_score})")

    loco_folds = nested_loco(df, grid, cache)
    lobo_folds = nested_lobo(df, grid, cache)

    all_selected = [production_w] + [f["selected_weights"] for f in loco_folds] + [f["selected_weights"] for f in lobo_folds]
    fold_labels = (
        ["production (no holdout)"]
        + [f"LOCO held-out={f['outer_test_cohort']}" for f in loco_folds]
        + [f"LOBO held-out={f['outer_test_biomarker']}" for f in lobo_folds]
    )

    per_fold_analysis = []
    for label, w in zip(fold_labels, all_selected):
        w_c, w_f, w_p = w
        is_boundary = (w_c == 0.0 or w_f == 0.0 or w_p >= 0.95)
        matches_reported = all(abs(a - b) < 1e-9 for a, b in zip(w, REPORTED_WEIGHTS))
        floor_would_have_bound = w_c < BINDING_THRESHOLD - 1e-9 or w_f < BINDING_THRESHOLD - 1e-9
        per_fold_analysis.append({
            "fold": label, "unconstrained_weights": list(w),
            "matches_reported_0.15_0.15_0.70": matches_reported,
            "is_boundary_solution": is_boundary,
            "floor_0.15_would_have_been_binding": floor_would_have_bound,
        })

    n_matching_reported = sum(a["matches_reported_0.15_0.15_0.70"] for a in per_fold_analysis)
    n_boundary = sum(a["is_boundary_solution"] for a in per_fold_analysis)
    n_floor_binding = sum(a["floor_0.15_would_have_been_binding"] for a in per_fold_analysis)

    if n_matching_reported == len(per_fold_analysis):
        floor_conclusion = (
            "The floor was NOT binding: every fold's unconstrained optimum lands EXACTLY at "
            "(0.15, 0.15, 0.70) even with the floor removed and weights free to range down to 0. "
            "This means 0.15/0.15/0.70 is the genuine interior optimum of the underlying objective "
            "(biomarker_metric), not an artifact of the floor constraint -- the search would have "
            "landed there regardless of whether the floor existed."
        )
    elif n_floor_binding == len(per_fold_analysis):
        floor_conclusion = (
            "The floor WAS actively binding on every fold: the unconstrained optimum consistently "
            "wants w_C and/or w_F below 0.15 (often at or near 0), pushing toward a boundary solution "
            "with w_P higher than 0.70. The reported (0.15, 0.15, 0.70) reflects the floor constraint, "
            "not the underlying objective's true unconstrained preference."
        )
    else:
        floor_conclusion = (
            f"Mixed result: {n_floor_binding}/{len(per_fold_analysis)} folds show the floor would have "
            f"been binding (unconstrained w_C and/or w_F below 0.15), while the rest land at or above "
            f"0.15 anyway. See per_fold_analysis for which folds diverge and how."
        )

    result = {
        "item": 14,
        "description": "Unconstrained (0-to-1, no 0.15 floor) LOCO/LOBO weight re-fit, same grid-search/nested-CV machinery as item 5",
        "methodology": {
            "floor_used": UNCONSTRAINED_FLOOR, "grid_step": GRID_STEP,
            "grid_size": len(grid),
            "reused_code": "item5_loco_lobo_refit.py's generate_weight_grid/select_best_weights/nested_loco/nested_lobo, "
                            "imported directly (this audit's own code, not GEM1's scoring modules) -- only the floor parameter changes",
        },
        "reported_weights": list(REPORTED_WEIGHTS),
        "unconstrained_production_weights": list(production_w),
        "unconstrained_production_score": production_score,
        "loco_folds": loco_folds,
        "lobo_folds": lobo_folds,
        "per_fold_analysis": per_fold_analysis,
        "n_folds_matching_reported": n_matching_reported,
        "n_folds_total": len(per_fold_analysis),
        "n_folds_boundary_solution": n_boundary,
        "n_folds_where_015_floor_would_bind": n_floor_binding,
        "floor_binding_conclusion": floor_conclusion,
        "status": "RESOLVED",
        "max_abs_discrepancy": max(
            max(abs(a - b) for a, b in zip(w, REPORTED_WEIGHTS)) for w in all_selected
        ),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item14_unconstrained_weight_refit.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}\n")
    for a in per_fold_analysis:
        print(f"  {a['fold']:35s} -> {a['unconstrained_weights']}  "
              f"matches_reported={a['matches_reported_0.15_0.15_0.70']}  boundary={a['is_boundary_solution']}")
    print(f"\n{n_matching_reported}/{len(per_fold_analysis)} folds match reported weights exactly.")
    print(f"\nCONCLUSION: {floor_conclusion}")
    return result


if __name__ == "__main__":
    run()
