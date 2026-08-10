"""
Item 5 -- Independently re-run the LOCO/LOBO weight optimization from raw
per-reaction ranks (rank_consensus, rank_flux_confidence, rank_perturbation
columns of all_groups_calibrated_with_metadata.csv -- an explicitly allowed
raw input) and biomarker labels (docs/STEP11_BIOMARKER_MAPPING.md, transcribed
into audit_common.BIOMARKERS_FULL), to confirm an independent grid-search /
nested-CV fit converges to (0.15, 0.15, 0.70) -- not just that the stored
formula is internally arithmetically consistent (already confirmed R^2=1.0
elsewhere).

This does NOT import scripts/09_calibration.py. The grid search, nested CV
loop structure, and weighted-rank combination formula are re-derived from the
plain-English docstrings traced during audit setup, then reimplemented here
from scratch.

Fast (~seconds); run synchronously, no detached process needed.
"""

import itertools
import json
import os
import time

import numpy as np
import pandas as pd

from audit_common import (
    ALL_GROUPS, BIOMARKERS_FULL, CALIBRATION_WEIGHTS_JSON, DETAILS_DIR,
    biomarker_metric, evaluation_percentile, load_calibrated_metadata,
)

WEIGHT_FLOOR = 0.15
GRID_STEP = 0.05
REPORTED_WEIGHTS = (0.15, 0.15, 0.70)


def generate_weight_grid(floor=WEIGHT_FLOOR, step=GRID_STEP):
    grid = []
    n_steps = round((1 - 3 * floor) / step) + 1
    for i in range(n_steps + 1):
        w1 = round(floor + i * step, 10)
        for j in range(n_steps + 1):
            w2 = round(floor + j * step, 10)
            w3 = round(1 - w1 - w2, 10)
            if w3 >= floor - 1e-9 and w1 + w2 <= 1 - floor + 1e-9:
                grid.append((w1, w2, round(w3, 10)))
    return grid


def combined_score_col(df, w1, w2, w3):
    score = (
        w1 * df["rank_consensus"].fillna(0)
        + w2 * df["rank_flux_confidence"].fillna(0)
        + w3 * df["rank_perturbation"].fillna(0)
    )
    return np.where(df["is_active"], score, 0.0)


def eval_pct_for_weights(df, w1, w2, w3, cache):
    key = (round(w1, 6), round(w2, 6), round(w3, 6))
    if key in cache:
        return cache[key]
    df = df.copy()
    df["_score"] = combined_score_col(df, w1, w2, w3)
    ep = evaluation_percentile(df, "_score")
    df["_ep"] = ep
    cache[key] = df
    return df


def select_best_weights(df, groups, biomarker_ids, grid, cache):
    best_w, best_score = None, -np.inf
    for (w1, w2, w3) in grid:
        scored = eval_pct_for_weights(df, w1, w2, w3, cache)
        m = biomarker_metric(scored, "_ep", biomarker_ids, groups)
        if m > best_score:
            best_score, best_w = m, (w1, w2, w3)
    return best_w, best_score


def nested_loco(df, grid, cache):
    cohorts = sorted({c for c, g in ALL_GROUPS})
    folds = []
    for test_cohort in cohorts:
        remaining = [c for c in cohorts if c != test_cohort]
        inner_scores = {}
        for w in grid:
            fold_scores = []
            for inner_test in remaining:
                inner_train = [c for c in remaining if c != inner_test]
                train_groups = [(c, g) for c, g in ALL_GROUPS if c in inner_train]
                if not train_groups:
                    continue
                scored = eval_pct_for_weights(df, *w, cache)
                fold_scores.append(biomarker_metric(scored, "_ep", BIOMARKERS_FULL, train_groups))
            inner_scores[w] = np.mean(fold_scores) if fold_scores else -np.inf
        best_w = max(inner_scores, key=inner_scores.get)
        test_groups = [(c, g) for c, g in ALL_GROUPS if c == test_cohort]
        scored = eval_pct_for_weights(df, *best_w, cache)
        outer_score = biomarker_metric(scored, "_ep", BIOMARKERS_FULL, test_groups)
        folds.append({
            "outer_test_cohort": test_cohort,
            "selected_weights": list(best_w),
            "inner_selection_score": inner_scores[best_w],
            "outer_held_out_score": outer_score,
        })
    return folds


def nested_lobo(df, grid, cache):
    biomarkers = sorted(BIOMARKERS_FULL.keys())
    folds = []
    for test_bm in biomarkers:
        remaining = [b for b in biomarkers if b != test_bm]
        inner_scores = {}
        for w in grid:
            fold_scores = []
            for inner_test in remaining:
                inner_train_bms = [b for b in remaining if b != inner_test]
                train_bm_ids = {b: BIOMARKERS_FULL[b] for b in inner_train_bms}
                if not train_bm_ids:
                    continue
                scored = eval_pct_for_weights(df, *w, cache)
                fold_scores.append(biomarker_metric(scored, "_ep", train_bm_ids, ALL_GROUPS))
            inner_scores[w] = np.mean(fold_scores) if fold_scores else -np.inf
        best_w = max(inner_scores, key=inner_scores.get)
        test_bm_ids = {test_bm: BIOMARKERS_FULL[test_bm]}
        scored = eval_pct_for_weights(df, *best_w, cache)
        outer_score = biomarker_metric(scored, "_ep", test_bm_ids, ALL_GROUPS)
        folds.append({
            "outer_test_biomarker": test_bm,
            "selected_weights": list(best_w),
            "inner_selection_score": inner_scores[best_w],
            "outer_held_out_score": outer_score,
        })
    return folds


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    grid = generate_weight_grid()
    print(f"Weight grid size: {len(grid)}")
    cache = {}

    production_w, production_score = select_best_weights(df, ALL_GROUPS, BIOMARKERS_FULL, grid, cache)
    print(f"Production weights (independent grid search, no holdout): {production_w} (score={production_score})")

    loco_folds = nested_loco(df, grid, cache)
    lobo_folds = nested_lobo(df, grid, cache)

    with open(CALIBRATION_WEIGHTS_JSON) as f:
        stored = json.load(f)
    stored_production = tuple(stored["production_weights_used"].values())

    def w_close(a, b, tol=1e-9):
        return all(abs(x - y) < tol for x, y in zip(a, b))

    all_selected = [production_w] + [f["selected_weights"] for f in loco_folds] + [f["selected_weights"] for f in lobo_folds]
    max_abs_discrepancy = max(
        max(abs(a - b) for a, b in zip(w, REPORTED_WEIGHTS)) for w in all_selected
    )
    all_match_reported = all(w_close(w, REPORTED_WEIGHTS, tol=GRID_STEP / 2 + 1e-9) for w in all_selected)

    result = {
        "item": 5,
        "description": "LOCO/LOBO weight optimization independently re-run from raw ranks + biomarker labels",
        "weight_grid_size": len(grid),
        "reported_production_weights": list(REPORTED_WEIGHTS),
        "stored_json_production_weights": list(stored_production),
        "independent_production_weights": list(production_w),
        "independent_production_score": production_score,
        "loco_folds": loco_folds,
        "lobo_folds": lobo_folds,
        "all_folds_converge_to_reported_weights": all_match_reported,
        "max_abs_weight_discrepancy": max_abs_discrepancy,
        "status": "PASS" if all_match_reported else "FAIL",
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item5_loco_lobo_refit.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"Status: {result['status']}, max_abs_weight_discrepancy={max_abs_discrepancy}")
    return result


if __name__ == "__main__":
    run()
