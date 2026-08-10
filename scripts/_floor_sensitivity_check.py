"""
One-off diagnostic (not part of the production pipeline): re-run the nested
LOCO/LOBO validation for the full_isoform biomarker set at a lower weight
floor, to check whether 0.15 is artificially restraining the fit or close to
a near-optimal decomposable solution. Reuses 09_calibration.py's functions
unmodified; only the grid's floor differs.

2026-07-30 update (Figure 5 reproducibility gap, docs/MANUSCRIPT_FIGURES.md):
this originally only printed to stdout, so its numbers lived only as ad hoc
console output, not in any regenerable structured file. Now also
writes data/calibration/floor_sensitivity_results.json. Also fixes a real bug
found while doing this: the original "production weights at this floor" line
called m.select_best_weights(), which internally calls
m.generate_weight_grid() with NO floor argument -- always defaulting to
09_calibration.py's WEIGHT_FLOOR (0.15) regardless of which floor this loop
was actually testing, so that one line printed the identical (0.15,0.15,0.7)
at every floor. Fixed below by building the floor-specific grid explicitly
and reusing 09_calibration.py's own selection logic (max over
biomarker_metric) against it, instead of calling the un-parameterized helper.
"""
import os
import sys
import json
import importlib
import numpy as np

sys.path.insert(0, "scripts")
m = importlib.import_module("09_calibration")

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
OUT_PATH = os.path.join(BASE_DIR, "data", "calibration", "floor_sensitivity_results.json")


def select_best_weights_at_floor(df, biomarker_mapping, floor):
    """Same selection rule as 09_calibration.py's select_best_weights (grid
    search maximizing biomarker_metric over ALL_GROUPS), but with an explicit
    floor -- fixes the bug described in the module docstring above."""
    grid = generate_grid(floor)
    best_weights, best_score = None, -np.inf
    for w in grid:
        df_pct = m.evaluation_percentile(df, w)
        score = m.biomarker_metric(df_pct, biomarker_mapping, m.ALL_GROUPS)
        if score > best_score:
            best_score, best_weights = score, w
    return best_weights, best_score


def generate_grid(floor, step=0.05):
    return m.generate_weight_grid(step=step, floor=floor)


def nested_loco_custom(df, biomarker_mapping, floor):
    grid = generate_grid(floor)
    cohorts = list(m.GROUPS.keys())
    outer_results = []
    for outer_test_cohort in cohorts:
        train_cohorts = [c for c in cohorts if c != outer_test_cohort]
        outer_test_groups = [(c, g) for c, g in m.ALL_GROUPS if c == outer_test_cohort]

        inner_scores = {}
        for w in grid:
            df_pct = m.evaluation_percentile(df, w)
            vals = []
            for inner_val_cohort in train_cohorts:
                inner_val_groups = [(c, g) for c, g in m.ALL_GROUPS if c == inner_val_cohort]
                vals.append(m.biomarker_metric(df_pct, biomarker_mapping, inner_val_groups))
            inner_scores[w] = float(np.mean(vals))
        best_w = max(inner_scores, key=inner_scores.get)

        df_pct_outer = m.evaluation_percentile(df, best_w)
        outer_score = m.biomarker_metric(df_pct_outer, biomarker_mapping, outer_test_groups)
        outer_results.append({"outer_test_cohort": outer_test_cohort, "weights": best_w, "outer_score": outer_score})
    return outer_results, float(np.mean([r["outer_score"] for r in outer_results]))


def nested_lobo_custom(df, biomarker_mapping, floor):
    grid = generate_grid(floor)
    names = list(biomarker_mapping.keys())
    outer_results = []
    for outer_test_bio in names:
        train_bios = [b for b in names if b != outer_test_bio]
        inner_scores = {}
        for w in grid:
            df_pct = m.evaluation_percentile(df, w)
            vals = []
            for inner_val_bio in train_bios:
                vals.append(m.biomarker_metric(df_pct, biomarker_mapping, m.ALL_GROUPS, [inner_val_bio]))
            inner_scores[w] = float(np.mean(vals))
        best_w = max(inner_scores, key=inner_scores.get)

        df_pct_outer = m.evaluation_percentile(df, best_w)
        outer_score = m.biomarker_metric(df_pct_outer, biomarker_mapping, m.ALL_GROUPS, [outer_test_bio])
        outer_results.append({"outer_test_biomarker": outer_test_bio, "weights": best_w, "outer_score": outer_score})
    return outer_results, float(np.mean([r["outer_score"] for r in outer_results]))


def main():
    growth_support = m.load_growth_support()
    rows = [m.load_group_raw_scores(c, g, growth_support) for c, g in m.ALL_GROUPS]
    import pandas as pd
    df = pd.concat(rows, ignore_index=True)
    df = m.add_rank_normalization(df)

    results = {}
    for floor in [0.15, 0.05, 0.0]:
        print(f"\n=== floor={floor} ===")
        loco_folds, loco_mean = nested_loco_custom(df, m.BIOMARKERS_FULL, floor)
        for f in loco_folds:
            print(f"  LOCO outer_test={f['outer_test_cohort']:10} weights={f['weights']} outer_score={f['outer_score']:.4f}")
        print(f"  LOCO mean outer-held-out: {loco_mean:.4f}")

        lobo_folds, lobo_mean = nested_lobo_custom(df, m.BIOMARKERS_FULL, floor)
        for f in lobo_folds:
            print(f"  LOBO outer_test={f['outer_test_biomarker']:24} weights={f['weights']} outer_score={f['outer_score']:.4f}")
        print(f"  LOBO mean outer-held-out: {lobo_mean:.4f}")

        prod_w, prod_score = select_best_weights_at_floor(df, m.BIOMARKERS_FULL, floor)
        print(f"  production weights at this floor (bug-fixed, floor correctly applied): "
              f"{prod_w} (score {prod_score:.4f})")

        results[str(floor)] = {
            "loco_folds": loco_folds,
            "loco_mean": loco_mean,
            "lobo_folds": lobo_folds,
            "lobo_mean": lobo_mean,
            "production_weights_at_this_floor": prod_w,
            "production_selection_score_at_this_floor": prod_score,
        }

    with open(OUT_PATH, "w") as f:
        json.dump({
            "floors_tested": [0.15, 0.05, 0.0],
            "biomarker_set": "full_isoform",
            "results_by_floor": results,
            "note": (
                "Persisted 2026-07-30 to close the Figure 5 reproducibility gap "
                "(docs/MANUSCRIPT_FIGURES.md) -- previously these numbers existed "
                "only as ad hoc console output, not a persisted file. Regenerate "
                "by re-running this script; values should match "
                "docs/MANUSCRIPT_FIGURES.md's Figure 5 floor-sensitivity table."
            ),
        }, f, indent=2, default=str)
    print(f"\nWrote {OUT_PATH}")


if __name__ == "__main__":
    main()
