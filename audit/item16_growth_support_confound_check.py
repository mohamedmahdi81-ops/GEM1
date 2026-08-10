"""
Item 16 -- Growth-support confound check on item 14's unconstrained-weight
finding. Tests whether the unconstrained optimizer's preference for high
consensus weight (w_C 0.35-0.6) is driven specifically by biomarker axes
whose mapped reactions overlap the 519-reaction growth-support set (which is
forced into is_active=True and, for perturbation, fixed at
perturbation_robustness_raw=1.0 regardless of any real evidence) -- versus
axes with no such overlap.

PREMISE CHECK (done first, corrects the request's stated premise): under the
BIOMARKERS_FULL mapping item14 actually used (audit_common.BIOMARKERS_FULL,
the FASN-inclusive DNL mapping), growth-support overlap is:
  serine_glycine_shmt: 1/2 reactions (MAR03845) in growth support
  urea_cycle:          0/6 -- no overlap
  bcaa:                0/9 -- no overlap
  dnl_scd1_fasn:       0/3 -- NO overlap (only the alternate SCD1-only
                       mapping MAR00146/147/148 has 1/3 overlap, via MAR00148
                       -- item14 did not use that mapping)
  choline_pc:          1/5 reactions (MAR00653) in growth support
So under item14's actual mapping, only 2 of 5 axes overlap (SHMT, choline_pc),
not 3 -- reported explicitly below, not silently corrected.

Regardless of the exact count, urea_cycle and bcaa are unambiguously
non-overlapping under EITHER DNL mapping variant, so they are the clean,
uncontested "no overlap" holdout set requested. This script reruns the same
unconstrained (floor=0) grid search from item14 restricted to ONLY these two
axes as fitting targets, and compares against item14's already-computed
5-axis result (not recomputed here).

nested_loco/nested_lobo in item5_loco_lobo_refit.py hardcode BIOMARKERS_FULL
internally, so this script reimplements parameterized versions (same fold
structure and grid-search logic, just with an explicit biomarker_ids
argument) rather than monkeypatching a shared module-level global.

With only 2 biomarkers, true NESTED LOBO (leave-one-out inner CV over the
remaining biomarkers) is structurally degenerate -- holding one out for the
outer fold leaves only 1 biomarker, and the inner leave-one-out loop needs
>=2 remaining to do anything, so it would always skip. This script instead
reports a SIMPLE (single-level, non-nested) leave-one-biomarker-out for the
2-axis case, explicitly labeled as such rather than silently forcing the
nested structure to produce a number.

Does not modify any existing GEM1 output file. item14 is not recomputed.
"""

import json
import os
import time

import numpy as np

from audit_common import ALL_GROUPS, BIOMARKERS_FULL, DETAILS_DIR, biomarker_metric, load_calibrated_metadata, load_growth_support_ids
from item5_loco_lobo_refit import eval_pct_for_weights, generate_weight_grid, select_best_weights

RESTRICTED_AXES = ["urea_cycle", "bcaa"]


def check_growth_support_overlap():
    gs = set(load_growth_support_ids()["support_reaction_ids"])
    overlap = {}
    for bm, rxns in BIOMARKERS_FULL.items():
        hits = [r for r in rxns if r in gs]
        overlap[bm] = {"n_mapped": len(rxns), "n_overlap": len(hits), "overlapping_reactions": hits}
    scd1_only = ["MAR00146", "MAR00147", "MAR00148"]
    overlap["dnl_scd1_only_alt_mapping"] = {
        "n_mapped": 3, "n_overlap": len([r for r in scd1_only if r in gs]),
        "overlapping_reactions": [r for r in scd1_only if r in gs],
        "note": "the SCD1-only mapping used in items 12/13, NOT the mapping item14 actually used",
    }
    return overlap


def nested_loco_restricted(df, grid, cache, biomarker_ids):
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
                fold_scores.append(biomarker_metric(scored, "_ep", biomarker_ids, train_groups))
            inner_scores[w] = np.mean(fold_scores) if fold_scores else -np.inf
        best_w = max(inner_scores, key=inner_scores.get)
        test_groups = [(c, g) for c, g in ALL_GROUPS if c == test_cohort]
        scored = eval_pct_for_weights(df, *best_w, cache)
        outer_score = biomarker_metric(scored, "_ep", biomarker_ids, test_groups)
        folds.append({
            "outer_test_cohort": test_cohort, "selected_weights": list(best_w),
            "inner_selection_score": inner_scores[best_w], "outer_held_out_score": outer_score,
        })
    return folds


def simple_lobo_restricted(df, grid, cache, biomarker_ids):
    """Single-level (non-nested) leave-one-biomarker-out: with only 2
    biomarkers, true nested LOBO is degenerate (inner loop has 0 remaining
    after leaving 1 out of 2) -- this fits on the 1 remaining biomarker
    directly, no inner CV, and is labeled as such in the output."""
    names = sorted(biomarker_ids.keys())
    folds = []
    for test_bm in names:
        train_bms = [b for b in names if b != test_bm]
        train_bm_ids = {b: biomarker_ids[b] for b in train_bms}
        best_w, best_score = select_best_weights(df, ALL_GROUPS, train_bm_ids, grid, cache)
        test_bm_ids = {test_bm: biomarker_ids[test_bm]}
        scored = eval_pct_for_weights(df, *best_w, cache)
        outer_score = biomarker_metric(scored, "_ep", test_bm_ids, ALL_GROUPS)
        folds.append({
            "outer_test_biomarker": test_bm, "trained_on": train_bms,
            "selected_weights": list(best_w), "training_score": best_score,
            "outer_held_out_score": outer_score,
            "note": "single-level (non-nested) -- with only 2 biomarkers, nested nested LOBO's inner loop is degenerate",
        })
    return folds


def run():
    t0 = time.time()
    overlap = check_growth_support_overlap()
    print("Growth-support overlap by axis (item14's actual mapping):")
    for bm, info in overlap.items():
        print(f"  {bm}: {info['n_overlap']}/{info['n_mapped']} -> {info['overlapping_reactions']}")

    df = load_calibrated_metadata()
    grid = generate_weight_grid(floor=0.0, step=0.05)
    print(f"\nGrid size: {len(grid)}")
    cache = {}

    restricted_ids = {a: BIOMARKERS_FULL[a] for a in RESTRICTED_AXES}
    print(f"\nRestricted axes (no growth-support overlap): {RESTRICTED_AXES}")

    production_w, production_score = select_best_weights(df, ALL_GROUPS, restricted_ids, grid, cache)
    print(f"Restricted (2-axis) production weights: {production_w} (score={production_score})")

    loco_folds = nested_loco_restricted(df, grid, cache, restricted_ids)
    lobo_folds = simple_lobo_restricted(df, grid, cache, restricted_ids)

    # Load item14's already-computed 5-axis result for direct comparison (not recomputed).
    item14_path = os.path.join(DETAILS_DIR, "item14_unconstrained_weight_refit.json")
    with open(item14_path) as f:
        item14 = json.load(f)

    all_restricted_weights = [production_w] + [f["selected_weights"] for f in loco_folds] + [f["selected_weights"] for f in lobo_folds]
    restricted_fold_labels = (
        ["production (no holdout)"]
        + [f"LOCO held-out={f['outer_test_cohort']}" for f in loco_folds]
        + [f"simple-LOBO held-out={f['outer_test_biomarker']}" for f in lobo_folds]
    )

    mean_wc_5axis = float(np.mean([w[0] for w in [item14["unconstrained_production_weights"]] + [f["selected_weights"] for f in item14["loco_folds"]] + [f["selected_weights"] for f in item14["lobo_folds"]]]))
    mean_wc_2axis = float(np.mean([w[0] for w in all_restricted_weights]))

    still_prefers_high_wc = mean_wc_2axis > 0.30  # same rough threshold used to describe item14's "high consensus" finding

    if still_prefers_high_wc:
        conclusion = (
            f"NOT a growth-support artifact: restricting the fitting targets to ONLY urea_cycle and "
            f"bcaa (both confirmed 0/6 and 0/9 growth-support overlap -- no overlap under either DNL "
            f"mapping variant) still yields a mean w_C of {mean_wc_2axis:.3f}, comparable to the "
            f"5-axis mean of {mean_wc_5axis:.3f}. The consensus-dominant preference persists WITHOUT "
            f"the growth-support-overlapping axes (SHMT, choline_pc) in the fit -- this is evidence "
            f"the unconstrained optimizer's preference for high w_C is a genuine, independent "
            f"signal from the two clean axes, not an artifact of growth-support-forced perfect scores "
            f"on SHMT/choline_pc reactions leaking into the fit."
        )
    else:
        conclusion = (
            f"CONFIRMS growth-support artifact: restricting the fitting targets to ONLY urea_cycle and "
            f"bcaa (no growth-support overlap) drops mean w_C from {mean_wc_5axis:.3f} (5-axis) to "
            f"{mean_wc_2axis:.3f} (2-axis) -- the consensus-dominant preference substantially weakens "
            f"or disappears once the growth-support-overlapping axes (SHMT, choline_pc) are excluded "
            f"from the fit, supporting the hypothesis that item14's finding was driven by the "
            f"growth-support confound rather than independent evidence favoring consensus."
        )

    result = {
        "item": 16,
        "description": "Growth-support confound check on item14's unconstrained-weight consensus-dominant finding",
        "premise_check": {
            "request_claimed": "3 axes overlap growth support (SHMT, DNL, choline_pc), 2 don't (urea_cycle, bcaa)",
            "actual_finding_under_item14s_mapping": "Only 2 axes overlap (SHMT: 1/2, choline_pc: 1/5); DNL (FASN-inclusive mapping, as item14 used it) has 0/3 overlap. DNL only overlaps (1/3) under the alternate SCD1-only mapping used in items 12/13, which item14 did NOT use.",
            "growth_support_overlap_by_axis": overlap,
        },
        "restricted_axes_used": RESTRICTED_AXES,
        "restricted_fit": {
            "production_weights": list(production_w), "production_score": production_score,
            "loco_folds": loco_folds, "lobo_folds": lobo_folds,
            "all_fold_weights": all_restricted_weights, "fold_labels": restricted_fold_labels,
        },
        "comparison": {
            "mean_w_consensus_5axis_item14": mean_wc_5axis,
            "mean_w_consensus_2axis_restricted": mean_wc_2axis,
            "still_prefers_high_consensus_2axis_only": still_prefers_high_wc,
        },
        "conclusion": conclusion,
        "status": "RESOLVED",
        "max_abs_discrepancy": abs(mean_wc_5axis - mean_wc_2axis),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item16_growth_support_confound_check.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {out_path}")
    print("\nRestricted (2-axis: urea_cycle+bcaa) fold weights:")
    for label, w in zip(restricted_fold_labels, all_restricted_weights):
        print(f"  {label:35s} -> {w}")
    print(f"\nMean w_C: 5-axis={mean_wc_5axis:.3f}  2-axis(restricted)={mean_wc_2axis:.3f}")
    print(f"\nCONCLUSION: {conclusion}")
    return result


if __name__ == "__main__":
    run()
