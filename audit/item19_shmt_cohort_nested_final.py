"""
Item 19 -- Final attempt to reconcile the original SHMT p=2.1e-5 claim via
the two-level nested (per-group -> per-cohort Fisher -> across-cohort
Fisher) combination, requested explicitly as "distinct from the two
variants already tried (flat 10-group, and one-level-across-all-groups)".

IMPORTANT: this exact procedure (per-group empirical p -> Fisher-combine
WITHIN each of the 3 cohorts -> Fisher-combine the 3 cohort-level p-values
again) is NOT a new variant -- it is precisely what
item11_calibrated_enrichment_audit.py's nested_test() function already
computed and reported as the "two-level nested" design (item11's per-axis
"nested_design" field), at N=20,000 permutations:
    raw combined p = 8.898e-06
    BH-adjusted p  = 4.449e-05
This script does not silently reuse that number -- it recomputes the
identical procedure fresh, from scratch, at N=100,000 permutations (matching
the higher precision used in item13's frozen definitive procedure) as a
clean, dedicated, directly-citable confirmation, and reports both the
N=20,000 (item11, cited) and N=100,000 (this script, fresh) results side by
side.

Does not modify any existing GEM1 output file. Fast (~seconds); synchronous.
"""

import json
import os
import time

import numpy as np
from scipy import stats

from audit_common import ALL_GROUPS, DETAILS_DIR, evaluation_percentile, load_calibrated_metadata

SEED = 42
N_PERMUTATIONS = 100_000
SHMT_REACTIONS = ["MAR03845", "MAR04792"]
REPORTED_VALUE = 2.1e-5
ORDER_OF_MAGNITUDE_TOL = 10.0  # "within an order of magnitude" = ratio <= 10x in either direction

COHORT_GROUP_INDICES = {}
for _i, (_c, _g) in enumerate(ALL_GROUPS):
    COHORT_GROUP_INDICES.setdefault(_c, []).append(_i)


def fisher_combined(pvals):
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1.0)
    stat = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)
    p_combined = 1 - stats.chi2.cdf(stat, df)
    return float(stat), float(p_combined), df


def build_score_matrix(df):
    df = df.copy()
    df["_ep"] = evaluation_percentile(df, "calibrated_confidence_score")
    reaction_ids = sorted(df["reaction_id"].unique())
    r_to_idx = {rid: i for i, rid in enumerate(reaction_ids)}
    n_groups, n_rxn = len(ALL_GROUPS), len(reaction_ids)
    mat = np.zeros((n_groups, n_rxn))
    for gi, (cohort, group) in enumerate(ALL_GROUPS):
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        idx = sub["reaction_id"].map(r_to_idx).values
        mat[gi, idx] = sub["_ep"].values
    return mat, reaction_ids, r_to_idx


def two_level_nested_test(mat, r_to_idx, rxn_ids, n_reactions, rng, n_perm):
    idx = [r_to_idx[r] for r in rxn_ids if r in r_to_idx]
    k = len(idx)
    per_group = []
    for gi, (cohort, group) in enumerate(ALL_GROUPS):
        observed_g = float(mat[gi, idx].max())
        random_idx = rng.integers(0, n_reactions, size=(n_perm, k))
        null_g = mat[gi, random_idx].max(axis=1)
        n_ge = int(np.sum(null_g >= observed_g))
        p_g = (1 + n_ge) / (1 + n_perm)
        per_group.append({"cohort": cohort, "group": group, "observed_statistic": observed_g,
                           "n_null_ge_observed": n_ge, "p_value": float(p_g)})

    per_cohort = {}
    for cohort, gis in COHORT_GROUP_INDICES.items():
        group_ps = [per_group[gi]["p_value"] for gi in gis]
        stat, p_c, df = fisher_combined(group_ps)
        per_cohort[cohort] = {"n_groups": len(gis), "group_p_values": group_ps,
                               "fisher_statistic": stat, "df": df, "p_value": p_c}

    cohort_ps = [per_cohort[c]["p_value"] for c in per_cohort]
    stat_final, p_final, df_final = fisher_combined(cohort_ps)

    return {
        "per_group": per_group, "per_cohort": per_cohort,
        "final_fisher_statistic": stat_final, "final_df": df_final, "p_value": p_final,
    }


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    mat, reaction_ids, r_to_idx = build_score_matrix(df)
    n_reactions = len(reaction_ids)

    rng = np.random.default_rng(SEED)
    fresh_100k = two_level_nested_test(mat, r_to_idx, SHMT_REACTIONS, n_reactions, rng, N_PERMUTATIONS)

    item11_path = os.path.join(DETAILS_DIR, "item11_calibrated_enrichment_audit.json")
    with open(item11_path) as f:
        item11 = json.load(f)
    prior_20k = item11["axes"]["serine_glycine_shmt"]["nested_design"]

    def ratio_within_order_of_magnitude(p):
        r = max(p / REPORTED_VALUE, REPORTED_VALUE / p)
        return r, r <= ORDER_OF_MAGNITUDE_TOL

    ratio_100k, within_100k = ratio_within_order_of_magnitude(fresh_100k["p_value"])
    ratio_20k, within_20k = ratio_within_order_of_magnitude(prior_20k["p_value"])

    already_tried_note = (
        "This exact procedure (per-group p -> Fisher-combine within each cohort -> Fisher-combine "
        "across cohorts) is NOT a new, previously-untried variant -- it is precisely "
        "item11_calibrated_enrichment_audit.py's 'nested_design' (two-level nested), already computed "
        "and reported at N=20,000: raw p=8.898e-06. The 'flat 10-group' and 'one-level-across-all-"
        "groups' variants referenced in the request are DIFFERENT designs (item13's frozen flat "
        "procedure, and item11/item12's nested_test_one_level respectively) -- this two-level "
        "cohort-then-overall design was already the closest match found in item11, not an untested "
        "possibility."
    )

    if within_100k:
        conclusion = (
            f"MATCHES within an order of magnitude -- DO NOT retract. The two-level nested "
            f"(cohort-then-overall) Fisher combination, recomputed fresh at N=100,000 permutations, "
            f"gives a raw combined p of {fresh_100k['p_value']:.4g}, a ratio of {ratio_100k:.2f}x from "
            f"the reported 2.1e-5 (well within the 10x order-of-magnitude threshold; for comparison, "
            f"the N=20,000 result -- already reported in item11 -- was {ratio_20k:.2f}x off). This is "
            f"the best-matching mechanism found across every variant tested in this audit (flat, "
            f"one-level nested, two-level nested), and it lands close enough to 2.1e-5 that a nested "
            f"Fisher combination of this specific structure is a credible, defensible explanation for "
            f"the original figure. It is not an exact bit-for-bit reproduction (no source script exists "
            f"to confirm the precise seed/draw count/tie-breaking used originally), but 'unrecoverable, "
            f"retract' is NOT the conclusion this evidence supports -- 'plausible mechanism identified, "
            f"exact reproduction not possible without source code' is the accurate summary."
        )
    else:
        conclusion = (
            f"Does NOT match within an order of magnitude even under this two-level nested design at "
            f"N=100,000 (ratio {ratio_100k:.2f}x from 2.1e-5). Per instruction, this value should be "
            f"treated as UNRECOVERABLE from any procedure tested in this audit, and formally retracted "
            f"or corrected in the manuscript rather than pursued further."
        )

    result = {
        "item": 19,
        "description": "Final SHMT p=2.1e-5 reconciliation attempt via two-level nested (cohort-then-overall) Fisher combination",
        "already_tried_note": already_tried_note,
        "fresh_100k_result": fresh_100k,
        "prior_20k_result_cited_from_item11": {
            "raw_p": prior_20k["p_value"], "bh_adjusted_p": prior_20k["bh_adjusted_p"],
            "per_cohort_p": {c: v["p_value"] for c, v in prior_20k["per_cohort"].items()},
        },
        "reported_value": REPORTED_VALUE,
        "ratio_100k_vs_reported": ratio_100k, "within_order_of_magnitude_100k": within_100k,
        "ratio_20k_vs_reported": ratio_20k, "within_order_of_magnitude_20k": within_20k,
        "conclusion": conclusion,
        "status": "MATCHES_WITHIN_ORDER_OF_MAGNITUDE -- NOT unrecoverable" if within_100k else "UNRECOVERABLE -- RECOMMEND RETRACTION",
        "max_abs_discrepancy": abs(fresh_100k["p_value"] - REPORTED_VALUE),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item19_shmt_cohort_nested_final.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"\n{already_tried_note}\n")
    print(f"Fresh N=100,000 two-level nested raw p: {fresh_100k['p_value']:.6g}")
    print(f"Prior N=20,000 (item11) raw p: {prior_20k['p_value']:.6g}, BH-adjusted: {prior_20k['bh_adjusted_p']:.6g}")
    print(f"Reported: {REPORTED_VALUE:.6g}")
    print(f"Ratio (100k): {ratio_100k:.3f}x -- within order of magnitude: {within_100k}")
    print(f"\nCONCLUSION: {conclusion}")
    return result


if __name__ == "__main__":
    run()
