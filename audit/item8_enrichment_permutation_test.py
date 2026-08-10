"""
Item 8 (priority 2) -- Enrichment recovery counts and p-values for the
Flux Uncertainty-Only and Uncalibrated Equal Weighting proxy architectures.

IMPORTANT PROVENANCE NOTE, checked before writing this script: the request
that produced this item cites "the existing Section 3.8 biomarker enrichment
test" as the methodology to mirror. An exhaustive grep of every .py/.ipynb/
.R/.md file in the repository for "enrichment", "Fisher", "Benjamini",
"permutation", "20000" found exactly ONE hit outside this audit folder:
docs/MANUSCRIPT_FIGURES.md, in a section that explicitly states "Primary
statistical test: Not applicable yet -- no extended analysis has been built"
and "this entire figure is unbuilt". There is NO existing Section 3.8
enrichment test anywhere in this codebase to mirror -- this is the same class
of finding as the missing Table 3/4 sources in the prior audit round.

Given that, this script does NOT claim to reproduce an established procedure.
It implements Fisher's method + Benjamini-Hochberg correction fresh, from
first principles, exactly as specified in the request (20,000-draw Monte
Carlo permutation test per biomarker axis, BH correction across the 5 axes,
Fisher's combined p-value per architecture), and reports whatever the honest
result is.

Score construction: for each architecture (flux_uncertainty_only w=(0,1,0),
uncalibrated_equal_weighting w=(1/3,1/3,1/3)), the per-reaction score is built
as: weighted rank-sum of rank_consensus/rank_flux_confidence/rank_perturbation,
zero-filled for inactive reactions, THEN evaluation-percentile'd (rank among
ALL 12931 reactions per group, method='min'). This is not an arbitrary pick
among item2's 8-cell grid -- it is the ONE representation that matches how
GEM1's own real pipeline actually constructs biomarker_ranking_by_group.csv
from calibrated_confidence_score (itself an all-rows-zero-filled raw score,
eval-percentile'd afterward in 09_calibration.py's evaluation_percentile()) --
i.e. the only population/representation choice with an actual in-repo
precedent, rather than a choice made to fit a target number.

Permutation null: for each of 20,000 draws, each biomarker's reaction set is
replaced by a same-size random set of reaction IDs drawn (with replacement,
negligible collision probability at n=12931 reactions vs k<=9 draws) from the
full 12931-reaction universe, evaluated with the same biomarker_metric
(mean over 10 groups of max score among the drawn reactions) used for the
real biomarker sets. One-sided p-value = (1 + #null >= observed) / (1 + 20000).

Fast (~seconds, pure vectorized numpy on already-existing data -- no new
solver calls); run synchronously, no detached process needed.
"""

import json
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

from audit_common import (
    ALL_GROUPS, BIOMARKERS_FULL, DETAILS_DIR, evaluation_percentile,
    load_calibrated_metadata,
)

N_PERMUTATIONS = 20000
SEED = 42
ALPHA = 0.05

PROXIES = {
    "flux_uncertainty_only": (0.0, 1.0, 0.0),
    "uncalibrated_equal_weighting": (1.0 / 3, 1.0 / 3, 1.0 / 3),
}

REPORTED = {
    "flux_uncertainty_only": {"axes_recovered": 1, "p_value": 0.018},
    "uncalibrated_equal_weighting": {"axes_recovered": 1, "p_value": 8.4e-4},
}


def weighted_rank_score_zero_filled(df, w1, w2, w3):
    score = (
        w1 * df["rank_consensus"]
        + w2 * df["rank_flux_confidence"]
        + w3 * df["rank_perturbation"]
    )
    return np.where(df["is_active"], score.fillna(0), 0.0)


def build_score_matrix(df, w1, w2, w3):
    """Returns (score_matrix[n_groups, n_reactions], reaction_ids[n_reactions])."""
    df = df.copy()
    df["_raw"] = weighted_rank_score_zero_filled(df, w1, w2, w3)
    df["_ep"] = evaluation_percentile(df, "_raw")
    reaction_ids = sorted(df["reaction_id"].unique())
    r_to_idx = {rid: i for i, rid in enumerate(reaction_ids)}
    n_groups, n_rxn = len(ALL_GROUPS), len(reaction_ids)
    mat = np.zeros((n_groups, n_rxn))
    for gi, (cohort, group) in enumerate(ALL_GROUPS):
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        idx = sub["reaction_id"].map(r_to_idx).values
        mat[gi, idx] = sub["_ep"].values
    return mat, reaction_ids, r_to_idx


def observed_biomarker_metric(mat, r_to_idx, rxn_ids):
    idx = [r_to_idx[r] for r in rxn_ids if r in r_to_idx]
    if not idx:
        return float("nan")
    return float(mat[:, idx].max(axis=1).mean())


def null_distribution(mat, k, n_reactions, rng, n_perm=N_PERMUTATIONS):
    random_idx = rng.integers(0, n_reactions, size=(n_perm, k))
    drawn = mat[:, random_idx]  # shape (n_groups, n_perm, k)
    per_perm_max = drawn.max(axis=2)  # (n_groups, n_perm)
    return per_perm_max.mean(axis=0)  # (n_perm,)


def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


def fisher_combined(pvals):
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1.0)
    stat = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)
    p_combined = 1 - stats.chi2.cdf(stat, df)
    return float(stat), float(p_combined)


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    n_reactions_total = df["reaction_id"].nunique()

    biomarker_names = sorted(BIOMARKERS_FULL.keys())
    results = {}

    for arch_name, (w1, w2, w3) in PROXIES.items():
        rng = np.random.default_rng(SEED)
        mat, reaction_ids, r_to_idx = build_score_matrix(df, w1, w2, w3)
        n_rxn = len(reaction_ids)

        per_axis = {}
        raw_pvals = []
        for bm in biomarker_names:
            rxn_ids = BIOMARKERS_FULL[bm]
            observed = observed_biomarker_metric(mat, r_to_idx, rxn_ids)
            k = len([r for r in rxn_ids if r in r_to_idx])
            null = null_distribution(mat, k, n_rxn, rng)
            p = (1 + np.sum(null >= observed)) / (1 + N_PERMUTATIONS)
            per_axis[bm] = {
                "observed_statistic": observed,
                "n_reactions_in_biomarker": k,
                "null_mean": float(null.mean()),
                "null_std": float(null.std()),
                "p_value_raw": float(p),
            }
            raw_pvals.append(p)

        bh_adjusted = benjamini_hochberg(raw_pvals)
        for bm, q in zip(biomarker_names, bh_adjusted):
            per_axis[bm]["p_value_bh_adjusted"] = float(q)
            per_axis[bm]["recovered_at_alpha_0.05"] = bool(q < ALPHA)

        n_recovered = int(sum(per_axis[bm]["recovered_at_alpha_0.05"] for bm in biomarker_names))
        fisher_stat, fisher_p = fisher_combined(raw_pvals)

        results[arch_name] = {
            "per_axis": per_axis,
            "n_axes_recovered_of_5": n_recovered,
            "fisher_combined_statistic": fisher_stat,
            "fisher_combined_p_value": fisher_p,
            "reported_axes_recovered": REPORTED[arch_name]["axes_recovered"],
            "reported_p_value": REPORTED[arch_name]["p_value"],
            "axes_recovered_discrepancy": n_recovered - REPORTED[arch_name]["axes_recovered"],
            "p_value_abs_discrepancy": abs(fisher_p - REPORTED[arch_name]["p_value"]),
        }

    max_abs_discrepancy = max(
        max(results[a]["p_value_abs_discrepancy"], abs(results[a]["axes_recovered_discrepancy"]))
        for a in PROXIES
    )

    summary = {
        "item": 8,
        "description": "Enrichment recovery counts and p-values for Flux Uncertainty-Only and Uncalibrated Equal Weighting, independently computed (no in-repo Section 3.8 methodology exists to mirror -- confirmed by exhaustive search)",
        "ground_truth_search": {
            "searched_patterns": ["enrichment", "Fisher", "Benjamini", "permutation", "20000", "Section 3."],
            "finding": (
                "No existing enrichment-test script or notebook anywhere in the repository. The only "
                "Benjamini-Hochberg mention in the entire codebase is in docs/MANUSCRIPT_FIGURES.md, "
                "explicitly describing a figure whose 'Primary statistical test' is 'Not applicable yet "
                "-- no extended analysis has been built' and which is flagged 'this entire figure is "
                "unbuilt'. This audit's enrichment test is a first-source, independent implementation "
                "of the specified methodology (Fisher's method + BH correction, 20000-draw permutation), "
                "not a reproduction of an established in-repo procedure."
            ),
        },
        "methodology": {
            "n_permutations": N_PERMUTATIONS, "seed": SEED, "alpha": ALPHA,
            "score_construction": "weighted rank-sum, zero-filled for inactive, then evaluation-percentile "
                                   "(rank among all 12931 reactions per group) -- matches the one population/"
                                   "representation choice with an actual in-repo precedent (how "
                                   "calibrated_confidence_score -> biomarker_ranking_by_group.csv works)",
            "null_construction": "same-size random reaction sets drawn with replacement from all 12931 "
                                  "reactions, evaluated with the same max-over-isoforms-then-mean-over-groups "
                                  "statistic as the real biomarker sets",
        },
        "results_by_architecture": results,
        "max_abs_discrepancy": max_abs_discrepancy,
        "status": "PASS" if all(
            results[a]["p_value_abs_discrepancy"] < 0.05 and results[a]["axes_recovered_discrepancy"] == 0
            for a in PROXIES
        ) else "FIRST_SOURCE_NO_PRIOR_TO_RECONCILE",
        "status_meaning": (
            "No traceable prior computation exists for the reported values (see ground_truth_search) -- "
            "status reflects whether this audit's independently-computed values happen to agree with the "
            "previously reported figures, not whether a reconciliation against a verified source succeeded."
        ),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item8_enrichment_permutation_test.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {out_path}")
    for arch, r in results.items():
        print(f"\n{arch}: axes_recovered={r['n_axes_recovered_of_5']}/5 (reported {r['reported_axes_recovered']}), "
              f"fisher_p={r['fisher_combined_p_value']:.6g} (reported {r['reported_p_value']:.6g})")
    print(f"\nStatus: {summary['status']}")
    return summary


if __name__ == "__main__":
    run()
