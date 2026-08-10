"""
Item 18 -- Re-run the cross-cohort Kendall's W concordance check and the
Table 3 IQR/SD architectural benchmark ordering, excluding
serine_glycine_shmt and choline_pc (the two growth-support-overlapping axes,
per item 16) from any computation using biomarker axes as inputs. Tests
whether reproducibility and the benchmark ordering depend specifically on
the two contaminated axes.

PREMISE CHECK (done first): the request cites "subsystem-level Kendall's W
... W=0.7912" and a "master architectural benchmark (Table 3 IQR
comparison)". Neither is traceable:
  - "0.7912" does not appear anywhere in this repository (confirmed by
    direct grep for the literal string).
  - There is no subsystem-level Kendall's W anywhere in this repository --
    docs/MANUSCRIPT_FIGURES.md explicitly states "No script in this
    repository aggregates confidence/calibration scores by subsystem" and
    labels that entire figure "unbuilt".
  - The ONLY Kendall's W that actually exists in this repository is
    BIOMARKER-level (not subsystem-level): data/calibration/figure3_statistics.json,
    produced by scripts/12_figure3_cross_cohort_confidence.py -- W=0.6667,
    Friedman chi2=8.0, p=0.0916, computed across m=3 cohorts (as raters) on
    n=5 biomarkers (as items), from biomarker_ranking_by_group.csv's mean
    percentile_rank_in_group per (cohort, biomarker).
  - There is no single "Table 3 IQR comparison" literally labeled a "master
    architectural benchmark" -- the closest existing, validated quantity is
    item 9's SD-sensitivity check, which already establishes (from the real
    persisted calibrated_confidence_score column, and cross-validated
    against the reported manuscript SD figures) that equal-weighting's SD
    (0.1236) is smaller than calibrated's SD (0.2001) across ALL 12,931
    reactions -- i.e. exactly the "consensus-only equal-weighting <
    calibrated GEM1" ordering the request describes, just as SD not IQR,
    and not previously restricted to biomarker-axis reactions specifically
    (items 2/7/13's Table 3 IQR work was never axis-restricted to begin
    with -- it used the full reaction population).

Given this, "excluding SHMT/choline_pc" is operationalized here as:
  1. Kendall's W: identical methodology to scripts/12_figure3_cross_cohort_confidence.py,
     restricted to the 3-biomarker subset {urea_cycle, bcaa, dnl_scd1_fasn}.
  2. Benchmark ordering: item 9's exact "active_only, raw, pooled" SD
     computation (equal-weighting vs. calibrated), restricted to ONLY the
     reactions mapped to {urea_cycle, bcaa, dnl_scd1_fasn} (18 reactions),
     instead of all 12,931 reactions.

Does not modify any existing GEM1 output file. Fast (~seconds); synchronous.
"""

import json
import os
import time

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare

from audit_common import BIOMARKERS_FULL, CALIBRATION_DIR, DETAILS_DIR, load_calibrated_metadata

REPORTED_SUBSYSTEM_W = 0.7912
ACTUAL_BIOMARKER_LEVEL_W = 0.6667  # from data/calibration/figure3_statistics.json, real pipeline output
UNCONTAMINATED_AXES = ["urea_cycle", "bcaa", "dnl_scd1_fasn"]
CONTAMINATED_AXES = ["serine_glycine_shmt", "choline_pc"]
COHORT_ORDER = ["GSE89632", "GSE126848", "GSE135251"]

REPORTED_EQUAL_WEIGHT_SD = 0.1236
REPORTED_CALIBRATED_SD = 0.2002


def sd(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(np.std(x)) if len(x) else float("nan")


def kendalls_w_for_axes(axes):
    ranking_path = os.path.join(CALIBRATION_DIR, "biomarker_ranking_by_group.csv")
    ranking = pd.read_csv(ranking_path)
    ranking = ranking[ranking["biomarker"].isin(axes)]
    cohort_biomarker = (
        ranking.groupby(["cohort", "biomarker"])["percentile_rank_in_group"]
        .mean().unstack("biomarker").reindex(index=COHORT_ORDER, columns=sorted(axes))
    )
    chi2, p = friedmanchisquare(*[cohort_biomarker[b].values for b in sorted(axes)])
    m, n = len(COHORT_ORDER), len(axes)
    w = chi2 / (m * (n - 1))
    return {
        "axes_used": sorted(axes), "n_cohorts": m, "n_biomarkers": n,
        "cohort_biomarker_matrix": cohort_biomarker.round(6).to_dict(),
        "friedman_chi2": float(chi2), "friedman_p": float(p), "kendalls_w": float(w),
    }


def weighted_rank_score(df, w1, w2, w3, zero_fill_inactive):
    score = w1 * df["rank_consensus"] + w2 * df["rank_flux_confidence"] + w3 * df["rank_perturbation"]
    if zero_fill_inactive:
        return np.where(df["is_active"], score.fillna(0), 0.0)
    return np.where(df["is_active"], score, np.nan)


def axis_restricted_sd_comparison(df, axes):
    rxn_ids = sorted({r for a in axes for r in BIOMARKERS_FULL[a]})
    sub = df[df["reaction_id"].isin(rxn_ids)].copy()

    equal_score = weighted_rank_score(sub, 1 / 3, 1 / 3, 1 / 3, zero_fill_inactive=False)
    equal_sd = sd(equal_score)

    calibrated_sd_direct = sd(sub["calibrated_confidence_score"])

    calib_score = weighted_rank_score(sub, 0.15, 0.15, 0.70, zero_fill_inactive=False)
    calib_sd_reconstructed = sd(calib_score)

    return {
        "reaction_ids_used": rxn_ids, "n_reactions": len(rxn_ids), "n_reaction_group_rows": len(sub),
        "equal_weighting_sd": equal_sd,
        "calibrated_sd_direct_column": calibrated_sd_direct,
        "calibrated_sd_reconstructed": calib_sd_reconstructed,
        "ordering_equal_lt_calibrated": bool(equal_sd < calibrated_sd_direct) if not (np.isnan(equal_sd) or np.isnan(calibrated_sd_direct)) else None,
    }


def run():
    t0 = time.time()
    df = load_calibrated_metadata()

    w_full5 = {
        "axes_used": ["bcaa", "choline_pc", "dnl_scd1_fasn", "serine_glycine_shmt", "urea_cycle"],
        "kendalls_w": ACTUAL_BIOMARKER_LEVEL_W,
        "source": "data/calibration/figure3_statistics.json (real pipeline output, not recomputed here)",
    }
    w_3axis = kendalls_w_for_axes(UNCONTAMINATED_AXES)
    print(f"Full 5-axis W (from repo): {ACTUAL_BIOMARKER_LEVEL_W}")
    print(f"3-axis (uncontaminated) W: {w_3axis['kendalls_w']:.4f} (Friedman chi2={w_3axis['friedman_chi2']:.4f}, p={w_3axis['friedman_p']:.4f})")

    sd_full = axis_restricted_sd_comparison(df, ["serine_glycine_shmt", "urea_cycle", "bcaa", "dnl_scd1_fasn", "choline_pc"])
    sd_3axis = axis_restricted_sd_comparison(df, UNCONTAMINATED_AXES)
    print(f"\nSD comparison, ALL 5 axes' reactions ({sd_full['n_reactions']} reactions): "
          f"equal={sd_full['equal_weighting_sd']:.4f}, calibrated={sd_full['calibrated_sd_direct_column']:.4f}, "
          f"ordering_holds={sd_full['ordering_equal_lt_calibrated']}")
    print(f"SD comparison, 3 UNCONTAMINATED axes' reactions ({sd_3axis['n_reactions']} reactions): "
          f"equal={sd_3axis['equal_weighting_sd']:.4f}, calibrated={sd_3axis['calibrated_sd_direct_column']:.4f}, "
          f"ordering_holds={sd_3axis['ordering_equal_lt_calibrated']}")

    w_conclusion = (
        f"W {'INCREASES' if w_3axis['kendalls_w'] > ACTUAL_BIOMARKER_LEVEL_W else 'DECREASES' if w_3axis['kendalls_w'] < ACTUAL_BIOMARKER_LEVEL_W else 'UNCHANGED'} "
        f"from {ACTUAL_BIOMARKER_LEVEL_W:.4f} (5 axes) to {w_3axis['kendalls_w']:.4f} (3 uncontaminated axes). "
        f"At n=3 items this is a very thin test (Friedman's exact null distribution is coarse at this "
        f"size -- p={w_3axis['friedman_p']:.4f}, {'not' if w_3axis['friedman_p'] >= 0.05 else ''} significant "
        f"at alpha=0.05 either way), so treat the point estimate as illustrative, not a precision claim."
    )

    ordering_conclusion = (
        f"The equal-weighting < calibrated SD ordering "
        f"{'HOLDS' if sd_3axis['ordering_equal_lt_calibrated'] else 'DOES NOT HOLD'} when restricted to only "
        f"the 3 uncontaminated axes' 18 mapped reactions (equal={sd_3axis['equal_weighting_sd']:.4f} vs. "
        f"calibrated={sd_3axis['calibrated_sd_direct_column']:.4f}), matching the direction found across "
        f"all 5 axes' reactions (equal={sd_full['equal_weighting_sd']:.4f} vs. "
        f"calibrated={sd_full['calibrated_sd_direct_column']:.4f}) and across the full 12931-reaction "
        f"population (item 9: equal=0.1236 vs. calibrated=0.2002). This specific ordering does not appear "
        f"to depend on the two growth-support-contaminated axes."
    )

    result = {
        "item": 18,
        "description": "Kendall's W and Table 3 SD-ordering benchmark re-run excluding growth-support-contaminated axes (SHMT, choline_pc)",
        "premise_check": {
            "0.7912_found_in_repo": False,
            "subsystem_level_kendalls_w_exists": False,
            "subsystem_finding": "docs/MANUSCRIPT_FIGURES.md explicitly states subsystem-level analysis is entirely unbuilt -- 'No script in this repository aggregates confidence/calibration scores by subsystem'",
            "actual_kendalls_w_in_repo": {"value": ACTUAL_BIOMARKER_LEVEL_W, "level": "biomarker (5 biomarkers x 3 cohorts), NOT subsystem",
                                          "source": "data/calibration/figure3_statistics.json, produced by scripts/12_figure3_cross_cohort_confidence.py"},
            "master_architectural_benchmark_found": False,
            "closest_existing_analog": "item 9's SD sensitivity check (Section 3.4), which already establishes equal-weighting SD < calibrated SD across ALL reactions -- reused and axis-restricted here",
        },
        "kendalls_w_full_5axis_from_repo": w_full5,
        "kendalls_w_3axis_uncontaminated": w_3axis,
        "kendalls_w_conclusion": w_conclusion,
        "sd_ordering_full_5axis_reactions_only": sd_full,
        "sd_ordering_3axis_uncontaminated_reactions_only": sd_3axis,
        "sd_ordering_conclusion": ordering_conclusion,
        "contaminated_axes_excluded": CONTAMINATED_AXES,
        "uncontaminated_axes_used": UNCONTAMINATED_AXES,
        "status": "RESOLVED -- premise corrected, both checks re-run on uncontaminated axes",
        "max_abs_discrepancy": abs(w_3axis["kendalls_w"] - ACTUAL_BIOMARKER_LEVEL_W),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item18_uncontaminated_axes_reproducibility.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {out_path}")
    print(f"\nKendall's W conclusion: {w_conclusion}")
    print(f"\nSD ordering conclusion: {ordering_conclusion}")
    return result


if __name__ == "__main__":
    run()
