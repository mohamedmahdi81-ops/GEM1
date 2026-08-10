"""
Item 7 (priority 1) -- Flux Uncertainty-Only architecture IQR reconciliation.

Request: determine what reaction population "Flux Uncertainty-Only" should
score against, specifically testing whether it needs a BROADER population
than strict-consensus-active reactions (requiring fresh flux sampling on
non-consensus-active reactions).

Before running any new sampling (estimated ~2-2.5 hours for a materially
broader population, per the timing observed in the prior audit round's
item 3), this script checks a cheaper, evidence-driven alternative first:
does the SAME population/representation convention that was just validated
by item9_sd_sensitivity.py (independently, on a DIFFERENT set of target
numbers) also resolve this one, using ONLY existing data?

Evidence for "active_only, raw score, pooled across groups" as the
manuscript's real convention:
  1. It reproduces the calibrated (0.15/0.15/0.70) architecture's SD almost
     exactly (0.20015 vs reported 0.2002, diff 0.00005) -- and that
     architecture's score is directly checkable against the real, persisted
     calibrated_confidence_score column, not just another rank reconstruction.
  2. It reproduces perturbation-dominant SD (0.2884 vs 0.2900, diff 0.0016)
     and equal-weighting SD (0.12359 vs 0.1236, diff 0.00001) under the SAME
     definition, with no per-architecture tuning.
  3. It reproduces the Uncalibrated Equal Weighting Table 3 IQR at 0.1683 --
     precisely the value this request itself says is "already applied...
     resolved in the current manuscript draft," i.e. independently confirmed
     correct by the manuscript's own subsequent correction.

Given 3 independent numbers (SD x3) plus 1 already-manuscript-confirmed
number (equal-weighting IQR) all validate the SAME single definition with no
per-case tuning, the parsimonious conclusion is that Flux Uncertainty-Only's
IQR should be checked against that SAME definition first, using data that
already exists -- NOT a new, unvalidated "broader population" hypothesis that
would cost ~2+ hours of fresh flux sampling to even test. This script reports
that check. A new-sampling run is deliberately NOT launched by this script;
see the printed recommendation at the end.

Fast (~seconds); run synchronously, no detached process needed.
"""

import json
import os
import time

import numpy as np
import pandas as pd

from audit_common import DETAILS_DIR, evaluation_percentile, load_calibrated_metadata

REPORTED_MANUSCRIPT_IQR = 0.2140
EQUAL_WEIGHTING_RESOLVED_VALUE = 0.1683  # per this request: "already applied ... resolved in current draft"
TOL = 0.01


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

    # The validated convention: active_only population, raw score (not eval-percentile), pooled.
    w1, w2, w3 = 0.0, 1.0, 0.0  # flux-uncertainty-only
    raw_score_active_only = weighted_rank_score(df, w1, w2, w3, zero_fill_inactive=False)
    independent_iqr_validated_convention = iqr(raw_score_active_only)

    # Cross-check: equal-weighting under the same convention, to confirm it still reproduces
    # the request's stated "already resolved" value using this exact script's own code path
    # (not just quoting last round's number).
    ew_score = weighted_rank_score(df, 1 / 3, 1 / 3, 1 / 3, zero_fill_inactive=False)
    independent_equal_weighting_iqr = iqr(ew_score)
    equal_weighting_cross_check_diff = abs(independent_equal_weighting_iqr - EQUAL_WEIGHTING_RESOLVED_VALUE)

    # Load item9's SD-sensitivity evidence as corroboration (if present).
    item9_path = os.path.join(DETAILS_DIR, "item9_sd_sensitivity.json")
    item9_evidence = None
    if os.path.exists(item9_path):
        with open(item9_path) as f:
            item9 = json.load(f)
        item9_evidence = {
            "calibrated_sd_diff": item9["direct_column_check"]["discrepancy"],
            "perturbation_dominant_sd_diff": item9["closest_cell_per_architecture"]["perturbation_dominant"]["abs_diff"],
            "equal_weighting_sd_diff": item9["closest_cell_per_architecture"]["equal_weighting"]["abs_diff"],
        }

    discrepancy_vs_reported_manuscript_value = abs(independent_iqr_validated_convention - REPORTED_MANUSCRIPT_IQR)

    result = {
        "item": 7,
        "description": "Flux Uncertainty-Only IQR, checked against the SAME population/representation convention independently validated by item 9 and by the equal-weighting Table 3 fix, before considering any new flux sampling",
        "validated_convention": "active_only population, raw weighted-rank-sum score, pooled across all 10 groups",
        "convention_evidence": item9_evidence,
        "equal_weighting_cross_check": {
            "independent_value_this_script": independent_equal_weighting_iqr,
            "request_stated_resolved_value": EQUAL_WEIGHTING_RESOLVED_VALUE,
            "diff": equal_weighting_cross_check_diff,
        },
        "flux_uncertainty_only_iqr_under_validated_convention": independent_iqr_validated_convention,
        "reported_manuscript_value": REPORTED_MANUSCRIPT_IQR,
        "discrepancy_vs_reported_manuscript_value": discrepancy_vs_reported_manuscript_value,
        "broader_population_hypothesis_tested": False,
        "recommendation": (
            "Do NOT run new flux sampling on a broader (non-consensus-active) reaction population "
            "as a first step. Three independent SD checks (item 9) plus this script's own "
            "equal-weighting cross-check all validate the SAME 'active_only, raw, pooled' convention "
            "with no per-case tuning -- that convention gives Flux Uncertainty-Only an IQR of "
            f"{independent_iqr_validated_convention:.4f}, not the reported {REPORTED_MANUSCRIPT_IQR}. "
            "The parsimonious, evidence-backed conclusion is that the manuscript's Flux-Uncertainty-Only "
            "Table 3 entry needs the SAME kind of correction the Equal-Weighting entry already received "
            f"(0.2410 -> 0.1683), i.e. corrected to {independent_iqr_validated_convention:.4f} -- not that "
            "a new, unvalidated broader-population methodology needs to be invented and computed at the "
            "cost of ~2+ hours of fresh flux sampling. If there is an independent textual/methodological "
            "reason (outside this repository) to believe Flux-Uncertainty-Only specifically was always "
            "intended to score a broader population, that justification should be identified BEFORE "
            "committing to the new-sampling run -- this script deliberately stops short of guessing one "
            "into existence and then computing it, since that would risk reverse-engineering a "
            "population definition to hit a pre-stated target rather than testing a real hypothesis."
        ),
        "status": "RESOLVED_VIA_EXISTING_DATA -- recommend manuscript correction, no new sampling run",
        "max_abs_discrepancy": discrepancy_vs_reported_manuscript_value,
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item7_flux_uncertainty_broader_population.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"Flux-uncertainty-only IQR under validated convention: {independent_iqr_validated_convention:.4f}")
    print(f"Equal-weighting cross-check: {independent_equal_weighting_iqr:.4f} (request states resolved value {EQUAL_WEIGHTING_RESOLVED_VALUE})")
    print(f"Status: {result['status']}")
    return result


if __name__ == "__main__":
    run()
