"""
Item 10 (priority 4) -- Section 3.7 per-group BCAA severity means.

Manuscript claims: mean percentile ranks for the BCAA axis across
healthy/steatosis/NASH = 0.5066, 0.9044, 0.8430.

Independently recomputed directly from data/calibration/biomarker_ranking_by_group.csv
(an existing raw pipeline output, not a scoring/calibration *module* -- this file is
data, not code; no scoring/calibration script is imported here), by taking the
`percentile_rank_in_group` column for biomarker=="bcaa", grouped by disease severity
(healthy/steatosis/NASH) and averaged across the 3 cohorts that have that group label.
GSE126848's 4th group ("obese_no_NAFLD") is excluded from the 3-point severity ladder,
since it isn't part of the healthy->steatosis->NASH progression the other two cohorts
share -- same population choice implied by the manuscript's "healthy/steatosis/NASH"
framing (3 categories, not 4).

Fast (~seconds); run synchronously, no detached process needed.
"""

import json
import os
import time

import numpy as np
import pandas as pd

from audit_common import CALIBRATION_DIR, DETAILS_DIR

BIOMARKER_RANKING_CSV = os.path.join(CALIBRATION_DIR, "biomarker_ranking_by_group.csv")

REPORTED_MEANS = {"healthy": 0.5066, "steatosis": 0.9044, "NASH": 0.8430}
SEVERITY_GROUPS = ["healthy", "steatosis", "NASH"]
TOL = 0.001


def run():
    t0 = time.time()
    df = pd.read_csv(BIOMARKER_RANKING_CSV)
    bcaa = df[df["biomarker"] == "bcaa"]

    per_group_detail = {}
    independent_means = {}
    for sev in SEVERITY_GROUPS:
        sub = bcaa[bcaa["group"] == sev]
        per_group_detail[sev] = sub[["cohort", "percentile_rank_in_group"]].to_dict("records")
        independent_means[sev] = float(sub["percentile_rank_in_group"].mean())

    discrepancies = {sev: abs(independent_means[sev] - REPORTED_MEANS[sev]) for sev in SEVERITY_GROUPS}
    max_abs_discrepancy = max(discrepancies.values())

    result = {
        "item": 10,
        "description": "Section 3.7 per-group BCAA severity means, independently recomputed from biomarker_ranking_by_group.csv",
        "population": "biomarker=='bcaa' rows, grouped by disease severity (healthy/steatosis/NASH), "
                       "averaged across the 3 cohorts sharing that group label; GSE126848's "
                       "'obese_no_NAFLD' 4th group excluded (not part of the 3-point severity ladder)",
        "per_group_detail": per_group_detail,
        "independent_means": independent_means,
        "reported_means": REPORTED_MEANS,
        "discrepancies": discrepancies,
        "max_abs_discrepancy": max_abs_discrepancy,
        "pass_threshold": TOL,
        "status": "PASS" if max_abs_discrepancy < TOL else "FAIL",
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item10_bcaa_severity_means.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"Independent means: {independent_means}")
    print(f"Status: {result['status']}, max_abs_discrepancy={max_abs_discrepancy:.6f}")
    return result


if __name__ == "__main__":
    run()
