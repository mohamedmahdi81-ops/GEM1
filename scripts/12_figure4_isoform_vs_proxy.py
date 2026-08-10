"""
GEM1 - Step 12: Figure 4 - Full-Isoform vs. Single-Proxy Validation.

Per docs/MANUSCRIPT_FIGURES.md Figure 4.

Panel A: grouped bar chart, LOCO scores (3 cohort folds), full-isoform vs.
single-proxy.
Panel B: grouped bar chart, LOBO scores (5 biomarker folds), full-isoform
vs. single-proxy, with the choline_pc single-proxy failure (~0.0001)
labeled explicitly per the spec's "not hidden" instruction.

Statistical analysis -- the spec states these were "computed fresh this
session ... not yet a saved script -- worth adding as a permanent check if
this figure is finalized." This script IS that permanent check: every
number below is computed here, directly from
data/calibration/calibration_weights_and_validation.json's per-fold arrays,
not hardcoded from the spec text.

RESOLVED DISCREPANCY (2026-08-05): this computation originally surfaced a
mismatch against the frozen spec's stated "Effect size: mean paired
difference = +0.176" -- this script's direct computation gives +0.336
(median +0.297). The spec has since been corrected to +0.336
(docs/MANUSCRIPT_FIGURES.md Figure 4, updated 2026-08-05) after tracing the
original number: this is not a git repository and no earlier draft of the
spec exists on disk, so the +0.176 derivation itself is unrecoverable, but
the Wilcoxon W/p and sign-test p reported alongside it all reproduce
exactly from the same fold data (ruling out stale/different underlying
data), and multiple alternative interpretations of "mean paired difference"
were tried and none reproduce +0.176 either -- by elimination, a one-off
arithmetic/transcription slip at ad hoc computation time, not a methodology
difference. This script is now the permanent, saved source of record for
these numbers, superseding that ad hoc computation.
"""

import os
import json
import numpy as np
from scipy.stats import wilcoxon, binomtest
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CAL_DIR = os.path.join(BASE_DIR, "data", "calibration")
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(CAL_DIR, "calibration_weights_and_validation.json")) as f:
    calib = json.load(f)

full = calib["results_by_biomarker_set"]["full_isoform"]
single = calib["results_by_biomarker_set"]["single_proxy"]

loco_cohorts = [fold["outer_test_cohort"] for fold in full["loco"]["folds"]]
loco_full = [fold["outer_held_out_score"] for fold in full["loco"]["folds"]]
loco_single = [fold["outer_held_out_score"] for fold in single["loco"]["folds"]]

lobo_biomarkers = [fold["outer_test_biomarker"] for fold in full["lobo"]["folds"]]
lobo_full = [fold["outer_held_out_score"] for fold in full["lobo"]["folds"]]
lobo_single = [fold["outer_held_out_score"] for fold in single["lobo"]["folds"]]

all_full = np.array(loco_full + lobo_full)
all_single = np.array(loco_single + lobo_single)
all_labels = loco_cohorts + lobo_biomarkers
diffs = all_full - all_single

w_stat, w_p = wilcoxon(all_full, all_single)
nonzero = diffs[diffs != 0]
n_nonzero = len(nonzero)
n_pos = int(np.sum(nonzero > 0))
sign_result = binomtest(n_pos, n_nonzero, 0.5, alternative="greater")

mean_diff = float(np.mean(diffs))
median_diff = float(np.median(diffs))

print("=== Figure 4 statistics (computed fresh from calibration_weights_and_validation.json) ===")
print(f"Matched folds (n=8): {all_labels}")
print(f"full_isoform scores: {np.round(all_full, 4)}")
print(f"single_proxy scores: {np.round(all_single, 4)}")
print(f"diffs (full - single): {np.round(diffs, 4)}")
print(f"Wilcoxon signed-rank: W={w_stat}, p={w_p} (n_nonzero={n_nonzero})")
print(f"Sign test (one-sided, n={n_nonzero}, n_pos={n_pos}): p={sign_result.pvalue:.4f}")
print(f"Mean paired difference: {mean_diff:.4f}  (spec corrected 2026-08-05 to match this value)")
print(f"Median paired difference: {median_diff:.4f}  (spec corrected 2026-08-05 to match this value)")

stats_out = {
    "matched_fold_labels": all_labels,
    "full_isoform_scores": all_full.tolist(),
    "single_proxy_scores": all_single.tolist(),
    "wilcoxon_W": float(w_stat),
    "wilcoxon_p": float(w_p),
    "n_nonzero_pairs": n_nonzero,
    "sign_test_p_one_sided": float(sign_result.pvalue),
    "mean_paired_difference": mean_diff,
    "median_paired_difference": median_diff,
    "spec_correction_note": (
        "docs/MANUSCRIPT_FIGURES.md Figure 4 originally stated mean paired difference = +0.176 "
        "(median +0.31). This script's computation from the same source JSON gave "
        f"{mean_diff:.4f} (median {median_diff:.4f}) instead. Traced 2026-08-05: no git history "
        "or earlier draft exists to recover the original derivation, but the W/p/sign-test "
        "numbers reported alongside +0.176 all reproduce exactly, ruling out stale data; several "
        "alternative interpretations of 'mean paired difference' were tried and none reproduce "
        "+0.176. Concluded a one-off arithmetic/transcription slip, not a methodology "
        "difference. The spec has been corrected to match this script's output."
    ),
}
with open(os.path.join(CAL_DIR, "figure4_statistics.json"), "w") as f:
    json.dump(stats_out, f, indent=2)

# ---------------------------------------------------------------------------
# Panel A: LOCO grouped bars
# ---------------------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6))

x = np.arange(len(loco_cohorts))
width = 0.35
axA.bar(x - width / 2, loco_full, width, label="Full-isoform", color="#5f9e6e")
axA.bar(x + width / 2, loco_single, width, label="Single-proxy", color="#c0783c")
axA.set_xticks(x)
axA.set_xticklabels(loco_cohorts)
axA.set_ylabel("LOCO outer-held-out score")
axA.set_ylim(0, 1.32)
axA.axhline(np.mean(loco_full), color="#5f9e6e", linestyle="--", linewidth=1, alpha=0.6)
axA.axhline(np.mean(loco_single), color="#c0783c", linestyle="--", linewidth=1, alpha=0.6)
axA.set_title("A", loc="left", fontsize=15, fontweight="bold")
axA.legend(loc="upper left")
axA.text(0.98, 0.99,
          f"LOCO mean:\nfull={np.mean(loco_full):.3f}\nsingle={np.mean(loco_single):.3f}",
          transform=axA.transAxes, ha="right", va="top", fontsize=9)

# ---------------------------------------------------------------------------
# Panel B: LOBO grouped bars, choline_pc outlier callout
# ---------------------------------------------------------------------------
xb = np.arange(len(lobo_biomarkers))
axB.bar(xb - width / 2, lobo_full, width, label="Full-isoform", color="#5f9e6e")
axB.bar(xb + width / 2, lobo_single, width, label="Single-proxy", color="#c0783c")
axB.set_xticks(xb)
axB.set_xticklabels(lobo_biomarkers, rotation=20, ha="right")
axB.set_ylabel("LOBO outer-held-out score")
axB.set_ylim(0, 1.32)
axB.set_title("B", loc="left", fontsize=15, fontweight="bold")
axB.legend(loc="upper left")

choline_idx = lobo_biomarkers.index("choline_pc")
axB.annotate(
    f"single-proxy catastrophic\nfailure: {lobo_single[choline_idx]:.4f}",
    xy=(choline_idx + width / 2, lobo_single[choline_idx]),
    xytext=(choline_idx + width / 2 - 0.9, 0.45),
    fontsize=8.5, color="#8b3a2a",
    arrowprops=dict(arrowstyle="-|>", color="#8b3a2a"),
)
axB.text(0.98, 0.99,
          f"LOBO mean:\nfull={np.mean(lobo_full):.3f}\nsingle={np.mean(lobo_single):.3f}",
          transform=axB.transAxes, ha="right", va="top", fontsize=9)

caption = (
    "Figure 4. Held-out validation comparing full-isoform aggregation (MAX over all mapped "
    f"reactions per biomarker) against a single-best-reaction proxy. (A) LOCO: {np.mean(loco_full):.3f} "
    f"vs. {np.mean(loco_single):.3f}. (B) LOBO: {np.mean(lobo_full):.3f} vs. {np.mean(lobo_single):.3f}, "
    "per biomarker, highlighting a catastrophic single-proxy failure for choline/PC "
    f"({lobo_single[choline_idx]:.4f}) absent under full-isoform aggregation. Wilcoxon signed-rank "
    f"W={w_stat}, p={w_p:.5f} across 8 matched folds (2 exact ties dropped). Mean paired "
    f"difference = {mean_diff:+.3f} (see figure4_statistics.json; this corrected an earlier "
    "spec value of +0.176, traced to an unrecoverable ad hoc transcription slip)."
)
fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=8.3, wrap=True)

out_path = os.path.join(OUT_DIR, "figure4_isoform_vs_proxy.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Wrote {out_path}")
