"""
GEM1 - Step 12: Figure 5 - Calibration and Validation.

Per docs/MANUSCRIPT_FIGURES.md Figure 5.

Panel A: schematic of the nested CV structure (LOCO / LOBO, outer / inner
folds), reusing the same box style as Figure 1 Panel B for visual
continuity, per the spec.
Panel B: LOCO/LOBO mean scores at floor=0.15/0.05/0.0, with error bars
showing between-fold SD, illustrating the "narrow performance range"
relative to ordinary fold-to-fold variability.
Panel C: final production weight composition (0.15/0.15/0.70), shown as a
pie breakdown annotated with the interpretability rationale -- explicitly
NOT framed as the empirical optimum (it measurably isn't, per Panel B / the
floor=0.05 numbers).

Statistical analysis -- same situation as Figure 4: the spec states these
numbers were "computed fresh this session ... the Wilcoxon test itself
computed ad hoc ... not yet a saved script." This script is that permanent,
reproducible computation, run directly against
data/calibration/floor_sensitivity_results.json (floor=0.15 vs floor=0.05,
paired by matched fold identity, n=8, no ties this time).

RESOLVED DISCREPANCY (2026-08-05): this computation reproduces the spec's
effect-size numbers exactly (mean paired difference 0.0036, SD 0.0059,
between-fold SD 0.037 LOCO / 0.123 LOBO) but originally did not reproduce
the spec's stated Wilcoxon p-value of 0.1797 -- scipy.stats.wilcoxon on
these exact 8 matched pairs (W=8.0, matching the spec's W) gives p=0.1953
under its 'auto'/'exact' methods. Investigated before correcting: neither
scripts/09_calibration.py nor scripts/_floor_sensitivity_check.py contains
any Wilcoxon call at all, so the original number was not produced by a
different script using different parameters -- like Figure 4's number, it
was computed ad hoc outside any saved script. Explicitly checked
scipy.stats.wilcoxon's 'exact', 'approx' (with and without continuity
correction) methods against these 8 pairs: exact/auto=0.1953,
approx-no-correction=0.1614, approx-with-correction=0.1834 -- none
reproduce 0.1797, ruling out "different standard method" as the
explanation. Concluded a one-off arithmetic slip, same class as Figure 4's
effect-size error. The spec has been corrected to p=0.1953
(docs/MANUSCRIPT_FIGURES.md Figure 5, updated 2026-08-05); the qualitative
conclusion ("not statistically significant, directly supporting the
'statistically indistinguishable' claim") is unchanged, since both the old
and corrected p-value clear "not significant at alpha=0.05" identically.
This script is now the permanent, saved source of record, and the project
standard going forward is scipy.stats.wilcoxon(..., method='auto') (the
default), which resolves to the exact distribution here (n=8, no ties/
zeros).
"""

import os
import json
import numpy as np
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CAL_DIR = os.path.join(BASE_DIR, "data", "calibration")
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(CAL_DIR, "floor_sensitivity_results.json")) as f:
    floor_data = json.load(f)

r15 = floor_data["results_by_floor"]["0.15"]
r05 = floor_data["results_by_floor"]["0.05"]
r00 = floor_data["results_by_floor"]["0.0"]

def fold_scores(block):
    loco = [f["outer_score"] for f in block["loco_folds"]]
    lobo = [f["outer_score"] for f in block["lobo_folds"]]
    return loco, lobo

loco15, lobo15 = fold_scores(r15)
loco05, lobo05 = fold_scores(r05)
loco00, lobo00 = fold_scores(r00)

all15 = np.array(loco15 + lobo15)
all05 = np.array(loco05 + lobo05)
diffs = all05 - all15

w_stat, w_p = wilcoxon(all15, all05)
mean_diff = float(np.mean(diffs))
sd_diff = float(np.std(diffs, ddof=1))
loco_sd_15 = float(np.std(loco15, ddof=1))
lobo_sd_15 = float(np.std(lobo15, ddof=1))

print("=== Figure 5 statistics (computed fresh from floor_sensitivity_results.json) ===")
print(f"floor=0.15 scores (LOCO+LOBO, n=8): {np.round(all15, 4)}")
print(f"floor=0.05 scores (LOCO+LOBO, n=8): {np.round(all05, 4)}")
print(f"Wilcoxon: W={w_stat}, p={w_p}  (spec corrected 2026-08-05 to match this value; method=scipy default 'auto')")
print(f"Mean paired diff (0.05-0.15): {mean_diff:.4f}  (spec: 0.0036, matches)")
print(f"SD paired diff: {sd_diff:.4f}  (spec: 0.0059, matches)")
print(f"Between-fold SD @floor=0.15: LOCO={loco_sd_15:.4f} (spec 0.037), LOBO={lobo_sd_15:.4f} (spec 0.123)")

stats_out = {
    "wilcoxon_W": float(w_stat),
    "wilcoxon_p": float(w_p),
    "wilcoxon_method": "auto (scipy default; resolves to exact for n=8, no ties/zeros)",
    "mean_paired_difference_005_minus_015": mean_diff,
    "sd_paired_difference": sd_diff,
    "between_fold_sd_loco_floor015": loco_sd_15,
    "between_fold_sd_lobo_floor015": lobo_sd_15,
    "spec_correction_note": (
        "docs/MANUSCRIPT_FIGURES.md Figure 5 originally stated Wilcoxon p=0.1797 for this "
        f"comparison. This script's computation from the same source data gives p={w_p:.4f} "
        f"(W={w_stat} matched the original exactly, as did every effect-size number). Traced "
        "2026-08-05: neither 09_calibration.py nor _floor_sensitivity_check.py contains a "
        "Wilcoxon call, so the original was computed ad hoc outside any script; scipy's "
        "'exact'/'approx' (+/- continuity correction) methods were all checked and none "
        "reproduce 0.1797. Concluded a one-off arithmetic slip, not a methodology difference. "
        "The spec has been corrected to match this script's output; the qualitative conclusion "
        "('not statistically significant') is unchanged at alpha=0.05 either way."
    ),
}
with open(os.path.join(CAL_DIR, "figure5_statistics.json"), "w") as f:
    json.dump(stats_out, f, indent=2)

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], hspace=0.4, wspace=0.3)

# --- Panel A: nested CV schematic (same box style as Figure 1 Panel B) ---
axA = fig.add_subplot(gs[0, :])
axA.set_xlim(0, 10)
axA.set_ylim(0, 3.4)
axA.axis("off")
axA.set_title("A", loc="left", fontsize=15, fontweight="bold")


def draw_box(ax, xy, w, h, text, facecolor, fontsize=8.5, dashed=False, textcolor="white"):
    box = FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.02,rounding_size=0.02",
                          linewidth=1.3, edgecolor="black", facecolor=facecolor,
                          linestyle="dashed" if dashed else "solid", zorder=2)
    ax.add_patch(box)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center",
             fontsize=fontsize, color=textcolor, zorder=3)


def draw_arrow(ax, start, end, color="black"):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color=color, lw=1.4, zorder=1))


# LOCO row
draw_box(axA, (0.2, 2.0), 1.7, 1.0, "LOCO outer fold\n(1 of 3 cohorts\nheld out)", "#3b6ea5")
draw_arrow(axA, (1.9, 2.5), (2.5, 2.5))
draw_box(axA, (2.5, 2.0), 2.0, 1.0, "Inner LOCO\n(2 folds over\nremaining cohorts)\nselects weights", "#8a8a8a")
draw_arrow(axA, (4.5, 2.5), (5.1, 2.5))
draw_box(axA, (5.1, 2.0), 2.1, 1.0, "Apply fixed weights\nONCE to held-out\ncohort", "#5f9e6e")
draw_arrow(axA, (7.2, 2.5), (7.8, 2.5))
draw_box(axA, (7.8, 2.0), 1.9, 1.0, "outer_held_out_score\n(generalization\nestimate)", "#5f9e6e")

# LOBO row
draw_box(axA, (0.2, 0.4), 1.7, 1.0, "LOBO outer fold\n(1 of 5 biomarkers\nheld out)", "#3b6ea5")
draw_arrow(axA, (1.9, 0.9), (2.5, 0.9))
draw_box(axA, (2.5, 0.4), 2.0, 1.0, "Inner LOBO\n(4 folds over\nremaining biomarkers)\nselects weights", "#8a8a8a")
draw_arrow(axA, (4.5, 0.9), (5.1, 0.9))
draw_box(axA, (5.1, 0.4), 2.1, 1.0, "Apply fixed weights\nONCE to held-out\nbiomarker", "#5f9e6e")
draw_arrow(axA, (7.2, 0.9), (7.8, 0.9))
draw_box(axA, (7.8, 0.4), 1.9, 1.0, "outer_held_out_score\n(generalization\nestimate)", "#5f9e6e")

axA.text(0.05, 3.15, "LOCO (3 outer folds)", fontsize=9, fontweight="bold")
axA.text(0.05, 1.55, "LOBO (5 outer folds)", fontsize=9, fontweight="bold")

# --- Panel B: floor sensitivity, mean +/- SD ---
axB = fig.add_subplot(gs[1, 0])
floors = ["0.15", "0.05", "0.0"]
loco_means = [np.mean(loco15), np.mean(loco05), np.mean(loco00)]
loco_sds = [np.std(loco15, ddof=1), np.std(loco05, ddof=1), np.std(loco00, ddof=1)]
lobo_means = [np.mean(lobo15), np.mean(lobo05), np.mean(lobo00)]
lobo_sds = [np.std(lobo15, ddof=1), np.std(lobo05, ddof=1), np.std(lobo00, ddof=1)]

x = np.arange(len(floors))
width = 0.35
axB.bar(x - width / 2, loco_means, width, yerr=loco_sds, capsize=4, label="LOCO", color="#3b6ea5")
axB.bar(x + width / 2, lobo_means, width, yerr=lobo_sds, capsize=4, label="LOBO", color="#c0783c")
axB.set_xticks(x)
axB.set_xticklabels([f"floor={f}" for f in floors])
axB.set_ylabel("Mean outer-held-out score\n(error bars = between-fold SD)")
axB.set_ylim(0, 1.55)
axB.set_title("B", loc="left", fontsize=15, fontweight="bold")
axB.legend(loc="upper left")
axB.text(0.98, 0.99,
          f"Wilcoxon (0.15 vs 0.05, n=8 matched folds): W={w_stat}, p={w_p:.3f}\n"
          f"paired-diff SD ({sd_diff:.4f}) is {loco_sd_15/sd_diff:.0f}-{lobo_sd_15/sd_diff:.0f}x "
          "smaller than between-fold SD",
          transform=axB.transAxes, ha="right", va="top", fontsize=7.8,
          bbox=dict(boxstyle="round", facecolor="#f5f5f5", edgecolor="#999999"))

# --- Panel C: production weight composition ---
axC = fig.add_subplot(gs[1, 1])
weights = [0.15, 0.15, 0.70]
labels = ["w1: consensus\n(0.15)", "w2: flux confidence\n(0.15)", "w3: perturbation\nrobustness (0.70)"]
colors = ["#3b6ea5", "#c0783c", "#5f9e6e"]
axC.pie(weights, labels=labels, colors=colors, autopct="%1.0f%%", startangle=90,
        wedgeprops=dict(edgecolor="white", linewidth=1.5), textprops=dict(fontsize=9))
axC.set_title("C", loc="left", fontsize=15, fontweight="bold")
axC.text(0, -1.35,
          "Retained for decomposability/interpretability,\nNOT because it is the numerical optimum\n"
          "(floor=0.05/0.0 select (0.5, 0.05, 0.45), score +0.0036 higher)",
          ha="center", va="top", fontsize=8, style="italic")

caption = (
    "Figure 5. Nested cross-validation for empirical weight calibration. (A) LOCO/LOBO nested "
    f"structure. (B) Performance is a near-tie across weight-floor settings (Wilcoxon W={w_stat}, "
    f"p={w_p:.3f}, paired by fold), differing by roughly {loco_sd_15/sd_diff:.0f}-{lobo_sd_15/sd_diff:.0f}x "
    "less than ordinary fold-to-fold variability. (C) Production weights (0.15/0.15/0.70) were "
    "retained for decomposability, not because they are the numerical optimum -- a marginally "
    "higher-scoring, consensus-dominant alternative exists and is disclosed. See "
    "figure5_statistics.json (this corrected an earlier spec value of p=0.1797, traced to an "
    "unrecoverable ad hoc arithmetic slip; same qualitative conclusion either way)."
)
fig.text(0.5, 0.0, caption, ha="center", va="top", fontsize=8.2, wrap=True)

out_path = os.path.join(OUT_DIR, "figure5_calibration_validation.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Wrote {out_path}")
