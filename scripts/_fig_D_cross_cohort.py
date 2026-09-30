"""
GEM1 - Figure D: cross-cohort matched-disease-state reproducibility. Two
panels: Spearman correlation restricted to the "either cohort detected"
subset (addresses the large tied/null mass directly, per instruction), and
Jaccard overlap of detected reaction sets, for all 9 same-state cohort pairs
(healthy x3 cohorts, steatosis x3, NASH x3; obese_no_NAFLD has no cross-
cohort match and is excluded here by construction, not by omission).
Reads only 04_cross_cohort_reproducibility.csv (already computed, read-only).
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
TABLES = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "tables")
FIGS = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "figures")
os.makedirs(FIGS, exist_ok=True)

STATE_COLORS = {"healthy": "#0072B2", "steatosis": "#E69F00", "NASH": "#D55E00"}

t4 = pd.read_csv(os.path.join(TABLES, "04_cross_cohort_reproducibility.csv"))
t4.to_csv(os.path.join(FIGS, "figD_source_data.csv"), index=False)
t4 = t4.sort_values(["state", "cohort_A", "cohort_B"]).reset_index(drop=True)
t4["pair_label"] = t4["cohort_A"] + "\nvs\n" + t4["cohort_B"]
colors = [STATE_COLORS[s] for s in t4["state"]]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
x = range(len(t4))

ax1.bar(x, t4["spearman_r_detected_subset"], color=colors)
ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_xticks(list(x))
ax1.set_xticklabels(t4["pair_label"], fontsize=7.5)
ax1.set_ylabel("Spearman r (detected-in-either subset only)")
ax1.set_title("D1. Rank concordance, detected subset\n(excludes the joint-zero tied mass)")

ax2.bar(x, t4["jaccard_detected_overlap"], color=colors)
for i, (v, n) in enumerate(zip(t4["jaccard_detected_overlap"], t4["n_either_detected"])):
    ax2.text(i, v + 0.01, f"n={n}", ha="center", fontsize=7)
ax2.set_xticks(list(x))
ax2.set_xticklabels(t4["pair_label"], fontsize=7.5)
ax2.set_ylim(0, max(t4["jaccard_detected_overlap"]) * 1.25)
ax2.set_ylabel("Jaccard overlap of detected reaction sets")
ax2.set_title("D2. Detected-set overlap (both/either detected)")

handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in STATE_COLORS.values()]
fig.legend(handles, STATE_COLORS.keys(), loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.04))
fig.suptitle("Figure D - Cross-cohort matched-state reproducibility (9 same-state cohort pairs)", y=1.1, fontsize=12)
fig.tight_layout()

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figD_cross_cohort.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("Wrote figD_cross_cohort.{pdf,svg,png} + figD_source_data.csv")
