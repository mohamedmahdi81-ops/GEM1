"""
GEM1 - Step 12: Figure 3 - Cross-Cohort Confidence Analysis.

Per docs/MANUSCRIPT_FIGURES.md Figure 3.

Panel A: distribution (violin) of calibrated_confidence_score among ACTIVE
reactions only (is_active == True), one violin per cohort, pooling that
cohort's groups. Restricting to active reactions matters: inactive
reactions are fixed at 0 by construction (Step 11's "Layer 1 is
foundational" rule), so including them would just show a huge spike at 0
that has nothing to do with cross-cohort consistency.

Panel B: biomarker x group grid (5 biomarkers x 10 groups), color-coded by
each biomarker's max-aggregated percentile rank in that group -- read
directly from data/calibration/biomarker_ranking_by_group.csv (Open Item 2
of the spec), not recomputed here.

Statistical analysis -- per the spec's own explicit methodological note:
a naive Kruskal-Wallis/ANOVA directly on calibrated_confidence_score across
cohorts is REJECTED as circular (rank-normalization forces every group's
active-reaction mean to ~0.5 by construction; CLAUDE.md's own
calibration_manifest.csv confirms this: all 10 groups sit at 0.5001). This
script does not compute that test.

- Primary statistic: nested LOCO cross-validation, already computed in Step
  11 (scripts/09_calibration.py) -- read here, not recomputed, from
  data/calibration/calibration_weights_and_validation.json
  (full_isoform.loco.mean_outer_held_out_score and per-fold range).
- Secondary statistic: the spec explicitly flags this as "not yet computed
  ... should be added at Step 12, not fabricated here" -- Kendall's
  coefficient of concordance (W) across the 3 cohorts on the 5 biomarkers'
  rank order, computed fresh in this script from
  biomarker_ranking_by_group.csv (mean percentile per cohort per biomarker,
  averaged over that cohort's own groups, then ranked). Implemented via
  scipy.stats.friedmanchisquare (m=3 raters/cohorts as blocks, n=5
  treatments/biomarkers), using the identity W = chi2 / (m*(n-1)).
  This is a genuinely new computation -- reported honestly regardless of
  outcome, not tuned to a pre-stated expectation (none was stated for this
  secondary test).
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CAL_DIR = os.path.join(BASE_DIR, "data", "calibration")
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

COHORT_ORDER = ["GSE89632", "GSE126848", "GSE135251"]
GROUP_ORDER = {
    "GSE89632": ["healthy", "steatosis", "NASH"],
    "GSE126848": ["healthy", "obese_no_NAFLD", "steatosis", "NASH"],
    "GSE135251": ["healthy", "steatosis", "NASH"],
}
BIOMARKER_ORDER = ["serine_glycine_shmt", "urea_cycle", "bcaa", "dnl_scd1_fasn", "choline_pc"]
BIOMARKER_LABELS = {
    "serine_glycine_shmt": "Serine/glycine\n(SHMT)",
    "urea_cycle": "Urea cycle",
    "bcaa": "BCAA\ncatabolism",
    "dnl_scd1_fasn": "De novo\nlipogenesis",
    "choline_pc": "Choline/PC",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
scores = pd.read_csv(os.path.join(CAL_DIR, "all_groups_calibrated.csv"),
                      usecols=["cohort", "group", "is_active", "calibrated_confidence_score"])
active = scores[scores["is_active"]].copy()

ranking = pd.read_csv(os.path.join(CAL_DIR, "biomarker_ranking_by_group.csv"))

with open(os.path.join(CAL_DIR, "calibration_weights_and_validation.json")) as f:
    calib_json = json.load(f)
loco = calib_json["results_by_biomarker_set"]["full_isoform"]["loco"]
loco_mean = loco["mean_outer_held_out_score"]
loco_fold_scores = [fold["outer_held_out_score"] for fold in loco["folds"]]
loco_range = (min(loco_fold_scores), max(loco_fold_scores))

# ---------------------------------------------------------------------------
# Secondary statistic: Kendall's W across cohorts on biomarker rank order
# ---------------------------------------------------------------------------
# Mean percentile per (cohort, biomarker), averaged over that cohort's groups.
cohort_biomarker = (
    ranking.groupby(["cohort", "biomarker"])["percentile_rank_in_group"]
    .mean()
    .unstack("biomarker")
    .reindex(index=COHORT_ORDER, columns=BIOMARKER_ORDER)
)
# friedmanchisquare(*treatments): each argument = one biomarker's values across the 3 cohort-blocks.
chi2, friedman_p = friedmanchisquare(*[cohort_biomarker[b].values for b in BIOMARKER_ORDER])
m_raters, n_items = len(COHORT_ORDER), len(BIOMARKER_ORDER)
kendalls_w = chi2 / (m_raters * (n_items - 1))

print("=== Figure 3 statistics ===")
print(f"Primary (Step 11 LOCO, read not recomputed): mean={loco_mean:.4f}, "
      f"range=({loco_range[0]:.4f}, {loco_range[1]:.4f})")
print(f"Secondary (Kendall's W, computed fresh this script): "
      f"W={kendalls_w:.4f}, Friedman chi2={chi2:.4f}, p={friedman_p:.4f} "
      f"(m={m_raters} cohorts, n={n_items} biomarkers)")
print("Cohort x biomarker mean-percentile matrix used for Kendall's W:")
print(cohort_biomarker.round(3).to_string())

stats_out = {
    "primary_loco_mean_full_isoform": loco_mean,
    "primary_loco_fold_range": list(loco_range),
    "secondary_kendalls_w": kendalls_w,
    "secondary_friedman_chi2": chi2,
    "secondary_friedman_p": friedman_p,
    "secondary_n_cohorts": m_raters,
    "secondary_n_biomarkers": n_items,
    "secondary_note": (
        "Kendall's W across cohorts on biomarker rank order, newly computed at Step 12 per "
        "the spec's explicit instruction (docs/MANUSCRIPT_FIGURES.md Figure 3 'Secondary "
        "statistical test' field). Not previously computed anywhere in the pipeline."
    ),
}
with open(os.path.join(CAL_DIR, "figure3_statistics.json"), "w") as f:
    json.dump(stats_out, f, indent=2)

# ---------------------------------------------------------------------------
# Panel A: violin of calibrated_confidence_score among active reactions, per cohort
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.3], hspace=0.4)

axA = fig.add_subplot(gs[0])
violin_data = [active.loc[active["cohort"] == c, "calibrated_confidence_score"].values for c in COHORT_ORDER]
parts = axA.violinplot(violin_data, showmeans=True, showmedians=True)
for pc in parts["bodies"]:
    pc.set_facecolor("#3b6ea5")
    pc.set_alpha(0.6)
axA.set_xticks(range(1, len(COHORT_ORDER) + 1))
axA.set_xticklabels([f"{c}\n(n={len(violin_data[i]):,})" for i, c in enumerate(COHORT_ORDER)])
axA.set_ylabel("calibrated_confidence_score\n(active reactions only)")
axA.set_title("A", loc="left", fontsize=15, fontweight="bold")
axA.annotate(
    f"Primary evidence for reproducibility: nested LOCO cross-validation\n"
    f"mean held-out score = {loco_mean:.3f} (range {loco_range[0]:.3f}-{loco_range[1]:.3f})\n"
    f"[not a test on this panel's raw distributions -- see caption]",
    xy=(0.98, 0.97), xycoords="axes fraction", ha="right", va="top", fontsize=8,
    bbox=dict(boxstyle="round", facecolor="#fff7e6", edgecolor="#c0783c"),
)

# ---------------------------------------------------------------------------
# Panel B: biomarker x group heatmap of percentile rank
# ---------------------------------------------------------------------------
axB = fig.add_subplot(gs[1])
group_cols = [(c, g) for c in COHORT_ORDER for g in GROUP_ORDER[c]]
col_labels = [f"{c}\n{g}" for c, g in group_cols]

grid = np.full((len(BIOMARKER_ORDER), len(group_cols)), np.nan)
for i, biomarker in enumerate(BIOMARKER_ORDER):
    for j, (c, g) in enumerate(group_cols):
        row = ranking[(ranking["biomarker"] == biomarker) & (ranking["cohort"] == c) & (ranking["group"] == g)]
        if len(row):
            grid[i, j] = row["percentile_rank_in_group"].values[0]

im = axB.imshow(grid, aspect="auto", cmap="viridis", vmin=0, vmax=1)
axB.set_xticks(range(len(group_cols)))
axB.set_xticklabels(col_labels, rotation=90, fontsize=7.5)
axB.set_yticks(range(len(BIOMARKER_ORDER)))
axB.set_yticklabels([BIOMARKER_LABELS[b] for b in BIOMARKER_ORDER], fontsize=8.5)
axB.set_title("B", loc="left", fontsize=15, fontweight="bold")
for i in range(len(BIOMARKER_ORDER)):
    for j in range(len(group_cols)):
        val = grid[i, j]
        if not np.isnan(val):
            axB.text(j, i, f"{val:.2f}", ha="center", va="center",
                      fontsize=6.5, color="white" if val < 0.6 else "black")
cbar = fig.colorbar(im, ax=axB, fraction=0.025, pad=0.01)
cbar.set_label("Max-aggregated percentile rank\n(within group, among all 12,931 reactions)", fontsize=8)

caption = (
    "Figure 3. Calibrated confidence scores across three independent NAFLD/MASLD cohorts "
    "(GSE89632, GSE126848, GSE135251). (A) Score distributions per cohort (active reactions "
    "only). (B) Ranking of the 5 literature-curated calibration biomarkers across all 10 disease "
    f"groups, supported quantitatively by nested leave-one-cohort-out cross-validation (mean "
    f"held-out score {loco_mean:.3f}). Secondary check: Kendall's W = {kendalls_w:.3f} "
    f"(Friedman chi2={chi2:.2f}, p={friedman_p:.3f}) for agreement among the 3 cohorts on "
    "biomarker rank order -- newly computed at Step 12, see figure3_statistics.json."
)
fig.text(0.5, 0.005, caption, ha="center", va="bottom", fontsize=8.3, wrap=True)

out_path = os.path.join(OUT_DIR, "figure3_cross_cohort_confidence.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Wrote {out_path}")
