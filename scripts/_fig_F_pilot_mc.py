"""
GEM1 - Figure F (supplementary): (F1) Bland-Altman magnitude-agreement plot,
production vs. pilot PR_objective, for the 2 groups with pilot data; (F2)
Monte Carlo precision -- distribution of PR_objective_se / |PR_objective|
(coefficient-of-variation-like ratio) among detected reactions, per group.
Reads 07b_bland_altman_source.csv (already computed) and reaction_scores.csv
per ACCEPTED group for the CV panel. Read-only.
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROD_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_production")
TABLES = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "tables")
FIGS = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "figures")
os.makedirs(FIGS, exist_ok=True)

GROUPS = [
    "GSE126848__NASH", "GSE126848__healthy", "GSE126848__obese_no_NAFLD",
    "GSE126848__steatosis", "GSE135251__NASH", "GSE135251__healthy",
    "GSE135251__steatosis", "GSE89632__NASH", "GSE89632__healthy",
    "GSE89632__steatosis",
]

ba = pd.read_csv(os.path.join(TABLES, "07b_bland_altman_source.csv"))
corr = pd.read_csv(os.path.join(TABLES, "07a_pilot_vs_production_correlation.csv"))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

colors = {"GSE126848__obese_no_NAFLD": "#0072B2", "GSE135251__healthy": "#E69F00"}
for g, sub in ba.groupby("group"):
    ax1.scatter(sub["mean_prod_pilot"], sub["diff_prod_minus_pilot"], s=14, alpha=0.6,
                color=colors.get(g, "gray"), label=g.replace("__", "/"))
    row = corr[corr["group"] == g].iloc[0]
    ax1.axhline(row["bland_altman_bias"], color=colors.get(g, "gray"), linestyle="-", linewidth=1)
    ax1.axhline(row["bland_altman_LoA_lo"], color=colors.get(g, "gray"), linestyle="--", linewidth=0.8)
    ax1.axhline(row["bland_altman_LoA_hi"], color=colors.get(g, "gray"), linestyle="--", linewidth=0.8)
ax1.axhline(0, color="black", linewidth=0.6)
ax1.set_xlabel("Mean of (production, pilot) PR_objective")
ax1.set_ylabel("Production - pilot PR_objective")
ax1.set_title("F1. Bland-Altman magnitude agreement\n(solid=bias, dashed=95% limits of agreement, per group)")
ax1.legend(fontsize=8, frameon=False)

cv_rows = []
for g in GROUPS:
    with open(os.path.join(PROD_DIR, g, "manifest.json")) as f:
        manifest = json.load(f)
    scores = pd.read_csv(os.path.join(PROD_DIR, g, "reaction_scores.csv"))
    det = scores[scores["occ_obj_k"] > 0].copy()
    det["cv"] = det["PR_objective_se"] / det["PR_objective"].abs()
    cv_rows.append(det["cv"].values)

bp = ax2.boxplot(cv_rows, tick_labels=[g.replace("__", "\n") for g in GROUPS], showfliers=False, patch_artist=True)
for patch in bp["boxes"]:
    patch.set_facecolor("#56B4E9")
    patch.set_alpha(0.7)
ax2.set_yscale("log")
ax2.set_ylabel("SE(PR_objective) / |PR_objective|  (detected reactions only, log scale)")
ax2.set_title("F2. Monte Carlo precision by group\n(M=800 backgrounds; background-sampling variance only)")
ax2.tick_params(axis="x", labelsize=7, rotation=90)

fig.suptitle("Figure F (supplementary) - Pilot/production agreement and Monte Carlo precision", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figF_pilot_mc_supplement.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)

source = ba.copy()
source.to_csv(os.path.join(FIGS, "figF_source_data.csv"), index=False)
print("Wrote figF_pilot_mc_supplement.{pdf,svg,png} + figF_source_data.csv")
