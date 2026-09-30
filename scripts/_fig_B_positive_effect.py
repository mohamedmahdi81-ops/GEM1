"""
GEM1 - Figure B: ECDF of PR_objective among detected reactions (occ_obj_k>0)
per group, log-x axis so the large zero/near-zero mass (excluded here by
construction) does not obscure the positive-effect distribution. Reads only
reaction_scores.csv per ACCEPTED group (read-only).
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROD_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_production")
FIGS = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "figures")
os.makedirs(FIGS, exist_ok=True)

GROUPS = [
    "GSE126848__NASH", "GSE126848__healthy", "GSE126848__obese_no_NAFLD",
    "GSE126848__steatosis", "GSE135251__NASH", "GSE135251__healthy",
    "GSE135251__steatosis", "GSE89632__NASH", "GSE89632__healthy",
    "GSE89632__steatosis",
]
CMAP = plt.get_cmap("tab10")

fig, ax = plt.subplots(figsize=(8, 6))
source_rows = []
for i, g in enumerate(GROUPS):
    with open(os.path.join(PROD_DIR, g, "manifest.json")) as f:
        manifest = json.load(f)
    assert manifest["status"] == "ACCEPTED"
    scores = pd.read_csv(os.path.join(PROD_DIR, g, "reaction_scores.csv"))
    det = scores.loc[scores["occ_obj_k"] > 0, ["reaction_id", "PR_objective"]].copy()
    det["group"] = g
    source_rows.append(det)
    x = np.sort(det["PR_objective"].values)
    y = np.arange(1, len(x) + 1) / len(x)
    ax.step(x, y, where="post", label=g.replace("__", "/"), color=CMAP(i % 10), linewidth=1.4)

ax.set_xscale("log")
ax.set_xlabel("PR_objective (log scale), reactions with occ_obj_k > 0 only")
ax.set_ylabel("Empirical CDF")
ax.set_title("Figure B - Positive-effect PR_objective distribution (ECDF, per group)\n"
              "Zero/sub-TAU mass (91.7-96.1% of eligible reactions per group) excluded by construction")
ax.legend(loc="lower right", fontsize=7.5, ncol=1, frameon=False)
ax.grid(alpha=0.25, which="both")
fig.tight_layout()

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figB_positive_effect_ecdf.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)

source = pd.concat(source_rows, ignore_index=True)
source.to_csv(os.path.join(FIGS, "figB_source_data.csv"), index=False)
print(f"Wrote figB_positive_effect_ecdf.{{pdf,svg,png}} + figB_source_data.csv ({len(source)} rows)")
