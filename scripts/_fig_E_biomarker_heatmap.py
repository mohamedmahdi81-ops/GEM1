"""
GEM1 - Figure E: 25 curated calibration/biomarker reactions x 10 groups,
categorical heatmap preserving NE (not eligible) / evaluated-zero / evaluated-
nonzero as three distinct visual states, with nonzero cells additionally
colored by PR_objective magnitude (log10, normalized within this panel's
nonzero cells only -- noted in the source data / caption, not a network-wide
scale). Reads reaction_scores.csv per ACCEPTED group and imports the frozen
BIOMARKERS_FULL set from 09_calibration.py (not redefined). Read-only.
"""
import importlib.util
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

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

spec = importlib.util.spec_from_file_location("step09_calibration", os.path.join(BASE_DIR, "scripts", "09_calibration.py"))
_calib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_calib)
BIOMARKERS_FULL = _calib.BIOMARKERS_FULL
biomarker_rxns = sorted({rid for lst in BIOMARKERS_FULL.values() for rid in lst})
rxn_to_set = {}
for name, lst in BIOMARKERS_FULL.items():
    for rid in lst:
        rxn_to_set.setdefault(rid, []).append(name)

status = pd.DataFrame(index=biomarker_rxns, columns=GROUPS, dtype=object)
pr_val = pd.DataFrame(index=biomarker_rxns, columns=GROUPS, dtype=float)

for g in GROUPS:
    with open(os.path.join(PROD_DIR, g, "manifest.json")) as f:
        manifest = json.load(f)
    assert manifest["status"] == "ACCEPTED"
    scores = pd.read_csv(os.path.join(PROD_DIR, g, "reaction_scores.csv")).set_index("reaction_id")
    for rid in biomarker_rxns:
        if rid not in scores.index:
            status.loc[rid, g] = "NE"
        elif scores.loc[rid, "occ_obj_k"] > 0:
            status.loc[rid, g] = "nonzero"
            pr_val.loc[rid, g] = scores.loc[rid, "PR_objective"]
        else:
            status.loc[rid, g] = "zero"

long_rows = []
for rid in biomarker_rxns:
    for g in GROUPS:
        long_rows.append({"reaction_id": rid, "biomarker_set": "+".join(rxn_to_set[rid]), "group": g,
                           "status": status.loc[rid, g], "PR_objective": pr_val.loc[rid, g]})
pd.DataFrame(long_rows).to_csv(os.path.join(FIGS, "figE_source_data.csv"), index=False)

nonzero_vals = pr_val.values[status.values == "nonzero"]
log_vals = np.log10(np.clip(nonzero_vals.astype(float), 1e-12, None))
norm = Normalize(vmin=log_vals.min(), vmax=log_vals.max())
cmap = plt.get_cmap("viridis")

n_r, n_g = len(biomarker_rxns), len(GROUPS)
rgba = np.zeros((n_r, n_g, 4))
NE_COLOR = (0.55, 0.55, 0.55, 1.0)
ZERO_COLOR = (0.92, 0.92, 0.92, 1.0)
for i, rid in enumerate(biomarker_rxns):
    for j, g in enumerate(GROUPS):
        s = status.loc[rid, g]
        if s == "NE":
            rgba[i, j] = NE_COLOR
        elif s == "zero":
            rgba[i, j] = ZERO_COLOR
        else:
            v = pr_val.loc[rid, g]
            logv = np.log10(max(v, 1e-12))
            rgba[i, j] = cmap(norm(logv))

fig, ax = plt.subplots(figsize=(9, 9))
ax.imshow(rgba, aspect="auto")
ax.set_xticks(range(n_g))
ax.set_xticklabels([g.replace("__", "/") for g in GROUPS], rotation=90, fontsize=8)
ax.set_yticks(range(n_r))
ax.set_yticklabels([f"{rid} ({'+'.join(rxn_to_set[rid])})" for rid in biomarker_rxns], fontsize=7)
ax.set_title("Figure E - 25 curated biomarker reactions x 10 groups\n"
             "gray=NE (not dropout-eligible), light=evaluated-zero, colored=evaluated-nonzero (log10 PR_objective, panel-relative scale)")

sm = ScalarMappable(norm=norm, cmap=cmap)
cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
cbar.set_label("log10(PR_objective) | nonzero cells only, this panel's range")

legend_handles = [plt.Rectangle((0, 0), 1, 1, color=NE_COLOR, label="NE"),
                   plt.Rectangle((0, 0), 1, 1, color=ZERO_COLOR, label="evaluated-zero")]
ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.25, 1.0), frameon=False, fontsize=8)

fig.tight_layout()
for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figE_biomarker_heatmap.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("Wrote figE_biomarker_heatmap.{pdf,svg,png} + figE_source_data.csv")
