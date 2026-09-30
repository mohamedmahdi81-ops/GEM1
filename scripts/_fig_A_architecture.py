"""
GEM1 - Figure A: per-group NE / evaluated-zero / evaluated-nonzero architecture,
with feasibility-flip fraction annotated separately (not assumed nested inside
nonzero without verification). Reads only 01_sparsity_occurrence_by_group.csv
(already computed, read-only w.r.t. production). Vector (PDF+SVG) + PNG.
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

GRAY, LGRAY, ORANGE, BLUE = "#595959", "#BFBFBF", "#E69F00", "#0072B2"

t1 = pd.read_csv(os.path.join(TABLES, "01_sparsity_occurrence_by_group.csv"))
t1 = t1.sort_values("group").reset_index(drop=True)
t1.to_csv(os.path.join(FIGS, "figA_source_data.csv"), index=False)

labels = [g.replace("__", " / ") for g in t1["group"]]
ne = t1["n_NE"].values
zero = t1["n_evaluated_zero"].values
nonzero = t1["n_evaluated_nonzero_delta_gt_TAU"].values
flip_frac_of_eval = t1["frac_flip_of_evaluated"].values * 100

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.4, 1]})

y = range(len(labels))
ax1.barh(y, ne, color=GRAY, label="NE (not dropout-eligible)")
ax1.barh(y, zero, left=ne, color=LGRAY, label="Evaluated, $\\delta_{obj} \\leq \\tau$ (all M=800 bg)")
ax1.barh(y, nonzero, left=ne + zero, color=ORANGE, label="Evaluated, $\\delta_{obj} > \\tau$ in $\\geq$1 bg")
ax1.set_yticks(list(y))
ax1.set_yticklabels(labels, fontsize=9)
ax1.set_xlabel("Reactions (of 5,718 full-network union)")
ax1.set_title("A1. NE / evaluated-zero / evaluated-nonzero architecture")
ax1.legend(loc="lower right", fontsize=7.5, frameon=False)
ax1.invert_yaxis()

ax2.barh(y, flip_frac_of_eval, color=BLUE)
for i, v in enumerate(flip_frac_of_eval):
    ax2.text(v + 0.02, i, f"{v:.2f}%", va="center", fontsize=8)
ax2.set_yticks(list(y))
ax2.set_yticklabels(labels, fontsize=9)
ax2.set_xlabel("% of evaluated reactions with a feasibility flip")
ax2.set_title("A2. Feasibility-flip fraction\n(not assumed nested in A1's nonzero segment)")
ax2.invert_yaxis()

fig.suptitle("Figure A - Per-group perturbation-response architecture (TAU=1e-6, M=800)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figA_architecture.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("Wrote figA_architecture.{pdf,svg,png} + figA_source_data.csv")
