"""
GEM1 - Figure C: observed cross-group recurrence (reactions detected in k of
10 groups) vs. an eligibility-aware permutation null (5000 reps, per-group
random subset of size n_detected(g) drawn from that group's own eligible
pool, exactly preserving each group's eligibility set and detection count).
Reads only 03b_recurrence_null_summary.csv (already computed, read-only).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
TABLES = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "tables")
FIGS = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "figures")
os.makedirs(FIGS, exist_ok=True)

ORANGE, BLUE = "#E69F00", "#0072B2"

t3b = pd.read_csv(os.path.join(TABLES, "03b_recurrence_null_summary.csv"))
t3b.to_csv(os.path.join(FIGS, "figC_source_data.csv"), index=False)

k = t3b["n_groups_detected"].values
obs = t3b["observed_n_reactions"].values
null_mean = t3b["null_mean"].values
null_lo = t3b["null_p2_5"].values
null_hi = t3b["null_p97_5"].values

fig, ax = plt.subplots(figsize=(8, 5.5))
width = 0.38
ax.bar(k - width / 2, obs, width=width, color=ORANGE, label="Observed")
yerr_lo = np.clip(null_mean - null_lo, 0, None)
yerr_hi = np.clip(null_hi - null_mean, 0, None)
ax.bar(k + width / 2, null_mean, width=width, color=BLUE, alpha=0.85,
       yerr=[yerr_lo, yerr_hi], capsize=3,
       error_kw={"linewidth": 1}, label="Null mean (95% empirical band)")
ax.set_yscale("symlog", linthresh=1)
ax.set_xlabel("Number of groups (of 10) in which a reaction is detected ($\\delta_{obj}>\\tau$ in $\\geq$1 background)")
ax.set_ylabel("Number of reactions (symlog scale)")
ax.set_xticks(list(k))
ax.set_title("Figure C - Cross-group recurrence vs. eligibility-aware permutation null\n"
             "(N=5000 permutations; per-group shuffle within that group's own eligible pool)")
ax.legend(frameon=False)
ax.grid(alpha=0.25, axis="y")
fig.tight_layout()

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(FIGS, f"figC_recurrence_null.{ext}"),
                dpi=300 if ext == "png" else None, bbox_inches="tight")
plt.close(fig)
print("Wrote figC_recurrence_null.{pdf,svg,png} + figC_source_data.csv")
