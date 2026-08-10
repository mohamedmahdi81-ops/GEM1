"""
GEM1 - Step 12: Figure 2 - Positioning GEM1 Relative to Existing Frameworks.

Per docs/MANUSCRIPT_FIGURES.md Figure 2: a qualitative literature-comparison
table, NOT derived from GEM1's own computed data. The spec is explicit:
"populate the table only from this literature-verified specification. Do not
independently reinterpret or extend literature claims." Every cell, the two
footnotes, and the citation list below are therefore transcribed verbatim
from the frozen spec -- nothing here is computed, inferred, or rephrased.

No statistical analysis applies (spec: N/A throughout).
"""

import os
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

COLUMNS = ["Criterion", "Mardinoglu 2014", "Mardinoglu 2016", "BayFlux", "GEMsembler", "GEM1"]

ROWS = [
    ["Context-specific reconstruction from transcriptomics", "Yes (INIT)", "Yes (iMAT)", "No", "No", "Yes"],
    ["Multiple reconstruction algorithms combined (1)", "No", "No", "No", "Yes", "Yes"],
    ["Flux uncertainty quantification (2)", "No", "Partial (FVA bounds)", "Yes", "No", "Yes"],
    ["Cross-algorithm consensus scoring", "No", "No", "No", "Yes", "Yes"],
    ["Hierarchical/decomposable confidence score", "No", "No", "No", "No", "Yes"],
    ["Empirically calibrated via cross-validated weight fitting", "No", "No", "No", "No", "Yes (nested LOCO/LOBO)"],
    ["Disease biomarker application", "Yes (NAFLD)", "Yes (NAFLD)", "No", "No", "Yes (NAFLD)"],
    ["Cross-cohort validation", "Partial (1 validation cohort)", "No", "No", "No", "Yes (3 cohorts)"],
]

FOOTNOTE_1 = (
    "(1) GEMsembler integrates outputs from different whole-model reconstruction tools across "
    "microbial genomes (demonstrated on Lactiplantibacillus plantarum and Escherichia coli, "
    "outperforming gold-standard models in auxotrophy and gene essentiality predictions). GEM1 "
    "integrates outputs from different context-specific extraction algorithms (iMAT/GIMME/"
    "FASTCORE/tINIT) applied to transcriptomic data against one shared human base model. "
    "Different consensus paradigms, not the same capability applied to different organisms."
)
FOOTNOTE_2 = (
    "(2) BayFlux uses Bayesian inference and MCMC to quantify metabolic flux uncertainty from "
    "13C-labeling data. GEM1 quantifies uncertainty in context-specific metabolic predictions "
    "through flux sampling and Monte Carlo perturbation robustness, evaluated against "
    "transcriptomics-derived models, not isotope-labeling data. Different uncertainty "
    "frameworks answering different questions."
)

CAPTION = (
    "Figure 2. Feature comparison of GEM1 against the two foundational NAFLD GSMM studies and "
    "two methodologically related but non-overlapping tools (BayFlux: flux uncertainty from "
    "isotope tracing; GEMsembler: cross-tool structural consensus in microbial GEMs). GEM1's "
    "contribution is the combination, not unique ownership of any individual capability (see "
    "footnotes)."
)

fig, ax = plt.subplots(figsize=(13, 5.2))
ax.axis("off")

table = ax.table(
    cellText=ROWS,
    colLabels=COLUMNS,
    cellLoc="center",
    loc="center",
    colWidths=[0.34, 0.13, 0.14, 0.10, 0.11, 0.16],
)
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.9)

for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor("#999999")
    if row == 0:
        cell.set_facecolor("#3b6ea5")
        cell.set_text_props(color="white", fontweight="bold")
    else:
        cell.set_facecolor("#eef3f8" if row % 2 == 0 else "white")
        if col == 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.02
        if col == 5:  # GEM1 column
            cell.set_facecolor("#d9ead3" if row % 2 == 0 else "#e9f5e6")
            cell.set_text_props(fontweight="bold")

fig.text(0.06, 0.10, FOOTNOTE_1, ha="left", va="top", fontsize=7.3, wrap=True)
fig.text(0.06, 0.055, FOOTNOTE_2, ha="left", va="top", fontsize=7.3, wrap=True)
fig.text(0.5, 0.005, CAPTION, ha="center", va="bottom", fontsize=8.5, wrap=True)

out_path = os.path.join(OUT_DIR, "figure2_novelty_table.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Wrote {out_path}")
