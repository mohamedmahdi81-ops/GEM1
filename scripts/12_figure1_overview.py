"""
GEM1 - Step 12: Figure 1 - GEM1 Framework Overview.

Conceptual/schematic figure (per docs/MANUSCRIPT_FIGURES.md Figure 1). No
statistical test applies -- the spec explicitly marks every "Statistical
analysis" field N/A for this figure, since it depicts pipeline architecture,
not a quantitative result. This script therefore draws a schematic rather
than plotting data, but it DOES verify (not assume) that every pipeline
output the schematic claims to depict actually exists on disk, and fails
loudly if any is missing, so the figure is never generated from a pipeline
state it misrepresents.

Panel A: horizontal 4-stage pipeline schematic (Steps 1-11), matching the
four named components of the Scientific Claim:
  (1) data/model acquisition
  (2) context-specific reconstruction (iMAT/GIMME/FASTCORE/tINIT -> one base model)
  (3) uncertainty analysis (flux sampling + perturbation testing, parallel)
  (4) consensus + hierarchical confidence engine + empirical calibration -> calibrated_confidence_score

Panel B (inset): the confidence engine's 4 evidence layers (Layer 0 raw
calls -> Layer 1 consensus -> Layer 2 GIMME corroboration, annotation-only,
drawn dashed to show it does not feed the score -> Layer 3 calibrated
combination). Reused visually (same box style) in Figure 5 Panel A for the
nested-CV schematic, per the spec's "reused visually in later figures for
continuity" instruction.

Caption text below is transcribed verbatim from the spec's "Caption outline".
"""

import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Verify required pipeline outputs exist (per the spec's "Verified pipeline
# support" table) before drawing a schematic that claims they do.
# ---------------------------------------------------------------------------
REQUIRED_PATHS = [
    os.path.join(BASE_DIR, "models", "Human-GEM.json"),
    os.path.join(BASE_DIR, "models", "Human-GEM.xml"),
    os.path.join(BASE_DIR, "data", "geo_cohorts"),
    os.path.join(BASE_DIR, "data", "context_specific_models", "extraction_manifest.csv"),
    os.path.join(BASE_DIR, "data", "flux_analysis", "flux_manifest.csv"),
    os.path.join(BASE_DIR, "data", "flux_sampling", "flux_sampling_manifest.csv"),
    os.path.join(BASE_DIR, "data", "perturbation_testing", "perturbation_manifest.csv"),
    os.path.join(BASE_DIR, "data", "consensus_scores", "consensus_manifest.csv"),
    os.path.join(BASE_DIR, "data", "confidence_scores", "confidence_manifest.csv"),
    os.path.join(BASE_DIR, "data", "calibration", "calibration_manifest.csv"),
    os.path.join(BASE_DIR, "data", "calibration", "calibration_weights_and_validation.json"),
]
missing = [p for p in REQUIRED_PATHS if not os.path.exists(p)]
if missing:
    raise FileNotFoundError(
        "Figure 1 depicts these pipeline stages as complete; the following "
        f"required outputs are missing, so the schematic would misrepresent "
        f"pipeline state: {missing}"
    )

STAGE_COLOR = "#3b6ea5"
UNCERTAINTY_COLOR = "#c0783c"
CONFIDENCE_COLOR = "#5f9e6e"
LAYER_COLORS = ["#8a8a8a", "#3b6ea5", "#b0b0b0", "#5f9e6e"]


def draw_box(ax, xy, w, h, text, facecolor, fontsize=9, dashed=False, textcolor="white"):
    box = FancyBboxPatch(
        xy, w, h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.4,
        edgecolor="black",
        facecolor=facecolor,
        linestyle="dashed" if dashed else "solid",
        zorder=2,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + w / 2, xy[1] + h / 2, text,
        ha="center", va="center", fontsize=fontsize, color=textcolor,
        wrap=True, zorder=3,
    )


def draw_arrow(ax, start, end, style="-|>", color="black", lw=1.6, connectionstyle="arc3,rad=0.0"):
    arrow = FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=14,
        color=color, lw=lw, connectionstyle=connectionstyle, zorder=1,
    )
    ax.add_patch(arrow)


fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 1, height_ratios=[2.1, 1.3], hspace=0.35)

# ---------------------------------------------------------------------------
# Panel A: 4-stage horizontal pipeline schematic
# ---------------------------------------------------------------------------
axA = fig.add_subplot(gs[0])
axA.set_xlim(0, 10)
axA.set_ylim(0, 4.2)
axA.axis("off")
axA.set_title("A", loc="left", fontsize=15, fontweight="bold")

# Stage 1: data/model acquisition
draw_box(axA, (0.1, 2.6), 1.9, 1.1, "Human-GEM\nbase model\n(12,931 rxns)", STAGE_COLOR)
draw_box(axA, (0.1, 1.1), 1.9, 1.1,
         "3 GEO cohorts\nGSE89632/126848/135251\n(NAFLD/MASLD)", STAGE_COLOR)

# Stage 2: context-specific reconstruction (4 algorithms -> 1 base model box)
algo_x = 2.55
algo_labels = ["iMAT", "GIMME", "FASTCORE", "tINIT"]
for i, lab in enumerate(algo_labels):
    y = 3.55 - i * 0.85
    draw_box(axA, (algo_x, y), 1.5, 0.65, lab, STAGE_COLOR, fontsize=9)
    draw_arrow(axA, (2.0, 3.15), (algo_x, y + 0.325), color="#555555", lw=1.0)
draw_box(axA, (4.4, 1.5), 1.55, 1.4, "Context-specific\nmodels\n(10 groups x 4 algos)",
         STAGE_COLOR, fontsize=8.5)
for i in range(4):
    y = 3.55 - i * 0.85
    draw_arrow(axA, (algo_x + 1.5, y + 0.325), (4.4, 2.2), color="#555555", lw=1.0)

# Stage 3: uncertainty analysis (parallel branches)
draw_box(axA, (6.3, 2.6), 1.7, 1.0, "Flux sampling\n(Monte Carlo)", UNCERTAINTY_COLOR, fontsize=8.5)
draw_box(axA, (6.3, 1.2), 1.7, 1.0, "Perturbation\ntesting (MC)", UNCERTAINTY_COLOR, fontsize=8.5)
draw_arrow(axA, (5.95, 2.2), (6.3, 3.1), color="black")
draw_arrow(axA, (5.95, 2.2), (6.3, 1.7), color="black")

# Stage 4: consensus + confidence engine + calibration -> output
draw_box(axA, (8.35, 1.5), 1.55, 1.4,
         "Consensus +\nconfidence engine +\nempirical calibration",
         CONFIDENCE_COLOR, fontsize=8)
draw_arrow(axA, (8.0, 3.1), (8.35, 2.55), color="black", connectionstyle="arc3,rad=-0.2")
draw_arrow(axA, (8.0, 1.7), (8.35, 1.95), color="black", connectionstyle="arc3,rad=0.2")

axA.annotate(
    "calibrated_confidence_score",
    xy=(9.9, 2.2), xytext=(9.95, 2.2),
    ha="left", va="center", fontsize=9.5, fontweight="bold", color="#2a5a35",
    arrowprops=dict(arrowstyle="-|>", color="#2a5a35", lw=1.8),
)

stage_labels = [
    "1. Data / model\nacquisition",
    "2. Context-specific\nreconstruction",
    "3. Uncertainty\nanalysis",
    "4. Consensus +\nconfidence + calibration",
]
stage_x = [1.0, 4.6, 7.1, 9.0]
for x, lab in zip(stage_x, stage_labels):
    axA.text(x, 4.0, lab, ha="center", va="bottom", fontsize=9.5, fontweight="bold")

# ---------------------------------------------------------------------------
# Panel B: confidence engine's 4 evidence layers (inset, reused in Fig 5)
# ---------------------------------------------------------------------------
axB = fig.add_subplot(gs[1])
axB.set_xlim(0, 10)
axB.set_ylim(0, 2.2)
axB.axis("off")
axB.set_title("B", loc="left", fontsize=15, fontweight="bold")

layer_texts = [
    "Layer 0\nRaw per-algorithm\nactive/inactive calls",
    "Layer 1\nprimary_consensus_score\n(FASTCORE+iMAT+tINIT)",
    "Layer 2\nGIMME corroboration\n(annotation only)",
    "Layer 3\ncalibrated_confidence_score\n(rank-weighted combination)",
]
box_w, gap = 2.05, 0.35
x0 = 0.3
for i, (txt, color) in enumerate(zip(layer_texts, LAYER_COLORS)):
    x = x0 + i * (box_w + gap)
    dashed = (i == 2)
    draw_box(axB, (x, 0.7), box_w, 1.0, txt, color, fontsize=8, dashed=dashed)
    if i < 3:
        draw_arrow(axB, (x + box_w, 1.2), (x + box_w + gap, 1.2), color="black")

axB.text(
    x0 + 2 * (box_w + gap) + box_w / 2, 0.45,
    "dashed = does not feed the score\n(independent robustness check only)",
    ha="center", va="top", fontsize=7.5, style="italic", color="#555555",
)

caption = (
    "Figure 1. Overview of the GEM1 pipeline. (A) Four-stage workflow from raw transcriptomic "
    "cohorts through context-specific reconstruction (four independent algorithms), uncertainty "
    "quantification (flux sampling and perturbation robustness), and a hierarchical, "
    "empirically-calibrated confidence engine, applied here to three independent NAFLD/MASLD "
    "cohorts. (B) Structure of the confidence engine's evidence layers."
)
fig.text(0.5, 0.01, caption, ha="center", va="bottom", fontsize=9, wrap=True)

out_path = os.path.join(OUT_DIR, "figure1_overview.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Wrote {out_path}")
