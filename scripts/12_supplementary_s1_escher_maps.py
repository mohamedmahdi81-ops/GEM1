"""
GEM1 - Step 12: Supplementary Figure S1 - Escher Pathway Confidence Maps.

Per docs/MANUSCRIPT_FIGURES.md Supplementary S1. Uses the Escher library
itself (not a matplotlib imitation), per this task's explicit instruction.

The spec's own "Verified support" section already established that no
pre-built Human-GEM Escher map ships with the installed `escher` package
(only RECON1/iJO1366/iMM904 maps are bundled, none of which use Human-GEM's
MAR-prefixed reaction/metabolite IDs) -- so this script builds two small,
hand-laid-out, schema-valid Escher maps from scratch, one per representative
pathway named in the spec's "Panel layout": urea cycle and serine/glycine
one-carbon metabolism (SHMT). Reaction/metabolite IDs, names, and full
stoichiometry are read directly from models/Human-GEM.json (not guessed),
using exactly the reaction set already locked in for these two biomarkers in
scripts/09_calibration.py's BIOMARKERS_FULL dict (serine_glycine_shmt,
urea_cycle).

SIMPLIFICATION, disclosed rather than hidden: to keep each map a readable
schematic, only each reaction's primary substrate/product (and a couple of
key named cofactors, e.g. THF/5,10-methylene-THF, NH3/CO2) are drawn as
nodes; pure "energy currency" cofactors (ATP, ADP, Pi, PPi, AMP, H+) are
omitted from the map entirely. This affects only the VISUAL layout, not any
scientific claim -- the full stoichiometry for every reaction is on record
in models/Human-GEM.json and scripts/09_calibration.py's biomarker mapping,
and the omitted species are standard secondary cofactors that most published
pathway cartoons (including Escher's own convention of a
`hide_secondary_metabolites` option) also suppress for readability.

Data overlaid: calibrated_confidence_score (Step 11 output) for a
representative cohort/group pair. GSE135251 was chosen as the representative
cohort because it is this project's largest, most clinically-graded cohort
(spans healthy -> steatosis -> NASH -> fibrosis/cirrhosis), and healthy vs.
NASH as the representative group
pair because that is the clearest disease-contrast the spec's own framing
("high-confidence metabolic alterations", "biologically meaningful
pathways") calls for. This is a visualization choice, not a new scientific
claim -- the spec itself leaves the specific representative pathway/group
choice open ("the natural candidates").

Rendering: escher.Builder generates a standalone, self-contained HTML file
per (pathway, group) via save_html() -- this is real, unmodified Escher
output, referencing escher.min.js from unpkg.com (Escher's own documented
standalone-HTML mechanism; requires network access, confirmed available on
this machine). Each HTML is then screenshotted headlessly via Selenium +
Microsoft Edge (confirmed working on this machine) to produce a static PNG,
since Escher itself has no headless/selenium-free static image export. The
four screenshots (2 pathways x 2 groups) are assembled into one composite
supplementary figure with matplotlib, on a shared 0-1 confidence color
scale.

No statistical test applies (spec: N/A throughout).
"""

import os
import io
import json
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
CAL_DIR = os.path.join(BASE_DIR, "data", "calibration")
MODEL_PATH = os.path.join(BASE_DIR, "models", "Human-GEM.json")
OUT_DIR = os.path.join(BASE_DIR, "results", "figures")
HTML_DIR = os.path.join(OUT_DIR, "escher_html")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(HTML_DIR, exist_ok=True)

REPRESENTATIVE_COHORT = "GSE135251"
REPRESENTATIVE_GROUPS = ["healthy", "NASH"]

# Same reaction sets as scripts/09_calibration.py's BIOMARKERS_FULL.
SHMT_REACTIONS = ["MAR03845", "MAR04792"]
UREA_CYCLE_REACTIONS = ["MAR03873", "MAR03809", "MAR03811", "MAR03813", "MAR03816", "MAR08426"]

COLOR_LOW = "#d73027"
COLOR_MID = "#fee08b"
COLOR_HIGH = "#1a9850"
CMAP = LinearSegmentedColormap.from_list("confidence", [COLOR_LOW, COLOR_MID, COLOR_HIGH])

# ---------------------------------------------------------------------------
# Load reaction metadata from Human-GEM directly (not guessed)
# ---------------------------------------------------------------------------
with open(MODEL_PATH) as f:
    model = json.load(f)
MET_BY_ID = {m["id"]: m for m in model["metabolites"]}
RXN_BY_ID = {r["id"]: r for r in model["reactions"]}


def rxn_meta(rxn_id):
    r = RXN_BY_ID[rxn_id]
    genes = r.get("gene_reaction_rule", "") or ""
    return {
        "name": r.get("name", rxn_id),
        "gene_reaction_rule": genes,
        "metabolites": r["metabolites"],
    }


# ---------------------------------------------------------------------------
# Minimal Escher map builder (schema verified against a real downloaded
# Escher map, e_coli_core.Core metabolism, this session -- node_type values
# 'metabolite'/'midmarker', reaction fields name/bigg_id/reversibility/
# metabolites/segments/label_x/label_y all confirmed against that reference).
# ---------------------------------------------------------------------------
class MapBuilder:
    def __init__(self, map_name):
        self.map_name = map_name
        self.nodes = {}
        self.reactions = {}
        self._next_id = 1

    def _new_id(self):
        nid = str(self._next_id)
        self._next_id += 1
        return nid

    def add_metabolite(self, x, y, name, bigg_id, primary=True, label_dx=25, label_dy=0):
        nid = self._new_id()
        self.nodes[nid] = {
            "node_type": "metabolite",
            "x": x, "y": y,
            "bigg_id": bigg_id, "name": name,
            "label_x": x + label_dx, "label_y": y + label_dy,
            "node_is_primary": primary,
        }
        return nid

    def add_reaction(self, rxn_id, name, gene_reaction_rule, mid_x, mid_y, edges, reversibility=False):
        """edges: list of (node_id, coefficient, bigg_id) tuples -- one
        segment per edge, all radiating from a single midmarker node."""
        mid_id = self._new_id()
        self.nodes[mid_id] = {"node_type": "midmarker", "x": mid_x, "y": mid_y}
        segments = {}
        metabolites = []
        for node_id, coeff, bigg_id in edges:
            seg_id = self._new_id()
            if coeff < 0:
                frm, to = node_id, mid_id
            else:
                frm, to = mid_id, node_id
            segments[seg_id] = {"from_node_id": frm, "to_node_id": to, "b1": None, "b2": None}
            metabolites.append({"coefficient": coeff, "bigg_id": bigg_id})
        genes = []
        for tok in gene_reaction_rule.replace("(", "").replace(")", "").split(" or "):
            tok = tok.strip()
            if tok and tok.lower() != "and":
                genes.append({"bigg_id": tok, "name": tok})
        self.reactions[rxn_id] = {
            "name": name, "bigg_id": rxn_id, "reversibility": reversibility,
            "label_x": mid_x, "label_y": mid_y - 30,
            "gene_reaction_rule": gene_reaction_rule,
            "genes": genes,
            "metabolites": metabolites,
            "segments": segments,
        }

    def to_map_json(self, width, height):
        meta = {
            "map_name": self.map_name,
            "map_id": self.map_name.replace(" ", "_"),
            "map_description": f"GEM1 Supplementary Figure S1 -- {self.map_name} (custom layout, "
                                 "built from models/Human-GEM.json; not a published Escher map).",
            "homepage": "https://escher.github.io",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        }
        data = {
            "reactions": self.reactions,
            "nodes": self.nodes,
            "canvas": {"x": -50, "y": -50, "width": width, "height": height},
            "text_labels": {},
        }
        return [meta, data]


def build_shmt_map():
    b = MapBuilder("Serine/glycine one-carbon metabolism (SHMT)")

    # Row 1: MAR03845 (SHMT1, cytosol) -- substrates serine(c), THF(c); products glycine(c), 5,10-mTHF(c) [H2O omitted]
    ser_c = b.add_metabolite(150, 200, "L-serine (c)", "MAM02896c", primary=True)
    gly_c = b.add_metabolite(750, 200, "Glycine (c)", "MAM01986c", primary=True)
    thf_c = b.add_metabolite(450, 60, "THF (c)", "MAM02980c", primary=False, label_dy=-15)
    mthf_c = b.add_metabolite(450, 340, "5,10-methylene-THF (c)", "MAM01045c", primary=False, label_dy=15)
    m1 = rxn_meta("MAR03845")
    b.add_reaction(
        "MAR03845", "SHMT1 (cytosol): " + m1["name"], m1["gene_reaction_rule"], 450, 200,
        edges=[(ser_c, -1.0, "MAM02896c"), (thf_c, -1.0, "MAM02980c"),
               (gly_c, 1.0, "MAM01986c"), (mthf_c, 1.0, "MAM01045c")],
    )

    # Row 2: MAR04792 (SHMT2, mitochondria) -- substrates 5,10-mTHF(m), glycine(m); products serine(m), THF(m) [H2O omitted]
    ser_m = b.add_metabolite(150, 600, "L-serine (m)", "MAM02896m", primary=True)
    gly_m = b.add_metabolite(750, 600, "Glycine (m)", "MAM01986m", primary=True)
    thf_m = b.add_metabolite(450, 460, "THF (m)", "MAM02980m", primary=False, label_dy=-15)
    mthf_m = b.add_metabolite(450, 740, "5,10-methylene-THF (m)", "MAM01045m", primary=False, label_dy=15)
    m2 = rxn_meta("MAR04792")
    b.add_reaction(
        "MAR04792", "SHMT2 (mitochondria): " + m2["name"], m2["gene_reaction_rule"], 450, 600,
        edges=[(gly_m, -1.0, "MAM01986m"), (mthf_m, -1.0, "MAM01045m"),
               (ser_m, 1.0, "MAM02896m"), (thf_m, 1.0, "MAM02980m")],
    )

    return b.to_map_json(900, 800), SHMT_REACTIONS


def build_urea_cycle_map():
    b = MapBuilder("Urea cycle")
    cx, cy, R = 500, 500, 320

    def pt(angle_deg, radius=R):
        rad = np.radians(angle_deg)
        return cx + radius * np.cos(rad), cy - radius * np.sin(rad)

    x, y = pt(90)
    nh3 = b.add_metabolite(x, y, "NH3 (m)", "MAM02578m", primary=True)
    x, y = pt(110, R + 80)
    co2 = b.add_metabolite(x, y, "CO2 (m)", "MAM01596m", primary=False)

    x, y = pt(30)
    carbp = b.add_metabolite(x, y, "Carbamoyl-P (m)", "MAM01420m", primary=True)

    x, y = pt(-30)
    citr = b.add_metabolite(x, y, "Citrulline (c)", "MAM01588c", primary=True)

    x, y = pt(-90)
    argsucc = b.add_metabolite(x, y, "Argininosuccinate (c)", "MAM01366c", primary=True)
    x, y = pt(-100, R + 90)
    asp = b.add_metabolite(x, y, "Aspartate (c)", "MAM01370c", primary=False)

    x, y = pt(-150)
    arg = b.add_metabolite(x, y, "Arginine (c)", "MAM01365c", primary=True)
    x, y = pt(-125, R + 90)
    fum = b.add_metabolite(x, y, "Fumarate (c)", "MAM01862c", primary=False)

    x, y = pt(150)
    orn = b.add_metabolite(x, y, "Ornithine", "MAM02658c", primary=True)

    x, y = pt(180, R - 150)
    urea = b.add_metabolite(x, y, "Urea", "MAM03121c", primary=False)

    # R1: MAR03873 (CPS1) -- NH3 + CO2 -> carbamoyl-P
    mx, my = pt(60)
    m = rxn_meta("MAR03873")
    b.add_reaction("MAR03873", "CPS1: " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(nh3, -1.0, "MAM02578m"), (co2, -1.0, "MAM01596m"), (carbp, 1.0, "MAM01420m")])

    # R2: MAR03809 (OTC) -- carbamoyl-P + ornithine -> citrulline
    mx, my = pt(0)
    m = rxn_meta("MAR03809")
    b.add_reaction("MAR03809", "OTC: " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(carbp, -1.0, "MAM01420m"), (orn, -1.0, "MAM02658c"), (citr, 1.0, "MAM01588c")])

    # R3: MAR03811 (ASS1) -- citrulline + aspartate -> argininosuccinate
    mx, my = pt(-60)
    m = rxn_meta("MAR03811")
    b.add_reaction("MAR03811", "ASS1: " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(citr, -1.0, "MAM01588c"), (asp, -1.0, "MAM01370c"), (argsucc, 1.0, "MAM01366c")])

    # R4: MAR03813 (ASL) -- argininosuccinate -> arginine + fumarate
    mx, my = pt(-120)
    m = rxn_meta("MAR03813")
    b.add_reaction("MAR03813", "ASL: " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(argsucc, -1.0, "MAM01366c"), (arg, 1.0, "MAM01365c"), (fum, 1.0, "MAM01862c")])

    # R5: MAR03816 (ARG1, cytosol) -- arginine -> ornithine + urea
    mx, my = pt(180, R + 60)
    m = rxn_meta("MAR03816")
    b.add_reaction("MAR03816", "ARG1 (cytosol): " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(arg, -1.0, "MAM01365c"), (orn, 1.0, "MAM02658c"), (urea, 1.0, "MAM03121c")])

    # R6: MAR08426 (ARG2, mitochondria) -- arginine -> ornithine + urea (parallel isoform)
    mx, my = pt(165, R - 80)
    m = rxn_meta("MAR08426")
    b.add_reaction("MAR08426", "ARG2 (mitochondria): " + m["name"], m["gene_reaction_rule"], mx, my,
                    edges=[(arg, -1.0, "MAM01365m"), (orn, 1.0, "MAM02658m"), (urea, 1.0, "MAM03121m")])

    return b.to_map_json(1200, 1200), UREA_CYCLE_REACTIONS


# ---------------------------------------------------------------------------
# Load calibrated_confidence_score for the representative cohort/groups
# ---------------------------------------------------------------------------
def load_scores(cohort, group, reaction_ids):
    path = os.path.join(CAL_DIR, f"{cohort}__{group}__calibrated.csv")
    df = pd.read_csv(path, usecols=["reaction_id", "calibrated_confidence_score"])
    df = df.set_index("reaction_id")
    scores = {}
    for rid in reaction_ids:
        if rid not in df.index:
            raise KeyError(f"{rid} not found in {path}")
        scores[rid] = float(df.loc[rid, "calibrated_confidence_score"])
    return scores


# ---------------------------------------------------------------------------
# Render one map+data combo to a standalone HTML via Escher, then screenshot
# it headlessly via Selenium + Edge to get a static PNG.
# ---------------------------------------------------------------------------
def render_map_to_png(map_json, reaction_data, out_png_path, out_html_path, window_size=(1000, 900)):
    from escher import Builder

    builder = Builder(
        map_json=json.dumps(map_json),
        reaction_data=reaction_data,
        reaction_styles=["color", "size", "text"],
        # Fixed absolute domain [0,1], NOT 'min'/'max' types (those auto-scale
        # to each panel's own data range, which would silently break the
        # "shared 0-1 scale across panels" claim in the composite caption --
        # caught by visually inspecting the first render, where each panel's
        # auto-generated Escher legend showed a different max, e.g. 0.56 vs
        # 0.79).
        reaction_scale=[
            {"type": "value", "value": 0.0, "color": COLOR_LOW, "size": 6},
            {"type": "value", "value": 0.5, "color": COLOR_MID, "size": 12},
            {"type": "value", "value": 1.0, "color": COLOR_HIGH, "size": 20},
        ],
        reaction_no_data_color="#cccccc",
        reaction_no_data_size=4,
        hide_secondary_metabolites=False,
        identifiers_on_map="bigg_id",
        never_ask_before_quit=True,
        menu="none",
        scroll_behavior="none",
    )
    builder.save_html(out_html_path)

    from selenium import webdriver
    from selenium.webdriver.edge.options import Options

    opts = Options()
    opts.add_argument("--headless=new")
    opts.add_argument(f"--window-size={window_size[0]},{window_size[1]}")
    driver = webdriver.Edge(options=opts)
    try:
        driver.get("file:///" + os.path.abspath(out_html_path).replace("\\", "/"))
        time.sleep(4)  # allow escher.min.js (from unpkg) to load and render the SVG
        driver.save_screenshot(out_png_path)
    finally:
        driver.quit()


def main():
    pathways = {
        "urea_cycle": build_urea_cycle_map(),
        "one_carbon_shmt": build_shmt_map(),
    }

    panel_pngs = {}  # (pathway_key, group) -> png path
    for pathway_key, (map_json, reaction_ids) in pathways.items():
        for group in REPRESENTATIVE_GROUPS:
            scores = load_scores(REPRESENTATIVE_COHORT, group, reaction_ids)
            print(f"{pathway_key} / {REPRESENTATIVE_COHORT}/{group}: {scores}")
            html_path = os.path.join(HTML_DIR, f"{pathway_key}__{REPRESENTATIVE_COHORT}__{group}.html")
            png_path = os.path.join(HTML_DIR, f"{pathway_key}__{REPRESENTATIVE_COHORT}__{group}.png")
            render_map_to_png(map_json, scores, png_path, html_path)
            panel_pngs[(pathway_key, group)] = png_path

    # ---------------------------------------------------------------------
    # Assemble composite supplementary figure
    # ---------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))
    pathway_titles = {"urea_cycle": "Urea cycle", "one_carbon_shmt": "Serine/glycine one-carbon metabolism (SHMT)"}
    for row, pathway_key in enumerate(["urea_cycle", "one_carbon_shmt"]):
        for col, group in enumerate(REPRESENTATIVE_GROUPS):
            ax = axes[row, col]
            img = mpimg.imread(panel_pngs[(pathway_key, group)])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(f"{pathway_titles[pathway_key]} -- {REPRESENTATIVE_COHORT}/{group}", fontsize=10)

    sm = ScalarMappable(cmap=CMAP)
    sm.set_array([0, 1])
    cbar = fig.colorbar(sm, ax=axes, orientation="horizontal", fraction=0.03, pad=0.04, aspect=40)
    cbar.set_label("calibrated_confidence_score", fontsize=9)

    caption = (
        "Supplementary Figure S1. Escher pathway maps showing calibrated confidence scores "
        f"overlaid on urea cycle and one-carbon (SHMT) metabolism, for {REPRESENTATIVE_COHORT} "
        "healthy vs. NASH groups, illustrating that GEM1's confidence estimates map onto "
        "biologically coherent, literature-recognized pathways rather than scattered, "
        "uninterpretable reactions. Reaction arrow color/width encode "
        "calibrated_confidence_score (0-1 scale, shared across panels); pure energy-currency "
        "cofactors (ATP/ADP/Pi/H+) are omitted from the layout for readability (see script "
        "docstring)."
    )
    fig.text(0.5, 0.01, caption, ha="center", va="bottom", fontsize=8.3, wrap=True)

    out_path = os.path.join(OUT_DIR, "supplementary_s1_escher_maps.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Wrote {out_path}")
    print(f"Standalone interactive Escher HTML files (real, unmodified Escher output) in {HTML_DIR}")


if __name__ == "__main__":
    main()
