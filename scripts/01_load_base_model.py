"""
GEM1 - Step 2: Load and verify the base genome-scale metabolic model.

Base model: Human-GEM (SysBioChalmers/Metabolic Atlas)
Rationale (see GEM1-novelty-dossier.md in the project):
  - Human-GEM descends from the HMR2.0 / iHepatocytes2322 lineage -- the exact
    model family used in the foundational NAFLD GSMM papers (Mardinoglu et al.
    2014, 2016) that GEM1's calibration set is built from. That shared lineage
    makes reaction/gene-level comparison to prior findings more direct.
  - Actively maintained, cleanly curated compartmentalization.
  - Distributed as SBML L3V1 FBCv3 -- loads cleanly into COBRApy and is the
    format troppo's extraction algorithms expect.

Run this on a machine with internet access (this script downloads the model).
Requires: cobra  (pip install cobra)
"""

import os
import urllib.request
import cobra

MODEL_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "models")
MODEL_PATH = os.path.join(MODEL_DIR, "Human-GEM.xml")
MODEL_URL = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.xml"

# Expected specs as of the current SysBioChalmers/Human-GEM main branch
# (verified via the repo README on 2026-07-10). If your downloaded version's
# counts differ substantially, the model has been updated upstream --
# note the new counts in your lab notebook, it isn't a bug.
EXPECTED_APPROX = {"reactions": 12931, "metabolites": 8461, "genes": 2848}


def download_model():
    os.makedirs(MODEL_DIR, exist_ok=True)
    if os.path.exists(MODEL_PATH):
        print(f"Model already present at {MODEL_PATH}, skipping download.")
        return
    print(f"Downloading Human-GEM from {MODEL_URL} ...")
    urllib.request.urlretrieve(MODEL_URL, MODEL_PATH)
    print("Download complete.")


def load_and_verify():
    model = cobra.io.read_sbml_model(MODEL_PATH)
    print(f"Loaded model: {model.id}")
    print(f"  Reactions:   {len(model.reactions)} (expected ~{EXPECTED_APPROX['reactions']})")
    print(f"  Metabolites: {len(model.metabolites)} (expected ~{EXPECTED_APPROX['metabolites']})")
    print(f"  Genes:       {len(model.genes)} (expected ~{EXPECTED_APPROX['genes']})")

    # Sanity check: model should grow biomass under default (unconstrained) bounds
    solution = model.optimize()
    print(f"  Default FBA objective value: {solution.objective_value:.4f} "
          f"(status: {solution.status})")
    if solution.status != "optimal" or solution.objective_value <= 0:
        print("  WARNING: model did not produce a positive optimal objective. "
              "Check the biomass/objective reaction before proceeding to Step 5.")
    return model


if __name__ == "__main__":
    download_model()
    model = load_and_verify()
    cobra.io.save_json_model(model, os.path.join(MODEL_DIR, "Human-GEM.json"))
    print(f"Saved JSON copy to {os.path.join(MODEL_DIR, 'Human-GEM.json')} "
          "(faster to reload than SBML in later steps).")
