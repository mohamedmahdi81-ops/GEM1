"""
GEM1 - Step 4: Pull the three independent NAFLD/MASLD liver-biopsy cohorts
selected for the first case study, and produce a clean gene x sample
expression matrix for each -- ready to feed into Step 5's context-specific
extraction.

Cohorts selected (see GEM1-novelty-dossier.md for the selection rationale --
all three appear in the independent MASLD-progression meta-analysis by
Life Science Alliance 2024, https://doi.org/10.26508/lsa.202302517, which
gives external validation that these are field-standard, high-quality series):

  1. GSE89632  - Illumina microarray, 63 samples: healthy (24) / simple
                 steatosis (20) / NASH (19). Wruck et al.
  2. GSE126848 - RNA-seq, 57 samples: healthy normal-weight (14) / obese (12)
                 / NAFL (15) / NASH (16). Suppli et al.
  3. GSE135251 - RNA-seq, large multi-center European cohort spanning the
                 full steatosis -> NASH -> fibrosis -> cirrhosis spectrum.
                 Govaere et al., Sci Transl Med 2020. Exact current sample
                 count is confirmed at download time below -- print it and
                 record it in the dossier once you've pulled the data.

Two platforms (microarray + RNA-seq) and three genuinely independent patient
populations were chosen deliberately -- this is what "dataset agreement"
(30% of the FluxConfidence-derived GEM1 score) is measured across.

Requires: GEOparse, pandas   (pip install GEOparse pandas)
Run on a machine with internet access.
"""

import os
import GEOparse
import pandas as pd

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "geo_cohorts")
os.makedirs(OUT_DIR, exist_ok=True)

COHORTS = ["GSE89632", "GSE126848", "GSE135251"]


def fetch_and_summarize(gse_id: str):
    print(f"\n=== {gse_id} ===")
    gse = GEOparse.get_GEO(geo=gse_id, destdir=OUT_DIR)

    # Build a samples-x-metadata table so you can define disease-group
    # labels (healthy / steatosis / NASH / fibrosis) for Step 5's
    # context-specific extraction per group.
    meta_rows = []
    for gsm_name, gsm in gse.gsms.items():
        meta_rows.append({
            "sample": gsm_name,
            "title": gsm.metadata.get("title", [""])[0],
            "characteristics": " | ".join(gsm.metadata.get("characteristics_ch1", [])),
        })
    meta_df = pd.DataFrame(meta_rows)
    meta_path = os.path.join(OUT_DIR, f"{gse_id}_sample_metadata.csv")
    meta_df.to_csv(meta_path, index=False)
    print(f"  Samples: {len(meta_df)}")
    print(f"  Sample metadata saved to {meta_path}")
    print("  --> Open this file and assign each sample to a disease group "
          "(healthy/steatosis/NASH/fibrosis) before Step 5.")

    # Expression matrix (probes/genes x samples). For RNA-seq series GEOparse
    # pulls the supplementary processed matrix if GEO provides one; for
    # microarray series it assembles from the platform table + GSM values.
    try:
        expr_df = gse.pivot_samples("VALUE")
        expr_path = os.path.join(OUT_DIR, f"{gse_id}_expression_matrix.csv")
        expr_df.to_csv(expr_path)
        print(f"  Expression matrix ({expr_df.shape[0]} features x "
              f"{expr_df.shape[1]} samples) saved to {expr_path}")
    except Exception as e:
        print(f"  Could not auto-build an expression matrix for {gse_id} "
              f"from GSM VALUE fields ({e}). This is common for RNA-seq "
              "series that publish counts only as a supplementary file -- "
              "check gse.metadata['supplementary_file'] and download that "
              "directly instead.")
        print(f"  Supplementary files listed: "
              f"{gse.metadata.get('supplementary_file', 'none')}")

    return gse, meta_df


if __name__ == "__main__":
    results = {}
    for gse_id in COHORTS:
        results[gse_id] = fetch_and_summarize(gse_id)

    print("\nAll cohorts fetched. Next: label each sample's disease group in "
          "the *_sample_metadata.csv files, then run 03_build_context_specific_models.py")
