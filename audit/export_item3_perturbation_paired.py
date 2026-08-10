"""
Export the full paired reaction-level perturbation_robustness_raw dataset
(original pipeline run vs. item 3's independent re-execution) for building a
real scatter plot -- no subsampling, no aggregation.

Reads only audit/details/item3_per_reaction_diffs.csv (already computed by
item3_flux_sampling_perturbation.py's own comparison step -- not recomputed
here) and data/context_specific_models/minimal_growth_support.json (raw
pipeline output) for the is_growth_support flag. Does not modify any
existing GEM1 output file; writes only to audit/details/.
"""

import os

import pandas as pd

from audit_common import DETAILS_DIR, load_growth_support_ids

diffs_path = os.path.join(DETAILS_DIR, "item3_per_reaction_diffs.csv")
out_path = os.path.join(DETAILS_DIR, "item3_perturbation_paired_values.csv")

df = pd.read_csv(diffs_path)
pert = df[df["metric"] == "perturbation_robustness_raw"].copy()

growth_support_ids = set(load_growth_support_ids()["support_reaction_ids"])
pert["is_growth_support"] = pert["reaction_id"].isin(growth_support_ids)

out = pert.rename(columns={
    "independent": "perturbation_robustness_raw_reexecution",
    "existing": "perturbation_robustness_raw_original",
})[[
    "reaction_id", "cohort", "group",
    "perturbation_robustness_raw_original", "perturbation_robustness_raw_reexecution",
    "is_growth_support",
]]

out.to_csv(out_path, index=False)

print(f"Wrote {out_path}")
print(f"Rows: {len(out)}")
print(f"Growth-support rows: {out['is_growth_support'].sum()}")
print(f"Non-growth-support rows: {(~out['is_growth_support']).sum()}")
print(out.head())
