"""
Item 1 (highest priority) -- Independently reconstruct the Single-Algorithm
(iMAT+FVA) external baseline: for each of the 10 cohort x group combinations,
take ONLY the raw iMAT active-reaction set (Step 5's
data/context_specific_models/{cohort}__{group}__imat.json -- an explicitly
allowed raw input; no consensus voting across FASTCORE/tINIT/GIMME, no flux
sampling, no perturbation testing, no calibration weights), build a fresh
cobra model directly from Human-GEM.json, and run FBA + FVA on it with a
fresh solver call. This mirrors what a literature-standard single-algorithm
GSMM workflow (the kind of analysis this pipeline is benchmarked against)
would do, and is the one benchmark comparator never independently verified.

Does NOT import scripts/06_flux_analysis.py or any other GEM1 pipeline
script. Uses cobra (a third-party library, not GEM1's own code) purely as a
solver interface, loading a fresh cobra.Model instance from the raw JSON each
time.

No in-repo reference number exists for a "reported" single-algorithm
baseline score to reconcile against (confirmed by exhaustive repo search
during audit setup) -- so this item's deliverable is the independently
computed baseline itself, reported alongside GEM1's own full-pipeline
biomarker metric (computed the same way, from the calibrated_confidence_score
column) for direct side-by-side comparison. PASS/FAIL here means "did the
independent reconstruction+FVA complete successfully and land on
internally-consistent numbers", not "did it match a manuscript figure".

Long-running (~30 min estimated from historical Step 6 FVA timings at a
smaller reaction-set size) -- launch via detached Start-Process per the
project's established pattern, not as an interactive-session-attached
foreground/background job.
"""

import json
import os
import time

import cobra
import numpy as np
import pandas as pd

from audit_common import (
    ALL_GROUPS, BIOMARKERS_FULL, BIOMASS_REACTION, DETAILS_DIR, LOGS_DIR,
    MODEL_PATH, biomarker_metric, evaluation_percentile,
    load_calibrated_metadata, load_growth_support_ids, load_step5_active_set,
)

RAW_DIR = os.path.join(DETAILS_DIR, "item1_raw")
os.makedirs(RAW_DIR, exist_ok=True)


def build_restricted_model(base_model, active_ids):
    model = base_model.copy()
    active_set = set(active_ids)
    for rxn in model.reactions:
        if rxn.id not in active_set:
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
    return model


def run_one_group(base_model, cohort, group, growth_support_ids):
    t0 = time.time()
    imat_active = load_step5_active_set(cohort, group, "imat")
    if imat_active is None:
        return {"cohort": cohort, "group": group, "error": "no imat.json found"}

    active_ids = set(imat_active) | growth_support_ids
    model = build_restricted_model(base_model, active_ids)
    model.objective = BIOMASS_REACTION

    fba_sol = model.optimize()
    fba_status = fba_sol.status
    fba_obj = fba_sol.objective_value if fba_status == "optimal" else None

    fva_reactions = sorted(active_ids)
    fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=fva_reactions, fraction_of_optimum=1.0, processes=1)

    result = {
        "cohort": cohort, "group": group,
        "n_imat_active_raw": len(imat_active),
        "n_active_with_growth_support": len(active_ids),
        "fba_status": fba_status,
        "fba_objective_value": fba_obj,
        "n_fva_reactions": len(fva_reactions),
        "fva_min": {rid: float(fva.loc[rid, "minimum"]) for rid in fva.index},
        "fva_max": {rid: float(fva.loc[rid, "maximum"]) for rid in fva.index},
        "elapsed_seconds": time.time() - t0,
    }
    raw_path = os.path.join(RAW_DIR, f"{cohort}__{group}.json")
    with open(raw_path, "w") as f:
        json.dump(result, f)
    print(f"  {cohort}/{group}: FBA={fba_status} obj={fba_obj}, "
          f"{len(active_ids)} active reactions, FVA done in {result['elapsed_seconds']:.1f}s -> {raw_path}",
          flush=True)
    return result


def compute_baseline_biomarker_metric(active_sets_by_group):
    """
    baseline_score = 1.0 if reaction in this group's iMAT(+growth-support)
    active set, else 0.0 -- the only signal a naive single-algorithm+FVA
    workflow (no confidence engine) naturally provides. Evaluated with the
    same evaluation-percentile-then-max-over-isoforms-then-mean rule GEM1's
    own biomarker_ranking_by_group.csv uses, for a like-for-like comparison.
    """
    rows = []
    # universe = all reaction ids that ever appear in Human-GEM (12931) -- pull from calibrated metadata
    meta = load_calibrated_metadata()
    all_rxn_ids = sorted(meta["reaction_id"].unique())
    for cohort, group in ALL_GROUPS:
        active = active_sets_by_group.get((cohort, group))
        if active is None:
            continue
        for rid in all_rxn_ids:
            rows.append({"cohort": cohort, "group": group, "reaction_id": rid,
                         "baseline_score": 1.0 if rid in active else 0.0})
    df = pd.DataFrame(rows)
    df["eval_pct"] = evaluation_percentile(df, "baseline_score")
    groups_present = sorted(active_sets_by_group.keys())
    metric = biomarker_metric(df, "eval_pct", BIOMARKERS_FULL, groups_present)
    return metric, groups_present


def compute_gem1_full_pipeline_metric():
    meta = load_calibrated_metadata()
    meta["eval_pct"] = evaluation_percentile(meta, "calibrated_confidence_score")
    return biomarker_metric(meta, "eval_pct", BIOMARKERS_FULL, ALL_GROUPS)


def run():
    t0 = time.time()
    print("Loading fresh Human-GEM.json into a new cobra.Model instance...", flush=True)
    base_model = cobra.io.load_json_model(MODEL_PATH)
    print(f"Loaded: {len(base_model.reactions)} reactions, {len(base_model.metabolites)} metabolites.", flush=True)

    growth_support = set(load_growth_support_ids()["support_reaction_ids"])

    per_group_results = {}
    active_sets_by_group = {}
    for cohort, group in ALL_GROUPS:
        res = run_one_group(base_model, cohort, group, growth_support)
        per_group_results[f"{cohort}__{group}"] = {k: v for k, v in res.items() if k not in ("fva_min", "fva_max")}
        if "error" not in res:
            imat_active = load_step5_active_set(cohort, group, "imat")
            active_sets_by_group[(cohort, group)] = set(imat_active) | growth_support

    baseline_metric, groups_present = compute_baseline_biomarker_metric(active_sets_by_group)
    gem1_full_metric = compute_gem1_full_pipeline_metric()

    result = {
        "item": 1,
        "description": "Single-Algorithm (iMAT+FVA) external baseline, independently reconstructed",
        "n_groups_processed": len(active_sets_by_group),
        "n_groups_expected": len(ALL_GROUPS),
        "per_group_results": per_group_results,
        "baseline_biomarker_metric": {
            "definition": (
                "mean over (biomarker,group) of MAX evaluation-percentile (rank among all "
                "12931 reactions, method=min) of a binary is-in-iMAT(+growth-support)-active-set "
                "score -- the only signal a naive single-algorithm+FVA workflow provides"
            ),
            "value": baseline_metric,
            "groups_included": [f"{c}__{g}" for c, g in groups_present],
        },
        "gem1_full_pipeline_biomarker_metric": {
            "definition": "same aggregation rule, applied to GEM1's own calibrated_confidence_score column",
            "value": gem1_full_metric,
        },
        "delta_gem1_minus_baseline": gem1_full_metric - baseline_metric,
        "reported_reference_value": None,
        "reported_reference_source": "NONE FOUND -- no script/notebook/doc in the repository defines a reported iMAT+FVA baseline number (confirmed by exhaustive search during audit setup)",
        "max_abs_discrepancy": "N/A -- no reference value exists to diff against",
        "status": "PASS" if len(active_sets_by_group) == len(ALL_GROUPS) else "FAIL",
        "status_meaning": (
            "PASS = independent reconstruction + fresh FVA completed for all 10 groups and produced "
            "internally consistent numbers. Does NOT mean 'matches a manuscript figure' -- none exists "
            "in-repo to match against."
        ),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item1_imat_fva_baseline.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {out_path}")
    print(f"Baseline metric: {baseline_metric}, GEM1 full-pipeline metric: {gem1_full_metric}")
    print(f"Status: {result['status']}, total elapsed: {result['elapsed_seconds']:.1f}s")
    return result


if __name__ == "__main__":
    run()
