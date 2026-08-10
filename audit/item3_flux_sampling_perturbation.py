"""
Item 3 -- Independently re-execute flux sampling and perturbation robustness
on the consensus-active reaction sets, and compare against flux_confidence_raw
/ perturbation_robustness_raw in all_groups_calibrated_with_metadata.csv using
distribution/CI-based tolerance (not exact-match -- both steps are stochastic).

Consensus active sets are recomputed independently from the raw Step 5 JSONs
(audit_common.compute_consensus_from_raw: set intersection over
fastcore/imat/tinit, mirroring 04_consensus_scoring.py's documented
*definition*, not its code) unioned with the growth-support set -- NOT read
from data/consensus_scores/*.json.

Formulas reimplemented from scratch here (re-derived from the plain-English
trace of 09_calibration.py's docstring during audit setup, not copy-pasted):
  flux_confidence_raw = clip(1 - (p95-p5)/(fva_max-fva_min), 0, 1); 1.0 if
    fva_max==fva_min; NaN if not active.
  perturbation_robustness_score(R) = P(feasible | R present) - P(feasible | R absent)
  perturbation_robustness_raw = (score+1)/2, except growth-support reactions
    fixed at 1.0 (never eligible for removal); NaN if not active.

Methodology parameters (n_samples, thinning, seed, obj_frac, n_iterations,
magnitude, feasibility obj_frac) are reused from the documented protocol --
these define WHAT was run, not GEM1's own scoring arithmetic, so reusing them
is what makes this a check of the *pipeline's execution*, not a different
experiment.

Long-running (~75 min estimated from historical Step 7+8 timings) -- launch
via detached Start-Process, not a Claude-Code-attached job.
"""

import json
import os
import time

import cobra
import numpy as np
import pandas as pd
from cobra.sampling import OptGPSampler

from audit_common import (
    ALL_GROUPS, DETAILS_DIR, MODEL_PATH, BIOMASS_REACTION,
    compute_consensus_from_raw, load_calibrated_metadata, load_growth_support_ids,
)

RAW_DIR = os.path.join(DETAILS_DIR, "item3_raw")
os.makedirs(RAW_DIR, exist_ok=True)

N_SAMPLES = 500
THINNING = 100
SEED = 42
SAMPLING_OBJ_FRAC = 0.9
N_PERTURBATION_ITER = 200
PERTURBATION_MAGNITUDE = 0.3
PERTURBATION_OBJ_FRAC = 0.99


def build_bounds_zeroed_model(base_model, active_ids):
    model = base_model.copy()
    active_set = set(active_ids)
    for rxn in model.reactions:
        if rxn.id not in active_set:
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
    return model


def build_reduced_model(base_model, active_ids):
    model = base_model.copy()
    active_set = set(active_ids)
    to_remove = [r for r in model.reactions if r.id not in active_set]
    model.remove_reactions(to_remove, remove_orphans=True)
    return model


def run_fva_fraction1(model, reaction_ids):
    model.solver = "gurobi"
    return cobra.flux_analysis.flux_variability_analysis(model, reaction_list=sorted(reaction_ids), fraction_of_optimum=1.0, processes=1)


def run_sampling(reduced_model, fba_optimum):
    reduced_model.solver = "glpk"
    biomass_rxn = reduced_model.reactions.get_by_id(BIOMASS_REACTION)
    floor = SAMPLING_OBJ_FRAC * fba_optimum
    biomass_rxn.lower_bound = max(biomass_rxn.lower_bound, floor)
    sampler = OptGPSampler(reduced_model, processes=1, thinning=THINNING, seed=SEED)
    samples = sampler.sample(N_SAMPLES)
    return samples


def summarize_samples(samples: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "mean": samples.mean(), "std": samples.std(), "min": samples.min(),
        "p5": samples.quantile(0.05), "median": samples.median(),
        "p95": samples.quantile(0.95), "max": samples.max(),
    })


def run_perturbation(bz_model, candidate_ids, target_biomass, rng):
    bz_model.solver = "glpk"
    candidate_ids = sorted(candidate_ids)
    n = len(candidate_ids)
    removed_mask = rng.random((N_PERTURBATION_ITER, n)) < PERTURBATION_MAGNITUDE
    feasible = np.zeros(N_PERTURBATION_ITER, dtype=bool)

    original_bounds = {rid: (bz_model.reactions.get_by_id(rid).lower_bound,
                              bz_model.reactions.get_by_id(rid).upper_bound) for rid in candidate_ids}

    for it in range(N_PERTURBATION_ITER):
        removed = [candidate_ids[i] for i in range(n) if removed_mask[it, i]]
        for rid in removed:
            rxn = bz_model.reactions.get_by_id(rid)
            rxn.lower_bound, rxn.upper_bound = 0.0, 0.0
        sol = bz_model.optimize()
        feasible[it] = sol.status == "optimal" and sol.objective_value is not None and sol.objective_value >= target_biomass
        for rid in removed:
            rxn = bz_model.reactions.get_by_id(rid)
            rxn.lower_bound, rxn.upper_bound = original_bounds[rid]

    scores = {}
    for i, rid in enumerate(candidate_ids):
        present_mask = ~removed_mask[:, i]
        absent_mask = removed_mask[:, i]
        rate_present = feasible[present_mask].mean() if present_mask.any() else np.nan
        rate_absent = feasible[absent_mask].mean() if absent_mask.any() else np.nan
        scores[rid] = {
            "n_trials_present": int(present_mask.sum()), "n_trials_absent": int(absent_mask.sum()),
            "feasibility_rate_when_present": float(rate_present) if not np.isnan(rate_present) else None,
            "feasibility_rate_when_absent": float(rate_absent) if not np.isnan(rate_absent) else None,
            "perturbation_robustness_score": float(rate_present - rate_absent) if not (np.isnan(rate_present) or np.isnan(rate_absent)) else None,
        }
    return scores, feasible


def run_one_group(base_model, cohort, group, growth_support_ids, seed_offset):
    t0 = time.time()
    consensus = compute_consensus_from_raw(cohort, group)
    strict = consensus["strict_consensus"]
    active_ids = strict | growth_support_ids
    print(f"  {cohort}/{group}: strict_consensus={len(strict)}, active(+growth support)={len(active_ids)}", flush=True)

    # fresh FVA (fraction_of_optimum=1.0) over bounds-zeroed model, mirrors Step 6's active-set convention
    bz_model = build_bounds_zeroed_model(base_model, active_ids)
    bz_model.objective = BIOMASS_REACTION
    fva = run_fva_fraction1(bz_model, active_ids)

    # fresh FBA optimum on the reduced model, for the sampling floor
    reduced_model = build_reduced_model(base_model, active_ids)
    reduced_model.objective = BIOMASS_REACTION
    reduced_model.solver = "glpk"
    fba_sol = reduced_model.optimize()
    fba_optimum = fba_sol.objective_value
    print(f"    fresh FBA optimum: {fba_optimum}", flush=True)

    samples = run_sampling(reduced_model, fba_optimum)
    summary = summarize_samples(samples)
    print(f"    sampling done: {samples.shape[0]} samples x {samples.shape[1]} reactions in {time.time()-t0:.1f}s", flush=True)

    # perturbation: candidates = strict consensus minus growth support (never eligible for removal)
    candidates = strict - growth_support_ids
    target_biomass = PERTURBATION_OBJ_FRAC * fba_optimum
    rng = np.random.default_rng(SEED + seed_offset)
    bz_model2 = build_bounds_zeroed_model(base_model, active_ids)
    bz_model2.objective = BIOMASS_REACTION
    pert_scores, feasible = run_perturbation(bz_model2, candidates, target_biomass, rng)
    print(f"    perturbation done: {len(candidates)} candidates, "
          f"{feasible.mean():.3f} overall feasibility rate, total {time.time()-t0:.1f}s", flush=True)

    # --- derive flux_confidence_raw and perturbation_robustness_raw ---
    flux_confidence_raw = {}
    for rid in active_ids:
        if rid not in fva.index or rid not in summary.index:
            continue
        fva_min, fva_max = fva.loc[rid, "minimum"], fva.loc[rid, "maximum"]
        if fva_max == fva_min:
            flux_confidence_raw[rid] = 1.0
        else:
            p5, p95 = summary.loc[rid, "p5"], summary.loc[rid, "p95"]
            val = 1.0 - (p95 - p5) / (fva_max - fva_min)
            flux_confidence_raw[rid] = float(np.clip(val, 0.0, 1.0))

    perturbation_robustness_raw = {}
    for rid in active_ids:
        if rid in growth_support_ids:
            perturbation_robustness_raw[rid] = 1.0
        elif rid in pert_scores and pert_scores[rid]["perturbation_robustness_score"] is not None:
            perturbation_robustness_raw[rid] = (pert_scores[rid]["perturbation_robustness_score"] + 1.0) / 2.0

    result = {
        "cohort": cohort, "group": group,
        "n_strict_consensus": len(strict), "n_active": len(active_ids),
        "fba_optimum": fba_optimum,
        "flux_confidence_raw": flux_confidence_raw,
        "perturbation_robustness_raw": perturbation_robustness_raw,
        "perturbation_overall_feasibility_rate": float(feasible.mean()),
        "elapsed_seconds": time.time() - t0,
    }
    raw_path = os.path.join(RAW_DIR, f"{cohort}__{group}.json")
    with open(raw_path, "w") as f:
        json.dump(result, f)
    print(f"  -> {raw_path}", flush=True)
    return result


def compare_to_existing(all_results):
    existing = load_calibrated_metadata()
    rows = []
    for key, res in all_results.items():
        cohort, group = res["cohort"], res["group"]
        sub = existing[(existing["cohort"] == cohort) & (existing["group"] == group)].set_index("reaction_id")
        for rid, new_val in res["flux_confidence_raw"].items():
            if rid in sub.index:
                old_val = sub.loc[rid, "flux_confidence_raw"]
                if pd.notna(old_val):
                    rows.append({"cohort": cohort, "group": group, "reaction_id": rid, "metric": "flux_confidence_raw",
                                 "independent": new_val, "existing": float(old_val), "diff": new_val - float(old_val)})
        for rid, new_val in res["perturbation_robustness_raw"].items():
            if rid in sub.index:
                old_val = sub.loc[rid, "perturbation_robustness_raw"]
                if pd.notna(old_val):
                    rows.append({"cohort": cohort, "group": group, "reaction_id": rid, "metric": "perturbation_robustness_raw",
                                 "independent": new_val, "existing": float(old_val), "diff": new_val - float(old_val)})
    return pd.DataFrame(rows)


def tolerance_summary(diffs: pd.DataFrame):
    out = {}
    for metric, sub in diffs.groupby("metric"):
        d = sub["diff"].values
        corr = np.corrcoef(sub["independent"], sub["existing"])[0, 1] if len(sub) > 1 else float("nan")
        spearman = pd.Series(sub["independent"]).corr(pd.Series(sub["existing"]), method="spearman")
        out[metric] = {
            "n": len(sub),
            "mean_diff": float(np.mean(d)),
            "std_diff": float(np.std(d)),
            "mean_abs_diff": float(np.mean(np.abs(d))),
            "max_abs_diff": float(np.max(np.abs(d))),
            "pct_within_0.15": float(np.mean(np.abs(d) <= 0.15)),
            "pct_within_0.25": float(np.mean(np.abs(d) <= 0.25)),
            "pearson_r": float(corr),
            "spearman_r": float(spearman),
            "ci95_mean_diff": [
                float(np.mean(d) - 1.96 * np.std(d) / np.sqrt(len(d))),
                float(np.mean(d) + 1.96 * np.std(d) / np.sqrt(len(d))),
            ],
        }
    return out


def run():
    t0 = time.time()
    print("Loading fresh Human-GEM.json into a new cobra.Model instance...", flush=True)
    base_model = cobra.io.load_json_model(MODEL_PATH)
    growth_support = set(load_growth_support_ids()["support_reaction_ids"])

    all_results = {}
    for i, (cohort, group) in enumerate(ALL_GROUPS):
        try:
            res = run_one_group(base_model, cohort, group, growth_support, seed_offset=i)
            all_results[f"{cohort}__{group}"] = res
        except Exception as e:
            print(f"  {cohort}/{group} FAILED: {e}", flush=True)
            all_results[f"{cohort}__{group}"] = {"cohort": cohort, "group": group, "error": str(e)}

    ok_results = {k: v for k, v in all_results.items() if "error" not in v}
    diffs = compare_to_existing(ok_results)
    tol_summary = tolerance_summary(diffs) if len(diffs) else {}

    # PASS criterion: distribution-based, not exact-match -- Spearman r >= 0.5 and
    # mean_abs_diff <= 0.25 for both metrics (generous stochastic-reproduction bar,
    # since sampling/perturbation are seeded but not bit-identical across independent code).
    status = "PASS"
    for metric in ("flux_confidence_raw", "perturbation_robustness_raw"):
        s = tol_summary.get(metric)
        if s is None or s["spearman_r"] < 0.5 or s["mean_abs_diff"] > 0.25:
            status = "FAIL"

    result = {
        "item": 3,
        "description": "Flux sampling + perturbation robustness independently re-executed on independently-recomputed consensus-active reaction sets",
        "methodology_params": {
            "n_samples": N_SAMPLES, "thinning": THINNING, "seed": SEED,
            "sampling_obj_frac": SAMPLING_OBJ_FRAC, "n_perturbation_iter": N_PERTURBATION_ITER,
            "perturbation_magnitude": PERTURBATION_MAGNITUDE, "perturbation_obj_frac": PERTURBATION_OBJ_FRAC,
        },
        "n_groups_processed": len(ok_results), "n_groups_expected": len(ALL_GROUPS),
        "n_groups_failed": len(all_results) - len(ok_results),
        "tolerance_summary": tol_summary,
        "tolerance_criteria": "PASS requires Spearman r >= 0.5 AND mean_abs_diff <= 0.25 for both flux_confidence_raw and perturbation_robustness_raw (distribution/CI-based, not exact-match)",
        "max_abs_discrepancy": max((s["max_abs_diff"] for s in tol_summary.values()), default=None),
        "status": status,
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item3_flux_sampling_perturbation.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    diffs.to_csv(os.path.join(DETAILS_DIR, "item3_per_reaction_diffs.csv"), index=False)
    print(f"\nWrote {out_path}")
    print(json.dumps(tol_summary, indent=2))
    print(f"Status: {status}, total elapsed: {result['elapsed_seconds']:.1f}s")
    return result


if __name__ == "__main__":
    run()
