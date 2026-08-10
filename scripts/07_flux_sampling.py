"""
GEM1 - Step 7: Probabilistic flux sampling.

Purpose: characterize the *distribution* of feasible flux states for each
group's context-specific model, not just a single FBA/pFBA point solution
(Step 6) -- this is the uncertainty-quantification layer the confidence
engine's reserved `flux_sampling_uncertainty_score` (Step 10, Layer 3) will
eventually consume. Does not touch Step 9 (structural consensus) or Step 10
(confidence engine) in any way; both remain exactly as already computed and
locked.

Methodology (per Mohammed's instructions, 2026-07-18):

1. Model scaffold reused verbatim from Step 6, not reimplemented. This
   script imports `load_base_model`, `load_growth_support`,
   `load_strict_consensus`, and `build_constrained_model` directly from
   `06_flux_analysis.py` (via importlib, since the module name starts with a
   digit) rather than re-deriving the active-reaction-set logic -- so the
   exact same reaction set (Step 9's strict_consensus_reactions unioned with
   the 519-reaction growth-support set, see CLAUDE.md's Step 6 methodology
   section) is used here. Step 6 is a functional scaffold consumed by this
   step, not re-derived or reinterpreted.

2. Biomass floor, not a single point. Sampling the model exactly at its FBA
   optimum would (generically) collapse the feasible region to a single
   point or a thin degenerate slice, since MAR13082 would be pinned at
   exactly one value. Instead this script re-solves FBA on the Step 6
   scaffold (cross-checked against Step 6's own stored fba_objective_value,
   see validate_step6_inputs() below) and imposes `MAR13082 >= obj_frac *
   fba_optimum` as a lower-bound constraint before sampling -- the same
   obj_frac=0.9 convention already used throughout this pipeline (GIMME's
   own obj_frac, the minimal_growth_support derivation). This opens up a
   genuine, non-degenerate polytope of "reasonably growing" flux states to
   sample from, while still respecting the group's context-specific active-
   reaction set.

3. Sampler: cobra's OptGPSampler (hit-and-run, parallelizable, cobra's
   recommended sampler for genome-scale models over ACHR). Samples are drawn
   in two sequential batches (`sampler.batch(n_samples // 2, 2)`) from the
   SAME running chain -- not two independent chains -- so comparing batch 1
   vs batch 2 is a first-half-vs-second-half convergence/stationarity check
   in the standard MCMC sense (see validate_samples() below), not a
   between-chain agreement check.

4. Reproducibility: every run is seeded (`--seed`, default 42) and every
   sampling parameter (n_samples, thinning, processes, obj_frac) is a CLI
   flag with a documented default, not hardcoded.

Validation performed on every group before results are trusted (see
validate_samples()):
  a. Structural feasibility: cobra's own `sampler.validate()` on every drawn
     sample (checks S.v = 0 equality and bound feasibility per sample).
  b. Range containment: every sampled flux, for every reaction, falls within
     that reaction's FVA-derived [min, max] range under the SAME obj_frac
     constraint used for sampling (FVA re-run fresh here, not reused from
     Step 6, since Step 6's FVA was under a different constraint --
     fraction_of_optimum=1.0 pinning biomass near its exact optimum -- not
     this step's looser obj_frac=0.9 floor).
  c. Convergence: two-sample Kolmogorov-Smirnov test per reaction, first
     batch vs second batch of the same chain; reports the fraction of
     reactions with a significant (p < 0.01) distributional shift as a
     stationarity diagnostic. Some fraction of false positives is expected
     by chance at this alpha across thousands of reactions -- this is a
     diagnostic to review, not an automatic pass/fail gate (documented
     explicitly in the output rather than silently thresholded).

Outputs per group in data/flux_sampling/:
  {cohort}__{group}__samples.csv.gz -- full raw sample matrix (n_samples x
    n_active_reactions), for reproducibility and any future re-analysis.
  {cohort}__{group}__summary.csv -- per-reaction mean/std/median/min/max/
    p5/p95 across samples -- the primary artifact downstream steps should
    consume (Step 8 perturbation baseline, Step 10 Layer 3 eventually).
  {cohort}__{group}__validation.json -- the three validation results above.
flux_sampling_manifest.csv -- one row per group: parameters used, runtime,
  and a rollup of the validation results.

Run from gem1-main (GEM1 conda env) -- standard cobra (0.31.1), Gurobi as
the LP solver (same as Step 6, no MILP in this step).
"""

import os
import sys
import json
import time
import argparse
import importlib.util
import numpy as np
import pandas as pd
import cobra
from cobra.sampling import OptGPSampler
from scipy.stats import ks_2samp

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.dirname(__file__)
CONSENSUS_DIR = os.path.join(BASE_DIR, "data", "consensus_scores")
FLUX_DIR = os.path.join(BASE_DIR, "data", "flux_analysis")
OUT_DIR = os.path.join(BASE_DIR, "data", "flux_sampling")
os.makedirs(OUT_DIR, exist_ok=True)

BIOMASS_ID = "MAR13082"


def _import_step06():
    """Load 06_flux_analysis.py as a module (name starts with a digit, so a
    plain `import` isn't possible) to reuse its scaffold functions verbatim
    rather than duplicating/re-deriving them."""
    path = os.path.join(SCRIPTS_DIR, "06_flux_analysis.py")
    spec = importlib.util.spec_from_file_location("step06_flux_analysis", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_step6_inputs(step06, cohort, group, model, active_set):
    """Cross-check this script's freshly-solved FBA against Step 6's own
    stored fba_objective_value for the same reaction set, before sampling --
    confirms Step 6's outputs are still valid/reproducible inputs, per
    Mohammed's instruction to verify inputs before implementing Step 7."""
    flux_json_path = os.path.join(FLUX_DIR, f"{cohort}__{group}__flux.json")
    with open(flux_json_path) as f:
        step6_result = json.load(f)

    with model:
        model.objective = BIOMASS_ID
        sol = model.optimize()
        fresh_fba = sol.objective_value if sol.status == "optimal" else None

    stored_fba = step6_result["fba_objective_value"]
    match = fresh_fba is not None and abs(fresh_fba - stored_fba) < 1e-4
    stored_active = step6_result["n_active_reactions"]
    active_match = len(active_set) == stored_active

    return {
        "step6_fba_objective_value": stored_fba,
        "freshly_solved_fba_objective_value": fresh_fba,
        "fba_reproduced": match,
        "step6_n_active_reactions": stored_active,
        "this_step_n_active_reactions": len(active_set),
        "active_set_size_matches": active_match,
    }, fresh_fba


def build_sampling_model(base_model, active_set, growth_support, fba_optimum, obj_frac, step06):
    """
    2026-07-19: solver switched to GLPK here (applies to both the FVA and
    the sampling that follow), not just inside the sampler -- see
    draw_samples()'s docstring for the full diagnosis. Every LP solve this
    model is used for downstream (FVA re-solves, OptGP's warmup points and
    hit-and-run steps) is pure LP, so GLPK avoids Gurobi's per-solve WLS
    network handshake entirely without any loss of correctness.

    Superseded for the actual sampling model by build_reduced_sampling_model()
    below -- kept here only because validate_step6_inputs() still needs a
    bounds-zeroed (not reaction-removed) model to directly reproduce Step 6's
    own stored FBA value under Step 6's own exact representation.
    """
    model, n_blocked, unblocked = step06.build_constrained_model(base_model, active_set, growth_support)
    model.solver = "glpk"
    model.objective = BIOMASS_ID
    floor = obj_frac * fba_optimum
    biomass_rxn = model.reactions.get_by_id(BIOMASS_ID)
    biomass_rxn.lower_bound = max(biomass_rxn.lower_bound, floor)
    return model, unblocked, floor


def build_reduced_sampling_model(base_model, active_set, growth_support, fba_optimum, obj_frac):
    """
    2026-07-19, performance fix (confirmed with Mohammed before implementing,
    same day): a smoke test sampling the bounds-zeroed 12,931-reaction model
    (Step 6's representation, reused verbatim) took ~67 minutes for just 20
    samples at thinning=10 -- at the real parameters (500 samples,
    thinning=100, x10 groups) this would not finish in any practical time.
    OptGP's warmup-point generation and hit-and-run stepping cost scale with
    the FULL model's variable/constraint count regardless of how many
    reactions are pinned to (0, 0); Step 6 deliberately zeroed bounds rather
    than removing reactions to preserve stoichiometric structure for
    FBA/pFBA/FVA (cheap regardless of model size there), but that same
    property makes sampling pay full genome-scale cost to explore what is
    actually only an ~4,700-dimensional space.

    Fix, scoped to the sampling model only: reactions NOT in the active set
    (strict consensus | growth support) are REMOVED via cobra's
    `remove_reactions(..., remove_orphans=True)` rather than bounds-zeroed,
    for this step's own model only -- Step 6's stored outputs, its own
    model-building function, and every other locked step are untouched.
    This does not change the feasible region for any surviving reaction:
    a reaction already pinned to (0, 0) contributes nothing but an
    always-satisfiable zero-valued variable to the LP, so removing it cannot
    change the optimal value or feasible range of any reaction that remains.
    Confirmed empirically in process_group(), which cross-checks FBA on this
    reduced model against the bounds-zeroed model's freshly-solved value
    (itself already checked against Step 6's stored value) before sampling
    proceeds -- if the two representations disagree, that's a real
    correctness bug, not an expected consequence of this optimization, and
    the run stops rather than sampling on unreconciled inputs.
    """
    unblocked = active_set | growth_support
    model = base_model.copy()
    to_remove = [r for r in model.reactions if r.id not in unblocked]
    model.remove_reactions(to_remove, remove_orphans=True)
    model.solver = "glpk"
    model.objective = BIOMASS_ID
    floor = obj_frac * fba_optimum
    biomass_rxn = model.reactions.get_by_id(BIOMASS_ID)
    biomass_rxn.lower_bound = max(biomass_rxn.lower_bound, floor)
    return model, unblocked, floor


def run_fva_for_sampling(model, unblocked):
    """
    2026-07-19, bug found via the FVA-range-violation validation check
    itself (not by inspection): cobra's flux_variability_analysis() defaults
    to fraction_of_optimum=1.0, which adds its OWN constraint pinning the
    model's objective (MAR13082) near its exact optimum -- on top of the
    manual `biomass_rxn.lower_bound = floor` (obj_frac=0.9 of optimum)
    already set in build_reduced_sampling_model(). That stacked FVA under a
    much tighter box (biomass ~124.84) than what the sampler actually
    explores (biomass anywhere in [0.9*optimum, optimum] = [~112.35,
    ~124.84]), producing an internally-inconsistent comparison: FVA reported
    near-zero ranges for reactions that the wider sampled region legitimately
    uses. Confirmed directly -- several "violations" had a sampled value of
    636 against an FVA range of essentially [0, 0], which is only possible
    if FVA and the sampler were solving under different constraints, not a
    sampler defect. Fix: fraction_of_optimum=0.0 tells FVA not to add any
    additional objective-fraction constraint beyond the model's own existing
    bounds, so it respects exactly the manually-set floor -- the same
    feasible region the sampler operates in.
    """
    from cobra.flux_analysis import flux_variability_analysis
    return flux_variability_analysis(model, reaction_list=sorted(unblocked), fraction_of_optimum=0.0)


def draw_samples(model, n_samples, thinning, processes, seed):
    """
    2026-07-19: OptGP's warmup-point generation plus every hit-and-run step
    each issue their own LP solve. With Gurobi as the model's solver (as
    inherited from Step 6's scaffold, needed there for MILP-adjacent
    performance), each of those solves re-authenticates a fresh Gurobi WLS
    (Web License Service) environment over the network -- confirmed directly
    from a smoke-test log showing dozens of repeated "Read LP format model
    from file" + "Set parameter WLSAccessID/WLSSecret" + "Academic license"
    blocks before even the first batch of 20 samples completed, reproduced
    identically with processes=1 (ruling out multiprocessing/WLS-concurrency
    contention as the cause -- it's the per-solve WLS handshake itself).
    OptGP's own solves are pure LP (hit-and-run within a fixed polytope, no
    MILP), so Gurobi's MILP strength is not needed here. The model's solver
    is switched to GLPK in build_sampling_model() (applies here too, since
    OptGPSampler.__init__ does model.copy() which preserves whatever solver
    the model already has) -- a solver-backend choice, not a change to the
    sampling methodology, reaction set, or any already-locked step's output.
    """
    sampler = OptGPSampler(model, thinning=thinning, processes=processes, seed=seed)
    batch_size = n_samples // 2
    batch1 = sampler.batch(batch_size, 1).__next__()
    batch2 = sampler.batch(n_samples - batch_size, 1).__next__()
    return sampler, batch1, batch2


def validate_samples(sampler, batch1, batch2, fva_df, ks_alpha=0.01):
    all_samples = pd.concat([batch1, batch2], ignore_index=True)

    validity_codes = sampler.validate(all_samples.values)
    n_valid = int(np.sum(validity_codes == "v"))
    pct_valid = 100.0 * n_valid / len(validity_codes)

    fva_min = fva_df["minimum"]
    fva_max = fva_df["maximum"]
    tol = 1e-6
    out_of_range_counts = {}
    for rxn in all_samples.columns:
        lo, hi = fva_min.get(rxn), fva_max.get(rxn)
        if lo is None or hi is None:
            continue
        col = all_samples[rxn]
        n_out = int(((col < lo - tol) | (col > hi + tol)).sum())
        if n_out:
            out_of_range_counts[rxn] = n_out
    n_reactions_with_violations = len(out_of_range_counts)
    total_violations = sum(out_of_range_counts.values())
    pct_range_ok = 100.0 * (1 - n_reactions_with_violations / all_samples.shape[1])

    # 2026-07-19: report effect size (KS D statistic, and a Cohen's-d-like
    # standardized mean shift) alongside the raw significance count, not
    # instead of it. A smoke test showed the raw "% significant at p<0.01"
    # figure is misleading on its own -- it went UP (59.5% -> 75.8%) when
    # n_samples per batch increased 10x (20 -> 200), purely because the KS
    # test gains statistical power to detect small differences as batches
    # grow, not because convergence got worse. Effect size is what actually
    # tells you whether a "significant" shift is practically meaningful.
    ks_pvalues, ks_stats, mean_shifts = {}, {}, {}
    for rxn in all_samples.columns:
        a, b = batch1[rxn].values, batch2[rxn].values
        if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.isclose(a[0], b[0]):
            continue
        d, p = ks_2samp(a, b)
        ks_pvalues[rxn] = p
        ks_stats[rxn] = d
        pooled_std = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        mean_shifts[rxn] = abs(a.mean() - b.mean()) / pooled_std if pooled_std > 1e-9 else 0.0
    n_tested = len(ks_pvalues)
    n_significant = sum(1 for p in ks_pvalues.values() if p < ks_alpha)
    pct_significant = 100.0 * n_significant / n_tested if n_tested else 0.0
    ks_d_values = np.array(list(ks_stats.values())) if ks_stats else np.array([0.0])
    shift_values = np.array(list(mean_shifts.values())) if mean_shifts else np.array([0.0])
    n_large_shift = int(np.sum(shift_values > 0.5))  # >0.5 pooled-SD shift: a conventional "medium" effect size

    return {
        "n_samples_total": len(all_samples),
        "pct_samples_structurally_valid": round(pct_valid, 4),
        "n_reactions_with_range_violations": n_reactions_with_violations,
        "total_range_violation_count": total_violations,
        "ks_d_statistic_mean": round(float(ks_d_values.mean()), 4),
        "ks_d_statistic_median": round(float(np.median(ks_d_values)), 4),
        "ks_d_statistic_max": round(float(ks_d_values.max()), 4),
        "standardized_mean_shift_median": round(float(np.median(shift_values)), 4),
        "n_reactions_with_medium_or_larger_shift": n_large_shift,
        "pct_reactions_with_medium_or_larger_shift": round(100.0 * n_large_shift / len(shift_values), 4),
        "pct_reactions_within_fva_range": round(pct_range_ok, 4),
        "ks_alpha": ks_alpha,
        "n_reactions_ks_tested": n_tested,
        "n_reactions_ks_significant": n_significant,
        "pct_reactions_ks_significant": round(pct_significant, 4),
    }, all_samples


def summarize_samples(all_samples):
    desc = all_samples.describe(percentiles=[0.05, 0.5, 0.95]).T
    desc = desc.rename(columns={"5%": "p5", "50%": "median", "95%": "p95"})
    desc.index.name = "reaction_id"
    return desc.reset_index()[["reaction_id", "mean", "std", "min", "p5", "median", "p95", "max"]]


def process_group(step06, cohort, group, base_model, growth_support, args):
    t0 = time.time()
    active_set, consensus_def = step06.load_strict_consensus(cohort, group)
    model, n_blocked, unblocked = step06.build_constrained_model(base_model.copy(), active_set, growth_support)

    input_check, fresh_fba = validate_step6_inputs(step06, cohort, group, model, active_set)
    if not input_check["fba_reproduced"]:
        raise RuntimeError(
            f"{cohort}/{group}: freshly-solved FBA ({fresh_fba}) does not reproduce "
            f"Step 6's stored value ({input_check['step6_fba_objective_value']}) -- "
            "Step 6 inputs are not reproducible, stopping rather than sampling on "
            "possibly-stale/inconsistent data."
        )

    sampling_model, unblocked2, floor = build_reduced_sampling_model(
        base_model.copy(), active_set, growth_support, fresh_fba, args.obj_frac
    )

    with sampling_model:
        sampling_model.objective = BIOMASS_ID
        reduced_fba_sol = sampling_model.optimize()
    reduced_fba_value = reduced_fba_sol.objective_value if reduced_fba_sol.status == "optimal" else None
    reduced_matches_bounds_zeroed = (
        reduced_fba_value is not None and abs(reduced_fba_value - fresh_fba) < 1e-4
    )
    input_check["reduced_model_fba_value"] = reduced_fba_value
    input_check["reduced_model_matches_bounds_zeroed_model"] = reduced_matches_bounds_zeroed
    if not reduced_matches_bounds_zeroed:
        raise RuntimeError(
            f"{cohort}/{group}: reduced (reaction-removed) sampling model's FBA "
            f"({reduced_fba_value}) does not match the bounds-zeroed model's FBA "
            f"({fresh_fba}) at unconstrained-biomass-objective -- these representations "
            "should be mathematically equivalent; stopping rather than sampling on a "
            "model that disagrees with the already-validated Step 6 scaffold."
        )

    fva_df = run_fva_for_sampling(sampling_model, unblocked2)

    sampler, batch1, batch2 = draw_samples(
        sampling_model, args.n_samples, args.thinning, args.processes, args.seed
    )
    validation, all_samples = validate_samples(sampler, batch1, batch2, fva_df, args.ks_alpha)
    summary_df = summarize_samples(all_samples)
    runtime_s = time.time() - t0

    samples_path = os.path.join(OUT_DIR, f"{cohort}__{group}__samples.csv.gz")
    all_samples.to_csv(samples_path, index=False, compression="gzip")
    summary_path = os.path.join(OUT_DIR, f"{cohort}__{group}__summary.csv")
    summary_df.to_csv(summary_path, index=False)
    validation_payload = {
        "cohort": cohort, "group": group,
        "consensus_definition_used": consensus_def,
        "n_active_reactions": len(unblocked2),
        "obj_frac": args.obj_frac,
        "biomass_floor": floor,
        "fba_optimum_used": fresh_fba,
        "n_samples": args.n_samples,
        "thinning": args.thinning,
        "processes": args.processes,
        "seed": args.seed,
        "runtime_seconds": round(runtime_s, 2),
        "step6_input_check": input_check,
        "validation": validation,
    }
    validation_path = os.path.join(OUT_DIR, f"{cohort}__{group}__validation.json")
    with open(validation_path, "w") as f:
        json.dump(validation_payload, f, indent=2)

    return validation_payload, samples_path, summary_path, validation_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", type=str, default=None,
                         help="Restrict to a single 'cohort:group' combo, e.g. "
                              "GSE89632:healthy -- for testing before a full run.")
    parser.add_argument("--n-samples", type=int, default=500,
                         help="Total retained samples per group (split into 2 "
                              "sequential batches for the convergence check). Default 500.")
    parser.add_argument("--thinning", type=int, default=100,
                         help="OptGP thinning factor (cobra default 100).")
    parser.add_argument("--processes", type=int, default=4,
                         help="Parallel processes for OptGP (machine has 12 logical cores).")
    parser.add_argument("--seed", type=int, default=42,
                         help="Random seed, for reproducibility.")
    parser.add_argument("--obj-frac", type=float, default=0.9,
                         help="Biomass floor as a fraction of each group's FBA optimum "
                              "(same 0.9 convention used elsewhere in this pipeline).")
    parser.add_argument("--ks-alpha", type=float, default=0.01,
                         help="Significance threshold for the per-reaction split-chain "
                              "KS convergence diagnostic.")
    args = parser.parse_args()

    step06 = _import_step06()

    print("Loading Human-GEM base model...")
    base_model = step06.load_base_model()
    growth_support = step06.load_growth_support()
    print(f"  {len(base_model.reactions)} reactions. Growth-support set: {len(growth_support)} reactions.")
    print(f"  Params: n_samples={args.n_samples} thinning={args.thinning} "
          f"processes={args.processes} seed={args.seed} obj_frac={args.obj_frac}")

    consensus_manifest = pd.read_csv(os.path.join(CONSENSUS_DIR, "consensus_manifest.csv"))
    combos = consensus_manifest[["cohort", "group"]].drop_duplicates().values.tolist()

    if args.only:
        cohort_f, group_f = args.only.split(":")
        combos = [(c, g) for c, g in combos if c == cohort_f and g == group_f]
        if not combos:
            print(f"No matching combo for --only {args.only}")
            sys.exit(1)

    manifest_rows = []
    for cohort, group in combos:
        print(f"\n{cohort}/{group}:")
        payload, samples_path, summary_path, validation_path = process_group(
            step06, cohort, group, base_model, growth_support, args
        )
        v = payload["validation"]
        print(f"  fba_reproduced={payload['step6_input_check']['fba_reproduced']} "
              f"biomass_floor={payload['biomass_floor']:.4f}")
        print(f"  valid_samples={v['pct_samples_structurally_valid']}% "
              f"reactions_in_fva_range={v['pct_reactions_within_fva_range']}% "
              f"ks_significant={v['n_reactions_ks_significant']}/{v['n_reactions_ks_tested']} "
              f"({v['pct_reactions_ks_significant']}%)")
        print(f"  effect size: ks_D_median={v['ks_d_statistic_median']} ks_D_max={v['ks_d_statistic_max']} "
              f"medium+_shift_reactions={v['n_reactions_with_medium_or_larger_shift']} "
              f"({v['pct_reactions_with_medium_or_larger_shift']}%)")
        print(f"  runtime={payload['runtime_seconds']}s -> {samples_path}")

        manifest_rows.append({
            "cohort": cohort, "group": group,
            "consensus_definition_used": payload["consensus_definition_used"],
            "n_active_reactions": payload["n_active_reactions"],
            "n_samples": payload["n_samples"],
            "thinning": payload["thinning"],
            "seed": payload["seed"],
            "obj_frac": payload["obj_frac"],
            "biomass_floor": payload["biomass_floor"],
            "fba_reproduced": payload["step6_input_check"]["fba_reproduced"],
            "pct_samples_structurally_valid": v["pct_samples_structurally_valid"],
            "pct_reactions_within_fva_range": v["pct_reactions_within_fva_range"],
            "pct_reactions_ks_significant": v["pct_reactions_ks_significant"],
            "ks_d_statistic_median": v["ks_d_statistic_median"],
            "pct_reactions_with_medium_or_larger_shift": v["pct_reactions_with_medium_or_larger_shift"],
            "runtime_seconds": payload["runtime_seconds"],
            "samples_path": samples_path,
            "summary_path": summary_path,
            "validation_path": validation_path,
        })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(OUT_DIR, "flux_sampling_manifest.csv")
    if args.only:
        print(f"\n(--only was used; not overwriting the full {manifest_path})")
    else:
        manifest_df.to_csv(manifest_path, index=False)
        print(f"\nManifest saved to {manifest_path}")


if __name__ == "__main__":
    main()
