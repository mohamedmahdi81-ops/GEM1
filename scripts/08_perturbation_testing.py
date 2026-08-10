"""
GEM1 - Step 8: Monte Carlo perturbation testing (structural robustness).

Purpose: test how robust each group's STRUCTURAL consensus network (Step 9's
strict_consensus_reactions) is to reconstruction/algorithm-disagreement
noise -- a distinct axis from Step 7 (which characterizes flux variability
WITHIN a fixed network). Does not touch Step 6/7/9/10 in any way; both
remain exactly as already computed and locked.

Methodology (design confirmed with Mohammed 2026-07-19, before implementing):
each Monte Carlo trial independently removes every strict_consensus_reactions
reaction with probability `magnitude` (the growth-support set is NEVER
eligible for removal -- it's a deliberately forced hard-required scaffold,
not part of the "consensus signal" being stress-tested; dropping it would
trivially and uninformatively break every trial), then re-solves FBA on the
perturbed model and records feasibility (biomass >= obj_frac * that group's
unperturbed FBA optimum). obj_frac defaults to 0.99, NOT this pipeline's
usual 0.9 convention -- found via smoke-testing that 0.9 makes this step
tautologically always-feasible: the growth-support set was itself derived
(Step 5) to guarantee exactly 0.9 of max biomass on its own, so with 0.9 as
the bar, structural_robustness_index is 1.0 by construction regardless of
what else is removed (confirmed empirically: still 1.0 even removing 50% of
the ~3,600 non-growth-support consensus reactions per trial). 0.99 requires
near-full network contribution beyond the protected backbone, so the test is
actually sensitive to what's removed. See "Known open items" in CLAUDE.md
for the full diagnosis.

Per-reaction robustness score: since each reaction is independently present
or absent across the N trials, this script estimates each reaction's
marginal association with feasibility:
    perturbation_robustness_score(R) =
        P(feasible | R present) - P(feasible | R absent)
A reaction structurally required for growth will show a large positive
value (feasible far more often when present). A redundant/replaceable
reaction will show a value near zero. This is an approximate, correlational
marginal-effect estimate (multiple reactions are perturbed simultaneously
per trial, so this is not a controlled, isolated causal effect for any
single reaction) -- documented as such, not oversold as exact.

Group-level structural_robustness_index = fraction of all N trials that
remained feasible at the given magnitude -- the headline cross-group/
cross-cohort comparable metric, reported with a Wilson-score binomial
confidence interval (quantifying the Monte Carlo estimate's own
uncertainty, not just a bare point estimate).

Model scaffold reused from Step 6/7 (via importlib, module names start with
digits): same active-reaction-set definitions, same base model. For speed
(potentially thousands of FBA solves across all groups), the REDUCED
(reaction-removed, not bounds-zeroed) model is built ONCE per group -- same
performance fix already validated in Step 7 -- then each trial temporarily
bounds-zeroes the randomly-selected removed subset inside a `with model:`
context (auto-reverting cobra context manager), avoiding a full model
rebuild per trial.

Reproducibility: every run is seeded (`--seed`, default 42) via a single
numpy Generator advanced sequentially across trials (not re-seeded per
trial), and every parameter (`--magnitude`, `--n-iterations`, `--obj-frac`)
is a CLI flag with a documented default.

Validation performed on every group before results are trusted:
  a. Baseline (zero reactions removed) FBA reproduces Step 6's stored
     fba_objective_value -- confirms Step 6/7's scaffold is still a valid,
     reproducible input before perturbing it.
  b. Growth-support reactions are never present in any trial's removed set
     (a hard structural invariant, checked directly against every trial's
     log, not assumed from the sampling code alone).
  c. Feasibility rate is monotonically non-increasing (within Monte Carlo
     noise) as a function of how many reactions were removed in a trial --
     binned by n_removed decile, reported for inspection rather than
     silently assumed.
  d. Split-half stability: structural_robustness_index computed from the
     first half of trials vs the second half, with each half's own Wilson
     CI -- large, non-overlapping CIs would indicate n_iterations is too
     small for a stable estimate.

Outputs per group in data/perturbation_testing/:
  {cohort}__{group}__reactions.csv -- per-reaction robustness scores.
  {cohort}__{group}__trials.csv.gz -- complete per-trial log (trial index,
    n_removed, feasible, objective_value, removed reaction IDs) -- full
    audit trail, not just aggregates.
  {cohort}__{group}__validation.json -- all four validation checks above.
perturbation_manifest.csv -- one row per group: parameters, headline
  structural_robustness_index + CI, runtime, validation rollup.

Run from gem1-main (GEM1 conda env) -- standard cobra (0.31.1), GLPK solver
(pure LP, no MILP in this step -- same reasoning as Step 7's solver fix).
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

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.dirname(__file__)
CONSENSUS_DIR = os.path.join(BASE_DIR, "data", "consensus_scores")
FLUX_DIR = os.path.join(BASE_DIR, "data", "flux_analysis")
OUT_DIR = os.path.join(BASE_DIR, "data", "perturbation_testing")
os.makedirs(OUT_DIR, exist_ok=True)

BIOMASS_ID = "MAR13082"


def _import_module(name, filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_group_model(step06, step07, base_model, active_set, growth_support):
    """
    One reduced (reaction-removed) model per group, built once and reused
    across all trials -- same performance-critical technique validated in
    Step 7 (removing a reaction pinned to (0,0) cannot change any surviving
    reaction's feasible range, so this is equivalent to bounds-zeroing, just
    far cheaper for the many solves this step needs).
    """
    model = base_model.copy()
    unblocked = active_set | growth_support
    to_remove = [r for r in model.reactions if r.id not in unblocked]
    model.remove_reactions(to_remove, remove_orphans=True)
    model.solver = "glpk"
    model.objective = BIOMASS_ID
    return model, unblocked


def run_trials(model, dropout_candidates, growth_support, n_iterations, magnitude, target_biomass, rng):
    """
    dropout_candidates: sorted list of strict_consensus reaction IDs eligible
    for random removal (growth_support reactions are structurally excluded
    from this list by the caller, never eligible).

    Returns (trials_df, removed_mask) where removed_mask is a boolean
    (n_iterations x n_candidates) array -- used directly by
    compute_reaction_scores() for vectorized per-reaction scoring, instead
    of re-parsing the human-readable `removed_reactions` string column
    (kept in trials_df purely for the audit-trail log). The original
    string-regex-based scoring approach took ~110s for just 20 trials on
    one group -- almost entirely regex overhead, not the FBA solves
    themselves -- found via a smoke test before committing to the full run.
    """
    candidates = np.array(dropout_candidates)
    n_candidates = len(candidates)
    removed_mask = np.zeros((n_iterations, n_candidates), dtype=bool)
    trials = []
    for trial_idx in range(n_iterations):
        mask = rng.random(n_candidates) < magnitude
        removed_mask[trial_idx] = mask
        removed = candidates[mask]
        assert not (set(removed) & growth_support), (
            f"Trial {trial_idx}: growth-support reaction(s) ended up in the removed set -- "
            "structural invariant violated, this must never happen."
        )
        with model:
            for rid in removed:
                model.reactions.get_by_id(rid).bounds = (0, 0)
            sol = model.optimize()
            obj = sol.objective_value if sol.status == "optimal" else None
        feasible = obj is not None and obj >= target_biomass - 1e-6
        trials.append({
            "trial": trial_idx,
            "n_removed": int(mask.sum()),
            "feasible": feasible,
            "objective_value": obj,
            "removed_reactions": ";".join(removed),
        })
    return pd.DataFrame(trials), removed_mask


def wilson_ci(n_success, n_total, z=1.96):
    if n_total == 0:
        return (0.0, 0.0, 0.0)
    p = n_success / n_total
    denom = 1 + z ** 2 / n_total
    center = (p + z ** 2 / (2 * n_total)) / denom
    half = (z * np.sqrt(p * (1 - p) / n_total + z ** 2 / (4 * n_total ** 2))) / denom
    return (p, max(0.0, center - half), min(1.0, center + half))


def compute_reaction_scores(trials_df, dropout_candidates, removed_mask):
    """
    Vectorized via the boolean removed_mask (n_iterations x n_candidates)
    returned by run_trials(), rather than re-parsing the removed_reactions
    string per reaction (see run_trials()'s docstring for why that was
    replaced -- it was the dominant cost of this whole step).
    """
    presence = ~removed_mask  # True where reaction was present (not removed) that trial
    feasible = trials_df["feasible"].values.astype(float)
    n_iterations = len(feasible)

    n_present = presence.sum(axis=0)  # per-candidate count across trials
    n_absent = n_iterations - n_present
    feasible_and_present = presence.T @ feasible  # per-candidate sum of feasible trials where present
    total_feasible = feasible.sum()
    feasible_and_absent = total_feasible - feasible_and_present

    with np.errstate(invalid="ignore", divide="ignore"):
        rate_present = np.where(n_present > 0, feasible_and_present / np.maximum(n_present, 1), np.nan)
        rate_absent = np.where(n_absent > 0, feasible_and_absent / np.maximum(n_absent, 1), np.nan)
    score = rate_present - rate_absent

    return pd.DataFrame({
        "reaction_id": dropout_candidates,
        "n_trials_present": n_present.astype(int),
        "n_trials_absent": n_absent.astype(int),
        "feasibility_rate_when_present": rate_present,
        "feasibility_rate_when_absent": rate_absent,
        "perturbation_robustness_score": score,
    })


def validate_group(step06, cohort, group, model, trials_df, growth_support, target_biomass):
    flux_json_path = os.path.join(FLUX_DIR, f"{cohort}__{group}__flux.json")
    with open(flux_json_path) as f:
        step6_result = json.load(f)

    with model:
        model.objective = BIOMASS_ID
        baseline_sol = model.optimize()
        baseline_obj = baseline_sol.objective_value if baseline_sol.status == "optimal" else None
    baseline_matches = baseline_obj is not None and abs(baseline_obj - step6_result["fba_objective_value"]) < 1e-4

    all_removed = set()
    for s in trials_df["removed_reactions"]:
        if s:
            all_removed.update(s.split(";"))
    growth_support_leak = len(all_removed & growth_support)

    trials_df = trials_df.copy()
    trials_df["n_removed_decile"] = pd.qcut(trials_df["n_removed"], q=min(10, trials_df["n_removed"].nunique()), duplicates="drop")
    monotonicity = trials_df.groupby("n_removed_decile", observed=True)["feasible"].mean().to_dict()
    monotonicity = {str(k): round(float(v), 4) for k, v in monotonicity.items()}

    n = len(trials_df)
    half = n // 2
    first, second = trials_df.iloc[:half], trials_df.iloc[half:]
    p1, lo1, hi1 = wilson_ci(int(first["feasible"].sum()), len(first))
    p2, lo2, hi2 = wilson_ci(int(second["feasible"].sum()), len(second))

    return {
        "step6_fba_objective_value": step6_result["fba_objective_value"],
        "baseline_fba_objective_value": baseline_obj,
        "baseline_reproduces_step6": baseline_matches,
        "target_biomass": target_biomass,
        "growth_support_reactions_ever_removed": growth_support_leak,
        "feasibility_rate_by_n_removed_decile": monotonicity,
        "split_half_first": {"p": round(p1, 4), "ci_lo": round(lo1, 4), "ci_hi": round(hi1, 4), "n": len(first)},
        "split_half_second": {"p": round(p2, 4), "ci_lo": round(lo2, 4), "ci_hi": round(hi2, 4), "n": len(second)},
        "split_half_cis_overlap": not (hi1 < lo2 or hi2 < lo1),
    }


def process_group(step06, step07, cohort, group, base_model, growth_support, args, rng):
    t0 = time.time()
    active_set, consensus_def = step06.load_strict_consensus(cohort, group)
    model, unblocked = build_group_model(step06, step07, base_model.copy(), active_set, growth_support)

    with open(os.path.join(FLUX_DIR, f"{cohort}__{group}__flux.json")) as f:
        step6_fba = json.load(f)["fba_objective_value"]
    target_biomass = args.obj_frac * step6_fba

    dropout_candidates = sorted(active_set - growth_support)
    trials_df, removed_mask = run_trials(model, dropout_candidates, growth_support, args.n_iterations, args.magnitude, target_biomass, rng)

    reaction_scores = compute_reaction_scores(trials_df, dropout_candidates, removed_mask)
    validation = validate_group(step06, cohort, group, model, trials_df, growth_support, target_biomass)

    n_feasible = int(trials_df["feasible"].sum())
    p, ci_lo, ci_hi = wilson_ci(n_feasible, len(trials_df))
    runtime_s = time.time() - t0

    reactions_path = os.path.join(OUT_DIR, f"{cohort}__{group}__reactions.csv")
    reaction_scores.to_csv(reactions_path, index=False)
    trials_path = os.path.join(OUT_DIR, f"{cohort}__{group}__trials.csv.gz")
    trials_df.drop(columns=["n_removed_decile"], errors="ignore").to_csv(trials_path, index=False, compression="gzip")

    payload = {
        "cohort": cohort, "group": group,
        "consensus_definition_used": consensus_def,
        "n_strict_consensus_reactions": len(active_set - growth_support),
        "n_growth_support_reactions": len(growth_support),  # always the full 519, unioned in regardless of overlap with active_set
        "n_dropout_candidates": len(dropout_candidates),
        "magnitude": args.magnitude,
        "n_iterations": args.n_iterations,
        "seed": args.seed,
        "obj_frac": args.obj_frac,
        "target_biomass": target_biomass,
        "structural_robustness_index": p,
        "structural_robustness_ci_lo": ci_lo,
        "structural_robustness_ci_hi": ci_hi,
        "runtime_seconds": round(runtime_s, 2),
        "validation": validation,
    }
    validation_path = os.path.join(OUT_DIR, f"{cohort}__{group}__validation.json")
    with open(validation_path, "w") as f:
        json.dump(payload, f, indent=2)

    return payload, reactions_path, trials_path, validation_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", type=str, default=None,
                         help="Restrict to a single 'cohort:group' combo, e.g. "
                              "GSE89632:healthy -- for testing before a full run.")
    parser.add_argument("--magnitude", type=float, default=0.3,
                         help="Per-reaction independent removal probability per trial. Default 0.3 (30%%), "
                              "calibrated empirically 2026-07-19 (scripts/_08 smoke tests): 0.05-0.15 showed "
                              "no discrimination at obj_frac=0.99 (100%% robust), 0.5 was heavily saturated "
                              "toward infeasible (~40%% robust); 0.3 gives a well-centered, non-saturated "
                              "structural_robustness_index (~0.6-0.8 range observed) informative enough to "
                              "differentiate groups.")
    parser.add_argument("--n-iterations", type=int, default=200,
                         help="Number of Monte Carlo trials per group. Default 200.")
    parser.add_argument("--seed", type=int, default=42,
                         help="Random seed (single Generator advanced sequentially across all trials).")
    parser.add_argument("--obj-frac", type=float, default=0.99,
                         help="Feasibility threshold as a fraction of each group's unperturbed FBA "
                              "optimum. Deliberately NOT the pipeline's usual 0.9 convention -- "
                              "found via smoke-testing (2026-07-19) that 0.9 makes this step "
                              "tautologically always-feasible, since the permanently-protected "
                              "519-reaction growth-support set was itself derived to guarantee "
                              "exactly 0.9 of max biomass on its own, independent of anything else "
                              "removed. 0.99 requires near-full network contribution, so the test is "
                              "actually sensitive to which consensus-only reactions are removed.")
    args = parser.parse_args()

    step06 = _import_module("step06_flux_analysis", "06_flux_analysis.py")
    step07 = _import_module("step07_flux_sampling", "07_flux_sampling.py")

    print("Loading Human-GEM base model...")
    base_model = step06.load_base_model()
    growth_support = step06.load_growth_support()
    print(f"  {len(base_model.reactions)} reactions. Growth-support set: {len(growth_support)} reactions "
          f"(never eligible for removal).")
    print(f"  Params: magnitude={args.magnitude} n_iterations={args.n_iterations} "
          f"seed={args.seed} obj_frac={args.obj_frac}")

    consensus_manifest = pd.read_csv(os.path.join(CONSENSUS_DIR, "consensus_manifest.csv"))
    combos = consensus_manifest[["cohort", "group"]].drop_duplicates().values.tolist()

    if args.only:
        cohort_f, group_f = args.only.split(":")
        combos = [(c, g) for c, g in combos if c == cohort_f and g == group_f]
        if not combos:
            print(f"No matching combo for --only {args.only}")
            sys.exit(1)

    rng = np.random.default_rng(args.seed)

    manifest_rows = []
    for cohort, group in combos:
        print(f"\n{cohort}/{group}:")
        payload, reactions_path, trials_path, validation_path = process_group(
            step06, step07, cohort, group, base_model, growth_support, args, rng
        )
        v = payload["validation"]
        print(f"  baseline_reproduces_step6={v['baseline_reproduces_step6']} "
              f"target_biomass={payload['target_biomass']:.4f}")
        print(f"  structural_robustness_index={payload['structural_robustness_index']:.4f} "
              f"CI=[{payload['structural_robustness_ci_lo']:.4f}, {payload['structural_robustness_ci_hi']:.4f}] "
              f"growth_support_leak={v['growth_support_reactions_ever_removed']}")
        print(f"  split_half: first={v['split_half_first']['p']:.4f} second={v['split_half_second']['p']:.4f} "
              f"CIs_overlap={v['split_half_cis_overlap']}")
        print(f"  runtime={payload['runtime_seconds']}s -> {reactions_path}")

        manifest_rows.append({
            "cohort": cohort, "group": group,
            "consensus_definition_used": payload["consensus_definition_used"],
            "n_strict_consensus_reactions": payload["n_strict_consensus_reactions"],
            "n_dropout_candidates": payload["n_dropout_candidates"],
            "magnitude": payload["magnitude"],
            "n_iterations": payload["n_iterations"],
            "seed": payload["seed"],
            "structural_robustness_index": payload["structural_robustness_index"],
            "structural_robustness_ci_lo": payload["structural_robustness_ci_lo"],
            "structural_robustness_ci_hi": payload["structural_robustness_ci_hi"],
            "baseline_reproduces_step6": v["baseline_reproduces_step6"],
            "growth_support_reactions_ever_removed": v["growth_support_reactions_ever_removed"],
            "split_half_cis_overlap": v["split_half_cis_overlap"],
            "runtime_seconds": payload["runtime_seconds"],
            "reactions_path": reactions_path,
            "trials_path": trials_path,
            "validation_path": validation_path,
        })

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = os.path.join(OUT_DIR, "perturbation_manifest.csv")
    if args.only:
        print(f"\n(--only was used; not overwriting the full {manifest_path})")
    else:
        manifest_df.to_csv(manifest_path, index=False)
        print(f"\nManifest saved to {manifest_path}")


if __name__ == "__main__":
    main()
