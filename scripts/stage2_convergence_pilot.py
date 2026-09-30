"""
GEM1 - Stage 2 convergence pilot: how many perturbation Monte Carlo trials
(N) does perturbation_robustness_score need before it is a stable per-reaction
estimator?

Motivation: the 2026-08 manuscript audit (audit/item3_flux_sampling_perturbation.py,
see audit/details/item3_flux_sampling_perturbation.json) independently re-ran
Step 8's protocol at its production n_iterations=200 and found that, while
aggregate statistics matched well (Pearson r~0.956 across all reactions), the
per-reaction score among non-growth-support reactions was essentially
uncorrelated between the two runs (Spearman ~0.005) -- i.e. N=200 is not
enough for a stable per-reaction estimate, only a stable per-group/aggregate
one. This pilot measures inter-replicate reliability (Pearson r between two
independent replicates at the same N, same group) across a small N ladder
for 3 representative groups, to find what N a corrected Stage 2 rerun should
actually use.

Reuses build_group_model, run_trials, compute_reaction_scores, and wilson_ci
directly from scripts/08_perturbation_testing.py (imported via importlib,
same pattern 08 itself uses for 06_flux_analysis.py / 07_flux_sampling.py,
since the filenames start with a digit) -- none of that logic is
reimplemented here, so pilot results are directly comparable to production
Step 8 behavior. active_set and growth_support are loaded from the stored
Step 6/9 consensus files via step06.load_strict_consensus() /
step06.load_growth_support(), not recomputed from raw Step 5 output.

Does not modify 08_perturbation_testing.py, any other existing script, or
any existing data file (Step 6/7/8/9 outputs untouched). Writes only new
files under data/perturbation_pilot/.

Run from the GEM1 conda env, detached (same Start-Process pattern used for
the Step 5/7/8 reruns -- see logs/step8_rerun.log etc. and audit/README.md's
"How to run (long-running items)" section):

Start-Process -FilePath "<path-to-gem1-conda-env>\python.exe" `
  -ArgumentList "scripts\\stage2_convergence_pilot.py" `
  -WorkingDirectory "<repo-root>" `
  -RedirectStandardOutput "<repo-root>\\logs\\stage2_convergence_pilot.log" `
  -RedirectStandardError "<repo-root>\\logs\\stage2_convergence_pilot_err.log" `
  -WindowStyle Hidden
"""

import os
import json
import time
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
SCRIPTS_DIR = os.path.dirname(__file__)
FLUX_DIR = os.path.join(BASE_DIR, "data", "flux_analysis")
OUT_DIR = os.path.join(BASE_DIR, "data", "perturbation_pilot")
os.makedirs(OUT_DIR, exist_ok=True)

# Matches Step 8's production defaults exactly (--magnitude 0.3, --obj-frac 0.99).
MAGNITUDE = 0.3
OBJ_FRAC = 0.99

N_TIERS = [200, 500, 1000, 2000]
N_EXTRA_TIER = 5000
N_REPLICATES = 2

GROUPS = [
    ("GSE135251", "healthy"),
    ("GSE135251", "steatosis"),
    ("GSE126848", "obese_no_NAFLD"),
]

# Decision thresholds for whether N=5000 is worth running (see main() below
# for the reasoning): meaningful absolute reliability gain from N=200 to
# N=2000, and the most recent doubling (1000->2000) still delivering at
# least half the gain the previous doubling (500->1000) delivered -- i.e.
# not yet clearly plateaued -- and not already saturated near r=1.
MIN_MEANINGFUL_IMPROVEMENT = 0.10
MIN_STILL_RISING_RATIO = 0.5
SATURATION_R = 0.97


def _import_module(name, filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def seed_for(n, replicate_idx):
    # Unique across the entire pilot: N tiers are spaced by >=200*10=2000,
    # far more than the +0/+1 replicate offset, so no collisions are possible
    # across any (N, replicate) pair in N_TIERS + [N_EXTRA_TIER].
    return 900000 + n * 10 + replicate_idx


def reliability(scores_a, scores_b):
    """
    Pearson r between two replicates' perturbation_robustness_score, on
    reactions present in both. Note: run_trials()/compute_reaction_scores()
    only ever score dropout_candidates (active_set - growth_support) -- the
    growth-support reactions the audit found deterministically tied at 1.0
    are never part of this script's output in the first place, so no
    explicit tie-exclusion is needed here (unlike the audit's comparison
    against the original run's full-reaction-set output).
    """
    merged = scores_a.merge(scores_b, on="reaction_id", suffixes=("_a", "_b"))
    merged = merged.dropna(subset=["perturbation_robustness_score_a", "perturbation_robustness_score_b"])
    if len(merged) < 2:
        return float("nan"), len(merged)
    r, _ = pearsonr(merged["perturbation_robustness_score_a"], merged["perturbation_robustness_score_b"])
    return r, len(merged)


def run_one_replicate(step08, model, dropout_candidates, growth_support, n, target_biomass, cohort, group, replicate_idx):
    seed = seed_for(n, replicate_idx)
    rng = np.random.default_rng(seed)

    t0 = time.time()
    trials_df, removed_mask = step08.run_trials(
        model, dropout_candidates, growth_support, n, MAGNITUDE, target_biomass, rng
    )
    scores = step08.compute_reaction_scores(trials_df, dropout_candidates, removed_mask)
    runtime_s = time.time() - t0
    feasibility_rate = float(trials_df["feasible"].mean())

    out_path = os.path.join(OUT_DIR, f"{cohort}__{group}__N{n}__rep{replicate_idx}.csv")
    scores.to_csv(out_path, index=False)

    print(
        f"  {cohort}/{group} N={n} rep={replicate_idx} seed={seed} "
        f"runtime={runtime_s:.1f}s overall_feasibility_rate={feasibility_rate:.4f} -> {out_path}",
        flush=True,
    )
    return scores


def run_n_tier(step08, model, dropout_candidates, growth_support, n, target_biomass, cohort, group):
    reps = [
        run_one_replicate(step08, model, dropout_candidates, growth_support, n, target_biomass, cohort, group, i)
        for i in range(N_REPLICATES)
    ]
    r, n_pairs = reliability(reps[0], reps[1])
    print(f"  {cohort}/{group} N={n} inter-replicate Pearson r (n_reactions={n_pairs}) = {r:.4f}", flush=True)
    return r


def main():
    step06 = _import_module("step06_flux_analysis", "06_flux_analysis.py")
    step07 = _import_module("step07_flux_sampling", "07_flux_sampling.py")
    step08 = _import_module("step08_perturbation_testing", "08_perturbation_testing.py")

    print("Loading Human-GEM base model...", flush=True)
    base_model = step06.load_base_model()
    growth_support = step06.load_growth_support()
    print(f"  {len(base_model.reactions)} reactions. Growth-support set: {len(growth_support)} reactions.", flush=True)
    print(f"  Params: magnitude={MAGNITUDE} obj_frac={OBJ_FRAC} N_tiers={N_TIERS} replicates={N_REPLICATES}", flush=True)

    group_state = {}
    reliability_by_n = {n: {} for n in N_TIERS}

    for cohort, group in GROUPS:
        print(f"\n=== {cohort}/{group} ===", flush=True)
        active_set, consensus_def = step06.load_strict_consensus(cohort, group)
        model, unblocked = step08.build_group_model(step06, step07, base_model.copy(), active_set, growth_support)
        dropout_candidates = sorted(active_set - growth_support)

        with open(os.path.join(FLUX_DIR, f"{cohort}__{group}__flux.json")) as f:
            step6_fba = json.load(f)["fba_objective_value"]
        target_biomass = OBJ_FRAC * step6_fba
        print(
            f"  consensus_definition_used={consensus_def} n_dropout_candidates={len(dropout_candidates)} "
            f"target_biomass={target_biomass:.6f} (0.99 * stored Step 6 fba_objective_value={step6_fba:.6f})",
            flush=True,
        )

        group_state[(cohort, group)] = (model, dropout_candidates, target_biomass)

        for n in N_TIERS:
            r = run_n_tier(step08, model, dropout_candidates, growth_support, n, target_biomass, cohort, group)
            reliability_by_n[n][(cohort, group)] = r

    print("\n=== N=200 -> N=2000 reliability summary ===", flush=True)
    improvements = []
    for cohort, group in GROUPS:
        r200 = reliability_by_n[200][(cohort, group)]
        r2000 = reliability_by_n[2000][(cohort, group)]
        improvements.append(r2000 - r200)
        print(f"  {cohort}/{group}: r(N=200)={r200:.4f}  r(N=2000)={r2000:.4f}  improvement={r2000 - r200:+.4f}", flush=True)

    mean_r200 = np.mean([reliability_by_n[200][g] for g in GROUPS])
    mean_r500 = np.mean([reliability_by_n[500][g] for g in GROUPS])
    mean_r1000 = np.mean([reliability_by_n[1000][g] for g in GROUPS])
    mean_r2000 = np.mean([reliability_by_n[2000][g] for g in GROUPS])

    improvement_200_to_2000 = mean_r2000 - mean_r200
    gain_500_to_1000 = mean_r1000 - mean_r500
    gain_1000_to_2000 = mean_r2000 - mean_r1000
    still_rising = gain_1000_to_2000 >= MIN_STILL_RISING_RATIO * gain_500_to_1000 if gain_500_to_1000 > 0 else gain_1000_to_2000 > 0
    meaningful = improvement_200_to_2000 >= MIN_MEANINGFUL_IMPROVEMENT
    saturated = mean_r2000 >= SATURATION_R

    run_5000 = meaningful and still_rising and not saturated

    print(
        f"\nDecision inputs: mean_r200={mean_r200:.4f} mean_r500={mean_r500:.4f} "
        f"mean_r1000={mean_r1000:.4f} mean_r2000={mean_r2000:.4f}",
        flush=True,
    )
    print(
        f"  improvement(200->2000)={improvement_200_to_2000:+.4f} (threshold >= {MIN_MEANINGFUL_IMPROVEMENT}) -> meaningful={meaningful}",
        flush=True,
    )
    print(
        f"  gain(500->1000)={gain_500_to_1000:+.4f} gain(1000->2000)={gain_1000_to_2000:+.4f} "
        f"(need latter >= {MIN_STILL_RISING_RATIO} x former) -> still_rising={still_rising}",
        flush=True,
    )
    print(f"  saturated (mean_r2000 >= {SATURATION_R})={saturated}", flush=True)
    print(f"  => {'RUNNING N=5000' if run_5000 else 'SKIPPING N=5000'}", flush=True)

    if run_5000:
        print(f"\n=== N={N_EXTRA_TIER} (reliability still rising meaningfully at N=2000) ===", flush=True)
        for cohort, group in GROUPS:
            model, dropout_candidates, target_biomass = group_state[(cohort, group)]
            r = run_n_tier(step08, model, dropout_candidates, growth_support, N_EXTRA_TIER, target_biomass, cohort, group)
            print(f"  {cohort}/{group}: r(N={N_EXTRA_TIER})={r:.4f}", flush=True)

    print("\nStage 2 convergence pilot complete.", flush=True)


if __name__ == "__main__":
    main()
