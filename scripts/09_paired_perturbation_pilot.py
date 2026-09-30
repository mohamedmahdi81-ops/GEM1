"""
GEM1 - Step 9 (pilot): Paired-counterfactual perturbation redesign.

STATUS: design/spec only, not yet run. Requires the actual GEM1 environment
(cobra, Human-GEM model, 06_flux_analysis.py / 07_flux_sampling.py /
08_perturbation_testing.py present in scripts/, consensus_scores + flux_analysis
data on disk) -- none of which exist in the diagnostics sandbox this was written
in. This file is a precise, ready-to-run specification for Claude Code / Mohammed
to execute in the actual GEM1 environment, per the paired-counterfactual redesign
agreed 2026-09-18 (see claude/GEM1-phase1-paired-counterfactual-design.md), revised
the same day per review comments (reaction-subset diversity, preregistered primary
endpoint, relative-objective statistic, soft violation logging -- see that doc's
"Revisions after review" section for the full rationale behind each change below).

Why this exists (one paragraph, see the project doc for full reasoning):
08_perturbation_testing.py's per-reaction score, rate_present - rate_absent,
compares DIFFERENT random trials (different ~1,000-reaction backgrounds) to
each other. Phase-1 regression screening (ridge/logistic, then a stricter
disjoint-half test) showed that cross-trial background variation swamps any
individual reaction's marginal signal: coefficient Spearman correlation
between independent disjoint halves of the same N=200 data was ~0 in all
three groups tested. This script instead evaluates THE SAME background
perturbation B twice -- once with reaction R present, once with R additionally
removed -- so the only thing that differs between the paired arms is R itself.
This is the paired-difference / common-random-numbers design (equivalently,
the standard unbiased estimator of a coordinate's "influence" in the analysis
of Boolean functions, or a Morris-type elementary-effects design in global
sensitivity analysis) applied to this problem. It keeps GEM1's intended
concept -- R's contribution to feasibility under an uncertain, heavily
perturbed network background -- rather than collapsing to a single-reaction
knockout (which would fix background = the intact network, a different and
much weaker biological question).

Compute-efficiency design (the reason this isn't just "run 08 twice per
reaction"): background draws are generated ONCE per group and reused across
EVERY reaction in the pilot subset S. For a given background B, the
"baseline" (R present) arm is identical for every R that happens to not be
in B already (~(1 - magnitude) of S on average) -- so only ONE solve is
needed per background for all of those reactions' baseline arms, not one
per (reaction, background) pair. Only the "R additionally removed" arm needs
a new solve per (reaction, background) pair, and it is done via a NESTED
`with model:` context on top of the already-applied background bounds --
i.e. it is a single incremental bound change (R's bounds -> (0,0)) from an
already-solved basis, which GLPK/cobra should be able to warm-start from
(measure this empirically first -- see `--timing-only`, do not assume it).

Total solves for M backgrounds x |S| reactions at removal probability m:
    M * (1 + |S| * (1 - m))
Example: M=120, |S|=44, m=0.3 -> 120 * (1 + 44*0.7) = 120 * 31.8 = 3,816 solves.
Each reaction is eligible in ~M*(1-m) ~= 84 paired backgrounds per seed at
M=120 -- print this explicitly at startup so the expected precision is visible
before the run, not discovered afterward.

M=120 IS AN INITIAL PILOT DEPTH, NOT A FINAL STATISTICAL DEPTH. Do not
auto-escalate it. The intended sequence is: --timing-only -> M=120 at two
independent seeds -> inspect signal (scatter, Spearman/Pearson, top-K overlap,
per-reaction SE, whether the same reactions repeatedly show nonzero effects)
-> increase M only if that inspection justifies it. ~84 paired observations
per reaction is probably enough to see whether PR_objective (continuous) has
a real signal; it may well be too few for PR_binary if actual feasibility
flips are rare for most reactions in the subset (a reaction with 2-3 flips
out of 84 backgrounds has a very imprecise flip probability regardless of how
good the pairing is) -- this is expected going in, not a failure of the pilot,
and is exactly why PR_objective is preregistered as primary (see below).

Structural invariant (worth asserting, not assuming): removing a reaction
can only shrink the LP's feasible flux space, so for any (B, R) pair,
objective_value(B + R removed) <= objective_value(B) always, and a network
feasible at baseline can only stay feasible or become infeasible when R is
additionally removed -- never the reverse. VIOLATIONS ARE LOGGED, NOT
IMMEDIATELY FATAL: a bare `assert` that kills the run on the first violation
would discard everything computed so far over what might be pure numerical
noise (LP solver tolerance is typically ~1e-9 to 1e-6). Instead, every
violation is recorded with its magnitude, the run continues, and a bucketed
violation summary (numerical-tolerance-scale vs. clearly-a-bug-scale) is
printed and saved at the end; the run only hard-fails (nonzero exit, after
writing whatever was computed) if any single violation exceeds
`--violation-hard-fail-threshold` (default 1e-3, well above ordinary solver
tolerance -- a violation that large means investigate a real bug, e.g. a
growth-support leak, not "average away small numerical noise").

PREREGISTERED interpretation (decided now, not chosen post hoc after seeing
which variant looks best):
  PRIMARY:    PR_objective (and PR_relative_objective, see below) -- the
              continuous objective drop is the highest-power variant: it
              captures partial biomass degradation that never crosses the
              0.99-of-optimum feasibility threshold, which PR_binary and
              PR_hybrid both discard entirely (e.g. baseline optimum 115,
              with-R-removed optimum 113 -> PR_binary sees delta=0, i.e. two
              reactions with very different true effects can look identical
              to it, while PR_objective correctly sees a difference of 2).
  SECONDARY:  PR_binary -- direct paired analogue of the current production
              score, useful for comparability, but expected to be noisy at
              M=120 if flips are rare (see above).
  DIAGNOSTIC: PR_hybrid = P(flip) * E[delta_obj | flip], and its two-part
              decomposition (p_flip, mean_obj_drop_given_flip) -- useful for
              interpreting *why* a reaction scores as it does (rare-but-severe
              vs. reliably-mild), not a replacement for PR_objective, whose
              unconditional mean already folds frequency and magnitude
              together without an ad hoc combination rule.
Also computed (zero extra solver calls, purely derived from data already
collected): PR_relative_objective = mean over eligible backgrounds of
delta_obj / max(baseline_objective, epsilon) -- relative biomass loss.
Doesn't matter much for a single-group pilot with one baseline capacity
scale, but becomes important once/if this expands across cohorts/groups
with different baseline biomass capacities, where raw PR_objective isn't
directly comparable. CAUTION: when baseline_objective is itself near zero
(the background alone is already nearly infeasible), this ratio can blow up
-- report both mean and median for this statistic, and inspect the
denominator distribution before trusting it, rather than assuming it behaves
like a bounded percentage.

Reaction subset (see proposed_pilot_reaction_subset_v2.csv, and the design
doc's "Revisions after review" section for full construction detail):
NOT selected solely by the old (unstable) production score anymore -- that
would make the new method's validation partly dependent on the estimator
being replaced. Built from four independent sources instead: (A) old-score
extremes, 6 reactions, for comparison/continuity only; (B) old-score
near-zero, 6 reactions, likewise for comparison; (C) score-blind random draw
from the full candidate pool, 28 reactions, selected with NO reference to the
old score at all -- the component that actually tests whether the paired
statistic discovers reaction-specific effects rather than reproducing the old
one; (D) subsystem diversity -- NOT YET INCORPORATED, the subsystem mapping
file (data/calibration/all_groups_calibrated_with_metadata.csv) isn't
available in the diagnostics session this was designed in; recommend checking
component C's actual subsystem spread in the real environment (where that
file exists) before freezing the list, and topping up manually if it's
lopsided; (E) biologically interpretable / calibration reactions, pulled from
GEM1-step11-biomarker-mapping.md's 5 locked calibration axes (serine/glycine,
urea cycle, BCAA, lipogenesis, choline/PC) -- only 5 of the 14 candidate
biomarker reactions turned out to be in GSE126848/obese_no_NAFLD's
dropout-eligible universe at all (the other 9 are apparently growth-support-
protected or outside this group's active_set -- worth independently
confirming which, and worth noting in the manuscript's limitations either
way: several of GEM1's own literature-locked calibration reactions currently
cannot even be perturbation-tested in this group under the existing
protocol). Total 44 reactions, close to the originally proposed 45-reaction
scale.

Validation before trusting pilot results (mirrors 08's own validation
philosophy -- verify against source data, don't assume):
  a. Per-pair monotonicity, logged not asserted-fatal (see above).
  b. Growth-support reactions never touched (inherited invariant from 08,
     re-checked here since this script separately manipulates bounds).
  c. Two independent seeds (independent background draws, same reaction
     subset S, same M) -- report Pearson/Spearman of each variant's
     per-reaction score between the two seeds, and top-K overlap at
     K=10/20 (|S|=44 is too small for the K=50/100 used in the earlier
     diagnostics) -- using the SAME reproducibility metrics as
     GEM1-perturbation-reproducibility-diagnostics.md and
     GEM1-phase1-regression-screening.md for direct comparability.
     NO PRE-SET PASS/FAIL THRESHOLD (e.g. "Spearman must exceed 0.8") --
     with only 44 reactions that would be arbitrarily rigid. Instead compare
     directly against what's already on file: the old method produced
     rho~=0 (regression-adjusted) / rho~=0.005 (raw score) for independent
     ranking. The paired method needs to show a large, consistent
     improvement across MULTIPLE lines of evidence together -- score
     scatter, Spearman, Pearson, top-10/20 overlap, per-reaction SE, and
     whether the same reactions repeatedly show nonzero effects -- not just
     one number clearing a bar. Spearman 0.7-0.9 with coherent top-K overlap
     would be compelling; 0.1-0.2 is not a win just because it beats 0.005.
  d. Reports observed wall-clock time for cold vs. warm-started solves,
     so the M x |S| budget for a full run (if the pilot succeeds) can be
     estimated from real numbers, not the guess above.

STAGED VALIDATION PLAN (do not skip stages): (1) --timing-only on
GSE126848/obese_no_NAFLD -- the group with the only real predictive signal
in the Phase-1 regression screen, an intentionally FAVORABLE test case for
"does the concept work at all". (2) If (1)'s timing is sane, run M=120 at
two seeds on the same favorable group. (3) If (2) shows a large, coherent
improvement over the rho~=0 baseline, repeat the SAME reaction subset and
M on one of the difficult groups (healthy or steatosis) with independent
seeds again, to test whether the improvement generalizes rather than being
specific to the one group that was already the best case. (4) Only after
both a favorable AND a difficult group show real improvement should
full-scale (~3,500+ reaction) computation even be considered -- a single
group's success is not sufficient justification for that jump.

Run from gem1-main (GEM1 conda env), same as 08_perturbation_testing.py.
"""

import os
import sys
import json
import time
import argparse
import importlib.util
import numpy as np
import pandas as pd

SCRIPTS_DIR = os.path.dirname(__file__)
BASE_DIR = os.path.dirname(SCRIPTS_DIR)
OUT_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_pilot")
os.makedirs(OUT_DIR, exist_ok=True)

VIOLATION_NUMERICAL_TOL = 1e-6  # below this, not even logged -- ordinary LP solver tolerance


def _import_module(name, filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_paired_backgrounds(model, dropout_candidates, priority_reactions, target_biomass,
                            n_backgrounds, magnitude, rng, rel_epsilon=1e-6, timing_only=False):
    """
    dropout_candidates: full sorted list of a group's dropout-eligible reactions
        (identical to 08_perturbation_testing.py's `dropout_candidates` -- same
        pool the background is drawn from, so PR estimates stay comparable to
        the production score's universe).
    priority_reactions: the small pilot subset S (list of reaction IDs, subset
        of dropout_candidates) to compute paired scores for.

    Returns (background_log_df, pair_log_df, violations_df, timing_dict).
    Violations are logged, not fatal -- see module docstring. Caller decides
    whether to hard-fail based on violations_df's worst magnitude.
    """
    candidates = np.array(dropout_candidates)
    n_candidates = len(candidates)
    priority_set = set(priority_reactions)
    assert priority_set.issubset(set(dropout_candidates)), \
        "priority_reactions must be a subset of dropout_candidates"

    background_rows = []
    pair_rows = []
    violation_rows = []
    cold_times, warm_times = [], []

    for b_idx in range(n_backgrounds):
        mask = rng.random(n_candidates) < magnitude
        removed = candidates[mask]
        removed_set = set(removed)

        t0 = time.time()
        with model:
            for rid in removed:
                model.reactions.get_by_id(rid).bounds = (0, 0)
            sol = model.optimize()
            base_obj = sol.objective_value if sol.status == "optimal" else 0.0
            base_feasible = sol.status == "optimal" and base_obj >= target_biomass - 1e-6
            cold_times.append(time.time() - t0)

            background_rows.append({
                "background": b_idx, "n_removed": int(mask.sum()),
                "baseline_objective": base_obj, "baseline_feasible": bool(base_feasible),
            })

            eligible = [r for r in priority_reactions if r not in removed_set]
            for rid in eligible:
                t1 = time.time()
                with model:
                    model.reactions.get_by_id(rid).bounds = (0, 0)
                    sol_r = model.optimize()
                    obj_r = sol_r.objective_value if sol_r.status == "optimal" else 0.0
                    feasible_r = sol_r.status == "optimal" and obj_r >= target_biomass - 1e-6
                warm_times.append(time.time() - t1)

                delta_obj = base_obj - obj_r
                delta_bin = int(base_feasible) - int(feasible_r)
                delta_rel = delta_obj / max(base_obj, rel_epsilon)

                # structural invariant: removing a reaction cannot help feasibility.
                # LOGGED, not asserted-fatal -- see module docstring.
                if delta_obj < -VIOLATION_NUMERICAL_TOL or delta_bin == -1:
                    violation_rows.append({
                        "background": b_idx, "reaction_id": rid,
                        "baseline_objective": base_obj, "with_r_removed_objective": obj_r,
                        "delta_objective": delta_obj, "delta_binary": delta_bin,
                    })

                pair_rows.append({
                    "background": b_idx, "reaction_id": rid,
                    "baseline_objective": base_obj, "with_r_removed_objective": obj_r,
                    "delta_objective": delta_obj, "delta_relative_objective": delta_rel,
                    "baseline_feasible": bool(base_feasible), "with_r_removed_feasible": bool(feasible_r),
                    "delta_binary": delta_bin,
                })

        if timing_only and b_idx >= 4:
            break

    timing = {
        "n_cold_solves": len(cold_times),
        "mean_cold_solve_s": float(np.mean(cold_times)) if cold_times else None,
        "n_warm_solves": len(warm_times),
        "mean_warm_solve_s": float(np.mean(warm_times)) if warm_times else None,
        "warm_speedup_x": (float(np.mean(cold_times)) / float(np.mean(warm_times))
                            if cold_times and warm_times and np.mean(warm_times) > 0 else None),
    }
    return pd.DataFrame(background_rows), pd.DataFrame(pair_rows), pd.DataFrame(violation_rows), timing


def summarize_violations(violations_df):
    """Bucket violations by magnitude so a 1e-8 numerical blip and a 0.5 real bug
    are never conflated. Returns (bucket_counts_dict, worst_magnitude)."""
    if len(violations_df) == 0:
        return {"1e-6_to_1e-4": 0, "1e-4_to_1e-2": 0, "gt_1e-2": 0}, 0.0
    mags = violations_df["delta_objective"].abs().to_numpy()
    buckets = {
        "1e-6_to_1e-4": int(((mags >= 1e-6) & (mags < 1e-4)).sum()),
        "1e-4_to_1e-2": int(((mags >= 1e-4) & (mags < 1e-2)).sum()),
        "gt_1e-2": int((mags >= 1e-2).sum()),
    }
    return buckets, float(mags.max())


def summarize_reaction_scores(pair_log_df):
    """PR_objective is the preregistered PRIMARY endpoint (sorted on below).
    PR_binary is secondary. PR_hybrid + its decomposition (p_flip,
    mean_obj_drop_given_flip) is diagnostic. PR_relative_objective is reported
    alongside PR_objective for cross-group comparability later (mean AND
    median, since it can be heavy-tailed when baseline_objective is near
    zero -- see module docstring)."""
    rows = []
    for rid, g in pair_log_df.groupby("reaction_id"):
        n = len(g)
        delta_obj = g["delta_objective"].to_numpy()
        delta_rel = g["delta_relative_objective"].to_numpy()
        delta_bin = g["delta_binary"].to_numpy()
        flip = delta_bin == 1
        p_flip = flip.mean()
        mag_given_flip = delta_obj[flip].mean() if flip.any() else 0.0
        rows.append({
            "reaction_id": rid,
            "n_eligible_backgrounds": n,
            "PR_objective": delta_obj.mean(),
            "PR_objective_se": delta_obj.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
            "PR_relative_objective_mean": delta_rel.mean(),
            "PR_relative_objective_median": float(np.median(delta_rel)),
            "PR_binary": delta_bin.mean(),
            "PR_binary_se": delta_bin.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
            "p_flip": p_flip,
            "mean_obj_drop_given_flip": mag_given_flip,
            "PR_hybrid": p_flip * mag_given_flip,
        })
    # sorted by the preregistered PRIMARY endpoint, not whichever looks best
    return pd.DataFrame(rows).sort_values("PR_objective", ascending=False)


def process_group_pilot(step06, step08, cohort, group, base_model, growth_support, priority_reactions,
                         n_backgrounds, magnitude, obj_frac, seed, rel_epsilon=1e-6, timing_only=False):
    active_set, consensus_def = step06.load_strict_consensus(cohort, group)
    model, unblocked = step08.build_group_model(step06, None, base_model.copy(), active_set, growth_support)

    flux_dir = os.path.join(BASE_DIR, "data", "flux_analysis")
    with open(os.path.join(flux_dir, f"{cohort}__{group}__flux.json")) as f:
        step6_fba = json.load(f)["fba_objective_value"]
    target_biomass = obj_frac * step6_fba

    dropout_candidates = sorted(active_set - growth_support)
    missing = set(priority_reactions) - set(dropout_candidates)
    if missing:
        print(f"  WARNING: {len(missing)} priority reactions not in this group's dropout_candidates "
              f"(likely growth-support-protected or outside this group's active_set) -- skipping them: "
              f"{sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}")
        priority_reactions = [r for r in priority_reactions if r not in missing]

    rng = np.random.default_rng(seed)
    background_log, pair_log, violations, timing = run_paired_backgrounds(
        model, dropout_candidates, priority_reactions, target_biomass,
        n_backgrounds, magnitude, rng, rel_epsilon=rel_epsilon, timing_only=timing_only,
    )
    scores = summarize_reaction_scores(pair_log) if len(pair_log) else pd.DataFrame()
    return background_log, pair_log, violations, scores, timing, target_biomass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", type=str, required=True,
                         help="'cohort:group', e.g. GSE126848:obese_no_NAFLD -- proposed representative "
                              "FAVORABLE pilot group (the one cohort with real, if weak, CV predictive "
                              "signal in the Phase-1 regression screen). Per the staged validation plan, "
                              "run this group first; only if it shows real improvement over the rho~=0 "
                              "baseline should the SAME subset be repeated on a DIFFICULT group (healthy "
                              "or steatosis) before considering anything larger.")
    parser.add_argument("--reaction-subset-csv", type=str, required=True,
                         help="CSV with a reaction_id column -- the small pilot subset S. See "
                              "proposed_pilot_reaction_subset_v2.csv: 44 reactions built from four "
                              "independent sources (NOT selected solely by the old unstable score, so "
                              "the new method isn't validated against reactions the old method already "
                              "picked out) -- old-score extremes (6, comparison only), old-score "
                              "near-zero (6, comparison only), a score-blind random draw (28, the actual "
                              "test of whether this discovers new signal), and biologically interpretable "
                              "calibration reactions from GEM1-step11-biomarker-mapping.md's 5 locked "
                              "axes (subsystem diversity not yet incorporated -- check and top up in the "
                              "real environment where the subsystem mapping file exists).")
    parser.add_argument("--n-backgrounds", type=int, default=120,
                         help="M, number of independent background draws, shared across all reactions in "
                              "the subset. INITIAL PILOT DEPTH, NOT FINAL STATISTICAL DEPTH -- do not "
                              "auto-escalate. Start here, inspect the two-seed comparison, and only "
                              "increase if that inspection justifies it.")
    parser.add_argument("--magnitude", type=float, default=0.3, help="inherited from 08 for comparability")
    parser.add_argument("--obj-frac", type=float, default=0.99, help="inherited from 08 for comparability")
    parser.add_argument("--rel-epsilon", type=float, default=1e-6,
                         help="Floor for PR_relative_objective's denominator (max(baseline_objective, "
                              "this)), to avoid divide-by-near-zero blowup when a background is already "
                              "nearly infeasible on its own.")
    parser.add_argument("--seed", type=int, default=101,
                         help="Background-draw seed. Run twice with two different seeds (e.g. 101 and 202) "
                              "for the cross-seed reproducibility check -- this is the actual validation "
                              "criterion, not a single run.")
    parser.add_argument("--violation-hard-fail-threshold", type=float, default=1e-3,
                         help="If any monotonicity-invariant violation's |delta_objective| exceeds this, "
                              "the run hard-fails AFTER writing all output and printing the full violation "
                              "summary (not silently, and not mid-run before results are saved). Default "
                              "1e-3 is well above ordinary LP solver tolerance (~1e-9 to 1e-6) -- a "
                              "violation this large means investigate a real bug, not average it away.")
    parser.add_argument("--timing-only", action="store_true",
                         help="Run only 5 backgrounds, report cold vs. warm-start solve timing, and exit "
                              "-- do this FIRST, before committing to a full pilot's compute budget.")
    args = parser.parse_args()

    step06 = _import_module("step06_flux_analysis", "06_flux_analysis.py")
    step08 = _import_module("step08_perturbation_testing", "08_perturbation_testing.py")

    cohort, group = args.group.split(":")
    priority_reactions = pd.read_csv(args.reaction_subset_csv)["reaction_id"].tolist()
    expected_eligible = args.n_backgrounds * (1 - args.magnitude)

    print("Loading Human-GEM base model...")
    base_model = step06.load_base_model()
    growth_support = step06.load_growth_support()
    print(f"  Pilot group: {cohort}/{group}. Reaction subset size |S|={len(priority_reactions)}. "
          f"M={args.n_backgrounds} backgrounds. seed={args.seed}. "
          f"Expected eligible backgrounds/reaction ~= {expected_eligible:.0f} (M*(1-magnitude)). "
          f"{'TIMING-ONLY (5 backgrounds)' if args.timing_only else ''}")

    background_log, pair_log, violations, scores, timing, target_biomass = process_group_pilot(
        step06, step08, cohort, group, base_model, growth_support, priority_reactions,
        args.n_backgrounds, args.magnitude, args.obj_frac, args.seed,
        rel_epsilon=args.rel_epsilon, timing_only=args.timing_only,
    )

    print(f"\nTiming: {json.dumps(timing, indent=2)}")
    if args.timing_only:
        implied_total_solves = args.n_backgrounds * (1 + len(priority_reactions) * (1 - args.magnitude))
        implied_wall_s = (timing["mean_cold_solve_s"] or 0) * args.n_backgrounds + \
                          (timing["mean_warm_solve_s"] or 0) * (implied_total_solves - args.n_backgrounds)
        print(f"Implied full-pilot solve count at M={args.n_backgrounds}, |S|={len(priority_reactions)}: "
              f"{implied_total_solves:.0f} solves, ~{implied_wall_s:.1f}s wall-clock "
              f"(extrapolated from this timing sample -- re-check on a larger timing sample if this "
              f"matters for a go/no-go compute decision).")
        return

    tag = f"{cohort}__{group}__seed{args.seed}"
    background_log.to_csv(os.path.join(OUT_DIR, f"{tag}__backgrounds.csv"), index=False)
    pair_log.to_csv(os.path.join(OUT_DIR, f"{tag}__pairs.csv.gz"), index=False, compression="gzip")
    scores.to_csv(os.path.join(OUT_DIR, f"{tag}__reaction_scores.csv"), index=False)
    violations.to_csv(os.path.join(OUT_DIR, f"{tag}__violations.csv"), index=False)
    with open(os.path.join(OUT_DIR, f"{tag}__timing.json"), "w") as f:
        json.dump({"timing": timing, "target_biomass": target_biomass,
                    "n_backgrounds": args.n_backgrounds, "magnitude": args.magnitude,
                    "obj_frac": args.obj_frac, "seed": args.seed,
                    "priority_reactions": priority_reactions}, f, indent=2)

    bucket_counts, worst_mag = summarize_violations(violations)
    print(f"\nMonotonicity-invariant violations (logged, see {tag}__violations.csv): "
          f"total={len(violations)}, by magnitude bucket={bucket_counts}, worst={worst_mag:.2e}")

    print(f"\nWrote scores for {len(scores)} reactions to {tag}__reaction_scores.csv "
          f"(sorted by PR_objective, the preregistered PRIMARY endpoint)")
    print(scores.to_string(index=False))
    print(f"\nRun again with a different --seed (e.g. 202) on the same --reaction-subset-csv, then compare "
          f"PR_objective (primary) / PR_binary (secondary) / PR_hybrid (diagnostic) across the two seed "
          f"runs -- Pearson+Spearman, top-10/20 overlap, per-reaction SE, and whether the same reactions "
          f"repeatedly show nonzero effects -- against the rho~=0 baseline already on file, before "
          f"trusting any of these numbers or considering a difficult-group repeat.")

    if worst_mag > args.violation_hard_fail_threshold:
        print(f"\nHARD FAIL: worst violation magnitude {worst_mag:.2e} exceeds "
              f"--violation-hard-fail-threshold={args.violation_hard_fail_threshold:.2e}. "
              f"All output above was written to disk for post-mortem, but DO NOT TRUST these results "
              f"until the violation is investigated (likely a growth-support leak or model-construction "
              f"bug, not real biology).")
        sys.exit(1)


if __name__ == "__main__":
    main()
