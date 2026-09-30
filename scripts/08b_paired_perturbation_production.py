"""
GEM1 - Production paired-counterfactual perturbation scoring.

FROZEN SPECIFICATION v1.0 (see accompanying audit; nothing below is tunable from
the CLI except --group and --workers-note -- M, magnitude, tau, and obj_frac are
hardcoded constants, deliberately not exposed as flags, so a production run
cannot silently drift from the frozen decision):

    TAU = 1e-6      (numerical noise floor; theory/numerics-first, from LP solver
                      tolerance ~1e-9-1e-6 -- NOT derived from pilot data. The
                      pilot's observed real-vs-noise gap is supporting evidence,
                      not the derivation.)
    M = 800         (frozen production background budget)
    MAGNITUDE = 0.30 (shared Bernoulli(q=0.30) background draw, same distribution
                      the pilot and 08_perturbation_testing.py both use)
    OBJ_FRAC = 0.99  (feasibility threshold, unchanged from pilot/production)

Scores computed: PR_objective (expected impact, PRIMARY), PR_relative_objective
(mean+median, secondary/comparability), PR_binary (occurrence, empirical k/n with
Jeffreys 95% CI, SECONDARY), PR_hybrid + its decomposition p_flip /
mean_obj_drop_given_flip (conditional severity, SECONDARY/diagnostic). Legacy
0.15/0.15/0.70 combination weights remain SUSPENDED -- not used, not
recalibrated, not referenced anywhere in this script. This script does not
write to, or read from, any consensus/flux/confidence-score file -- the new
score is not integrated with anything yet.

Reaction subset: ALL dropout-eligible reactions for the group (active_set -
growth_support), i.e. the full network, not a curated subset -- matches the
scale validated by the engineering check (../_engineering_scaling_check.py).

Monte Carlo uncertainty scope (stated explicitly, also written into every
group's manifest.json): SE/CI reported here reflect background-sampling
(Monte Carlo) variance CONDITIONAL ON the already-constructed, fixed
context-specific group model (fixed active_set, fixed growth_support, fixed
consensus algorithm output). They do NOT capture model-construction
uncertainty, consensus-algorithm uncertainty, or cross-cohort uncertainty.

Status handling (the reason this is a new script, not a patch to the pilot):
  - "optimal"   -> raw objective_value used as-is.
  - "infeasible" -> a real, EXPECTED LP outcome (this IS part of the occurrence
                    facet's signal, not a failure). raw objective_value is None
                    (never fabricated). For SCORING ONLY, a separate
                    scoring_objective field applies the frozen 0.0-floor
                    convention (zero growth capacity) -- this convention is
                    inherited unchanged from the pilot's PR_objective
                    definition, applied explicitly and logged, never silently.
  - anything else (unbounded, infeasible_or_unbounded, time-limit, solver
                    exception, etc.) -> "unexpected". raw AND scoring objective
                    are both None. The pair/background is EXCLUDED from every
                    statistic (not floored to 0, not counted as occurrence-negative
                    or occurrence-positive) and logged to unexpected_status.csv.
                    A baseline (cold-solve) unexpected status excludes the ENTIRE
                    background from every reaction that run, not just one pair.

Checkpoint/resume: each background's draw is generated from an independent
child SeedSequence (np.random.SeedSequence(group_seed).spawn(M)[b_idx]), so
background b_idx's random draw is IDENTICAL regardless of run/resume history --
resuming after a crash reproduces exactly the same backgrounds already on disk
plus the remaining ones, deterministically. A background is "done" iff both its
background shard and pair shard exist and are non-empty; the run scans for
missing indices and only (re)computes those.

Output layout (data/paired_perturbation_production/{cohort}__{group}/):
  manifest.json                 - provenance, config, status (RUNNING/ACCEPTED/REJECTED)
  checkpoint.json                - completed background indices
  backgrounds/bg_{idx:04d}.csv   - 1 row: this background's raw+scoring cold-solve result
  pairs/bg_{idx:04d}.csv.gz      - 1 row per eligible reaction tested against this background
  unexpected_status.csv          - every non-{optimal,infeasible} status observed, with context
  violations.csv                 - monotonicity-invariant violations (diagnostic, logged not fatal
                                    below --violation-hard-fail-threshold)
  reaction_scores.csv            - FINAL derived summary, written only after all M backgrounds
                                    complete AND integrity checks pass (see run_integrity_checks)
  timing_summary.json             - aggregate + per-background timing telemetry

All per-background/per-pair writes are atomic (write to .tmp, then os.replace) so
an interruption mid-write never leaves a corrupt shard that looks complete.
"""
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import beta

SCRIPTS_DIR = os.path.dirname(__file__)
BASE_DIR = os.path.dirname(SCRIPTS_DIR)
OUT_ROOT = os.path.join(BASE_DIR, "data", "paired_perturbation_production")

# ---- FROZEN CONSTANTS (not CLI-configurable, by design) --------------------
TAU = 1e-6
M = 800
MAGNITUDE = 0.30
OBJ_FRAC = 0.99
EXPECTED_STATUSES = {"optimal", "infeasible"}
VIOLATION_HARD_FAIL_THRESHOLD = 1e-3
UNEXPECTED_STATUS_HARD_FAIL_RATE = 1e-3  # >0.1% unexpected statuses -> reject group
ELIGIBILITY_TOLERANCE = 0.05  # |observed - (1-MAGNITUDE)| must stay within this

# Fixed, disjoint from pilot seeds (101/202/303/404) and the engineering seed
# (999001). Alphabetical by cohort then group -- documented here, never derived
# from a non-reproducible source like Python's hash().
GROUP_SEED = {
    ("GSE126848", "NASH"): 800001,
    ("GSE126848", "healthy"): 800002,
    ("GSE126848", "obese_no_NAFLD"): 800003,
    ("GSE126848", "steatosis"): 800004,
    ("GSE135251", "NASH"): 800005,
    ("GSE135251", "healthy"): 800006,
    ("GSE135251", "steatosis"): 800007,
    ("GSE89632", "NASH"): 800008,
    ("GSE89632", "healthy"): 800009,
    ("GSE89632", "steatosis"): 800010,
}


def _import_module(name, filename):
    path = os.path.join(SCRIPTS_DIR, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _atomic_write_csv(df, path, **to_csv_kwargs):
    tmp = path + ".tmp"
    df.to_csv(tmp, index=False, **to_csv_kwargs)
    os.replace(tmp, path)


def _atomic_write_json(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2, default=str)
    os.replace(tmp, path)


def _validate_shard(b_idx, bg_path, pair_path):
    """Read and cross-check a background's shard pair BEFORE trusting it as
    complete during resume -- existence alone is not validation. Returns False
    (triggering deterministic recomputation of this background) on any
    inconsistency, including a corrupted/truncated file."""
    try:
        bg_row = pd.read_csv(bg_path)
        pairs = pd.read_csv(pair_path)
    except Exception:
        return False
    if len(bg_row) != 1:
        return False
    bg_row = bg_row.iloc[0]
    if int(bg_row["background"]) != b_idx:
        return False
    if len(pairs) != int(bg_row["n_pairs_computed"]):
        return False
    return True


def jeffreys_ci(x, n, alpha=0.05):
    lo = 0.0 if x <= 0 else beta.ppf(alpha / 2, x + 0.5, n - x + 0.5)
    hi = 1.0 if x >= n else beta.ppf(1 - alpha / 2, x + 0.5, n - x + 0.5)
    return lo, hi


def get_git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=BASE_DIR, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return None


def build_manifest(cohort, group, seed, dropout_candidates, active_set_path, growth_support_path):
    return {
        "cohort": cohort, "group": group, "seed": seed,
        "config": {"TAU": TAU, "M": M, "magnitude": MAGNITUDE, "obj_frac": OBJ_FRAC,
                   "legacy_weights_suspended": True},
        "n_dropout_candidates": len(dropout_candidates),
        "provenance": {
            "git_commit": get_git_commit(),
            "script": os.path.relpath(__file__, BASE_DIR),
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
            "active_set_source": active_set_path, "growth_support_source": growth_support_path,
        },
        "monte_carlo_uncertainty_scope": (
            "SE/CI reflect background-sampling (Monte Carlo) variance conditional on the "
            "fixed, already-constructed context-specific group model (fixed active_set, "
            "fixed growth_support, fixed consensus algorithm output). Does NOT capture "
            "model-construction, consensus-algorithm, or cross-cohort uncertainty."
        ),
        "status": "RUNNING",
        "rejection_reasons": [],
    }


def solve_and_classify(model, target_biomass):
    """Returns (status, raw_objective, scoring_objective, feasible).
    scoring_objective applies the frozen 0.0-floor ONLY for status=='infeasible'
    (an expected outcome). Any other non-optimal status yields
    scoring_objective=None and must be excluded from every statistic upstream."""
    try:
        sol = model.optimize()
        status = sol.status
    except Exception as e:
        return f"exception:{type(e).__name__}", None, None, False

    if status == "optimal":
        raw_obj = sol.objective_value
        feasible = raw_obj >= target_biomass - 1e-6
        return status, raw_obj, raw_obj, feasible
    elif status == "infeasible":
        return status, None, 0.0, False
    else:
        return status, None, None, False


def run_group(cohort, group, resume=True):
    step06 = _import_module("step06_flux_analysis", "06_flux_analysis.py")
    step08 = _import_module("step08_perturbation_testing", "08_perturbation_testing.py")

    out_dir = os.path.join(OUT_ROOT, f"{cohort}__{group}")
    bg_dir = os.path.join(out_dir, "backgrounds")
    pair_dir = os.path.join(out_dir, "pairs")
    os.makedirs(bg_dir, exist_ok=True)
    os.makedirs(pair_dir, exist_ok=True)

    manifest_path = os.path.join(out_dir, "manifest.json")
    checkpoint_path = os.path.join(out_dir, "checkpoint.json")
    unexpected_path = os.path.join(out_dir, "unexpected_status.csv")
    violations_path = os.path.join(out_dir, "violations.csv")

    print(f"Loading base model + growth support for {cohort}/{group}...", flush=True)
    base_model = step06.load_base_model()
    growth_support = step06.load_growth_support()
    active_set, consensus_def = step06.load_strict_consensus(cohort, group)
    model, unblocked = step08.build_group_model(step06, None, base_model.copy(), active_set, growth_support)
    dropout_candidates = sorted(active_set - growth_support)
    n_candidates = len(dropout_candidates)
    candidates_arr = np.array(dropout_candidates)

    flux_dir = os.path.join(BASE_DIR, "data", "flux_analysis")
    with open(os.path.join(flux_dir, f"{cohort}__{group}__flux.json")) as f:
        step6_fba = json.load(f)["fba_objective_value"]
    target_biomass = OBJ_FRAC * step6_fba

    seed = GROUP_SEED[(cohort, group)]
    child_seeds = np.random.SeedSequence(seed).spawn(M)

    if resume and os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
        if manifest["n_dropout_candidates"] != n_candidates:
            raise SystemExit(
                f"Config drift detected for {cohort}/{group}: manifest recorded "
                f"{manifest['n_dropout_candidates']} dropout candidates, current environment "
                f"computes {n_candidates}. Refusing to resume with a mismatched candidate pool "
                f"(likely growth_support/consensus file changed since this run started)."
            )
        print(f"  Resuming existing run (status={manifest['status']}).", flush=True)
    else:
        manifest = build_manifest(
            cohort, group, seed, dropout_candidates,
            active_set_path=f"data/consensus_scores/{cohort}__{group}__consensus.json",
            growth_support_path="data/context_specific_models/minimal_growth_support.json",
        )
        _atomic_write_json(manifest, manifest_path)

    if resume and os.path.exists(checkpoint_path):
        with open(checkpoint_path) as f:
            checkpoint = json.load(f)
    else:
        checkpoint = {"completed_backgrounds": []}
    completed = set(checkpoint["completed_backgrounds"])

    unexpected_rows = [] if not os.path.exists(unexpected_path) else pd.read_csv(unexpected_path).to_dict("records")
    violation_rows = [] if not os.path.exists(violations_path) else pd.read_csv(violations_path).to_dict("records")

    print(f"  dropout_candidates={n_candidates}. M={M}. seed={seed}. "
          f"Resuming with {len(completed)}/{M} backgrounds already complete.", flush=True)

    run_wall_start = time.time()
    for b_idx in range(M):
        bg_path = os.path.join(bg_dir, f"bg_{b_idx:04d}.csv")
        pair_path = os.path.join(pair_dir, f"bg_{b_idx:04d}.csv.gz")
        if (b_idx in completed and os.path.exists(bg_path) and os.path.exists(pair_path)
                and _validate_shard(b_idx, bg_path, pair_path)):
            continue
        if b_idx in completed:
            print(f"  WARNING: background {b_idx} was marked complete but failed shard "
                  f"validation -- recomputing deterministically (same seed, so this is safe).",
                  flush=True)

        rng = np.random.default_rng(child_seeds[b_idx])
        mask = rng.random(n_candidates) < MAGNITUDE
        removed = candidates_arr[mask]
        removed_set = set(removed)
        eligible = [r for r in dropout_candidates if r not in removed_set]

        bg_wall_start = time.time()
        with model:
            for rid in removed:
                model.reactions.get_by_id(rid).bounds = (0, 0)

            t0 = time.time()
            base_status, base_raw_obj, base_scoring_obj, base_feasible = solve_and_classify(model, target_biomass)
            t_cold = time.time() - t0

            if base_status not in EXPECTED_STATUSES:
                unexpected_rows.append({
                    "background": b_idx, "reaction_id": None, "context": "baseline",
                    "status": base_status, "n_removed": int(mask.sum()),
                })
                # Baseline itself untrustworthy -> exclude this ENTIRE background,
                # not just one pair. Still write shards (empty pairs) so the
                # background index is marked done and resume logic is simple.
                pair_rows = []
            else:
                pair_rows = []
                for position, rid in enumerate(eligible):
                    with model:
                        model.reactions.get_by_id(rid).bounds = (0, 0)
                        t1 = time.time()
                        status_r, raw_obj_r, scoring_obj_r, feasible_r = solve_and_classify(model, target_biomass)
                        t_warm = time.time() - t1

                    if status_r not in EXPECTED_STATUSES:
                        unexpected_rows.append({
                            "background": b_idx, "reaction_id": rid, "context": "warm",
                            "status": status_r, "n_removed": int(mask.sum()),
                        })
                        pair_rows.append({
                            "background": b_idx, "position": position, "reaction_id": rid,
                            "status": status_r, "raw_objective": None, "scoring_objective": None,
                            "delta_objective": None, "delta_relative_objective": None,
                            "delta_binary": None, "monotonicity_violation": False,
                            "solve_time_s": t_warm, "excluded_unexpected": True,
                        })
                        continue

                    delta_obj = base_scoring_obj - scoring_obj_r
                    delta_bin = int(base_feasible) - int(feasible_r)
                    delta_rel = delta_obj / max(base_scoring_obj, 1e-6)
                    violation = delta_obj < -TAU or delta_bin == -1
                    if violation:
                        violation_rows.append({
                            "background": b_idx, "reaction_id": rid,
                            "baseline_scoring_objective": base_scoring_obj,
                            "with_r_removed_scoring_objective": scoring_obj_r,
                            "delta_objective": delta_obj, "delta_binary": delta_bin,
                        })
                    pair_rows.append({
                        "background": b_idx, "position": position, "reaction_id": rid,
                        "status": status_r, "raw_objective": raw_obj_r, "scoring_objective": scoring_obj_r,
                        "delta_objective": delta_obj, "delta_relative_objective": delta_rel,
                        "delta_binary": delta_bin, "monotonicity_violation": violation,
                        "solve_time_s": t_warm, "excluded_unexpected": False,
                    })

        bg_wall_s = time.time() - bg_wall_start
        bg_row = pd.DataFrame([{
            "background": b_idx, "seed_entropy": str(child_seeds[b_idx].entropy),
            "n_removed": int(mask.sum()), "n_eligible": len(eligible),
            "eligibility_fraction": len(eligible) / n_candidates,
            "baseline_status": base_status, "baseline_raw_objective": base_raw_obj,
            "baseline_scoring_objective": base_scoring_obj, "baseline_feasible": base_feasible,
            "cold_solve_time_s": t_cold, "wall_clock_s": bg_wall_s,
            "n_pairs_computed": len(pair_rows),
        }])
        _atomic_write_csv(bg_row, bg_path)
        _atomic_write_csv(pd.DataFrame(pair_rows), pair_path, compression="gzip")

        completed.add(b_idx)
        checkpoint["completed_backgrounds"] = sorted(completed)
        _atomic_write_json(checkpoint, checkpoint_path)
        if unexpected_rows:
            _atomic_write_csv(pd.DataFrame(unexpected_rows), unexpected_path)
        if violation_rows:
            _atomic_write_csv(pd.DataFrame(violation_rows), violations_path)

        print(f"  [{time.strftime('%H:%M:%S')}] background {b_idx}/{M} done "
              f"(status={base_status}, n_eligible={len(eligible)}, wall={bg_wall_s:.1f}s, "
              f"{len(completed)}/{M} complete)", flush=True)

    total_wall_s = time.time() - run_wall_start
    print(f"\nAll backgrounds attempted. Total wall this session: {total_wall_s:.1f}s. "
          f"Running integrity checks...", flush=True)

    accepted, reasons, scores = run_integrity_checks_and_summarize(
        out_dir, manifest, checkpoint, n_candidates,
    )
    manifest["status"] = "ACCEPTED" if accepted else "REJECTED"
    manifest["rejection_reasons"] = reasons
    manifest["completed_at_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    _atomic_write_json(manifest, manifest_path)

    if accepted:
        _atomic_write_csv(scores, os.path.join(out_dir, "reaction_scores.csv"))
        print(f"ACCEPTED. reaction_scores.csv written ({len(scores)} reactions).")
    else:
        print(f"REJECTED. Reasons: {reasons}. reaction_scores.csv NOT written. "
              f"Raw background/pair shards preserved for post-mortem.")


def run_integrity_checks_and_summarize(out_dir, manifest, checkpoint, n_candidates):
    """Reads back every shard (never trusts in-memory state alone), recomputes
    final statistics from raw data, and decides ACCEPT/REJECT. Returns
    (accepted: bool, reasons: list[str], scores_df_or_None)."""
    reasons = []
    bg_dir = os.path.join(out_dir, "backgrounds")
    pair_dir = os.path.join(out_dir, "pairs")

    completed = sorted(checkpoint["completed_backgrounds"])
    if len(completed) != M or completed != list(range(M)):
        missing = sorted(set(range(M)) - set(completed))
        reasons.append(f"incomplete run: {len(completed)}/{M} backgrounds, missing {missing[:10]}...")
        return False, reasons, None

    bg_frames, pair_frames = [], []
    for b_idx in range(M):
        bp = os.path.join(bg_dir, f"bg_{b_idx:04d}.csv")
        pp = os.path.join(pair_dir, f"bg_{b_idx:04d}.csv.gz")
        if not (os.path.exists(bp) and os.path.exists(pp)):
            reasons.append(f"missing shard for background {b_idx}")
            continue
        bg_frames.append(pd.read_csv(bp))
        pair_frames.append(pd.read_csv(pp))
    if reasons:
        return False, reasons, None

    bg_df = pd.concat(bg_frames, ignore_index=True)
    pair_df = pd.concat(pair_frames, ignore_index=True) if pair_frames else pd.DataFrame()

    mean_elig = bg_df["eligibility_fraction"].mean()
    if abs(mean_elig - (1 - MAGNITUDE)) > ELIGIBILITY_TOLERANCE:
        reasons.append(f"eligibility fraction {mean_elig:.4f} deviates from expected "
                        f"{1 - MAGNITUDE:.4f} by more than {ELIGIBILITY_TOLERANCE}")

    # NOTE: sum n_pairs_computed, not n_eligible -- these differ whenever a
    # background's baseline was itself "unexpected" (n_eligible reflects the
    # computed eligible list, but zero pair rows are written in that branch).
    total_pairs_expected = int(bg_df["n_pairs_computed"].sum())
    total_pairs_actual = len(pair_df)
    if total_pairs_actual != total_pairs_expected:
        reasons.append(f"pair row count mismatch: expected {total_pairs_expected}, got {total_pairs_actual}")

    n_baseline_unexpected = (~bg_df["baseline_status"].isin(EXPECTED_STATUSES)).sum()
    n_pair_unexpected = pair_df["excluded_unexpected"].sum() if "excluded_unexpected" in pair_df else 0
    total_attempted = total_pairs_actual + n_baseline_unexpected
    unexpected_rate = (n_baseline_unexpected + n_pair_unexpected) / max(total_attempted, 1)
    if unexpected_rate > UNEXPECTED_STATUS_HARD_FAIL_RATE:
        reasons.append(f"unexpected-status rate {unexpected_rate:.4%} exceeds threshold "
                        f"{UNEXPECTED_STATUS_HARD_FAIL_RATE:.4%}")

    valid_pairs = pair_df[~pair_df["excluded_unexpected"]] if len(pair_df) else pair_df
    if len(valid_pairs):
        worst_violation = valid_pairs.loc[valid_pairs["monotonicity_violation"], "delta_objective"].abs().max()
        worst_violation = 0.0 if pd.isna(worst_violation) else worst_violation
        if worst_violation > VIOLATION_HARD_FAIL_THRESHOLD:
            reasons.append(f"worst monotonicity violation {worst_violation:.2e} exceeds "
                            f"hard-fail threshold {VIOLATION_HARD_FAIL_THRESHOLD:.2e}")

    if reasons:
        return False, reasons, None

    rows = []
    for rid, g in valid_pairs.groupby("reaction_id"):
        n = len(g)
        delta_obj = g["delta_objective"].to_numpy()
        delta_rel = g["delta_relative_objective"].to_numpy()
        delta_bin = g["delta_binary"].to_numpy()

        # --- Secondary structure (1): CONTINUOUS above-tau occurrence.
        # Distinct from the binary-flip structure below -- do not conflate.
        # This is what "healthy" needed: sparse real effects that never cross
        # the hard feasibility threshold are invisible to PR_binary/p_flip but
        # visible here.
        above_tau = delta_obj > TAU
        k_obj = int(above_tau.sum())
        p_obj = k_obj / n
        occ_obj_lo, occ_obj_hi = jeffreys_ci(k_obj, n)
        m_obj = float(delta_obj[above_tau].mean()) if above_tau.any() else 0.0

        # --- Secondary structure (2): BINARY feasibility-flip occurrence,
        # plus its conditional-severity diagnostic. Distinct from (1).
        flip = delta_bin == 1
        k_flip = int(flip.sum())
        p_flip = k_flip / n
        occ_flip_lo, occ_flip_hi = jeffreys_ci(k_flip, n)
        mag_given_flip = delta_obj[flip].mean() if flip.any() else 0.0

        rows.append({
            "reaction_id": rid, "n_eligible": n,
            # PRIMARY: raw, unthresholded mean -- exact frozen v1.0 definition.
            "PR_objective": delta_obj.mean(),
            "PR_objective_se": delta_obj.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
            "PR_relative_objective_mean": delta_rel.mean(),
            "PR_relative_objective_median": float(np.median(delta_rel)),
            # secondary (1): continuous above-tau occurrence
            "occ_obj_k": k_obj, "occ_obj_p": p_obj,
            "occ_obj_jeffreys_lo": occ_obj_lo, "occ_obj_jeffreys_hi": occ_obj_hi,
            "occ_obj_m": m_obj,
            # secondary (2): binary feasibility-flip occurrence + conditional severity
            "occ_flip_k": k_flip, "occ_flip_p": p_flip,
            "occ_flip_jeffreys_lo": occ_flip_lo, "occ_flip_jeffreys_hi": occ_flip_hi,
            "occ_flip_mean_obj_drop_given_flip": mag_given_flip,
            "PR_hybrid": p_flip * mag_given_flip,
        })
    scores = pd.DataFrame(rows).sort_values("PR_objective", ascending=False)
    return True, [], scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", required=True, help="'cohort:group', e.g. GSE126848:obese_no_NAFLD")
    parser.add_argument("--no-resume", action="store_true", help="ignore any existing checkpoint (DESTRUCTIVE -- confirm before use)")
    args = parser.parse_args()
    cohort, group = args.group.split(":")
    if (cohort, group) not in GROUP_SEED:
        raise SystemExit(f"Unknown group {cohort}:{group}. Known: {list(GROUP_SEED)}")
    run_group(cohort, group, resume=not args.no_resume)


if __name__ == "__main__":
    main()
