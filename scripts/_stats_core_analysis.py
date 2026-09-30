"""
GEM1 - post-production advanced statistics core (read-only).

Consumes only: data/paired_perturbation_production/*/{manifest,checkpoint}.json
and reaction_scores.csv (10 ACCEPTED groups), data/paired_perturbation_pilot/
(2 groups), and the frozen BIOMARKERS_FULL set (imported, not redefined) from
09_calibration.py. Reuses the production script's own jeffreys_ci() and frozen
constants (TAU, M, GROUP_SEED) via import -- not re-derived. Reaction-level
data only (no pair-shard decompression across all 10 groups: ~335MB gzip /
~20M rows, judged not worth the runtime for this pass -- see limitations in
the written report. This is a documented scope decision, not a hidden one).

Writes ONLY to data/paired_perturbation_statistics/tables/. Never writes to
data/paired_perturbation_production/, data/paired_perturbation_pilot/, or
data/paired_perturbation_synthesis/.

No weight recalibration, no consensus/flux integration, no TAU change, no
resimulation. Descriptive + randomization/permutation inference only -- no
biological interpretation asserted.
"""
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import beta, pearsonr, spearmanr, kendalltau, mannwhitneyu, hypergeom, fisher_exact

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROD_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_production")
PILOT_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_pilot")
OUT_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_statistics", "tables")
os.makedirs(OUT_DIR, exist_ok=True)

RNG = np.random.default_rng(20260929)  # fixed seed, documented, for all randomization below
N_PERM = 5000

GROUPS = [
    "GSE126848__NASH", "GSE126848__healthy", "GSE126848__obese_no_NAFLD",
    "GSE126848__steatosis", "GSE135251__NASH", "GSE135251__healthy",
    "GSE135251__steatosis", "GSE89632__NASH", "GSE89632__healthy",
    "GSE89632__steatosis",
]


def _import(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BASE_DIR, "scripts", filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_prod = _import("step08b", "08b_paired_perturbation_production.py")
_calib = _import("step09_calibration", "09_calibration.py")
jeffreys_ci = _prod.jeffreys_ci
TAU = _prod.TAU
M = _prod.M
BIOMARKERS_FULL = _calib.BIOMARKERS_FULL


def load_group(g):
    gdir = os.path.join(PROD_DIR, g)
    with open(os.path.join(gdir, "manifest.json")) as f:
        manifest = json.load(f)
    assert manifest["status"] == "ACCEPTED"
    scores = pd.read_csv(os.path.join(gdir, "reaction_scores.csv"))
    return manifest, scores


DATA = {g: load_group(g) for g in GROUPS}
print(f"Loaded {len(DATA)} ACCEPTED groups.")

ALL_REACTIONS = sorted(set().union(*[set(s["reaction_id"]) for _, s in DATA.values()]))
N_FULL_NETWORK = len(ALL_REACTIONS)
RXN_IDX = {r: i for i, r in enumerate(ALL_REACTIONS)}
N_R = len(ALL_REACTIONS)
N_G = len(GROUPS)
print(f"Full-network union: {N_R} reactions x {N_G} groups")

# Eligibility (E) and detection (D) boolean matrices, reactions x groups
E = np.zeros((N_R, N_G), dtype=bool)
D = np.zeros((N_R, N_G), dtype=bool)
PR = np.full((N_R, N_G), np.nan)
for gi, g in enumerate(GROUPS):
    _, scores = DATA[g]
    idx = scores["reaction_id"].map(RXN_IDX).values
    E[idx, gi] = True
    D[idx, gi] = scores["occ_obj_k"].values > 0
    PR[idx, gi] = scores["PR_objective"].values

# ===========================================================================
# 1. SPARSITY / OCCURRENCE
# ===========================================================================
rows = []
for gi, g in enumerate(GROUPS):
    manifest, scores = DATA[g]
    n_eval = int(E[:, gi].sum())
    n_ne = N_R - n_eval
    n_flip = int((scores["occ_flip_k"] > 0).sum())
    n_nonzero = int(D[:, gi].sum())
    n_zero = n_eval - n_nonzero

    def jci(k, n):
        lo, hi = jeffreys_ci(k, n)
        return k / n, lo, hi

    p_ne, lo_ne, hi_ne = jci(n_ne, N_R)
    p_nonzero, lo_nonzero, hi_nonzero = jci(n_nonzero, n_eval)
    p_flip, lo_flip, hi_flip = jci(n_flip, n_eval)
    rows.append({
        "group": g, "n_full_network": N_R, "n_evaluated_eligible": n_eval, "n_NE": n_ne,
        "frac_NE_of_full_network": p_ne, "frac_NE_jeffreys_lo": lo_ne, "frac_NE_jeffreys_hi": hi_ne,
        "n_evaluated_zero": n_zero, "n_evaluated_nonzero_delta_gt_TAU": n_nonzero,
        "frac_nonzero_of_evaluated": p_nonzero,
        "frac_nonzero_jeffreys_lo": lo_nonzero, "frac_nonzero_jeffreys_hi": hi_nonzero,
        "n_feasibility_flip": n_flip, "frac_flip_of_evaluated": p_flip,
        "frac_flip_jeffreys_lo": lo_flip, "frac_flip_jeffreys_hi": hi_flip,
    })
t1 = pd.DataFrame(rows)
t1.to_csv(os.path.join(OUT_DIR, "01_sparsity_occurrence_by_group.csv"), index=False)
print("Wrote 01_sparsity_occurrence_by_group.csv")

# ===========================================================================
# 2. POSITIVE-EFFECT MAGNITUDE (PR_objective | occ_obj_k>0), reaction-level
# ===========================================================================
def bootstrap_ci(x, stat_fn, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    boots = np.empty(n_boot)
    n = len(x)
    for b in range(n_boot):
        boots[b] = stat_fn(rng.choice(x, size=n, replace=True))
    return np.percentile(boots, 2.5), np.percentile(boots, 97.5)

rows = []
detected_pr_by_group = {}
for gi, g in enumerate(GROUPS):
    _, scores = DATA[g]
    detected = scores.loc[scores["occ_obj_k"] > 0, "PR_objective"].values
    detected_pr_by_group[g] = detected
    if len(detected) == 0:
        continue
    med_lo, med_hi = bootstrap_ci(detected, np.median, seed=hash(g) % (2**31))
    p90 = np.percentile(detected, 90)
    p90_lo, p90_hi = bootstrap_ci(detected, lambda a: np.percentile(a, 90), seed=(hash(g) + 1) % (2**31))
    rows.append({
        "group": g, "n_detected": len(detected),
        "median": np.median(detected), "median_boot_lo": med_lo, "median_boot_hi": med_hi,
        "q25": np.percentile(detected, 25), "q75": np.percentile(detected, 75),
        "IQR": np.percentile(detected, 75) - np.percentile(detected, 25),
        "p90": p90, "p90_boot_lo": p90_lo, "p90_boot_hi": p90_hi,
        "p95": np.percentile(detected, 95), "max": detected.max(),
    })
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(OUT_DIR, "02_positive_effect_magnitude_by_group.csv"), index=False)
print("Wrote 02_positive_effect_magnitude_by_group.csv")

# ===========================================================================
# 3. CROSS-GROUP RECURRENCE + eligibility-aware permutation null
# ===========================================================================
elig_freq = E.sum(axis=1)
detect_freq = D.sum(axis=1)
t3a = pd.DataFrame({
    "reaction_id": ALL_REACTIONS, "n_groups_eligible": elig_freq, "n_groups_detected": detect_freq,
})
t3a.to_csv(os.path.join(OUT_DIR, "03a_reaction_eligibility_and_detection_frequency.csv"), index=False)

obs_hist = np.bincount(detect_freq, minlength=N_G + 1)

n_evaluated_g = E.sum(axis=0)
n_detected_g = D.sum(axis=0)
eligible_idx_by_group = [np.flatnonzero(E[:, gi]) for gi in range(N_G)]

null_hists = np.zeros((N_PERM, N_G + 1), dtype=int)
for p in range(N_PERM):
    detect_count = np.zeros(N_R, dtype=int)
    for gi in range(N_G):
        pool = eligible_idx_by_group[gi]
        k = n_detected_g[gi]
        chosen = RNG.choice(pool, size=k, replace=False)
        detect_count[chosen] += 1
    null_hists[p] = np.bincount(detect_count, minlength=N_G + 1)

null_mean = null_hists.mean(axis=0)
null_sd = null_hists.std(axis=0, ddof=1)
rows = []
for k in range(N_G + 1):
    obs = obs_hist[k]
    z = (obs - null_mean[k]) / null_sd[k] if null_sd[k] > 0 else np.nan
    # two-sided empirical p: fraction of null draws at least as extreme as observed
    if obs >= null_mean[k]:
        p_emp = 2 * min(np.mean(null_hists[:, k] >= obs), 0.5)
    else:
        p_emp = 2 * min(np.mean(null_hists[:, k] <= obs), 0.5)
    rows.append({
        "n_groups_detected": k, "observed_n_reactions": int(obs),
        "null_mean": null_mean[k], "null_sd": null_sd[k],
        "null_p2_5": np.percentile(null_hists[:, k], 2.5),
        "null_p97_5": np.percentile(null_hists[:, k], 97.5),
        "effect_size_z": z, "empirical_p_two_sided": p_emp,
    })
t3b = pd.DataFrame(rows)
t3b.to_csv(os.path.join(OUT_DIR, "03b_recurrence_null_summary.csv"), index=False)
print(f"Wrote 03a/03b recurrence tables (N_PERM={N_PERM}, per-group eligible-pool "
      f"shuffle preserving each group's exact n_evaluated and n_detected).")

# ===========================================================================
# 4. CROSS-COHORT REPRODUCIBILITY (matched disease states)
# ===========================================================================
STATE_TO_GROUPS = {}
for g in GROUPS:
    cohort, state = g.split("__", 1)
    STATE_TO_GROUPS.setdefault(state, []).append((cohort, g))

rows = []
for state, lst in STATE_TO_GROUPS.items():
    if len(lst) < 2:
        continue
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            (c1, g1), (c2, g2) = lst[i], lst[j]
            gi1, gi2 = GROUPS.index(g1), GROUPS.index(g2)
            joint_elig = E[:, gi1] & E[:, gi2]
            n_joint = int(joint_elig.sum())
            pr1 = PR[joint_elig, gi1]
            pr2 = PR[joint_elig, gi2]
            r_p, p_p = pearsonr(pr1, pr2)
            r_s, p_s = spearmanr(pr1, pr2)
            tau, p_tau = kendalltau(pr1, pr2)

            det1 = D[joint_elig, gi1]
            det2 = D[joint_elig, gi2]
            either = det1 | det2
            n_either = int(either.sum())
            if n_either >= 3:
                r_p_d, p_p_d = pearsonr(pr1[either], pr2[either])
                r_s_d, p_s_d = spearmanr(pr1[either], pr2[either])
            else:
                r_p_d = p_p_d = r_s_d = p_s_d = np.nan

            both = int((det1 & det2).sum())
            only1 = int((det1 & ~det2).sum())
            only2 = int((~det1 & det2).sum())
            neither = n_joint - both - only1 - only2
            jaccard = both / (both + only1 + only2) if (both + only1 + only2) > 0 else np.nan
            table = [[both, only1], [only2, neither]]
            odds_ratio, fisher_p = fisher_exact(table)
            rows.append({
                "state": state, "cohort_A": c1, "cohort_B": c2,
                "n_joint_eligible": n_joint,
                "pearson_r_all_joint": r_p, "pearson_p_all_joint": p_p,
                "spearman_r_all_joint": r_s, "spearman_p_all_joint": p_s,
                "kendall_tau_all_joint": tau, "kendall_p_all_joint": p_tau,
                "n_either_detected": n_either,
                "pearson_r_detected_subset": r_p_d, "spearman_r_detected_subset": r_s_d,
                "n_both_detected": both, "n_only_A_detected": only1, "n_only_B_detected": only2,
                "n_neither_detected": neither, "jaccard_detected_overlap": jaccard,
                "odds_ratio_fisher": odds_ratio, "fisher_exact_p": fisher_p,
            })
t4 = pd.DataFrame(rows)
t4.to_csv(os.path.join(OUT_DIR, "04_cross_cohort_reproducibility.csv"), index=False)
print("Wrote 04_cross_cohort_reproducibility.csv")

# ===========================================================================
# 5. BIOMARKER / CALIBRATION REACTIONS
# ===========================================================================
biomarker_rxns = sorted({rid for lst in BIOMARKERS_FULL.values() for rid in lst})
assert len(biomarker_rxns) == 25
rxn_to_set = {}
for name, lst in BIOMARKERS_FULL.items():
    for rid in lst:
        rxn_to_set.setdefault(rid, []).append(name)

rows = []
for rid in biomarker_rxns:
    ridx = RXN_IDX.get(rid)
    if ridx is None:
        # Not a dropout candidate (NE) in ANY of the 10 groups -- e.g. always
        # inside growth_support / never in any group's active_set. Recorded
        # explicitly rather than silently dropped.
        rows.append({
            "reaction_id": rid, "biomarker_set": "+".join(rxn_to_set[rid]),
            "n_groups_eligible": 0, "n_groups_detected": 0,
            "mean_PR_objective_when_eligible": None, "median_PR_objective_when_eligible": None,
        })
        continue
    n_elig = int(E[ridx].sum())
    n_det = int(D[ridx].sum())
    pr_vals = PR[ridx, E[ridx]]
    rows.append({
        "reaction_id": rid, "biomarker_set": "+".join(rxn_to_set[rid]),
        "n_groups_eligible": n_elig, "n_groups_detected": n_det,
        "mean_PR_objective_when_eligible": np.mean(pr_vals) if n_elig else None,
        "median_PR_objective_when_eligible": np.median(pr_vals) if n_elig else None,
    })
t5a = pd.DataFrame(rows)
t5a.to_csv(os.path.join(OUT_DIR, "05a_biomarker_per_reaction_summary.csv"), index=False)

# Per-group eligibility-aware permutation: draw m_g random reactions from that
# group's OWN eligible pool (not the biomarker set), 5000 reps, compare mean
# PR_objective and detection count of the real biomarker-eligible subset.
rows = []
mw_rows = []
for gi, g in enumerate(GROUPS):
    _, scores = DATA[g]
    elig_mask = E[:, gi]
    bio_mask = np.zeros(N_R, dtype=bool)
    for rid in biomarker_rxns:
        ridx = RXN_IDX.get(rid)
        if ridx is not None:
            bio_mask[ridx] = True
    bio_elig_mask = bio_mask & elig_mask
    m_g = int(bio_elig_mask.sum())
    if m_g == 0:
        continue
    pool_idx = np.flatnonzero(elig_mask)
    obs_mean = PR[bio_elig_mask, gi].mean()
    obs_ndet = int(D[bio_elig_mask, gi].sum())

    null_means = np.empty(N_PERM)
    null_ndets = np.empty(N_PERM, dtype=int)
    pr_pool = PR[pool_idx, gi]
    det_pool = D[pool_idx, gi]
    for p in range(N_PERM):
        chosen = RNG.choice(len(pool_idx), size=m_g, replace=False)
        null_means[p] = pr_pool[chosen].mean()
        null_ndets[p] = det_pool[chosen].sum()

    z_mag = (obs_mean - null_means.mean()) / null_means.std(ddof=1)
    p_mag = 2 * min(np.mean(null_means >= obs_mean), np.mean(null_means <= obs_mean))
    p_mag = min(p_mag, 1.0)

    n_pool = len(pool_idx)
    K_pool_detected = int(det_pool.sum())
    p_hyper = hypergeom.sf(obs_ndet - 1, n_pool, K_pool_detected, m_g)

    rows.append({
        "group": g, "m_g_biomarker_eligible": m_g, "n_pool_eligible": n_pool,
        "observed_mean_PR_objective": obs_mean, "null_mean_of_means": null_means.mean(),
        "null_sd_of_means": null_means.std(ddof=1), "effect_size_z_magnitude": z_mag,
        "permutation_p_magnitude_two_sided": p_mag,
        "observed_n_detected": obs_ndet, "pool_n_detected": K_pool_detected,
        "hypergeom_p_detection_upper_tail": p_hyper,
    })

    full_pr = scores["PR_objective"].values
    bio_pr = PR[bio_elig_mask, gi]
    if m_g >= 2:
        stat, p_mw = mannwhitneyu(bio_pr, full_pr, alternative="two-sided")
        mw_rows.append({"group": g, "mannwhitney_U": stat, "mannwhitney_p_uncorrected": p_mw})

t5b = pd.DataFrame(rows)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked_mono = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked_mono, 0, 1)
    return q

t5b["BH_q_magnitude"] = bh_fdr(t5b["permutation_p_magnitude_two_sided"].values)
t5b["BH_q_detection"] = bh_fdr(t5b["hypergeom_p_detection_upper_tail"].values)
t5b.to_csv(os.path.join(OUT_DIR, "05b_biomarker_permutation_by_group.csv"), index=False)

t5_mw = pd.DataFrame(mw_rows)
t5_mw.to_csv(os.path.join(OUT_DIR, "05b_exploratory_mannwhitney_only.csv"), index=False)

# Stouffer combination across groups (equal weight)
z_vals = t5b["effect_size_z_magnitude"].dropna().values
from scipy.stats import norm
combined_z = z_vals.sum() / np.sqrt(len(z_vals))
combined_p = 2 * (1 - norm.cdf(abs(combined_z)))
pd.DataFrame([{
    "method": "Stouffer combination of per-group eligibility-aware permutation z-scores (equal weight)",
    "n_groups_combined": len(z_vals), "combined_Z": combined_z, "combined_p_two_sided": combined_p,
}]).to_csv(os.path.join(OUT_DIR, "05c_biomarker_combined_stouffer.csv"), index=False)
print("Wrote 05a/05b/05c biomarker tables (permutation-based primary inference; "
      "Mann-Whitney retained as exploratory-only in 05b_exploratory_mannwhitney_only.csv)")

# ===========================================================================
# 6. MONTE CARLO PRECISION
# ===========================================================================
rows = []
for gi, g in enumerate(GROUPS):
    _, scores = DATA[g]
    det = scores[scores["occ_obj_k"] > 0].copy()
    det["cv"] = det["PR_objective_se"] / det["PR_objective"].abs()
    jw = det["occ_obj_jeffreys_hi"] - det["occ_obj_jeffreys_lo"]
    rare = det["occ_obj_k"] <= 5
    well = det["occ_obj_k"] >= 50
    rows.append({
        "group": g, "n_detected": len(det),
        "median_CV_SE_over_|PR_objective|": det["cv"].median(),
        "q25_CV": det["cv"].quantile(0.25), "q75_CV": det["cv"].quantile(0.75),
        "n_rare_event_k_le_5": int(rare.sum()),
        "median_jeffreys_width_rare_k_le_5": jw[rare].median() if rare.any() else None,
        "n_well_characterized_k_ge_50": int(well.sum()),
        "median_jeffreys_width_well_k_ge_50": jw[well].median() if well.any() else None,
        "M_backgrounds": M,
    })
t6 = pd.DataFrame(rows)
t6.to_csv(os.path.join(OUT_DIR, "06_monte_carlo_precision.csv"), index=False)
print("Wrote 06_monte_carlo_precision.csv")

# ===========================================================================
# 7. PILOT -> PRODUCTION VALIDATION
# ===========================================================================
PILOT_GROUPS = {
    "GSE126848__obese_no_NAFLD": [101, 202],
    "GSE135251__healthy": [101, 202, 303, 404],
}
corr_rows = []
ba_rows = []
for g, seeds in PILOT_GROUPS.items():
    _, prod_scores = DATA[g]
    frames = []
    for s in seeds:
        p = os.path.join(PILOT_DIR, f"{g}__seed{s}__reaction_scores.csv")
        df = pd.read_csv(p)
        frames.append(df[["reaction_id", "PR_objective"]])
    pilot_pooled = (pd.concat(frames, ignore_index=True).groupby("reaction_id")["PR_objective"]
                     .mean().reset_index().rename(columns={"PR_objective": "PR_objective_pilot"}))
    merged = prod_scores.merge(pilot_pooled, on="reaction_id", how="inner")
    a, b = merged["PR_objective"].values, merged["PR_objective_pilot"].values
    r_p, p_p = pearsonr(a, b)
    r_s, p_s = spearmanr(a, b)
    tau, p_tau = kendalltau(a, b)
    diff = a - b
    mean_ab = (a + b) / 2
    bias = diff.mean()
    sd_diff = diff.std(ddof=1)
    loa_lo, loa_hi = bias - 1.96 * sd_diff, bias + 1.96 * sd_diff
    corr_rows.append({
        "group": g, "n_overlap": len(merged),
        "pearson_r": r_p, "pearson_p": p_p, "spearman_r": r_s, "spearman_p": p_s,
        "kendall_tau": tau, "kendall_p": p_tau,
        "mean_abs_diff": np.abs(diff).mean(), "median_abs_diff": np.median(np.abs(diff)),
        "bland_altman_bias": bias, "bland_altman_sd_diff": sd_diff,
        "bland_altman_LoA_lo": loa_lo, "bland_altman_LoA_hi": loa_hi,
    })
    for rid, m, d, pv, ppv in zip(merged["reaction_id"], mean_ab, diff, a, b):
        ba_rows.append({"group": g, "reaction_id": rid, "mean_prod_pilot": m, "diff_prod_minus_pilot": d,
                         "PR_objective_prod": pv, "PR_objective_pilot": ppv})
t7a = pd.DataFrame(corr_rows)
t7a.to_csv(os.path.join(OUT_DIR, "07a_pilot_vs_production_correlation.csv"), index=False)
t7b = pd.DataFrame(ba_rows)
t7b.to_csv(os.path.join(OUT_DIR, "07b_bland_altman_source.csv"), index=False)
print("Wrote 07a_pilot_vs_production_correlation.csv, 07b_bland_altman_source.csv")

print("\nAll core statistics tables written to:", OUT_DIR)
print("N_PERM =", N_PERM, "| RNG seed = 20260929 | fixed & documented")
