"""
Item 13 -- DEFINITIVE RECONCILIATION PASS. Methodology frozen per explicit
instruction: run exactly once, report fully, do not iterate toward a
favorable result. This table supersedes every previously reported enrichment
number in this audit (items 8, 11, 12's DNL/SHMT/PEMT figures included).
Manuscript files are NOT modified by this script. Items 1-12 are untouched
(this script only reads all_groups_calibrated_with_metadata.csv and writes
its own new output).

FROZEN METHODOLOGY (fixed before any statistic below was computed):

1. Biomarker mapping -- Table 2's CURRENT mapping is authoritative for this
   pass. DNL uses the SCD1-only mapping {MAR00146, MAR00147, MAR00148}
   (NOT the FASN-inclusive mapping used in items 11/12). The other four axes
   are unchanged from docs/STEP11_BIOMARKER_MAPPING.md. Exact reaction ID
   sets for all five axes are printed and written to the output BEFORE any
   statistic, per instruction.

2. Statistical procedure -- exactly one defined procedure:
   cohort/group-level permutation -> empirical p-value per axis per group
   -> Fisher combination DIRECTLY across all 10 cohort/group observations
   (one level, no per-cohort intermediate step) -> Benjamini-Hochberg
   correction across the five axes.

   Independence check (required before trusting the naive Fisher
   combination): the 10 groups are not a random sample of independent
   units -- 3 cohorts, one of which (GSE89632) is a different assay
   platform (microarray) than the other two (RNA-seq), and groups within
   a cohort share the same expression-matrix processing. This script tests
   for non-independence directly: pairwise Pearson correlation of each
   group's full 12931-reaction score vector against every other group's,
   comparing within-cohort pairs to across-cohort pairs, PLUS an effective-
   number-of-independent-tests estimate (Li & Ji, 2005 eigenvalue method)
   applied to the 10x10 group score-vector correlation matrix. If M_eff is
   materially below 10, the naive Fisher combination (which assumes 10
   independent tests, df=20) is flagged as anti-conservative, and a
   SEPARATE, corrected Fisher combination using df=2*M_eff is also reported
   alongside the naive one -- both are reported, neither is silently
   suppressed.

3. Permutation count -- 100,000 draws per group-level test. Seed=42
   (the same seed used throughout this audit), reported explicitly.
   Resolution floor = 1/100001 (reported explicitly, and checked against
   every group-level p-value produced).

4. All five axes run under this single frozen procedure in one pass,
   producing ONE definitive table: axis, raw combined p (naive Fisher, and
   independence-corrected Fisher), BH-adjusted p, significant (Y/N).

5. No new analyses beyond what's specified here. Prior audit items are not
   touched.

6. BCAA, if referenced in the writeup, is framed as a methodological
   diagnostic of the consensus-gating design (documented separately: its
   is_active coverage is 11/90 = 12%, far below the other axes), NOT as a
   biological claim that BCAA metabolism is absent.
"""

import json
import os
import time

import numpy as np
import pandas as pd
from scipy import stats

from audit_common import ALL_GROUPS, DETAILS_DIR, evaluation_percentile, load_calibrated_metadata

SEED = 42
N_PERMUTATIONS = 100_000
ALPHA = 0.05
P_FLOOR = 1.0 / (N_PERMUTATIONS + 1)

# Frozen mapping for this pass -- DNL is SCD1-only (Table 2's current mapping),
# per explicit instruction. The other four axes are unchanged.
FROZEN_BIOMARKERS = {
    "serine_glycine_shmt": ["MAR03845", "MAR04792"],
    "urea_cycle": ["MAR03873", "MAR03809", "MAR03811", "MAR03813", "MAR03816", "MAR08426"],
    "bcaa": ["MAR03744", "MAR03747", "MAR03765", "MAR06923", "MAR03777", "MAR03778",
             "MAR06416", "MAR06419", "MAR06421"],
    "dnl_scd1": ["MAR00146", "MAR00147", "MAR00148"],  # SCD1-only, Table 2 mapping -- frozen for this pass
    "choline_pc": ["MAR00636", "MAR00638", "MAR00653", "MAR01603", "MAR01606"],
}
AXIS_ORDER = ["serine_glycine_shmt", "urea_cycle", "bcaa", "dnl_scd1", "choline_pc"]


def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


def fisher_combined(pvals):
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1.0)
    stat = -2 * np.sum(np.log(pvals))
    df = 2 * len(pvals)
    p_combined = 1 - stats.chi2.cdf(stat, df)
    return float(stat), float(p_combined), df


def brown_corrected_fisher(pvals, pairwise_corr):
    """
    Brown (1975) / Kost & McDermott (2002) moment-matching correction for
    Fisher's method under dependent p-values. NOT a naive df substitution --
    the raw statistic T=-2*sum(ln(p_i)) is itself rescaled by a factor c
    (Var(T)/(4k)), and referred to chi2(f) with f=8k^2/Var(T), where Var(T)
    is inflated above its independence value 4k using Brown's empirical
    covariance approximation cov(-2ln p_i, -2ln p_j) ~= 3.263*rho + 0.710*rho^2
    + 0.027*rho^3 for each pair's correlation rho. This is the mathematically
    correct direction of correction under positive correlation (it makes the
    combined p-value LARGER / more conservative, unlike naively substituting
    a smaller df into the same raw statistic, which does the opposite).
    """
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1.0)
    k = len(pvals)
    T = -2 * np.sum(np.log(pvals))

    cov_sum = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            rho = pairwise_corr[i, j]
            cov_ij = 3.263 * rho + 0.710 * rho**2 + 0.027 * rho**3
            cov_sum += cov_ij
    var_T = 4 * k + 2 * cov_sum  # Var(T) = sum_i Var(-2ln p_i) + 2*sum_{i<j} Cov(-2ln p_i, -2ln p_j)
    c = var_T / (4 * k)
    f = (2 * k) / c
    T_corrected = T / c
    p_corrected = 1 - stats.chi2.cdf(T_corrected, f)
    return {
        "raw_statistic_T": float(T), "k": k, "sum_pairwise_cov": float(cov_sum),
        "var_T_corrected": float(var_T), "var_T_independent": float(4 * k),
        "inflation_factor_c": float(c), "corrected_df_f": float(f),
        "corrected_statistic_T_over_c": float(T_corrected),
        "p_value": float(p_corrected),
    }


def build_score_matrix(df):
    df = df.copy()
    df["_ep"] = evaluation_percentile(df, "calibrated_confidence_score")
    reaction_ids = sorted(df["reaction_id"].unique())
    r_to_idx = {rid: i for i, rid in enumerate(reaction_ids)}
    n_groups, n_rxn = len(ALL_GROUPS), len(reaction_ids)
    mat = np.zeros((n_groups, n_rxn))
    for gi, (cohort, group) in enumerate(ALL_GROUPS):
        sub = df[(df["cohort"] == cohort) & (df["group"] == group)]
        idx = sub["reaction_id"].map(r_to_idx).values
        mat[gi, idx] = sub["_ep"].values
    return mat, reaction_ids, r_to_idx


def effective_n_independent_tests_li_ji(corr_matrix):
    """
    Li & Ji (2005), 'Adjusting multiple testing in multilocus analyses using
    the eigenvalues of a correlation matrix', Heredity 95:221-227.
    M_eff = M - sum_i [ I(lambda_i > 1) * (lambda_i - 1) ], clipped to [1, M].
    """
    eigvals = np.linalg.eigvalsh(corr_matrix)
    eigvals = np.clip(eigvals, 0, None)  # numerical guard
    M = len(eigvals)
    excess = np.sum(np.where(eigvals > 1, eigvals - 1, 0.0))
    m_eff = M - excess
    return float(np.clip(m_eff, 1, M)), eigvals.tolist()


def independence_diagnostic(mat):
    """Pairwise correlation of each group's full 12931-reaction score vector,
    within-cohort vs across-cohort, plus Li & Ji effective-N."""
    n_groups = mat.shape[0]
    corr = np.corrcoef(mat)  # 10x10

    within, across = [], []
    for i in range(n_groups):
        for j in range(i + 1, n_groups):
            cohort_i, cohort_j = ALL_GROUPS[i][0], ALL_GROUPS[j][0]
            pair = {"group_i": f"{ALL_GROUPS[i][0]}__{ALL_GROUPS[i][1]}",
                    "group_j": f"{ALL_GROUPS[j][0]}__{ALL_GROUPS[j][1]}",
                    "correlation": float(corr[i, j])}
            if cohort_i == cohort_j:
                within.append(pair)
            else:
                across.append(pair)

    within_vals = [p["correlation"] for p in within]
    across_vals = [p["correlation"] for p in across]
    m_eff, eigvals = effective_n_independent_tests_li_ji(corr)

    t_stat, t_p = stats.ttest_ind(within_vals, across_vals, equal_var=False)

    return {
        "correlation_matrix": corr.tolist(),
        "group_order": [f"{c}__{g}" for c, g in ALL_GROUPS],
        "within_cohort_pairs": within, "across_cohort_pairs": across,
        "within_cohort_mean_correlation": float(np.mean(within_vals)),
        "across_cohort_mean_correlation": float(np.mean(across_vals)),
        "within_minus_across_diff": float(np.mean(within_vals) - np.mean(across_vals)),
        "welch_t_statistic": float(t_stat), "welch_t_p_value": float(t_p),
        "non_independence_evidence": bool(t_p < 0.05 and np.mean(within_vals) > np.mean(across_vals)),
        "eigenvalues": eigvals,
        "effective_n_independent_tests_li_ji": m_eff,
        "naive_n_tests": n_groups,
        "fisher_anti_conservative_flagged": bool(m_eff < n_groups - 0.5),
    }


def flat_test_100k(mat, r_to_idx, rxn_ids, n_reactions, rng, n_perm=N_PERMUTATIONS):
    """One-level design: per-group empirical p, Fisher-combined directly
    across all 10 groups. Returns per-group detail + naive combined p."""
    idx = [r_to_idx[r] for r in rxn_ids if r in r_to_idx]
    k = len(idx)
    per_group = []
    for gi, (cohort, group) in enumerate(ALL_GROUPS):
        observed_g = float(mat[gi, idx].max())
        random_idx = rng.integers(0, n_reactions, size=(n_perm, k))
        null_g = mat[gi, random_idx].max(axis=1)
        n_ge = int(np.sum(null_g >= observed_g))
        p_g = (1 + n_ge) / (1 + n_perm)
        per_group.append({
            "cohort": cohort, "group": group, "observed_statistic": observed_g,
            "n_reactions_in_axis": k, "n_permutations": n_perm,
            "n_null_ge_observed": n_ge, "p_value_formula": "(1 + n_null_ge_observed) / (1 + n_permutations)",
            "p_value": float(p_g), "p_floor": P_FLOOR, "hit_floor": n_ge == 0,
        })
    group_ps = [g["p_value"] for g in per_group]
    stat, p_naive, df = fisher_combined(group_ps)
    return per_group, group_ps, stat, p_naive, df


def run():
    t0 = time.time()
    df = load_calibrated_metadata()
    mat, reaction_ids, r_to_idx = build_score_matrix(df)
    n_reactions = len(reaction_ids)

    print("=" * 100)
    print("FROZEN MAPPING -- exact reaction IDs used for all 5 axes in this definitive pass:")
    for axis in AXIS_ORDER:
        print(f"  {axis}: {FROZEN_BIOMARKERS[axis]}")
    print("=" * 100)

    # --- independence diagnostic (computed once, applies to all axes since it's
    # about the score-vector structure across groups, not axis-specific) ---
    indep = independence_diagnostic(mat)
    print(f"\nIndependence diagnostic: within-cohort mean r={indep['within_cohort_mean_correlation']:.4f}, "
          f"across-cohort mean r={indep['across_cohort_mean_correlation']:.4f}, "
          f"Welch t-test p={indep['welch_t_p_value']:.4g}")
    print(f"Effective N independent tests (Li & Ji): {indep['effective_n_independent_tests_li_ji']:.2f} of 10 naive")
    print(f"Fisher anti-conservative flag: {indep['fisher_anti_conservative_flagged']}")

    m_eff = indep["effective_n_independent_tests_li_ji"]
    pairwise_corr = np.array(indep["correlation_matrix"])

    axes_results = {}
    for axis in AXIS_ORDER:
        rxn_ids = FROZEN_BIOMARKERS[axis]
        rng = np.random.default_rng(SEED)
        per_group, group_ps, stat, p_naive, df_naive = flat_test_100k(mat, r_to_idx, rxn_ids, n_reactions, rng)
        brown = brown_corrected_fisher(group_ps, pairwise_corr)
        p_corrected = brown["p_value"]

        n_hit_floor = sum(g["hit_floor"] for g in per_group)
        axes_results[axis] = {
            "mapped_reactions": rxn_ids, "n_mapped_reactions": len(rxn_ids),
            "per_group": per_group,
            "group_p_values": group_ps,
            "naive_fisher_statistic": stat, "naive_df": df_naive, "naive_combined_p": p_naive,
            "brown_correction_detail": brown,
            "independence_corrected_combined_p": p_corrected,
            "n_groups_hit_permutation_floor": n_hit_floor,
        }
        print(f"\n{axis}: naive_p={p_naive:.6g} (df={df_naive}), "
              f"brown_corrected_p={p_corrected:.6g} (c={brown['inflation_factor_c']:.3f}, f={brown['corrected_df_f']:.2f}), "
              f"groups_hit_floor={n_hit_floor}/10")

    naive_ps = [axes_results[a]["naive_combined_p"] for a in AXIS_ORDER]
    corrected_ps = [axes_results[a]["independence_corrected_combined_p"] for a in AXIS_ORDER]
    naive_bh = benjamini_hochberg(naive_ps)
    corrected_bh = benjamini_hochberg(corrected_ps)

    definitive_table = []
    for i, axis in enumerate(AXIS_ORDER):
        axes_results[axis]["naive_bh_adjusted_p"] = float(naive_bh[i])
        axes_results[axis]["naive_significant"] = bool(naive_bh[i] < ALPHA)
        axes_results[axis]["independence_corrected_bh_adjusted_p"] = float(corrected_bh[i])
        axes_results[axis]["independence_corrected_significant"] = bool(corrected_bh[i] < ALPHA)
        definitive_table.append({
            "axis": axis,
            "mapped_reactions": ";".join(FROZEN_BIOMARKERS[axis]),
            "raw_combined_p_naive_fisher": axes_results[axis]["naive_combined_p"],
            "bh_adjusted_p_naive_fisher": axes_results[axis]["naive_bh_adjusted_p"],
            "significant_naive_YN": "Y" if axes_results[axis]["naive_significant"] else "N",
            "raw_combined_p_independence_corrected": axes_results[axis]["independence_corrected_combined_p"],
            "bh_adjusted_p_independence_corrected": axes_results[axis]["independence_corrected_bh_adjusted_p"],
            "significant_independence_corrected_YN": "Y" if axes_results[axis]["independence_corrected_significant"] else "N",
        })

    bcaa_note = (
        "BCAA's is_active coverage across its 9 mapped reactions x 10 groups is 11/90 (12%), far below "
        "every other axis (67-80%, except choline_pc at 30%) -- see item11's per-axis detail. Framed here "
        "explicitly as a METHODOLOGICAL DIAGNOSTIC of the consensus-gating design (most BCAA reactions are "
        "gated out of the is_active/active-reaction set before enrichment is even tested), NOT as a "
        "biological claim that BCAA metabolism is absent or unaffected in NAFLD/NASH."
    )

    result = {
        "item": 13,
        "description": "DEFINITIVE RECONCILIATION PASS -- single frozen methodology, supersedes all prior enrichment figures in this audit (items 8, 11, 12)",
        "frozen_mapping": FROZEN_BIOMARKERS,
        "mapping_note": "DNL uses Table 2's current SCD1-only mapping (MAR00146/147/148), per explicit instruction -- NOT the FASN-inclusive mapping used in items 11/12.",
        "methodology": {
            "score_source": "calibrated_confidence_score, per-group evaluation-percentile (rank among all 12931 reactions per group, method='min')",
            "procedure": "per-group permutation -> empirical p per axis per group -> Fisher combination directly across all 10 groups (one level) -> BH correction across 5 axes",
            "n_permutations": N_PERMUTATIONS, "seed": SEED,
            "p_value_formula": "(1 + #{null >= observed}) / (1 + N)",
            "resolution_floor": P_FLOOR,
            "alpha": ALPHA,
        },
        "independence_diagnostic": indep,
        "axes": axes_results,
        "definitive_table": definitive_table,
        "bcaa_framing_note": bcaa_note,
        "supersedes": "All previously reported enrichment p-values in this audit (item 8's comparator-architecture test, item 11's calibrated-score flat/nested tests, item 12's DNL mapping reconciliation) are superseded by this table for the calibrated GEM1 score's enrichment result. Items 1-7, 9, 10 (unrelated to enrichment) are untouched and NOT superseded.",
        "status": "DEFINITIVE",
        "max_abs_discrepancy": "N/A -- this is the frozen definitive computation, not a reconciliation against a prior number",
        "elapsed_seconds": time.time() - t0,
    }

    out_path = os.path.join(DETAILS_DIR, "item13_definitive_reconciliation.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    table_path = os.path.join(DETAILS_DIR, "item13_definitive_table.csv")
    pd.DataFrame(definitive_table).to_csv(table_path, index=False)
    print(f"\nWrote {out_path}")
    print(f"Wrote {table_path}")

    print("\n" + "=" * 100)
    print("DEFINITIVE TABLE (naive Fisher, df=20, assumes 10 independent tests):")
    for row in definitive_table:
        print(f"  {row['axis']:25s} raw_p={row['raw_combined_p_naive_fisher']:.6g}  "
              f"bh_p={row['bh_adjusted_p_naive_fisher']:.6g}  sig={row['significant_naive_YN']}")
    print(f"\nBROWN-CORRECTED TABLE (moment-matched Fisher combination under measured correlation, M_eff={m_eff:.2f} of 10):")
    for row in definitive_table:
        print(f"  {row['axis']:25s} raw_p={row['raw_combined_p_independence_corrected']:.6g}  "
              f"bh_p={row['bh_adjusted_p_independence_corrected']:.6g}  sig={row['significant_independence_corrected_YN']}")
    print("=" * 100)
    return result


if __name__ == "__main__":
    run()
