"""
GEM1 - ad hoc: cross-group numerical synthesis of the 10 ACCEPTED paired-
counterfactual production groups (data/paired_perturbation_production/).

READ-ONLY with respect to production: only reads manifest.json/checkpoint.json/
reaction_scores.csv from each of the 10 groups, plus the pilot's
reaction_scores.csv and the frozen BIOMARKERS_FULL reaction set imported
directly from 09_calibration.py (not redefined here). Writes new tables only
to data/paired_perturbation_synthesis/ -- never touches
data/paired_perturbation_production/ or data/paired_perturbation_pilot/.

No weight recalibration, no consensus/flux combination, no TAU change, no
rerun of any simulation. Numeric/statistical synthesis only -- no biological
interpretation is asserted anywhere in this script or its outputs.
"""
import glob
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, mannwhitneyu

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROD_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_production")
PILOT_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_pilot")
OUT_DIR = os.path.join(BASE_DIR, "data", "paired_perturbation_synthesis")
os.makedirs(OUT_DIR, exist_ok=True)

TAU = 1e-6

GROUPS = [
    "GSE126848__NASH", "GSE126848__healthy", "GSE126848__obese_no_NAFLD",
    "GSE126848__steatosis", "GSE135251__NASH", "GSE135251__healthy",
    "GSE135251__steatosis", "GSE89632__NASH", "GSE89632__healthy",
    "GSE89632__steatosis",
]

# Import the frozen biomarker reaction set directly from the existing
# calibration script rather than re-typing it (single source of truth).
_spec = importlib.util.spec_from_file_location(
    "step09_calibration", os.path.join(BASE_DIR, "scripts", "09_calibration.py"))
_calib = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_calib)
BIOMARKERS_FULL = _calib.BIOMARKERS_FULL

PILOT_GROUPS = {
    "GSE126848__obese_no_NAFLD": ("GSE126848__obese_no_NAFLD", [101, 202]),
    "GSE135251__healthy": ("GSE135251__healthy", [101, 202, 303, 404]),
}


def load_group(g):
    gdir = os.path.join(PROD_DIR, g)
    with open(os.path.join(gdir, "manifest.json")) as f:
        manifest = json.load(f)
    assert manifest["status"] == "ACCEPTED", f"{g} is not ACCEPTED -- refusing to synthesize"
    scores = pd.read_csv(os.path.join(gdir, "reaction_scores.csv"))
    return manifest, scores


def main():
    data = {g: load_group(g) for g in GROUPS}
    print(f"Loaded {len(data)} ACCEPTED groups.")

    # ------------------------------------------------------------------
    # 1. Cross-group distribution of raw PR_objective
    # ------------------------------------------------------------------
    rows = []
    pooled_pr = []
    for g, (manifest, scores) in data.items():
        pr = scores["PR_objective"]
        pooled_pr.append(pr)
        rows.append({
            "group": g, "n_reactions": len(pr),
            "mean": pr.mean(), "std": pr.std(), "min": pr.min(),
            "p10": pr.quantile(0.10), "p25": pr.quantile(0.25),
            "median": pr.median(), "p75": pr.quantile(0.75),
            "p90": pr.quantile(0.90), "max": pr.max(),
        })
    t1 = pd.DataFrame(rows).sort_values("group")
    t1.to_csv(os.path.join(OUT_DIR, "01_cross_group_PR_objective_distribution.csv"), index=False)

    pooled = pd.concat(pooled_pr, ignore_index=True)
    with open(os.path.join(OUT_DIR, "01b_pooled_PR_objective_distribution.json"), "w") as f:
        json.dump({
            "n_reaction_rows_pooled_across_10_groups": int(len(pooled)),
            "mean": float(pooled.mean()), "std": float(pooled.std()),
            "min": float(pooled.min()), "p10": float(pooled.quantile(0.10)),
            "p25": float(pooled.quantile(0.25)), "median": float(pooled.median()),
            "p75": float(pooled.quantile(0.75)), "p90": float(pooled.quantile(0.90)),
            "max": float(pooled.max()),
            "note": "Pooled = every (reaction, group) row concatenated, unweighted "
                    "by group; NOT deduplicated by reaction_id.",
        }, f, indent=2)
    print("Wrote 01_cross_group_PR_objective_distribution.csv, 01b_pooled_*.json")

    # ------------------------------------------------------------------
    # 2. Per-group candidate count, above-TAU count, occ_obj/occ_flip sparsity
    # ------------------------------------------------------------------
    rows = []
    for g, (manifest, scores) in data.items():
        n = len(scores)
        n_above_tau_obj = int((scores["occ_obj_k"] > 0).sum())
        n_above_tau_flip = int((scores["occ_flip_k"] > 0).sum())
        rows.append({
            "group": g,
            "n_dropout_candidates_manifest": manifest["n_dropout_candidates"],
            "n_reactions_scored": n,
            "n_occ_obj_nonzero": n_above_tau_obj,
            "occ_obj_sparsity_frac_zero": 1 - n_above_tau_obj / n,
            "n_occ_flip_nonzero": n_above_tau_flip,
            "occ_flip_sparsity_frac_zero": 1 - n_above_tau_flip / n,
        })
    t2 = pd.DataFrame(rows).sort_values("group")
    t2.to_csv(os.path.join(OUT_DIR, "02_group_candidate_and_sparsity_summary.csv"), index=False)
    print("Wrote 02_group_candidate_and_sparsity_summary.csv")

    # ------------------------------------------------------------------
    # 3 & 4. Cross-group recurrence, preserving NE vs evaluated-zero
    # ------------------------------------------------------------------
    all_reactions = sorted(set().union(*[set(s["reaction_id"]) for _, s in data.values()]))
    print(f"Full-network union across all 10 groups: {len(all_reactions)} distinct reaction_ids")

    per_group_status = {}  # g -> {reaction_id: 'NE'|'zero'|'nonzero'}
    per_group_pr = {}      # g -> {reaction_id: PR_objective}
    for g, (manifest, scores) in data.items():
        present = set(scores["reaction_id"])
        status = {}
        pr_map = dict(zip(scores["reaction_id"], scores["PR_objective"]))
        k_map = dict(zip(scores["reaction_id"], scores["occ_obj_k"]))
        for rid in all_reactions:
            if rid not in present:
                status[rid] = "NE"
            elif k_map[rid] > 0:
                status[rid] = "nonzero"
            else:
                status[rid] = "zero"
        per_group_status[g] = status
        per_group_pr[g] = pr_map

    matrix_rows = []
    for rid in all_reactions:
        row = {"reaction_id": rid}
        n_ne = n_zero = n_nonzero = 0
        for g in GROUPS:
            s = per_group_status[g][rid]
            row[g] = s
            if s == "NE":
                n_ne += 1
            elif s == "zero":
                n_zero += 1
            else:
                n_nonzero += 1
        row["n_groups_NE"] = n_ne
        row["n_groups_evaluated_zero"] = n_zero
        row["n_groups_evaluated_nonzero"] = n_nonzero
        row["recurrence_count_nonzero"] = n_nonzero  # reactions with delta_objective>TAU in >=1 background, per group
        matrix_rows.append(row)
    t34 = pd.DataFrame(matrix_rows).sort_values("recurrence_count_nonzero", ascending=False)
    t34.to_csv(os.path.join(OUT_DIR, "03_04_reaction_by_group_status_matrix.csv"), index=False)

    hist = t34["recurrence_count_nonzero"].value_counts().sort_index()
    hist_df = hist.rename_axis("n_groups_with_nonzero_occ_obj").reset_index(name="n_reactions")
    hist_df.to_csv(os.path.join(OUT_DIR, "03b_recurrence_histogram.csv"), index=False)
    print("Wrote 03_04_reaction_by_group_status_matrix.csv, 03b_recurrence_histogram.csv")

    # ------------------------------------------------------------------
    # 5. 25 calibration/biomarker reactions vs full-network distribution
    # ------------------------------------------------------------------
    biomarker_rxns = sorted({rid for lst in BIOMARKERS_FULL.values() for rid in lst})
    assert len(biomarker_rxns) == 25, f"expected 25 biomarker reactions, got {len(biomarker_rxns)}"
    rxn_to_biomarker = {}
    for name, lst in BIOMARKERS_FULL.items():
        for rid in lst:
            rxn_to_biomarker.setdefault(rid, []).append(name)

    detail_rows = []
    for g, (manifest, scores) in data.items():
        s = scores.copy()
        s["pct_rank_PR_objective"] = s["PR_objective"].rank(pct=True) * 100
        s_idx = s.set_index("reaction_id")
        for rid in biomarker_rxns:
            if rid in s_idx.index:
                r = s_idx.loc[rid]
                detail_rows.append({
                    "reaction_id": rid, "biomarker_set": "+".join(rxn_to_biomarker[rid]),
                    "group": g, "status": "evaluated",
                    "PR_objective": r["PR_objective"],
                    "pct_rank_PR_objective": r["pct_rank_PR_objective"],
                    "occ_obj_k": r["occ_obj_k"], "occ_obj_p": r["occ_obj_p"],
                    "n_eligible": r["n_eligible"],
                })
            else:
                detail_rows.append({
                    "reaction_id": rid, "biomarker_set": "+".join(rxn_to_biomarker[rid]),
                    "group": g, "status": "NE",
                    "PR_objective": None, "pct_rank_PR_objective": None,
                    "occ_obj_k": None, "occ_obj_p": None, "n_eligible": None,
                })
    t5 = pd.DataFrame(detail_rows)
    t5.to_csv(os.path.join(OUT_DIR, "05_biomarker_reactions_detail.csv"), index=False)

    summary_rows = []
    for g, (manifest, scores) in data.items():
        full_pr = scores["PR_objective"]
        sub = t5[(t5["group"] == g) & (t5["status"] == "evaluated")]
        bio_pr = sub["PR_objective"].astype(float)
        row = {
            "group": g,
            "n_biomarker_reactions_evaluated": len(bio_pr),
            "n_biomarker_reactions_NE": 25 - len(bio_pr),
            "biomarker_mean_PR_objective": bio_pr.mean() if len(bio_pr) else None,
            "biomarker_median_PR_objective": bio_pr.median() if len(bio_pr) else None,
            "biomarker_mean_pct_rank": sub["pct_rank_PR_objective"].astype(float).mean() if len(bio_pr) else None,
            "full_network_mean_PR_objective": full_pr.mean(),
            "full_network_median_PR_objective": full_pr.median(),
            "full_network_n": len(full_pr),
        }
        if len(bio_pr) >= 2:
            stat, p = mannwhitneyu(bio_pr, full_pr, alternative="two-sided")
            row["mannwhitney_U"] = stat
            row["mannwhitney_p_two_sided"] = p
        else:
            row["mannwhitney_U"] = None
            row["mannwhitney_p_two_sided"] = None
        summary_rows.append(row)
    t5b = pd.DataFrame(summary_rows).sort_values("group")
    t5b.to_csv(os.path.join(OUT_DIR, "05b_biomarker_vs_full_network_summary.csv"), index=False)
    print("Wrote 05_biomarker_reactions_detail.csv, 05b_biomarker_vs_full_network_summary.csv")

    # ------------------------------------------------------------------
    # 6. Monte Carlo uncertainty quantification (already-frozen outputs)
    # ------------------------------------------------------------------
    rows = []
    for g, (manifest, scores) in data.items():
        jw_obj = scores["occ_obj_jeffreys_hi"] - scores["occ_obj_jeffreys_lo"]
        jw_flip = scores["occ_flip_jeffreys_hi"] - scores["occ_flip_jeffreys_lo"]
        rows.append({
            "group": g,
            "PR_objective_se_mean": scores["PR_objective_se"].mean(),
            "PR_objective_se_median": scores["PR_objective_se"].median(),
            "PR_objective_se_max": scores["PR_objective_se"].max(),
            "occ_obj_jeffreys_width_mean": jw_obj.mean(),
            "occ_obj_jeffreys_width_median": jw_obj.median(),
            "occ_flip_jeffreys_width_mean": jw_flip.mean(),
            "occ_flip_jeffreys_width_median": jw_flip.median(),
            "monte_carlo_uncertainty_scope": manifest["monte_carlo_uncertainty_scope"],
        })
    t6 = pd.DataFrame(rows).sort_values("group")
    t6.to_csv(os.path.join(OUT_DIR, "06_monte_carlo_uncertainty_summary.csv"), index=False)
    print("Wrote 06_monte_carlo_uncertainty_summary.csv")

    # ------------------------------------------------------------------
    # 7. Production vs pilot comparison, where pilot data exists
    # ------------------------------------------------------------------
    corr_rows = []
    merged_all = []
    for g, (pilot_prefix, seeds) in PILOT_GROUPS.items():
        manifest, prod_scores = data[g]
        pilot_frames = []
        for seed in seeds:
            p = os.path.join(PILOT_DIR, f"{pilot_prefix}__seed{seed}__reaction_scores.csv")
            df = pd.read_csv(p)
            df["seed"] = seed
            pilot_frames.append(df)
        pilot_pooled = (pd.concat(pilot_frames, ignore_index=True)
                         .groupby("reaction_id")["PR_objective"].mean()
                         .reset_index().rename(columns={"PR_objective": "PR_objective_pilot"}))
        merged = prod_scores.merge(pilot_pooled, on="reaction_id", how="inner")
        merged["group"] = g
        merged_all.append(merged[["group", "reaction_id", "PR_objective", "PR_objective_pilot"]])
        r_p, p_p = pearsonr(merged["PR_objective"], merged["PR_objective_pilot"])
        r_s, p_s = spearmanr(merged["PR_objective"], merged["PR_objective_pilot"])
        corr_rows.append({
            "group": g, "pilot_seeds": seeds, "pilot_M_effective": 120 * len(seeds),
            "n_overlapping_reactions": len(merged),
            "pearson_r": r_p, "pearson_p": p_p,
            "spearman_r": r_s, "spearman_p": p_s,
        })
    t7a = pd.DataFrame(corr_rows)
    t7a.to_csv(os.path.join(OUT_DIR, "07a_pilot_vs_production_correlation.csv"), index=False)
    t7b = pd.concat(merged_all, ignore_index=True)
    t7b.to_csv(os.path.join(OUT_DIR, "07b_pilot_vs_production_merged_values.csv"), index=False)
    print("Wrote 07a_pilot_vs_production_correlation.csv, 07b_pilot_vs_production_merged_values.csv")
    print(f"\nNote: pilot data only exists for {list(PILOT_GROUPS.keys())} "
          f"(2 of 10 groups) -- comparison is 'where applicable' as instructed, "
          f"the other 8 groups have no pilot counterpart.")

    print("\nAll synthesis tables written to:", OUT_DIR)
    print("Raw production outputs under data/paired_perturbation_production/ were only read, never modified.")


if __name__ == "__main__":
    main()
