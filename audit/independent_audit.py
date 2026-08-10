"""
GEM1 -- Independent Arithmetic and Provenance Audit

Orchestrator / aggregator. This script itself does not import any of GEM1's
own scoring/calibration modules (scripts/04-11); it reads only raw source
outputs and results already produced by this audit's own item scripts
(item1_imat_fva_baseline.py, item2_iqr_reconciliation.py,
item3_flux_sampling_perturbation.py, item4_growth_support_lp.py,
item5_loco_lobo_refit.py, item6_troppo_reextract.py).

Usage:
  python independent_audit.py --run-fast     # runs items 2, 4, 5 (seconds, synchronous)
  python independent_audit.py --finalize     # aggregates whatever is in audit/details/
                                              # into audit_report.csv / audit_summary.json / README.md

Items 1, 3, and 6 are long-running (~30/~75/~70 min respectively) and must be
launched as detached Start-Process jobs in the correct conda env -- see
README.md "How to run" for the exact commands. Re-run --finalize after they
complete to fold their results in.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from audit_common import AUDIT_DIR, DETAILS_DIR

ITEM_FILES = {
    1: "item1_imat_fva_baseline.json",
    2: "item2_iqr_reconciliation.json",
    3: "item3_flux_sampling_perturbation.json",
    4: "item4_growth_support_lp.json",
    5: "item5_loco_lobo_refit.json",
    6: "item6_troppo_reextract.json",
    7: "item7_flux_uncertainty_broader_population.json",
    8: "item8_enrichment_permutation_test.json",
    9: "item9_sd_sensitivity.json",
    10: "item10_bcaa_severity_means.json",
    11: "item11_calibrated_enrichment_audit.json",
    12: "item12_dnl_mapping_reconciliation.json",
    13: "item13_definitive_reconciliation.json",
    14: "item14_unconstrained_weight_refit.json",
    15: "item15_null_sampling_verification.json",
    16: "item16_growth_support_confound_check.json",
    17: "item17_per_axis_weight_preferences.json",
    18: "item18_uncontaminated_axes_reproducibility.json",
    19: "item19_shmt_cohort_nested_final.json",
}

ITEM_TITLES = {
    1: "Single-Algorithm (iMAT+FVA) external baseline",
    2: "Flux Uncertainty-Only / Uncalibrated Equal Weighting IQR reconciliation (Table 3)",
    3: "Flux sampling + perturbation robustness re-execution",
    4: "Growth-support LP rebuild (519 reactions, 90% biomass floor)",
    5: "LOCO/LOBO weight optimization refit",
    6: "Active-reaction counts via independent troppo reconstruction",
    7: "Flux Uncertainty-Only IQR reconciliation (broader-population question)",
    8: "Enrichment recovery counts and p-values (Flux Uncertainty-Only / Equal Weighting)",
    9: "Component-level SD sensitivity (Section 3.4)",
    10: "Per-group BCAA severity means (Section 3.7)",
    11: "Calibrated-score enrichment audit (Section 3.8 discrepancy investigation)",
    12: "DNL enrichment mapping reconciliation (2/5 vs 3/5 significant axes)",
    13: "DEFINITIVE reconciliation pass (frozen methodology, 100k draws, supersedes all prior enrichment figures)",
    14: "Unconstrained (0-1) LOCO/LOBO weight re-fit -- was the 0.15 floor binding?",
    15: "Null-sampling verification (reaction-count-matched draws, code citation)",
    16: "Growth-support confound check on item14's consensus-dominant finding",
    17: "Per-axis unconstrained weight preferences (all 5 biomarker axes individually)",
    18: "Kendall's W and SD-ordering benchmark, excluding growth-support-contaminated axes",
    19: "Final SHMT p=2.1e-5 reconciliation attempt (two-level nested, cohort-then-overall)",
}


def run_fast_items():
    import item2_iqr_reconciliation
    import item4_growth_support_lp
    import item5_loco_lobo_refit

    print("=== Item 4: growth-support LP ===")
    item4_growth_support_lp.run()
    print("\n=== Item 5: LOCO/LOBO refit ===")
    item5_loco_lobo_refit.run()
    print("\n=== Item 2: IQR reconciliation ===")
    item2_iqr_reconciliation.run()


def load_item_result(item_num):
    path = os.path.join(DETAILS_DIR, ITEM_FILES[item_num])
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def finalize():
    import csv

    rows = []
    summary = {"generated_at": datetime.now(timezone.utc).isoformat(), "items": {}}

    for item_num in range(1, 7):
        result = load_item_result(item_num)
        title = ITEM_TITLES[item_num]
        if result is None:
            rows.append({
                "item": item_num, "title": title, "status": "NOT_RUN",
                "max_abs_discrepancy": "", "notes": "Detail file not found -- long-running item not yet launched or still in progress.",
            })
            summary["items"][item_num] = {"status": "NOT_RUN", "title": title}
            continue

        status = result.get("status", "UNKNOWN")
        max_abs = result.get("max_abs_discrepancy", "")
        notes = result.get("note") or result.get("status_meaning") or result.get("description", "")
        rows.append({
            "item": item_num, "title": title, "status": status,
            "max_abs_discrepancy": max_abs, "notes": notes,
        })
        summary["items"][item_num] = result

    report_path = os.path.join(AUDIT_DIR, "audit_report.csv")
    with open(report_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["item", "title", "status", "max_abs_discrepancy", "notes"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {report_path}")

    summary_path = os.path.join(AUDIT_DIR, "audit_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {summary_path}")

    write_readme(rows, summary)
    return rows, summary


def finalize_append(new_item_nums):
    """
    Append-only finalize: preserves existing audit_report.csv / audit_summary.json
    rows for items not in new_item_nums, adds/updates rows only for the given
    item numbers, and appends a new dated section to README.md rather than
    regenerating the whole file. Used for a second (or later) audit run that
    continues past items already reported by a prior --finalize run.
    """
    import csv

    report_path = os.path.join(AUDIT_DIR, "audit_report.csv")
    summary_path = os.path.join(AUDIT_DIR, "audit_summary.json")

    existing_rows = []
    if os.path.exists(report_path):
        with open(report_path, newline="") as f:
            existing_rows = list(csv.DictReader(f))
    existing_by_item = {int(r["item"]): r for r in existing_rows}

    existing_summary = {"items": {}, "run_history": []}
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            existing_summary = json.load(f)
        existing_summary.setdefault("run_history", [])

    new_rows = []
    for item_num in new_item_nums:
        result = load_item_result(item_num)
        title = ITEM_TITLES[item_num]
        if result is None:
            row = {"item": item_num, "title": title, "status": "NOT_RUN",
                   "max_abs_discrepancy": "", "notes": "Detail file not found."}
            existing_summary["items"][str(item_num)] = {"status": "NOT_RUN", "title": title}
        else:
            status = result.get("status", "UNKNOWN")
            max_abs = result.get("max_abs_discrepancy", "")
            notes = result.get("note") or result.get("status_meaning") or result.get("description", "")
            row = {"item": item_num, "title": title, "status": status,
                   "max_abs_discrepancy": max_abs, "notes": notes}
            existing_summary["items"][str(item_num)] = result
        existing_by_item[item_num] = row
        new_rows.append(row)

    all_rows = [existing_by_item[k] for k in sorted(existing_by_item.keys())]
    with open(report_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["item", "title", "status", "max_abs_discrepancy", "notes"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Appended to {report_path} (items {new_item_nums})")

    run_timestamp = datetime.now(timezone.utc).isoformat()
    existing_summary["run_history"].append({"timestamp": run_timestamp, "items_added_or_updated": new_item_nums})
    existing_summary["generated_at"] = run_timestamp
    with open(summary_path, "w") as f:
        json.dump(existing_summary, f, indent=2)
    print(f"Appended to {summary_path}")

    append_readme_section(new_item_nums, existing_summary, run_timestamp)
    return new_rows, existing_summary


def append_readme_section(new_item_nums, summary, run_timestamp):
    readme_path = os.path.join(AUDIT_DIR, "README.md")
    lines = []
    lines.append(f"\n---\n\n## Run 2 -- {run_timestamp}\n")
    lines.append(
        "Continuation of the audit above: items 7-10, covering Table 3's remaining Flux "
        "Uncertainty-Only reconciliation, Section 3.4 (component SD sensitivity), Section 3.7 "
        "(BCAA severity means), and Section 3.8 (enrichment recovery/p-values). Same rules as "
        "Run 1: no GEM1 scoring/calibration module imported, no existing GEM1 output file "
        "modified. Rows below are appended to audit_report.csv / audit_summary.json -- Run 1's "
        "items 1-6 are unchanged.\n"
    )
    lines.append("### Run 2 summary\n")
    lines.append("| # | Item | Status | Max abs. discrepancy |")
    lines.append("|---|------|--------|----------------------|")
    for item_num in new_item_nums:
        r = summary["items"].get(str(item_num)) or summary["items"].get(item_num)
        title = ITEM_TITLES[item_num]
        status = r.get("status", "NOT_RUN") if isinstance(r, dict) else "NOT_RUN"
        max_abs = r.get("max_abs_discrepancy", "") if isinstance(r, dict) else ""
        lines.append(f"| {item_num} | {title} | **{status}** | {max_abs} |")
    lines.append("")

    r7 = summary["items"].get("7") or summary["items"].get(7)
    if isinstance(r7, dict) and "flux_uncertainty_only_iqr_under_validated_convention" in r7:
        lines.append("### 7. Flux Uncertainty-Only IQR reconciliation (broader-population question)\n")
        lines.append(
            f"Before running any new flux sampling, checked whether the SAME population/representation "
            f"convention independently validated by item 9 (three separate SD checks, one against a "
            f"directly-persisted column) also resolves this one. It does: under `{r7['validated_convention']}`, "
            f"the equal-weighting cross-check reproduces **{r7['equal_weighting_cross_check']['independent_value_this_script']:.4f}** "
            f"(request-stated resolved value: {r7['equal_weighting_cross_check']['request_stated_resolved_value']}, "
            f"diff {r7['equal_weighting_cross_check']['diff']:.2e}), and Flux-Uncertainty-Only's IQR under the "
            f"same convention is **{r7['flux_uncertainty_only_iqr_under_validated_convention']:.4f}** "
            f"(reported manuscript value: {r7['reported_manuscript_value']}).\n\n"
            f"**No new flux sampling was run.** Recommendation: {r7['recommendation']}\n\n"
        )

    r8 = summary["items"].get("8") or summary["items"].get(8)
    if isinstance(r8, dict) and "results_by_architecture" in r8:
        lines.append("### 8. Enrichment recovery counts and p-values\n")
        gt = r8.get("ground_truth_search", {})
        lines.append(f"**Ground-truth search finding**: {gt.get('finding', '')}\n\n")
        for arch, r in r8["results_by_architecture"].items():
            lines.append(
                f"- **{arch}**: axes recovered = {r['n_axes_recovered_of_5']}/5 (reported: "
                f"{r['reported_axes_recovered']}/5), Fisher combined p = {r['fisher_combined_p_value']:.4g} "
                f"(reported: {r['reported_p_value']:.4g})\n"
            )
        lines.append("")

    r9 = summary["items"].get("9") or summary["items"].get(9)
    if isinstance(r9, dict) and "closest_cell_per_architecture" in r9:
        lines.append("### 9. Component-level SD sensitivity (Section 3.4)\n")
        for arch, c in r9["closest_cell_per_architecture"].items():
            lines.append(
                f"- **{arch}**: independent SD = {c['value']:.4f} (reported: {r9['reported_sd'][arch]}, "
                f"diff {c['abs_diff']:.2e}) via `{c['cell']}`\n"
            )
        lines.append("")

    r10 = summary["items"].get("10") or summary["items"].get(10)
    if isinstance(r10, dict) and "independent_means" in r10:
        lines.append("### 10. Per-group BCAA severity means (Section 3.7)\n")
        for sev, val in r10["independent_means"].items():
            lines.append(f"- **{sev}**: independent = {val:.4f} (reported: {r10['reported_means'][sev]})\n")
        lines.append("")

    r11 = summary["items"].get("11") or summary["items"].get(11)
    if isinstance(r11, dict) and "axes" in r11:
        lines.append("### 11. Calibrated-score enrichment audit (Section 3.8 discrepancy investigation)\n")
        lines.append(
            f"Manuscript files were NOT modified. Full per-axis, per-group detail: "
            f"`audit/details/item11_axis_table.csv` and `audit/details/item11_calibrated_enrichment_audit.json`.\n\n"
            f"**Ground truth search**: {r11['ground_truth_search']['finding']}\n\n"
            f"**Flat design** (one permutation test per axis, all 10 groups pooled into one statistic): "
            f"axes significant at BH<0.05 = {r11['flat_design_overall']['axes_significant_bh_0.05'] or 'none'}.\n\n"
            f"**Nested design** (per-group p -> Fisher-combine per cohort -> Fisher-combine across cohorts): "
            f"axes significant at BH<0.05 = {r11['nested_design_overall']['axes_significant_bh_0.05'] or 'none'} "
            f"-- matches the previously-reported significant set (SHMT, choline_pc) exactly.\n\n"
            f"**SHMT discrepancy**: {r11['shmt_discrepancy_source_determination']['conclusion']}\n\n"
            f"**PEMT/choline_pc seed-variance check** ({len(r11['pemt_choline_pc_seed_variance_check']['seeds_used'])} seeds): "
            f"{r11['pemt_choline_pc_seed_variance_check']['interpretation']}\n\n"
            f"**Newly significant axes vs. previously reported**: flat={r11['newly_significant_axes']['flat'] or 'none'}, "
            f"nested={r11['newly_significant_axes']['nested'] or 'none'}.\n"
        )
        if r11.get("bcaa_flag"):
            lines.append(f"\n> **[BCAA FLAG]** {r11['bcaa_flag']}\n")

    r12 = summary["items"].get("12") or summary["items"].get(12)
    if isinstance(r12, dict) and "item11_dnl_results" in r12:
        lines.append("### 12. DNL enrichment mapping reconciliation (2/5 vs 3/5 significant axes)\n")
        lines.append(
            f"Manuscript files were NOT modified. Full detail: `audit/details/item12_dnl_mapping_reconciliation.json`.\n\n"
            f"- item11 mapping: `{r12['item11_mapping_used']}` -- flat bh_p={r12['item11_dnl_results']['flat_bh_p_recomputed']:.4g}, "
            f"nested bh_p={r12['item11_dnl_results']['nested_bh_p_recomputed']:.4g} (both non-significant)\n"
            f"- Alternative mapping (reportedly used outside this repo): `{r12['alternative_mapping_used']}` -- "
            f"flat raw_p={r12['alternative_mapping_results']['flat_raw_p']:.4g}, bh_p={r12['alternative_mapping_results']['flat_bh_p_recomputed']:.4g}; "
            f"nested raw_p={r12['alternative_mapping_results']['nested_raw_p']:.4g}, bh_p={r12['alternative_mapping_results']['nested_bh_p_recomputed']:.4g} "
            f"(nested significant)\n"
            f"- GPR check: alternative mapping's 3 reactions confirmed to share an identical GPR "
            f"({r12['gpr_check']['alternative_mapping']['summary']['distinct_gprs']}), independently verified "
            f"from Human-GEM.json, consistent with docs/STEP11_BIOMARKER_MAPPING.md's claim.\n\n"
            f"{r12['conclusion']}\n"
        )

    r13 = summary["items"].get("13") or summary["items"].get(13)
    if isinstance(r13, dict) and "definitive_table" in r13:
        lines.append("### 13. DEFINITIVE reconciliation pass (frozen methodology)\n")
        lines.append(
            f"**Manuscript files were NOT modified.** This table supersedes all previously reported "
            f"enrichment p-values in this audit (items 8, 11, 12). Full detail: "
            f"`audit/details/item13_definitive_reconciliation.json` and `item13_definitive_table.csv`.\n\n"
            f"**Frozen mapping** (exact reaction IDs, all 5 axes):\n\n"
        )
        for axis, rxns in r13["frozen_mapping"].items():
            lines.append(f"- `{axis}`: {rxns}\n")
        lines.append(
            f"\n{r13['mapping_note']}\n\n"
            f"**Methodology**: {r13['methodology']['procedure']}. "
            f"N={r13['methodology']['n_permutations']:,} permutations per group-level test, seed={r13['methodology']['seed']}, "
            f"resolution floor=1/{r13['methodology']['n_permutations']+1}={r13['methodology']['resolution_floor']:.3g}, "
            f"empirical-p formula: `{r13['methodology']['p_value_formula']}`.\n\n"
            f"**Independence diagnostic**: within-cohort mean correlation={r13['independence_diagnostic']['within_cohort_mean_correlation']:.4f}, "
            f"across-cohort mean correlation={r13['independence_diagnostic']['across_cohort_mean_correlation']:.4f} "
            f"(Welch t-test p={r13['independence_diagnostic']['welch_t_p_value']:.3g}). "
            f"Effective number of independent tests (Li & Ji 2005 eigenvalue method): "
            f"**{r13['independence_diagnostic']['effective_n_independent_tests_li_ji']:.2f} of 10 naive**. "
            f"Fisher's method flagged as anti-conservative: **{r13['independence_diagnostic']['fisher_anti_conservative_flagged']}**. "
            f"Both the naive Fisher result AND a Brown (1975)/Kost-McDermott moment-matched correction "
            f"(which rescales the combined statistic itself under the measured correlation structure, not "
            f"just the degrees of freedom) are reported below, side by side.\n\n"
            f"**DEFINITIVE TABLE — naive Fisher (df=20, assumes 10 independent tests):**\n\n"
            f"| Axis | Raw combined p | BH-adjusted p | Significant |\n|---|---|---|---|\n"
        )
        for row in r13["definitive_table"]:
            lines.append(
                f"| {row['axis']} | {row['raw_combined_p_naive_fisher']:.4g} | "
                f"{row['bh_adjusted_p_naive_fisher']:.4g} | {row['significant_naive_YN']} |\n"
            )
        lines.append(
            f"\n**Brown-corrected table (accounts for measured non-independence):**\n\n"
            f"| Axis | Raw combined p | BH-adjusted p | Significant |\n|---|---|---|---|\n"
        )
        for row in r13["definitive_table"]:
            lines.append(
                f"| {row['axis']} | {row['raw_combined_p_independence_corrected']:.4g} | "
                f"{row['bh_adjusted_p_independence_corrected']:.4g} | {row['significant_independence_corrected_YN']} |\n"
            )
        lines.append(f"\n{r13['bcaa_framing_note']}\n")
        lines.append(f"\n{r13['supersedes']}\n")

    r14 = summary["items"].get("14") or summary["items"].get(14)
    if isinstance(r14, dict) and "per_fold_analysis" in r14:
        lines.append("### 14. Unconstrained (0-1) LOCO/LOBO weight re-fit\n")
        lines.append(
            f"Manuscript files were NOT modified. Same grid-search/nested-CV code as item 5, floor "
            f"lowered from 0.15 to 0.0 (grid size {r14['methodology']['grid_size']} vs. item 5's 78). "
            f"Full detail: `audit/details/item14_unconstrained_weight_refit.json`.\n\n"
            f"**{r14['n_folds_matching_reported']}/{r14['n_folds_total']} folds match the reported "
            f"(0.15, 0.15, 0.70)** -- none do. Every fold instead lands with w_C (consensus) MUCH "
            f"higher than 0.15 (0.35-0.6), w_F (flux confidence) at the grid's smallest step (0.05, "
            f"i.e. wanting to go below the old floor), and w_P (perturbation) correspondingly lower "
            f"than 0.70 (0.35-0.6).\n\n"
            f"| Fold | Unconstrained (w_C, w_F, w_P) |\n|---|---|\n"
        )
        for a in r14["per_fold_analysis"]:
            lines.append(f"| {a['fold']} | {a['unconstrained_weights']} |\n")
        lines.append(
            f"\n**Precision on \"was the floor binding\"**: the floor bound asymmetrically, not "
            f"uniformly. It bound HARD on w_F (flux confidence) -- the unconstrained search "
            f"consistently wants w_F near 0 (0.05, the smallest nonzero grid step), well below the "
            f"old 0.15 floor. It did NOT bind on w_C in the direction of holding it down -- the "
            f"unconstrained search wants w_C substantially ABOVE 0.15 (0.35-0.6), i.e. the floor's "
            f"reported value of 0.15 for consensus was actually much LOWER than what the objective "
            f"prefers. What happened under the floor-constrained search: forcing w_F up to its 0.15 "
            f"minimum (which the objective doesn't want at all) ate into the weight budget, and the "
            f"remaining 0.85 split between w_C and w_P within the floor's restricted region landed at "
            f"0.15/0.70 as the best AVAILABLE combination -- not because w_C wanted to be low, but "
            f"because the floor on w_F left less room and the grid's floor-constrained search space "
            f"didn't include the true unconstrained preference (high w_C, near-zero w_F, moderate w_P). "
            f"**Conclusion: the floor was actively constraining the result, and the true unconstrained "
            f"preference is a substantially different weighting (consensus-dominant, not "
            f"perturbation-dominant) than what's reported.**\n"
        )

    r15 = summary["items"].get("15") or summary["items"].get(15)
    if isinstance(r15, dict) and "code_citation_frozen_procedure" in r15:
        lines.append("### 15. Null-sampling verification\n")
        cit = r15["code_citation_frozen_procedure"]
        lines.append(
            f"**{r15['answer']}**\n\n"
            f"Code citation (`{cit['file']}`, function `{cit['function']}`):\n\n"
            f"```python\n"
        )
        for line_no, code in cit["lines"].items():
            lines.append(f"{line_no}: {code}\n")
        lines.append(f"```\n\nPer-axis k used in item13's definitive pass: {r15['per_axis_k_used_in_item13']}\n")

    r16 = summary["items"].get("16") or summary["items"].get(16)
    if isinstance(r16, dict) and "restricted_fit" in r16:
        lines.append("### 16. Growth-support confound check on item14's consensus-dominant finding\n")
        pc = r16["premise_check"]
        lines.append(
            f"Manuscript files were NOT modified.\n\n"
            f"**Premise check**: the request's claim was 3 axes overlap growth-support (SHMT, DNL, "
            f"choline_pc) vs. 2 that don't (urea_cycle, bcaa). Actual finding under item14's mapping: "
            f"only 2 axes overlap (SHMT 1/2 via MAR03845, choline_pc 1/5 via MAR00653); DNL (the "
            f"FASN-inclusive mapping item14 used) has 0/3 overlap -- it only overlaps (1/3, MAR00148) "
            f"under the alternate SCD1-only mapping used in items 12/13, which item14 did not use. "
            f"urea_cycle and bcaa are unambiguously non-overlapping under either DNL variant, so they "
            f"remain the clean holdout set.\n\n"
            f"**Comparison** (mean weights, w_C/w_F/w_P):\n\n"
            f"| Fit | w_C | w_F | w_P |\n|---|---|---|---|\n"
            f"| 5-axis (item14) | {r16['comparison']['mean_weights_5axis_item14_wC_wF_wP'][0]} | "
            f"{r16['comparison']['mean_weights_5axis_item14_wC_wF_wP'][1]} | "
            f"{r16['comparison']['mean_weights_5axis_item14_wC_wF_wP'][2]} |\n"
            f"| 2-axis (urea_cycle+bcaa only) | {r16['comparison']['mean_weights_2axis_restricted_wC_wF_wP'][0]} | "
            f"{r16['comparison']['mean_weights_2axis_restricted_wC_wF_wP'][1]} | "
            f"{r16['comparison']['mean_weights_2axis_restricted_wC_wF_wP'][2]} |\n\n"
            f"{r16['comparison']['note_on_prior_auto_generated_conclusion']}\n\n"
            f"{r16['conclusion']}\n"
        )

    r17 = summary["items"].get("17") or summary["items"].get(17)
    if isinstance(r17, dict) and "per_axis_weights" in r17:
        lines.append("### 17. Per-axis unconstrained weight preferences\n")
        lines.append(
            f"Manuscript files were NOT modified. Single-axis unconstrained (floor=0) production fits "
            f"for each of the 5 biomarker axes individually. urea_cycle and bcaa carried forward from "
            f"item 16's single-axis legs and independently cross-checked here (both reproduce exactly). "
            f"DNL uses the same FASN-inclusive mapping as items 14/16 (not the SCD1-only mapping from "
            f"items 12/13).\n\n"
            f"| Axis | w_C | w_F | w_P | Dominant |\n|---|---|---|---|---|\n"
        )
        for axis, info in r17["per_axis_weights"].items():
            w = info["weights"]
            dom = r17["dominant_component_by_axis"][axis]
            lines.append(f"| {axis} | {w[0]} | {w[1]} | {w[2]} | {dom} |\n")
        lines.append(
            f"\n**Pattern**: serine_glycine_shmt and choline_pc -- the two axes with a growth-support-"
            f"forced reaction in their mapped set (item 16) -- land on IDENTICAL weights (0.10, 0.05, "
            f"0.85), both strongly perturbation-dominant. DNL (FASN-inclusive, no growth-support "
            f"overlap) is strongly consensus-dominant (0.70, 0.00, 0.30) -- the opposite pattern. BCAA "
            f"is flux-confidence-dominant (0.25, 0.70, 0.05). urea_cycle is the only axis with no single "
            f"dominant component (0.45, 0.30, 0.25), the closest thing to a \"balanced\" individual "
            f"signal among the five. In short: 5 axes produce (at most) 4 distinct preference profiles, "
            f"and the only pair that agrees (SHMT/choline_pc) does so for a plausible mechanical reason "
            f"(shared growth-support-forced reaction) rather than a shared biological signal. No single "
            f"weight vector is preferred by more than 2 of the 5 axes -- this is consistent with, and "
            f"extends, items 14/16's finding that the reported production weights are an artifact of "
            f"averaging together axes with genuinely incompatible individual preferences, not a "
            f"consensus any individual axis actually supports.\n"
        )

    r18 = summary["items"].get("18") or summary["items"].get(18)
    if isinstance(r18, dict) and "kendalls_w_3axis_uncontaminated" in r18:
        lines.append("### 18. Kendall's W and SD-ordering benchmark, excluding contaminated axes\n")
        pc = r18["premise_check"]
        lines.append(
            f"Manuscript files were NOT modified.\n\n"
            f"**Premise check**: \"0.7912\" does not appear anywhere in this repository (direct grep, "
            f"confirmed). There is no subsystem-level Kendall's W anywhere in this repository -- "
            f"`docs/MANUSCRIPT_FIGURES.md` explicitly states that entire figure is unbuilt (\"No script "
            f"in this repository aggregates confidence/calibration scores by subsystem\"). The only "
            f"Kendall's W that actually exists is **biomarker-level** (5 biomarkers x 3 cohorts): "
            f"**W={pc['actual_kendalls_w_in_repo']['value']}**, from `data/calibration/figure3_statistics.json` "
            f"(produced by `scripts/12_figure3_cross_cohort_confidence.py`). There is also no single "
            f"\"master architectural benchmark (Table 3 IQR comparison)\" literally in the repo -- the "
            f"closest existing, validated analog is item 9's SD-sensitivity check.\n\n"
            f"**Kendall's W, restricted to the 3 uncontaminated axes** (urea_cycle, bcaa, dnl_scd1_fasn): "
            f"**W={r18['kendalls_w_3axis_uncontaminated']['kendalls_w']:.4f}** "
            f"(Friedman chi2={r18['kendalls_w_3axis_uncontaminated']['friedman_chi2']:.4f}, "
            f"p={r18['kendalls_w_3axis_uncontaminated']['friedman_p']:.4f}), down from 0.6667 with all 5 "
            f"axes. {r18['kendalls_w_conclusion']}\n\n"
            f"**SD ordering (equal-weighting vs. calibrated), restricted to the 3 uncontaminated axes' 18 "
            f"mapped reactions**: equal={r18['sd_ordering_3axis_uncontaminated_reactions_only']['equal_weighting_sd']:.4f}, "
            f"calibrated={r18['sd_ordering_3axis_uncontaminated_reactions_only']['calibrated_sd_direct_column']:.4f}. "
            f"{r18['sd_ordering_conclusion']}\n"
        )

    r19 = summary["items"].get("19") or summary["items"].get(19)
    if isinstance(r19, dict) and "fresh_100k_result" in r19:
        lines.append("### 19. Final SHMT p=2.1e-5 reconciliation attempt\n")
        lines.append(
            f"Manuscript files were NOT modified.\n\n"
            f"**{r19['already_tried_note']}**\n\n"
            f"Recomputed fresh at N=100,000 (matching item 13's frozen-procedure precision) as a clean, "
            f"dedicated citation:\n\n"
            f"- Fresh N=100,000 raw combined p: **{r19['fresh_100k_result']['p_value']:.4g}**\n"
            f"- Prior N=20,000 (item 11) raw p: {r19['prior_20k_result_cited_from_item11']['raw_p']:.4g}, "
            f"BH-adjusted: {r19['prior_20k_result_cited_from_item11']['bh_adjusted_p']:.4g}\n"
            f"- Reported value: {r19['reported_value']:.4g}\n"
            f"- Ratio (N=100,000 vs. reported): **{r19['ratio_100k_vs_reported']:.2f}x** "
            f"(within order-of-magnitude threshold: {r19['within_order_of_magnitude_100k']})\n\n"
            f"{r19['conclusion']}\n"
        )

    with open(readme_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Appended new section to {readme_path}")


def write_readme(rows, summary):
    lines = []
    lines.append("# GEM1 Independent Arithmetic and Provenance Audit\n")
    lines.append(f"Generated: {summary['generated_at']}\n")
    lines.append(
        "This audit reads only raw source outputs (Human-GEM.json, Step 5 troppo JSONs, "
        "GEO expression matrices, all_groups_calibrated_with_metadata.csv) and independently "
        "recomputes every value from scratch. It never imports GEM1's own scoring/calibration "
        "modules (scripts/04_consensus_scoring.py, 05_confidence_engine.py, 09_calibration.py, "
        "10_join_calibration_metadata.py, 11_biomarker_ranking_by_group.py). No existing GEM1 "
        "output file was modified by this audit.\n"
    )
    lines.append("## Summary\n")
    lines.append("| # | Item | Status | Max abs. discrepancy |")
    lines.append("|---|------|--------|----------------------|")
    for r in rows:
        lines.append(f"| {r['item']} | {r['title']} | **{r['status']}** | {r['max_abs_discrepancy']} |")
    lines.append("")

    lines.append("## Item details\n")

    # Item 1
    lines.append("### 1. Single-Algorithm (iMAT+FVA) external baseline\n")
    r1 = summary["items"].get(1)
    if isinstance(r1, dict) and "baseline_biomarker_metric" in r1:
        lines.append(
            f"Independently reconstructed from each group's raw iMAT-only Step 5 output "
            f"(`data/context_specific_models/*__imat.json`), plus growth-support union, with a "
            f"fresh cobra/Gurobi FBA+FVA call (no `06_flux_analysis.py`). "
            f"**No manuscript/reported reference value exists anywhere in this repository for this "
            f"baseline** (confirmed by exhaustive search) -- so PASS/FAIL here reflects whether the "
            f"independent reconstruction completed successfully on all {r1.get('n_groups_expected','?')} "
            f"groups, not agreement with a manuscript figure.\n\n"
            f"- Baseline biomarker metric (single-algorithm, presence-only signal): "
            f"**{r1['baseline_biomarker_metric']['value']:.4f}**\n"
            f"- GEM1 full-pipeline biomarker metric (same aggregation rule, `calibrated_confidence_score`): "
            f"**{r1['gem1_full_pipeline_biomarker_metric']['value']:.4f}**\n"
            f"- Delta (GEM1 - baseline): **{r1['delta_gem1_minus_baseline']:.4f}**\n"
        )
    else:
        lines.append("Not yet run -- long-running (~30 min), launch as a detached process (see below).\n")

    # Item 2
    lines.append("### 2. Table 3 IQR reconciliation (Flux Uncertainty-Only / Uncalibrated Equal Weighting)\n")
    r2 = summary["items"].get(2)
    if isinstance(r2, dict):
        gt = r2.get("ground_truth_search", {})
        lines.append(
            f"**Ground-truth search finding**: grepped every `.py`/`.ipynb`/`.R` file in the repo for "
            f"`\"0.2140\"`, `\"0.2410\"`, `\"iMAT+FVA\"`, `\"single-algorithm\"`, `\"master benchmark\"`, "
            f"`\"structural consensus-only\"` -- **{gt.get('matches_found', 0)} matches**. "
            f"{gt.get('finding', '')}\n\n"
            f"Since no script defines which reaction population / score representation the two "
            f"ablations refer to, computed IQR under every combination of proxy definition "
            f"(flux-uncertainty-only `w=(0,1,0)` vs equal-weight `w=(1/3,1/3,1/3)`), score "
            f"representation (raw weighted-rank-sum vs evaluation-percentile), and reaction "
            f"population (active-only vs all-rows-zero-filled, pooled vs per-group-mean) -- "
            f"see `audit/details/item2_iqr_reconciliation.json` for the full 16-cell grid.\n\n"
            f"- Reported (manuscript text, unverifiable in-repo): flux-only=**{r2['reported_manuscript_values']['flux_uncertainty_only_iqr']}**, "
            f"equal-weight=**{r2['reported_manuscript_values']['uncalibrated_equal_weighting_iqr']}**\n"
            f"- Prior independent reconstruction: flux-only=**{r2['prior_independent_reconstruction_values']['flux_uncertainty_only_iqr']}**, "
            f"equal-weight=**{r2['prior_independent_reconstruction_values']['uncalibrated_equal_weighting_iqr']}**\n"
            f"- Grid cells matching reported flux-only value (±0.01): {list(r2['cells_matching_reported_flux_only_iqr'].keys()) or 'none'}\n"
            f"- Grid cells matching reported equal-weight value (±0.01): {list(r2['cells_matching_reported_equal_weight_iqr'].keys()) or 'none'}\n"
            f"- Closest cell to reported flux-only ({r2['reported_manuscript_values']['flux_uncertainty_only_iqr']}): "
            f"`{r2['closest_cell_to_reported_flux_only_iqr']['cell']}` = "
            f"{r2['closest_cell_to_reported_flux_only_iqr']['value']:.4f} "
            f"(abs diff {r2['closest_cell_to_reported_flux_only_iqr']['abs_diff']:.4f})\n"
            f"- Closest cell to reported equal-weight ({r2['reported_manuscript_values']['uncalibrated_equal_weighting_iqr']}): "
            f"`{r2['closest_cell_to_reported_equal_weight_iqr']['cell']}` = "
            f"{r2['closest_cell_to_reported_equal_weight_iqr']['value']:.4f} "
            f"(abs diff {r2['closest_cell_to_reported_equal_weight_iqr']['abs_diff']:.4f})\n"
            f"- Cells matching the *prior* independent reconstruction (0.4765 / 0.1683), confirming that "
            f"earlier attempt used 'active-only reactions, raw weighted-rank score, pooled across groups': "
            f"flux-only={list(r2['cells_matching_prior_flux_only_iqr'].keys())}, "
            f"equal-weight={list(r2['cells_matching_prior_equal_weight_iqr'].keys())}\n"
        )
        if r2.get("anomaly_note"):
            lines.append(
                f"\n> **[ANOMALY FLAGGED]**: {r2['anomaly_note']}\n"
            )

    # Item 3
    lines.append("### 3. Flux sampling + perturbation robustness re-execution\n")
    r3 = summary["items"].get(3)
    if isinstance(r3, dict) and "tolerance_summary" in r3:
        lines.append(
            f"Consensus-active reaction sets independently recomputed from raw Step 5 JSONs "
            f"(not read from `data/consensus_scores/*.json`); OptGP sampling + Monte Carlo "
            f"perturbation re-executed fresh (`{r3['n_groups_processed']}/{r3['n_groups_expected']}` "
            f"groups completed). Compared against `flux_confidence_raw`/`perturbation_robustness_raw` "
            f"using distribution/CI-based tolerance ({r3['tolerance_criteria']}).\n\n"
        )
        for metric, s in r3["tolerance_summary"].items():
            lines.append(
                f"- **{metric}** (n={s['n']}): mean diff={s['mean_diff']:.4f}, "
                f"mean |diff|={s['mean_abs_diff']:.4f}, max |diff|={s['max_abs_diff']:.4f}, "
                f"Spearman r={s['spearman_r']:.3f}, 95% CI on mean diff={s['ci95_mean_diff']}\n"
            )
        if "status_detail" in r3:
            lines.append("\n**Per-metric verdicts** (a single PASS/FAIL would hide a real finding here):\n")
            for metric, verdict in r3["status_detail"].items():
                lines.append(f"- *{metric}*: {verdict}\n")
        td = r3.get("tolerance_summary", {}).get("perturbation_robustness_raw", {}).get("tie_decomposition")
        if td:
            lines.append(
                f"\n> **[FINDING]** {td['finding']}\n"
                f">\n> Tie decomposition: {td['n_tied_at_1.0_both']} of "
                f"{td['n_tied_at_1.0_both'] + td['n_non_tied']} rows are growth-support reactions "
                f"deterministically fixed at 1.0 in both runs (exact ties by construction). Among the "
                f"remaining {td['n_non_tied']} reactions: Pearson r={td['non_tied_pearson_r']:.4f}, "
                f"Spearman r={td['non_tied_spearman_r']:.4f}, mean |diff|={td['non_tied_mean_abs_diff']:.4f}.\n"
            )
    else:
        lines.append("Not yet run -- long-running (~75 min), launch as a detached process (see below).\n")

    # Item 4
    lines.append("### 4. Growth-support LP rebuild\n")
    r4 = summary["items"].get(4)
    if isinstance(r4, dict):
        lines.append(
            f"Rebuilt directly from `models/Human-GEM.json` via a fresh `scipy.optimize.linprog` "
            f"(HiGHS) call -- no cobra, no troppo.\n\n"
            f"- Unconstrained max biomass: independent=**{r4['independent_max_biomass']}**, "
            f"reported=**{r4['reported_max_biomass']}**, diff={r4['discrepancy_max_biomass']:.2e}\n"
            f"- Achieved biomass @ 90% floor: independent=**{r4['independent_achieved_biomass']}**, "
            f"reported=**{r4['reported_achieved_biomass']}**, diff={r4['discrepancy_achieved_biomass']:.2e}\n"
            f"- Support-set size: independent=**{r4['independent_n_support']}**, reported=**{r4['reported_n_support']}**, "
            f"Jaccard overlap={r4['support_set_jaccard_overlap']:.3f} (informational -- L1-min solutions aren't unique)\n"
        )

    # Item 5
    lines.append("### 5. LOCO/LOBO weight optimization refit\n")
    r5 = summary["items"].get(5)
    if isinstance(r5, dict):
        lines.append(
            f"Independent grid search (floor=0.15, step=0.05, {r5['weight_grid_size']} candidate weight "
            f"triples) + nested LOCO (3 outer folds) / LOBO (5 outer folds) CV, from raw rank columns "
            f"+ biomarker labels only.\n\n"
            f"- Independent production weights (no holdout): **{r5['independent_production_weights']}** "
            f"(reported: {r5['reported_production_weights']})\n"
            f"- All LOCO+LOBO folds converge to reported weights: **{r5['all_folds_converge_to_reported_weights']}**\n"
            f"- Max abs weight discrepancy: {r5['max_abs_weight_discrepancy']:.4f}\n"
        )

    # Item 6
    lines.append("### 6. Active-reaction counts via independent troppo reconstruction\n")
    r6 = summary["items"].get(6)
    if isinstance(r6, dict) and "n_combinations_reconstructed" in r6:
        lines.append(
            f"Fresh troppo-based extraction (own code, not importing `03_build_context_specific_models.py`), "
            f"from raw GEO expression matrices + Human-GEM.json, with a freshly-recomputed FastCC "
            f"consistent set ({r6.get('fastcc_consistent_reactions','?')} reactions).\n\n"
            f"- Combinations reconstructed: {r6['n_combinations_reconstructed']}/{r6['n_combinations_total']}\n"
            f"- Within tolerance ({r6['tolerances_used']}): {r6['n_within_tolerance']}/{r6['n_combinations_total']}\n"
            f"- Max abs discrepancy (active-reaction count): {r6['max_abs_discrepancy']}\n"
            f"- Full per-combination comparison: `audit/details/item6_comparison.csv`\n"
        )
    else:
        lines.append("Not yet run -- long-running (~70-90 min), launch as a detached process (see below).\n")

    lines.append("## How to run (long-running items)\n")
    lines.append(
        "Paths below are relative to the repository root and the conda env names "
        "used during development (`GEM1` for the main pipeline env, `gem1-troppo` "
        "for the Step 5/extraction env) -- substitute `<path-to-conda-envs>` for "
        "your local conda installation's env directory (typically found via "
        "`conda env list`).\n"
    )
    lines.append("```powershell")
    lines.append('# Item 1 (GEM1/gem1-main env, ~30 min)')
    lines.append('Start-Process -FilePath "<path-to-conda-envs>\\GEM1\\python.exe" '
                  '-ArgumentList "item1_imat_fva_baseline.py" -WorkingDirectory "audit" '
                  '-RedirectStandardOutput "audit\\logs\\item1.log" '
                  '-RedirectStandardError "audit\\logs\\item1_err.log" -WindowStyle Hidden')
    lines.append('')
    lines.append('# Item 3 (GEM1/gem1-main env, ~75 min)')
    lines.append('Start-Process -FilePath "<path-to-conda-envs>\\GEM1\\python.exe" '
                  '-ArgumentList "item3_flux_sampling_perturbation.py" -WorkingDirectory "audit" '
                  '-RedirectStandardOutput "audit\\logs\\item3.log" '
                  '-RedirectStandardError "audit\\logs\\item3_err.log" -WindowStyle Hidden')
    lines.append('')
    lines.append('# Item 6 (gem1-troppo env, ~70-90 min)')
    lines.append('Start-Process -FilePath "<path-to-conda-envs>\\gem1-troppo\\python.exe" '
                  '-ArgumentList "item6_troppo_reextract.py" -WorkingDirectory "audit" '
                  '-RedirectStandardOutput "audit\\logs\\item6.log" '
                  '-RedirectStandardError "audit\\logs\\item6_err.log" -WindowStyle Hidden')
    lines.append("```")
    lines.append("\n(Run these from the repository root, or adjust `-WorkingDirectory` accordingly.) "
                  "Check progress with `Get-Content .\\audit\\logs\\item1.log -Tail 30`, confirm alive with "
                  "`Get-Process python`. Re-run `python independent_audit.py --finalize` after each completes "
                  "to fold its result into `audit_report.csv`/`audit_summary.json`/this README.\n")

    with open(os.path.join(AUDIT_DIR, "README.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote {os.path.join(AUDIT_DIR, 'README.md')}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-fast", action="store_true", help="Run items 2, 4, 5 (synchronous)")
    parser.add_argument("--finalize", action="store_true", help="Aggregate audit/details/*.json into the final report (overwrites, items 1-6)")
    parser.add_argument("--finalize-append", type=str, default=None,
                         help="Comma-separated item numbers to append to the existing report (e.g. 7,8,9,10)")
    args = parser.parse_args()

    if not args.run_fast and not args.finalize and not args.finalize_append:
        parser.print_help()
        return

    if args.run_fast:
        run_fast_items()
    if args.finalize:
        finalize()
    if args.finalize_append:
        item_nums = [int(x) for x in args.finalize_append.split(",")]
        finalize_append(item_nums)


if __name__ == "__main__":
    main()
