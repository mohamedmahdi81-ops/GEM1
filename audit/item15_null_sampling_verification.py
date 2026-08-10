"""
Item 15 -- Null-sampling verification for the frozen enrichment procedure
(item13_definitive_reconciliation.py, the S8.2 / definitive naive-vs-Brown-
corrected 5-axis table). Quick code-reading confirmation, NOT a rerun: does
the null distribution for each biomarker axis use reaction-count-matched
sampling (draw size = that axis's own number of mapped reactions), or some
other fixed/different sample size?

Writes only its own small detail JSON to audit/details/; does not modify any
existing GEM1 output file or re-execute any computation.
"""

import json
import os
import time

from audit_common import DETAILS_DIR

CITATION = {
    "file": "item13_definitive_reconciliation.py",
    "function": "flat_test_100k",
    "lines": {
        216: 'idx = [r_to_idx[r] for r in rxn_ids if r in r_to_idx]',
        217: 'k = len(idx)',
        221: 'random_idx = rng.integers(0, n_reactions, size=(n_perm, k))',
    },
}

CROSS_CHECK_CITATION = {
    "file": "item11_calibrated_enrichment_audit.py",
    "functions": ["flat_test (lines 140-143)", "nested_test (lines 158-163)"],
    "same_pattern": True,
    "note": "item11 (the earlier exploratory flat/nested comparison) uses the identical "
            "idx=[...]; k=len(idx); rng.integers(0, n_reactions, size=(n_perm, k)) pattern -- "
            "confirms this was not a one-off in the frozen item13 script, it's this audit's "
            "consistent convention throughout.",
}

PER_AXIS_K = {
    "serine_glycine_shmt": 2, "urea_cycle": 6, "bcaa": 9, "dnl_scd1": 3, "choline_pc": 5,
}


def run():
    t0 = time.time()
    result = {
        "item": 15,
        "description": "Null-sampling verification for the frozen enrichment procedure (item13) -- confirms reaction-count-matched null draws",
        "question": "Is the null distribution for each biomarker axis drawn as a reaction-count-matched sample (N mapped reactions -> N-reaction null draws), or some other way?",
        "answer": "CONFIRMED: reaction-count-matched. k (the null draw's reaction-set size) is set to len(idx), where idx is built directly from that axis's own mapped reaction IDs (rxn_ids) via r_to_idx lookup -- not a fixed constant, not a different pool-derived size.",
        "code_citation_frozen_procedure": CITATION,
        "cross_check_citation": CROSS_CHECK_CITATION,
        "per_axis_k_used_in_item13": PER_AXIS_K,
        "verification_method": "Direct source-code reading (grep + manual line inspection) of the actual .py files used to produce item13's definitive table -- not a rerun, not a recomputation.",
        "status": "CONFIRMED",
        "max_abs_discrepancy": "N/A -- code-citation confirmation, not a numerical reconciliation",
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item15_null_sampling_verification.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"\nAnswer: {result['answer']}")
    print(f"\nCitation ({CITATION['file']}, function {CITATION['function']}):")
    for line_no, code in CITATION["lines"].items():
        print(f"  {line_no}: {code}")
    return result


if __name__ == "__main__":
    run()
