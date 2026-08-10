"""
Item 4 -- Independently re-derive the growth-support LP (519 reactions, 90%
biomass floor) directly from models/Human-GEM.json using a fresh scipy/HiGHS
solver call. No troppo, no cobra, no GEM1 pipeline code -- builds the
stoichiometric matrix and both LPs by hand, exactly mirroring (but not
importing) the derivation documented in docs/GEM1-roadmap-schedule.md.

Fast (~seconds); run synchronously, no detached process needed.
"""

import json
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix, hstack, eye

from audit_common import load_human_gem_raw, load_growth_support_ids, BIOMASS_REACTION, DETAILS_DIR
import os

REPORTED_MAX_BIOMASS = 124.8681483774457
REPORTED_ACHIEVED_BIOMASS = 112.38133353970113
REPORTED_N_SUPPORT = 519
OBJ_FRAC = 0.9
SUPPORT_EPS = 1e-6


def build_S(model_json):
    mets = model_json["metabolites"]
    rxns = model_json["reactions"]
    met_idx = {m["id"]: i for i, m in enumerate(mets)}
    n_met, n_rxn = len(mets), len(rxns)
    S = lil_matrix((n_met, n_rxn))
    lb = np.zeros(n_rxn)
    ub = np.zeros(n_rxn)
    rxn_ids = []
    for j, r in enumerate(rxns):
        rxn_ids.append(r["id"])
        lb[j] = r.get("lower_bound", 0.0)
        ub[j] = r.get("upper_bound", 0.0)
        for met_id, coeff in r["metabolites"].items():
            S[met_idx[met_id], j] = coeff
    return csr_matrix(S), lb, ub, rxn_ids


def run():
    t0 = time.time()
    model_json = load_human_gem_raw()
    S, lb, ub, rxn_ids = build_S(model_json)
    n_met, n_rxn = S.shape
    obj_idx = rxn_ids.index(BIOMASS_REACTION)
    print(f"Built S: {n_met} metabolites x {n_rxn} reactions, {S.nnz} nonzeros.")

    # --- Step A: unconstrained FBA maximizing biomass ---
    c = np.zeros(n_rxn)
    c[obj_idx] = -1.0  # linprog minimizes
    res_fba = linprog(c, A_eq=S, b_eq=np.zeros(n_met), bounds=list(zip(lb, ub)), method="highs")
    if not res_fba.success:
        raise RuntimeError(f"Unconstrained FBA failed: {res_fba.message}")
    max_biomass = -res_fba.fun
    print(f"Unconstrained max biomass (fresh LP): {max_biomass}")

    # --- Step B: L1-minimal flux subject to biomass >= OBJ_FRAC * max_biomass ---
    target = OBJ_FRAC * max_biomass
    # split-variable trick: v = vp - vn, vp,vn >= 0
    vp_ub = np.clip(ub, 0, None)
    vn_ub = np.clip(-lb, 0, None)
    bounds2 = [(0, u) for u in vp_ub] + [(0, u) for u in vn_ub]
    c2 = np.ones(2 * n_rxn)  # minimize sum(vp) + sum(vn)

    Sp = S
    Sn = -S
    A_eq2 = hstack([Sp, Sn])
    b_eq2 = np.zeros(n_met)

    # biomass row: vp[obj] - vn[obj] >= target  ->  -vp[obj] + vn[obj] <= -target
    row = lil_matrix((1, 2 * n_rxn))
    row[0, obj_idx] = -1.0
    row[0, n_rxn + obj_idx] = 1.0
    A_ub2 = csr_matrix(row)
    b_ub2 = np.array([-target])

    res_l1 = linprog(c2, A_eq=A_eq2, b_eq=b_eq2, A_ub=A_ub2, b_ub=b_ub2, bounds=bounds2, method="highs")
    if not res_l1.success:
        raise RuntimeError(f"L1-min growth-support LP failed: {res_l1.message}")

    vp = res_l1.x[:n_rxn]
    vn = res_l1.x[n_rxn:]
    v = vp - vn
    achieved_biomass = v[obj_idx]
    support_mask = (vp + vn) > SUPPORT_EPS
    support_ids = sorted(rxn_ids[i] for i in np.where(support_mask)[0])
    n_support = len(support_ids)

    print(f"Achieved biomass at 90% floor: {achieved_biomass}")
    print(f"Independent support-set size: {n_support} (reported: {REPORTED_N_SUPPORT})")

    stored = load_growth_support_ids()
    stored_ids = set(stored["support_reaction_ids"])
    my_ids = set(support_ids)
    jaccard = len(stored_ids & my_ids) / len(stored_ids | my_ids)

    disc_max_biomass = abs(max_biomass - REPORTED_MAX_BIOMASS)
    disc_achieved = abs(achieved_biomass - REPORTED_ACHIEVED_BIOMASS)
    max_abs_discrepancy = max(disc_max_biomass, disc_achieved)

    result = {
        "item": 4,
        "description": "Growth-support LP (519 reactions, 90% biomass floor) rebuilt from Human-GEM.json via fresh scipy/HiGHS LP",
        "independent_max_biomass": max_biomass,
        "reported_max_biomass": REPORTED_MAX_BIOMASS,
        "discrepancy_max_biomass": disc_max_biomass,
        "independent_achieved_biomass": achieved_biomass,
        "reported_achieved_biomass": REPORTED_ACHIEVED_BIOMASS,
        "discrepancy_achieved_biomass": disc_achieved,
        "independent_n_support": n_support,
        "reported_n_support": REPORTED_N_SUPPORT,
        "n_support_diff": n_support - REPORTED_N_SUPPORT,
        "support_set_jaccard_overlap": jaccard,
        "max_abs_discrepancy": max_abs_discrepancy,
        "pass_threshold": 1e-3,
        "status": "PASS" if max_abs_discrepancy < 1e-3 else "FAIL",
        "note": (
            "L1-minimal support sets are not guaranteed unique under degenerate "
            "alternate optima -- exact reaction-count/identity match is informational, "
            "not part of the PASS/FAIL gate. PASS/FAIL is gated on the two LP objective "
            "values (unconstrained max biomass, achieved 90%-floor biomass), which ARE "
            "unique for a well-posed LP."
        ),
        "elapsed_seconds": time.time() - t0,
    }
    out_path = os.path.join(DETAILS_DIR, "item4_growth_support_lp.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {out_path}")
    print(json.dumps({k: v for k, v in result.items() if k not in ("note",)}, indent=2))
    return result


if __name__ == "__main__":
    run()
