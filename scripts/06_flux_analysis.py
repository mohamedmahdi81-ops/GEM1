"""
GEM1 - Step 6: Flux analyses (FBA, pFBA, FVA) on each group's strict-consensus
context-specific model.

Design decisions (confirmed with Mohammed before implementation, 2026-07-17):

1. Constrained model, not reaction deletion. Reactions NOT in a group's
   strict consensus set have both bounds set to (0, 0) on a copy of the full
   Human-GEM base model, rather than being removed from the model. Deleting
   reactions can silently orphan metabolites that only participate in
   deleted reactions and mask mass-balance structure; zeroing bounds keeps
   the full stoichiometric network (all 8,461 metabolites) intact and only
   blocks flux through reactions the consensus doesn't support. Reactions
   that stay active keep their original Human-GEM bounds unchanged.

2. One model per GROUP (10 total), built from
   data/consensus_scores/{cohort}__{group}__consensus.json's
   strict_consensus_reactions -- the exact reaction set Step 9/10 already
   scored (final_confidence_score == 1.0 reactions), not a per-algorithm
   reconstruction. For GSE126848/healthy this file's "strict" set is already
   the reduced (iMAT+tINIT only) definition from Step 9 -- reused as-is, and
   consensus_definition_used is carried into this step's output so the
   reduced group stays flagged all the way through the pipeline.

Objective: Human-GEM's own built-in objective, MAR13082 ("generic human cell
biomass reaction") -- not overridden.

FVA: cobra's own default (fraction_of_optimum=1.0), run over the full
unblocked reaction set (strict consensus plus the growth-support additions
below) -- scoring flux ranges for reactions already forced to zero is not
informative.

Run from gem1-main -- standard cobra (0.31.1 in this env), Gurobi already the
default solver here (fast for LP; no MILP in this step, so no Gurobi
gap/time-limit patch is needed, unlike Step 5).

2026-07-18 growth-support fix (confirmed with Mohammed before implementation,
same day): before running Step 6 for real, directly tested FBA on each
group's strict_consensus_reactions set as this script builds it -- 9 of 10
groups came back with biomass flux exactly 0.0 (only GSE126848/healthy, the
2-algorithm reduced group, sustained growth). Root-caused, NOT a Step 5 bug:
FASTCORE's own individual reaction set (untouched by the Step 5
growth-capability fix, which only ever applied to iMAT/tINIT) independently
sustains 123-125 biomass flux on its own, for all 9 groups it has a model
for -- so every one of the three core algorithms CAN individually sustain
growth. The strict 3-way intersection breaks growth anyway, because each
algorithm routes flux to biomass through a different specific combination of
reactions; intersecting drops the connecting reactions each one needed, even
though none of the three algorithms is individually deficient. This is a
Step 6 modeling problem, not a Step 5/9/10 defect -- those steps are locked
and were not touched.

Fix, scoped to this script only: `build_constrained_model` now also keeps
the same 519-reaction minimal growth-support set from
`data/context_specific_models/minimal_growth_support.json` unblocked,
alongside each group's strict-consensus set, purely so the LP this step
solves is growth-capable. This does not change Step 9/10's stored
consensus/confidence definitions in any way -- `n_strict_consensus_reactions`
in this step's own output still reports the original Step 9 count
unmodified; the additional unblocked reactions are reported separately
(`n_growth_support_added`) for transparency, and FVA is run over the full
unblocked set (consensus + additions) so the forced-open reactions' flux
ranges are visible too, not hidden.
"""

import os
import sys
import json
import argparse
import pandas as pd
import cobra
from cobra.flux_analysis import pfba, flux_variability_analysis

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "models", "Human-GEM.json")
CONSENSUS_DIR = os.path.join(BASE_DIR, "data", "consensus_scores")
GROWTH_SUPPORT_PATH = os.path.join(
    BASE_DIR, "data", "context_specific_models", "minimal_growth_support.json"
)
OUT_DIR = os.path.join(BASE_DIR, "data", "flux_analysis")
os.makedirs(OUT_DIR, exist_ok=True)


def load_base_model():
    return cobra.io.load_json_model(MODEL_PATH)


def load_growth_support():
    with open(GROWTH_SUPPORT_PATH) as f:
        return set(json.load(f)["support_reaction_ids"])


def load_strict_consensus(cohort, group):
    path = os.path.join(CONSENSUS_DIR, f"{cohort}__{group}__consensus.json")
    with open(path) as f:
        data = json.load(f)
    return set(data["strict_consensus_reactions"]), data["consensus_definition_used"]


def build_constrained_model(base_model, active_set, growth_support):
    unblocked = active_set | growth_support
    model = base_model.copy()
    n_blocked = 0
    for rxn in model.reactions:
        if rxn.id not in unblocked:
            rxn.bounds = (0, 0)
            n_blocked += 1
    return model, n_blocked, unblocked


def analyze_group(cohort, group, base_model, growth_support, run_fva=True):
    active_set, consensus_def = load_strict_consensus(cohort, group)
    model, n_blocked, unblocked = build_constrained_model(base_model, active_set, growth_support)
    n_growth_support_added = len(unblocked) - len(active_set)

    objective_rxn_ids = [r.id for r in model.reactions if r.objective_coefficient != 0]

    result = {
        "cohort": cohort,
        "group": group,
        "consensus_definition_used": consensus_def,
        "n_strict_consensus_reactions": len(active_set),
        "n_growth_support_added": n_growth_support_added,
        "n_active_reactions": len(unblocked),
        "n_blocked_reactions": n_blocked,
        "objective_reaction": ",".join(objective_rxn_ids),
    }

    fba_sol = model.optimize()
    result["fba_status"] = fba_sol.status
    result["fba_objective_value"] = (
        float(fba_sol.objective_value) if fba_sol.status == "optimal" else None
    )

    if fba_sol.status != "optimal" or not fba_sol.objective_value:
        print(f"  ** {cohort}/{group}: FBA status = {fba_sol.status!r}, "
              f"objective value = {fba_sol.objective_value!r} -- "
              f"STOPPING before pFBA/FVA for this group, not guessing around it.")
        result["pfba_status"] = None
        result["pfba_objective_value"] = None
        result["pfba_total_flux"] = None
        result["n_fva_reactions"] = 0
        return result, None, fba_sol

    pfba_sol = pfba(model)
    result["pfba_status"] = pfba_sol.status
    result["pfba_objective_value"] = float(pfba_sol.objective_value)
    result["pfba_total_flux"] = float(pfba_sol.fluxes.abs().sum())

    fva_df = None
    if run_fva:
        fva_df = flux_variability_analysis(model, reaction_list=sorted(unblocked))
        result["n_fva_reactions"] = len(fva_df)

    return result, fva_df, fba_sol


def write_group_output(cohort, group, result, fva_df, fba_sol):
    out_path = os.path.join(OUT_DIR, f"{cohort}__{group}__flux.json")
    payload = dict(result)
    if fva_df is not None:
        payload["fva"] = {
            rid: {"minimum": float(row["minimum"]), "maximum": float(row["maximum"])}
            for rid, row in fva_df.iterrows()
        }
    else:
        payload["fva"] = {}
    with open(out_path, "w") as f:
        json.dump(payload, f)
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", type=str, default=None,
                         help="Restrict to a single 'cohort:group' combo, e.g. "
                              "GSE89632:healthy -- for testing before a full run.")
    parser.add_argument("--no-fva", action="store_true",
                         help="Skip FVA (useful for a fast FBA/pFBA-only smoke test).")
    args = parser.parse_args()

    print("Loading Human-GEM base model...")
    base_model = load_base_model()
    print(f"  {len(base_model.reactions)} reactions, {len(base_model.metabolites)} metabolites, "
          f"objective = {base_model.objective.expression}")

    growth_support = load_growth_support()
    print(f"  Loaded {len(growth_support)}-reaction growth-support set "
          f"(kept unblocked alongside each group's strict consensus set).")

    consensus_manifest = pd.read_csv(os.path.join(CONSENSUS_DIR, "consensus_manifest.csv"))
    combos = consensus_manifest[["cohort", "group"]].drop_duplicates().values.tolist()

    if args.only:
        cohort_f, group_f = args.only.split(":")
        combos = [(c, g) for c, g in combos if c == cohort_f and g == group_f]
        if not combos:
            print(f"No matching combo for --only {args.only}")
            sys.exit(1)

    manifest_rows = []
    for cohort, group in combos:
        print(f"\n{cohort}/{group}:")
        result, fva_df, fba_sol = analyze_group(cohort, group, base_model, growth_support, run_fva=not args.no_fva)
        print(f"  consensus_definition_used={result['consensus_definition_used']} "
              f"n_strict_consensus={result['n_strict_consensus_reactions']} "
              f"n_growth_support_added={result['n_growth_support_added']} "
              f"n_active={result['n_active_reactions']} n_blocked={result['n_blocked_reactions']}")
        print(f"  FBA: status={result['fba_status']} objective_value={result['fba_objective_value']}")
        if result["pfba_status"] is not None:
            print(f"  pFBA: status={result['pfba_status']} "
                  f"objective_value={result['pfba_objective_value']} "
                  f"total_flux={result['pfba_total_flux']}")
            if fva_df is not None:
                print(f"  FVA: {result['n_fva_reactions']} reactions, "
                      f"range of minima [{fva_df['minimum'].min():.4g}, {fva_df['minimum'].max():.4g}], "
                      f"range of maxima [{fva_df['maximum'].min():.4g}, {fva_df['maximum'].max():.4g}]")

        out_path = write_group_output(cohort, group, result, fva_df, fba_sol)
        result["output"] = out_path
        manifest_rows.append(result)
        print(f"  -> {out_path}")

    manifest_df = pd.DataFrame(manifest_rows).drop(columns=["objective_reaction"], errors="ignore")
    manifest_path = os.path.join(OUT_DIR, "flux_manifest.csv")
    if args.only:
        print(f"\n(--only was used; not overwriting the full {manifest_path})")
    else:
        manifest_df.to_csv(manifest_path, index=False)
        print(f"\nManifest saved to {manifest_path}")


if __name__ == "__main__":
    main()
