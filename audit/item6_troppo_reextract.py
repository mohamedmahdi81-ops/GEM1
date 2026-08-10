"""
Item 6 -- Independently reconstruct active-reaction sets for all 39/40
cohort x group x algorithm combinations directly from raw expression data
(GEO expression matrices) and Human-GEM.json via troppo, and compare
active-reaction counts against data/context_specific_models/extraction_manifest.csv.

Does NOT import scripts/03_build_context_specific_models.py. This is a fresh
implementation, written independently against troppo's own API (methods
directly verified from troppo's source during original pipeline development,
documented in that script's docstrings -- reused here as "this is what the
algorithms are specified to do", not as borrowed scoring/calibration logic;
03_build_context_specific_models.py is a reconstruction/extraction script,
not one of the scoring/calibration modules (04/05/09/10/11) this audit is
required to avoid).

Reuses the *documented protocol parameters* (percentile thresholds, GIMME
log1p scaling + flux_threshold=0.25, IMAT exp_thresholds=(25,75) percentile,
tINIT median-centering, growth-support-set union into iMAT/tINIT, Gurobi
MIPGap=0.02/TimeLimit=180 caps) because those define WHAT Step 5 is specified
to compute -- an independence check has to run the same specified procedure
on the same raw inputs, or it isn't testing reproducibility at all.

FastCC (needed for FASTCORE's core-reaction filtering) is recomputed fresh
here rather than reusing data/context_specific_models/_fastcc_consistent_set.json,
for maximum independence (model-only, ~15-20 min one-time cost).

MUST run under the gem1-troppo conda env (Python 3.10, has troppo/cobamp/
cobra 0.24/optlang 1.8.3/GEOparse/gurobipy) -- NOT the GEM1/gem1-main env.

Long-running (~70-90 min including fresh FastCC) -- launch via detached
Start-Process per the project's established pattern.
"""

import json
import os
import time

import numpy as np
import pandas as pd
import cobra
import GEOparse

from troppo.methods_wrappers import ReconstructionWrapper
from troppo.methods.reconstruction.fastcore import FastcoreProperties
from troppo.methods.reconstruction.gimme import GIMMEProperties
from troppo.methods.reconstruction.imat import IMATProperties
from troppo.methods.reconstruction.tINIT import tINITProperties
from troppo.methods.gapfill.fastcc import FastCC, FastCCProperties

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(REPO_ROOT, "models", "Human-GEM.json")
GEO_DIR = os.path.join(REPO_ROOT, "data", "geo_cohorts")
CONTEXT_MODELS_DIR = os.path.join(REPO_ROOT, "data", "context_specific_models")
GROWTH_SUPPORT_PATH = os.path.join(CONTEXT_MODELS_DIR, "minimal_growth_support.json")
AUDIT_DIR = os.path.join(REPO_ROOT, "audit")
DETAILS_DIR = os.path.join(AUDIT_DIR, "details")
RAW_DIR = os.path.join(DETAILS_DIR, "item6_raw")
os.makedirs(RAW_DIR, exist_ok=True)

COHORTS = ["GSE89632", "GSE126848", "GSE135251"]
ALGORITHMS = ["fastcore", "gimme", "imat", "tinit"]
BIOMASS_REACTION = "MAR13082"


def patch_gurobi_solve_limits(mip_gap=0.02, time_limit=180):
    from cobamp.core.optimization import LinearSystemOptimizer, Solution
    original_optimize = LinearSystemOptimizer.optimize

    def patched_optimize(self):
        is_gurobi = getattr(self, "solver", None) == "GUROBI"
        if is_gurobi:
            try:
                self.model.problem.Params.MIPGap = mip_gap
                self.model.problem.Params.TimeLimit = time_limit
            except Exception:
                pass
        solution = original_optimize(self)
        if is_gurobi and solution.status() != "optimal":
            try:
                if self.model.problem.SolCount > 0:
                    solution = Solution(solution.var_values(), "optimal", objective_value=solution.objective_value())
            except Exception:
                pass
        return solution

    LinearSystemOptimizer.optimize = patched_optimize


def patch_gimme_reaction_activity_bug():
    from troppo.methods.reconstruction.gimme import GIMMESolution

    def fixed_get_reaction_activity(self, flux_threshold):
        gimme_fluxes = np.array([kv[1] for kv in self.var_values().items()])
        exp_arr = np.array(self.exp_vector, dtype=float)
        activity = np.zeros(gimme_fluxes.shape)
        ones = (exp_arr > flux_threshold) | (exp_arr == -1)
        twos = gimme_fluxes > 0
        activity[ones] = 1
        activity[twos & ~ones] = 2
        return [idx for idx, val in enumerate(activity) if val != 0]

    GIMMESolution.get_reaction_activity = fixed_get_reaction_activity


def build_entrez_to_model_gene(model):
    mapping = {}
    for g in model.genes:
        entrez = g.annotation.get("ncbigene")
        if entrez:
            mapping[str(entrez)] = g.id
    return mapping


def load_gene_level_expression(gse_id, model, entrez_to_model_gene):
    expr_path = os.path.join(GEO_DIR, f"{gse_id}_expression_matrix.csv")
    expr = pd.read_csv(expr_path, index_col=0)
    if gse_id == "GSE89632":
        gse = GEOparse.get_GEO(geo=gse_id, destdir=GEO_DIR, silent=True)
        gpl = gse.gpls["GPL14951"]
        probe_to_entrez = gpl.table.set_index("ID")["Entrez_Gene_ID"].dropna().astype(int).astype(str).to_dict()
        probe_to_gene = {p: entrez_to_model_gene[e] for p, e in probe_to_entrez.items() if e in entrez_to_model_gene}
        mapped = expr.loc[expr.index.intersection(probe_to_gene.keys())].copy()
        mapped.index = [probe_to_gene[p] for p in mapped.index]
        gene_expr = mapped.groupby(mapped.index).mean()
    else:
        model_gene_ids = set(g.id for g in model.genes)
        gene_expr = expr.loc[expr.index.intersection(model_gene_ids)]
    return gene_expr


def load_cohort_meta(gse_id):
    meta = pd.read_csv(os.path.join(GEO_DIR, f"{gse_id}_sample_metadata.csv"))
    return meta


def group_mean_expression(expr, meta, group):
    samples = meta.loc[meta["disease_group"] == group, "sample"]
    samples = [s for s in samples if s in expr.columns]
    return expr[samples].mean(axis=1)


def compute_reaction_scores(rw, gene_scores):
    scores = rw.model_reader.get_reaction_scores(gene_scores.to_dict(), and_fx=min, or_fx=max, as_vector=True)
    return np.array([np.nan if v is None else v for v in scores], dtype=float)


def get_fastcc_consistent_set_fresh(rw):
    fastcc_props = FastCCProperties(flux_threshold=1e-4, method="original", solver="GLPK")
    fastcc = FastCC(rw.S, rw.lb, rw.ub, fastcc_props)
    consistent_idx, *_ = fastcc.run()
    return set(int(i) for i in consistent_idx)


def build_no_gpr_mask(model):
    return np.array([len(r.genes) == 0 for r in model.reactions])


def run_fastcore(rw, reaction_scores, consistent_set):
    threshold = np.nanpercentile(reaction_scores, 75)
    core_idx = [i for i, v in enumerate(reaction_scores) if not np.isnan(v) and v >= threshold]
    core_idx = [i for i in core_idx if i in consistent_set]
    if not core_idx:
        raise ValueError("No flux-consistent core reactions above threshold for FASTCORE.")
    return list(rw.run(FastcoreProperties(core=core_idx, solver="GLPK")))


def run_gimme(rw, reaction_scores, biomass_idx, no_gpr_mask):
    finite_mask = ~np.isnan(reaction_scores)
    finite = reaction_scores[finite_mask]
    log_finite = np.log1p(np.clip(finite, a_min=0, a_max=None))
    lo, hi = log_finite.min(), log_finite.max()
    scaled_finite = (log_finite - lo) / (hi - lo) if hi > lo else np.zeros_like(log_finite)
    scaled = np.zeros_like(reaction_scores)
    scaled[finite_mask] = scaled_finite
    nan_fallback = np.where(no_gpr_mask, -1.0, 0.0)
    exp_vector = np.where(finite_mask, scaled, nan_fallback)
    props = GIMMEProperties(exp_vector=list(exp_vector), objectives=[{biomass_idx: 1.0}], obj_frac=0.9,
                             flux_threshold=0.25, solver="GUROBI",
                             reaction_ids=rw.model_reader.r_ids, metabolite_ids=rw.model_reader.m_ids)
    return list(rw.run(props))


def run_imat(rw, reaction_scores, growth_support_idx):
    finite = reaction_scores[~np.isnan(reaction_scores)]
    lo_thr, hi_thr = np.percentile(finite, [25, 75])
    exp_vector = np.where(np.isnan(reaction_scores), -1.0, reaction_scores)
    props = IMATProperties(exp_vector=exp_vector, exp_thresholds=(float(lo_thr), float(hi_thr)),
                            core=list(int(i) for i in growth_support_idx))
    idx = set(int(i) for i in rw.run(props)) | set(int(i) for i in growth_support_idx)
    return list(idx)


def run_tinit(rw, reaction_scores, growth_support_idx):
    finite = reaction_scores[~np.isnan(reaction_scores)]
    center = np.median(finite)
    centered = np.where(np.isnan(reaction_scores), 0.0, reaction_scores - center)
    props = tINITProperties(reactions_scores=list(centered), essential_reactions=list(int(i) for i in growth_support_idx), solver="GUROBI")
    idx = set(int(i) for i in rw.run(props)) | set(int(i) for i in growth_support_idx)
    return list(idx)


def run():
    t0 = time.time()
    patch_gurobi_solve_limits()
    patch_gimme_reaction_activity_bug()

    base_model = cobra.io.load_json_model(MODEL_PATH)
    print(f"Base model: {len(base_model.reactions)} reactions.", flush=True)

    rw = ReconstructionWrapper(model=base_model, ttg_ratio=9999)
    biomass_idx = rw.model_reader.r_ids.index(BIOMASS_REACTION)
    no_gpr_mask = build_no_gpr_mask(base_model)

    with open(GROWTH_SUPPORT_PATH) as f:
        growth_support_ids = set(json.load(f)["support_reaction_ids"])
    r_id_to_idx = {rid: i for i, rid in enumerate(rw.model_reader.r_ids)}
    growth_support_idx = np.array(sorted(r_id_to_idx[rid] for rid in growth_support_ids if rid in r_id_to_idx))

    print("Running fresh FastCC (independent of the pipeline's cached _fastcc_consistent_set.json)...", flush=True)
    tcc0 = time.time()
    consistent_set = get_fastcc_consistent_set_fresh(rw)
    print(f"FastCC: {len(consistent_set)}/{len(rw.model_reader.r_ids)} consistent, {time.time()-tcc0:.0f}s", flush=True)

    entrez_to_model_gene = build_entrez_to_model_gene(base_model)

    existing_manifest = pd.read_csv(os.path.join(CONTEXT_MODELS_DIR, "extraction_manifest.csv"))

    results = []
    for gse_id in COHORTS:
        expr = load_gene_level_expression(gse_id, base_model, entrez_to_model_gene)
        meta = load_cohort_meta(gse_id)
        groups = sorted(meta["disease_group"].dropna().unique())
        for group in groups:
            gene_scores = group_mean_expression(expr, meta, group)
            reaction_scores = compute_reaction_scores(rw, gene_scores)
            for algo in ALGORITHMS:
                tA = time.time()
                try:
                    if algo == "fastcore":
                        idx = run_fastcore(rw, reaction_scores, consistent_set)
                    elif algo == "gimme":
                        idx = run_gimme(rw, reaction_scores, biomass_idx, no_gpr_mask)
                    elif algo == "imat":
                        idx = run_imat(rw, reaction_scores, growth_support_idx)
                    elif algo == "tinit":
                        idx = run_tinit(rw, reaction_scores, growth_support_idx)
                    active_ids = [rw.model_reader.r_ids[i] for i in idx]
                    n_active = len(active_ids)
                    error = None
                except Exception as e:
                    active_ids, n_active, error = None, None, str(e)

                elapsed = time.time() - tA
                existing_row = existing_manifest[
                    (existing_manifest["cohort"] == gse_id) & (existing_manifest["group"] == group) &
                    (existing_manifest["algorithm"] == algo)
                ]
                existing_n = existing_row["n_active_reactions"].iloc[0] if len(existing_row) else None
                existing_n = None if pd.isna(existing_n) else existing_n

                row = {
                    "cohort": gse_id, "group": group, "algorithm": algo,
                    "independent_n_active": n_active, "existing_n_active": existing_n,
                    "abs_diff": (abs(n_active - existing_n) if (n_active is not None and existing_n is not None) else None),
                    "rel_diff": (abs(n_active - existing_n) / existing_n if (n_active is not None and existing_n and existing_n > 0) else None),
                    "error": error, "elapsed_seconds": elapsed,
                }
                results.append(row)
                print(f"  {gse_id}/{group}/{algo}: independent={n_active} existing={existing_n} ({elapsed:.0f}s)", flush=True)

                if active_ids is not None:
                    with open(os.path.join(RAW_DIR, f"{gse_id}__{group}__{algo}.json"), "w") as f:
                        json.dump({"active_reactions": active_ids}, f)

    df = pd.DataFrame(results)
    df.to_csv(os.path.join(DETAILS_DIR, "item6_comparison.csv"), index=False)

    # Tolerance: FASTCORE is deterministic LP given same threshold/consistent-set -> expect exact
    # or near-exact match. iMAT/tINIT are MILP under time/gap caps -> allow relative tolerance.
    # GIMME is QP-based, near-deterministic given fixed exp_vector -> tight tolerance.
    TOLERANCES = {"fastcore": 0.02, "gimme": 0.05, "imat": 0.15, "tinit": 0.15}
    df["within_tolerance"] = df.apply(
        lambda r: (r["rel_diff"] is not None and r["rel_diff"] <= TOLERANCES.get(r["algorithm"], 0.15))
        or (r["error"] is not None and pd.isna(r["existing_n_active"])),
        axis=1,
    )
    n_expected = 39  # 40 combos - 1 documented legitimate FASTCORE failure (GSE126848/healthy)
    n_reconstructed = int((df["independent_n_active"].notna()).sum())
    n_within_tol = int(df["within_tolerance"].sum())

    max_abs_discrepancy = float(df["abs_diff"].dropna().max()) if df["abs_diff"].notna().any() else None

    summary = {
        "item": 6,
        "description": "Active-reaction counts for all cohort x group x algorithm combinations, independently reconstructed via fresh troppo code",
        "n_combinations_total": len(df),
        "n_combinations_reconstructed": n_reconstructed,
        "n_combinations_expected_success": n_expected,
        "n_within_tolerance": n_within_tol,
        "tolerances_used": TOLERANCES,
        "fastcc_consistent_reactions": len(consistent_set),
        "max_abs_discrepancy": max_abs_discrepancy,
        "status": "PASS" if n_within_tol >= n_expected else "FAIL",
        "elapsed_seconds": time.time() - t0,
    }
    with open(os.path.join(DETAILS_DIR, "item6_troppo_reextract.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nStatus: {summary['status']}, {n_within_tol}/{len(df)} within tolerance, total elapsed {summary['elapsed_seconds']:.0f}s", flush=True)
    return summary


if __name__ == "__main__":
    run()
