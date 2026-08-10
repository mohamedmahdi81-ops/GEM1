"""
GEM1 - Step 5: Build context-specific models with all four reconstruction
algorithms (FASTCORE, GIMME, IMAT, tINIT) for each disease group in each
cohort. This is the reconstruction-algorithm-consensus half of the
FluxConfidence score (20% weight) -- the same disease group, extracted four
different ways, should largely agree if a candidate signature is real rather
than an artifact of the extraction method.

Gene-ID mapping (per cohort):
  - GSE89632 (Illumina microarray, GPL14951): raw expression is indexed by
    probe ID. Bridge probe -> Entrez Gene ID (GEO's own platform table) ->
    Human-GEM model gene ID (Human-GEM's `gene.annotation['ncbigene']`
    cross-reference). Where >1 probe maps to the same model gene, the probes
    are mean-aggregated (standard practice for multi-probe/gene microarray
    data) -- flagged here per Mohammed's standing instruction to note rather
    than pause on this kind of default.
  - GSE126848 / GSE135251 (RNA-seq): 02b_build_rnaseq_expression_matrices.py
    already produced CPM-normalized matrices indexed directly by Ensembl
    gene ID, which is Human-GEM's own gene namespace (confirmed directly:
    every `model.genes[i].id` is an ENSG id) -- no bridge needed, just
    restrict to genes the model actually has.

troppo API note: `ReconstructionWrapper.run_fastcore`/`run_gimme`/`run_imat`/
`run_tinit` referenced by an earlier draft of this script DO NOT EXIST in the
installed troppo 0.0.7 (verified via `dir(ReconstructionWrapper)` -- only
`run(properties)` and `run_from_omics(...)` exist). `run_from_omics`'s
integration-strategy path is unverified and troppo's own source carries a
"TODO change later, idk why core_idx is here" comment in
`FastcoreProperties.from_integrated_scores`, so this script instead builds
each algorithm's Properties object directly and calls `rw.run(properties)` --
the same low-level pattern already verified end-to-end for FASTCORE in
`_test_troppo_one_combo.py` (FastCC: 11641/12931 consistent, FASTCORE: 6300
reactions retained, submodel FBA optimal). GIMME/IMAT/tINIT's Properties
kwargs (`objectives` as a list of {reaction_index: coefficient} dicts for
GIMME, `exp_thresholds` as a (low, high) percentile pair for IMAT, centered
`reactions_scores` for tINIT) were derived by reading troppo's own algorithm
source (`optimize_gimme`, `IMAT.run_imat`), not from example docs -- these
are reasonable, standard defaults, not tuned values; re-check before trusting
outputs for anything beyond a first consensus pass.

Requires: cobra, troppo, GEOparse, pandas, numpy   (gem1-troppo env)
Run from gem1-troppo (needs troppo); Step 2 model + Step 4 cohorts must
already be in place, incl. 02b_build_rnaseq_expression_matrices.py's output.
"""

import os
import json
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

MODEL_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "models", "Human-GEM.json")
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "geo_cohorts")
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "context_specific_models")
os.makedirs(OUT_DIR, exist_ok=True)

COHORTS = ["GSE89632", "GSE126848", "GSE135251"]
ALGORITHMS = ["fastcore", "gimme", "imat", "tinit"]
BIOMASS_REACTION = "MAR13082"  # Human-GEM's "Generic human cell biomass reaction"

# 2026-07-17: Step 6 (FBA on the strict-consensus reaction sets) found iMAT
# and tINIT's extracted networks cannot sustain biomass flux -- their Step 5
# extraction never protected biomass or its supporting pathway, unlike
# FASTCORE (which happens to retain a growth-capable network) and GIMME
# (permissive enough that it does too). Forcing only the biomass reaction's
# own evidence score up is not sufficient -- the algorithm can mark biomass
# "active" while still excluding the upstream precursor-supply reactions it
# needs to actually carry flux. Fix: computed the actual minimal
# growth-supporting reaction set via a standalone LP (scipy.optimize.linprog,
# HiGHS) directly on Human-GEM.json's stoichiometry -- max-biomass FBA on the
# full unconstrained model reproduces the exact 124.8681 objective already
# recorded from the original Step 5 smoke test (confirms the LP is built
# correctly), then a parsimonious (L1-minimal) solution sustaining 90% of that
# (matching this script's existing obj_frac=0.9 convention for GIMME) has
# support on 519 reactions -- this is the set forced into iMAT/tINIT below,
# not just the single biomass reaction. See data/context_specific_models/
# minimal_growth_support.json for the full list and derivation.
#
# 2026-07-18 methodology check (before running the fix, not after): is
# forcing a growth-support set into the reconstruction itself defensible, or
# should this be a separate post-reconstruction gap-fill instead? Checked
# against the primary literature rather than assuming either way:
#   - Original INIT (Agren et al. 2012, PLOS Comp Biol) explicitly does NOT
#     use a biomass equation or guarantee growth capability -- not the
#     relevant precedent here.
#   - tINIT, the task-driven successor actually used for Human-GEM-family
#     models (Agren et al. 2014, Mol Syst Biol; see also the Human-GEM/
#     ftINIT guide), works in two stages: (1) identify reactions that are
#     hard-required for a predefined set of essential metabolic tasks
#     (including growth-related tasks) and force them into the model as
#     mandatory BEFORE/DURING the main optimization -- these "cannot be
#     removed during optimization" -- then (2) a post-hoc gap-fill pass only
#     for whichever individual tasks still fail after that. So forcing an
#     essential/growth-supporting reaction set in at reconstruction time is
#     not a deviation from tINIT's design -- it IS tINIT's design; running
#     tINIT without it (as the original Step 5 pass did) is what was
#     methodologically incomplete.
#   - troppo's own tINITProperties confirms this isn't an invented
#     workaround: it exposes a first-class `essential_reactions` parameter
#     (list of reaction indices) built exactly for this -- those reactions
#     get their lower bound floored to a nonzero minimum and are excluded
#     from the removable-candidate set, then unconditionally appended to the
#     result. Likewise troppo's IMATProperties exposes a `core` parameter
#     (reactions unioned directly into the "should be high/active" set,
#     independent of exp_vector/exp_thresholds) -- a purpose-built mechanism
#     for forcing specific reactions into iMAT's result, standard practice in
#     applied iMAT studies since plain iMAT (Zur/Shlomi 2010) has no built-in
#     growth guarantee either.
#   Conclusion: use these first-class parameters (see run_imat/run_tinit
#   below) rather than the score-inflation approach originally implemented --
#   inflating exp_vector/reactions_scores only biases the optimizer toward
#   keeping the growth-support set, it doesn't guarantee it, and it isn't
#   what either algorithm's own literature-documented mechanism does.
GROWTH_SUPPORT_PATH = os.path.join(OUT_DIR, "minimal_growth_support.json")
# 2026-07-18: tINIT's essential_reactions genuinely guarantees growth
# capability (validated 10/10 via scripts/_validate_growth_fix.py) and is
# unaffected by the iMAT core-parameter correction above, so it's excluded
# from this re-run to avoid re-solving 10 already-correct MILPs (~20-30 min).
# Restore to {"imat", "tinit"} if tINIT ever needs regenerating again.
FORCE_REGENERATE = {"imat", "tinit"}


def patch_gurobi_solve_limits(mip_gap=0.02, time_limit=180):
    """
    IMAT (troppo/methods/reconstruction/imat.py) sets no MIPGap/TimeLimit at
    all -- on Human-GEM-scale problems it ran >20 min without converging
    (see scripts/_test_all_algos_gurobi.log). tINIT (.../tINIT.py,
    solve_problem) hardcodes MIPGap=1e-3 with no TimeLimit, also too tight to
    converge quickly at this scale. Verified fix (scripts/_test_tinit_capped.log
    and the IMAT B&B trace in _test_all_algos_gurobi.log: first feasible
    incumbent already has a <1.4% gap by ~130s): MIPGap=2%, TimeLimit=180s
    converges both algorithms in 2-3 minutes with a stable solution.
    cobamp.core.optimization.LinearSystemOptimizer.optimize is the one choke
    point every algorithm's solve passes through (runs after tINIT's own
    hardcoded MIPGap=1e-3 line, so this overrides it), so patching it here
    avoids touching the installed troppo/cobamp source.

    Harder instances (e.g. GSE89632/NASH, GSE126848/NASH,
    GSE126848/obese_no_NAFLD) don't reach the 2% gap within 180s -- Gurobi
    stops at "Time limit reached" with a *feasible* incumbent (gap ~11%,
    SolCount > 0) but optlang's status for that isn't 'optimal', so tINIT's
    own `solution.status() != 'optimal'` check discards it as infeasible
    (`return` -> None), which then crashes downstream in
    `np.sort(None)` -> "axis -1 is out of bounds for array of dimension 0"
    (observed in scripts/03_run_full.log, first full-matrix run 2026-07-17).
    Fix: if Gurobi found >=1 feasible incumbent (SolCount > 0) even though it
    didn't prove the gap, treat it as usable -- rebuild the Solution with
    status forced to 'optimal' rather than raising. This is a looser bar than
    proven-optimal, consistent with this script's documented posture
    (reasonable defaults for a first consensus pass, not tuned/final values).
    """
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
                    solution = Solution(
                        solution.var_values(), "optimal",
                        objective_value=solution.objective_value(),
                    )
            except Exception:
                pass
        return solution

    LinearSystemOptimizer.optimize = patched_optimize


def patch_gimme_reaction_activity_bug():
    """
    troppo 0.0.7's GIMMESolution.get_reaction_activity (gimme.py:158) computes
    `ones = (np.array(self.exp_vector) > flux_threshold) | (self.exp_vector == -1)`.
    self.exp_vector is stored as a plain Python list (never cast to ndarray),
    so `self.exp_vector == -1` is a single Python `False` (list-to-scalar
    comparison), not an elementwise mask -- the documented rule "reactions
    with no usable expression score (-1) are always kept active" silently
    never fires, for every cohort. Confirmed via diagnostic (2026-07-17,
    CLAUDE.md gotcha): this alone doesn't explain the observed ~5x GIMME
    active-reaction gap between microarray and RNA-seq cohorts (the bug
    removes the same ~5,150 NaN-scored reactions from `ones` regardless of
    platform) -- see the min-max scaling fix in run_gimme() below for the
    other half of the interaction. Fixed here by replacing the method with a
    corrected version (np.array(self.exp_vector) cast before comparing)
    rather than editing installed troppo/cobamp source, consistent with
    patch_gurobi_solve_limits() above.
    """
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


def load_base_model():
    return cobra.io.load_json_model(MODEL_PATH)


def build_entrez_to_model_gene(model):
    mapping = {}
    for g in model.genes:
        entrez = g.annotation.get("ncbigene")
        if entrez:
            mapping[str(entrez)] = g.id
    return mapping


def load_gene_level_expression(gse_id, model, entrez_to_model_gene):
    """Returns a (model gene ID x GSM sample) expression DataFrame."""
    expr_path = os.path.join(DATA_DIR, f"{gse_id}_expression_matrix.csv")
    expr = pd.read_csv(expr_path, index_col=0)

    if gse_id == "GSE89632":
        gse = GEOparse.get_GEO(geo=gse_id, destdir=DATA_DIR, silent=True)
        gpl = gse.gpls["GPL14951"]
        probe_to_entrez = (
            gpl.table.set_index("ID")["Entrez_Gene_ID"].dropna().astype(int).astype(str).to_dict()
        )
        probe_to_gene = {
            probe: entrez_to_model_gene[entrez]
            for probe, entrez in probe_to_entrez.items()
            if entrez in entrez_to_model_gene
        }
        mapped = expr.loc[expr.index.intersection(probe_to_gene.keys())].copy()
        mapped.index = [probe_to_gene[p] for p in mapped.index]
        gene_expr = mapped.groupby(mapped.index).mean()
        print(f"  {gse_id}: mapped {len(probe_to_gene)} probes -> "
              f"{gene_expr.shape[0]} model genes (mean-aggregated where >1 probe/gene)")
    else:
        model_gene_ids = set(g.id for g in model.genes)
        gene_expr = expr.loc[expr.index.intersection(model_gene_ids)]
        print(f"  {gse_id}: {gene_expr.shape[0]} / {len(model_gene_ids)} model genes "
              f"found directly (already Ensembl-indexed, CPM-normalized)")

    return gene_expr


def load_cohort_meta(gse_id):
    meta_path = os.path.join(DATA_DIR, f"{gse_id}_sample_metadata.csv")
    meta = pd.read_csv(meta_path)
    if "disease_group" not in meta.columns:
        raise ValueError(
            f"{meta_path} has no 'disease_group' column yet. Open the CSV, "
            "add a 'disease_group' column (e.g. healthy / steatosis / NASH), "
            "and re-run."
        )
    return meta


def group_mean_expression(expr: pd.DataFrame, meta: pd.DataFrame, group: str) -> pd.Series:
    samples = meta.loc[meta["disease_group"] == group, "sample"]
    samples = [s for s in samples if s in expr.columns]
    if not samples:
        raise ValueError(f"No expression columns matched group '{group}'.")
    return expr[samples].mean(axis=1)


def compute_reaction_scores(rw, gene_scores: pd.Series) -> np.ndarray:
    scores = rw.model_reader.get_reaction_scores(
        gene_scores.to_dict(), and_fx=min, or_fx=max, as_vector=True
    )
    return np.array([np.nan if v is None else v for v in scores], dtype=float)


FASTCC_CACHE_PATH = os.path.join(OUT_DIR, "_fastcc_consistent_set.json")


def get_fastcc_consistent_set(rw) -> set:
    """
    Model-only, same every run -- cached to disk because this run has been
    killed mid-execution by the host process exiting multiple times
    (2026-07-17, see CLAUDE.md gotchas), and FastCC alone takes ~15-20 min on
    Human-GEM's 12,931 reactions, so re-deriving it on every restart wastes
    the bulk of the time between kills.
    """
    if os.path.exists(FASTCC_CACHE_PATH):
        with open(FASTCC_CACHE_PATH) as f:
            return set(json.load(f))
    fastcc_props = FastCCProperties(flux_threshold=1e-4, method="original", solver="GLPK")
    fastcc = FastCC(rw.S, rw.lb, rw.ub, fastcc_props)
    consistent_idx, *_ = fastcc.run()
    consistent_set = set(int(i) for i in consistent_idx)
    with open(FASTCC_CACHE_PATH, "w") as f:
        json.dump(sorted(consistent_set), f)
    return consistent_set


def run_fastcore(rw, reaction_scores, consistent_set):
    threshold = np.nanpercentile(reaction_scores, 75)
    core_idx = [i for i, v in enumerate(reaction_scores) if not np.isnan(v) and v >= threshold]
    core_idx = [i for i in core_idx if i in consistent_set]
    if not core_idx:
        raise ValueError("No flux-consistent core reactions above threshold for FASTCORE.")
    props = FastcoreProperties(core=core_idx, solver="GLPK")
    return list(rw.run(props))


def build_no_gpr_mask(model, r_ids) -> np.ndarray:
    """
    True for reactions with NO gene association at all (empty GPR -- transport,
    exchange, spontaneous reactions, a normal feature of genome-scale
    reconstructions: 5,149/12,931 in Human-GEM). These are the reactions
    GIMME's exp_vector==-1 rule was originally meant for (no expression
    evidence can ever apply to them, so they're exempt from the low-expression
    penalty). Distinct from reactions that DO have a GPR but whose score came
    back NaN because a mapped gene's expression wasn't measured/found (a
    mapping artifact, not biology) -- see run_gimme(). Assumes r_ids is in
    the same order as model.reactions (verified 2026-07-17 diagnostic).
    """
    return np.array([len(r.genes) == 0 for r in model.reactions])


def run_gimme(rw, reaction_scores, biomass_idx, no_gpr_mask):
    """
    Two fixes applied here 2026-07-17 (see CLAUDE.md gotcha for the full
    diagnostic and before/after numbers):

    1. Scaling: plain min-max scaling of raw reaction_scores is scale-
    sensitive. GSE89632 (microarray) scores are already log2 intensities
    (range ~7.4-15.8, compact). GSE126848/GSE135251 (RNA-seq) scores are
    linear CPM (range 0 to ~175,000-243,000, driven by a handful of extreme
    outlier genes e.g. albumin/apolipoproteins in liver tissue) -- plain
    min-max on that range crushed all but 2 of 12,931 RNA-seq reaction scores
    below the flux_threshold=0.25 cutoff, versus ~6,900 for microarray. log1p
    before min-max compresses that outlier-driven linear range while being a
    near-no-op in relative ordering for microarray's already-log2 values.
    Scoped to GIMME only -- FASTCORE/IMAT use percentile thresholds directly
    on reaction_scores (invariant under any monotonic transform), tINIT is
    untouched; the shared reaction_scores array itself is left raw here, only
    GIMME's local exp_vector is transformed.

    2. NaN-reaction handling: NaN reaction scores were previously all set to
    -1 (auto-active), conflating genuinely gene-less reactions (no_gpr_mask,
    5,149/12,931 -- legitimately exempt from expression evidence, matches
    GIMME's original intent) with the much smaller set of reactions that DO
    have a GPR but scored NaN anyway because a mapped gene's expression
    wasn't found (3-27 reactions per cohort, confirmed 2026-07-17 -- almost
    all NaN reactions turned out to be the former, not the latter). Only
    no_gpr_mask reactions get -1 now; GPR-present-but-unmapped reactions get
    0 (i.e. treated as the lowest possible expression -- penalized in GIMME's
    QP like any other low-expressed reaction, included only if the network
    needs their flux, same as everything else with weak/no evidence).
    """
    finite_mask = ~np.isnan(reaction_scores)
    finite = reaction_scores[finite_mask]
    log_finite = np.log1p(np.clip(finite, a_min=0, a_max=None))
    lo, hi = log_finite.min(), log_finite.max()
    scaled_finite = (log_finite - lo) / (hi - lo) if hi > lo else np.zeros_like(log_finite)
    scaled = np.zeros_like(reaction_scores)
    scaled[finite_mask] = scaled_finite
    nan_fallback = np.where(no_gpr_mask, -1.0, 0.0)
    exp_vector = np.where(finite_mask, scaled, nan_fallback)
    props = GIMMEProperties(
        exp_vector=list(exp_vector),
        objectives=[{biomass_idx: 1.0}],
        obj_frac=0.9,
        flux_threshold=0.25,
        solver="GUROBI",
        reaction_ids=rw.model_reader.r_ids,
        metabolite_ids=rw.model_reader.m_ids,
    )
    return list(rw.run(props))


def run_imat(rw, reaction_scores, growth_support_idx):
    """
    2026-07-18, corrected same day after validation caught it: IMATProperties'
    `core` parameter does NOT guarantee inclusion, contrary to what the
    methodology note above GROWTH_SUPPORT_PATH assumed before this was run.
    Verified directly against troppo's own source (methods/reconstruction/
    imat.py, run_imat(), line ~124-128): `core` reactions are only unioned
    into `high_idx`, the "should be active" classification -- the MILP then
    REWARDS (via binary Hpos/Hneg indicator variables in its objective) but
    does not REQUIRE those reactions to end up in `to_keep`. This is a soft
    preference, mechanically identical in strength to the score-inflation
    approach this function used before the 2026-07-18 methodology check, not
    the hard constraint that check concluded troppo provided. Confirmed
    empirically: the growth-fix validation (scripts/_validate_growth_fix.py)
    found all 10 iMAT groups extracted with `core` set to the growth-support
    set still returned 0.0 biomass flux when restricted to iMAT's own result
    (tINIT's `essential_reactions`, by contrast, passed 10/10 -- it really
    does floor lower bounds, exclude from the removable set, and
    unconditionally append, as documented).

    Fix: apply the same "real constraint, not a bias" principle the
    methodology note already endorsed, but the only way to get a real
    constraint out of iMAT without directional flux data (needed to safely
    floor individual reaction lower bounds -- not saved by the growth-support
    LP) is a post-hoc hard union of growth_support_idx into iMAT's own
    `to_keep`, mirroring tINIT's "unconditionally appended" semantics
    directly rather than routing through iMAT's classification machinery.
    This is guaranteed to restore growth capability by construction: the
    519-reaction growth-support set alone already sustains >=90% max biomass
    (see minimal_growth_support.json), and taking the union with iMAT's own
    result can only enlarge the feasible flux space (added reactions can
    always be left at zero flux if iMAT's own result already provides
    biomass), so the known-feasible growth-support solution remains feasible
    in the union.
    """
    finite = reaction_scores[~np.isnan(reaction_scores)]
    lo_thr, hi_thr = np.percentile(finite, [25, 75])
    exp_vector = np.where(np.isnan(reaction_scores), -1.0, reaction_scores)
    props = IMATProperties(
        exp_vector=exp_vector,
        exp_thresholds=(float(lo_thr), float(hi_thr)),
        core=list(int(i) for i in growth_support_idx),
    )
    idx = set(int(i) for i in rw.run(props))
    idx |= set(int(i) for i in growth_support_idx)
    return list(idx)


def run_tinit(rw, reaction_scores, growth_support_idx):
    """
    2026-07-18: uses tINITProperties' first-class `essential_reactions`
    parameter (verified against troppo's own source,
    methods/reconstruction/tINIT.py) rather than biasing reactions_scores --
    essential reactions have their lower bound floored to a nonzero minimum,
    are excluded from the removable-candidate set during preprocessing, and
    are unconditionally appended to the returned result. This mirrors the
    published tINIT algorithm's own two-stage design (hard-require reactions
    needed for essential metabolic tasks, then gap-fill only what's left) --
    see the methodology note above GROWTH_SUPPORT_PATH -- and is a stronger
    guarantee than the score-inflation this function used before this edit.
    """
    finite = reaction_scores[~np.isnan(reaction_scores)]
    center = np.median(finite)
    centered = np.where(np.isnan(reaction_scores), 0.0, reaction_scores - center)
    props = tINITProperties(
        reactions_scores=list(centered),
        essential_reactions=list(int(i) for i in growth_support_idx),
        solver="GUROBI",
    )
    idx = set(int(i) for i in rw.run(props))
    idx |= set(int(i) for i in growth_support_idx)
    return list(idx)


def extract_one(rw, reaction_scores, consistent_set, biomass_idx, algorithm, no_gpr_mask, growth_support_idx):
    if algorithm == "fastcore":
        idx = run_fastcore(rw, reaction_scores, consistent_set)
    elif algorithm == "gimme":
        idx = run_gimme(rw, reaction_scores, biomass_idx, no_gpr_mask)
    elif algorithm == "imat":
        idx = run_imat(rw, reaction_scores, growth_support_idx)
    elif algorithm == "tinit":
        idx = run_tinit(rw, reaction_scores, growth_support_idx)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    return [rw.model_reader.r_ids[i] for i in idx]


def main():
    patch_gurobi_solve_limits()
    patch_gimme_reaction_activity_bug()
    base_model = load_base_model()
    print(f"Base model loaded: {len(base_model.reactions)} reactions.")

    rw = ReconstructionWrapper(model=base_model, ttg_ratio=9999)
    biomass_idx = rw.model_reader.r_ids.index(BIOMASS_REACTION)
    no_gpr_mask = build_no_gpr_mask(base_model, rw.model_reader.r_ids)

    with open(GROWTH_SUPPORT_PATH) as f:
        growth_support_ids = set(json.load(f)["support_reaction_ids"])
    r_id_to_idx = {rid: i for i, rid in enumerate(rw.model_reader.r_ids)}
    growth_support_idx = np.array(sorted(r_id_to_idx[rid] for rid in growth_support_ids if rid in r_id_to_idx))
    print(f"Loaded minimal growth-supporting set: {len(growth_support_idx)}/{len(growth_support_ids)} "
          f"reactions matched in model_reader.r_ids (forced into iMAT/tINIT below).")

    print("Running FastCC once (model-only, reused across every cohort/group/algorithm)...")
    consistent_set = get_fastcc_consistent_set(rw)
    print(f"FastCC: {len(consistent_set)}/{len(rw.model_reader.r_ids)} reactions flux-consistent.")

    entrez_to_model_gene = build_entrez_to_model_gene(base_model)

    manifest = []
    for gse_id in COHORTS:
        try:
            expr = load_gene_level_expression(gse_id, base_model, entrez_to_model_gene)
            meta = load_cohort_meta(gse_id)
        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping {gse_id}: {e}")
            continue

        groups = sorted(meta["disease_group"].dropna().unique())
        print(f"\n{gse_id}: groups = {groups}")

        for group in groups:
            gene_scores = group_mean_expression(expr, meta, group)
            reaction_scores = compute_reaction_scores(rw, gene_scores)
            n_scored = int(np.sum(~np.isnan(reaction_scores)))
            print(f"  {group}: {n_scored}/{len(reaction_scores)} reactions scored")

            for algo in ALGORITHMS:
                out_path = os.path.join(OUT_DIR, f"{gse_id}__{group}__{algo}.json")
                if os.path.exists(out_path) and algo not in FORCE_REGENERATE:
                    with open(out_path) as f:
                        n_existing = len(json.load(f)["active_reactions"])
                    print(f"  Skipping {gse_id} / {group} / {algo} (already extracted, "
                          f"{n_existing} active reactions -> {out_path})")
                    manifest.append({
                        "cohort": gse_id, "group": group, "algorithm": algo,
                        "n_active_reactions": n_existing, "output": out_path,
                    })
                    continue
                action = "Re-extracting (growth-capability fix)" if os.path.exists(out_path) else "Extracting"
                print(f"  {action} {gse_id} / {group} / {algo} ...")
                try:
                    result = extract_one(rw, reaction_scores, consistent_set, biomass_idx, algo, no_gpr_mask, growth_support_idx)
                    with open(out_path, "w") as f:
                        json.dump({"active_reactions": list(result)}, f)
                    manifest.append({
                        "cohort": gse_id, "group": group, "algorithm": algo,
                        "n_active_reactions": len(result), "output": out_path,
                    })
                    print(f"    -> {len(result)} active reactions -> {out_path}")
                except Exception as e:
                    print(f"    FAILED ({algo} on {gse_id}/{group}): {e}")
                    manifest.append({
                        "cohort": gse_id, "group": group, "algorithm": algo,
                        "n_active_reactions": None, "output": None, "error": str(e),
                    })

    manifest_path = os.path.join(OUT_DIR, "extraction_manifest.csv")
    pd.DataFrame(manifest).to_csv(manifest_path, index=False)
    print(f"\nManifest of all extracted models saved to {manifest_path}")
    print("Next: Step 9 (reconstruction-algorithm consensus scoring) compares "
          "the active-reaction sets across the 4 algorithms per cohort/group.")


if __name__ == "__main__":
    main()
