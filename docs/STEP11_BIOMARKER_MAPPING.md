# GEM1 — Step 11 Calibration Biomarker → Reaction Mapping

Status: **DRAFT, for review. Not yet used in any code.** This document maps the 5 calibration
biomarkers locked in under "Decisions locked in" in `CLAUDE.md` to specific Human-GEM
(`models/Human-GEM.json`) reaction IDs, so Step 11's empirical weight calibration has a fixed,
version-controlled, auditable ground truth to calibrate against — rather than reaction choices
made ad hoc inside calibration code.

All reaction IDs, names, compartments, and gene-reaction-rules (GPR, Ensembl gene IDs) below were
read directly from `models/Human-GEM.json` (searched by subsystem, reaction name, and metabolite
name — not guessed or recalled from memory). Gene **symbols** for the items originally flagged as
uncertain (BCKDH complex subunits, CHKA/CHKB, PCYT1A/1B) have since been verified via a live
MyGene.info lookup (2026-07-21) — see the "Remaining uncertainties" section for the confirmed
mapping. The FASN reaction choice and the SCD1 triplicate have also been resolved (biomarker 4,
below) after tracing the actual reaction chemistry in the model file.

No Step 11 code has been written. This document is the deliverable; implementation is a separate,
later step pending review of everything below.

---

## 1. Serine/glycine depletion (one-carbon metabolism)

**Biological rationale**: Serine and glycine are interconverted by serine hydroxymethyltransferase
(SHMT), which also feeds one-carbon units (as 5,10-methylene-THF) into folate metabolism, purine
synthesis, and the methionine/SAM cycle. Reduced serine/glycine availability is proposed to reflect
increased diversion of one-carbon units to support the proliferative/lipogenic demands of a
steatotic/NASH liver.

**Literature rationale** (per `CLAUDE.md`'s "Decisions locked in"): Mardinoglu et al. 2014 (GSMM
study) + an independent 2025 metabolomics vote-counting meta-analysis + a 2023 SHMT2 mechanistic
study specifically.

**Reactions selected**:

| Reaction ID | Name | Compartment | Gene (Ensembl) | Gene symbol (best-confidence, unverified) |
|---|---|---|---|---|
| `MAR03845` | 5,10-Methylenetetrahydrofolate:glycine hydroxymethyltransferase | cytosol (`c`) | ENSG00000176974 | SHMT1 |
| `MAR04792` | 5,10-Methylenetetrahydrofolate:glycine hydroxymethyltransferase | mitochondria (`m`) | ENSG00000182199 | SHMT2 |

**Why selected**: These are the only two reactions in the "Glycine, serine and threonine
metabolism" subsystem whose name and metabolites (serine, glycine, 5,10-methylene-THF, H2O) exactly
match the canonical SHMT reaction. They are compartment-specific isoforms of the same
interconversion, not alternative pathways.

**Both compartments included, not just SHMT2**: even though the specific literature citation above
is about SHMT2 (mitochondrial), SHMT1 (cytosolic) catalyzes the same reaction and is also
NAFLD-relevant in the broader one-carbon-metabolism literature. Including both captures pathway-level
impairment rather than assuming the mitochondrial isoform alone carries the whole signal; if a
single-reaction proxy is preferred instead, SHMT2 (`MAR04792`) is the better-justified choice given
the specific citation.

**Alternatives considered and rejected**: `MAR03901` (glycine N-methyltransferase, GNMT) and
`MAR04582` (glycine amidinotransferase) also involve glycine and one-carbon/methyl-donor chemistry,
but neither directly interconverts serine and glycine — GNMT consumes glycine as a SAM-methyl sink
(regulates SAM:SAH ratio, a different axis), and the amidinotransferase reaction is about creatine
synthesis. Excluded as off-target for this specific biomarker.

---

## 2. Ureagenesis / arginine-cycle impairment

**Biological rationale**: The urea cycle disposes of ammonia via a 5-enzyme pathway
(CPS1 → OTC → ASS1 → ASL → ARG1/ARG2). Impaired ureagenesis is proposed as a NAFLD/NASH marker
reflecting broader hepatocyte metabolic dysfunction, not a single-enzyme defect.

**Literature rationale**: Mardinoglu et al. 2016 (GSMM) + the same 2025 meta-analysis, with arginine
independently reported down.

**Note on Human-GEM's own "Urea cycle" subsystem label**: Human-GEM's subsystem literally named
"Urea cycle" contains exactly one reaction (`MAR02527`, polyamine-related, not a urea cycle enzyme)
— a model/annotation quirk, not a missing pathway. All 5 real urea cycle enzymes are present in the
model but filed under "Arginine and proline metabolism" (4 of them) and "Alanine, aspartate and
glutamate metabolism" (CPS1). Found by searching reaction names/metabolites directly rather than
trusting the subsystem label.

**Reactions selected** (all 5 canonical urea cycle steps):

| Step | Reaction ID | Name | Compartment | Gene (Ensembl) | Gene symbol (best-confidence, unverified) |
|---|---|---|---|---|---|
| 1. CPS1 | `MAR03873` | Carbon-dioxide:ammonia ligase (ADP-forming, carbamate-phosphorylating) | mitochondria (`m`) | ENSG00000021826 | CPS1 |
| 2. OTC | `MAR03809` | Carbamoyl-phosphate:L-ornithine carbamoyltransferase | mitochondria (`m`) | ENSG00000036473 | OTC |
| 3. ASS1 | `MAR03811` | L-Citrulline:L-aspartate ligase (AMP-forming) | cytosol (`c`) | ENSG00000130707 | ASS1 |
| 4. ASL | `MAR03813` | 2-(N-omega-L-arginino)succinate arginine-lyase (fumarate-forming) | cytosol (`c`) | ENSG00000126522 | ASL |
| 5a. ARG1 | `MAR03816` | L-Arginine amidinohydrolase | cytosol (`c`) | ENSG00000118520 | ARG1 |
| 5b. ARG2 | `MAR08426` | L-Arginine amidinohydrolase | mitochondria (`m`) | ENSG00000081181 | ARG2 |

**Why selected**: Each reaction's substrate/product metabolite list matches the textbook urea cycle
step exactly (e.g. `MAR03873`: NH3 + CO2 + 2 ATP → carbamoyl-phosphate; `MAR03816`/`MAR08426`:
arginine + H2O → ornithine + urea). CPS1 was the one step not found in this session's first pass
(reported to Mohammed as "found OTC and arginase, still need CPS1/ASS1/ASL") — a second, more
targeted search (reactions containing the metabolite "carbamoyl-phosphate" as a product, filtered to
the mitochondrial, ammonia-consuming variant rather than the cytosolic glutamine-consuming
pyrimidine-synthesis variant) found it at `MAR03873`. ASS1/ASL were found the same way, searching
for "argininosuccinate" as a metabolite rather than by reaction name.

**Alternatives considered and rejected**: `MAR04034` ("HCO3-:L-glutamine amido-ligase
(ADP-forming, carbamate-...)") also produces carbamoyl-phosphate but is CPS2/CAD, the cytosolic
enzyme feeding pyrimidine (not urea) synthesis — wrong pathway despite the similar product. Excluded.
`MAR09919` ("Exchange of argininosuccinate") and `MAR09911` ("Transport of L-Arginosuccinic Acid")
are boundary/transport reactions, not the synthesis/lysis steps — excluded as not representing
enzyme activity.

**ARG1 vs. ARG2, both included**: ARG1 (cytosolic, `MAR03816`) is the liver-predominant,
canonical urea-cycle arginase; ARG2 (mitochondrial, `MAR08426`) is generally considered more of an
extrahepatic/broadly-expressed isoform in the literature. Both are included here for completeness
and because Human-GEM models them as distinct reactions, but if a single-reaction proxy is wanted,
ARG1 is the better-justified choice for a liver-focused NAFLD model.

---

## 3. Elevated BCAAs (valine, leucine, isoleucine)

**Biological rationale**: Elevated circulating branched-chain amino acids are a well-established
NAFLD/insulin-resistance marker, reflecting impaired BCAA catabolism (not simply increased intake).
The pathway has two stages for all three BCAAs: (1) a shared-specificity transaminase (BCAT) step,
then (2) the rate-limiting, shared branched-chain ketoacid dehydrogenase (BCKDH) complex.

**Literature rationale**: Multiple independent reviews/cohorts, plus the same 2025 meta-analysis.

**Reactions selected**:

Transaminase (BCAT) step — one gene per compartment handles all three BCAAs (not three separate
enzymes per amino acid):

| Reaction ID | Amino acid | Compartment | Gene (Ensembl) | Gene symbol (best-confidence, unverified) |
|---|---|---|---|---|
| `MAR03744` | valine | mitochondria (`m`) | ENSG00000105552 | BCAT2 |
| `MAR03747` | valine | cytosol (`c`) | ENSG00000060982 | BCAT1 |
| `MAR03765` | leucine | mitochondria (`m`) | ENSG00000105552 | BCAT2 |
| `MAR06923` | leucine | cytosol (`c`) | ENSG00000060982 | BCAT1 |
| `MAR03777` | isoleucine | cytosol (`c`) | ENSG00000060982 | BCAT1 |
| `MAR03778` | isoleucine | mitochondria (`m`) | ENSG00000105552 | BCAT2 |

Branched-chain ketoacid dehydrogenase (BCKDH) complex — rate-limiting step, one multi-subunit
complex handles all three ketoacid substrates:

| Reaction ID | Substrate (from) | Compartment | Gene (Ensembl, 3-subunit complex) |
|---|---|---|---|
| `MAR06416` | valine-derived ketoacid | mitochondria (`m`) | ENSG00000083123 and ENSG00000137992 and ENSG00000248098 |
| `MAR06419` | isoleucine-derived ketoacid | mitochondria (`m`) | ENSG00000083123 and ENSG00000137992 and ENSG00000248098 |
| `MAR06421` | leucine-derived ketoacid | mitochondria (`m`) | ENSG00000083123 and ENSG00000137992 and ENSG00000248098 |

**Gene symbols for the BCKDH complex, confirmed via MyGene.info**: `ENSG00000083123` = `BCKDHB`
(E1-beta), `ENSG00000137992` = `DBT` (E2), `ENSG00000248098` = `BCKDHA` (E1-alpha) — the 3-gene AND
rule is exactly the expected E1-alpha/E1-beta/E2 subunit composition. DLD (E3) is not part of the
rule, which is expected since it's shared across BCKDH/PDH/OGDH and not BCKDH-specific.

**Why selected, and why both BCAT and BCKDH rather than just one**: user (Mohammed) explicitly noted
that including both the transaminase and the dehydrogenase step is stronger than picking one,
because elevated BCAAs in NAFLD reflect pathway-level impairment, not a single-enzyme lesion — a
weight-calibration signal built from only the first step would miss impairment specifically at the
rate-limiting BCKDH step (the step most commonly implicated as the actual bottleneck in the human
BCAA/insulin-resistance literature). Both cytosolic (BCAT1) and mitochondrial (BCAT2) transaminase
reactions are included for the same completeness reasoning as the SHMT and arginase cases above.

**Alternatives considered and rejected**: none found that were a closer match — the transaminase and
BCKDH steps are the two literature-standard control points for BCAA catabolism, so no other reaction
in the "Valine, leucine, and isoleucine metabolism" subsystem was a plausible alternative for this
specific biomarker (the subsystem also contains further downstream degradation reactions toward
acetyl-CoA/propionyl-CoA, which are further from the "elevated BCAA" signal itself and were not
included).

---

## 4. Increased de novo lipogenesis (SCD1/FASN flux)

**Biological rationale**: NAFLD is characterized by increased hepatic de novo lipogenesis — new
fatty acid synthesis from acetyl-CoA/malonyl-CoA (FASN) and subsequent desaturation to
monounsaturated fatty acids (SCD1) — contributing to hepatic triglyceride accumulation.

**Literature rationale**: Multiple independent stable-isotope tracer studies (per `CLAUDE.md`).

**SCD1 — reaction selected**:

| Reaction ID | Name | Product | Compartment | Gene (Ensembl) | Gene symbol |
|---|---|---|---|---|---|
| `MAR00146` | Stearoyl Coenzyme A Desaturase (C18:0-CoA → C18:1-CoA) | oleoyl-CoA (Δ9-cis) | cytosol (`c`) | ENSG00000099194 | SCD |

**Resolved**: Human-GEM contains two further reactions with identical gene/substrate/cofactors —
`MAR00147` (product: "11-Octadecenoyl Coenzyme A", a Δ11 positional isomer) and `MAR00148`
(product: "(2E)-octadecenoyl-CoA", a Δ2 geometric/positional isomer). These are not duplicate
encodings of the same reaction as first thought — they are SCD acting on the same substrate but
producing three different C18:1-CoA isomers. `MAR00146`'s product, oleoyl-CoA (Δ9-cis, oleic acid),
is the canonical SCD1 product referenced throughout the NAFLD/DNL literature (the standard
"desaturation index" is oleate/stearate or palmitoleate/palmitate) — `MAR00147`/`MAR00148` are
excluded as alternate positional isomers, not the isomer this biomarker refers to. Since all three
share an identical single-gene GPR, this choice doesn't change what the confidence engine would
score (they'd get identical expression-driven activity calls) — it only matters for biological
accuracy of the biomarker identity, which favors `MAR00146` unambiguously.

**FASN — resolved**:

Tracing the actual model chemistry (metabolite names/stoichiometry, not just gene annotations)
resolves this differently than either originally-considered candidate:

- `MAR01471` ("B-Ketoacyl Synthetase (Palmitate, N-C16:0)": acetyl-CoA + 7 malonyl-CoA + 14 NADPH →
  palmitate + 8 CoA + 7 CO2 + 6 H2O + 14 NADP+) has the textbook-exact net stoichiometry but a
  genuinely blank `gene_reaction_rule` — a lumped pseudo-reaction bypassing the ACP-linked
  mechanism. Beyond lacking a symbol, this is a structural problem for calibration: with no GPR it
  cannot receive an expression-driven activity call in Step 5's troppo extraction (every algorithm
  keys off gene expression via the GPR), so it can't respond to differential FASN expression across
  NAFLD/healthy/NASH groups. **Excluded** as unusable for this purpose despite the clean biology.
- `MAR02323` / `MAR02325` / `MAR02327` ("Fatty Acid Synthase Polyunsaturated Fatty Acid
  Biosynthesis") do carry `ENSG00000169710` (FASN), but tracing their metabolites shows they are
  **not** the primary palmitate-producing pathway: they form the linear chain hexadecanoyl-ACP →
  3-oxostearoyl-ACP (`MAR02322`) → 3-hydroxystearoyl-ACP (`MAR02323`) → octadecenoyl-ACP
  (`MAR02325`) → stearoyl-ACP (`MAR02327`), i.e. **C16→C18 chain elongation**, not palmitate
  synthesis. The chain's first (committing) step, `MAR02322`, is annotated to `ENSG00000151093` =
  **OXSM** ("3-oxoacyl-ACP synthase, mitochondrial") — a distinct mitochondrial-type
  fatty-acid-synthesis enzyme, not FASN — even though everything here is modeled in the cytosolic
  compartment (confirming the earlier-flagged "Mitochondrial" name/compartment-label
  inconsistency is a real annotation quirk, inherited from the model's source nomenclature, not a
  faithful representation of the core cytosolic FASN cycle). **Excluded** as an off-target match,
  the same way GNMT/creatine reactions were excluded from the SHMT biomarker above.
- **The correct reaction is `MAR02182`** ("hexadecanoyl-[acyl-carrier protein] hydrolase":
  hexadecanoyl-ACP + H2O → palmitate + ACP + H+), the terminal, product-releasing step of a
  complete, mechanistically faithful, gene-annotated FASN cycle in subsystem "Fatty acid
  biosynthesis (even-chain)" (`MAR02150` through `MAR02182`, ~33 reactions: sequential
  loading/condensation/reduction/dehydration/reduction steps building the acyl chain from
  acetyl-CoA up to hexadecanoyl-ACP). `MAR02182`'s GPR is `ENSG00000152463 (OLAH) or
  ENSG00000169710 (FASN)` — it simultaneously matches the correct biology (produces palmitate, the
  literature's DNL readout) and carries a real, expression-responsive FASN GPR. `MAR02150` (first
  committed step, acetyl-CoA loading onto ACP, GPR = FASN alone) is also included for pathway-entry
  coverage, matching the "first step + rate-limiting/terminal step" pattern already used for
  BCAT+BCKDH and choline-kinase+CCT above.

| Reaction ID | Name | Role | Compartment | Gene (Ensembl) | Gene symbol |
|---|---|---|---|---|---|
| `MAR02150` | acetyl-CoA:[acyl-carrier-protein] S-acetyltransferase | first committed step | cytosol (`c`) | ENSG00000169710 | FASN |
| `MAR02182` | hexadecanoyl-[acyl-carrier protein] hydrolase | terminal, palmitate-releasing step | cytosol (`c`) | ENSG00000152463 or ENSG00000169710 | OLAH or FASN |

---

## 5. Choline / phosphatidylcholine (PC) depletion

**Biological rationale**: Choline is required for PC synthesis via the Kennedy pathway; PC can also
be synthesized choline-independently via sequential methylation of phosphatidylethanolamine (PE) by
PEMT, using SAM as the methyl donor. Choline/PC depletion is implicated in NAFLD both directly
(impaired VLDL/triglyceride export, which requires PC) and via its connection to the same one-carbon
methyl-donor pool as the serine/glycine axis (biomarker 1) — a genetic, cohort, and RCT evidence
base per `CLAUDE.md`.

**Reactions selected**:

Kennedy pathway (choline-dependent):

| Reaction ID | Name | Compartment | Gene (Ensembl) | Gene symbol (best-confidence, unverified) |
|---|---|---|---|---|
| `MAR00636` | ATP:choline phosphotransferase (choline kinase) | cytosol (`c`) | ENSG00000100288 or ENSG00000110721 | CHKB or CHKA (confirmed via MyGene.info) |
| `MAR00638` | CTP:choline-phosphate cytidylyltransferase (rate-limiting step, CCT) | cytosol (`c`) | ENSG00000102230 or ENSG00000161217 | PCYT1B or PCYT1A (confirmed via MyGene.info) |

PEMT pathway (choline-independent, compensatory):

| Reaction ID | Name | Compartment | Gene (Ensembl) | Gene symbol (best-confidence, unverified) |
|---|---|---|---|---|
| `MAR00653` | SAM:phosphatidylethanolamine N-methyltransferase | cytosol (`c`) | ENSG00000133027 | PEMT |
| `MAR01603` | Phosphatidylethanolamine N-Methyltransferase | mitochondria (`m`) | ENSG00000133027 | PEMT |
| `MAR01606` | Phosphatidylethanolamine N-Methyltransferase | ER (`r`) | ENSG00000133027 | PEMT |

**Why selected**: `MAR00638` (CCT) is the textbook rate-limiting step of PC synthesis via the Kennedy
pathway — the single best single-reaction proxy if a minimal set is wanted. `MAR00636` (choline
kinase) is the first committed step and included for pathway-level coverage, same reasoning as
biomarkers 1-3 above. PEMT (`MAR00653`/`01603`/`01606`, same gene, 3 compartment variants) is
included specifically because it represents the choline-independent compensatory route — its
inclusion means the biomarker set captures both "not enough choline coming in" and "not enough
methyl-donor-driven compensation," which is closer to the actual mechanistic picture than choline
availability alone.

**Alternatives considered and rejected**: `MAR00629`/`00630`/`00632`/`00633`/`00634`/`00635` (PC
hydrolysis/remodeling/esterification reactions, e.g. phospholipase-type activity, lecithin-cholesterol
acyltransferase-type activity) all involve PC as a substrate or product but represent PC *turnover/
remodeling*, not net *synthesis* — excluded as not directly representing the "depletion" mechanism
this biomarker is about. `MAR00640`/`00641` (acetylcholine synthesis/hydrolysis) are a distinct,
neurotransmitter-related choline use, not the PC-synthesis axis — excluded.

---

## 9. Remaining uncertainties (summary, consolidated from all 5 sections above)

- **Gene symbols verified via MyGene.info (2026-07-21)** — all confirmed, no discrepancies from
  this document's earlier best-confidence guesses:
  - BCKDH complex: `ENSG00000083123`=`BCKDHB`, `ENSG00000137992`=`DBT`, `ENSG00000248098`=`BCKDHA`.
  - `ENSG00000100288`=`CHKB`, `ENSG00000110721`=`CHKA`.
  - `ENSG00000102230`=`PCYT1B`, `ENSG00000161217`=`PCYT1A`.
- **FASN reaction choice — resolved** (biomarker 4): neither original candidate was right.
  `MAR01471` is excluded (no GPR, can't carry an expression signal); `MAR02323`/`02325`/`02327` are
  excluded (traced to a C16→C18 elongation pathway initiated by OXSM, not the primary FASN/DNL
  pathway). Selected instead: `MAR02182` (palmitate-releasing terminal step, GPR includes FASN) and
  `MAR02150` (first committed step, GPR = FASN alone) — see biomarker 4 above for full reasoning.
- **SCD1's three near-duplicate reactions — resolved**: they produce three different C18:1-CoA
  positional isomers, not true duplicates. `MAR00146` (oleoyl-CoA, the literature-referenced SCD1
  product) is selected; `MAR00147`/`MAR00148` (non-canonical positional isomers) are excluded.
- **Urea cycle's CPS1/ASS1/ASL reactions were not found in this session's first pass** (reported
  earlier as "only found OTC and arginase") **and were subsequently found** in a second, more
  targeted search of this same session, documented above — flagging this so the earlier
  under-complete report isn't mistaken for a still-open gap.
- **ARG1 vs. ARG2, SHMT1 vs. SHMT2, BCAT1 vs. BCAT2, and PEMT's 3 compartments**: this document
  currently recommends including all compartment/isoform variants per biomarker for completeness,
  but flags that a single-reaction proxy is possible and may be preferable for calibration
  simplicity — this is a methodology choice for Mohammed to confirm, not decided here.
- **No Step 11 weighting/calibration methodology has been decided yet** — this document only fixes
  *which reactions*, not *how* the eventual `calibrated_confidence_score` combines Layer 1
  (consensus), Step 7 (flux sampling uncertainty), and Step 8 (perturbation robustness) into a score
  that ranks these reactions/pathways appropriately. That's the next design conversation, after this
  mapping is reviewed.
