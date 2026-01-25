# Spatial Research Pipeline — Meta‑Architecture and Didactic Core Model

Optimised Precipitation Network Design — Burgwald

This project is not a collection of scripts.
It is a structured research pipeline that transforms spatial observations into a justified decision space.

This project is not about running scripts.
It is about understanding how spatial knowledge is constructed.

Every step transforms reality into a more abstract representation.
Every transformation introduces assumptions, simplifications and loss of information.

The primary learning goal is conceptual modelling and epistemic reasoning, not programming.

---

## Conceptual Layers and Epistemic Roles

Each layer transforms information and introduces abstraction and uncertainty.

| Layer | Name                               | Epistemic Role                                                    |
| ----: | ---------------------------------- | ----------------------------------------------------------------- |
|    S0 | Problem Space                      | Research question, scale logic, decision targets                  |
|    S1 | Observation Space                  | Reproducible construction of spatial observations                 |
|    S2 | Feature Space                      | Process proxies and abstract descriptors                          |
|    S3 | Structural Space                   | Segmentation, stratification, typification                        |
|    S4 | Signature Space                    | Deterministic, interpretable segment descriptors                  |
|   S4L | Learned Signature Space (optional) | Learned representations, embeddings, regime assignments           |
|    S5 | Decision Space                     | Candidate construction, constraint handling, trade‑off resolution |
|    S6 | Validation Space                   | Stability, sensitivity, robustness, transferability               |

> Note: S4L is an explicit optional extension. The core pipeline remains S0–S5.

---

## How to Read This Document (Didactic Reader Guidance)

This document is not a workflow manual and not a software specification. It is a conceptual map of how spatial knowledge is constructed, constrained and justified in this project. The modelling spaces are not sequential processing steps but epistemic commitments: each space defines what counts as an object, a variable, a structure, a candidate and a valid argument.

When reading, focus on three questions:

1. What is being abstracted away at this stage?
2. Which assumptions become fixed and no longer negotiable downstream?
3. Which kinds of errors become structural rather than accidental?

The purpose is not to memorise the pipeline but to understand where meaning, uncertainty and responsibility enter the system. The concrete scripts are merely implementations of these commitments.

---

## The Modelling Spaces

The pipeline is organised as a sequence of modelling spaces. Each space represents a deliberate epistemic transformation: data are reduced, structured, abstracted and normalised in order to become operable for reasoning and decision making. None of these spaces is merely technical. Each introduces assumptions about scale, relevance, causality and admissible uncertainty.

The spaces are not simply sequential processing stages. They represent different epistemic commitments about what counts as an object, a signal, a pattern or a candidate. Errors or mismatches introduced early propagate and cannot be repaired by downstream sophistication.

### S0 — Problem Space

The Problem Space defines what the system is actually about. It specifies the scientific or operational question, the spatial and temporal scales at which the question is meaningful, and the criteria by which later results will be judged. This includes explicit decisions about what constitutes success, what trade‑offs are acceptable, and which uncertainties can be tolerated.

S0 is intentionally not computational. It is a modelling commitment. If S0 is ill-posed — for example by mixing incompatible scales, vague objectives, or hidden normative assumptions — no amount of downstream optimisation can compensate for this.

### S1 — Observation Space

The Observation Space defines how reality is operationalised as data. It specifies which sensors, archives and products are admitted as representations of the physical world, and how these heterogeneous sources are harmonised into a coherent spatial reference frame.

In the Burgwald pipeline this includes digital elevation models, Sentinel‑2 spectral stacks, land‑cover products and DWD station archives. Preprocessing steps such as mosaicking, reprojection, resampling, temporal selection and masking are not technical details; they are epistemic filters that determine which aspects of reality remain visible and which are suppressed.

Uncertainty enters explicitly at this stage through sensor limitations, temporal mismatches, spatial aggregation and preprocessing artefacts. S1 produces reproducible but still heterogeneous observations, not yet interpretable process variables.

### S2 — Feature Space

The Feature Space transforms raw observations into abstract descriptors that are assumed to be informative about relevant processes. Terrain derivatives, spectral indices, texture measures, hydrological predictors or vegetation proxies are not measurements of processes themselves; they are hypotheses encoded as variables.

Aggregation windows, scaling choices, normalization strategies and proxy definitions embed strong modelling assumptions. S2 therefore marks the first major abstraction step: continuous physical signals are compressed into low‑dimensional descriptors that trade fidelity for tractability and comparability.

Errors introduced here tend to be structural rather than random, because the feature definitions constrain what kinds of patterns can ever become visible downstream.

### S3 — Structural Space

The Structural Space discretises continuous fields into explicit spatial units. Segmentation, stratification and typification convert raster fields into polygons, objects or strata that become the atomic carriers of later signatures and decisions.

In this project, segmentation is performed using MeanShift with stability screening via local parameter perturbation and ARI-based agreement analysis. Minimum size constraints and polygon extraction enforce geometric regularity and computational tractability.

Crucially, structure is not discovered in a neutral sense. It is imposed by algorithmic choices, scale parameters and stability thresholds. S3 therefore defines the geometry of the decision space itself: what counts as a unit, what counts as neighbourhood, and what kinds of heterogeneity can even be expressed.

The output of S3 is still primarily geometric. Segments exist as spatial containers, not yet as comparable entities in an abstract attribute space.

### Transition: From Structure to Signature (S3 → S4)

The transition from S3 to S4 is the decisive epistemic pivot of the pipeline. Geometry is converted into vectors of attributes. Spatial objects become elements of a metric space in which similarity, distance, clustering and ranking become formally meaningful operations.

At this point, spatial relations are no longer represented primarily through topology or adjacency, but through numerical descriptors derived from aggregated measurements. This conversion enables comparison, fusion and multi-domain reasoning — but it also collapses spatial context, internal heterogeneity and geometric nuance into summary statistics.

Errors introduced here are particularly consequential: once encoded as signatures, downstream processes cannot recover lost spatial structure. The choice of metrics therefore determines which kinds of explanations and decisions remain expressible.

### S4 — Signature Space

The Signature Space assembles interpretable, deterministic descriptors for each structural unit. All downstream reasoning operates exclusively on this signature table.

The productive artefact is the consolidated attribute stack:

```
layer0_segments_attrstack_metrics (GPKG)
```

This table fuses multiple metric domains: information-theoretic measures (entropy, mutual information), physiographic descriptors (elevation, slope, southness), hydrological indicators (flow accumulation, stream distance, network order), biostructural metrics (canopy height, canopy fraction, structural variability) and land-cover composition fractions.

S4 deliberately contains no decision logic, optimisation or ranking. It defines the representational state space in which all later reasoning takes place. Signatures are intended to be interpretable, reproducible and stable under controlled perturbations of upstream processing.

The quality and internal consistency of S4 determine the expressive power and epistemic limits of the entire pipeline.

---

## Optional Extension — S4L: Learned Signatures

S4L introduces representation learning in a controlled and auditable way. Instead of using machine learning as an implicit decision engine, learning is confined to the construction of additional representations: embeddings, regime memberships or latent typologies.

These learned artifacts extend the signature space but remain frozen inputs to decision logic. Training provenance, versioning, diagnostics and stability assessments are treated as first‑class metadata. This separation prevents the entanglement of representation learning with decision rules and preserves interpretability and auditability at the decision layer.

S4L is optional and does not modify the core S0–S5 logic.

S4L introduces representation learning without collapsing decision logic into ML.

Typical artifacts:

* embeddings
* regime memberships
* latent typologies

Constraints:

* versioned outputs
* explicit training provenance
* frozen representations
* diagnostics attached to each artifact

S5 may consume S4 and S4L jointly — but never learn representations internally.

---

## S5 — Decision and Candidate Selection Space

The Decision Space operates on signatures and produces explicit candidate sets. Its role is not to deliver a single optimal solution, but to make the construction of admissible alternatives transparent, controllable and reviewable.

S5 formalises judgement rather than replacing it. Constraints, coverage logic, ranking rules and spacing criteria are expressed explicitly so that their effects can be inspected, challenged and modified. The resulting candidates are therefore not answers, but structured proposals whose rationale can be traced back through the pipeline.

Normatively, S5 avoids two common traps: treating scoring as objective truth, and collapsing decision making into optimisation. Scores are only relational ordering devices; they do not carry intrinsic meaning. Decisions remain situated, contingent and accountable.

The Decision Space transforms signatures into explicit candidate sets. Unlike typical optimisation pipelines, S5 does not attempt to compute a single optimal solution. Instead, it constructs admissible and interpretable candidate spaces under explicitly stated constraints.

S5 operates on the assumption that spatial decision problems are fundamentally underdetermined: multiple structurally different solutions may be equally defensible depending on weighting, coverage logic and spatial spacing constraints. The role of S5 is therefore to make these trade‑offs explicit, reproducible and inspectable.

Conceptually, S5 integrates three coupled mechanisms.

First, the admissible design space is restricted using hard constraints. Filter expressions define which segments are even allowed to participate in candidate construction (for example forest-only constraints, watershed masks or minimum data completeness thresholds).

Second, the admissible segments are structured into domain-specific representational spaces. Segments are clustered within selected metric subspaces, and distance to the respective cluster center operationalises representativeness. This produces explicit strata and quantitative proximity measures rather than vague similarity notions.

Third, candidates are selected using deterministic ranking and spacing logic. Ranking resolves multi-criteria trade-offs either lexicographically or via explicit MCDA scoring. Spatial thinning enforces geometric independence and prevents spatial clustering artefacts.

The output of S5 is never a surface or a prediction map. It is a finite, explicitly constructed set of candidate locations and associated justification metadata.

S5 transforms signatures into explicit candidate sets.

S5 is not optimisation in the abstract sense. It is *controlled construction of admissible candidate spaces*.

S5 fulfils three coupled roles:

1. **Design space restriction**
   Hard filters define admissible segments (e.g. forest‑only, watershed masks).

2. **Candidate space structuring**
   Segments are clustered into strata in multiple domain‑specific spaces. Distance to cluster center formalises representativeness.

3. **Selection and trade‑off resolution**
   Candidates are ranked using lexicographic ordering or MCDA. Optional spatial spacing enforces geometric independence.

Outputs are explicit candidate sets (segments + representative points), never continuous fields.

---

## Domain‑Specific Decision Spaces in S5

Each decision space operates on a subset of S4 signatures and produces:

* stratum assignment per segment
* distance‑to‑center metric
* representative candidate segments
* representative candidate points

### Information‑Theoretic Space (IT)

Input variables:

* H_norm (or H)
* U

Method:

* scaling
* silhouette‑based k selection
* k‑means clustering
* Euclidean distance to cluster center

Purpose:

Capture landscape structural heterogeneity and compositional complexity.

---

### Physiographic Space

Input variables:

* elev_mean
* slope_mean_deg
* southness_mean

Purpose:

Capture macro‑topographic process regimes (radiation exposure, drainage context, orographic gradients).

---

### Hydrological Space

Input variables:

* strahler_max
* flowacc_p90
* dist_stream_mean

Purpose:

Capture hydrological connectivity and flow regime differentiation.

---

### Biostructural Space

Input variables:

* chm_p95_mean
* canopy_fraction_mean
* chm_sd_mean

Purpose:

Capture vegetation structure and canopy complexity relevant for interception and microclimate.

---

### Coverage Space

Input variables:

* cov_forest
* cov_agri
* cov_grass
* cov_built
* cov_water

Purpose:

Capture land‑cover composition and surface heterogeneity.

---

### MCDA Fusion Space

Inputs:

* per‑domain distance‑to‑center metrics

Processing:

* min‑max normalization
* weighted linear aggregation
* global ranking

Purpose:

Construct cross‑domain representativeness rankings.

MCDA is strictly a scoring model — not a learning model.

---

## Selector Layer — Constraint‑Based Candidate Construction

A separate selector module applies scenario‑specific constraints:

Each scenario defines:

1. Hard filter expression
2. Target coverage variable
3. Quota per target stratum
4. Ranking mode

   * lexicographic
   * MCDA
5. Minimum spatial spacing

Outputs are scenario‑specific candidate sets stored in temporary workspaces for rapid iteration.

This layer supports systematic exploration of design hypotheses without contaminating productive outputs.

---

## Design Principles for S5

* No geometry modification in S5
* No implicit optimisation
* All filters explicit
* All weights explicit
* Deterministic reproducibility
* Clear separation between scoring and selection

---

## Repository Mapping

| Conceptual Layer | Folder             |
| ---------------- | ------------------ |
| Infrastructure   | src/_core          |
| S1 Observation   | src/S1_observation |
| S2 Features      | src/S2_features    |
| S3 Structure     | src/S3_structure   |
| S4 Signatures    | src/S4_signatures  |
| S5 Decisions     | src/S5_decisions   |
| Libraries        | src/lib            |
| Tools            | src/tools          |
| Orchestration    | src/run            |

---

## Typical Failure Modes

* Confusing signatures with decisions
* Treating clusters as ontological truth
* Hiding filters inside scripts
* Mixing learning and decision logic
* Over‑tuning MCDA weights
* Ignoring spatial spacing effects

---

## Analyst Responsibility

The analyst is responsible for maintaining epistemic clarity across the pipeline. This includes knowing which modelling space is currently active, which assumptions have already been fixed upstream, and which degrees of freedom remain legitimately open.

Responsibility also means making abstraction choices explicit, documenting why particular representations were chosen, and being able to explain how alternative abstractions would plausibly change the outcome. Code alone is never sufficient documentation of reasoning.

---

## Success Criterion

A result is successful when its origin can be reconstructed: why it looks as it does, which modelling commitments produced it, and under which alternative assumptions it would differ. Reproducibility without interpretability is insufficient.

---

## What This Project Does NOT Provide

This project does not promise a universal optimum, automated judgement or guaranteed transferability. Its purpose is to make spatial reasoning explicit, inspectable and defensible — not to eliminate uncertainty or replace expert responsibility.
