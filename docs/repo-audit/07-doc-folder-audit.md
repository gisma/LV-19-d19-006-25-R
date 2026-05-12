# Documentation-Folder Audit

## Required Inventory

### Documentation Folders Found

```text
src/S1_observation/doc
src/S2_features/doc
src/S3_structure/doc
src/S4_signatures/doc
src/S5-decision/doc
```

### Documentation Files Found

```text
src/README_1st.md
src/S1_observation/doc/S1-meta.md
src/S2_features/doc/S2-meta.qmd
src/S3_structure/archive/meta.md
src/S3_structure/archive/segmentation-meta.md
src/S3_structure/doc/S3-meta.qmd
src/S4_signatures/doc/04-3_signatures_RF-segs_classification.qmd
src/S4_signatures/doc/S4-meta.qmd
src/S5-decision/05-9_decision_selector.qmd
src/S5-decision/doc/S5-meta.qmd
src/_core/meta.md
src/r-libs/meta.md
src/run/meta.md
src/tools/meta.md
```

## Documentation Files Explicitly Read

- FACT: All files listed in the inventory above were read for this targeted documentation-folder audit.
- FACT: These files were treated as primary sources for workflow, conceptual structure, layer semantics, and intended use.

## Per-Document Audit

### `src/README_1st.md`

- FACT: Layer/module: top-level `src` meta-architecture.
- FACT: Main purpose: Defines a spatial research pipeline for optimized precipitation network design in Burgwald.
- FACT: Workflow statements: Describes the transformation chain from S1 observations through S2 features, S3 structure, S4 signatures, optional S4L learned signatures, S5 decisions, and S6 validation.
- FACT: Conceptual statements: The project is about how spatial knowledge is constructed; every step abstracts reality and introduces assumptions, simplifications, and information loss.
- FACT: Links to scripts, outputs, or registry keys: Names `layer0_segments_attrstack_metrics (GPKG)` as the productive S4 artefact and describes S5 domain variables and MCDA inputs.
- FACT: Contradictions or gaps compared with the previous audit: The documented repository mapping lists `src/S5_decisions` and `src/lib`, while the actual repository uses `src/S5-decision` and `src/r-libs`.

### `src/S1_observation/doc/S1-meta.md`

- FACT: Layer/module: S1 Observation.
- FACT: Main purpose: Defines S1 as primary observations and externally provided datasets.
- FACT: Workflow statements: Acquisition is via explicit retrieval scripts; provider preprocessing is limited to mosaicking, clipping, reprojection, and format conversion.
- FACT: Conceptual statements: S1 preserves original spatial and temporal support and avoids semantic enrichment beyond the source data.
- FACT: Links to scripts, outputs, or registry keys: Mentions base geodata, Sentinel-2 observations, and DWD meteorological observations as product families.
- FACT: Contradictions or gaps compared with the previous audit: The file sharpens the boundary that derived indices, interpolations, classifications, aggregations, segmentation objects, and decisions are not S1.

### `src/S2_features/doc/S2-meta.qmd`

- FACT: Layer/module: S2 Features.
- FACT: Main purpose: Explains S2 as the first layer where data become representation.
- FACT: Workflow statements: Defines five pipelines: terrain/hydrology, vegetation biostructure, atmospheric context, temporal disturbance signals, and semantic anchoring.
- FACT: Conceptual statements: S2 representations are not predictions, classifications, or semantic labels; they encode assumptions about spatial structure, scale, noise, and relevance.
- FACT: Links to scripts, outputs, or registry keys: Describes DEM-derived 10 m grids, dCHM aggregation, wind tables/figures, kNDVI/bfast descriptors, and segment purity/training data.
- FACT: Contradictions or gaps compared with the previous audit: Adds that S2-3 wind context intentionally remains orthogonal to the spatial modelling chain and is not rasterized into the predictor stack.

### `src/S3_structure/doc/S3-meta.qmd`

- FACT: Layer/module: S3 Structure.
- FACT: Main purpose: Defines segmentation and structural spatial units.
- FACT: Workflow statements: Predictor bands are z-scored, PCA is applied, OTB LargeScaleMeanShift segments joint spatial-feature space, ARI perturbation screening controls stability, and postconditions enforce structural constraints.
- FACT: Conceptual statements: S3 segments are operational structural entities, not semantic objects.
- FACT: Links to scripts, outputs, or registry keys: Names MeanShift parameters `spatialr`, `ranger`, and `minsize`; aligns with `src/S3_structure/03-1_analysis_base-segmentation.R`.
- FACT: Contradictions or gaps compared with the previous audit: Clarifies that S-layer assignment is not identical to execution order.

### `src/S3_structure/archive/segmentation-meta.md`

- FACT: Layer/module: S3 archive documentation.
- FACT: Main purpose: Detailed operational documentation for multi-scale MeanShift segmentation with perturbation-based stability screening.
- FACT: Workflow statements: Documents z-score, PCA, MeanShift, ARI_prev, scale selection, CHM-driven scale construction, fallback grid, and knee-based size filtering.
- FACT: Conceptual statements: Stability is not validation and does not imply ecological correctness; it enforces operational consistency.
- FACT: Links to scripts, outputs, or registry keys: Names OTB BandMathX, DimensionalityReduction, LargeScaleMeanShift and parameters such as `sample_fact`, `spatialr`, `ranger`, `minsize`.
- FACT: Contradictions or gaps compared with the previous audit: Adds important fallback-grid and CHM-driven scale semantics that the first audit underrepresented.

### `src/S3_structure/archive/meta.md`

- FACT: Layer/module: S3 archive.
- FACT: Main purpose: Defines archive scripts as experimental, legacy, or redundant.
- FACT: Workflow statements: Archive scripts are not part of the standard pipeline.
- FACT: Conceptual statements: Archive material is retained for methodological comparison and historical reference.
- FACT: Links to scripts, outputs, or registry keys: No specific keys.
- FACT: Contradictions or gaps compared with the previous audit: Adds an explicit rule that archive scripts should have no production dependency unless stated.

### `src/S4_signatures/doc/S4-meta.qmd`

- FACT: Layer/module: S4 Signatures.
- FACT: Main purpose: Defines signatures as explicit structural representations attached to stable segments.
- FACT: Workflow statements: S4 consumes segment geometry and derived raster fields, then joins metric tables into an attribute stack without changing geometry.
- FACT: Conceptual statements: A signature is a deterministic representational contract, not merely a performance-oriented machine-learning feature vector.
- FACT: Links to scripts, outputs, or registry keys: Names the attribute stack as the canonical S4 artefact and describes information-theoretic, physiographic, hydrological, and biostructural domains.
- FACT: Contradictions or gaps compared with the previous audit: Adds symbolic tokens/controlled compression as an intended S4 concept; implementation and registry location remain unclear.

### `src/S4_signatures/doc/04-3_signatures_RF-segs_classification.qmd`

- FACT: Layer/module: S4/S4L transition documentation.
- FACT: Main purpose: Documents segment-based CLC supervision and RF classification on Sentinel signatures.
- FACT: Workflow statements: CLC is aggregated at object level; segment statistics are used rather than pixels; caret preprocessing removes near-zero variance, linear dependencies, and high correlations; spatial validation is mandatory.
- FACT: Conceptual statements: The objective is scale stabilization under semantic control, not land-cover mapping as a cartographic product.
- FACT: Links to scripts, outputs, or registry keys: Links conceptually to the segment layer, Sentinel predictor stack, CLC raster, supervised segment layer, predicted segment layer, and persistent training table.
- FACT: Contradictions or gaps compared with the previous audit: The document says learning operates instrumentally in S4L, while the implementation is located under `src/S4_signatures`.

### `src/S5-decision/doc/S5-meta.qmd`

- FACT: Layer/module: S5 Decision.
- FACT: Main purpose: Defines S5 as decision and candidate selection space.
- FACT: Workflow statements: S5 separates metric harmonisation, evaluation, and selection.
- FACT: Conceptual statements: S5 is not modelling or discovery; weights, thresholds, and selector rules encode value judgements.
- FACT: Links to scripts, outputs, or registry keys: Mentions candidate-generating logic, strata, distance-to-centre measures, MCDA, selector rules, quotas, and spacing.
- FACT: Contradictions or gaps compared with the previous audit: The title says "Validation and Selection Space", but validation is otherwise defined as S6; it also says candidate-generating logic is completed before S5, while other docs place candidate construction in S5.

### `src/S5-decision/05-9_decision_selector.qmd`

- FACT: Layer/module: S5 selector tutorial.
- FACT: Main purpose: Explains S4 signatures -> S5 decisions -> S6 validation.
- FACT: Workflow statements: S4 provides `layer0_segments_attrstack_metrics.gpkg`; S5 derives per-domain and selector-based candidates; S6 compares scenario outputs.
- FACT: Conceptual statements: MCDA is a scoring module; selector is the design engine combining constraints, coverage, and spacing.
- FACT: Links to scripts, outputs, or registry keys: Names `05-9_decisions_mcda_fuse_metrics.R`, `05-9_decisions_selector_settings.R`, `mcda_add_score()` in `metrics-fun.R`, `layer0_segments_attrstack_metrics.gpkg`, and expected distance/stratum/filter columns.
- FACT: Contradictions or gaps compared with the previous audit: Adds explicit scenario settings A/B/C and S6 validation expectations that were underrepresented before.

### `src/_core/meta.md`

- FACT: Layer/module: Core Infrastructure.
- FACT: Main purpose: Stable project infrastructure for setup, paths, toolchain registration, global options, and constants.
- FACT: Workflow statements: Core should contain no data processing or modelling logic.
- FACT: Conceptual statements: Core is not part of the scientific model.
- FACT: Links to scripts, outputs, or registry keys: Implies setup and path registry responsibilities.
- FACT: Contradictions or gaps compared with the previous audit: Raises a gap where code may use hard-coded local tool paths despite the rule against hard-coded absolute paths.

### `src/r-libs/meta.md`

- FACT: Layer/module: Library layer.
- FACT: Main purpose: Reusable functional building blocks.
- FACT: Workflow statements: Functions should be stateless, explicit, path-independent, side-effect-free, and unit-agnostic where possible.
- FACT: Conceptual statements: Library layer is infrastructure only.
- FACT: Links to scripts, outputs, or registry keys: No direct keys.
- FACT: Contradictions or gaps compared with the previous audit: Creates an audit question for helper functions that download, write, or depend on project-specific paths.

### `src/run/meta.md`

- FACT: Layer/module: Run layer.
- FACT: Main purpose: Explicit entry points and standard workflow orchestration.
- FACT: Workflow statements: Typical usage examples are `run_S1_all.R`, `run_S2_all.R`, and `run_S3_all.R`.
- FACT: Conceptual statements: Execution logic only; no modelling logic.
- FACT: Links to scripts, outputs, or registry keys: Names expected run scripts.
- FACT: Contradictions or gaps compared with the previous audit: Inventory found only `src/run/meta.md`, not the named run scripts.

### `src/tools/meta.md`

- FACT: Layer/module: Tool layer.
- FACT: Main purpose: External tooling and glue scripts.
- FACT: Workflow statements: Includes JavaScript helpers, shell scripts, API helpers, and external engine wrappers.
- FACT: Conceptual statements: Pure infrastructure; changes here must not alter scientific assumptions.
- FACT: Links to scripts, outputs, or registry keys: Aligns with `src/tools/*.js` vegetation-index helpers.
- FACT: Contradictions or gaps compared with the previous audit: No major contradiction; it clarifies tool-layer epistemic status.

