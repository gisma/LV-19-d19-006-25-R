# Architecture

## Documentation-Folder Audit Update

- FACT (`src/README_1st.md`): The architecture is an epistemic S-layer model, not only a directory convention. S0 defines the problem, S1 observations, S2 representations/features, S3 structural units, S4 deterministic signatures, optional S4L learned representations, S5 decisions, and S6 validation.
- FACT (`src/S1_observation/doc/S1-meta.md`): S1 contains primary observations and externally provided datasets only; it explicitly excludes derived features, indices, interpolations, classifications, aggregations, segmentation objects, and decisions.
- FACT (`src/S2_features/doc/S2-meta.qmd`): S2 is the first layer where data become representation. It contains five pipelines: terrain/hydrology, vegetation biostructure, atmospheric context, temporal disturbance signals, and semantic anchoring.
- FACT (`src/S3_structure/doc/S3-meta.qmd`; `src/S3_structure/archive/segmentation-meta.md`): S3 converts continuous raster fields into operational structural units through z-score standardisation, PCA, OTB MeanShift segmentation, ARI-based stability screening, and size-filtered polygonisation.
- FACT (`src/S4_signatures/doc/S4-meta.qmd`): S4 turns stable segment geometry into explicit representable objects by attaching deterministic signatures; the attribute stack is the canonical S4 interface.
- FACT (`src/S5-decision/doc/S5-meta.qmd`): S5 separates metric harmonisation, evaluation, and selection; weights, thresholds, and selector rules encode explicit value judgements.
- FACT (`src/S5-decision/05-9_decision_selector.qmd`): The practical S5 architecture uses `layer0_segments_attrstack_metrics.gpkg` as the truth source; MCDA is the scoring module and the selector is the design engine combining constraints, coverage, and spacing.
- QUESTION (`src/S5-decision/doc/S5-meta.qmd` vs. `src/README_1st.md` and `src/S5-decision/05-9_decision_selector.qmd`): Is domain stratification conceptually completed before S5, or is it part of S5 candidate construction? The documentation uses both framings.
- QUESTION (`src/run/meta.md`): The run layer is documented as the place for `run_S1_all.R`, `run_S2_all.R`, and `run_S3_all.R`, but those files were not found in the documentation inventory.

