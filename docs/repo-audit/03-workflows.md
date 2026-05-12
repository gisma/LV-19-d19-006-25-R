# Workflows

## Documentation-Folder Audit Update

- FACT (`src/S1_observation/doc/S1-meta.md`): S1 workflows acquire primary observations through explicit retrieval scripts and allow only technical provider preprocessing such as mosaicking, clipping, reprojection, and format conversion.
- FACT (`src/S2_features/doc/S2-meta.qmd`): S2 workflows are fivefold: S2-1 terrain/hydrology on a canonical 10 m grid, S2-2 DOM-DEM biostructure aggregation, S2-3 DWD wind context as tables/figures, S2-4 kNDVI/bfast temporal signals, and S2-5 CLC-to-segment semantic anchoring.
- FACT (`src/S3_structure/doc/S3-meta.qmd`): S3 workflow standardises predictor bands, applies PCA, segments with OTB LargeScaleMeanShift, evaluates parameter stability with ARI, and produces structural segment polygons.
- FACT (`src/S3_structure/archive/segmentation-meta.md`): S3 scale construction has two documented modes: CHM-driven candidate scales when CHM/DSM is available, and a fallback fixed scale grid; size filtering uses a distribution knee before polygonisation.
- FACT (`src/S4_signatures/doc/S4-meta.qmd`): S4 workflow joins deterministic metric domains onto stable segment geometry without geometry modification.
- FACT (`src/S4_signatures/doc/04-3_signatures_RF-segs_classification.qmd`): The RF/CLC workflow aggregates CLC semantics at object level, constructs segment signatures, applies numerical feature filtering, and requires spatial validation.
- FACT (`src/S5-decision/05-9_decision_selector.qmd`): S5 has two ranking workflows: standalone MCDA for global representativeness and selector settings for constrained design with hard filters, target coverage, ranking, and spacing.
- FACT (`src/S5-decision/05-9_decision_selector.qmd`): S6 is expected to evaluate scenario outputs by coverage completeness, spacing distribution, overlap with the MCDA baseline, and per-domain representativeness summaries.
- QUESTION: The documentation names S6 validation as the next step, but no S6 source documentation folder or implementation was found by the requested `src/` documentation inventory.

