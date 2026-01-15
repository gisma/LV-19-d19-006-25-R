# S1 — Observation Space

## Purpose
Construct a reproducible spatial observation baseline.

This layer creates harmonised representations of the physical world.

## Typical Inputs
- Raw geodata (DEM, land cover, Sentinel, DWD)
- External APIs and archives
- AOI definitions

## Typical Outputs
- Harmonised rasters and vectors
- Reference geometry
- Cached raw products

## Abstraction
Physical reality → Discrete measurable layers

## Gains
- Reproducibility
- Spatial consistency
- Controlled uncertainty

## Losses
- Sensor uncertainty
- Temporal generalisation
- Resolution limits

## Special Note
Multiple Sentinel pipelines exist intentionally to expose abstraction trade-offs
(server-side processing, cloud cubes, fully manual workflows).

## Not in Scope
- Interpretation
- Feature engineering
- Segmentation
- Decision logic
