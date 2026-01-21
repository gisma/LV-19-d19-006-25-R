# S1: Observations and primary data products  
Burgwald decision stack: meta-level description

## Conceptual position of S1

S1 denotes the layer of primary observations and externally provided datasets on
which the entire decision stack is built. S1 artefacts represent measurements,
surveys, or authoritative reference products that are not derived within the
project itself. They constitute the empirical substrate of the system.

S1 products preserve the spatial and temporal support defined by the original
measurement process. For remote sensing data this is the sensor pixel and its
acquisition time; for in-situ networks it is the station and its sampling
interval; for authoritative geodata it is the provider-defined spatial grid or
vector geometry.

S1 addresses a single foundational question:

> What has been observed, measured, or authoritatively provided, prior to any
> project-specific transformation or interpretation?

S1 deliberately avoids any semantic enrichment beyond what is inherent in the
original data source.

---

## Boundary conditions: what S1 explicitly excludes

S1 does not contain derived features, indices, interpolations, classifications,
or aggregations. These belong to S2. It does not contain segmentation objects,
which belong to S3, nor any higher-level abstractions, representations, or
decisions.

Whenever a product reflects a transformation, combination, or interpretation of
measurements performed within the project, it has left the S1 domain.

---

## Functional role of S1 in the decision stack

S1 fulfils three tightly coupled functions.

First, it provides the authoritative geometric and thematic reference for the
project. Digital terrain and surface models define the geometric substrate.
Administrative and land-cover reference layers define spatial context.
OpenStreetMap layers provide infrastructure and land-use context. These datasets
anchor all subsequent processing steps in a consistent spatial frame.

Second, S1 provides raw observational input for quantitative modelling. Satellite
imagery delivers multi-spectral measurements with known radiometric and geometric
properties. In-situ networks deliver meteorological time series with documented
sampling characteristics and uncertainties. These observations remain unaltered
at S1 and serve as the basis for all derived quantities.

Third, S1 provides traceability and auditability. Because S1 products are either
external authoritative datasets or raw measurements, they enable downstream
products to be traced back to their empirical origin. This is essential for
reproducibility, error attribution, and methodological transparency.

---

## Acquisition and handling principles

S1 data are acquired via explicit, version-controlled retrieval scripts. The
scripts document source endpoints, spatial extent, temporal coverage, and access
constraints. All downloads are deterministic with respect to the chosen
parameters and are stored under canonical project paths.

Provider-specific preprocessing is limited to operations that are required to
make the data technically usable, such as mosaicking tiles, clipping to the area
of interest, reprojection to a common coordinate reference system, or format
conversion. No semantic transformation or analytical interpretation is performed
at this stage.

Temporal coverage and update logic are controlled explicitly by the retrieval
scripts. Changes in upstream provider data therefore remain visible and
auditable.

---

## Product families in the Burgwald context

Base geodata comprise digital terrain models, digital surface models, reference
land-cover layers, and selected OpenStreetMap extracts. These datasets define the
geometric and thematic baseline of the study area.

Remote sensing observations comprise Sentinel-2 Level-2A scenes retrieved via
STAC-based interfaces or cloud-native data cubes. The data preserve original
radiometric measurements and acquisition geometry.

Meteorological observations comprise hourly and sub-hourly station data for wind
and precipitation retrieved from the national data provider. These data preserve
the original temporal resolution and station geometry.

All these products remain observational and are not analytically modified within
S1.

---

## Semantic versus organisational ordering

S1 defines the semantic origin of data, not necessarily their position in the
execution pipeline. Some S1 retrieval steps may be triggered conditionally or
lazily, depending on downstream requirements. This does not change their role as
primary observations.

The defining criterion is whether a product represents an external measurement
or authoritative dataset rather than a project-internal transformation.

---

## Contract statement

S1 comprises all primary observations and externally provided datasets that enter
the system without project-specific analytical transformation. These artefacts
preserve the spatial and temporal support of the original measurement process and
form the empirical foundation for all subsequent processing, representation, and
decision layers.
