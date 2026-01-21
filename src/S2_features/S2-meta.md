# S2: Derived raster and observation products

Burgwald decision stack: meta-level description

## Conceptual position of S2

S2 denotes the layer of spatially explicit products that are derived directly from observations (S1) and remain defined on the native spatial support of the measurement system. For raster data, this is pixels; for in-situ networks, it is stations.

S2 artefacts are transformations of measurements, not abstractions of spatial objects. They preserve physical or statistical interpretability and remain anchored in the observational domain. Typical S2 products include raster layers, per-pixel indices, pixel-wise classifications, probability fields and descriptive station summaries.

S2 addresses a specific yet fundamental question:

> What can be derived directly from measurements in space and time, before any segmentation, object formation or decision logic is introduced?

This restriction is intentional. It prevents premature semantic commitments and keeps uncertainty and scale effects explicit.

## Boundary conditions: what S2 explicitly excludes

S2 does not contain spatial objects, segmentations or region abstractions.
These belong to S3. It also excludes aggregated signatures or descriptors of objects, which belong to S4. Similarly, it does not contain learned representations or embeddings, which are assigned to S4L. It also excludes ranking, selection, optimisation and policy logic, which belong to S5.

A product leaves the S2 domain whenever it requires reasoning over spatial units rather than pixels or stations, or whenever it encodes semantic decisions rather than measurements.

##  Functional role of S2 in the decision stack

S2 fulfils three complementary functions.

First, it produces physically interpretable predictor layers. Raw measurements are transformed into quantities that express terrain structure, hydrological context, vegetation structure, spectral behaviour or temporal dynamics. These layers form the quantitative basis for segmentation, signature construction and representation learning, but are not semantic objects themselves.

Secondly, S2 supports pixel-level interpretation. Semantic information may appear at this stage, provided it is defined per pixel. Land-cover classifications, disturbance maps or probability surfaces are all legitimate S2 products. Pixel semantics do not imply object semantics; the unit of meaning remains the pixel.

Thirdly, S2 includes descriptive observational summaries. Aggregated statistics, diagnostic figures or station climatologies remain part of S2 as long as they are purely descriptive and do not introduce decision-making logic or optimisation criteria.

##  Design principles for S2 processing

All S2 products are subject to strict constraints regarding reproducibility and transparency.

Inputs and outputs are centrally registered via the project registry and must not be hard-coded. Transformations must be explicit; implicit thresholds, undocumented resampling, silent smoothing or hidden aggregation must be avoided.

Spatial support is controlled deliberately. The canonical working grid is a 10 m raster. Higher-resolution sources may be used internally, but they must be aggregated deterministically. Sub-pixel modelling is intentionally excluded at this level to avoid scale ambiguity.

S2 scripts may depend operationally on artefacts from other levels (for example, segmentation used as annotation infrastructure), but they must not modify or overwrite artefacts outside of S2.

## Product families in the Burgwald context

Terrain and hydrological derivatives transform digital elevation models into morphometric and flow-related predictors. These layers express physical structure and serve as quantitative inputs for later stages.

Biostructural derivatives transform surface and terrain models into canopy height statistics and structural indicators. The logic is strictly geometric: height differences are computed, negative artefacts are removed and blockwise statistics are derived on the canonical grid. No object assumptions are introduced.

Pixel-based classifications generate categorical or probabilistic land cover information from predictor stacks. Training data can be stabilised using segmentation as geometric scaffolding, but the prediction remains pixel-wise, so the result is still an S2 artefact.

Time-series analyses, such as change detection, operate on pixel stacks and produce per-pixel diagnostics of temporal behaviour.

Observational climate summaries condense station or gridded measurements into descriptive statistics and visualisations, without applying any decision logic.

##  Semantic versus organisational ordering

The semantic level of a product is defined by its meaning rather than the execution order in the processing pipeline. Some S2 products may depend on S3 artefacts for practical reasons; for example, when segmentation stabilises training region delineation. While this introduces an operational dependency, it does not change the semantic classification of the product as S2, provided the output remains pixel-based and non-object-centric.

This distinction between the semantic level and the orchestration order is essential in order to maintain a flexible architecture without blurring conceptual responsibilities.

## Contract statement

S2 comprises all spatially explicit products that are derived directly from observations, prior to any segmentation, abstraction or decision logic. These products preserve physical or statistical interpretability and form the quantitative basis for all higher-level representations and decisions.
