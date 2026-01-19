# Multi-Scale MeanShift Segmentation with Perturbation-Based Stability Screening (ARI_prev)

## Purpose

This workflow constructs spatial segmentation units from a Sentinel-derived predictor stack using OTB MeanShift segmentation.  
Instead of selecting segmentation parameters heuristically or visually, parameter choice is treated as an **operational control problem** and screened via **local stability under small perturbations**.

The goal is **not** to infer a “true physical scale”, but to identify parameter settings that:

- produce internally consistent segmentations,
- are locally robust against small parameter changes,
- remain computationally tractable and reproducible,
- form a stable structural basis for subsequent modelling stages.

---

## Conceptual Overview

### Pipeline Stages

1. **Z-score standardisation (BandMathX)**  
   Each predictor band is standardised to zero mean and unit variance:
   \[
   z = \frac{x - \mu}{\sigma}
   \]
   This removes scale dominance and ensures comparability across spectral, texture, and auxiliary features.

2. **PCA feature extraction (DimensionalityReduction)**  
   PCA reduces dimensionality and collinearity while preserving dominant variance structure.  
   This is a pragmatic feature transform, not a process model.

3. **MeanShift segmentation (LargeScaleMeanShift)**  
   Segmentation operates in joint spatial–feature space using three parameters:

   - **spatialr** – spatial bandwidth (pixel radius)
   - **ranger**   – feature-space bandwidth (similarity tolerance)
   - **minsize**  – minimum region size in region merging

4. **Stability screening via ARI_prev**  
   For each base parameter set:
   - Run baseline segmentation.
   - Generate a small set of local perturbations.
   - Re-run segmentation for each perturbation.
   - Compute Adjusted Rand Index (ARI) between baseline and perturbed labelings.
   - Aggregate:
     \[
     \text{ARI\_prev} = \text{mean}(\text{ARI}_i), \quad
     \sigma_{\text{ARI}} = \text{sd}(\text{ARI}_i)
     \]

5. **Scale selection**
   A composite robustness score is used:
   \[
   \text{score} = \text{ARI\_prev} - 0.5 \cdot \sigma_{\text{ARI}}
   \]
   The best scale maximises stability while penalising variability.

6. **Postcondition: size filtering**
   After final segmentation:
   - Label size distribution is analysed.
   - A jump / knee point in the distribution determines a minimum segment size.
   - Small fragments are masked before polygonisation.

---

## Why Stability Instead of “True Scale”?

Segmentation scale is inherently circular:

- The segmentation determines what structures become visible.
- Those structures would later be used to justify the segmentation scale.

To avoid implicit confirmation bias, scale selection is **not interpreted physically**, but treated as a **numerical stability property**:

> If small parameter perturbations produce similar partitions, the segmentation is locally robust.

This does **not** guarantee ecological correctness — it only enforces operational consistency.

---

## Adjusted Rand Index (ARI)

ARI measures similarity between two partitions corrected for chance agreement:

\[
\text{ARI} = \frac{\text{Index} - \text{Expected}}{\text{Max} - \text{Expected}}
\]

Properties:

- Range: \([-1, 1]\)
- 1 = identical partitions
- 0 = random agreement
- Negative = worse than random

In this workflow:

- ARI is computed on downsampled label rasters (nearest neighbour) for efficiency.
- Only paired valid cells are used.
- No full contingency tables are constructed (memory-safe implementation).

ARI here is **not** a validation metric — it is a *relative stability indicator*.

---

## Local Perturbation Strategy

Perturbations are generated adaptively around a base parameter set.

### Adaptive rules

Let the base parameters be:

- \( s = \text{spatialr} \)
- \( r = \text{ranger} \)
- \( m = \text{minsize} \)

Then:

- **Spatial perturbation**
  - If \( s \le 3 \): no spatial perturbation (relative change would be too large).
  - Else: \( s \pm 1 \)

- **Range perturbation**
  - \( \Delta r = \max(0.005,\; 0.10 \cdot r) \)

- **MinSize perturbation**
  - \( \Delta m = \max(5,\; 0.20 \cdot m) \)
  - Lower bound enforced: \( m_{\text{min}} \ge 0.8 \cdot m \)

- Candidate grid is constructed from all combinations.
- Baseline itself is removed.
- Deterministic subsampling limits to \( K \) perturbations (default \( K = 8 \)).

This ensures:

- Perturbations remain **local in relative parameter space**.
- No collapse to degenerate regimes (e.g. tiny minsize).
- Comparable stability semantics across scales.

---

## Sampling Strategy for ARI

To avoid memory explosions and integer overflows:

- Label rasters are optionally downsampled by a fixed factor (`sample_fact`) using nearest-neighbour resampling.
- ARI is computed on the resulting value vectors.
- No spatial coordinate indexing is used.
- No full raster-to-table conversions occur.

This keeps:

- Runtime bounded.
- Memory footprint stable.
- Numerical behaviour deterministic.

---

## Scale Construction

Two modes exist:

### A) CHM-driven (preferred when available)

If a high-resolution CHM/DSM exists:

1. Aggregate CHM to 10 m grid using maximum height.
2. Threshold canopy presence.
3. Label connected canopy patches.
4. Compute patch area distribution.
5. Convert patch areas to equivalent radii.
6. Use quantiles of radii as candidate spatial scales.
7. Derive `minsize` from patch area.
8. Derive `ranger` empirically from PCA feature-space distances.

This anchors segmentation scale in observed structural geometry.

### B) Fallback grid (current operational mode)

When no CHM is available, a fixed scale grid is used:

| scale_id | spatialr | ranger | minsize |
|-----------|-----------|---------|----------|
| s20m | 2 | 0.06 | 30 |
| s30m | 3 | 0.08 | 40 |
| s40m | 4 | 0.10 | 50 |
| s60m | 6 | 0.12 | 80 |
| s80m | 8 | 0.14 | 100 |
| s120m | 12 | 0.16 | 150 |

This grid is explicitly treated as provisional.

---

## Postcondition: Segment Size Filtering

Segmentation often produces many small fragments.

Instead of arbitrary thresholds:

1. Compute label size distribution (pixels per label).
2. Transform to log-scale.
3. Smooth (running median).
4. Detect maximum curvature (“knee”).
5. Use knee location as minimum pixel cutoff.

All labels smaller than this threshold are masked before polygonisation.

This enforces:

- Structural consistency of final objects.
- Suppression of unstable micro-fragments.
- Reproducible size logic.

---

## Interpretation of Results

Typical behaviour observed:

- Small scales → lower ARI_prev, high sensitivity to perturbations.
- Intermediate scales → increasing ARI_prev and stability.
- Large scales → ARI may saturate or slightly decrease due to over-smoothing.

Therefore:

- The optimum is usually a **stability plateau**, not necessarily the maximum scale.
- The composite score (mean − 0.5·sd) balances robustness and variability.

Final scale choice remains **operational**, not ontological.

---

## Limitations

- Stability does not imply semantic correctness.
- ARI is blind to spatial morphology and topology.
- Feature-space ranger calibration remains empirical.
- CHM-based scaling improves physical grounding but still depends on threshold choices.
- MeanShift segmentation remains sensitive to feature engineering.

These constraints are explicitly accepted.

---

## References

- Comaniciu, D., & Meer, P. (2002).  
  *Mean Shift: A robust approach toward feature space analysis.* IEEE TPAMI.

- Hubert, L., & Arabie, P. (1985).  
  *Comparing partitions.* Journal of Classification.

- Jolliffe, I., & Cadima, J. (2016).  
  *Principal component analysis: a review and recent developments.* Philosophical Transactions A.

---

