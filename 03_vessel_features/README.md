# Vessel features

Slide-level vascular features from the predicted vessel polygons and a tumour probability map
(patch-level tumour classifier; not included).

## Per-vessel morphology

Shape metrics are measured on each vessel polygon after removing the pixel staircase of the mask
contour (polygon simplification). With area $A$, perimeter $P$, convex-hull area $A_H$ and
convex-hull perimeter $P_H$:

| Metric | Formula |
|---|---|
| Eccentricity | $\sqrt{1 - \lambda_2 / \lambda_1}$, $\lambda_1 \ge \lambda_2$ eigenvalues of the second-order central area moments of the polygon (0 = circle) |
| Solidity | $A / A_H$ |
| Convexity | $P_H / P$ |
| Circularity | $4\pi A / P^2$ |
| Equivalent diameter | $2\sqrt{A / \pi}$ |

Polygons produced by patch stitching (very large polygons, or large polygons with a long right
angle) are removed as tile artefacts.

## Measurement region

Tumour cells of the probability map are thresholded. For resection specimens, the measurement
region is the deeper part of the tumour (Euclidean distance to the tumour boundary at or above a
slide-level percentile); biopsy specimens use the whole tumour region. A vessel belongs to the
region when its centroid lies in a region cell.

## Slide-level features

Shape cutoffs are fixed percentiles of each metric in an outcome-independent reference set of
slides and are applied unchanged to every cohort. For the morphology-eligible vessels of a slide:

| Feature | Definition |
|---|---|
| Elongated | proportion with eccentricity $\ge$ upper cutoff |
| Concave | proportion with solidity $\le$ lower cutoff |
| Irregular | proportion with convexity $\le$ lower cutoff |
| Round | proportion with circularity $\ge$ upper cutoff |
| Vessel calibre | upper percentile of the equivalent diameter (µm) |
| Vessel density | vessels per mm² of the measurement region |

Area thresholds, percentiles and slide inclusion rules are given in the Methods.

## Files

- `vessel_morphology.py` — shape metrics and tile-artefact filter
- `tumour_region.py` — measurement region from the tumour map
- `slide_features.py` — reference cutoffs and slide-level features
