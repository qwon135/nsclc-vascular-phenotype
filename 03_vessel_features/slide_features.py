"""Reference cutoffs and slide-level vascular features."""
import numpy as np
import pandas as pd

import tumour_region as TR
import vessel_morphology as VM

MIN_ELIGIBLE_AREA = ...    # morphology-eligible vessel area
MIN_ELIGIBLE_VESSELS = ... # slide inclusion
LOWER, UPPER = ...         # reference percentiles for the shape cutoffs
CALIBRE_QUANTILE = ...
PHENOTYPES = [("Elongated", "eccentricity", "high"), ("Concave", "solidity", "low"),
              ("Irregular", "convexity", "low"), ("Round", "circularity", "high")]


def vessel_table(vessel_polygons, prob, tissue_frac, cell_px, specimen):
    region = TR.measurement_region(prob, specimen)
    rows = []
    for q in vessel_polygons:
        cx, cy = q[:, 0].mean(), q[:, 1].mean()
        if TR.in_region(region, cell_px, cx, cy) and not VM.is_tile_artefact(q):
            rows.append(VM.vessel_shape(q))
    area_mm2 = TR.region_area_mm2(region, tissue_frac, cell_px, VM.MPP)
    return pd.DataFrame(rows), area_mm2


def reference_cutoffs(reference_vessels):
    """Fixed percentiles of each shape metric in the outcome-independent reference set."""
    v = reference_vessels[reference_vessels.area >= MIN_ELIGIBLE_AREA]
    return {m: dict(lower=np.percentile(v[m], LOWER), upper=np.percentile(v[m], UPPER))
            for _, m, _ in PHENOTYPES}


def slide_features(vessels, area_mm2, cutoffs):
    v = vessels[vessels.area >= MIN_ELIGIBLE_AREA]
    if len(v) < MIN_ELIGIBLE_VESSELS:
        return None
    out = {"Vessel density": len(vessels) / area_mm2}
    for name, metric, side in PHENOTYPES:
        c = cutoffs[metric]
        out[name] = (v[metric] >= c["upper"]).mean() if side == "high" else (v[metric] <= c["lower"]).mean()
    out["Vessel calibre"] = v.diameter.quantile(CALIBRE_QUANTILE)
    return out
