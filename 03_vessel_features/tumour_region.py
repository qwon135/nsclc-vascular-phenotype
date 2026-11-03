"""Measurement region from a tumour probability map (grid of cells over the slide)."""
import numpy as np
from scipy import ndimage as ndi

TUMOUR_PROB = ...      # tumour cell threshold
CORE_PERCENTILE = ...  # resection: deeper part of the tumour, see Methods


def tumour_cells(prob):
    return np.isfinite(prob) & (np.nan_to_num(prob) >= TUMOUR_PROB)


def tumour_core(prob):
    """Tumour cells whose distance to the nearest non-tumour cell is >= the slide-level percentile."""
    tum = tumour_cells(prob)
    dist = ndi.distance_transform_edt(tum)
    return tum & (dist >= np.percentile(dist[tum], CORE_PERCENTILE))


def measurement_region(prob, specimen):
    return tumour_core(prob) if specimen == "resection" else tumour_cells(prob)


def region_area_mm2(mask, tissue_frac, cell_px, mpp):
    return float((np.where(mask, tissue_frac, 0) * (cell_px * mpp / 1000) ** 2).sum())


def in_region(mask, cell_px, x, y):
    return mask[int(y // cell_px), int(x // cell_px)]
