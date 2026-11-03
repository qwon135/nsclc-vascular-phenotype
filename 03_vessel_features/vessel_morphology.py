"""Per-vessel morphology from predicted vessel polygons (pixel coordinates)."""
import numpy as np
from shapely.geometry import MultiPoint
from shapely.geometry import Polygon as ShapelyPolygon

MPP = ...              # um per pixel of the analysis grid
SIMPLIFY_TOL = ...     # px, removes the 1-px staircase of the mask contour
ARTEFACT_RULES = ...   # area limits and right-angle edge length, see Methods


def polygon_area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def perimeter(p):
    d = np.diff(np.vstack([p, p[:1]]), axis=0)
    return float(np.hypot(d[:, 0], d[:, 1]).sum())


def area_moments(p):
    """Area, centroid and central second moments (mxx, myy, mxy) of a closed polygon (Green's theorem)."""
    x, y = p[:, 0], p[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    a = cr.sum() / 2
    cx = ((x + x1) * cr).sum() / (6 * a)
    cy = ((y + y1) * cr).sum() / (6 * a)
    ixx = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12
    iyy = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12
    ixy = ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr).sum() / 24
    return abs(a), (cx, cy), (ixx / a - cx * cx, iyy / a - cy * cy, ixy / a - cx * cy)


def simplify(p, tol=SIMPLIFY_TOL):
    g = ShapelyPolygon(p).buffer(0)
    if g.geom_type != "Polygon":
        g = max(g.geoms, key=lambda a: a.area)
    return np.asarray(g.simplify(tol).exterior.coords)[:-1]


def vessel_shape(poly):
    p = simplify(np.asarray(poly, dtype=np.float64))
    area, _, (mxx, myy, mxy) = area_moments(p)
    per = perimeter(p)
    tr, det = mxx + myy, mxx * myy - mxy * mxy
    disc = np.sqrt(max(tr * tr / 4 - det, 0.0))
    l1, l2 = tr / 2 + disc, max(tr / 2 - disc, 0.0)
    hull = MultiPoint(p).convex_hull
    return dict(area=area,
                diameter=2 * np.sqrt(area / np.pi) * MPP,
                eccentricity=np.sqrt(1 - l2 / l1),
                solidity=area / hull.area,
                convexity=hull.length / per,
                circularity=4 * np.pi * area / per ** 2)


def is_tile_artefact(poly):
    """Straight edges / right angles from patch stitching: very large polygons, or large polygons
    with a right angle between two long edges (thresholds in ARTEFACT_RULES)."""
    ...
