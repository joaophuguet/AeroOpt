"""Analytically generated NACA four-digit coordinates, not measured datasets.

Equations: PDAS naca456thick4.html and naca456mean2.html, with thickness
normal to the piecewise-parabolic mean line. No imported contour is modified.
"""

import numpy as np
from numpy.typing import ArrayLike

from ._numeric import FloatPoints, immutable_points


def _parameters(designation: str) -> tuple[float, float, float]:
    if (not isinstance(designation, str) or len(designation) != 4
            or not designation.isascii() or not designation.isdigit()):
        raise ValueError("designation must be exactly four ASCII digits, e.g. '4412'")
    m, p, t = int(designation[0]) / 100, int(designation[1]) / 10, int(designation[2:]) / 100
    if t <= 0 or (m > 0 and p == 0) or (m == 0 and p != 0):
        raise ValueError("positive thickness and consistent camber/location digits are required")
    return m, p, t


def naca4_surfaces(designation: str, x: ArrayLike, *,
                   closed_te: bool = False) -> tuple[FloatPoints, FloatPoints]:
    """Return upper/lower coordinates at mean-line abscissas x in [0,1].

    x is a one-dimensional set of stations in the nominal unit-chord frame.
    Standard finite TE uses coefficient -0.1015; closed_te explicitly selects
    the modified -0.1036 coefficient. The nominal LE reference is (0,0), which
    need not be the minimum-x point of a cambered thickness-offset contour.
    """
    m, p, t = _parameters(designation)
    stations = np.asarray(x, dtype=np.float64)
    if (stations.ndim != 1 or len(stations) == 0 or not np.isfinite(stations).all()
            or np.any((stations < 0) | (stations > 1))):
        raise ValueError("x must be a nonempty finite 1D array within [0,1]")
    if not isinstance(closed_te, bool):
        raise ValueError("closed_te must be a bool")
    coefficient = -0.1036 if closed_te else -0.1015
    thickness = 5 * t * (0.2969 * np.sqrt(stations) - 0.1260 * stations
                        - 0.3516 * stations**2 + 0.2843 * stations**3
                        + coefficient * stations**4)
    if closed_te:
        # The polynomial is mathematically zero here; enforce its exact endpoint.
        thickness = np.where(stations == 1, 0., thickness)
    if m == 0:
        camber = np.zeros_like(stations)
        slope = np.zeros_like(stations)
    else:
        front = stations <= p
        camber = np.where(front, m / p**2 * (2 * p * stations - stations**2),
                          m / (1 - p)**2 * (1 - 2 * p + 2 * p * stations - stations**2))
        slope = np.where(front, 2 * m / p**2 * (p - stations),
                         2 * m / (1 - p)**2 * (p - stations))
    theta = np.arctan(slope)
    offset = np.column_stack((-thickness * np.sin(theta), thickness * np.cos(theta)))
    mean = np.column_stack((stations, camber))
    return immutable_points(mean + offset), immutable_points(mean - offset)


def naca4_points(designation: str, n_points: int = 161, *,
                 closed_te: bool = False, spacing: str = "cosine") -> FloatPoints:
    """Generate exactly odd n_points in upper TE -> LE -> lower TE order.

    Both TE endpoints are retained, and the nominal LE is included once at
    index (n_points-1)//2. For a closed TE the last row repeats the first.
    """
    if isinstance(n_points, bool) or not isinstance(n_points, int) or n_points < 5 or n_points % 2 == 0:
        raise ValueError("n_points must be an odd integer >= 5")
    if spacing not in ("cosine", "uniform"):
        raise ValueError("spacing must be 'cosine' or 'uniform'")
    count = (n_points + 1) // 2
    x = ((1 - np.cos(np.linspace(0, np.pi, count))) / 2
         if spacing == "cosine" else np.linspace(0, 1, count))
    upper, lower = naca4_surfaces(designation, x, closed_te=closed_te)
    return immutable_points(np.vstack((upper[::-1], lower[1:])))
