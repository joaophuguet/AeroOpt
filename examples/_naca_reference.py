"""Independent scalar continuous-geometry references for the validation study.

This module does NOT call the coordinate generator or derive references by
mesh extrapolation. Expanded PDAS coefficients, analytic normal-offset speed,
and Gauss-Legendre quadrature after x=u² provide independent area/perimeter.
"""

from functools import lru_cache
from math import hypot, sqrt

import numpy as np


def _local(code: str, u: float, closed: bool) -> tuple[float, float, float, float, float]:
    """Return thickness, d(thickness)/du, camber, slope, curvature numerator."""
    m, p, t = int(code[0]) * .01, int(code[1]) * .1, int(code[2:]) * .01
    x = u * u
    last = -.518 if closed else -.5075
    h = t * (1.4845 * u - .630 * u**2 - 1.758 * u**4 + 1.4215 * u**6 + last * u**8)
    dh = t * (1.4845 - 1.260 * u - 7.032 * u**3 + 8.529 * u**5 + 8 * last * u**7)
    if m == 0:
        return h, dh, 0., 0., 0.
    if x <= p:
        yc = m * (2 * x / p - (x / p)**2)
        g, gp = 2 * m * (p - x) / p**2, -2 * m / p**2
    else:
        yc = m * (1 - 2 * p + 2 * p * x - x*x) / (1 - p)**2
        g, gp = 2 * m * (p - x) / (1 - p)**2, -2 * m / (1 - p)**2
    return h, dh, yc, g, gp


@lru_cache(maxsize=None)
def _integrate(code: str, closed: bool, order: int) -> tuple[float, float]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    breaks = [0., 1.] if code[0] == "0" else [0., sqrt(int(code[1]) / 10), 1.]
    area, perimeter = 0., 0.
    for left, right in zip(breaks[:-1], breaks[1:], strict=True):
        for node, weight in zip(nodes, weights, strict=True):
            u = float(left + (node + 1) * (right - left) / 2)
            w = float(weight * (right - left) / 2)
            h, dh, _, g, gp = _local(code, u, closed)
            s = hypot(1., g)
            area += w * 4 * u * h * s
            perimeter += w * (hypot(2 * u * (s - h * gp / s**2), dh)
                              + hypot(2 * u * (s + h * gp / s**2), dh))
    perimeter += 2 * _local(code, 1., closed)[0] if not closed else 0.
    return area, perimeter


def continuous_metrics(code: str, closed_te: bool = False) -> dict:
    """Independent continuous area/length, with a quadrature-order diagnostic.

    The symmetric area is also evaluated by an exact polynomial integral.
    Order disagreement is an error indicator, not a rigorous error bound.
    """
    if code not in ("0012", "4412"):
        raise ValueError("independent references are validated only for 0012/4412")
    a128, p128 = _integrate(code, closed_te, 128)
    area, perimeter = _integrate(code, closed_te, 256)
    exact = None
    if code == "0012":
        # Integral 2*t*(A sqrt(x)+B*x+C*x²+D*x³+E*x⁴), from 0 to 1.
        exact = 2 * .12 * (1.4845 * 2/3 - .630/2 - 1.758/3 + 1.4215/4
                           - (.518 if closed_te else .5075)/5)
    return {"area": exact if exact is not None else area, "perimeter": perimeter,
            "area_reference": "exact integral" if exact is not None else "independent continuous quadrature",
            "quadrature_area_128_256_difference": abs(area - a128),
            "quadrature_perimeter_128_256_difference": abs(perimeter - p128),
            "quadrature_vs_exact_area": None if exact is None else abs(area - exact)}


def continuous_slot_distance(code: str, closed_te: bool = False) -> dict:
    """Exact cap-to-nose construction for this study, solved in continuous u.

    Main=finite-TE 0012 with x<=1, flap ratio=.3 at (1.01,0), no deflection.
    The global minimum flap x projects INSIDE the main TE cap. Thus the global
    boundary distance equals min(flap x)-1, without curve-to-curve meshing.
    For 4412 solve dx_upper/du=0 in the nose bracket [0,.25] by bisection.
    """
    if code not in ("0012", "4412"):
        raise ValueError("slot construction is specified only for 0012/4412")
    left, right = 0., .25

    def dx(u: float) -> float:
        h, dh, _, g, gp = _local(code, u, closed_te)
        s = hypot(1., g)
        return 2 * u - dh * g / s - h * (2 * u * gp) / s**3

    if code == "0012":
        u = 0.
    else:
        if not dx(left) < 0 < dx(right):
            raise ValueError("continuous nose derivative root is not bracketed")
        for _ in range(64):
            middle = (left + right) / 2
            if dx(middle) < 0:
                left = middle
            else:
                right = middle
        u = (left + right) / 2
    h, _, yc, g, _ = _local(code, u, closed_te)
    minimum_x, y = u*u - h*g/hypot(1., g), yc + h/hypot(1., g)
    if abs(.3 * y) >= .00126:
        raise ValueError("minimum-x point does not project inside the finite main TE cap")
    return {"distance": .01 + .3 * minimum_x, "flap_min_x": minimum_x,
            "flap_y_at_min_x": y, "nominal_x_at_min_x": u*u,
            "derivative_residual": abs(dx(u)),
            "reference": "continuous minimum flap x with projection inside finite main TE cap"}
