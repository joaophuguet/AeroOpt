"""Reproducible stress cases; expectations concern geometry only."""

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points


def validation_cases(*, symmetric: Airfoil | None = None,
                     cambered: Airfoil | None = None) -> list[tuple[str, MultiElementAirfoil, float, bool, set[str]]]:
    """Include analytic cap contact/distance and true containment, not only crossing."""
    if symmetric is None:
        symmetric = Airfoil.from_points(naca4_points("0012", 161), leading_edge_index=80,
                                       name="ANALYTIC_NACA_0012_161_FINITE_TE")
    if cambered is None:
        cambered = Airfoil.from_points(naca4_points("4412", 81), leading_edge_index=40,
                                      name="ANALYTIC_NACA_4412_81_FINITE_TE")
    definitions = [
        ("zero", symmetric, cambered, .3, (1.02, 0.), 0., .002, True, set()),
        ("positive", cambered, symmetric, .3, (1.02, -.04), 25., .002, True, set()),
        ("negative", cambered, symmetric, .2, (1.02, .04), -20., .002, True, set()),
        ("insufficient", symmetric, symmetric, .3, (1.00005, 0.), 0., .002, False,
         {"insufficient_clearance"}),
        ("contact", symmetric, symmetric, .3, (1., 0.), 0., 0., False,
         {"boundary_intersection"}),
        ("crossing", symmetric, cambered, .3, (.9, -.01), 10., 0., False,
         {"boundary_intersection", "solid_overlap"}),
        ("contained", symmetric, symmetric, .03, (.45, 0.), 0., 0., False, {"solid_overlap"}),
        ("above_tolerance", symmetric, symmetric, .3, (1. + 2e-10, 0.), 0., 0., True, set()),
        ("below_tolerance", symmetric, symmetric, .3, (1. + 5e-11, 0.), 0., 0., False,
         {"near_contact"}),
        ("small_resolved_flap", symmetric, symmetric, 1e-7, (1.02, 0.), 0., .01, True, set()),
    ]
    return [(label, MultiElementAirfoil(main, flap, FlapParameters(ratio, nominal, angle)),
             clearance, valid, codes)
            for label, main, flap, ratio, nominal, angle, clearance, valid, codes in definitions]
