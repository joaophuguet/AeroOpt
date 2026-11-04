"""Independent published ordinates, exact integrals and algebraic NACA identities."""

from decimal import Decimal, localcontext

import numpy as np
import pytest
from shapely.geometry import Polygon

from aeroopt.geometry import Airfoil, naca4_points, naca4_surfaces
from examples._naca_reference import continuous_metrics, continuous_slot_distance


@pytest.mark.parametrize("code", ["0012", "4412"])
@pytest.mark.parametrize("closed", [False, True])
@pytest.mark.parametrize("count", [41, 81, 161, 321])
def test_generated_shape_count_units_references_and_area(code, closed, count):
    points = naca4_points(code, count, closed_te=closed)
    foil = Airfoil.from_points(points, leading_edge_index=count // 2)
    assert len(points) == count
    assert foil.chord == pytest.approx(1., abs=1e-14)
    np.testing.assert_allclose(foil.leading_edge, [0., 0.], atol=1e-14)
    np.testing.assert_allclose(foil.trailing_edge, [1., 0.], atol=1e-14)
    assert foil.is_te_closed == closed
    reference = continuous_metrics(code, closed)
    # Error bounds concern THIS discretization study, not acceptance tolerances.
    area_error = abs(foil.normalized_area - reference["area"])
    assert area_error < .14 / ((count - 1) / 2)**2
    perimeter = Polygon(points).length
    assert abs(perimeter - reference["perimeter"]) < .5 / ((count - 1) / 2)**2


def test_symmetric_decimal_ordinate_independent_of_generator():
    with localcontext() as ctx:
        ctx.prec = 50
        x = Decimal('.3')
        half = Decimal('.6') * (Decimal('.2969') * x.sqrt() - Decimal('.126') * x
                               - Decimal('.3516') * x*x + Decimal('.2843') * x**3
                               - Decimal('.1015') * x**4)
    upper, lower = naca4_surfaces("0012", [.3])
    assert upper[0, 1] == pytest.approx(float(half), abs=2e-17)
    assert lower[0, 1] == pytest.approx(-float(half), abs=2e-17)


def test_published_4412_ordinates_at_zero_slope_station():
    # UIUC naca4412.dat at x=.4: upper=.0980, lower=-.0180 (4-place table).
    upper, lower = naca4_surfaces("4412", [.4])
    np.testing.assert_allclose(upper[0], [.4, .0980], atol=5e-5, rtol=0)
    np.testing.assert_allclose(lower[0], [.4, -.0180], atol=5e-5, rtol=0)
    np.testing.assert_allclose((upper + lower)[0] / 2, [.4, .04], atol=1e-15)


def test_normally_offset_camber_geometry():
    upper, lower = naca4_surfaces("4412", [.25])
    camber, slope = .034375, .075
    mean = (upper + lower)[0] / 2
    np.testing.assert_allclose(mean, [.25, camber], atol=1e-15)
    assert np.dot(upper[0] - lower[0], [1., slope]) == pytest.approx(0., abs=1e-15)


def test_symmetric_continuous_integral_is_exact_known_area():
    assert continuous_metrics("0012")["area"] == pytest.approx(.08221, abs=5e-17)
    assert continuous_metrics("0012", True)["area"] == pytest.approx(.081706, abs=5e-17)
    assert continuous_metrics("0012")["quadrature_vs_exact_area"] < 1e-15


def test_continuous_slot_reference_has_cap_projection():
    assert continuous_slot_distance("0012")["distance"] == .01
    cambered = continuous_slot_distance("4412")
    assert .00990 < cambered["distance"] < .00992
    assert abs(.3 * cambered["flap_y_at_min_x"]) < .00126
    assert cambered["derivative_residual"] < 1e-15


def test_nominal_le_and_farthest_sample_are_distinct_references():
    points = naca4_points("4412", 161)
    auto = Airfoil.from_points(points)
    nominal = Airfoil.from_points(points, leading_edge_index=80)
    assert auto.leading_edge_index == 79
    assert nominal.leading_edge_index == 80
    # Record a real limitation; do not "fix" the default by silently changing conventions.
    assert .0034 < np.linalg.norm(auto.leading_edge - nominal.leading_edge) < .0036
    assert .00029 < auto.chord - nominal.chord < .00031
    np.testing.assert_array_equal(auto.original_points, nominal.original_points)


@pytest.mark.parametrize("code", ["", "12", "abc4", "１２３４", "4012", "0112", "0000", 12])
def test_bad_designations(code):
    with pytest.raises(ValueError):
        naca4_points(code)


@pytest.mark.parametrize("count", [0, 4, 40, 4.5, True])
def test_bad_point_counts(count):
    with pytest.raises(ValueError):
        naca4_points("0012", count)


@pytest.mark.parametrize("stations", [[], [-.01], [1.01], [np.nan], [[.5]]])
def test_bad_stations(stations):
    with pytest.raises(ValueError):
        naca4_surfaces("0012", stations)


def test_bad_generation_options():
    with pytest.raises(ValueError):
        naca4_points("0012", spacing="unknown")
    with pytest.raises(ValueError):
        naca4_points("0012", closed_te=1)
