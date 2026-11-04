"""Analytic squares exercise boundary distance independently of overlap."""

import numpy as np
import pytest

from aeroopt.geometry import validate_elements

SQUARE = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])


def codes(report):
    return {issue.code for issue in report.issues}


def test_analytic_separation_and_clearance():
    report = validate_elements(SQUARE, SQUARE + [1.25, 0.], min_clearance=.25)
    assert report.valid
    assert report.min_distance == pytest.approx(.25, abs=1e-14)
    assert report.overlap_area == 0.
    assert report.to_dict()["aerodynamic_gap"] is None
    rejected = validate_elements(SQUARE, SQUARE + [1.25, 0.], min_clearance=.251)
    assert codes(rejected) == {"insufficient_clearance"}


def test_clearance_tolerance_boundary():
    assert validate_elements(SQUARE, SQUARE + [1.25, 0.], min_clearance=.25 + 5e-11).valid
    assert not validate_elements(SQUARE, SQUARE + [1.25, 0.], min_clearance=.25 + 2e-10).valid


def test_overlap_crossing():
    report = validate_elements(SQUARE, SQUARE + [.5, .5])
    assert codes(report) == {"boundary_intersection", "solid_overlap"}
    assert report.overlap_area == pytest.approx(.25, abs=1e-14)
    assert report.min_distance == 0.


def test_containment_has_nonzero_boundary_distance_but_is_invalid():
    report = validate_elements(SQUARE, .25 * SQUARE + [.25, .25])
    assert codes(report) == {"solid_overlap"}
    assert report.min_distance == pytest.approx(.25, abs=1e-14)
    assert report.overlap_area == pytest.approx(.0625, abs=1e-14)


@pytest.mark.parametrize("offset", [[1., 0.], [1., 1.]])
def test_tangent_contact_rejected_even_at_zero_clearance(offset):
    report = validate_elements(SQUARE, SQUARE + offset)
    assert codes(report) == {"boundary_intersection"}
    assert report.overlap_area == 0.


def test_disjoint_near_contact_rejected():
    report = validate_elements(SQUARE, SQUARE + [1. + 5e-11, 0.])
    assert codes(report) == {"near_contact"}


@pytest.mark.parametrize("points,expected", [
    ([[0, 0], [1, 1], [0, 1], [1, 0]], "self_intersection"),
    ([[0, 0], [1, 0], [2, 0]], "degenerate_contour"),
    ([[0, 0], [1, 0]], "degenerate_contour"),
    ([[0, 0], [1, 0], [1, 0], [0, 1]], "duplicate_points"),
    ([[0, 0], [1, 0], [np.inf, 1]], "non_finite"),
    ([[0, 0], [0, 0], [0, 0]], "degenerate_contour"),
])
def test_invalid_contours_report_structured_reasons(points, expected):
    report = validate_elements(points, SQUARE + [3., 0.])
    assert expected in codes(report)
    assert report.min_distance is None
    assert report.overlap_area is None
    assert report.issues[0].element == "main"
    assert report.to_dict()["issues"][0]["message"]


def test_open_contour_uses_closing_segment():
    # The closing edge x=0 is the nearest part, although absent from the input rows.
    report = validate_elements(SQUARE, SQUARE + [-1.2, 0.])
    assert report.valid
    assert report.min_distance == pytest.approx(.2, abs=1e-14)


@pytest.mark.parametrize("clearance", [-.1, np.inf, np.nan])
def test_invalid_clearance(clearance):
    with pytest.raises(ValueError):
        validate_elements(SQUARE, SQUARE + [3., 0.], min_clearance=clearance)

