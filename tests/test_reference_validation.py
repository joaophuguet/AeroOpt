"""Analytic airfoils under frame changes and physically distinct pair relations."""

import json

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, naca4_points, read_dat, validate_elements
from examples._naca_reference import continuous_metrics, continuous_slot_distance
from examples._validation_cases import validation_cases


@pytest.mark.parametrize("code", ["0012", "4412"])
@pytest.mark.parametrize("count", [41, 321])
@pytest.mark.parametrize("scale,angle", [(1e-9, -70.), (1., 130.), (1e9, 35.)])
def test_realistic_shape_normalization_preserves_metric_and_orientation(code, count, scale, angle):
    p = naca4_points(code, count)
    theta = np.deg2rad(angle)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    original = scale * (p @ rotation.T + [2.3, -1.7])
    untouched = original.copy()
    normalized = Airfoil.from_points(original, leading_edge_index=count // 2)
    np.testing.assert_array_equal(original, untouched)
    np.testing.assert_array_equal(normalized.original_points, original)
    np.testing.assert_allclose(normalized.normalized_points, p, atol=3e-14, rtol=0)
    before = np.linalg.norm(p[:, None] - p[None, :], axis=2)
    after = np.linalg.norm(normalized.normalized_points[:, None] - normalized.normalized_points[None, :], axis=2)
    np.testing.assert_allclose(after, before, atol=3e-14, rtol=0)
    area = Airfoil.from_points(p, leading_edge_index=count // 2).area
    assert normalized.area / scale**2 == pytest.approx(area, abs=2e-15)
    # Rotation must not turn the upper surface below the physical mean line.
    assert normalized.normalized_points[count // 4, 1] > 0


@pytest.mark.parametrize("case", validation_cases(), ids=lambda case: case[0])
def test_multielement_realistic_stress_cases(case, tmp_path):
    label, geometry, clearance, expected, required = case
    report = geometry.validate(min_clearance=clearance)
    assert report.valid == expected, report.to_dict()
    assert required.issubset({issue.code for issue in report.issues})
    if label == "contained":
        assert report.min_distance > 0 and report.overlap_area > 0
        assert not any(issue.code == "boundary_intersection" for issue in report.issues)
    if label == "contact":
        assert report.min_distance == 0 and report.overlap_area == 0
    if label == "insufficient":
        assert report.min_distance == pytest.approx(5e-5, abs=1e-14)
    if expected:
        path = geometry.export(tmp_path / label, configuration_id=label, min_clearance=clearance)
        _, main = read_dat(path / "main.dat")
        _, flap = read_dat(path / "flap.dat")
        np.testing.assert_array_equal(main, geometry.main_points)
        np.testing.assert_array_equal(flap, geometry.flap_points)
        after = validate_elements(main, flap, min_clearance=clearance)
        metadata = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
        assert metadata["configuration_id"] == label
        assert metadata["parameters"] == json.loads(json.dumps(geometry.parameters.to_dict()))
        assert metadata["files"] == {"main": "main.dat", "flap": "flap.dat"}
        assert after.min_distance == pytest.approx(report.min_distance, abs=1e-14)
        assert metadata["validation"] == json.loads(json.dumps(after.to_dict()))
        assert not geometry.main.is_te_closed  # Finite-TE round trips explicitly exercised.


@pytest.mark.parametrize("code", ["0012", "4412"])
@pytest.mark.parametrize("closed", [False, True])
def test_discretization_against_continuous_area_perimeter_and_distance(code, closed):
    from shapely.geometry import Polygon
    from aeroopt.geometry import FlapParameters, MultiElementAirfoil

    reference = continuous_metrics(code, closed)
    errors = []
    for count in [41, 81, 161, 321]:
        points = naca4_points(code, count, closed_te=closed)
        foil = Airfoil.from_points(points, leading_edge_index=count // 2)
        main = Airfoil.from_points(naca4_points("0012", count), leading_edge_index=count // 2)
        geometry = MultiElementAirfoil(main, foil, FlapParameters(.3, (1.01, 0.)))
        report = geometry.validate(min_clearance=.002)
        assert report.valid
        distance_error = abs(report.min_distance - continuous_slot_distance(code, closed)["distance"])
        assert distance_error < 9.1e-5
        errors.append((abs(foil.area-reference["area"]), abs(Polygon(points).length-reference["perimeter"])))
    # Absolute independent-reference bounds, in addition to convergence ratios.
    assert errors[-1][0] < 5.5e-6
    assert errors[-1][1] < 4e-5
    assert errors[-1][0] < errors[0][0] / 40
    assert errors[-1][1] < errors[0][1] / 40
