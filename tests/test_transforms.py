"""Analytic lengths, area, rotation direction, and immutable configuration state."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from aeroopt.geometry import FlapParameters, transform_flap


@pytest.mark.parametrize("angle", [-90., 0., 25., 90., 370.])
def test_rigid_distances_and_scaled_area(diamond, angle):
    from shapely.geometry import Polygon

    p = FlapParameters(.3, (1.1, -.2), angle, .02, -.03)
    transformed = transform_flap(diamond.normalized_points, p)
    before = np.linalg.norm(diamond.normalized_points[:, None] - diamond.normalized_points[None, :], axis=2)
    after = np.linalg.norm(transformed[:, None] - transformed[None, :], axis=2)
    np.testing.assert_allclose(after, .3 * before, atol=5e-15, rtol=0)
    assert Polygon(transformed).area == pytest.approx(.1 * .3**2, abs=1e-14)
    np.testing.assert_allclose(transformed[diamond.leading_edge_index], [1.12, -.23], atol=1e-14)


def test_positive_deflection_is_clockwise():
    points = [[0., 0.], [1., 0.], [0., 1.]]
    actual = transform_flap(points, FlapParameters(.5, (0., 0.), 90.))
    np.testing.assert_allclose(actual, [[0., 0.], [0., -.5], [.5, 0.]], atol=1e-15, rtol=0)


def test_absolute_configurations_are_reproducible(separated):
    original = separated.flap.original_points.copy()
    first = separated.configure(flap_deflection_deg=25., flap_dx=.02, flap_dy=-.03)
    other = first.configure(flap_deflection_deg=80., flap_dx=-.1)
    final = other.configure(flap_deflection_deg=25., flap_dx=.02, flap_dy=-.03)
    np.testing.assert_array_equal(first.flap_points, final.flap_points)
    np.testing.assert_array_equal(first.flap_points, first.flap_points)
    np.testing.assert_array_equal(separated.flap.original_points, original)
    assert separated.parameters.deflection_deg == 0.
    assert first.parameters.to_dict()["dx"] == .02
    assert other.parameters.dy == 0.
    with pytest.raises(FrozenInstanceError):
        first.parameters.dx = 8.


@pytest.mark.parametrize("kwargs", [
    {"chord_ratio": 0.}, {"chord_ratio": -1.}, {"chord_ratio": np.inf},
    {"deflection_deg": np.nan}, {"dx": np.inf}, {"dy": -np.inf},
    {"reference_position": (0., np.nan)}, {"reference_position": (0.,)},
    {"reference_position": (1e308, 0.), "dx": 1e308},
])
def test_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        FlapParameters(**kwargs)


def test_transform_rejects_overflow():
    with pytest.raises(ValueError, match="unrepresentable"):
        transform_flap([[1e308, 0.]], FlapParameters(10., (0., 0.)))

