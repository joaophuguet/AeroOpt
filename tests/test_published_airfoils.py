"""Optional pinned published coordinates; no network or unlicensed bundling in tests."""

import hashlib

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, read_dat, validate_elements


def test_published_points_references_and_original_bytes(published_reference):
    path, entry = published_reference
    raw = path.read_bytes()
    foil = Airfoil.from_dat(path, leading_edge_index=entry["reference_leading_edge_index"])
    auto = Airfoil.from_dat(path)
    assert len(foil.original_points) == entry["point_count"]
    assert foil.leading_edge_index == auto.leading_edge_index  # Observed, not assumed for all airfoils.
    np.testing.assert_allclose(foil.normalized_points[foil.leading_edge_index], [0., 0.], atol=1e-14)
    np.testing.assert_allclose((foil.normalized_points[0]+foil.normalized_points[-1])/2, [1., 0.], atol=1e-14)
    assert foil.normalized_area > 0
    if entry["name"] == "S1223":
        np.testing.assert_array_equal(foil.leading_edge, [.00005, .00178])
    else:
        np.testing.assert_array_equal(foil.leading_edge, [0., 0.])
        assert not foil.is_te_closed
    assert path.read_bytes() == raw
    assert hashlib.sha256(raw).hexdigest() == entry["sha256"]


@pytest.mark.parametrize("scale,angle", [(1e-6, -80.), (1., 135.), (1e6, 40.)])
def test_published_normalization_under_frame_changes(published_reference, scale, angle):
    from shapely.geometry import Polygon

    path, entry = published_reference
    foil = Airfoil.from_dat(path, leading_edge_index=entry["reference_leading_edge_index"])
    theta = np.deg2rad(angle)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    original = scale * (foil.original_points @ rotation.T + [2.3, -1.7])
    moved = Airfoil.from_points(original, leading_edge_index=foil.leading_edge_index)
    np.testing.assert_allclose(moved.normalized_points, foil.normalized_points, atol=2e-14, rtol=0)
    np.testing.assert_array_equal(moved.original_points, original)
    a = foil.normalized_points
    b = moved.normalized_points
    np.testing.assert_allclose(np.linalg.norm(a[:, None]-a[None,:],axis=2),
                                np.linalg.norm(b[:, None]-b[None,:],axis=2), atol=2e-14, rtol=0)
    assert moved.area / scale**2 == pytest.approx(foil.area, abs=1e-14)
    assert Polygon(a).exterior.is_ccw == Polygon(b).exterior.is_ccw


@pytest.mark.parametrize("angle,ratio", [(0., .3), (25., .2), (-20., .4)])
def test_published_pair_export_roundtrip(published_reference, angle, ratio, tmp_path):
    path, entry = published_reference
    foil = Airfoil.from_dat(path, leading_edge_index=entry["reference_leading_edge_index"])
    geometry = MultiElementAirfoil(foil, foil, FlapParameters(ratio, (1.1, -.04), angle))
    report = geometry.validate(min_clearance=.002)
    assert report.valid
    target = geometry.export(tmp_path / "published", configuration_id=entry["name"], min_clearance=.002)
    _, main = read_dat(target / "main.dat")
    _, flap = read_dat(target / "flap.dat")
    np.testing.assert_array_equal(main, geometry.main_points)
    np.testing.assert_array_equal(flap, geometry.flap_points)
    assert validate_elements(main, flap, min_clearance=.002).min_distance == pytest.approx(report.min_distance, abs=1e-14)
