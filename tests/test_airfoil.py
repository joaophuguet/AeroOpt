"""DAT contracts and analytical normalization/area checks."""

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, AirfoilFormatError, read_dat


@pytest.mark.parametrize("header", ["", "NACA 0012\n", "# synthetic\n", "\ufeffDemo profile\n"])
def test_valid_dat(tmp_path, header):
    path = tmp_path / "foil.dat"
    path.write_text(header + "\n1e0 0\n 5e-1\t1E-1\n0 0\n.5 -.1\n1 0\n", encoding="utf-8")
    foil = Airfoil.from_dat(path)
    assert foil.chord == pytest.approx(1., abs=1e-14)
    assert foil.area == pytest.approx(.1, abs=1e-14)
    assert foil.source == path.resolve()
    assert foil.is_te_closed


@pytest.mark.parametrize("text,match", [
    ("", "three coordinate"),
    ("title\n1 0\n0 0\n", "three coordinate"),
    ("title\nanother header\n1 0\n", "two numeric"),
    ("title\n1 0 2\n", "two numeric"),
    ("1 0\n0 broken\n", "two numeric"),
    ("1,0 0,0\n", "two numeric"),
    ("1 0\n0 0 # comment\n", "two numeric"),
    ("Demo\n3 3\n1 0\n", "count table"),
    ("NaN 0\n0 0\n1 0\n", "non-finite"),
    ("1 0\n0 inf\n.5 .1\n", "non-finite"),
    ("1 0\n0 -1e999\n.5 .1\n", "non-finite"),
])
def test_bad_dat_has_explanation(tmp_path, text, match):
    path = tmp_path / "bad.dat"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(AirfoilFormatError, match=match):
        Airfoil.from_dat(path)


@pytest.mark.parametrize("scale", [1e-6, .2, 1., 7., 1e6])
@pytest.mark.parametrize("angle", [0., 40., 120., 230.])
def test_normalization_is_uniform_rigid(diamond_points, scale, angle):
    theta = np.deg2rad(angle)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    # Translation proportional to scale avoids artificially destroying input precision.
    original = scale * (diamond_points @ rotation.T + [3.1, -1.7])
    foil = Airfoil.from_points(original)
    np.testing.assert_allclose(foil.normalized_points, diamond_points, atol=2e-14, rtol=0)
    assert foil.chord == pytest.approx(scale, rel=2e-14)
    assert foil.area == pytest.approx(.1 * scale**2, rel=2e-14)
    assert foil.leading_edge_index == 2
    # A rotated foil's minimum global x need not be its leading reference.
    if angle == 120.:
        assert np.argmin(original[:, 0]) != foil.leading_edge_index


def test_open_te_cap_and_midpoint():
    points = np.array([[1., .02], [.5, .1], [0., 0.], [.5, -.1], [1., -.02]])
    foil = Airfoil.from_points(points)
    assert not foil.is_te_closed
    assert foil.area == pytest.approx(.11, abs=1e-14)
    np.testing.assert_array_equal(foil.trailing_edge, [1., 0.])
    np.testing.assert_array_equal(foil.original_points, points)
    assert len(foil.normalized_points) == len(points)


def test_explicit_le_reference_for_undersampled_nose():
    # Here the farthest sample proxy differs slightly from the selected nose.
    points = [[1., .01], [.3, .1], [0., .02], [.001, -.1], [.5, -.08], [1., -.01]]
    auto = Airfoil.from_points(points)
    manual = Airfoil.from_points(points, leading_edge_index=2)
    assert auto.leading_edge_index == 3
    assert manual.leading_edge_index == 2
    assert manual.reference_method == "explicit_leading_edge_index"
    np.testing.assert_allclose(manual.normalized_points[2], [0., 0.], atol=1e-14)


@pytest.mark.parametrize("index", [0, 4, -1, 1.5, True])
def test_bad_le_indices(diamond_points, index):
    with pytest.raises(AirfoilFormatError, match="leading"):
        Airfoil.from_points(diamond_points, leading_edge_index=index)


@pytest.mark.parametrize("points,match", [
    ([[0., 0.], [0., 0.], [0., 0.]], "degenerate"),
    ([[1., 0.], [0., 0.], [.5, 0.]], "degenerate"),
    ([[1., 0.], [0., 0.]], "N >= 3"),
    ([[1., 0.], [np.nan, 0.], [0., 0.]], "non-finite"),
    ([[1., 0.], [.5, .1], [.5, .1], [0., 0.], [.5, -.1]], "duplicate"),
    ([[1., 0.], [.5, .1], [0., 0.], [.5, .1]], "duplicate"),
    ([[1., 0.], [.5, .1], [0., 0.], [1. + 1e-12, 0.]], "duplicate"),
])
def test_degenerate_and_duplicate_inputs(points, match):
    with pytest.raises(AirfoilFormatError, match=match):
        Airfoil.from_points(points)


def test_arrays_are_copied_and_cannot_be_unlocked(diamond_points):
    original = diamond_points.copy()
    foil = Airfoil.from_points(diamond_points)
    diamond_points[1] = [8., 9.]
    np.testing.assert_array_equal(foil.original_points, original)
    for array in (foil.original_points, foil.normalized_points, foil.leading_edge, foil.trailing_edge):
        with pytest.raises(ValueError):
            array.flat[0] = 9.
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_self_intersecting_input_not_repaired():
    points = [[1., 0.], [.3, .1], [0., 0.], [.7, .1], [.3, -.1], [1., 0.]]
    foil = Airfoil.from_points(points)
    assert foil.area is None
    np.testing.assert_array_equal(foil.original_points, points)


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_dat(tmp_path / "missing.dat")


def test_explicit_coordinate_header_disambiguates_translated_integer_pair(tmp_path, diamond_points):
    path = tmp_path / "translated.dat"
    np.savetxt(path, diamond_points + [2., 3.], fmt="%.17g",
               header="AeroOpt_COORDINATES_v0.1 translated", comments="")
    foil = Airfoil.from_dat(path)
    np.testing.assert_array_equal(foil.original_points, diamond_points + [2., 3.])
    np.testing.assert_allclose(foil.normalized_points, diamond_points, atol=1e-14, rtol=0)


def test_reversed_traversal_does_not_reflect(diamond_points):
    foil = Airfoil.from_points(diamond_points[::-1])
    np.testing.assert_array_equal(foil.normalized_points, diamond_points[::-1])


def test_invalid_dat_duplicate_rows_reported(tmp_path):
    path = tmp_path / "duplicate.dat"
    path.write_text("Demo\n1 0\n.5 .1\n0 0\n.5 .1\n1 0\n", encoding="utf-8")
    with pytest.raises(AirfoilFormatError, match="duplicate"):
        Airfoil.from_dat(path)

