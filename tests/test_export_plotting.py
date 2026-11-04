"""End-to-end global frame round trips, export guardrails, and figures."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from aeroopt.geometry import Airfoil, InvalidGeometryError, MultiElementAirfoil, read_dat, validate_elements
from aeroopt.geometry import FlapParameters
from aeroopt.geometry.plotting import compare_configurations, save_figure


def test_export_global_frame_roundtrip_and_reproducibility(separated, tmp_path):
    geometry = separated.configure(flap_deflection_deg=25., flap_dx=.02, flap_dy=-.03)
    target = geometry.export(tmp_path / "one", configuration_id="candidate_001", min_clearance=.002)
    metadata = json.loads((target / "metadata.json").read_text(encoding="utf-8"))
    _, main = read_dat(target / "main.dat")
    _, flap = read_dat(target / "flap.dat")
    np.testing.assert_array_equal(main, geometry.main_points)
    np.testing.assert_array_equal(flap, geometry.flap_points)
    np.testing.assert_array_equal(Airfoil.from_dat(target / "flap.dat").original_points, flap)
    report = validate_elements(main, flap, min_clearance=.002)
    assert metadata["validation"] == json.loads(json.dumps(report.to_dict()))
    assert metadata["parameters"]["deflection_deg"] == 25.
    assert metadata["units"]["coordinates"] == "x/c_main, y/c_main"
    assert metadata["configuration_id"] == "candidate_001"
    second = geometry.export(tmp_path / "two", configuration_id="candidate_001", min_clearance=.002)
    for filename in ("main.dat", "flap.dat", "metadata.json"):
        assert (target / filename).read_bytes() == (second / filename).read_bytes()


def test_invalid_geometry_never_exported(separated, tmp_path):
    geometry = separated.configure(flap_dx=-.7)
    target = tmp_path / "bad"
    with pytest.raises(InvalidGeometryError) as error:
        geometry.export(target)
    assert not error.value.report.valid
    assert not target.exists()


def test_export_clearance_must_be_passed_explicitly(separated, tmp_path):
    assert separated.validate().valid
    with pytest.raises(InvalidGeometryError):
        separated.export(tmp_path / "insufficient", min_clearance=.3)


def test_export_preserves_existing_files(separated, tmp_path):
    target = separated.export(tmp_path / "one")
    before = {file.name: file.read_bytes() for file in target.iterdir()}
    with pytest.raises(FileExistsError):
        separated.export(target)
    assert before == {file.name: file.read_bytes() for file in target.iterdir()}


def test_independent_point_counts_and_source_metadata(tmp_path):
    root = Path(__file__).resolve().parents[1]
    geometry = MultiElementAirfoil.from_dat(root / "examples/data/demo_main.dat",
        root / "examples/data/demo_flap.dat", flap_reference_position=(1.1, -.1))
    assert len(geometry.main_points) == 9
    assert len(geometry.flap_points) == 7
    target = geometry.export(tmp_path / "demo", min_clearance=.002)
    meta = json.loads((target / "metadata.json").read_text(encoding="utf-8"))
    assert meta["elements"]["main"]["source_file"] == str(geometry.main.source)
    assert meta["elements"]["flap"]["normalization"]["chord_original_units"] == 2.


def test_plot_labels_aspect_and_files(separated, tmp_path):
    ax = separated.plot()
    assert ax.get_aspect() == 1.
    assert ax.get_xlabel() and ax.get_ylabel()
    assert ax.get_legend() is not None
    assert len(ax.lines) == 2
    paths = save_figure(ax.figure, tmp_path / "geometry")
    assert paths[0].read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert "<svg" in paths[1].read_text(encoding="utf-8")
    plt.close(ax.figure)
    fig = compare_configurations([separated, separated.configure(flap_dx=-.7), separated],
                                  labels=["a", "b", "c"], min_clearance=.002)
    assert "REJEITADA" in fig.axes[1].get_title()
    assert not fig.axes[-1].get_visible()
    save_figure(fig, tmp_path / "comparison")
    plt.close(fig)


def test_invalid_comparison(separated):
    with pytest.raises(ValueError):
        compare_configurations([])
    with pytest.raises(ValueError):
        compare_configurations([separated], labels=[])


def test_global_integer_coordinates_export_reimport(diamond, tmp_path):
    # The first flap row is (3,3), which could otherwise resemble a count header.
    geometry = MultiElementAirfoil(diamond, diamond, FlapParameters(1., (2., 3.)))
    target = geometry.export(tmp_path / "translated")
    np.testing.assert_array_equal(Airfoil.from_dat(target / "flap.dat").original_points,
                                  geometry.flap_points)


def test_self_intersection_is_rejected_at_configuration_level(diamond, tmp_path):
    crossed = Airfoil.from_points([[1., 0.], [.3, .1], [0., 0.], [.7, .1], [.3, -.1], [1., 0.]])
    geometry = MultiElementAirfoil(diamond, crossed, FlapParameters(.3, (1.2, 0.)))
    report = geometry.validate()
    assert any(issue.code == "self_intersection" and issue.element == "flap" for issue in report.issues)
    with pytest.raises(InvalidGeometryError):
        geometry.export(tmp_path / "crossed")


@pytest.mark.parametrize("identifier", ["", "\ninvalid", "one\rtwo"])
def test_bad_configuration_id(separated, tmp_path, identifier):
    with pytest.raises(ValueError, match="configuration_id"):
        separated.export(tmp_path / "bad_id", configuration_id=identifier)
