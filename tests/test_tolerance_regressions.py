"""Regressions identified with realistic NACA shapes, not relaxed tolerances."""

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points, validate_elements


@pytest.mark.parametrize("ratio", [1e-5, 1e-6, 1e-7])
def test_small_resolved_flap_is_not_degenerate_due_to_its_scale(ratio):
    foil = Airfoil.from_points(naca4_points("0012", 81), leading_edge_index=40)
    geometry = MultiElementAirfoil(foil, foil, FlapParameters(ratio, (1.02, 0.)))
    report = geometry.validate(min_clearance=.01)
    assert report.valid, report.to_dict()
    assert report.min_distance == pytest.approx(.02, abs=1e-14)


def test_legitimate_nonuniform_sampling_does_not_change_validity():
    # Splitting a straight polygon segment changes only sampling, not its shape.
    main = Airfoil.from_points(naca4_points("0012", 81), leading_edge_index=40).normalized_points
    flap = main * 1e-5 + [1.02, 0.]
    dense = np.insert(flap, 1, flap[0] + 1e-4 * (flap[1] - flap[0]), axis=0)
    coarse_report = validate_elements(main, flap)
    dense_report = validate_elements(main, dense)
    assert coarse_report.valid, coarse_report.to_dict()
    assert dense_report.valid, dense_report.to_dict()
    assert coarse_report.min_distance == pytest.approx(dense_report.min_distance, abs=1e-14)
