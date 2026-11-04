"""Typed immutable inputs, finite JSON, exact scientific identity and result guards."""

import json
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil
from aeroopt.solvers import (AerodynamicCoefficients, ConditionParameter, FrozenJSON,
                            MockSolver, MomentReference, OperatingConditions, SolverCase)
from aeroopt.solvers import ConvergenceStatus, ExecutionEvidence, ParsedOutput
from aeroopt.solvers.serialization import canonical_json, read_json, scientific_hash, write_json


@pytest.mark.parametrize("field", ["alpha_deg", "reynolds", "mach"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), "1", True, None])
def test_operating_conditions_reject_bad_scalars(field, value):
    arguments = {"alpha_deg": 5., "reynolds": 200000., "mach": .05, field: value}
    with pytest.raises((TypeError, ValueError)):
        OperatingConditions(**arguments)


@pytest.mark.parametrize("reynolds,mach", [(0, 0), (-1, 0), (1, -1)])
def test_operating_domains(reynolds, mach):
    with pytest.raises(ValueError):
        OperatingConditions(0, reynolds, mach)


def test_no_arbitrary_solver_limits_and_no_physical_defaults():
    for alpha, re, mach in [(-1000., 1e-200, 0.), (10000., 1e100, 30.)]:
        item = OperatingConditions(alpha, re, mach)
        assert item.to_dict()["reynolds_reference"] == "main_chord"
        assert item.additional == ()
    with pytest.raises(TypeError):
        OperatingConditions()


def test_conditions_and_metadata_deep_immutability():
    params = [ConditionParameter("transition.main.upper", .3, "x/c_main")]
    conditions = OperatingConditions(5, 2e5, .05, params)
    params.clear()
    assert len(conditions.additional) == 1
    with pytest.raises(FrozenInstanceError):
        conditions.alpha_deg = 8
    data = {"nested": {"values": [1., 2.]}}
    frozen = FrozenJSON(data)
    data["nested"]["values"].append(3.)
    detached = frozen.to_dict()
    detached["nested"]["values"].append(4.)
    assert frozen.to_dict() == {"nested": {"values": [1., 2.]}}


def test_additional_fields_are_explicit_typed_and_unique():
    with pytest.raises(ValueError):
        ConditionParameter("n_crit", 9.)
    with pytest.raises(ValueError):
        ConditionParameter("n_crit", float("nan"), "dimensionless")
    with pytest.raises(TypeError):
        ConditionParameter("array", [1.], "dimensionless")
    parameter = ConditionParameter("n_crit", 9., "dimensionless")
    with pytest.raises(ValueError):
        OperatingConditions(0, 1e5, 0, [parameter, parameter])
    assert OperatingConditions.from_dict(OperatingConditions(3., 1e5, .1, [parameter]).to_dict()).additional == (parameter,)


def test_case_identity_order_path_independence_and_coordinate_preservation(separated, tmp_path):
    parameters = [ConditionParameter("z", "explicit"), ConditionParameter("a", .3, "dimensionless")]
    first = SolverCase(separated, OperatingConditions(5, 2e5, .05, parameters), MockSolver().identity, .002)
    other = MultiElementAirfoil(replace(separated.main, source=tmp_path / "elsewhere.dat"),
                               replace(separated.flap, name="different display name"), separated.parameters)
    second = SolverCase(other, OperatingConditions(5., 200000., .05, list(reversed(parameters))),
                        MockSolver().identity, .002)
    assert first.case_id == second.case_id
    for arr in [first.main_points, first.flap_points]:
        with pytest.raises(ValueError):
            arr.setflags(write=True)
    np.testing.assert_array_equal(first.main_points, separated.main_points)
    np.testing.assert_array_equal(first.flap_points, separated.flap_points)
    assert first.geometry is separated


@pytest.mark.parametrize("change", ["alpha", "reynolds", "mach", "additional", "coordinates", "version",
                                    "identifier", "settings", "schema", "clearance", "simulation"])
def test_hash_sensitive_to_relevant_inputs(separated, change):
    conditions = OperatingConditions(5., 2e5, .05)
    identity = MockSolver().identity
    original = SolverCase(separated, conditions, identity, .002)
    geometry, clearance, schema = separated, .002, original.schema_version
    if change == "alpha": conditions = replace(conditions, alpha_deg=5.01)
    if change == "reynolds": conditions = replace(conditions, reynolds=200001.)
    if change == "mach": conditions = replace(conditions, mach=.051)
    if change == "additional": conditions = replace(conditions, additional=(ConditionParameter("n_crit", 8., "dimensionless"),))
    if change == "coordinates": geometry = separated.configure(flap_dx=np.nextafter(.0, 1.) + .01)
    if change == "version": identity = replace(identity, version="next")
    if change == "identifier": identity = replace(identity, identifier="other")
    if change == "settings": identity = replace(identity, settings=FrozenJSON({"new": True}))
    if change == "schema": schema = "aeroopt.solver-input.future-test-only"
    if change == "clearance": clearance = .003
    if change == "simulation": identity = replace(identity, is_simulated=False)
    modified = SolverCase(geometry, conditions, identity, clearance, schema)
    assert modified.case_id != original.case_id


def test_hash_canonical_numbers_and_keys():
    assert scientific_hash({"b": 2., "a": -0.}) == scientific_hash({"a": 0., "b": 2.})
    assert scientific_hash(1.) != scientific_hash(np.nextafter(1., 2.).item())
    assert scientific_hash(1.) != scientific_hash(1)
    assert scientific_hash(1.) != scientific_hash({"float64_hex": "0x1.0000000000000p+0"})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_rejected_throughout_contract(value, tmp_path):
    with pytest.raises(ValueError): FrozenJSON({"nested": [value]})
    with pytest.raises(ValueError): canonical_json({"value": value})
    with pytest.raises(ValueError): scientific_hash({"value": value})
    with pytest.raises(ValueError): AerodynamicCoefficients(cl=value)
    with pytest.raises(ValueError): MomentReference(value, 0, "clockwise", "explicit test only")
    with pytest.raises(ValueError): write_json(tmp_path / "bad.json", {"value": value})
    assert not (tmp_path / "bad.json").exists()


@pytest.mark.parametrize("text", ['{"v":NaN}', '{"v":Infinity}', '{"v":1e999}', '{"v":1,"v":2}', '[]'])
def test_json_import_rejects_invalid_numbers_and_ambiguity(tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        read_json(path)


def test_missing_values_and_explicit_moment_reference():
    assert AerodynamicCoefficients().to_dict() == {"cl": None, "cd": None, "cm": None, "moment_reference": None}
    with pytest.raises(ValueError): AerodynamicCoefficients(cm=0.)
    # Artificial scalar used ONLY to verify serialization, never solver output.
    reference = MomentReference(.25, 0., "counterclockwise", "UNIT TEST ONLY: q_inf*c_main^2 per span")
    coefficient = AerodynamicCoefficients(cm=0., moment_reference=reference)
    assert coefficient.cm == 0.
    assert AerodynamicCoefficients.from_dict(coefficient.to_dict()) == coefficient


def test_json_atomic_publication_no_overwrite(tmp_path):
    path = tmp_path / "object.json"
    write_json(path, {"x": 1.})
    before = path.read_bytes()
    with pytest.raises(FileExistsError): write_json(path, {"x": 2.})
    assert path.read_bytes() == before
    assert not list(tmp_path.glob("*.pending"))


@pytest.mark.parametrize("field", ["capabilities", "messages", "command"])
def test_sequences_do_not_silently_split_strings(field):
    with pytest.raises(TypeError):
        if field == "capabilities": replace(MockSolver().identity, capabilities="not a tuple")
        if field == "messages": ParsedOutput(ConvergenceStatus.NOT_APPLICABLE, messages="not a tuple")
        if field == "command": ExecutionEvidence(0, 0., command="not argv")


def test_angle_convention_is_explicit_and_reimport_checked():
    data = OperatingConditions(5., 2e5, .05).to_dict()
    assert data["alpha_convention"] == "freestream_vector_(cos_alpha,sin_alpha)_in_global_xy"
    data["alpha_convention"] = "different sign convention"
    with pytest.raises(ValueError, match="angle-of-attack"):
        OperatingConditions.from_dict(data)
