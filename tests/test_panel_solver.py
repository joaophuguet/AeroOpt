"""Panel adapter lifecycle, scientific identity and refusal to publish failures."""

import hashlib
import json
import sys
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import subprocess

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points
from aeroopt.solvers import (ConditionParameter, ConvergenceStatus, MomentReference, OperatingConditions,
                            PanelSolver, RunStatus, SolverCase, SolverResult, SolverRunner)
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers._paths import safe_path
from aeroopt.solvers.serialization import read_json

REFERENCE = MomentReference(.25, 0, "clockwise", "q_inf*c_main^2")
CONDITIONS = OperatingConditions(6, 2e5, 0)
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def geometry():
    main = Airfoil.from_points(naca4_points("0012", 41, closed_te=True), leading_edge_index=20)
    flap = Airfoil.from_points(naca4_points("4412", 41, closed_te=True), leading_edge_index=20)
    return MultiElementAirfoil(main, flap, FlapParameters(.3, (1.02, -.04), 25))


def test_requires_explicit_moment_definition():
    with pytest.raises(TypeError):
        PanelSolver()
    with pytest.raises(TypeError):
        PanelSolver(None)
    with pytest.raises(ValueError, match="normalization"):
        PanelSolver(replace(REFERENCE, normalization="unknown"))
    solver = PanelSolver(REFERENCE)
    with pytest.raises(FrozenInstanceError):
        solver.moment_reference = replace(REFERENCE, x_over_c=.5)
    assert not solver.identity.is_simulated
    assert solver.identity.convergence_required
    assert "cd" not in solver.identity.capabilities


def test_lifecycle_cp_circulation_artifacts_and_result_roundtrip(geometry, tmp_path):
    runner = SolverRunner(tmp_path, "panel")
    original = (geometry.main_points.copy(), geometry.flap_points.copy())
    result = runner.run(PanelSolver(REFERENCE), geometry, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.SUCCESS, result.diagnostics
    assert result.convergence == ConvergenceStatus.CONVERGED
    assert result.physical_coefficients_available
    assert result.coefficients.cl > 0
    assert result.coefficients.cm < 0  # clockwise-positive, flap load aft
    assert result.coefficients.cd is None
    path = runner.execution_root / result.attempt_directory
    restored = SolverResult.load(path / "result.json")
    assert restored.to_dict() == result.to_dict()
    output = read_json(path / "panel_output.json")
    assert output["case_id"] == result.case_id
    assert output["moment_reference"] == REFERENCE.to_dict()
    assert len(output["elements"]) == 2
    assert output["diagnostics"]["unknown_count"] == 82
    assert output["diagnostics"]["impermeability_max"] < 1e-12
    for element in output["elements"]:
        assert len(element["cp"]) == element["panel_count"]
        assert element["circulation"] == pytest.approx(element["vortex_sheet_strength"] * sum(element["lengths"]))
    for artifact in result.artifacts:
        file = safe_path(path, artifact.path)
        assert hashlib.sha256(file.read_bytes()).hexdigest() == artifact.sha256
    np.testing.assert_array_equal(geometry.main_points, original[0])
    np.testing.assert_array_equal(geometry.flap_points, original[1])
    assert (path / "panel.log").is_file()


def test_repeated_attempt_identity_and_changed_settings(geometry, tmp_path):
    runner = SolverRunner(tmp_path, "panel")
    solver = PanelSolver(REFERENCE)
    first = runner.run(solver, geometry, CONDITIONS, min_clearance=.002)
    second = runner.run(solver, geometry, CONDITIONS, min_clearance=.002)
    assert first.case_id == second.case_id
    assert first.attempt_id != second.attempt_id
    changed = SolverCase(geometry, CONDITIONS, PanelSolver(replace(REFERENCE, x_over_c=.5)).identity, .002)
    assert changed.case_id != first.case_id
    changed_re = SolverCase(geometry, replace(CONDITIONS, reynolds=3e5), solver.identity, .002)
    assert changed_re.case_id != first.case_id  # preserved despite Re-independent physics


@pytest.mark.parametrize("conditions", [OperatingConditions(6, 2e5, .01),
    OperatingConditions(6, 2e5, 0, (ConditionParameter("transition", .5, "dimensionless"),))])
def test_unsupported_physics_rejected_before_prepare(geometry, conditions, tmp_path):
    runner = SolverRunner(tmp_path)
    result = runner.run(PanelSolver(REFERENCE), geometry, conditions, min_clearance=.002)
    assert result.status == RunStatus.INVALID_INPUT
    assert result.coefficients.empty
    assert not (runner.execution_root / result.attempt_directory / "panel_request.json").exists()


def test_tiny_clearance_and_invalid_geometry_are_not_repaired(tmp_path):
    p = naca4_points("0012", 41, closed_te=True)
    foil = Airfoil.from_points(p, leading_edge_index=20)
    runner = SolverRunner(tmp_path)
    for x, code in ((1 + 5e-7, "panel numerical clearance"), (.8, "geometry rejected")):
        geometry = MultiElementAirfoil(foil, foil, FlapParameters(.3, (x, 0)))
        before = geometry.flap_points.copy()
        result = runner.run(PanelSolver(REFERENCE), geometry, CONDITIONS, min_clearance=0)
        assert result.status == RunStatus.INVALID_INPUT
        assert code in result.diagnostics[0].message
        assert result.coefficients.empty
        np.testing.assert_array_equal(geometry.flap_points, before)


@pytest.mark.parametrize("failure", ["condition", "singular", "residual"])
def test_numerical_failure_is_not_converged_without_coefficients(geometry, tmp_path, monkeypatch, failure):
    if failure == "condition":
        monkeypatch.setattr(np.linalg, "cond", lambda a: 1e12)
    elif failure == "singular":
        def fail(a, b):
            raise np.linalg.LinAlgError("injected singular solve")
        monkeypatch.setattr(np.linalg, "solve", fail)
    else:
        monkeypatch.setattr(np.linalg, "solve", lambda a, b: np.zeros_like(b))
    runner = SolverRunner(tmp_path)
    result = runner.run(PanelSolver(REFERENCE), geometry, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.NOT_CONVERGED
    assert result.convergence == ConvergenceStatus.NOT_CONVERGED
    assert result.coefficients.empty and not result.physical_coefficients_available
    output = read_json(runner.execution_root / result.attempt_directory / "panel_output.json")
    assert not output["accepted"] and "cl" not in output


def test_execution_exception_preserves_evidence(geometry, tmp_path, monkeypatch):
    def fail(*args):
        raise RuntimeError("injected execution failure")
    monkeypatch.setattr(pm, "assemble", fail)
    runner = SolverRunner(tmp_path)
    result = runner.run(PanelSolver(REFERENCE), geometry, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.EXECUTION_FAILED
    assert result.coefficients.empty
    directory = runner.execution_root / result.attempt_directory
    assert (directory / "panel_request.json").is_file()
    assert (directory / "traceback.txt").is_file()


@pytest.mark.parametrize("defect", ["json", "case_id", "cp", "residual", "drag", "moment", "circulation", "mesh", "limitations"])
def test_corrupt_output_is_parsing_failure(geometry, tmp_path, monkeypatch, defect):
    original = PanelSolver.execute
    def corrupt(self, case, context):
        evidence = original(self, case, context)
        path = context.path("panel_output.json")
        output = read_json(path)
        if defect == "json":
            path.write_text("{bad", encoding="utf-8")
            return evidence
        if defect == "case_id":
            output["case_id"] = "0" * 64
        elif defect == "cp":
            output["elements"][0]["cp"][0] += 1
        elif defect == "residual":
            output["diagnostics"]["kutta_max"] = .1
        elif defect == "drag":
            output["cd"] = 0
        elif defect == "moment":
            output["cm"] += 1
        elif defect == "circulation":
            output["elements"][0]["circulation"] += 1
        elif defect == "mesh":
            output["elements"][0]["centers"][0][0] += 1
        else:
            output.pop("limitations")
        path.write_text(json.dumps(output), encoding="utf-8")
        return evidence
    monkeypatch.setattr(PanelSolver, "execute", corrupt)
    result = SolverRunner(tmp_path).run(PanelSolver(REFERENCE), geometry, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.PARSING_FAILED
    assert result.coefficients.empty
    assert any(a.path == "panel_output.json" for a in result.artifacts)


def test_cooperative_timeout_cannot_publish_coefficients(geometry, tmp_path):
    result = SolverRunner(tmp_path).run(PanelSolver(REFERENCE), geometry, CONDITIONS,
                                      min_clearance=.002, timeout_s=1e-9)
    assert result.status == RunStatus.TIMEOUT
    assert result.coefficients.empty


def test_demo_cli_and_refusal_to_overwrite(tmp_path):
    # Output is explicitly inside the repository, as required by the example.
    output = tmp_path / "demo"
    command = [sys.executable, str(ROOT / "examples/demo_panel_solver.py"), "--output", str(output)]
    first = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
    assert first.returncode == 0, first.stdout + first.stderr
    summary = read_json(output / "summary.json")
    assert summary["status"] == "success" and summary["coefficients"]["cd"] is None
    second = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
    assert second.returncode != 0 and "empty" in second.stderr


def test_verification_cli_persists_finite_analytical_and_runner_evidence(tmp_path):
    output = tmp_path / "verification"
    completed = subprocess.run([sys.executable, str(ROOT / "examples/verify_panel_solver.py"),
                               "--output", str(output)], cwd=ROOT, capture_output=True, text=True, timeout=90)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    report = read_json(output / "verification.json")
    assert len(report["joukowski"]) == 5
    assert report["joukowski"][-1]["cm_error"] < .003
    assert len(report["multielement_refinement"]) == 4
    assert len(report["finite_te_refinement"]) == 4
    assert (output / "refinement.png").is_file()
    assert "not viscous physical validation" in report["scope"]
