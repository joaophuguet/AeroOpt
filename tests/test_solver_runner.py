"""Staged failures, audit integrity, result reimport, retries and path boundaries."""

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from aeroopt.solvers import (AerodynamicCoefficients, AttemptContext, ConvergenceStatus,
                            MockScenario, MockSolver, OperatingConditions, ParsedOutput,
                            RunStatus, SolverResult, SolverRunner)
from aeroopt.solvers.serialization import FrozenJSON, read_json, write_json
from aeroopt.solvers._paths import native_path, safe_path

CONDITIONS = OperatingConditions(5., 2e5, .05)
EXPECTED = {MockScenario.SUCCESS: RunStatus.SUCCESS, MockScenario.INPUT_REJECTED: RunStatus.INVALID_INPUT,
            MockScenario.PREPARATION_FAILED: RunStatus.PREPARATION_FAILED,
            MockScenario.EXECUTION_FAILED: RunStatus.EXECUTION_FAILED, MockScenario.TIMEOUT: RunStatus.TIMEOUT,
            MockScenario.PARSING_FAILED: RunStatus.PARSING_FAILED, MockScenario.NOT_CONVERGED: RunStatus.NOT_CONVERGED}


@pytest.mark.parametrize("scenario", list(MockScenario))
def test_all_mock_states_are_persisted_and_nonphysical(separated, tmp_path, scenario):
    runner = SolverRunner(tmp_path, "cases")
    before = separated.flap.original_points.copy()
    result = runner.run(MockSolver(scenario), separated, CONDITIONS, min_clearance=.002, timeout_s=.25)
    assert result.status == EXPECTED[scenario], result.diagnostics
    assert result.coefficients.empty
    assert not result.physical_coefficients_available
    assert result.solver.is_simulated
    path = runner.execution_root / result.attempt_directory
    restored = SolverResult.load(path / "result.json")
    assert restored.to_dict() == result.to_dict()
    assert result.duration_s >= 0
    for artifact in result.artifacts:
        file = safe_path(path, artifact.path)
        assert file.is_relative_to(path)
        assert hashlib.sha256(file.read_bytes()).hexdigest() == artifact.sha256
        assert file.stat().st_size == artifact.size_bytes
    np.testing.assert_array_equal(separated.flap.original_points, before)
    if scenario == MockScenario.PREPARATION_FAILED:
        assert (path / "partial_prepare.txt").exists()
    if scenario == MockScenario.PARSING_FAILED:
        assert "deliberately invalid" in (path / "mock_output.json").read_text()
    if scenario == MockScenario.SUCCESS:
        assert result.convergence == ConvergenceStatus.NOT_APPLICABLE
    if scenario == MockScenario.NOT_CONVERGED:
        assert result.convergence == ConvergenceStatus.NOT_CONVERGED


def test_invalid_geometry_and_explicit_clearance_block_preparation(separated, tmp_path):
    runner = SolverRunner(tmp_path)
    for geometry, clearance in [(separated.configure(flap_dx=-.7), .0), (separated, .3)]:
        result = runner.run(MockSolver(), geometry, CONDITIONS, min_clearance=clearance)
        assert result.status == RunStatus.INVALID_INPUT
        assert not result.geometry_validation.to_dict()["valid"]
        assert result.case_id is not None
        path = runner.execution_root / result.attempt_directory
        assert not (path / "mock_request.json").exists()
        assert result.execution_duration_s is None


def test_attempts_and_execution_directories_do_not_change_analysis_identity(separated, tmp_path):
    runner = SolverRunner(tmp_path, "first")
    first = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=.002, timeout_s=1.)
    snapshot = {artifact.path: (runner.execution_root / first.attempt_directory / artifact.path).read_bytes()
                for artifact in first.artifacts}
    second = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=.002, timeout_s=2.)
    elsewhere = SolverRunner(tmp_path, "different")
    third = elsewhere.run(MockSolver(), separated, CONDITIONS, min_clearance=.002)
    assert first.case_id == second.case_id == third.case_id
    assert first.attempt_id != second.attempt_id
    assert first.attempt_directory != second.attempt_directory
    for name, raw in snapshot.items():
        assert (runner.execution_root / first.attempt_directory / name).read_bytes() == raw
    assert read_json(runner.execution_root / first.attempt_directory / "input.json") == read_json(
        elsewhere.execution_root / third.attempt_directory / "input.json")


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf"), True, "1"])
def test_invalid_timeout_has_finite_structured_result(separated, tmp_path, timeout):
    runner = SolverRunner(tmp_path)
    result = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=.002, timeout_s=timeout)
    assert result.status == RunStatus.INVALID_INPUT
    assert result.case_id is None
    SolverResult.load(runner.execution_root / result.attempt_directory / "result.json")


@pytest.mark.parametrize("clearance", [-1, float("nan"), float("inf"), True])
def test_invalid_clearance_is_not_prepared(separated, tmp_path, clearance):
    runner = SolverRunner(tmp_path)
    result = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=clearance)
    assert result.status == RunStatus.INVALID_INPUT
    assert result.case_id is None
    assert result.coefficients.empty


def test_unrepresentable_transform_retains_diagnostic(separated, tmp_path):
    from aeroopt.geometry import FlapParameters, MultiElementAirfoil

    geometry = MultiElementAirfoil(separated.main, separated.flap, FlapParameters(1e308, (1e308, 0)))
    runner = SolverRunner(tmp_path)
    result = runner.run(MockSolver(), geometry, CONDITIONS, min_clearance=0)
    assert result.status == RunStatus.INVALID_INPUT
    assert result.case_id is None
    assert "unrepresentable" in result.diagnostics[0].message


class CoefficientLeak(MockSolver):
    def parse(self, case, context, evidence):
        # Artificial value ONLY for a negative contract test; publication must fail.
        return ParsedOutput(ConvergenceStatus.NOT_APPLICABLE, AerodynamicCoefficients(cl=1.))


class NonfiniteParser(MockSolver):
    def parse(self, case, context, evidence):
        return ParsedOutput(ConvergenceStatus.NOT_APPLICABLE, AerodynamicCoefficients(cd=float("nan")))


@pytest.mark.parametrize("solver", [CoefficientLeak(), NonfiniteParser()])
def test_invalid_parsed_values_cannot_be_published(separated, tmp_path, solver):
    result = SolverRunner(tmp_path).run(solver, separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.PARSING_FAILED
    assert result.coefficients.empty
    assert not result.physical_coefficients_available


def test_result_rejects_inconsistent_or_forged_records(separated, tmp_path):
    result = SolverRunner(tmp_path).run(MockSolver(), separated, CONDITIONS, min_clearance=.002)
    with pytest.raises(ValueError): replace(result, coefficients=AerodynamicCoefficients(cd=0.))
    with pytest.raises(ValueError): replace(result, case_id="a" * 64)
    with pytest.raises(ValueError): replace(result, duration_s=float("nan"))
    with pytest.raises(ValueError): replace(result, provenance=FrozenJSON({"invalid": float("inf")}))
    with pytest.raises(ValueError): replace(result, convergence=ConvergenceStatus.NOT_CONVERGED)
    with pytest.raises(ValueError): replace(result, started_at_utc="2026-01-01T00:00:00")
    data = result.to_dict()
    data["physical_coefficients_available"] = True
    path = native_path(tmp_path / "tampered.json")
    write_json(path, data)
    with pytest.raises(ValueError, match="physical-result"):
        SolverResult.load(path)


@pytest.mark.parametrize("path", ["../escape", "a/../../b", "/outside", "C:/outside", "C:relative",
                                "//server/share/x", "x:stream", "CON", "logs/NUL.txt", "a/./b", "x.", "x "])
def test_paths_cannot_escape_or_use_windows_devices(tmp_path, path):
    with pytest.raises(ValueError): SolverRunner(tmp_path, path)
    with pytest.raises(ValueError): AttemptContext(tmp_path, None).path(path)


def test_long_windows_paths_spaces_and_unicode(separated, tmp_path):
    directory = "/".join(["long_" + "a" * 45] * 4 + ["espaços ç"])
    runner = SolverRunner(tmp_path, directory)
    result = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.SUCCESS, result.diagnostics
    assert len(str(runner.execution_root / result.attempt_directory / "result.json")) > 260
    assert SolverResult.load(runner.execution_root / result.attempt_directory / "result.json").case_id == result.case_id


def test_existing_result_and_adapter_files_not_overwritten(separated, tmp_path):
    runner = SolverRunner(tmp_path)
    result = runner.run(MockSolver(), separated, CONDITIONS, min_clearance=.002)
    path = runner.execution_root / result.attempt_directory
    before = (path / "result.json").read_bytes()
    with pytest.raises(FileExistsError): result.save(path / "result.json")
    with pytest.raises(FileExistsError): AttemptContext(path, None).write_text("adapter.log", "overwrite")
    assert (path / "result.json").read_bytes() == before


def test_junction_guard_without_requiring_symlink_privileges(tmp_path, monkeypatch):
    link = native_path(tmp_path / "junction")
    original = Path.is_junction
    monkeypatch.setattr(Path, "is_junction", lambda path: path == link or original(path))
    with pytest.raises(ValueError, match="junction"):
        AttemptContext(tmp_path, None).path("junction/outside.txt")


def test_identity_drift_rejected(separated, tmp_path):
    class Drift(MockSolver):
        def prepare(self, case, context):
            super().prepare(case, context)
            object.__setattr__(self, "scenario", MockScenario.NOT_CONVERGED)
    result = SolverRunner(tmp_path).run(Drift(), separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.PREPARATION_FAILED
    assert "identity changed" in result.diagnostics[-1].message


def test_python_library_system_exit_is_a_structured_failure(separated, tmp_path):
    class LibraryExit(MockSolver):
        def execute(self, case, context):
            context.write_text("partial_library.txt", "evidence before SystemExit")
            raise SystemExit(2)
    result = SolverRunner(tmp_path).run(LibraryExit(), separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.EXECUTION_FAILED
    assert result.diagnostics[-1].code == "SystemExit"
    assert any(a.path == "partial_library.txt" for a in result.artifacts)


def test_keyboard_interrupt_is_persisted_then_respected(separated, tmp_path):
    class Interrupted(MockSolver):
        def execute(self, case, context):
            raise KeyboardInterrupt()
    runner = SolverRunner(tmp_path)
    with pytest.raises(KeyboardInterrupt):
        runner.run(Interrupted(), separated, CONDITIONS, min_clearance=.002)
    paths = list(runner.execution_root.glob("attempt-*/result.json"))
    assert len(paths) == 1
    assert SolverResult.load(paths[0]).status == RunStatus.EXECUTION_FAILED


def test_adapter_cannot_overwrite_framework_input(separated, tmp_path):
    class Collision(MockSolver):
        def prepare(self, case, context):
            context.write_json("input.json", {"changed": True})
    result = SolverRunner(tmp_path).run(Collision(), separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.PREPARATION_FAILED
    assert result.diagnostics[-1].code == "FileExistsError"


def test_library_evidence_may_have_no_process_return_code(separated, tmp_path):
    class Library(MockSolver):
        def execute(self, case, context):
            return replace(super().execute(case, context), return_code=None)
    result = SolverRunner(tmp_path).run(Library(), separated, CONDITIONS, min_clearance=.002)
    assert result.status == RunStatus.SUCCESS


def test_demo_cli_saves_all_states_and_refuses_overwrite(tmp_path):
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    output = native_path(tmp_path / "demo")
    relative = output.relative_to(native_path(root)).as_posix()
    command = [sys.executable, str(root / "examples/demo_solver_interface.py"), "--output", relative]
    run = subprocess.run(command, cwd=root, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    data = read_json(output / "summary.json")
    assert len(data["runs"]) == 9
    assert data["runs"][0]["case_id"] == data["runs"][7]["case_id"]
    assert data["runs"][0]["attempt_id"] != data["runs"][7]["attempt_id"]
    assert all(row["coefficients"]["cl"] is None for row in data["runs"])
    assert set(row["status"] for row in data["runs"]) == {state.value for state in RunStatus}
    snapshot = (output / "summary.json").read_bytes()
    repeated = subprocess.run(command, cwd=root, capture_output=True, text=True)
    assert repeated.returncode != 0
    assert (output / "summary.json").read_bytes() == snapshot
