"""Real local subprocesses: timeout, process trees, logs and protected inputs."""

import sys
import time
from dataclasses import replace

import pytest

from aeroopt.solvers import (ConvergenceStatus, MockSolver, OperatingConditions, ParsedOutput,
                            RunStatus, SolverResult, SolverRunner)
from aeroopt.solvers.serialization import FrozenJSON, read_json
from solver_helpers import ProcessFixtureSolver, process_running

CONDITIONS = OperatingConditions(5., 2e5, .05)


def test_real_subprocess_logs_and_nonzero_status(separated, tmp_path):
    runner = SolverRunner(tmp_path)
    success = runner.run(ProcessFixtureSolver(), separated, CONDITIONS, min_clearance=.002, timeout_s=3.)
    assert success.status == RunStatus.SUCCESS, success.diagnostics
    path = runner.execution_root / success.attempt_directory
    assert "STDOUT" in (path / "process.stdout.log").read_text()
    assert "STDERR" in (path / "process.stderr.log").read_text()
    assert read_json(path / "execution.json")["return_code"] == 0
    failed = runner.run(ProcessFixtureSolver("exit7"), separated, CONDITIONS, min_clearance=.002, timeout_s=3.)
    assert failed.status == RunStatus.EXECUTION_FAILED
    assert "7" in failed.diagnostics[0].message


@pytest.mark.parametrize("mode,status", [("blocking_tree", RunStatus.TIMEOUT), ("orphan_on_success", RunStatus.SUCCESS)])
def test_blocking_tree_or_surviving_descendant_is_terminated(separated, tmp_path, mode, status):
    runner = SolverRunner(tmp_path)
    started = time.monotonic()
    result = runner.run(ProcessFixtureSolver(mode), separated, CONDITIONS, min_clearance=.002, timeout_s=1.5)
    assert result.status == status, result.diagnostics
    assert time.monotonic() - started < 10
    path = runner.execution_root / result.attempt_directory
    # A handshake ensures both actual parent and descendant existed BEFORE timeout.
    assert "CHILD_READY" in (path / "process.stdout.log").read_text()
    parent = int((path / "parent.pid").read_text())
    child = int((path / "child.pid").read_text())
    assert not process_running(parent)
    assert not process_running(child)
    evidence = read_json(path / "execution.json")
    assert evidence["timed_out"] == (status == RunStatus.TIMEOUT)
    assert "terminated" in evidence["tree_cleanup"]
    assert SolverResult.load(path / "result.json").status == status


def test_process_completion_does_not_imply_convergence(separated, tmp_path):
    class UnknownConvergence(ProcessFixtureSolver):
        @property
        def identity(self):
            # Contract fixture only: this declaration is not a physical solver.
            return replace(super().identity, convergence_required=True)
        def parse(self, case, context, evidence):
            return ParsedOutput(ConvergenceStatus.UNKNOWN)
    result = SolverRunner(tmp_path).run(UnknownConvergence(), separated, CONDITIONS,
                                        min_clearance=.002, timeout_s=3.)
    assert result.status == RunStatus.NOT_CONVERGED
    assert result.convergence == ConvergenceStatus.UNKNOWN


def test_subprocess_input_mutation_cannot_produce_success(separated, tmp_path):
    result = SolverRunner(tmp_path).run(ProcessFixtureSolver("mutate_input"), separated, CONDITIONS,
                                       min_clearance=.002, timeout_s=3.)
    assert result.status == RunStatus.PARSING_FAILED
    assert result.coefficients.empty
    assert not result.physical_coefficients_available


def test_failed_job_assignment_never_releases_target(separated, tmp_path, monkeypatch):
    import os
    if os.name != "nt":
        pytest.skip("Windows Job Object fail-closed test")
    from aeroopt.solvers._windows_job import WindowsJob

    def denied(self, pid):
        raise PermissionError("simulated job assignment denial")
    monkeypatch.setattr(WindowsJob, "assign", denied)
    result = SolverRunner(tmp_path).run(ProcessFixtureSolver("blocking_tree"), separated, CONDITIONS,
                                       min_clearance=.002, timeout_s=1.)
    assert result.status == RunStatus.EXECUTION_FAILED
    assert "assignment denial" in result.diagnostics[0].message
    path = SolverRunner(tmp_path).execution_root / result.attempt_directory
    assert not (path / "parent.pid").exists()
    assert not (path / "child.pid").exists()


def test_argv_is_not_a_shell_and_output_paths_are_guarded(separated, tmp_path):
    class Arguments(MockSolver):
        def execute(self, case, context):
            return context.run_process([sys.executable, "-c", "import sys; print(sys.argv[1])",
                                       "literal & echo unsafe > escaped.txt"])
    runner = SolverRunner(tmp_path)
    result = runner.run(Arguments(), separated, CONDITIONS, min_clearance=.002, timeout_s=3.)
    assert result.status == RunStatus.PARSING_FAILED  # No mock output file: execute is a fixture only.
    path = runner.execution_root / result.attempt_directory
    assert "literal & echo unsafe > escaped.txt" in (path / "process.stdout.log").read_text()
    assert not (path / "escaped.txt").exists()
