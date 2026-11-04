"""Sequential staged runner: immutable analyses, isolated auditable attempts."""

import hashlib
import platform
import sys
import time
import traceback
from dataclasses import replace
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import numpy as np

from aeroopt import __version__
from aeroopt.geometry import MultiElementAirfoil, read_dat, validate_elements

from ._paths import native_path, reject_links, safe_path
from .base import AerodynamicSolver, SolverCase, SolverIdentity
from .conditions import OperatingConditions, real_number
from .execution import AttemptContext, ExecutionEvidence, utc_now
from .results import (AerodynamicCoefficients, Artifact, ConvergenceStatus, Diagnostic,
                      ParsedOutput, RunStatus, SolverResult)
from .serialization import FrozenJSON, read_json, scientific_hash


def _artifact(root: Path, path: Path) -> Artifact:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return Artifact(path.relative_to(root).as_posix(), digest.hexdigest(), path.stat().st_size)


def _inventory(root: Path) -> tuple[tuple[Artifact, ...], tuple[Diagnostic, ...]]:
    artifacts, problems = [], []

    def visit(directory: Path) -> None:
        try:
            entries = sorted(directory.iterdir())
        except OSError as error:
            problems.append(Diagnostic("persist", "unreadable_directory", str(error)))
            return
        for path in entries:
            try:
                reject_links(path)
                if path.is_dir():
                    visit(path)
                elif path.is_file() and not path.name.endswith(".pending") and path.name != "result.json":
                    artifacts.append(_artifact(root, path))
            except (OSError, ValueError) as error:
                problems.append(Diagnostic("persist", "unsafe_or_unreadable_artifact", str(error)))

    visit(root)
    return tuple(artifacts), tuple(problems)


def _check_inputs(case: SolverCase, context: AttemptContext) -> None:
    if read_json(context.path("input.json")) != case.to_dict():
        raise ValueError("adapter modified the immutable scientific input manifest")
    for name, expected in (("main", case.main_points), ("flap", case.flap_points)):
        _, actual = read_dat(context.path(f"geometry/{name}.dat"))
        if not np.array_equal(actual, expected):
            raise ValueError(f"adapter modified the original global {name} coordinates")


def _geometry_provenance(geometry: MultiElementAirfoil) -> dict:
    return {"parameters": geometry.parameters.to_dict(), "elements": {
        name: {"name": foil.name, "source_path_display_only": None if foil.source is None else str(foil.source),
               "original_numeric_sha256": scientific_hash(foil.original_points.tolist()),
               "leading_edge_index": foil.leading_edge_index, "reference_method": foil.reference_method}
        for name, foil in (("main", geometry.main), ("flap", geometry.flap))}}


class SolverRunner:
    """One analysis at a time; every call allocates a new exclusive attempt.

    workspace is the explicit file boundary; directory must be relative to it.
    Filesystem failure preventing creation/publication raises OSError rather
    than falsely claiming that an audit record was persisted.
    """

    def __init__(self, workspace: str | Path, directory: str | Path = "outputs/solver_interface_v020") -> None:
        root = native_path(workspace)
        reject_links(root)
        if not root.is_dir():
            raise ValueError("workspace must be an existing directory")
        self.workspace = root.resolve()
        self.execution_root = safe_path(self.workspace, directory)
        self.execution_root.mkdir(parents=True, exist_ok=True)

    def run(self, solver: AerodynamicSolver, geometry: MultiElementAirfoil,
            conditions: OperatingConditions, *, min_clearance: float,
            timeout_s: float | None = None) -> SolverResult:
        """Validate -> prepare -> execute -> parse -> persist, without caching results.

        Timeout covers external execution, not arbitrary in-process Python code
        or disk preparation. An attempt with unrepresentable input has no case_id.
        """
        started_clock, started_at = time.monotonic(), utc_now()
        case, identity, input_error = None, None, None
        safe_timeout = None
        try:
            if timeout_s is not None:
                timeout_s = real_number(timeout_s, "timeout_s")
                if timeout_s <= 0:
                    raise ValueError("timeout_s must be positive")
                safe_timeout = timeout_s
            if not isinstance(solver, AerodynamicSolver):
                raise TypeError("solver must implement AerodynamicSolver")
            identity = solver.identity
            if not isinstance(identity, SolverIdentity):
                raise TypeError("adapter identity must be SolverIdentity")
            case = SolverCase(geometry, conditions, identity, min_clearance)
        except Exception as error:
            input_error = error
            if not isinstance(identity, SolverIdentity):
                identity = None
        timeout_s = safe_timeout
        attempt_id = uuid4().hex
        # The analysis ID lives in manifests, not in the path: a flat attempt
        # layout leaves room for legacy Windows solvers using relative paths.
        attempt_relative = f"attempt-{attempt_id}"
        attempt = safe_path(self.execution_root, attempt_relative)
        attempt.parent.mkdir(parents=True, exist_ok=True)
        attempt.mkdir(exist_ok=False)
        context = AttemptContext(attempt, timeout_s)
        context.write_text("events.jsonl", "")
        context.event("input", "attempt allocated", case_id=None if case is None else case.case_id,
                      attempt_id=attempt_id)
        status, convergence = RunStatus.INVALID_INPUT, ConvergenceStatus.UNKNOWN
        coefficients = AerodynamicCoefficients()
        diagnostics: list[Diagnostic] = []
        validation = None
        execution_duration = None
        stage = "input"
        execution_started = None
        interruption = None
        try:
            if input_error is not None:
                raise input_error
            stage = "prepare"
            context.write_json("input.json", case.to_dict())
            context.write_json("attempt.json", {"attempt_id": attempt_id, "case_id": case.case_id,
                "started_at_utc": started_at, "timeout_s": timeout_s, "initial_state": "allocated"})
            stage = "input"
            report = validate_elements(case.main_points, case.flap_points, min_clearance=case.min_clearance)
            validation = FrozenJSON(report.to_dict())
            stage = "prepare"
            context.write_json("geometry_validation.json", report.to_dict())
            stage = "input"
            if not report.valid:
                raise ValueError("geometry rejected: " + "; ".join(issue.code for issue in report.issues))
            reasons = solver.validate_input(geometry, conditions)
            if not isinstance(reasons, tuple) or not all(isinstance(item, str) for item in reasons):
                raise TypeError("adapter validate_input must return tuple[str, ...]")
            if reasons:
                raise ValueError("adapter rejected input: " + "; ".join(reasons))
            stage = "prepare"
            context.event(stage, "preparation started")
            for name, points in (("main", case.main_points), ("flap", case.flap_points)):
                path = context.path(f"geometry/{name}.dat")
                path.parent.mkdir(exist_ok=True)
                with path.open("x", encoding="utf-8", newline="\n") as stream:
                    np.savetxt(stream, points, fmt="%.17g",
                               header=f"AeroOpt_COORDINATES_v0.1 {name} {case.case_id}", comments="")
            solver.prepare(case, context)
            _check_inputs(case, context)
            if solver.identity != identity:
                raise ValueError("adapter identity changed during preparation")
            stage = "execute"
            context.event(stage, "execution started")
            execution_started = time.monotonic()
            context = replace(context, deadline=None if timeout_s is None else execution_started + timeout_s)
            evidence = solver.execute(case, context)
            execution_duration = time.monotonic() - execution_started
            if not isinstance(evidence, ExecutionEvidence):
                raise TypeError("execute must return ExecutionEvidence")
            context.write_json("execution.json", evidence.to_dict())
            if evidence.timed_out:
                status = RunStatus.TIMEOUT
                diagnostics.append(Diagnostic(stage, "process_timeout", "External execution deadline exceeded"))
            elif evidence.return_code not in (0, None):
                status = RunStatus.EXECUTION_FAILED
                diagnostics.append(Diagnostic(stage, "nonzero_exit", f"Process returned {evidence.return_code}"))
            else:
                stage = "parse"
                context.event(stage, "interpretation started")
                _check_inputs(case, context)
                if solver.identity != identity:
                    raise ValueError("adapter identity changed during execution")
                parsed = solver.parse(case, context, evidence)
                if not isinstance(parsed, ParsedOutput):
                    raise TypeError("parse must return ParsedOutput")
                if identity.is_simulated and not parsed.coefficients.empty:
                    raise ValueError("non-physical mock cannot report aerodynamic coefficients")
                convergence, coefficients = parsed.convergence, parsed.coefficients
                diagnostics.extend(Diagnostic(stage, "adapter_message", item) for item in parsed.messages)
                established = (convergence == ConvergenceStatus.CONVERGED or
                               (convergence == ConvergenceStatus.NOT_APPLICABLE and not identity.convergence_required))
                status = RunStatus.SUCCESS if established else RunStatus.NOT_CONVERGED
                if not established:
                    diagnostics.append(Diagnostic(stage, "convergence_not_established",
                                                  "Process completion does not establish aerodynamic convergence"))
        except (Exception, KeyboardInterrupt, SystemExit) as error:
            if isinstance(error, KeyboardInterrupt):
                interruption = error
            if stage == "execute" and execution_started is not None:
                execution_duration = time.monotonic() - execution_started
            status = {"input": RunStatus.INVALID_INPUT, "prepare": RunStatus.PREPARATION_FAILED,
                      "execute": RunStatus.EXECUTION_FAILED, "parse": RunStatus.PARSING_FAILED}[stage]
            if stage == "execute" and isinstance(error, TimeoutError):
                status = RunStatus.TIMEOUT
            convergence, coefficients = ConvergenceStatus.UNKNOWN, AerodynamicCoefficients()
            diagnostics.append(Diagnostic(stage, type(error).__name__, str(error)))
            context.write_text("traceback.txt", "".join(traceback.format_exception(error)))
            context.event(stage, "exception captured", exception_type=type(error).__name__, detail=str(error))
        context.event("persist", "adapter stages completed", stage_status=status.value)
        artifacts, inventory_problems = _inventory(attempt)
        diagnostics.extend(inventory_problems)
        if inventory_problems:
            status, convergence, coefficients = (RunStatus.EXECUTION_FAILED, ConvergenceStatus.UNKNOWN,
                                                 AerodynamicCoefficients())
        provenance = {"aeroopt_version": __version__, "python": platform.python_version(),
                      "python_executable_display_only": sys.executable, "platform": platform.platform(),
                      "libraries": {name: version(name) for name in ("numpy", "shapely")},
                      "execution_policy": {"timeout_s": timeout_s, "parallelism": "none"},
                      "geometry": None if case is None else _geometry_provenance(geometry)}
        result = SolverResult(
            case_id=None if case is None else case.case_id, attempt_id=attempt_id,
            attempt_directory=attempt_relative, solver=identity,
            conditions=conditions if isinstance(conditions, OperatingConditions) else None,
            status=status, convergence=convergence, coefficients=coefficients,
            diagnostics=tuple(diagnostics), duration_s=time.monotonic() - started_clock,
            execution_duration_s=execution_duration, started_at_utc=started_at, finished_at_utc=utc_now(),
            artifacts=artifacts, provenance=FrozenJSON(provenance),
            input_payload=None if case is None else FrozenJSON(case.to_dict()), geometry_validation=validation)
        result.save(context.path("result.json"))
        if interruption is not None:
            raise interruption
        return result
