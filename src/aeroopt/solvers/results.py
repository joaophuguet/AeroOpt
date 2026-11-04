"""Scientific result contracts: absence != zero; simulated coefficients forbidden."""

import re
from dataclasses import dataclass
from datetime import datetime, timedelta
from enum import StrEnum
from pathlib import Path

from ._paths import relative_path
from .base import SolverIdentity
from .conditions import OperatingConditions, real_number
from .serialization import RESULT_SCHEMA, FrozenJSON, read_json, scientific_hash, write_json


class RunStatus(StrEnum):
    SUCCESS = "success"
    INVALID_INPUT = "invalid_input"
    PREPARATION_FAILED = "preparation_failed"
    EXECUTION_FAILED = "execution_failed"
    TIMEOUT = "timeout"
    PARSING_FAILED = "parsing_failed"
    NOT_CONVERGED = "not_converged"


class ConvergenceStatus(StrEnum):
    UNKNOWN = "unknown"
    NOT_APPLICABLE = "not_applicable"
    CONVERGED = "converged"
    NOT_CONVERGED = "not_converged"


@dataclass(frozen=True)
class MomentReference:
    """Explicit global location, positive rotation and cm normalization definition."""

    x_over_c: float
    y_over_c: float
    positive_sense: str
    normalization: str

    def __post_init__(self) -> None:
        for name in ("x_over_c", "y_over_c"):
            object.__setattr__(self, name, real_number(getattr(self, name), name))
        if self.positive_sense not in ("clockwise", "counterclockwise"):
            raise ValueError("moment positive_sense must be explicit")
        if not isinstance(self.normalization, str) or not self.normalization.strip():
            raise ValueError("moment normalization must be explicit")

    def to_dict(self) -> dict:
        return {"x_over_c": self.x_over_c, "y_over_c": self.y_over_c,
                "positive_sense": self.positive_sense, "normalization": self.normalization}


@dataclass(frozen=True)
class AerodynamicCoefficients:
    """Optional finite values. cm requires an explicit reference, never a default."""

    cl: float | None = None
    cd: float | None = None
    cm: float | None = None
    moment_reference: MomentReference | None = None

    def __post_init__(self) -> None:
        for name in ("cl", "cd", "cm"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, real_number(getattr(self, name), name))
        if self.moment_reference is not None and not isinstance(self.moment_reference, MomentReference):
            raise TypeError("moment_reference must be MomentReference")
        if self.cm is not None and self.moment_reference is None:
            raise ValueError("cm requires an explicit moment reference and normalization")

    @property
    def empty(self) -> bool:
        return self.cl is None and self.cd is None and self.cm is None

    def to_dict(self) -> dict:
        return {"cl": self.cl, "cd": self.cd, "cm": self.cm,
                "moment_reference": None if self.moment_reference is None else self.moment_reference.to_dict()}

    @classmethod
    def from_dict(cls, data: dict) -> "AerodynamicCoefficients":
        if set(data) != {"cl", "cd", "cm", "moment_reference"}:
            raise ValueError("unexpected coefficients fields")
        reference = data["moment_reference"]
        return cls(data["cl"], data["cd"], data["cm"],
                   None if reference is None else MomentReference(**reference))


@dataclass(frozen=True)
class ParsedOutput:
    convergence: ConvergenceStatus
    coefficients: AerodynamicCoefficients = AerodynamicCoefficients()
    messages: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.convergence, ConvergenceStatus):
            raise TypeError("convergence must be ConvergenceStatus")
        if not isinstance(self.coefficients, AerodynamicCoefficients):
            raise TypeError("coefficients must be AerodynamicCoefficients")
        if isinstance(self.messages, (str, bytes)) or not all(isinstance(item, str) for item in self.messages):
            raise TypeError("messages must be text")
        object.__setattr__(self, "messages", tuple(self.messages))


@dataclass(frozen=True)
class Diagnostic:
    stage: str
    code: str
    message: str

    def __post_init__(self) -> None:
        if not isinstance(self.message, str):
            raise TypeError("diagnostic message must be text")
        if not all(isinstance(value, str) and value for value in (self.stage, self.code)):
            raise ValueError("diagnostic stage/code must be nonempty text")

    def to_dict(self) -> dict:
        return {"stage": self.stage, "code": self.code, "message": self.message}


@dataclass(frozen=True)
class Artifact:
    path: str
    sha256: str
    size_bytes: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", relative_path(self.path))
        if not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
            raise ValueError("artifact hash must be SHA-256 hex")
        if type(self.size_bytes) is not int or self.size_bytes < 0:
            raise ValueError("artifact size must be a nonnegative integer")

    def to_dict(self) -> dict:
        return {"path": self.path, "sha256": self.sha256, "size_bytes": self.size_bytes}


@dataclass(frozen=True)
class SolverResult:
    case_id: str | None
    attempt_id: str
    attempt_directory: str
    solver: SolverIdentity | None
    conditions: OperatingConditions | None
    status: RunStatus
    convergence: ConvergenceStatus
    coefficients: AerodynamicCoefficients
    diagnostics: tuple[Diagnostic, ...]
    duration_s: float
    execution_duration_s: float | None
    started_at_utc: str
    finished_at_utc: str
    artifacts: tuple[Artifact, ...]
    provenance: FrozenJSON
    input_payload: FrozenJSON | None
    geometry_validation: FrozenJSON | None

    def __post_init__(self) -> None:
        if self.case_id is not None and not re.fullmatch(r"[0-9a-f]{64}", self.case_id):
            raise ValueError("case_id must be SHA-256 hex or absent for unrepresentable input")
        if not re.fullmatch(r"[0-9a-f]{32}", self.attempt_id):
            raise ValueError("attempt_id must be uuid4 hex")
        object.__setattr__(self, "attempt_directory", relative_path(self.attempt_directory))
        if not isinstance(self.status, RunStatus) or not isinstance(self.convergence, ConvergenceStatus):
            raise TypeError("status/convergence must use their enums")
        if not isinstance(self.coefficients, AerodynamicCoefficients):
            raise TypeError("invalid coefficient record")
        for name in ("duration_s", "execution_duration_s"):
            value = getattr(self, name)
            if value is not None:
                value = real_number(value, name)
                if value < 0:
                    raise ValueError("durations must be nonnegative")
                object.__setattr__(self, name, value)
        if self.duration_s is None:
            raise ValueError("duration_s is required")
        if self.solver is not None and not isinstance(self.solver, SolverIdentity):
            raise TypeError("solver must be SolverIdentity")
        if self.conditions is not None and not isinstance(self.conditions, OperatingConditions):
            raise TypeError("conditions must be OperatingConditions")
        if self.solver is not None and self.solver.is_simulated and not self.coefficients.empty:
            raise ValueError("simulated adapters must never publish aerodynamic coefficients")
        if self.status not in (RunStatus.SUCCESS, RunStatus.NOT_CONVERGED) and not self.coefficients.empty:
            raise ValueError("failed executions cannot publish coefficients")
        if self.status == RunStatus.SUCCESS:
            if self.solver is None or self.convergence == ConvergenceStatus.UNKNOWN:
                raise ValueError("success requires a declared solver and known convergence policy")
            if self.convergence == ConvergenceStatus.NOT_CONVERGED:
                raise ValueError("process success is not convergence")
            if self.solver.convergence_required and self.convergence != ConvergenceStatus.CONVERGED:
                raise ValueError("this adapter requires confirmed convergence")
        if self.status == RunStatus.NOT_CONVERGED and self.convergence == ConvergenceStatus.CONVERGED:
            raise ValueError("not_converged status contradicts confirmed convergence")
        if self.status not in (RunStatus.SUCCESS, RunStatus.NOT_CONVERGED) and self.convergence != ConvergenceStatus.UNKNOWN:
            raise ValueError("failed stages cannot assert convergence")
        for name in ("started_at_utc", "finished_at_utc"):
            value = getattr(self, name)
            if not isinstance(value, str) or datetime.fromisoformat(value).utcoffset() != timedelta(0):
                raise ValueError("execution timestamps must explicitly use UTC")
        for name, cls in (("diagnostics", Diagnostic), ("artifacts", Artifact)):
            values = tuple(getattr(self, name))
            if not all(isinstance(item, cls) for item in values):
                raise TypeError(f"invalid {name}")
            object.__setattr__(self, name, values)
        if not isinstance(self.provenance, FrozenJSON):
            raise TypeError("provenance must be FrozenJSON")
        if self.geometry_validation is not None and not isinstance(self.geometry_validation, FrozenJSON):
            raise TypeError("geometry_validation must be FrozenJSON")
        if self.input_payload is not None:
            if not isinstance(self.input_payload, FrozenJSON):
                raise TypeError("input_payload must be FrozenJSON")
            data = self.input_payload.to_dict()
            if scientific_hash(data) != self.case_id:
                raise ValueError("case identity does not match input snapshot")
            if self.solver is None or data["solver"] != self.solver.to_dict():
                raise ValueError("solver identity does not match input snapshot")
            if self.conditions is None or data["conditions"] != self.conditions.to_dict():
                raise ValueError("conditions do not match input snapshot")
        elif self.case_id is not None:
            raise ValueError("identified case requires its input snapshot")

    @property
    def physical_coefficients_available(self) -> bool:
        return (self.status == RunStatus.SUCCESS and self.solver is not None
                and not self.solver.is_simulated and not self.coefficients.empty)

    def to_dict(self) -> dict:
        return {"schema_version": RESULT_SCHEMA, "case_id": self.case_id,
                "attempt_id": self.attempt_id, "attempt_directory": self.attempt_directory,
                "solver": None if self.solver is None else self.solver.to_dict(),
                "conditions": None if self.conditions is None else self.conditions.to_dict(),
                "status": self.status.value, "convergence": self.convergence.value,
                "coefficients": self.coefficients.to_dict(),
                "physical_coefficients_available": self.physical_coefficients_available,
                "diagnostics": [item.to_dict() for item in self.diagnostics],
                "duration_s": self.duration_s, "execution_duration_s": self.execution_duration_s,
                "started_at_utc": self.started_at_utc, "finished_at_utc": self.finished_at_utc,
                "artifacts": [item.to_dict() for item in self.artifacts],
                "provenance": self.provenance.to_dict(),
                "input_payload": None if self.input_payload is None else self.input_payload.to_dict(),
                "geometry_validation": None if self.geometry_validation is None else self.geometry_validation.to_dict()}

    def save(self, path: str | Path) -> None:
        write_json(path, self.to_dict())

    @classmethod
    def load(cls, path: str | Path) -> "SolverResult":
        data = read_json(path)
        if data.pop("schema_version") != RESULT_SCHEMA:
            raise ValueError("unsupported result schema")
        physical = data.pop("physical_coefficients_available")
        data["solver"] = None if data["solver"] is None else SolverIdentity.from_dict(data["solver"])
        data["conditions"] = None if data["conditions"] is None else OperatingConditions.from_dict(data["conditions"])
        data["status"] = RunStatus(data["status"])
        data["convergence"] = ConvergenceStatus(data["convergence"])
        data["coefficients"] = AerodynamicCoefficients.from_dict(data["coefficients"])
        data["diagnostics"] = tuple(Diagnostic(**item) for item in data["diagnostics"])
        data["artifacts"] = tuple(Artifact(**item) for item in data["artifacts"])
        data["provenance"] = FrozenJSON(data["provenance"])
        for name in ("input_payload", "geometry_validation"):
            data[name] = None if data[name] is None else FrozenJSON(data[name])
        result = cls(**data)
        if type(physical) is not bool or physical != result.physical_coefficients_available:
            raise ValueError("inconsistent physical-result label")
        return result
