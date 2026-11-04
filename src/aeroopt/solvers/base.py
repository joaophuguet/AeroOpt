"""Solver-neutral snapshots and staged adapter contract."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from aeroopt.geometry import MultiElementAirfoil
from aeroopt.geometry._numeric import FloatPoints, immutable_points

from .conditions import OperatingConditions, real_number
from .serialization import INPUT_SCHEMA, FrozenJSON, scientific_hash

if TYPE_CHECKING:
    from .execution import AttemptContext, ExecutionEvidence
    from .results import ParsedOutput


@dataclass(frozen=True)
class SolverIdentity:
    """Adapter identity, explicit simulation label, settings and capabilities.

    Real adapters must include binary/library version and relevant settings.
    Runtime working paths must never be placed in settings.
    """

    identifier: str
    version: str
    is_simulated: bool
    convergence_required: bool
    capabilities: tuple[str, ...]
    settings: FrozenJSON = FrozenJSON()

    def __post_init__(self) -> None:
        for name in ("identifier", "version"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"{name} must be nonempty")
        if type(self.is_simulated) is not bool or type(self.convergence_required) is not bool:
            raise TypeError("simulation/convergence flags must be bool")
        if not isinstance(self.settings, FrozenJSON):
            raise TypeError("settings must be FrozenJSON")
        if isinstance(self.capabilities, (str, bytes)):
            raise TypeError("capabilities must be a sequence of names, not one string")
        capabilities = tuple(self.capabilities)
        if not all(isinstance(item, str) and item.strip() for item in capabilities):
            raise ValueError("capabilities must contain nonempty text")
        object.__setattr__(self, "capabilities", tuple(sorted(set(capabilities))))

    def to_dict(self) -> dict:
        return {"identifier": self.identifier, "version": self.version,
                "is_simulated": self.is_simulated, "convergence_required": self.convergence_required,
                "capabilities": list(self.capabilities), "settings": self.settings.to_dict()}

    @classmethod
    def from_dict(cls, data: dict) -> "SolverIdentity":
        if set(data) != {"identifier", "version", "is_simulated", "convergence_required", "capabilities", "settings"}:
            raise ValueError("unexpected solver-identity fields")
        return cls(data["identifier"], data["version"], data["is_simulated"],
                   data["convergence_required"], tuple(data["capabilities"]), FrozenJSON(data["settings"]))


@dataclass(frozen=True, eq=False)
class SolverCase:
    """Immutable input snapshot, including global points actually hashed/executed."""

    geometry: MultiElementAirfoil
    conditions: OperatingConditions
    solver: SolverIdentity
    min_clearance: float
    schema_version: str = INPUT_SCHEMA
    main_points: FloatPoints = field(init=False, repr=False)
    flap_points: FloatPoints = field(init=False, repr=False)
    case_id: str = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.geometry, MultiElementAirfoil):
            raise TypeError("geometry must be MultiElementAirfoil")
        if not isinstance(self.conditions, OperatingConditions) or not isinstance(self.solver, SolverIdentity):
            raise TypeError("expected OperatingConditions and SolverIdentity")
        clearance = real_number(self.min_clearance, "min_clearance")
        if clearance < 0:
            raise ValueError("min_clearance must be nonnegative")
        if not isinstance(self.schema_version, str) or not self.schema_version:
            raise ValueError("schema_version must be nonempty")
        object.__setattr__(self, "min_clearance", clearance)
        object.__setattr__(self, "main_points", immutable_points(self.geometry.main_points))
        object.__setattr__(self, "flap_points", immutable_points(self.geometry.flap_points))
        object.__setattr__(self, "case_id", scientific_hash(self.to_dict()))

    def to_dict(self) -> dict:
        return {"schema_version": self.schema_version,
                "geometry": {"main": self.main_points.tolist(), "flap": self.flap_points.tolist(),
                             "frame": "global_main_chord", "closure": "straight_TE_cap"},
                "conditions": self.conditions.to_dict(), "solver": self.solver.to_dict(),
                "min_clearance": self.min_clearance}


class AerodynamicSolver(ABC):
    """Trusted adapters stage validation, preparation, execution and interpretation.

    Runner owns persistence, geometry rejection, attempt isolation and error
    classification. Use case.main_points/flap_points as the execution inputs.
    Python-library execution is cooperative; hard timeouts apply to run_process.
    """

    @property
    @abstractmethod
    def identity(self) -> SolverIdentity:
        """Return stable identity/settings for the complete analysis."""

    def validate_input(self, geometry: MultiElementAirfoil,
                       conditions: OperatingConditions) -> tuple[str, ...]:
        """Return adapter-specific rejection reasons; do not repair geometry."""
        return ()

    @abstractmethod
    def prepare(self, case: SolverCase, context: "AttemptContext") -> None:
        """Write solver-specific files without changing the geometry snapshot."""

    @abstractmethod
    def execute(self, case: SolverCase, context: "AttemptContext") -> "ExecutionEvidence":
        """Run an external process or a Python library; return process evidence."""

    @abstractmethod
    def parse(self, case: SolverCase, context: "AttemptContext",
              evidence: "ExecutionEvidence") -> "ParsedOutput":
        """Interpret convergence and available coefficients explicitly."""
