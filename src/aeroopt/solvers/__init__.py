"""Auditable solver infrastructure and the v0.3 inviscid panel adapter."""

from .base import AerodynamicSolver, SolverCase, SolverIdentity
from .conditions import ConditionParameter, OperatingConditions
from .execution import AttemptContext, ExecutionEvidence
from .mock import MockScenario, MockSolver
from .panel import PanelSolver
from .results import (AerodynamicCoefficients, Artifact, ConvergenceStatus, Diagnostic,
                      MomentReference, ParsedOutput, RunStatus, SolverResult)
from .runner import SolverRunner
from .serialization import FrozenJSON

__all__ = ["AerodynamicSolver", "SolverCase", "SolverIdentity", "OperatingConditions",
           "ConditionParameter", "AttemptContext", "ExecutionEvidence", "MockScenario", "MockSolver",
           "AerodynamicCoefficients", "Artifact", "ConvergenceStatus", "Diagnostic", "MomentReference",
           "ParsedOutput", "RunStatus", "SolverResult", "SolverRunner", "FrozenJSON", "PanelSolver"]
