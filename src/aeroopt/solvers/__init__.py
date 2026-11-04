"""Solver interface v0.2: execution infrastructure, no physical solver supplied."""

from .base import AerodynamicSolver, SolverCase, SolverIdentity
from .conditions import ConditionParameter, OperatingConditions
from .execution import AttemptContext, ExecutionEvidence
from .mock import MockScenario, MockSolver
from .results import (AerodynamicCoefficients, Artifact, ConvergenceStatus, Diagnostic,
                      MomentReference, ParsedOutput, RunStatus, SolverResult)
from .runner import SolverRunner
from .serialization import FrozenJSON

__all__ = ["AerodynamicSolver", "SolverCase", "SolverIdentity", "OperatingConditions",
           "ConditionParameter", "AttemptContext", "ExecutionEvidence", "MockScenario", "MockSolver",
           "AerodynamicCoefficients", "Artifact", "ConvergenceStatus", "Diagnostic", "MomentReference",
           "ParsedOutput", "RunStatus", "SolverResult", "SolverRunner", "FrozenJSON"]
