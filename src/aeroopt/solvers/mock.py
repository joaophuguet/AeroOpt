"""Deterministic NON-PHYSICAL adapter. All aerodynamic coefficients stay absent."""

import sys
import time
from dataclasses import dataclass
from enum import StrEnum

from .base import AerodynamicSolver, SolverCase, SolverIdentity
from .conditions import OperatingConditions
from .execution import AttemptContext, ExecutionEvidence
from .results import AerodynamicCoefficients, ConvergenceStatus, ParsedOutput
from .serialization import FrozenJSON, read_json


class MockScenario(StrEnum):
    SUCCESS = "success"
    INPUT_REJECTED = "input_rejected"
    PREPARATION_FAILED = "preparation_failed"
    EXECUTION_FAILED = "execution_failed"
    TIMEOUT = "timeout"
    PARSING_FAILED = "parsing_failed"
    NOT_CONVERGED = "not_converged"


@dataclass(frozen=True)
class MockSolver(AerodynamicSolver):
    """Infrastructure test scenarios, not an aerodynamic model or surrogate."""

    scenario: MockScenario = MockScenario.SUCCESS

    def __post_init__(self) -> None:
        object.__setattr__(self, "scenario", MockScenario(self.scenario))

    @property
    def identity(self) -> SolverIdentity:
        return SolverIdentity("aeroopt.non_physical_mock", "0.2.0", True, False,
            ("2d_multielement_input", "infrastructure_tests_only", "no_aerodynamic_predictions"),
            FrozenJSON({"scenario": self.scenario.value}))

    def validate_input(self, geometry, conditions: OperatingConditions) -> tuple[str, ...]:
        return ("NON-PHYSICAL mock input rejection",) if self.scenario == MockScenario.INPUT_REJECTED else ()

    def prepare(self, case: SolverCase, context: AttemptContext) -> None:
        context.write_json("mock_request.json", {"case_id": case.case_id, "scenario": self.scenario.value,
                                               "non_physical": True})
        if self.scenario == MockScenario.PREPARATION_FAILED:
            context.write_text("partial_prepare.txt", "NON-PHYSICAL preparation evidence\n")
            raise RuntimeError("NON-PHYSICAL simulated preparation failure")

    def execute(self, case: SolverCase, context: AttemptContext) -> ExecutionEvidence:
        started = time.monotonic()
        context.write_text("adapter.log", "NON-PHYSICAL MOCK: no aerodynamic calculations or coefficients.\n")
        if self.scenario == MockScenario.EXECUTION_FAILED:
            raise RuntimeError("NON-PHYSICAL simulated execution failure")
        if self.scenario == MockScenario.TIMEOUT:
            if context.timeout_s is None:
                raise ValueError("timeout mock scenario requires an explicit finite timeout_s")
            return context.run_process([sys.executable, "-u", "-c",
                "import time; print('NON_PHYSICAL_MOCK_BLOCK', flush=True); time.sleep(3600)"])
        if self.scenario == MockScenario.PARSING_FAILED:
            context.write_text("mock_output.json", "{deliberately invalid NON-PHYSICAL JSON")
        else:
            convergence = (ConvergenceStatus.NOT_CONVERGED if self.scenario == MockScenario.NOT_CONVERGED
                           else ConvergenceStatus.NOT_APPLICABLE)
            context.write_json("mock_output.json", {"non_physical": True, "convergence": convergence.value,
                                                   "coefficients": AerodynamicCoefficients().to_dict()})
        return ExecutionEvidence(0, time.monotonic() - started)

    def parse(self, case: SolverCase, context: AttemptContext, evidence: ExecutionEvidence) -> ParsedOutput:
        data = read_json(context.path("mock_output.json"))
        if data.get("non_physical") is not True:
            raise ValueError("missing NON-PHYSICAL mock label")
        return ParsedOutput(ConvergenceStatus(data["convergence"]),
                            AerodynamicCoefficients.from_dict(data["coefficients"]),
                            ("NON-PHYSICAL infrastructure test; no physical results.",))
