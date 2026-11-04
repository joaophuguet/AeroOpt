"""AeroOpt v0.2.0: NON-PHYSICAL solver execution and failure demonstration."""

import argparse
from pathlib import Path

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points
from aeroopt.solvers import MockScenario, MockSolver, OperatingConditions, RunStatus, SolverResult, SolverRunner
from aeroopt.solvers._paths import native_path, safe_path
from aeroopt.solvers.serialization import write_json

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="outputs/solver_interface_v020",
                        help="New/empty directory relative to the project workspace")
    parser.add_argument("--timeout", type=float, default=.5, help="External execution deadline, seconds")
    args = parser.parse_args()
    output = safe_path(native_path(ROOT), args.output)
    if output.exists() and any(output.iterdir()):
        parser.error("Output must be new/empty; previous evidence is never overwritten")
    runner = SolverRunner(ROOT, args.output)
    main = Airfoil.from_points(naca4_points("0012", 161), leading_edge_index=80,
                              name="ANALYTIC_NACA_0012")
    flap = Airfoil.from_points(naca4_points("4412", 81), leading_edge_index=40,
                              name="ANALYTIC_NACA_4412")
    geometry = MultiElementAirfoil(main, flap, FlapParameters(.3, (1.02, -.04), 25.))
    conditions = OperatingConditions(alpha_deg=6., reynolds=200000., mach=.05)
    print("NON-PHYSICAL infrastructure demonstration. cl, cd and cm remain absent.")
    expected = {MockScenario.SUCCESS: RunStatus.SUCCESS, MockScenario.INPUT_REJECTED: RunStatus.INVALID_INPUT,
                MockScenario.PREPARATION_FAILED: RunStatus.PREPARATION_FAILED,
                MockScenario.EXECUTION_FAILED: RunStatus.EXECUTION_FAILED, MockScenario.TIMEOUT: RunStatus.TIMEOUT,
                MockScenario.PARSING_FAILED: RunStatus.PARSING_FAILED, MockScenario.NOT_CONVERGED: RunStatus.NOT_CONVERGED}
    results = []
    for scenario in MockScenario:
        result = runner.run(MockSolver(scenario), geometry, conditions, min_clearance=.002, timeout_s=args.timeout)
        if result.status != expected[scenario]:
            raise AssertionError(f"Unexpected {scenario} result: {result.to_dict()}")
        results.append((scenario.value, result))
    retry = runner.run(MockSolver(), geometry, conditions, min_clearance=.002, timeout_s=args.timeout)
    assert retry.case_id == results[0][1].case_id and retry.attempt_id != results[0][1].attempt_id
    assert retry.status == RunStatus.SUCCESS
    results.append(("same_analysis_retry", retry))
    invalid = geometry.configure(flap_deflection_deg=0., flap_dx=-.6)
    rejected = runner.run(MockSolver(), invalid, conditions, min_clearance=.002, timeout_s=args.timeout)
    assert rejected.status == RunStatus.INVALID_INPUT
    results.append(("invalid_geometry", rejected))
    summary = []
    for label, result in results:
        relative = f"{result.attempt_directory}/result.json"
        loaded = SolverResult.load(safe_path(runner.execution_root, relative))
        assert loaded.to_dict() == result.to_dict()
        assert result.coefficients.empty and not result.physical_coefficients_available
        summary.append({"label": label, "case_id": result.case_id, "attempt_id": result.attempt_id,
                        "status": result.status.value, "convergence": result.convergence.value,
                        "coefficients": result.coefficients.to_dict(), "result_file": relative})
        print(f"{label:22} {result.status.value:20} analysis={result.case_id[:12]} attempt={result.attempt_id[:8]}")
    write_json(output / "summary.json", {"non_physical": True, "physical_coefficients_available": False,
                                        "conditions": conditions.to_dict(), "runs": summary})
    print(f"\n{len(results)} attempts persisted. Results and logs: {args.output}")


if __name__ == "__main__":
    main()
