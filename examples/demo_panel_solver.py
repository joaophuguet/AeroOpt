"""Run the existing coarse synthetic demonstration geometry without resampling."""

import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aeroopt.geometry import MultiElementAirfoil
from aeroopt.solvers import MomentReference, OperatingConditions, PanelSolver, RunStatus, SolverRunner
from aeroopt.solvers.serialization import read_json, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/panel_demo_v030")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        parser.error("Output must remain inside the AeroOpt workspace")
    if output.exists() and any(output.iterdir()):
        parser.error("Output directory must be empty; choose a new --output path")
    geometry = MultiElementAirfoil.from_dat(
        ROOT / "examples/data/demo_main.dat", ROOT / "examples/data/demo_flap.dat",
        flap_chord_ratio=.3, flap_reference_position=(1.10, -.05),
    ).configure(flap_deflection_deg=25)
    solver = PanelSolver(MomentReference(.25, 0, "clockwise", "q_inf*c_main^2"))
    runner = SolverRunner(ROOT, output.relative_to(ROOT))
    result = runner.run(solver, geometry, OperatingConditions(6, 200000, 0), min_clearance=.002)
    summary = {"status": result.status.value, "case_id": result.case_id,
               "attempt_directory": result.attempt_directory, "coefficients": result.coefficients.to_dict(),
               "warning": "Coarse synthetic contours: execution example, not a converged aerodynamic prediction."}
    if result.status != RunStatus.SUCCESS:
        write_json(output / "summary.json", summary)
        raise RuntimeError(str(result.diagnostics))
    panel_output = read_json(runner.execution_root / result.attempt_directory / "panel_output.json")
    summary["diagnostics"] = panel_output["diagnostics"]
    write_json(output / "summary.json", summary)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), layout="constrained")
    geometry.plot(ax=axes[0])
    axes[0].set_title("Existing synthetic contours: 9 main / 7 flap points, unchanged")
    for name, element in zip(("Main", "Flap"), panel_output["elements"], strict=True):
        centers = element["centers"]
        # Disconnected markers also display the finite-TE base panel explicitly.
        axes[1].plot([p[0] for p in centers], element["cp"], "o", label=f"{name} panel centers")
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Global x / main chord", ylabel="Cp", title="Inviscid pressure; coarse panels, no stall or drag model")
    axes[1].legend()
    axes[1].grid(alpha=.25)
    fig.savefig(output / "geometry_cp.png", dpi=160)
    plt.close(fig)
    print(f"{result.status}: CL={result.coefficients.cl:.8f}; CM(clockwise, quarter chord)={result.coefficients.cm:.8f}; CD=None")
    print(summary["warning"])
    print(f"Artifacts: {output}")


if __name__ == "__main__":
    main()
