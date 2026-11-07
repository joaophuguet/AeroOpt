"""Reproduce offline numerical verification; no viscous experimental validation."""

import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from aeroopt import __version__
from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points
from aeroopt.solvers import MomentReference, OperatingConditions, PanelSolver, RunStatus, SolverRunner
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers.serialization import read_json, write_json
from _panel_reference import joukowski_exact, joukowski_points


def cylinder_error(count):
    theta = np.linspace(0, 2*np.pi, count+1)
    points = np.column_stack((np.cos(theta), np.sin(theta)))
    delta = np.diff(points, axis=0)
    length = np.linalg.norm(delta, axis=1)
    tangent = delta / length[:,None]
    normal = np.column_stack((tangent[:,1], -tangent[:,0]))
    mesh = pm.PanelMesh(points[:-1], points[1:], (points[:-1]+points[1:])/2, tangent, normal, length, 1, count)
    source, _ = pm.influences((mesh,))
    an = np.einsum("ijk,ik->ij", source, normal)
    at = np.einsum("ijk,ik->ij", source, tangent)
    q = np.linalg.solve(an, -normal[:,0])
    cp = 1 - (at @ q + tangent[:,0])**2
    exact = 1 - 4*np.sin((theta[:-1]+theta[1:])/2)**2
    return {"panels": count, "cp_error_max": float(np.max(np.abs(cp-exact))),
            "normal_residual": float(np.max(np.abs(an @ q + normal[:,0])))}


def geometry(points, position=(1.02, -.04), deflection=25, flap_code="4412", closed=True):
    main = Airfoil.from_points(naca4_points("0012", points, closed_te=closed), leading_edge_index=(points-1)//2)
    flap = Airfoil.from_points(naca4_points(flap_code, points, closed_te=closed), leading_edge_index=(points-1)//2)
    return MultiElementAirfoil(main, flap, FlapParameters(.3, position, deflection))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/panel_verification_v030")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        parser.error("Output must remain inside the AeroOpt workspace")
    if output.exists() and any(output.iterdir()):
        parser.error("Output directory must be empty; choose a new --output path")
    output.mkdir(parents=True, exist_ok=True)
    report = {"aeroopt_version": __version__, "numpy_version": np.__version__,
              "scope": "Numerical verification of inviscid potential flow, not viscous physical validation",
              "cylinder": [cylinder_error(n) for n in (24, 48, 96)], "joukowski": [], "thin_airfoil": []}
    exact = joukowski_exact(4)
    report["joukowski_exact"] = {key: value for key, value in exact.items() if key != "cp"}
    for count in (80, 160, 320, 640, 1920):
        mesh = pm.make_mesh(joukowski_points(count), count//2)
        result = pm.solve(pm.assemble((mesh,), 4), .25, 0, "counterclockwise")
        if not result["accepted"]:
            raise RuntimeError(result["reason"])
        write_json(output / f"joukowski_{count}.json", result)
        report["joukowski"].append({"panels": count, "cl": result["cl"], "cm_ccw": result["cm"],
            "cl_error": abs(result["cl"] - exact["cl"]), "cm_error": abs(result["cm"] - exact["cm_ccw_quarter_chord"]),
            "circulation": result["elements"][0]["circulation"], "diagnostics": result["diagnostics"]})
    for count in (81, 161, 321):
        mesh = pm.make_mesh(naca4_points("0002", count, closed_te=True), (count-1)//2)
        plus = pm.solve(pm.assemble((mesh,), 2), .25, 0, "counterclockwise")
        minus = pm.solve(pm.assemble((mesh,), -2), .25, 0, "counterclockwise")
        slope = float((plus["cl"] - minus["cl"]) / np.deg2rad(4))
        report["thin_airfoil"].append({"panels": mesh.size, "slope_per_rad": slope,
                                      "relative_difference_from_2pi": slope/(2*np.pi)-1})
    runner = SolverRunner(ROOT, (output / "attempts").relative_to(ROOT))
    solver = PanelSolver(MomentReference(.25, 0, "clockwise", "q_inf*c_main^2"))
    conditions = OperatingConditions(6, 200000, 0)

    def run(label, shape, op=conditions):
        result = runner.run(solver, shape, op, min_clearance=.002)
        if result.status != RunStatus.SUCCESS:
            raise RuntimeError(str(result.diagnostics))
        data = read_json(runner.execution_root / result.attempt_directory / "panel_output.json")
        return {"label": label, "case_id": result.case_id, "attempt_directory": result.attempt_directory,
                "panels": sum(e["panel_count"] for e in data["elements"]),
                "cl": data["cl"], "cm_clockwise": data["cm"],
                "circulations_ccw": [e["circulation"] for e in data["elements"]], "diagnostics": data["diagnostics"]}

    report["multielement_refinement"] = [run(f"refinement_{n}", geometry(n)) for n in (41, 81, 161, 321)]
    report["finite_te_refinement"] = [run(f"finite_te_{n}", geometry(n, closed=False)) for n in (41, 81, 161, 321)]
    report["sensitivity"] = [run("straight", geometry(161, deflection=0)),
                             run("deflected_20", geometry(161, deflection=20)),
                             run("displaced", geometry(161, position=(1.2, -.1), deflection=20))]
    report["symmetric_zero"] = run("symmetric_zero", geometry(161, position=(1.05, 0), deflection=0, flap_code="0012"),
                                   OperatingConditions(0, 200000, 0))
    write_json(output / "verification.json", report)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    rows = report["joukowski"]
    axes[0].loglog([r["panels"] for r in rows], [r["cl_error"] for r in rows], "o-", label="CL error")
    axes[0].loglog([r["panels"] for r in rows], [r["cm_error"] for r in rows], "s-", label="CM error")
    axes[0].set(xlabel="Panels", ylabel="Absolute coefficient error", title="Exact Joukowski: slow cusp convergence")
    axes[0].legend()
    rows = report["multielement_refinement"]
    axes[1].loglog([r["panels"] for r in rows], [abs(r["diagnostics"]["pressure_circulation_cl_difference"]) for r in rows], "o-")
    axes[1].set(xlabel="Total main + flap panels", ylabel="|CL pressure - CL circulation|", title="Coupled NACA elements: consistency")
    axes[1].set_xticks([80, 160, 320, 640], labels=["80", "160", "320", "640"])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    for axis in axes:
        axis.grid(alpha=.25, which="both")
    fig.savefig(output / "refinement.png", dpi=160)
    plt.close(fig)
    print(f"Numerical verification written to {output / 'verification.json'}")
    print("No experimental viscous validation or stall/drag predictions.")


if __name__ == "__main__":
    main()
