"""Reproduce the v0.3.1 scientific review with explicit research variants.

Independent reference: linear nodal vorticity, not another Hess--Smith call.
Research variants never modify production inputs or replace production Kutta.
All arrays, provenance, residuals and performance measurements are preserved.
"""

import argparse
from decimal import Decimal, localcontext
import hashlib
import os
from pathlib import Path
import platform
import time
import tracemalloc

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from aeroopt import __version__
from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_surfaces
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers.serialization import scientific_hash, write_json
from _linear_vortex_reference import solve_linear_vortex
from _panel_reference import joukowski_exact
from _panel_review_tools import (assemble_dense_vectors, evaluate_strengths, extrapolated_kutta,
                                 integrate_pressure, joukowski_blasius, joukowski_sample, subdivide_cap)
from verify_panel_solver import geometry


def finite_geometry(points, te_power):
    """Explicit analytic station family x=1-(1-x_cosine)^p; no imported data."""
    x = (1-np.cos(np.linspace(0,np.pi,(points+1)//2)))/2
    x = 1-(1-x)**te_power
    foils = []
    for code in ("0012","4412"):
        upper,lower = naca4_surfaces(code,x,closed_te=False)
        p = np.vstack((upper[::-1],lower[1:]))
        foils.append(Airfoil.from_points(p,leading_edge_index=(points-1)//2))
    return MultiElementAirfoil(*foils,FlapParameters(.3,(1.02,-.04),25.))


def meshes_of(shape):
    return (pm.make_mesh(shape.main_points,shape.main.leading_edge_index),
            pm.make_mesh(shape.flap_points,shape.flap.leading_edge_index))


def solve_hs(meshes, alpha=6.):
    system = pm.assemble(meshes,alpha)
    result = pm.solve(system,.25,0.,"counterclockwise")
    if not result["accepted"]:
        raise RuntimeError(result)
    return system,result


def record_case(output,label,points,alpha,provenance,**results):
    inputs = {"coordinates": [np.asarray(p).tolist() for p in points], "alpha_deg": float(alpha),
              "moment_reference": {"xy": [.25,0.],"sense": "CCW","normalization": "q_inf*c_main^2"},
              "provenance": provenance}
    identifier = scientific_hash(inputs)
    write_json(output / "cases" / f"{label}.json", {"input": inputs,"input_sha256": identifier,**results})
    return identifier


def summary_hs(result):
    return {"cl": result["cl"], "cm_ccw": result["cm"],
            "circulations": [e["circulation"] for e in result["elements"]],
            "source_fluxes": [e["source_flux"] for e in result["elements"]], "diagnostics": result["diagnostics"]}


def compare_multielement(output,label,shape,alpha=6.):
    contours = (shape.main_points,shape.flap_points)
    _,hs = solve_hs(meshes_of(shape),alpha)
    lv = solve_linear_vortex(contours,alpha)
    report = {"label": label,"panels": sum(e["panel_count"] for e in hs["elements"]),
              "hess_smith": summary_hs(hs),
              "linear_vortex": {key:val for key,val in lv.items() if key != "elements"},
              "cp_rms_difference": [],"cp_length_weighted_rms_difference": [],
              "cl_difference": hs["cl"]-lv["cl"],"cm_difference": hs["cm"]-lv["cm_ccw"],
              "geometric_clearance_c_main": shape.validate(min_clearance=.002).min_distance,
              "parameters": shape.parameters.to_dict()}
    for a,b in zip(hs["elements"],lv["elements"],strict=True):
        diff = np.asarray(a["cp"])-b["cp"]
        length = np.asarray(a["lengths"])
        report["cp_rms_difference"].append(float(np.sqrt(np.mean(diff**2))))
        report["cp_length_weighted_rms_difference"].append(float(np.sqrt(np.sum(length*diff**2)/sum(length))))
    report["input_sha256"] = record_case(output,label,contours,alpha,
        {"generator": "analytic NACA, closed_te=True, explicit nominal LE", "parameters": shape.parameters.to_dict()},
        hess_smith=hs,linear_vortex=lv)
    return report


def far_field_review():
    rows = []
    for distance in (1e6,1e9,1e12):
        with localcontext() as context:
            context.prec = 60
            d = Decimal(int(distance))
            reference = float((d/(d-1)).ln())/(2*np.pi)
        old = float((np.log(distance)-np.log(distance-1))/(2*np.pi))
        new = float(pm.source_velocity(np.array([[distance,0.]]),np.array([[0.,0.]]),np.array([[1.,0.]]))[0,0,0])
        rows.append({"distance": distance,"reference": reference,"legacy_log_subtraction": old,"balanced_result": new,
                     "legacy_relative_error": abs(old/reference-1),"balanced_relative_error": abs(new/reference-1)})
    return rows


def benchmark(points,repetitions):
    meshes = meshes_of(geometry(points))
    rows = {}
    systems = {}
    for name,assemble in (("dense_vectors_v030",assemble_dense_vectors),("scalar_projection_v031",pm.assemble)):
        measurements = []
        for repeat in range(repetitions):
            start = time.perf_counter()
            system = assemble(meshes,6.)
            built = time.perf_counter()
            condition = float(np.linalg.cond(system.matrix))
            conditioned = time.perf_counter()
            strengths = np.linalg.solve(system.matrix,system.rhs)
            solved = time.perf_counter()
            measurements.append({"assembly_s": built-start,"svd_condition_s": conditioned-built,
                                 "linear_solve_s": solved-conditioned,"total_core_s": solved-start})
        # Profile allocations separately: tracing disproportionately slows the
        # original Python column loop and must not contaminate runtime claims.
        tracemalloc.start()
        profiled = assemble(meshes,6.)
        np.linalg.cond(profiled.matrix)
        np.linalg.solve(profiled.matrix,profiled.rhs)
        _,peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        systems[name] = system
        rows[name] = {"condition_2": condition,"measurements": measurements,
                      "peak_traced_allocation_bytes": peak,
                      "median": {key: float(np.median([r[key] for r in measurements])) for key in measurements[0]},
                      "coefficients": evaluate_strengths(system,strengths)}
    old,new = systems.values()
    rows["panel_count"] = sum(m.size for m in meshes)
    rows["matrix_max_abs_difference"] = float(np.max(np.abs(old.matrix-new.matrix)))
    rows["rhs_max_abs_difference"] = float(np.max(np.abs(old.rhs-new.rhs)))
    rows["scope"] = "core assembly + full SVD + solve; no runner, I/O or geometry validation"
    return rows


def inverse_bound_experiment(points):
    """Not adopted: conservative inverse-based bound, with a residual check."""
    system = pm.assemble(meshes_of(geometry(points)),6.)
    a = system.matrix
    start = time.perf_counter()
    inverse = np.linalg.solve(a,np.eye(len(a)))
    inverse_time = time.perf_counter()-start
    start = time.perf_counter()
    residual = float(np.linalg.norm(np.eye(len(a))-a@inverse,"fro"))
    bound = float(np.linalg.norm(a,"fro")*np.linalg.norm(inverse,"fro")/(1-residual)) if residual < 1 else None
    product_time = time.perf_counter()-start
    start = time.perf_counter()
    condition = float(np.linalg.cond(a))
    svd_time = time.perf_counter()-start
    return {"panels": len(a)-2,"inverse_s": inverse_time,"residual_and_bound_s": product_time,
            "svd_s": svd_time,"inverse_residual_fro": residual,"condition_2": condition,"frobenius_upper_bound": bound,
            "decision": "not adopted; extra dense storage/product, conservative bound, floating-point bound not interval-certified"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"outputs/panel_review_v031")
    parser.add_argument("--quick",action="store_true",help="Small smoke study, not the report's full evidence")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        parser.error("Output must remain inside the workspace")
    if output.exists() and any(output.iterdir()):
        parser.error("Output directory must be empty")
    (output/"cases").mkdir(parents=True,exist_ok=True)
    code_paths = (Path(__file__),Path(__file__).with_name("_linear_vortex_reference.py"),
                  Path(__file__).with_name("_panel_review_tools.py"),Path(pm.__file__),
                  Path(__file__).with_name("_panel_reference.py"),Path(__file__).with_name("verify_panel_solver.py"))
    (output/"sources").mkdir()
    for path in code_paths:
        with (output/"sources"/path.name).open("xb") as stream:
            stream.write(path.read_bytes())
    exact = joukowski_exact(4)
    report = {"version": __version__,"scope": "software and numerical verification; no experimental physical validation",
              "quick": args.quick,"python": platform.python_version(),"platform": platform.platform(),"numpy": np.__version__,
              "threads_environment": {k:os.environ.get(k) for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")},
              "code_sha256": {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in code_paths},
              "joukowski_exact": {key:val for key,val in exact.items() if key != "cp"},
              "joukowski_blasius": joukowski_blasius(4),"far_field_source": far_field_review(),
              "joukowski": [],"finite_te": [],"finite_te_density": [],"benchmarks": []}
    counts = (80,160) if args.quick else (80,160,320,640,1920)
    for skew in (0.,-.5):
        for count in counts:
            points,exact_cp = joukowski_sample(count,skew)
            mesh = pm.make_mesh(points,count//2)
            system,hs = solve_hs((mesh,),4.)
            try:
                lv = solve_linear_vortex((points,),4.)
            except ArithmeticError as error:
                # A reference is not exempt from conditioning checks. Preserve
                # its rejection, without inventing a value at finer cusp grids.
                lv = {"accepted": False,"reason": str(error)}
            analytic_cp_polygon = integrate_pressure((mesh,),exact_cp,4.)
            a,b = system.matrix.copy(),system.rhs.copy()
            a[-1] = 0
            a[-1,-1] = 1
            b[-1] = exact["circulation"]/sum(mesh.lengths)
            exact_gamma = evaluate_strengths(system,np.linalg.solve(a,b))
            extrapolated = extrapolated_kutta(system)
            label = f"joukowski_{count}_bias{int(-10*skew)}"
            identifier = record_case(output,label,(points,),4.,
                {"generator": "Joukowski circle", "epsilon": .08,"theta_mapping": "beta+skew*sin(beta)","skew": skew},
                hess_smith=hs,linear_vortex=lv,exact_cp_at_parameter_midpoints=exact_cp.tolist())
            report["joukowski"].append({"panels": count,"skew": skew,"input_sha256": identifier,
                "hess_smith": summary_hs(hs),"linear_vortex": {k:v for k,v in lv.items() if k != "elements"},
                "cp_comparison_scope": "exact Cp at parameter midpoints versus polygon collocation; cusp values not qualified",
                "hess_smith_cp_rms_difference": float(np.sqrt(np.mean((np.asarray(hs["elements"][0]["cp"])-exact_cp)**2))),
                "linear_vortex_cp_rms_difference": (float(np.sqrt(np.mean((np.asarray(lv["elements"][0]["cp"])-exact_cp)**2)))
                                                     if lv["accepted"] else None),
                "analytic_cp_polygon": analytic_cp_polygon,"forced_exact_circulation_diagnostic": exact_gamma,
                "extrapolated_kutta_diagnostic": extrapolated,
                "cl_error": abs(hs["cl"]-exact["cl"]),"cm_error": abs(hs["cm"]-exact["cm_ccw_quarter_chord"])})
    report["multielement"] = [compare_multielement(output,f"two_{n}",geometry(n)) for n in
                              ((41,81) if args.quick else (81,161,321,641))]
    sensitivity_points = 81 if args.quick else 321
    report["sensitivity"] = [compare_multielement(output,label,geometry(sensitivity_points,position=pos,deflection=delta))
        for label,pos,delta in (("base",(1.02,-.04),25),("less_streamwise_offset",(.98,-.04),25),
                               ("more_streamwise_offset",(1.10,-.04),25),("lower_flap",(1.02,-.07),25),
                               ("deflection_15",(1.02,-.04),15))]
    finite_counts = (41,81) if args.quick else (41,81,161,321,641)
    for n in finite_counts:
        shape = geometry(n,closed=False)
        original = meshes_of(shape)
        for parts in (1,4,8):
            meshes = tuple(subdivide_cap(m,parts) for m in original)
            system,hs = solve_hs(meshes)
            label = f"finite_{n}_cap{parts}"
            identifier = record_case(output,label,(shape.main_points,shape.flap_points),6.,
                {"generator": "NACA finite TE, cosine stations", "cap_subdivisions": parts},hess_smith=hs)
            report["finite_te"].append({"points_per_element": n,"cap_subdivisions": parts,"input_sha256": identifier,
                "hess_smith": summary_hs(hs),"terminal_length_over_cap": [float(m.lengths[0]/m.lengths[-1]) for m in original],
                "extrapolated_kutta_diagnostic": extrapolated_kutta(system)})
    for n in ((81,) if args.quick else (81,161,321)):
        for power in (.75,1.,1.25):
            shape = finite_geometry(n,power)
            _,hs = solve_hs(meshes_of(shape))
            identifier = record_case(output,f"finite_density_{n}_{int(power*100)}",(shape.main_points,shape.flap_points),6.,
                {"generator": "analytic NACA finite TE", "station_formula": "x=1-(1-x_cosine)^te_power", "te_power": power},
                hess_smith=hs)
            report["finite_te_density"].append({"points_per_element": n,"te_power": power,"input_sha256": identifier,
                                                "hess_smith": summary_hs(hs)})
    # Warm BLAS before timings; samples below are operational, not case identity.
    np.linalg.solve(np.eye(16),np.ones(16))
    for n in ((41,81) if args.quick else (41,81,161,321,641,961)):
        report["benchmarks"].append(benchmark(n,1 if args.quick else 3))
    report["condition_bound_experiment"] = [inverse_bound_experiment(n) for n in ((81,) if args.quick else (161,641))]
    write_json(output/"review.json",report)
    fig,axes = plt.subplots(1,3,figsize=(14,4),layout="constrained")
    for skew,label in ((0.,"HS uniform circle angle"),(-.5,"HS cusp-biased angle")):
        rows = [r for r in report["joukowski"] if r["skew"] == skew]
        axes[0].loglog([r["panels"] for r in rows],[r["cm_error"] for r in rows],"o-",label=label)
    rows = [r for r in report["joukowski"] if r["skew"] == 0 and r["linear_vortex"]["accepted"]]
    axes[0].loglog([r["panels"] for r in rows],[abs(r["linear_vortex"]["cm_ccw"]-exact["cm_ccw_quarter_chord"]) for r in rows],"s-",label="Linear vortex reference")
    axes[0].set(xlabel="Panels",ylabel="Absolute CM error",title="Exact Joukowski")
    rows = report["multielement"]
    axes[1].loglog([r["panels"] for r in rows],[abs(r["cl_difference"]) for r in rows],"o-",label="|CL HS - CL LV|")
    axes[1].loglog([r["panels"] for r in rows],[abs(r["cm_difference"]) for r in rows],"s-",label="|CM HS - CM LV|")
    axes[1].set(xlabel="Total panels",ylabel="Coefficient difference",title="Independent two-element comparison")
    axes[1].set_xticks([r["panels"] for r in rows],labels=[str(r["panels"]) for r in rows])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    for parts in (1,4,8):
        rows = [r for r in report["finite_te"] if r["cap_subdivisions"] == parts]
        axes[2].plot([2*(r["points_per_element"]-1) for r in rows],[r["hess_smith"]["cl"] for r in rows],"o-",label=f"{parts} base panels / element")
    axes[2].set(xlabel="Surface panels (two elements)",ylabel="CL",title="Finite TE: base sensitivity")
    for axis in axes:
        axis.grid(alpha=.25)
        axis.legend(fontsize=8)
    fig.savefig(output/"review.png",dpi=160)
    plt.close(fig)
    inventory = []
    for path in sorted(output.rglob("*")):
        if path.is_file():
            inventory.append({"path": path.relative_to(output).as_posix(),"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                              "size_bytes": path.stat().st_size})
    write_json(output/"artifacts.json",{"artifacts": inventory})
    print(f"Review evidence: {output/'review.json'}")
    print("Finite-TE closure remains unverified; inviscid results do not establish suitability for optimization.")


if __name__ == "__main__":
    main()
