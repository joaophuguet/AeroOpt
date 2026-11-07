"""Opt-in v0.3.2 analytic benchmark and conservative numerical ordering study.

Williams circle data are reproduced. The double-map interpretation is NOT
certified as the printed airfoil; its precisely declared analytic two-body
case is a fallback reference. No fitted coefficients, optimization or automatic
refinement inside PanelSolver. See docs/VALIDATION_V032.md.
"""

import argparse
from collections import Counter
import csv
import hashlib
from itertools import combinations
import os
from pathlib import Path
import platform
import time

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR",str(ROOT/".cache/matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from aeroopt import __version__
from aeroopt.geometry import validate_elements
from aeroopt.solvers import MomentReference, OperatingConditions, PanelSolver, RunStatus, SolverRunner
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers.serialization import read_json, scientific_hash, write_json
from aeroopt.verification import RefinementRecord, RefinementStudy, RankingPolicy, compare_pair
from _linear_vortex_reference import solve_linear_vortex
from _panel_reference import joukowski_exact, joukowski_points
from _panel_review_tools import subdivide_cap
from _pressure_quadrature import pressure_quadrature
from _williams_reference import PrintedParameterCase
from verify_panel_solver import geometry


def load_published(directory):
    catalog = read_json(Path(__file__).with_name("reference_catalog_v032.json"))
    for name,key in ((catalog["file"],"sha256"),(catalog["checks_file"],"checks_sha256")):
        data = (directory/name).read_bytes()
        if hashlib.sha256(data).hexdigest() != catalog[key]:
            raise ValueError(f"pinned Williams source/check hash mismatch: {name}")
    data = read_json(directory/catalog["checks_file"])
    if data["source_sha256"] != catalog["sha256"]:
        raise ValueError("transcription source mismatch")
    return catalog,data


def published_comparison(case,data):
    """Compare rounded values, retaining mismatch instead of tuning a map."""
    delta = np.angle(1.9-.4j)
    te = (np.asarray(case.trailing_points)-case.centers[0])*np.exp(-1j*delta)
    expected = np.asarray(data["stagnation_points_aligned_centers_frame"])
    xy = np.column_stack((te.real,te.imag))
    result = {"circle_stagnation_coordinate_max_error": float(np.max(abs(xy-expected))),
              "circle_circulation_max_error": float(np.max(abs(-case.gamma_raw-data["circle_circulations_clockwise"]))),
              "coordinate_projections": []}
    result["circle_table_reproduced_at_printed_precision"] = (
        result["circle_stagnation_coordinate_max_error"] <= 5e-5 and result["circle_circulation_max_error"] <= 5e-4)
    ref = case.reference(8192)
    result["airfoil_circulation_max_error"] = float(np.max(abs(-np.array(ref["element_circulations_ccw"])-data["airfoil_circulations_clockwise"])))
    result["total_cl_circulation_difference"] = ref["cl_circulation"]-data["total_cl_circulation"]
    result["element_cl_pressure_difference"] = (np.array(ref["element_cl"])-data["element_cl_pressure"]).tolist()
    # Projection is only a comparison of rounded printed coordinates. It is
    # never used to repair, resample or generate solver input geometry.
    for body,key in enumerate(("main_xy_cp","flap_xy_cp")):
        for x,y,cp in data.get(key,[]):
            target = x+1j*y
            c,tepoint = case.centers[body],case.trailing_points[body]

            def point(theta):
                theta = theta%(2*np.pi)
                return case.published_frame(c+(tepoint-c)*np.exp(1j*theta))

            if cp == 1.:
                theta = 0.
                calculated = case.published_frame(tepoint)
            else:
                grid = np.arange(2048)*2*np.pi/2048
                index = int(np.argmin(abs(point(grid)-target)))
                low,high = grid[index]-2*np.pi/2048,grid[index]+2*np.pi/2048
                ratio = (np.sqrt(5)-1)/2
                for step in range(64):
                    left,right = high-ratio*(high-low),low+ratio*(high-low)
                    if abs(point(left)-target) < abs(point(right)-target):
                        high = right
                    else:
                        low = left
                theta = ((low+high)/2)%(2*np.pi)
                calculated = point(theta)
            actual_cp = float(case.surface(body,np.array([theta]))[1][0])
            result["coordinate_projections"].append({"body": body,"published_xy": [x,y],
                "published_cp": cp,"nearest_curve_xy": [float(calculated.real),float(calculated.imag)],
                "coordinate_distance": float(abs(calculated-target)),"calculated_cp": actual_cp,
                "cp_difference": actual_cp-cp,"scope": "projection_diagnostic_not_input_repair"})
    max_distance = max((r["coordinate_distance"] for r in result["coordinate_projections"]),default=None)
    result["airfoil_coordinate_max_distance"] = max_distance
    result["airfoil_reproduction_certified"] = False
    result["limitation"] = "intermediate second-map construction and printed-frame reproduction remain unverified; no fitting"
    return result


def kernel_solution(points,alpha):
    validation = validate_elements(*points,min_clearance=0.)
    if not validation.valid or validation.min_distance < pm.MIN_CLEARANCE:
        raise ValueError("benchmark geometry outside panel numerical domain")
    meshes = tuple(pm.make_mesh(p,(len(p)-1)//2) for p in points)
    start = time.perf_counter()
    hs = pm.solve(pm.assemble(meshes,alpha),.25,0.,"counterclockwise")
    duration = time.perf_counter()-start
    return meshes,hs,duration


def benchmark(output,counts,published):
    rows,references = [],{}
    for name,direction in (("A_like",30.),("B_like",10.)):
        case = PrintedParameterCase(direction)
        start = time.perf_counter()
        coarse,fine = case.reference(8192),case.reference(16384)
        references[name] = {"reference": fine,"provenance": case.provenance(),
                            "analytic_flow_checks": case.flow_checks(),
                            "reference_evaluation_s": time.perf_counter()-start,
                            "quadrature_difference_cl": abs(coarse["cl"]-fine["cl"]),
                            "quadrature_difference_cm": abs(coarse["cm_ccw"]-fine["cm_ccw"])}
        if published is not None:
            references[name]["published_comparison"] = published_comparison(case,published[name[0]])
        for count in counts:
            points = case.contours(count)
            meshes,hs,hs_s = kernel_solution(points,case.alpha_deg)
            start = time.perf_counter()
            try:
                lv = solve_linear_vortex(points,case.alpha_deg)
            except ArithmeticError as error:
                lv = {"accepted": False,"reason": str(error)}
            lv_s = time.perf_counter()-start
            exact_cp = []
            location_errors = []
            for body in (0,1):
                theta = case.parameter_nodes(body,count)
                p,cp,_ = case.surface(body,(theta[:-1]+theta[1:])/2)
                exact_cp.append(cp)
                physical = np.column_stack((p.real,p.imag))
                location_errors.append(float(np.max(np.linalg.norm(physical-meshes[body].centers,axis=1))))
            inputs = {"coordinates": [p.tolist() for p in points],"alpha_deg": case.alpha_deg,
                      "provenance": case.provenance(),"moment_reference": [.25,0.,"CCW"]}
            row = {"case": name,"panels": 2*count,"hs_s": hs_s,"linear_vortex_s": lv_s,
                   "input_sha256": scientific_hash(inputs),"methods": {},"curve_to_panel_center_max_offset": location_errors}
            for method,result in (("hess_smith",hs),("linear_vortex",lv)):
                if not result["accepted"]:
                    row["methods"][method] = {"accepted": False,"reason": result["reason"]}
                    continue
                cl,cm = result["cl"],result["cm"] if method == "hess_smith" else result["cm_ccw"]
                gamma = [e["circulation"] for e in result["elements"]]
                row["methods"][method] = {"accepted": True,"cl": cl,"cm_ccw": cm,
                    "diagnostics": result["diagnostics"] if method == "hess_smith" else {
                        key: result[key] for key in ("condition_number_2","normal_residual","kutta_residual","cl_circulation")},
                    "cl_error": cl-fine["cl_circulation"],"cm_error": cm-fine["cm_blasius"],
                    "circulations": gamma,
                    "circulation_max_error": float(np.max(abs(np.array(gamma)-fine["element_circulations_ccw"]))),
                    "cp_rms_error": [float(np.sqrt(np.mean((np.array(e["cp"])-cp)**2)))
                                     for e,cp in zip(result["elements"],exact_cp,strict=True)]}
            write_json(output/"cases"/f"benchmark_{name}_{count}.json",{"input": inputs,"input_sha256": row["input_sha256"],
                "hess_smith": hs,"linear_vortex": lv,"exact_cp_curve_parameter_midpoints": [p.tolist() for p in exact_cp],
                "location_note": "exact curve parameter midpoints versus flat-panel collocation; geometry error included"})
            rows.append(row)
    return references,rows


def ranking_study(output,counts,quick):
    # Independent one-variable changes and very small changes around a fixed,
    # geometrically valid baseline. No optimization/search of design variables.
    candidates = [("base",25.,1.02,-.04,6.),("deflection_plus",26.,1.02,-.04,6.),
                  ("deflection_minus",24.,1.02,-.04,6.),("x_plus",25.,1.03,-.04,6.),
                  ("x_minus",25.,1.01,-.04,6.),("y_plus",25.,1.02,-.035,6.),
                  ("y_minus",25.,1.02,-.045,6.),("alpha_plus",25.,1.02,-.04,6.5),
                  ("alpha_minus",25.,1.02,-.04,5.5),("near_deflection",25.001,1.02,-.04,6.),
                  ("near_x",25.,1.02001,-.04,6.),("near_y",25.,1.02,-.03999,6.),
                  ("near_alpha",25.,1.02,-.04,6.0001)]
    if quick:
        candidates = [candidates[i] for i in (0,1,9,10,12)]
    runner = SolverRunner(ROOT,(output/"attempts").relative_to(ROOT))
    solver = PanelSolver(MomentReference(.25,0.,"counterclockwise","q_inf*c_main^2"))
    studies,metadata = [],[]
    for name,delta,x,y,alpha in candidates:
        records,attempts,times = [],[],[]
        start_study = time.perf_counter()
        for count in counts:
            shape = geometry(count+1,position=(x,y),deflection=delta)
            start = time.perf_counter()
            result = runner.run(solver,shape,OperatingConditions(alpha,200000,0),min_clearance=.002)
            hs_s = time.perf_counter()-start
            accepted = result.status == RunStatus.SUCCESS
            records.append(RefinementRecord("hess_smith",2*count,result.coefficients.cl,result.coefficients.cm,accepted))
            start = time.perf_counter()
            try:
                lv = solve_linear_vortex((shape.main_points,shape.flap_points),alpha)
            except ArithmeticError as error:
                lv = {"accepted": False,"reason": str(error)}
            lv_s = time.perf_counter()-start
            records.append(RefinementRecord("linear_vortex",2*count,lv.get("cl"),lv.get("cm_ccw"),lv["accepted"]))
            write_json(output/"cases"/f"ranking_{name}_{count}.json",{"case_id": result.case_id,
                "attempt_directory": result.attempt_directory,"linear_vortex": lv,
                "geometry_numeric_sha256": scientific_hash([shape.main_points.tolist(),shape.flap_points.tolist()])})
            attempts.append(result.attempt_directory)
            times.append({"panels": 2*count,"runner_hs_s": hs_s,"linear_vortex_s": lv_s})
        study = RefinementStudy(name,tuple(records))
        studies.append(study)
        metadata.append({"study": study.to_dict(),"configuration": {"deflection_deg": delta,"flap_le": [x,y],"alpha_deg": alpha},
                         "attempts": attempts,"timings": times,"total_history_s": time.perf_counter()-start_study,
                         "cl_diagnostics": study.analyze("cl"),"cm_diagnostics": study.analyze("cm")})
    policy = RankingPolicy(safety_factor=2.,comparison_floor=1e-12)
    pairs = [compare_pair(a,b,coefficient,policy=policy) for a,b in combinations(studies,2) for coefficient in ("cl","cm")]
    counts_result = {c:dict(Counter(p["classification"] for p in pairs if p["coefficient"] == c)) for c in ("cl","cm")}
    return {"studies": metadata,"pairs": pairs,"counts": counts_result,
            "scope": "local numerical ordering; heuristic indicators, no universal or physical error bounds"}


def remaining_limitations(output,quick):
    joukowski = []
    exact = joukowski_exact(4)
    for count in ((80,160) if quick else (160,640,1920)):
        meshes = (pm.make_mesh(joukowski_points(count),count//2),)
        solved = pm.solve(pm.assemble(meshes,4),.25,0.,"counterclockwise")
        variants = [pressure_quadrature(meshes,solved,4,order) for order in (2,4)]
        joukowski.append({"panels": count,"midpoint_cm": solved.get("cm"),"exact_cm": exact["cm_ccw_quarter_chord"],
                          "pressure_quadrature_diagnostics": variants,"production_improvement_adopted": False})
    finite = []
    for count in ((40,80) if quick else (80,160,320,640)):
        shape = geometry(count+1,closed=False)
        points = (shape.main_points,shape.flap_points)
        original = tuple(pm.make_mesh(p,count//2) for p in points)
        for parts in (1,4,8,16):
            meshes = tuple(subdivide_cap(m,parts) for m in original)
            solved = pm.solve(pm.assemble(meshes,6),.25,0.,"counterclockwise")
            finite.append({"surface_panels_per_element": count,"cap_parts": parts,"cl": solved.get("cl"),
                           "cm_ccw": solved.get("cm"),"diagnostics": solved["diagnostics"],
                           "limitation": "FINITE_TE_UNVERIFIED"})
    return {"joukowski": joukowski,"finite_te": finite,
            "note": "finite wedge angle of a sharp analytic TE does not validate a finite-thickness cap"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"outputs/verification_v032")
    parser.add_argument("--williams-data",type=Path,help="Local pinned PDF plus manually checked excerpts; no downloads")
    parser.add_argument("--quick",action="store_true",help="Smoke study only; not full scientific evidence")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT) or output.exists() and any(output.iterdir()):
        parser.error("Output must be a new/empty directory inside the workspace")
    (output/"cases").mkdir(parents=True,exist_ok=True)
    start = time.perf_counter()
    published = None
    if args.williams_data is not None:
        catalog,published = load_published(args.williams_data)
        (output/"originals").mkdir()
        for name in (catalog["file"],catalog["checks_file"]):
            with (output/"originals"/name).open("xb") as stream:
                stream.write((args.williams_data/name).read_bytes())
    counts = (40,80,160) if args.quick else (80,160,320,640)
    references,benchmark_rows = benchmark(output,counts,published)
    ranking = ranking_study(output,counts,args.quick)
    limitations = remaining_limitations(output,args.quick)
    report = {"aeroopt_version": __version__,"production_adapter_version": PanelSolver(MomentReference(.25,0.,"counterclockwise","q_inf*c_main^2")).identity.version,
              "python": platform.python_version(),"numpy": np.__version__,"quick": args.quick,
              "references": references,"benchmark": benchmark_rows,"ranking": ranking,"remaining_limitations": limitations,
              "total_study_s": time.perf_counter()-start,
              "acceptance": {"williams_published_airfoil_reconstruction_certified": False,
                             "independent_published_benchmark_acceptance_satisfied": False,
                             "analytic_two_body_fallback_available": True,"physical_validation": False,
                             "preliminary_automatic_optimization_ready": False}}
    write_json(output/"verification.json",report)
    with (output/"benchmark_errors.csv").open("x",encoding="utf-8",newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case","panels","method","cl_error","cm_error","max_circulation_error","cp_rms_main","cp_rms_flap"])
        for row in benchmark_rows:
            for method,data in row["methods"].items():
                if data["accepted"]:
                    writer.writerow([row["case"],row["panels"],method,data["cl_error"],data["cm_error"],
                                     data["circulation_max_error"],*data["cp_rms_error"]])
    with (output/"ranking_pairs.csv").open("x",encoding="utf-8",newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["a","b","coefficient","classification","relation","margin_indicator"])
        for row in ranking["pairs"]:
            writer.writerow([row[k] for k in ("a","b","coefficient","classification","relation","margin_indicator")])
    fig,axes = plt.subplots(1,3,figsize=(14,4),layout="constrained")
    for method in ("hess_smith","linear_vortex"):
        rows = [r for r in benchmark_rows if r["case"] == "A_like"]
        accepted = [r for r in rows if r["methods"][method]["accepted"]]
        axes[0].loglog([r["panels"] for r in accepted],[abs(r["methods"][method]["cl_error"]) for r in accepted],"o-",label=method)
        axes[1].loglog([r["panels"] for r in accepted],[abs(r["methods"][method]["cm_error"]) for r in accepted],"o-",label=method)
    axes[0].set(xlabel="Total panels",ylabel="Absolute CL error",title="Defined analytic two-body case")
    axes[1].set(xlabel="Total panels",ylabel="Absolute CM error",title="CCW quarter-chord moment")
    for axis in axes[:2]:
        axis.set_xticks([2*n for n in counts],labels=[str(2*n) for n in counts])
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.legend(fontsize=8)
        axis.grid(alpha=.25)
    categories = ("numerically_distinguishable","numerically_unresolved","inconsistent_across_refinement_or_methods")
    x = np.arange(3)
    axes[2].bar(x-.18,[ranking["counts"]["cl"].get(c,0) for c in categories],.36,label="CL")
    axes[2].bar(x+.18,[ranking["counts"]["cm"].get(c,0) for c in categories],.36,label="CM")
    axes[2].set_xticks(x,labels=["Distinguishable","Unresolved","Inconsistent"],rotation=15)
    axes[2].set(ylabel="Pairwise comparisons",title="Local numerical ordering indicators")
    axes[2].legend()
    fig.savefig(output/"verification.png",dpi=160)
    plt.close(fig)
    source_paths = (Path(__file__),Path(__file__).with_name("_williams_reference.py"),
                    Path(__file__).with_name("_linear_vortex_reference.py"),Path(__file__).with_name("_pressure_quadrature.py"),
                    Path(__file__).with_name("_panel_review_tools.py"),Path(__file__).with_name("verify_panel_solver.py"),
                    Path(pm.__file__),ROOT/"src/aeroopt/verification/refinement.py")
    (output/"sources").mkdir()
    for path in source_paths:
        with (output/"sources"/path.name).open("xb") as stream:
            stream.write(path.read_bytes())
    inventory = [{"path": p.relative_to(output).as_posix(),"sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                  "size_bytes": p.stat().st_size} for p in sorted(output.rglob("*")) if p.is_file()]
    write_json(output/"artifacts.json",{"artifacts": inventory})
    print(f"Evidence: {output/'verification.json'}")
    print("Published Williams airfoil reconstruction NOT certified; exact-flow fallback and local ordering diagnostics only.")


if __name__ == "__main__":
    main()
