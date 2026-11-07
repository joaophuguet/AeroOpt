"""Auditable v0.4.0 integration cases and costs; no download or physical validation."""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np

from aeroopt import __version__
from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points
from aeroopt.solvers import (AlulaSolver,AlulaViscousSettings,AttemptContext,MomentReference,
                            OperatingConditions,SolverRunner)
from aeroopt.solvers.alula import LOCK,RESOURCES,UPSTREAM_COMMIT,node_path
from aeroopt.solvers.serialization import read_json,scientific_hash,write_json

ROOT = Path(__file__).resolve().parents[1]
REF = MomentReference(.25,0,"counterclockwise","q_inf*c_main^2")


def assembly(panels=80,*,closed=True,flap_code="4412",position=(1.02,-.04),deflection=25):
    def foil(code):
        return Airfoil.from_points(naca4_points(code,panels+1,closed_te=closed),leading_edge_index=panels//2)
    return MultiElementAirfoil(foil("0012"),foil(flap_code),FlapParameters(.3,position,deflection))


def single_payload(solver,mode,alpha,*,closed=True,max_iterations=60):
    points = naca4_points("0012",81,closed_te=closed)
    return {"schema_version":"aeroopt.alula-input.v1","case_id":scientific_hash({"points":points.tolist(),"mode":mode,"alpha":alpha,"max_iterations":max_iterations}),
            "mode":mode,"elements":[{"name":"NACA0012 analytic","role":"main",
                "points":[{"x":float(x),"y":float(y)} for x,y in points],
                "closure":"sharp_repeated_endpoint" if closed else "straight_TE_cap"}],
            "conditions":{"alpha_deg":alpha,"reynolds":1e6,"mach":0},"moment_reference":REF.to_dict(),
            "settings":{"boundary_condition":"normal-velocity" if mode=="inviscid" and closed else "streamfunction",
                        "viscous":None if mode=="inviscid" else asdict(AlulaViscousSettings(9,(1,1),max_iterations=max_iterations))},
            "runtime":{"node_version":solver.identity.settings.to_dict()["node_version"],"upstream_commit":UPSTREAM_COMMIT}}


def standalone(solver,output,name,payload):
    payload=dict(payload)
    payload["case_id"]=scientific_hash({"input":{k:v for k,v in payload.items() if k!="case_id"},
        "bridge_sha256":solver.identity.settings.to_dict()["bridge_sha256"]})
    directory=output/name;directory.mkdir()
    context=AttemptContext(directory,60)
    context.write_json("input.json",payload)
    evidence=context.run_process([str(solver.node_executable),str(RESOURCES/"bridge.mjs"),
        "--upstream",str(solver.alula_directory),"--input",str(directory/"input.json"),
        "--output",str(directory/"output.json")],label="bridge")
    context.write_json("execution.json",evidence.to_dict())
    if evidence.return_code != 0 or evidence.timed_out:
        raise RuntimeError(f"Standalone bridge failed: {directory}")
    result=read_json(directory/"output.json")
    baseline=context.run_process([str(solver.node_executable),str(ROOT/"examples/alula_upstream_baseline.mjs"),
        str(solver.alula_directory),str(directory/"input.json"),str(directory/"upstream_baseline.json")],label="upstream")
    context.write_json("upstream_execution.json",baseline.to_dict())
    if baseline.return_code != 0 or baseline.timed_out:
        raise RuntimeError(f"Upstream baseline failed: {directory}")
    raw=read_json(directory/"upstream_baseline.json")
    differences={}
    if result["status"]["eligible"]:
        for key in ("cl","cd","cm"):
            actual,expected=result["coefficients"][key],raw.get(key)
            if key=="cm" and expected is not None: expected=-expected
            differences[key]=None if actual is None or expected is None else actual-expected
        if any(v is not None and abs(v)>1e-11 for v in differences.values()):
            raise ArithmeticError("Bridge/upstream numerical baseline mismatch")
    return {"name":name,"path":directory.relative_to(output).as_posix(),"status":result["status"],
            "coefficients":result["coefficients"],"baseline_differences":differences,
            "node_process_s":evidence.duration_s,"solver_s":result["metadata"]["timings"]["solver_ms"]/1000,
            "rss_bytes":result["metadata"]["memory"]["rss"],"warnings":result["warnings"],"limitations":result["limitations"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"outputs/alula_v040")
    parser.add_argument("--node",type=Path)
    parser.add_argument("--alula",type=Path)
    parser.add_argument("--repeats",type=int,default=3)
    args=parser.parse_args()
    output=args.output.resolve()
    if not output.is_relative_to(ROOT) or output.exists() and any(output.iterdir()) or not 1<=args.repeats<=10:
        parser.error("Use a new/empty workspace output and 1-10 repetitions")
    output.mkdir(parents=True,exist_ok=True)
    solver=AlulaSolver("inviscid",REF,node_executable=args.node,alula_directory=args.alula)
    context=AttemptContext(output,30)
    starts=[]
    for k in range(5):
        e=context.run_process([str(solver.node_executable),"--version"],label=f"startup_{k}")
        if e.return_code != 0 or e.timed_out: raise RuntimeError("Node startup probe failed")
        starts.append(e.duration_s)
    rows=[]
    rows.append(standalone(solver,output,"A_symmetric_zero",single_payload(solver,"inviscid",0)))
    for k in range(args.repeats):
        rows.append(standalone(solver,output,f"B_positive_{k}",single_payload(solver,"inviscid",4)))
        rows.append(standalone(solver,output,f"D_viscous_{k}",single_payload(solver,"viscous_incompressible",0,closed=False)))
    p=single_payload(solver,"viscous_incompressible",2,closed=False,max_iterations=1)
    p["settings"]["viscous"]["max_wake_iterations"]=0
    rows.append(standalone(solver,output,"E_nonconvergent",p))
    runner=SolverRunner(ROOT,(output/"attempts").relative_to(ROOT))
    assembly_rows=[]
    shape=assembly()
    for k in range(args.repeats):
        start=time.perf_counter()
        result=runner.run(solver,shape,OperatingConditions(6,1e6,0),min_clearance=.002,timeout_s=60)
        wall=time.perf_counter()-start
        path=runner.execution_root/result.attempt_directory
        raw=read_json(path/"alula_output.json")
        # Baseline logs and outputs stay outside the sealed runner attempt.
        baseline_path=output/f"C_upstream_{k}";baseline_path.mkdir()
        baseline=AttemptContext(baseline_path,60).run_process([str(solver.node_executable),str(ROOT/"examples/alula_upstream_baseline.mjs"),
            str(solver.alula_directory),node_path(path/"alula_input.json"),str(baseline_path/"upstream_baseline.json")],label="baseline")
        write_json(baseline_path/"execution.json",baseline.to_dict())
        if baseline.return_code != 0: raise RuntimeError("Assembly baseline failed")
        direct=read_json(baseline_path/"upstream_baseline.json")
        assembly_rows.append({"name":f"C_slotted_{k}","status":result.status.value,"case_id":result.case_id,
            "attempt_directory":result.attempt_directory,"coefficients":result.coefficients.to_dict(),"runner_s":wall,
            "node_process_s":read_json(path/"execution.json")["duration_s"],
            "solver_s":raw["metadata"]["timings"]["solver_ms"]/1000,"rss_bytes":raw["metadata"]["memory"]["rss"],
            "baseline_difference_cl":result.coefficients.cl-direct["cl"],
            "baseline_difference_cm":result.coefficients.cm+direct["cm"]})
    viscous=AlulaSolver("viscous_incompressible",REF,node_executable=solver.node_executable,
                       alula_directory=solver.alula_directory,viscous_settings=AlulaViscousSettings(9,(1,1)))
    start=time.perf_counter()
    viscous_result=runner.run(viscous,assembly(40,flap_code="0012",position=(1.2,-.15),deflection=0),
        OperatingConditions(0,1e6,0),min_clearance=.002,timeout_s=60)
    viscous_wall=time.perf_counter()-start
    finite=AlulaSolver("inviscid",REF,node_executable=solver.node_executable,alula_directory=solver.alula_directory,boundary_condition="streamfunction")
    finite_result=runner.run(finite,assembly(closed=False),OperatingConditions(6,1e6,0),min_clearance=.002,timeout_s=60)
    unsupported=runner.run(solver,shape,OperatingConditions(6,1e6,.1),min_clearance=.002,timeout_s=30)
    consistency={}
    for prefix in ("B_positive","D_viscous"):
        group=[r for r in rows if r["name"].startswith(prefix)]
        consistency[prefix]={key:None if group[0]["coefficients"][key] is None else
            max(r["coefficients"][key] for r in group)-min(r["coefficients"][key] for r in group) for key in ("cl","cd","cm")}
    consistency["C_slotted"]={key:max(r["coefficients"][key] for r in assembly_rows)-min(r["coefficients"][key] for r in assembly_rows) for key in ("cl","cm")}
    report={"aeroopt_version":__version__,"python":platform.python_version(),"platform":platform.platform(),
            "upstream_commit":UPSTREAM_COMMIT,"identity":solver.identity.to_dict(),"startup_supervised_s":starts,
            "startup_median_s":statistics.median(starts),"standalone_cases":rows,"assembly_cases":assembly_rows,
            "viscous_assembly":{"result":viscous_result.to_dict(),"runner_s":viscous_wall},
            "finite_base_assembly":finite_result.to_dict(),"unsupported":unsupported.to_dict(),"repeat_spreads":consistency,
            "physical_validation":False,"optimization_ready":False,
            "cost_scope":"Node process includes supervisor/startup/imports/integrity/JSON; runner includes Python audit/IO. RSS is a sample, not peak."}
    write_json(output/"verification.json",report)
    sources=output/"sources";sources.mkdir()
    for p in (Path(__file__),ROOT/"examples/alula_upstream_baseline.mjs",ROOT/"src/aeroopt/solvers/alula.py",
              RESOURCES/"bridge.mjs",RESOURCES/"upstream-lock.json",RESOURCES/"COPYING"):
        (sources/p.name).write_bytes(p.read_bytes())
    inventory=[{"path":p.relative_to(output).as_posix(),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"bytes":p.stat().st_size}
               for p in sorted(output.rglob("*")) if p.is_file()]
    write_json(output/"artifacts.json",{"artifacts":inventory})
    print(f"Integration evidence: {output/'verification.json'}")
    print("No physical validation or aerodynamic optimization qualification.")


if __name__=="__main__":
    main()
