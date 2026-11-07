"""Pinned real Node integration; dependency-dependent cases skip explicitly."""

from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil, naca4_points
from aeroopt.solvers import (AlulaSolver, AlulaViscousSettings, AttemptContext, ConvergenceStatus,
                            MomentReference, OperatingConditions, RunStatus, SolverCase, SolverRunner)
from aeroopt.solvers.alula import LOCK, RESOURCES, UPSTREAM_COMMIT, node_path, verify_installation
from aeroopt.solvers.serialization import read_json, scientific_hash
from solver_helpers import process_running

ROOT = Path(__file__).resolve().parents[1]
REF = MomentReference(.25,0,"counterclockwise","q_inf*c_main^2")
CONDITIONS = OperatingConditions(6,1e6,0)


def assembly(panels=80,*,closed=True,reverse=False,position=(1.02,-.04),deflection=25,flap_code="4412"):
    def foil(code):
        p = naca4_points(code,panels+1,closed_te=closed)
        return Airfoil.from_points(p[::-1] if reverse else p,leading_edge_index=panels//2)
    return MultiElementAirfoil(foil("0012"),foil(flap_code),FlapParameters(.3,position,deflection))


@pytest.fixture(scope="module")
def runtime():
    node = shutil.which("node")
    local = ROOT/".deps"/f"node-{LOCK['tested_node']['version']}-win-x64/node.exe"
    if node is None and local.is_file():
        node = str(local)
    source = ROOT/".deps/alula"/UPSTREAM_COMMIT
    if node is None or not source.is_dir():
        pytest.skip("Optional Node >=22 and pinned Alula installation required; run scripts/install_alula.py explicitly")
    verify_installation(source)  # A changed installation is an error, not a skip.
    return Path(node),source


@pytest.fixture(scope="module")
def solver(runtime):
    node,source = runtime
    return AlulaSolver("inviscid",REF,node_executable=node,alula_directory=source)


def bridge_input(solver,*,mode="inviscid",alpha=0,closed=True):
    points = naca4_points("0012",81,closed_te=closed)
    return {"schema_version":"aeroopt.alula-input.v1","case_id":"a"*64,"mode":mode,
            "elements":[{"name":"analytic NACA0012","role":"main","points":[{"x":float(x),"y":float(y)} for x,y in points],
                         "closure":"sharp_repeated_endpoint" if closed else "straight_TE_cap"}],
            "conditions":{"alpha_deg":alpha,"reynolds":1e6,"mach":0},"moment_reference":REF.to_dict(),
            "runtime":{"node_version":solver.identity.settings.to_dict()["node_version"],"upstream_commit":UPSTREAM_COMMIT},
            "settings":{"boundary_condition":"normal-velocity" if mode=="inviscid" and closed else "streamfunction",
                        "viscous":None if mode=="inviscid" else as_settings()}}


def as_settings(**changes):
    from dataclasses import asdict
    return json.loads(json.dumps(asdict(AlulaViscousSettings(9,(1,1),**changes))))


def run_bridge(solver,path,payload,*,raw_text=None):
    path.mkdir(parents=True,exist_ok=False)
    (path/"input.json").write_text(raw_text if raw_text is not None else json.dumps(payload),encoding="utf-8")
    context = AttemptContext(path,30)
    evidence = context.run_process([str(solver.node_executable),str(RESOURCES/"bridge.mjs"),"--upstream",str(solver.alula_directory),
        "--input",node_path(path/"input.json"),"--output",node_path(path/"output.json")],label="bridge")
    return evidence,read_json(path/"output.json") if (path/"output.json").exists() else None


@pytest.mark.parametrize("kwargs",[{"ncrit":15},{"trips":(0,1)},{"trips":(1,)},{"wake_count":3},
                                    {"max_iterations":0},{"wake_tolerance":0},{"max_wake_iterations":-1},
                                    {"tolerance":float("nan")},{"wake_relaxation":2}])
def test_viscous_settings_reject_invalid_or_implicit_physics(kwargs):
    values={"ncrit":9,"trips":(1,1),**kwargs}
    with pytest.raises((ValueError,TypeError)):
        AlulaViscousSettings(**values)


def test_missing_alula_has_actionable_diagnostic(tmp_path):
    with pytest.raises(FileNotFoundError,match="Missing pinned Alula"):
        AlulaSolver("inviscid",REF,alula_directory=tmp_path/"absent")


def test_missing_node_has_actionable_diagnostic(tmp_path,monkeypatch):
    monkeypatch.setattr("aeroopt.solvers.alula.verify_installation",lambda p:tmp_path)
    monkeypatch.setattr("aeroopt.solvers.alula.shutil.which",lambda p:None)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError,match="Node.js missing"):
        AlulaSolver("inviscid",REF)


@pytest.mark.parametrize("mode,settings",[("euler",None),("viscous_incompressible",None),
                                         ("inviscid",AlulaViscousSettings(9,(1,1)))])
def test_unsupported_modes_or_incompatible_settings_rejected(mode,settings):
    with pytest.raises(ValueError):
        AlulaSolver(mode,REF,viscous_settings=settings)


def test_identity_stable_across_installation_paths_and_attempts(solver,tmp_path):
    copied = tmp_path/"Alula cópia com espaços"
    shutil.copytree(solver.alula_directory,copied)
    other = AlulaSolver("inviscid",REF,node_executable=solver.node_executable,alula_directory=copied)
    assert solver.identity == other.identity
    assert str(copied) not in json.dumps(other.identity.to_dict())
    shape = assembly()
    assert SolverCase(shape,CONDITIONS,solver.identity,.002).case_id == SolverCase(shape,CONDITIONS,other.identity,.002).case_id
    changed = AlulaSolver("inviscid",replace(REF,x_over_c=.3),node_executable=solver.node_executable,alula_directory=copied)
    assert solver.identity != changed.identity
    assert solver.identity.settings.to_dict()["upstream_commit"] == UPSTREAM_COMMIT


def test_source_tampering_rejected(solver,tmp_path):
    copied = tmp_path/"source"
    shutil.copytree(solver.alula_directory,copied)
    with (copied/"src/inviscid/linear-vortex.js").open("ab") as stream:
        stream.write(b"\n// changed\n")
    with pytest.raises(ValueError,match="hash mismatch"):
        verify_installation(copied)


@pytest.mark.parametrize("alpha,expected_cl,expected_cm",[
    (0,0.,0.),(4,.4767885298826391,.006445853438382634)])
def test_single_element_against_direct_pinned_upstream_baseline(solver,tmp_path,alpha,expected_cl,expected_cm):
    # Values obtained by direct solveInviscid(naca4('0012',80), alpha), recorded
    # in VALIDATION_V040. The baseline is integration parity, not physical truth.
    evidence,output = run_bridge(solver,tmp_path/"single",bridge_input(solver,alpha=alpha))
    assert evidence.return_code == 0
    assert output["status"]["eligible"]
    assert output["coefficients"]["cl"] == pytest.approx(expected_cl,abs=2e-11)
    assert output["coefficients"]["cm"] == pytest.approx(expected_cm,abs=2e-11)
    assert output["coefficients"]["cd"] is None


def test_reliable_single_viscous_case_against_upstream_baseline(solver,tmp_path):
    evidence,output = run_bridge(solver,tmp_path/"viscous",bridge_input(solver,mode="viscous_incompressible",closed=False))
    assert evidence.return_code == 0
    assert output["status"]["eligible"]
    assert abs(output["coefficients"]["cl"]) < 2e-10
    assert abs(output["coefficients"]["cm"]) < 2e-10
    assert output["coefficients"]["cd"] == pytest.approx(.005344683098715396,abs=2e-10)
    assert output["upstream"]["diagnostics"]["equationResidual"] <= 1e-8
    assert "ALULA_FINITE_BASE_RESEARCH" in output["limitations"]
    assert output["upstream"]["boundaryLayer"]["surfaces"]


def test_nonconvergence_withholds_provisional_values(solver,tmp_path):
    payload = bridge_input(solver,mode="viscous_incompressible",alpha=2,closed=False)
    payload["settings"]["viscous"]["max_iterations"] = 1
    payload["settings"]["viscous"]["max_wake_iterations"] = 0
    evidence,output = run_bridge(solver,tmp_path/"unconverged",payload)
    assert evidence.return_code == 0  # Process completion is explicitly separate.
    assert not output["status"]["eligible"]
    assert output["status"]["coefficient_status"] in ("provisional","missing")
    assert all(v is None for v in output["coefficients"].values())
    assert output["upstream"]["status"] == "unconverged"
    assert output["nonfinite_diagnostics"]  # Infinite wake residual is labeled, never silently zeroed.


def test_finite_base_inviscid_retains_original_points_and_explicit_cap(solver,tmp_path):
    payload=bridge_input(solver,closed=False,alpha=4)
    evidence,output=run_bridge(solver,tmp_path/"finite",payload)
    assert evidence.return_code == 0
    assert output["status"]["eligible"]
    supplied=payload["elements"][0]["points"]
    assert output["geometry"][0]["points"] == supplied
    assert output["upstream"]["elements"][0]["points"] == supplied+[supplied[0]]
    assert "ALULA_FINITE_BASE_RESEARCH" in output["limitations"]
    assert output["coefficients"]["cd"] is None


def test_real_viscous_multielement_adapter_converges_without_qualifying_physics(solver,tmp_path):
    viscous=AlulaSolver("viscous_incompressible",REF,node_executable=solver.node_executable,
        alula_directory=solver.alula_directory,viscous_settings=AlulaViscousSettings(9,(1,1)))
    runner=SolverRunner(tmp_path,"attempts")
    shape=assembly(40,position=(1.2,-.15),deflection=0,flap_code="0012")
    result=runner.run(viscous,shape,OperatingConditions(0,1e6,0),min_clearance=.002,timeout_s=30)
    assert result.status == RunStatus.SUCCESS,result.diagnostics
    assert result.coefficients.cd is not None and result.coefficients.cd > 0
    output=read_json(runner.execution_root/result.attempt_directory/"alula_output.json")
    assert output["metadata"]["physical_validation"] is False
    assert output["metadata"]["aerodynamic_ranking_qualified"] is False
    assert len(output["upstream"]["boundaryLayer"]["wakes"]) == 2


@pytest.mark.parametrize("kind",["reynolds","alpha","few_panels","additional"])
def test_viscous_operating_limits_rejected(solver,kind):
    from aeroopt.solvers import ConditionParameter
    viscous=AlulaSolver("viscous_incompressible",REF,node_executable=solver.node_executable,
        alula_directory=solver.alula_directory,viscous_settings=AlulaViscousSettings(9,(1,1)))
    shape=assembly(40)
    c=OperatingConditions(0,1e6,0)
    if kind=="reynolds": c=OperatingConditions(0,49999,0)
    elif kind=="alpha": c=OperatingConditions(21,1e6,0)
    elif kind=="few_panels": shape=assembly(20)
    else: c=OperatingConditions(0,1e6,0,(ConditionParameter("unrecognized",1,"dimensionless"),))
    assert viscous.validate_input(shape,c)


@pytest.mark.parametrize("mutation",["mach","clockwise","nonfinite","too_few_panels","unknown_mode","invalid_json"])
def test_node_rejects_invalid_json_geometry_and_conditions(solver,tmp_path,mutation):
    p = bridge_input(solver,mode="viscous_incompressible")
    text = None
    if mutation == "mach": p["conditions"]["mach"] = .1
    elif mutation == "clockwise": p["elements"][0]["points"].reverse()
    elif mutation == "nonfinite": text = json.dumps(p).replace('"alpha_deg": 0','"alpha_deg": 1e999')
    elif mutation == "too_few_panels": p["elements"][0]["points"] = p["elements"][0]["points"][:10]
    elif mutation == "unknown_mode": p["mode"] = "euler"
    else: text = "{"
    evidence,output = run_bridge(solver,tmp_path/"bad",p,raw_text=text)
    assert evidence.return_code != 0
    assert output is None


@pytest.fixture(scope="module")
def real_run(solver,tmp_path_factory):
    path = tmp_path_factory.mktemp("alula_reference")
    shape = assembly()
    runner = SolverRunner(path,"attempts")
    result = runner.run(solver,shape,CONDITIONS,min_clearance=.002,timeout_s=30)
    assert result.status == RunStatus.SUCCESS,result.diagnostics
    return SolverCase(shape,CONDITIONS,solver.identity,.002),runner.execution_root/result.attempt_directory,result


def test_python_node_multielement_roundtrip_and_independent_force_signs(solver,real_run):
    case,path,result = real_run
    output = read_json(path/"alula_output.json")
    np.testing.assert_array_equal([[p['x'],p['y']] for p in output['geometry'][0]['points']],case.main_points)
    np.testing.assert_array_equal([[p['x'],p['y']] for p in output['geometry'][1]['points']],case.flap_points)
    force = np.zeros(2);moment = 0.
    for e in output["upstream"]["elements"]:
        p = np.array([[q['x'],q['y']] for q in e['points']])
        cp = np.array([q['cp'] for q in e['cp']])
        delta = np.diff(p,axis=0)
        df = cp[:,None]*np.column_stack((-delta[:,1],delta[:,0]))
        force += df.sum(axis=0)
        arm = (p[:-1]+p[1:])/2-[.25,0]
        moment += np.sum(arm[:,0]*df[:,1]-arm[:,1]*df[:,0])
    angle = np.deg2rad(CONDITIONS.alpha_deg)
    assert result.coefficients.cl == pytest.approx(force@[-np.sin(angle),np.cos(angle)],abs=2e-12)
    assert result.coefficients.cm == pytest.approx(moment,abs=2e-12)
    assert result.coefficients.cl == pytest.approx(3.1199434441091887,abs=2e-11)
    assert result.coefficients.cm == pytest.approx(.8066128645779502,abs=2e-11)
    assert result.coefficients.cd is None
    assert read_json(path/"execution.json")["tree_cleanup"] == "windows_job_terminated_and_empty" if sys.platform == "win32" else True


@pytest.mark.parametrize("kind",["invalid_json","nan","missing","provisional","diagnostic","geometry","case_id","input_mutated","limits"])
def test_parse_rejects_malformed_or_ineligible_output(solver,real_run,tmp_path,kind):
    case,path,_ = real_run
    context = AttemptContext(tmp_path,None)
    for name in ("alula_input.json","alula_output.json"):
        (tmp_path/name).write_bytes((path/name).read_bytes())
    output = read_json(tmp_path/"alula_output.json")
    if kind == "invalid_json": text = "{"
    else:
        if kind == "nan": output["coefficients"]["cl"] = float("nan")
        elif kind == "missing": output["coefficients"]["cm"] = None
        elif kind == "provisional": output["upstream"]["coefficientStatus"] = "unconverged"
        elif kind == "diagnostic": del output["upstream"]["diagnostics"]["linearResidual"]
        elif kind == "geometry": output["geometry"][0]["points"][1]["x"] += .01
        elif kind == "case_id": output["case_id"] = "b"*64
        elif kind == "limits": output["limitations"] = []
        elif kind == "input_mutated":
            changed = read_json(tmp_path/"alula_input.json")
            changed["conditions"]["alpha_deg"] = 7
            (tmp_path/"alula_input.json").write_text(json.dumps(changed))
        text = json.dumps(output)
    (tmp_path/"alula_output.json").write_text(text)
    with pytest.raises((ValueError,KeyError)):
        solver.parse(case,context,None)


def test_withheld_coefficients_are_never_zero(solver,real_run,tmp_path):
    case,path,_ = real_run
    for name in ("alula_input.json","alula_output.json"):
        (tmp_path/name).write_bytes((path/name).read_bytes())
    output = read_json(tmp_path/"alula_output.json")
    output["coefficients"] = dict.fromkeys(("cl","cd","cm"))
    output["status"]["eligible"] = False
    output["status"]["numerical_convergence"] = "not_converged"
    (tmp_path/"alula_output.json").write_text(json.dumps(output))
    parsed = solver.parse(case,AttemptContext(tmp_path,None),None)
    assert parsed.convergence == ConvergenceStatus.NOT_CONVERGED
    assert parsed.coefficients.empty


@pytest.mark.parametrize("mutation",["provisional","missing","missing_diagnostic","nonfinite","outside_model"])
def test_node_translation_never_qualifies_missing_provisional_or_outside_model(solver,real_run,tmp_path,mutation):
    _,path,_=real_run
    output=tmp_path/"translated.json"
    evidence=AttemptContext(tmp_path,15).run_process([str(solver.node_executable),str(ROOT/"tests/alula_translation_fixture.mjs"),
        str(RESOURCES/"bridge.mjs"),node_path(path/"alula_input.json"),node_path(path/"alula_output.json"),mutation,str(output)],label="translation")
    if mutation=="nonfinite":
        assert evidence.return_code != 0
        assert not output.exists()
    else:
        assert evidence.return_code == 0
        result=read_json(output)
        assert result["status"]["eligible"] is False
        assert all(v is None for v in result["coefficients"].values())
        if mutation=="outside_model":
            assert result["status"]["numerical_convergence"] == "converged"
            assert result["status"]["model_admissible"] is False


@pytest.mark.parametrize("kind",["mach","reverse","few","many","finite_needs_streamfunction"])
def test_adapter_rejects_unsupported_inputs_before_execution(solver,tmp_path,kind):
    shape,conditions = assembly(),CONDITIONS
    if kind == "mach": conditions = OperatingConditions(6,1e6,.1)
    elif kind == "reverse": shape = assembly(reverse=True)
    elif kind == "few":
        shape = assembly(panels=6)
    elif kind == "many": shape = assembly(panels=400)
    else: shape = assembly(closed=False)
    assert solver.validate_input(shape,conditions)
    r = SolverRunner(tmp_path,"attempts").run(solver,shape,conditions,min_clearance=.002,timeout_s=10)
    assert r.status == RunStatus.INVALID_INPUT
    assert r.coefficients.empty


def test_windows_paths_spaces_unicode_and_moment_reference(solver,real_run,tmp_path):
    relative = "espaços α"
    runner = SolverRunner(tmp_path,relative)
    other = AlulaSolver("inviscid",replace(REF,x_over_c=.4,y_over_c=.03,positive_sense="clockwise"),
                        node_executable=solver.node_executable,alula_directory=solver.alula_directory)
    result = runner.run(other,assembly(),CONDITIONS,min_clearance=.002,timeout_s=30)
    assert result.status == RunStatus.SUCCESS,result.diagnostics
    output = read_json(runner.execution_root/result.attempt_directory/"alula_output.json")
    assert result.coefficients.cm == output["upstream"]["coefficients"]["cm_clockwise"]
    assert output["moment_reference"] == other.moment_reference.to_dict()
    _,base_path,base_result = real_run
    force = np.array([sum(e[key] for e in read_json(base_path/"alula_output.json")["upstream"]["elements"])
                      for key in ("cx","cy")])
    expected_ccw = base_result.coefficients.cm-.15*force[1]+.03*force[0]
    assert result.coefficients.cm == pytest.approx(-expected_ccw,abs=3e-12)


@pytest.mark.skipif(sys.platform != "win32",reason="Windows process working-directory budget")
def test_windows_unsupported_long_cwd_rejected_explicitly(solver,tmp_path):
    runner = SolverRunner(tmp_path,"/".join(["long_path_"+"x"*40]*5))
    result = runner.run(solver,assembly(),CONDITIONS,min_clearance=.002,timeout_s=30)
    assert result.status == RunStatus.PREPARATION_FAILED
    assert any("240-character cwd budget" in d.message for d in result.diagnostics)
    assert result.coefficients.empty


@pytest.mark.parametrize("kind",["exit7","timeout_tree"])
def test_node_nonzero_exit_and_timeout_terminate_descendants(solver,tmp_path,kind):
    class Fault(AlulaSolver):
        def execute(self,case,context):
            if kind == "exit7":
                context.write_text("fault.mjs","process.exit(7);\n")
            else:
                context.write_text("child.mjs","console.log('CHILD_READY');setInterval(()=>{},1000);\n")
                context.write_text("fault.mjs", "import {spawn} from 'node:child_process';import {writeFileSync} from 'node:fs';import {fileURLToPath} from 'node:url';"
                    "const c=spawn(process.execPath,[fileURLToPath(new URL('./child.mjs',import.meta.url))],{stdio:['ignore','pipe','pipe']});"
                    "writeFileSync('pids.json',JSON.stringify([process.pid,c.pid]));c.stdout.pipe(process.stdout);setInterval(()=>{},1000);")
            return context.run_process([str(self.node_executable),node_path(context.path("fault.mjs"))],label="fault")
    faulty = Fault("inviscid",REF,node_executable=solver.node_executable,alula_directory=solver.alula_directory)
    runner = SolverRunner(tmp_path,"attempts")
    result = runner.run(faulty,assembly(),CONDITIONS,min_clearance=.002,timeout_s=1.5)
    assert result.status == (RunStatus.EXECUTION_FAILED if kind == "exit7" else RunStatus.TIMEOUT),result.diagnostics
    assert result.coefficients.empty
    path = runner.execution_root/result.attempt_directory
    evidence = read_json(path/"execution.json")
    assert "terminated" in evidence["tree_cleanup"]
    if kind == "timeout_tree":
        assert 'CHILD_READY' in (path/"fault.stdout.log").read_text()
        pids = json.loads((path/"pids.json").read_text())
        assert all(not process_running(pid) for pid in pids)
