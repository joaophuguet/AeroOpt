"""Independent mathematical checks and regressions from the scientific review."""

from decimal import Decimal, localcontext
import hashlib
from pathlib import Path
import sys

import numpy as np
import pytest

from aeroopt.geometry import MultiElementAirfoil
from aeroopt.solvers import (AttemptContext, MomentReference, OperatingConditions, PanelSolver,
                            RunStatus, SolverResult, SolverRunner)
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers.serialization import read_json

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"examples"))
from _linear_vortex_reference import linear_influence, solve_linear_vortex
from _panel_reference import joukowski_exact
from _panel_review_tools import (assemble_dense_vectors, integrate_pressure, joukowski_blasius,
                                 joukowski_sample, subdivide_cap)
from verify_panel_solver import geometry


@pytest.mark.parametrize("distance",[1e6,1e9,1e12])
def test_far_field_source_log_cancellation_regression(distance):
    with localcontext() as context:
        context.prec = 60
        d = Decimal(int(distance))
        exact = float((d/(d-1)).ln())/(2*np.pi)
    actual = pm.source_velocity(np.array([[distance,0.]]),np.array([[0.,0.]]),np.array([[1.,0.]]))[0,0,0]
    assert actual == pytest.approx(exact,rel=1e-13)


@pytest.mark.parametrize("target",[[1e-20,1e-12],[1-1e-12,1e-20]])
def test_balanced_radial_formula_retains_regular_near_endpoint_values(target):
    x,y = target
    actual = pm.source_velocity(np.array([target]),np.array([[0.,0.]]),np.array([[1.,0.]]))[0,0,0]
    exact = (np.log(np.hypot(x,y))-np.log(np.hypot(x-1,y)))/(2*np.pi)
    assert np.isfinite(actual)
    assert actual == pytest.approx(exact,abs=1e-14)


@pytest.mark.parametrize("reverse",[False,True])
@pytest.mark.parametrize("closed",[False,True])
def test_projected_assembly_preserves_all_dense_vector_equations(reverse,closed):
    shape = geometry(81,closed=closed)
    points = (shape.main_points,shape.flap_points)
    meshes = tuple(pm.make_mesh(p[::-1] if reverse else p,40) for p in points)
    old,new = assemble_dense_vectors(meshes,6),pm.assemble(meshes,6)
    for a,b in ((old.matrix,new.matrix),(old.normal,new.normal),(old.tangent,new.tangent),(old.rhs,new.rhs)):
        np.testing.assert_allclose(a,b,rtol=2e-13,atol=2e-14)
    a,b = pm.solve(old,.25,0,"counterclockwise"),pm.solve(new,.25,0,"counterclockwise")
    assert a["accepted"] and b["accepted"]
    np.testing.assert_allclose([a["cl"],a["cm"]],[b["cl"],b["cm"]],rtol=2e-12,atol=2e-12)
    assert np.shares_memory(new.matrix,new.normal)


@pytest.mark.parametrize("angle",[0.,.71,2.5])
def test_independent_linear_vortex_coefficients_against_point_vortex_quadrature(angle):
    start = .3-.2j
    end = start+1.2*np.exp(1j*angle)
    targets = np.array([-1+1j,2-1j,1+2j])
    left,right = linear_influence(targets,start,end)
    nodes,weights = np.polynomial.legendre.leggauss(128)
    s = (nodes+1)/2
    z = start+s*(end-start)
    point_kernel = abs(end-start)/(2j*np.pi*(targets[:,None]-z[None,:]))
    qleft = np.sum(point_kernel*(1-s)*weights/2,axis=1)
    qright = np.sum(point_kernel*s*weights/2,axis=1)
    np.testing.assert_allclose(left,qleft,atol=2e-14)
    np.testing.assert_allclose(right,qright,atol=2e-14)


@pytest.mark.parametrize("orientation",[-1,1])
def test_independent_linear_vortex_exterior_self_limit(orientation):
    a,b = linear_influence(np.array([.5+0j]),0j,1+0j,self_index=0,orientation=orientation)
    assert a[0] == pytest.approx(orientation*.25-1j/(2*np.pi),abs=1e-15)
    assert b[0] == pytest.approx(orientation*.25+1j/(2*np.pi),abs=1e-15)


def test_independent_linear_vortex_far_field_moments_do_not_cancel():
    targets = np.array([1e9+1j,1e12+2j])
    left,right = linear_influence(targets,0j,1+0j)
    nodes,weights = np.polynomial.legendre.leggauss(64)
    s = (nodes+1)/2
    kernel = 1/(2j*np.pi*(targets[:,None]-s))
    expected_left = np.sum(kernel*(1-s)*weights/2,axis=1)
    expected_right = np.sum(kernel*s*weights/2,axis=1)
    np.testing.assert_allclose(left,expected_left,rtol=1e-13,atol=0)
    np.testing.assert_allclose(right,expected_right,rtol=1e-13,atol=0)


@pytest.mark.parametrize("alpha",[-8.,-4.,0.,4.,8.])
def test_blasius_closed_form_independently_confirms_pressure_reference(alpha):
    blasius = joukowski_blasius(alpha)
    pressure = joukowski_exact(alpha)
    assert blasius["cl"] == pytest.approx(pressure["cl_pressure_quadrature"],abs=2e-13)
    assert blasius["cm_ccw"] == pytest.approx(pressure["cm_ccw_quarter_chord"],abs=2e-13)


def test_exact_pressure_on_same_polygon_resolves_moment_sign():
    points,cp = joukowski_sample(160)
    mesh = pm.make_mesh(points,80)
    exact = joukowski_blasius(4)
    polygon = integrate_pressure((mesh,),cp,4)
    solved = pm.solve(pm.assemble((mesh,),4),.25,0,"counterclockwise")
    assert polygon["cm_ccw"] > 0
    assert abs(polygon["cm_ccw"]-exact["cm_ccw"]) < 2e-5
    assert abs(solved["cm"]-exact["cm_ccw"]) > 100*abs(polygon["cm_ccw"]-exact["cm_ccw"])


def test_linear_reference_joukowski_refinement_and_independence(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("independent reference called the production kernel")
    monkeypatch.setattr(pm,"assemble",forbidden)
    monkeypatch.setattr(pm,"source_velocity",forbidden)
    exact = joukowski_blasius(4)
    errors = []
    for count in (80,160,320):
        points,_ = joukowski_sample(count)
        reference = solve_linear_vortex((points,),4)
        assert reference["normal_residual"] < 1e-10
        assert reference["kutta_residual"] < 1e-12
        assert reference["cm_ccw"] > 0
        errors.append(abs(reference["cm_ccw"]-exact["cm_ccw"]))
    assert np.all(np.diff(errors) < 0)
    assert errors[-1] < .00025
    assert abs(reference["cl"]-exact["cl"]) < .0014


def test_independent_reference_is_traversal_invariant():
    shape = geometry(81)
    p = (shape.main_points,shape.flap_points)
    forward = solve_linear_vortex(p,6)
    backward = solve_linear_vortex(tuple(a[::-1] for a in p),6)
    assert forward["cl"] == pytest.approx(backward["cl"],abs=2e-12)
    assert forward["cm_ccw"] == pytest.approx(backward["cm_ccw"],abs=2e-12)
    np.testing.assert_allclose(forward["circulations"],backward["circulations"],atol=2e-12)


def test_independent_reference_rejects_ill_conditioning(monkeypatch):
    points,_ = joukowski_sample(80)
    monkeypatch.setattr(np.linalg,"cond",lambda matrix: 1e12)
    with pytest.raises(ArithmeticError,match="cond_2"):
        solve_linear_vortex((points,),4)


def test_independent_multielement_comparison_converges_in_cp_lift_and_moment():
    errors = []
    for count in (81,161,321):
        shape = geometry(count)
        p = (shape.main_points,shape.flap_points)
        meshes = tuple(pm.make_mesh(a,(count-1)//2) for a in p)
        hs = pm.solve(pm.assemble(meshes,6),.25,0,"counterclockwise")
        lv = solve_linear_vortex(p,6)
        assert lv["normal_residual"] < 1e-10 and lv["kutta_residual"] < 1e-12
        cp_errors = [float(np.sqrt(np.mean((np.array(a["cp"])-b["cp"])**2)))
                     for a,b in zip(hs["elements"],lv["elements"],strict=True)]
        errors.append([abs(hs["cl"]-lv["cl"]),abs(hs["cm"]-lv["cm_ccw"]),*cp_errors])
    assert np.all(np.diff(np.array(errors),axis=0) < 0)
    assert errors[-1][0] < .006 and errors[-1][1] < .0031
    assert errors[-1][2] < .06 and errors[-1][3] < .03
    np.testing.assert_allclose([e["circulation"] for e in hs["elements"]],lv["circulations"],atol=.003)


def test_cap_subdivision_study_preserves_polygon_and_exposes_sensitivity():
    shape = geometry(81,closed=False)
    before = shape.main_points.copy()
    meshes = tuple(pm.make_mesh(p,40) for p in (shape.main_points,shape.flap_points))
    split = tuple(subdivide_cap(m,8) for m in meshes)
    for old,new in zip(meshes,split,strict=True):
        assert np.sum(old.lengths) == pytest.approx(np.sum(new.lengths),abs=1e-14)
        np.testing.assert_array_equal(old.starts[:old.surface_count],new.starts[:new.surface_count])
        np.testing.assert_array_equal(old.ends[:old.surface_count],new.ends[:new.surface_count])
    first = pm.solve(pm.assemble(meshes,6),.25,0,"counterclockwise")
    refined = pm.solve(pm.assemble(split,6),.25,0,"counterclockwise")
    assert first["accepted"] and refined["accepted"]
    assert abs(first["cl"]-refined["cl"]) > .01
    assert refined["diagnostics"]["condition_number_2"] < 200
    np.testing.assert_array_equal(shape.main_points,before)


def test_finite_te_explicit_limitation_survives_runner_and_roundtrip(tmp_path):
    shape = MultiElementAirfoil.from_dat(ROOT/"examples/data/demo_main.dat",ROOT/"examples/data/demo_flap.dat",
                                        flap_reference_position=(1.10,-.05)).configure(flap_deflection_deg=25)
    runner = SolverRunner(tmp_path)
    result = runner.run(PanelSolver(MomentReference(.25,0,"clockwise","q_inf*c_main^2")),shape,
                        OperatingConditions(6,200000,0),min_clearance=.002)
    assert result.status == RunStatus.SUCCESS
    assert any("FINITE_TE_UNVERIFIED" in d.message for d in result.diagnostics)
    directory = runner.execution_root/result.attempt_directory
    output = read_json(directory/"panel_output.json")
    assert output["limitations"] == ["finite_te_closure_unverified"]
    assert SolverResult.load(directory/"result.json").to_dict() == result.to_dict()


def test_review_cli_preserves_provenance_and_reference_artifacts(tmp_path):
    context = AttemptContext(tmp_path,60)
    output = tmp_path/"review"
    evidence = context.run_process([sys.executable,str(ROOT/"examples/review_panel_solver.py"),
                                    "--quick","--output",str(output)],label="review")
    assert evidence.return_code == 0 and not evidence.timed_out
    report = read_json(output/"review.json")
    assert report["quick"]
    assert len(report["joukowski"]) == 4 and len(report["multielement"]) == 2
    assert report["joukowski_blasius"]["cm_ccw"] > 0
    for artifact in read_json(output/"artifacts.json")["artifacts"]:
        data = (output/artifact["path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == artifact["sha256"]
        assert len(data) == artifact["size_bytes"]
