"""Independent analytic two-body verification; published airfoil remains uncertified."""

import hashlib
from pathlib import Path
import sys

import numpy as np
import pytest
from shapely.geometry import Point, Polygon

from aeroopt.geometry import validate_elements
from aeroopt.solvers import AttemptContext
from aeroopt.solvers import _panel_method as pm
from aeroopt.solvers.serialization import read_json, scientific_hash

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"examples"))
from _linear_vortex_reference import solve_linear_vortex
from _williams_reference import PrintedParameterCase, TwoCircleFlow, karman_trefftz
from verify_ranking_reliability import load_published, published_comparison
from _pressure_quadrature import pressure_quadrature


@pytest.mark.parametrize("direction,published_cw",[(30.,(5.169,1.778)),(10.,(3.154,.655))])
def test_circle_circulations_match_report_rounding_and_boundary_equations(direction,published_cw):
    case = PrintedParameterCase(direction)
    np.testing.assert_allclose(-case.gamma_raw,published_cw,rtol=0,atol=5e-4)
    checks = case.flow_checks()
    assert checks["circle_normal_velocity_max"] < 2e-13
    assert checks["stagnation_velocity_max"] < 2e-13
    assert checks["circulation_half_series_difference_max"] < 2e-13


def test_two_circle_solution_rotates_covariantly():
    case = PrintedParameterCase()
    rotation = np.exp(.73j)
    flow = TwoCircleFlow(tuple(c*rotation for c in case.centers),case.flow.radii)
    te = np.asarray(case.trailing_points)*rotation
    gamma = flow.circulations(te,.73)
    np.testing.assert_allclose(gamma,case.gamma_raw,atol=1e-13)
    z = np.array([3+2j,-3+1j])
    np.testing.assert_allclose(flow.velocity(z*rotation,.73,gamma),
                               case.flow.velocity(z,0,case.gamma_raw)/rotation,atol=1e-14)


@pytest.mark.parametrize("direction",[10.,30.])
def test_normalization_reference_and_branch_cuts_are_explicit(direction):
    case = PrintedParameterCase(direction)
    points = case.contours(160)
    assert validate_elements(*points,min_clearance=0.).valid
    np.testing.assert_allclose(points[0][80],[0,0],atol=2e-15)
    np.testing.assert_allclose(points[0][0],[1,0],atol=2e-15)
    for p in points:
        np.testing.assert_array_equal(p[0],p[-1])
    assert case.alpha_deg == pytest.approx(-np.rad2deg(np.angle(case.normalization[1])))
    # Both principal-log cuts must lie inside their respective solid contours.
    first = Polygon(np.column_stack((
        (case.centers[0]+1.096*np.exp(1j*np.linspace(0,2*np.pi,1001))).real,
        (case.centers[0]+1.096*np.exp(1j*np.linspace(0,2*np.pi,1001))).imag)))
    assert all(first.contains(Point(x,0)) for x in np.linspace(-.999,.999,50))
    theta = np.linspace(0,2*np.pi,2049)
    z = case.centers[1]+(case.trailing_points[1]-case.centers[1])*np.exp(1j*theta)
    mapped = karman_trefftz(z,1.,case.n)[0]
    flap = Polygon(np.column_stack((mapped.real,mapped.imag)))
    origin,c = case.second_map
    cut = origin+np.linspace(-.999,.999,50)*c*np.exp(-1j*np.deg2rad(direction))
    assert all(flap.contains(Point(p.real,p.imag)) for p in cut)


def test_conformal_derivative_and_critical_limits():
    case = PrintedParameterCase()
    z = np.array([3+2j,-3+1j,.2+2j])
    value,derivative = case.map_raw(z)
    numerical = (case.map_raw(z+1e-5)[0]-case.map_raw(z-1e-5)[0])/2e-5
    np.testing.assert_allclose(derivative,numerical,rtol=2e-9,atol=2e-10)
    value,derivative = karman_trefftz(np.array([-1,1]),1,case.n)
    np.testing.assert_allclose(value,[-case.n,case.n])
    np.testing.assert_array_equal(derivative,[0j,0j])
    for body in (0,1):
        _,cp,_ = case.surface(body,np.array([0.,2*np.pi]))
        np.testing.assert_array_equal(cp,[1,1])


@pytest.mark.parametrize("direction",[10.,30.])
def test_exact_pressure_circulation_and_independent_blasius_agree(direction,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("reference must not call panel formulations")
    monkeypatch.setattr(pm,"assemble",forbidden)
    monkeypatch.setattr(sys.modules["_linear_vortex_reference"],"solve_linear_vortex",forbidden)
    case = PrintedParameterCase(direction)
    coarse,fine = case.reference(4096),case.reference(8192)
    assert abs(fine["cl"]-fine["cl_circulation"]) < 3e-9
    assert abs(fine["cl_blasius"]-fine["cl_circulation"]) < 2e-13
    assert abs(fine["cm_ccw"]-fine["cm_blasius"]) < 3e-9
    assert abs(fine["cl"]-coarse["cl"]) < 7e-9
    assert abs(fine["cm_ccw"]-coarse["cm_ccw"]) < 7e-9
    assert fine["published_airfoil_reproduction_certified"] is False


@pytest.mark.parametrize("direction",[10.,30.])
def test_both_formulations_converge_to_declared_two_body_reference(direction):
    case = PrintedParameterCase(direction)
    exact = case.reference(8192)
    errors = {m:[] for m in ("hs","lv")}
    for count in (40,80,160):
        points = case.contours(count)
        meshes = tuple(pm.make_mesh(p,count//2) for p in points)
        hs = pm.solve(pm.assemble(meshes,case.alpha_deg),.25,0,"counterclockwise")
        lv = solve_linear_vortex(points,case.alpha_deg)
        assert hs["accepted"] and lv["accepted"]
        assert hs["diagnostics"]["kutta_max"] < 1e-9
        assert lv["kutta_residual"] < 1e-9
        for name,result in (("hs",hs),("lv",lv)):
            cp_errors = []
            for body in (0,1):
                nodes = case.parameter_nodes(body,count)
                cp = case.surface(body,(nodes[:-1]+nodes[1:])/2)[1]
                cp_errors.append(np.sqrt(np.mean((result["elements"][body]["cp"]-cp)**2)))
            gamma = np.array([e["circulation"] for e in result["elements"]])
            errors[name].append([abs(result["cl"]-exact["cl_circulation"]),
                                 abs(result["cm" if name == "hs" else "cm_ccw"]-exact["cm_blasius"]),
                                 max(abs(gamma-exact["element_circulations_ccw"])),*cp_errors])
    for history in errors.values():
        e = np.array(history)
        assert np.all(e[1:] < e[:-1])
        # Case-specific regression limits at 160 panels/body; no universal bound.
        assert e[-1,0] < .035
        assert e[-1,1] < .025
        assert e[-1,2] < .006
        assert max(e[-1,3:]) < .12


def test_research_quadrature_records_unsolved_boundary_and_preserves_solution():
    case = PrintedParameterCase()
    meshes = tuple(pm.make_mesh(p,40) for p in case.contours(80))
    solved = pm.solve(pm.assemble(meshes,case.alpha_deg),.25,0,"counterclockwise")
    original = scientific_hash(solved)
    diagnostic = pressure_quadrature(meshes,solved,case.alpha_deg,2)
    assert diagnostic["normal_velocity_max"] > 1e-3
    assert diagnostic["limitation"] == "research_only_no_impermeability_solution_at_quadrature_points"
    assert scientific_hash(solved) == original


@pytest.mark.parametrize("direction,key",[(30.,"A"),(10.,"B")])
def test_local_published_original_is_pinned_and_mismatch_not_hidden(request,direction,key):
    directory = request.config.getoption("--williams-data")
    if directory is None:
        pytest.skip("Williams originals/checks not redistributed; use --williams-data")
    _,data = load_published(directory)
    result = published_comparison(PrintedParameterCase(direction),data[key])
    assert result["circle_table_reproduced_at_printed_precision"]
    assert result["airfoil_circulation_max_error"] > 5e-5
    assert not result["airfoil_reproduction_certified"]


def test_quick_harness_preserves_inputs_hashes_and_uncertainty(tmp_path):
    context = AttemptContext(tmp_path,120)
    output = tmp_path/"evidence"
    argv = [sys.executable,str(ROOT/"examples/verify_ranking_reliability.py"),"--quick","--output",str(output)]
    execution = context.run_process(argv,label="quick")
    assert execution.return_code == 0, (tmp_path/"quick.stderr.log").read_text()
    report = read_json(output/"verification.json")
    assert report["quick"]
    assert not report["acceptance"]["independent_published_benchmark_acceptance_satisfied"]
    assert not report["acceptance"]["preliminary_automatic_optimization_ready"]
    assert len(report["benchmark"]) == 6
    assert len(report["ranking"]["pairs"]) == 20
    assert any(p["classification"] == "numerically_unresolved" for p in report["ranking"]["pairs"])
    for row in read_json(output/"artifacts.json")["artifacts"]:
        data = (output/row["path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == row["sha256"]
    benchmark = read_json(output/"cases/benchmark_A_like_40.json")
    assert scientific_hash(benchmark["input"]) == benchmark["input_sha256"]
    second = context.run_process(argv,label="overwrite")
    assert second.return_code != 0
    assert "new/empty" in (tmp_path/"overwrite.stderr.log").read_text()
