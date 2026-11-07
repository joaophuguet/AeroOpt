"""Numerical ordering policy, including counterexamples to overconfident ranks."""

import json

import pytest

from aeroopt.verification import RefinementRecord, RefinementStudy, RankingPolicy, compare_pair


def study(name,hs=(1.1,1.04,1.02),lv=(1.08,1.03,1.015),offset=0.,limitations=()):
    records = tuple(RefinementRecord(method,n,value+offset,-value-offset)
                    for method,values in (("hess_smith",hs),("linear_vortex",lv))
                    for n,value in zip((40,80,160),values,strict=True))
    return RefinementStudy(name,records,limitations)


def test_histories_steps_and_indicator_are_explicit_and_serializable():
    s = study("a")
    d = s.analyze("cl")
    assert d["indicator"] == pytest.approx(.02)
    assert d["finest_method_discrepancy"] == pytest.approx(.005)
    assert d["histories"]["hess_smith"][0]["delta_previous"] is None
    assert d["histories"]["hess_smith"][-1]["delta_previous"] == pytest.approx(-.02)
    assert d["indicator_kind"] == "heuristic_not_error_bound"
    assert not d["labels"]
    assert json.loads(json.dumps(s.to_dict()))["identifier"] == "a"


@pytest.mark.parametrize("coefficient",["cl","cm"])
def test_resolved_order_reverses_with_candidates_and_moment_is_arithmetic(coefficient):
    a,b = study("a"),study("b",offset=.2)
    result = compare_pair(a,b,coefficient)
    assert result["classification"] == "numerically_distinguishable"
    expected = "b_greater" if coefficient == "cl" else "a_greater"
    assert result["relation"] == expected
    reverse = compare_pair(b,a,coefficient)
    assert reverse["relation"] != expected
    assert reverse["margin_indicator"] == pytest.approx(result["margin_indicator"])


def test_nearby_values_and_common_discretization_bias_are_unresolved():
    a = study("a")
    r = compare_pair(a,study("b",offset=.001),"cl")
    assert r["classification"] == "numerically_unresolved"
    assert r["relation"] is None
    # Even perfect method agreement cannot defeat a large refinement increment.
    a,b = study("a",hs=(2,1.5,1.25),lv=(2,1.5,1.25)),study(
        "b",hs=(2,1.5,1.25),lv=(2,1.5,1.25),offset=.1)
    assert compare_pair(a,b,"cl")["classification"] == "numerically_unresolved"


@pytest.mark.parametrize("b",[
    study("b",hs=(1.11,1.05,1.03),lv=(1.07,1.02,1.005)),
    study("b",hs=(1.11,1.05,1.01),lv=(1.09,1.04,1.005)),
])
def test_method_or_mesh_ordering_crossover_is_inconsistent(b):
    result = compare_pair(study("a"),b,"cl")
    assert result["classification"] == "inconsistent_across_refinement_or_methods"
    assert "ordering_crossover" in result["labels"]


@pytest.mark.parametrize("hs,label",[
    ((1.,1.01,1.2),"noncontracting_refinement"),
    ((1.,1.1,1.05),"oscillatory_refinement"),
])
def test_bad_refinement_does_not_grant_confident_ordering(hs,label):
    a = study("a",hs=hs)
    result = compare_pair(a,study("b",offset=5.),"cl")
    assert result["classification"] == "numerically_unresolved"
    assert label in result["labels"]


@pytest.mark.parametrize("kind",["failed","missing_cm","missing_method","few_levels","mismatched_methods","mismatched_candidates"])
def test_incomplete_evidence_never_becomes_zero_or_resolved(kind):
    a,b = study("a"),study("b",offset=1.)
    rows = list(b.records)
    if kind == "failed":
        r = rows[-1]
        rows[-1] = RefinementRecord(r.method,r.panels,r.cl,r.cm,False)
    elif kind == "missing_cm":
        r = rows[-1]
        rows[-1] = RefinementRecord(r.method,r.panels,r.cl,None)
    elif kind == "missing_method":
        rows = rows[:3]
    elif kind == "few_levels":
        rows = [r for r in rows if r.panels != 40]
    elif kind == "mismatched_methods":
        r = rows[-1]
        rows[-1] = RefinementRecord(r.method,320,r.cl,r.cm)
    else:
        rows = [RefinementRecord(r.method,2*r.panels,r.cl,r.cm) for r in rows]
    b = RefinementStudy("b",tuple(rows))
    result = compare_pair(a,b,"cm")
    assert result["classification"] == "numerically_unresolved"
    assert result["relation"] is None
    if kind in ("failed","missing_cm"):
        assert b.analyze("cm")["histories"]["linear_vortex"][-1]["value"] is None


def test_unsupported_finite_te_has_no_indicator_even_with_smooth_histories():
    a = study("finite",limitations=("FINITE_TE_UNVERIFIED",))
    assert a.analyze("cl")["indicator"] is None
    result = compare_pair(a,study("other",offset=9),"cl")
    assert result["classification"] == "numerically_unresolved"
    assert "FINITE_TE_UNVERIFIED" in result["labels"]


@pytest.mark.parametrize("value",[float("nan"),float("inf"),True,"1"])
def test_nonfinite_or_nonreal_coefficients_rejected(value):
    with pytest.raises(ValueError):
        RefinementRecord("hess_smith",40,value,None)


@pytest.mark.parametrize("levels",[(40,40),(80,40)])
def test_levels_never_silently_sorted_or_deduplicated(levels):
    with pytest.raises(ValueError,match="strictly increasing"):
        RefinementStudy("a",tuple(RefinementRecord("hs",n,1.,None) for n in levels))


@pytest.mark.parametrize("kwargs",[{"safety_factor":.5},{"comparison_floor":0},{"safety_factor":float("inf")}])
def test_invalid_policy_rejected(kwargs):
    with pytest.raises(ValueError):
        RankingPolicy(**kwargs)


def test_coefficient_and_method_contract():
    with pytest.raises(ValueError):
        study("a").analyze("cd")
    with pytest.raises(ValueError):
        study("a").analyze("cl",methods=("hs","hs"))
