"""Conservative numerical ordering diagnostics, never mathematical error bounds.

No solver is executed here. Histories must be supplied explicitly, with the
same increasing mesh levels for both methods and candidates. A missing/failed
solution or an unsupported configuration never becomes a zero coefficient.
"""

from dataclasses import asdict, dataclass
import math
from numbers import Real


def _finite(value,name):
    if isinstance(value,bool) or not isinstance(value,Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


@dataclass(frozen=True)
class RefinementRecord:
    method: str
    panels: int
    cl: float | None
    cm: float | None
    valid: bool = True

    def __post_init__(self):
        if not isinstance(self.method,str) or not self.method:
            raise ValueError("method must be a nonempty name")
        if type(self.panels) is not int or self.panels < 1 or type(self.valid) is not bool:
            raise ValueError("positive integer panel count and boolean validity required")
        for name in ("cl","cm"):
            if getattr(self,name) is not None:
                object.__setattr__(self,name,_finite(getattr(self,name),name))


@dataclass(frozen=True)
class RankingPolicy:
    safety_factor: float = 2.
    comparison_floor: float = 1e-12

    def __post_init__(self):
        for name in ("safety_factor","comparison_floor"):
            object.__setattr__(self,name,_finite(getattr(self,name),name))
        if self.safety_factor < 1 or self.comparison_floor <= 0:
            raise ValueError("safety_factor >= 1 and comparison_floor > 0 required")


@dataclass(frozen=True)
class RefinementStudy:
    identifier: str
    records: tuple[RefinementRecord,...]
    limitations: tuple[str,...] = ()

    def __post_init__(self):
        if not isinstance(self.identifier,str) or not self.identifier:
            raise ValueError("study identifier must be nonempty")
        records = tuple(self.records)
        if not all(isinstance(r,RefinementRecord) for r in records):
            raise TypeError("records must contain RefinementRecord")
        if isinstance(self.limitations,str) or not all(isinstance(s,str) and s for s in self.limitations):
            raise ValueError("limitations must contain nonempty labels")
        for method in {r.method for r in records}:
            levels = [r.panels for r in records if r.method == method]
            if any(b <= a for a,b in zip(levels,levels[1:])):
                raise ValueError("mesh levels must be strictly increasing per method; no automatic sorting")
        object.__setattr__(self,"records",records)
        object.__setattr__(self,"limitations",tuple(self.limitations))

    def to_dict(self):
        return asdict(self)

    def analyze(self,coefficient,*,methods=("hess_smith","linear_vortex"),policy=RankingPolicy()):
        if coefficient not in ("cl","cm"):
            raise ValueError("only cl and cm histories are supported")
        if len(methods) < 2 or len(set(methods)) != len(methods):
            raise ValueError("at least two distinct methods required")
        histories,labels = {},list(self.limitations)
        for method in methods:
            history = []
            previous = None
            for record in (r for r in self.records if r.method == method):
                value = getattr(record,coefficient) if record.valid else None
                delta = None if value is None or previous is None else value-previous
                history.append({"panels": record.panels,"value": value,"delta_previous": delta})
                previous = value
            histories[method] = history
            if len(history) < 3:
                labels.append("insufficient_mesh_levels")
            if any(r["value"] is None for r in history):
                labels.append("failed_or_missing_coefficient")
        levels = [[r["panels"] for r in histories[m]] for m in methods]
        if any(level != levels[0] for level in levels[1:]):
            labels.append("mesh_levels_mismatch")
        indicator = discrepancy = None
        if not labels:
            last_steps = []
            for method in methods:
                history = histories[method]
                first,last = history[-2]["delta_previous"],history[-1]["delta_previous"]
                if abs(last) > abs(first)*(1+1e-8)+policy.comparison_floor:
                    labels.append("noncontracting_refinement")
                if last*first < -policy.comparison_floor**2:
                    labels.append("oscillatory_refinement")
                last_steps.append(abs(last))
            values = [histories[m][-1]["value"] for m in methods]
            discrepancy = max(values)-min(values)
            indicator = max(*last_steps,discrepancy)
        return {"coefficient": coefficient,"histories": histories,"indicator": indicator,
                "finest_method_discrepancy": discrepancy,"labels": sorted(set(labels)),
                "indicator_kind": "heuristic_not_error_bound","methods": list(methods)}


def compare_pair(a,b,coefficient,*,methods=("hess_smith","linear_vortex"),policy=RankingPolicy()):
    """Classify arithmetic ordering of one coefficient, not design desirability.

    U_c=max(last mesh increment in each method, finest method discrepancy).
    Distinguishable requires contracting nonoscillatory histories, common sign
    at both finest levels in both methods, and min|B-A| > factor*(U_A+U_B).
    Crossovers are inconsistent; absent/unsupported evidence is unresolved.
    This is a finite-study heuristic, not a universal bound or physical claim.
    """
    da,db = (s.analyze(coefficient,methods=methods,policy=policy) for s in (a,b))
    labels = set(da["labels"]+db["labels"])
    result = {"a": a.identifier,"b": b.identifier,"coefficient": coefficient,
              "classification": "numerically_unresolved","relation": None,"differences_b_minus_a": [],
              "margin_indicator": None,"a_indicator": da["indicator"],"b_indicator": db["indicator"],
              "policy": asdict(policy),"indicator_kind": "heuristic_not_error_bound",
              "scope": "numerical_ordering_only_no_physical_validation"}
    if da["indicator"] is None or db["indicator"] is None:
        result["labels"] = sorted(labels)
        return result
    for method in methods:
        ha,hb = da["histories"][method],db["histories"][method]
        if [r["panels"] for r in ha] != [r["panels"] for r in hb]:
            labels.add("candidate_mesh_levels_mismatch")
            result["labels"] = sorted(labels)
            return result
        for ra,rb in zip(ha[-2:],hb[-2:],strict=True):
            result["differences_b_minus_a"].append({"method": method,"panels": ra["panels"],
                                                   "difference": rb["value"]-ra["value"]})
    values = [r["difference"] for r in result["differences_b_minus_a"]]
    signs = {1 if v > policy.comparison_floor else -1 if v < -policy.comparison_floor else 0 for v in values}
    margin = policy.safety_factor*(da["indicator"]+db["indicator"])
    result["margin_indicator"] = margin
    if 1 in signs and -1 in signs:
        result["classification"] = "inconsistent_across_refinement_or_methods"
        labels.add("ordering_crossover")
    elif not labels and min(abs(v) for v in values) > max(margin,policy.comparison_floor):
        result["classification"] = "numerically_distinguishable"
        result["relation"] = "b_greater" if 1 in signs else "a_greater"
    else:
        labels.add("difference_not_resolved_by_indicator")
    result["labels"] = sorted(labels)
    return result
