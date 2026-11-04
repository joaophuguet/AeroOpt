"""Immutable explicit operating inputs; Reynolds is based on MAIN chord."""

import math
from dataclasses import dataclass
from numbers import Real

ALPHA_CONVENTION = "freestream_vector_(cos_alpha,sin_alpha)_in_global_xy"


def real_number(value: Real, name: str) -> float:
    """Accept real scalars, reject bool/text coercion and non-finite values."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, not bool or text")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} must be representable as finite float64") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class ConditionParameter:
    """Named scalar extension, e.g. transition location; no physical defaults.

    Numeric extensions require an explicit unit, including 'dimensionless'.
    Adapter validation assigns meaning to names and supported physical domains.
    """

    name: str
    value: float | int | str | bool
    unit: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("parameter name must be nonempty")
        if type(self.value) not in (float, int, str, bool):
            raise TypeError("parameter value must be float/int/str/bool")
        if type(self.value) in (float, int):
            real_number(self.value, self.name)
            if not isinstance(self.unit, str) or not self.unit.strip():
                raise ValueError("numeric extensions require an explicit unit")
        if self.unit is not None and (not isinstance(self.unit, str) or not self.unit.strip()):
            raise ValueError("unit must be nonempty when supplied")

    def to_dict(self) -> dict:
        return {"name": self.name, "value": self.value, "unit": self.unit}


@dataclass(frozen=True)
class OperatingConditions:
    """alpha_deg in degrees, Re_c_main > 0, Mach >= 0; no solver-specific limits."""

    alpha_deg: float
    reynolds: float
    mach: float
    additional: tuple[ConditionParameter, ...] = ()

    def __post_init__(self) -> None:
        for name in ("alpha_deg", "reynolds", "mach"):
            object.__setattr__(self, name, real_number(getattr(self, name), name))
        if self.reynolds <= 0 or self.mach < 0:
            raise ValueError("reynolds must be positive and mach nonnegative")
        params = tuple(self.additional)
        if not all(isinstance(item, ConditionParameter) for item in params):
            raise TypeError("additional must contain ConditionParameter objects")
        if len({item.name for item in params}) != len(params):
            raise ValueError("duplicate additional parameter names")
        object.__setattr__(self, "additional", tuple(sorted(params, key=lambda item: item.name)))

    def to_dict(self) -> dict:
        return {"alpha_deg": self.alpha_deg, "reynolds": self.reynolds, "mach": self.mach,
                "reynolds_reference": "main_chord",
                "alpha_convention": ALPHA_CONVENTION,
                "additional": [item.to_dict() for item in self.additional]}

    @classmethod
    def from_dict(cls, data: dict) -> "OperatingConditions":
        if set(data) != {"alpha_deg", "reynolds", "mach", "reynolds_reference", "alpha_convention", "additional"}:
            raise ValueError("unexpected operating-conditions fields")
        if data["reynolds_reference"] != "main_chord":
            raise ValueError("Reynolds reference must be main_chord")
        if data["alpha_convention"] != ALPHA_CONVENTION:
            raise ValueError("unsupported angle-of-attack convention")
        return cls(data["alpha_deg"], data["reynolds"], data["mach"],
                   tuple(ConditionParameter(**item) for item in data["additional"]))
