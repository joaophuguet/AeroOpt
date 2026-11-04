"""Deterministic, nonaccumulating flap placement in main-chord units."""

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike

from ._numeric import FloatPoints, finite_scalar, immutable_points


@dataclass(frozen=True)
class FlapParameters:
    """Angles in degrees; all coordinates/displacements in main chord c=1."""

    chord_ratio: float = 0.30
    reference_position: tuple[float, float] = (0.95, -0.02)
    deflection_deg: float = 0.0
    dx: float = 0.0
    dy: float = 0.0

    def __post_init__(self) -> None:
        for name in ("chord_ratio", "deflection_deg", "dx", "dy"):
            object.__setattr__(self, name, finite_scalar(getattr(self, name), name))
        if self.chord_ratio <= 0:
            raise ValueError("chord_ratio must be positive")
        if len(self.reference_position) != 2:
            raise ValueError("reference_position must contain two coordinates")
        object.__setattr__(self, "reference_position", tuple(
            finite_scalar(v, "reference_position") for v in self.reference_position
        ))
        if not np.isfinite(self.global_reference).all():
            raise ValueError("reference_position plus displacement must be finite")

    @property
    def global_reference(self) -> tuple[float, float]:
        """Actual global LE = nominal LE + (dx, dy)."""
        return (self.reference_position[0] + self.dx, self.reference_position[1] + self.dy)

    def to_dict(self) -> dict:
        """Return serializable parameters without losing units in field names."""
        return asdict(self)


def transform_flap(normalized_points: ArrayLike, parameters: FlapParameters) -> FloatPoints:
    """Scale about local LE=(0,0), rotate clockwise, then place the global LE.

    Always pass the original normalized flap, never a previously placed flap.
    """
    points = np.asarray(normalized_points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("expected finite (N, 2) normalized points")
    theta = np.deg2rad(parameters.deflection_deg % 360.0)
    cosine, sine = np.cos(theta), np.sin(theta)
    rotation = np.array([[cosine, sine], [-sine, cosine]])
    with np.errstate(over="ignore", invalid="ignore"):
        result = (points * parameters.chord_ratio) @ rotation.T + parameters.global_reference
    if not np.isfinite(result).all():
        raise ValueError("flap transform produced unrepresentable coordinates")
    return immutable_points(result)

