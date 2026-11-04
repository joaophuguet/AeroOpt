"""Small shared numerical contracts; all lengths are in the supplied frame."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatPoints = NDArray[np.float64]
LENGTH_TOL = 1e-10
AREA_TOL = 1e-12


def immutable_points(value: ArrayLike) -> FloatPoints:
    """Copy into an immutable bytes buffer, including protection from setflags."""
    array = np.asarray(value, dtype=np.float64)
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


def finite_scalar(value: float, name: str) -> float:
    """Return a finite float or raise a parameter-specific ValueError."""
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result

