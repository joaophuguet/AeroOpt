"""Strict single-contour DAT import and shape-preserving normalization."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike
from shapely.geometry import MultiPoint, Polygon

from ._numeric import AREA_TOL, LENGTH_TOL, FloatPoints, immutable_points

DAT_COORDINATE_HEADER = "AeroOpt_COORDINATES_v0.1"


class AirfoilFormatError(ValueError):
    """The input is not a supported finite, nondegenerate single contour."""


def read_dat(path: str | Path) -> tuple[str, FloatPoints]:
    """Read whitespace-separated x y rows and at most one initial text header.

    No count tables, comments within data, comma delimiters, or multiple blocks
    are interpreted. Values are preserved; no normalization occurs here.
    """
    path = Path(path)
    rows: list[tuple[float, float]] = []
    header: str | None = None
    for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        tokens = line.split()
        values: list[float] = []
        for token in tokens:
            try:
                values.append(float(token))
            except ValueError:
                pass
        try:
            float(tokens[0])
            first_is_number = True
        except ValueError:
            first_is_number = False
        if (not rows and header is None and not first_is_number
                and (tokens[0][0].isalpha() or tokens[0].startswith("#"))):
            header = line
            continue
        if len(tokens) != 2 or len(values) != 2:
            raise AirfoilFormatError(
                f"{path}:{line_number}: expected exactly two numeric coordinates"
            )
        if not np.isfinite(values).all():
            raise AirfoilFormatError(f"{path}:{line_number}: non-finite coordinate")
        explicitly_coordinates = header is not None and header.split()[0] == DAT_COORDINATE_HEADER
        if not rows and not explicitly_coordinates and all(v > 1 and v.is_integer() for v in values):
            raise AirfoilFormatError(
                f"{path}:{line_number}: ambiguous integer pair/count table; "
                "Lednicer format is unsupported (convert explicitly to TE-LE-TE); "
                f"for actual coordinates, declare header {DAT_COORDINATE_HEADER}"
            )
        rows.append((values[0], values[1]))
    if len(rows) < 3:
        raise AirfoilFormatError(f"{path}: expected at least three coordinate rows")
    return header or path.stem, immutable_points(rows)


@dataclass(frozen=True, eq=False)
class Airfoil:
    """Immutable sampled contour with explicit TE-LE-TE ordering contract.

    The trailing reference is the midpoint of first/last rows. The default
    leading reference is the interior sample farthest from that midpoint.
    Supply leading_edge_index to override this sampled proxy. Original units
    may be arbitrary but must be equal on both axes. No points are repaired.
    """

    original_points: FloatPoints
    name: str = "airfoil"
    source: Path | None = None
    leading_edge_index: int | None = None
    normalized_points: FloatPoints = field(init=False, repr=False)
    leading_edge: FloatPoints = field(init=False, repr=False)
    trailing_edge: FloatPoints = field(init=False, repr=False)
    chord: float = field(init=False)
    reference_method: str = field(init=False)

    def __post_init__(self) -> None:
        points = np.asarray(self.original_points, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
            raise AirfoilFormatError("expected an (N, 2) contour with N >= 3")
        if not np.isfinite(points).all():
            raise AirfoilFormatError("non-finite coordinate")
        span = float(np.linalg.norm(np.ptp(points, axis=0)))
        if not np.isfinite(span) or span <= 0:
            raise AirfoilFormatError("degenerate contour: zero or unrepresentable span")
        # Work in a relative frame so duplicate/degeneracy tests are scale independent.
        relative = (points - points[0]) / span
        for i in range(len(points) - 1):
            distances = np.linalg.norm(relative[i + 1:] - relative[i], axis=1)
            for offset in np.flatnonzero(distances <= LENGTH_TOL):
                j = i + 1 + int(offset)
                if i == 0 and j == len(points) - 1 and np.array_equal(points[0], points[-1]):
                    continue
                raise AirfoilFormatError(f"problematic duplicate/near-duplicate rows {i} and {j}")
        if MultiPoint(relative).convex_hull.area <= AREA_TOL:
            raise AirfoilFormatError("degenerate contour: collinear or negligible area")
        te = points[0] + (points[-1] - points[0]) / 2
        index = self.leading_edge_index
        method = "farthest_sample_from_te_midpoint" if index is None else "explicit_leading_edge_index"
        if index is None:
            index = int(np.argmax(np.linalg.norm((points - te) / span, axis=1)))
        if isinstance(index, bool) or not isinstance(index, (int, np.integer)):
            raise AirfoilFormatError("leading_edge_index must be an integer")
        if not 0 < index < len(points) - 1:
            raise AirfoilFormatError("leading edge must be an interior row; expected TE-LE-TE order")
        le = points[index]
        chord = float(np.linalg.norm(te - le))
        if chord / span <= LENGTH_TOL:
            raise AirfoilFormatError("degenerate chord")
        ex = (te - le) / chord
        ey = np.array([-ex[1], ex[0]])
        normalized = ((points - le) / chord) @ np.column_stack((ex, ey))
        if not np.isfinite(normalized).all():
            raise AirfoilFormatError("normalization produced unrepresentable coordinates")
        object.__setattr__(self, "original_points", immutable_points(points))
        object.__setattr__(self, "normalized_points", immutable_points(normalized))
        object.__setattr__(self, "leading_edge", immutable_points(le))
        object.__setattr__(self, "trailing_edge", immutable_points(te))
        object.__setattr__(self, "chord", chord)
        object.__setattr__(self, "leading_edge_index", int(index))
        object.__setattr__(self, "reference_method", method)
        if self.source is not None:
            object.__setattr__(self, "source", Path(self.source).resolve())

    @classmethod
    def from_dat(cls, path: str | Path, *, leading_edge_index: int | None = None) -> "Airfoil":
        """Import one supported DAT; preserve its source path and original samples."""
        name, points = read_dat(path)
        return cls(points, name=name, source=Path(path), leading_edge_index=leading_edge_index)

    @classmethod
    def from_points(cls, points: ArrayLike, *, name: str = "airfoil",
                    leading_edge_index: int | None = None) -> "Airfoil":
        """Build from an ordered contour, using the same checks as DAT import."""
        return cls(np.asarray(points, dtype=np.float64), name=name,
                   leading_edge_index=leading_edge_index)

    @property
    def is_te_closed(self) -> bool:
        """Whether the input explicitly repeats its first point at the end."""
        return bool(np.array_equal(self.original_points[0], self.original_points[-1]))

    @property
    def normalized_area(self) -> float | None:
        """Closed solid area in chord squared; None for invalid polygons."""
        polygon = Polygon(self.normalized_points)
        return float(polygon.area) if polygon.is_valid and polygon.area > AREA_TOL else None

    @property
    def area(self) -> float | None:
        """Area in original input units squared; open TE receives a straight cap."""
        area = self.normalized_area
        return None if area is None else area * self.chord**2

