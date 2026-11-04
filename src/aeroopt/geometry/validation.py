"""Structured contour/solid validation; no geometric repair is performed."""

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike
from shapely.geometry import LinearRing, Polygon
from shapely.validation import explain_validity

from ._numeric import AREA_TOL, LENGTH_TOL, finite_scalar


@dataclass(frozen=True)
class ValidationIssue:
    """Stable machine-readable code plus contextual human-readable detail."""

    code: str
    message: str
    element: str | None = None


@dataclass(frozen=True)
class ValidationReport:
    """Distances are in main chord units; absent metrics are undefined."""

    issues: tuple[ValidationIssue, ...]
    min_distance: float | None
    min_clearance: float
    overlap_area: float | None
    length_tolerance: float = LENGTH_TOL
    area_tolerance: float = AREA_TOL

    @property
    def valid(self) -> bool:
        return not self.issues

    def to_dict(self) -> dict:
        """Return a JSON-compatible report, with no NaN placeholders."""
        return {"valid": self.valid, **asdict(self), "units": "main_chord",
                "element_tolerance_frame": "own_bbox_diagonal",
                "area_tolerance_units": "own_bbox_diagonal_squared",
                "duplicate_tolerance_units": "own_bbox_diagonal",
                "contact_clearance_tolerance_units": "main_chord",
                "aerodynamic_gap": None, "aerodynamic_overlap": None,
                "gap_overlap_status": "pending_definition"}


def _polygon(points: ArrayLike, name: str,
             issues: list[ValidationIssue]) -> Polygon | None:
    array = np.asarray(points, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 2 or len(array) < 3:
        issues.append(ValidationIssue("degenerate_contour", "Expected N >= 3 coordinate pairs", name))
        return None
    if not np.isfinite(array).all():
        issues.append(ValidationIssue("non_finite", "Non-finite coordinates", name))
        return None
    span = float(np.hypot(*np.ptp(array, axis=0)))
    if not np.isfinite(span) or span <= 0:
        issues.append(ValidationIssue("degenerate_contour", "Zero or unrepresentable contour span", name))
        return None
    relative = (array - array[0]) / span
    # Exact final closure is allowed, but short edges/repeated internal vertices are not.
    unique = relative[:-1] if np.array_equal(array[0], array[-1]) else relative
    for i in range(len(unique) - 1):
        if np.any(np.linalg.norm(unique[i + 1:] - unique[i], axis=1) <= LENGTH_TOL):
            issues.append(ValidationIssue("duplicate_points", "Repeated/near-repeated contour vertices", name))
            break
    if len(unique) < 3:
        issues.append(ValidationIssue("degenerate_contour", "Fewer than three distinct vertices", name))
        return None
    polygon = Polygon(array)  # GEOS closes a finite-thickness TE with one straight segment.
    # Shape quality is intrinsic: scaling an otherwise resolved flap by r must
    # not multiply its degeneracy threshold by 1/r². Contact stays global.
    local_polygon = Polygon(relative)
    ring = LinearRing(relative)
    if not ring.is_simple:
        issues.append(ValidationIssue("self_intersection", explain_validity(local_polygon), name))
    if local_polygon.area <= AREA_TOL:
        issues.append(ValidationIssue("degenerate_contour", "Area / own span² <= area tolerance", name))
    if not polygon.is_valid and ring.is_simple:
        issues.append(ValidationIssue("invalid_contour", explain_validity(polygon), name))
    return polygon


def validate_elements(main_points: ArrayLike, flap_points: ArrayLike, *,
                      min_clearance: float = 0.0) -> ValidationReport:
    """Check closed contours and solids; contact is rejected even at zero clearance.

    Boundary distance remains distinct from solid overlap (including containment).
    Pair metrics are undefined when either element is invalid. No buffer(0) repair.
    """
    clearance = finite_scalar(min_clearance, "min_clearance")
    if clearance < 0:
        raise ValueError("min_clearance must be nonnegative")
    issues: list[ValidationIssue] = []
    main = _polygon(main_points, "main", issues)
    flap = _polygon(flap_points, "flap", issues)
    if issues or main is None or flap is None:
        return ValidationReport(tuple(issues), None, clearance, None)
    distance = float(main.boundary.distance(flap.boundary))
    overlap = float(main.intersection(flap).area)
    if main.boundary.intersects(flap.boundary):
        issues.append(ValidationIssue("boundary_intersection", "Element contours intersect or touch"))
    if overlap > 0.0:
        issues.append(ValidationIssue("solid_overlap", "Solid regions overlap, including containment"))
    if distance + LENGTH_TOL < clearance:
        issues.append(ValidationIssue("insufficient_clearance", "Boundary distance below requested clearance"))
    if not issues and distance <= LENGTH_TOL:
        issues.append(ValidationIssue("near_contact", "Boundary separation <= length tolerance"))
    return ValidationReport(tuple(issues), distance, clearance, overlap)

