"""Separate full-precision DAT contours and deterministic JSON metadata."""

import json
import platform
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import shapely

from aeroopt import __version__

from .airfoil import DAT_COORDINATE_HEADER, Airfoil
from .validation import ValidationReport

if TYPE_CHECKING:
    from .multielement import MultiElementAirfoil


class InvalidGeometryError(ValueError):
    """Export refused; report holds the complete geometric rejection reasons."""

    def __init__(self, report: ValidationReport) -> None:
        self.report = report
        super().__init__("Geometry rejected: " + "; ".join(issue.code for issue in report.issues))


def _source_metadata(airfoil: Airfoil) -> dict:
    return {
        "name": airfoil.name,
        "source_file": str(airfoil.source) if airfoil.source else None,
        "point_count": len(airfoil.original_points),
        "explicit_te_closure": airfoil.is_te_closed,
        "normalization": {
            "reference_method": airfoil.reference_method,
            "leading_edge_index": airfoil.leading_edge_index,
            "leading_edge_original": airfoil.leading_edge.tolist(),
            "trailing_edge_midpoint_original": airfoil.trailing_edge.tolist(),
            "chord_original_units": airfoil.chord,
            "input_units": "unspecified, same length unit on x and y",
        },
    }


def export_geometry(geometry: "MultiElementAirfoil", directory: str | Path, *,
                    configuration_id: str | None = None,
                    min_clearance: float = 0.0) -> Path:
    """Export valid geometry. Existing target files are never overwritten.

    All coordinates use the GLOBAL main-chord frame, including flap.dat.
    Reimport with read_dat/Airfoil.original_points to inspect that frame;
    Airfoil.normalized_points deliberately defines a new local frame.
    """
    directory = Path(directory)
    identifier = directory.name if configuration_id is None else configuration_id
    if not isinstance(identifier, str) or not identifier.strip() or any(c in identifier for c in "\r\n"):
        raise ValueError("configuration_id must be nonempty and single-line")
    # Freeze the exact arrays used for both validation and serialization.
    main, flap = geometry.main_points, geometry.flap_points
    from .validation import validate_elements

    report = validate_elements(main, flap, min_clearance=min_clearance)
    if not report.valid:
        raise InvalidGeometryError(report)
    targets = [directory / name for name in ("main.dat", "flap.dat", "metadata.json")]
    for target in targets:
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite {target}; choose a new output directory")
    metadata = {
        "schema_version": "aeroopt.geometry.v0.1",
        "configuration_id": identifier,
        "elements": {"main": _source_metadata(geometry.main), "flap": _source_metadata(geometry.flap)},
        "files": {"main": "main.dat", "flap": "flap.dat"},
        "parameters": geometry.parameters.to_dict(),
        "flap_global_leading_edge": list(geometry.parameters.global_reference),
        "units": {"coordinates": "x/c_main, y/c_main", "displacements": "c_main",
                  "nominal_position": "c_main", "deflection": "degree", "chord_ratio": "dimensionless",
                  "area": "c_main^2", "min_distance": "c_main"},
        "conventions": {
            "x_positive": "downstream", "y_positive": "up",
            "main_le": [0.0, 0.0], "main_te_midpoint": [1.0, 0.0],
            "positive_deflection": "clockwise",
            "flap_local_reference": "normalized leading edge (0,0)",
            "transform": "p_global = nominal + (dx,dy) + R(-deflection) @ (chord_ratio * p_normalized)",
            "te_closure": "straight segment between last and first sample for solid metrics",
            "dat_frame": "global normalized main-chord frame for BOTH files",
            "no_smoothing_or_resampling": True,
        },
        "validation": report.to_dict(),
        "software": {"aeroopt": __version__, "python": platform.python_version(),
                     "geos": shapely.geos_version_string,
                     **{name: version(name) for name in ("numpy", "shapely", "matplotlib")}},
    }
    serialized = json.dumps(metadata, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    directory.mkdir(parents=True, exist_ok=True)
    for target, element, points in zip(targets[:2], ("main", "flap"), (main, flap), strict=True):
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            np.savetxt(stream, points, fmt="%.17g",
                       header=f"{DAT_COORDINATE_HEADER} {element} {identifier}", comments="")
    with targets[2].open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(serialized)
    return directory

