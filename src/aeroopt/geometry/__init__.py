"""Public API for AeroOpt Geometry Engine v0.1."""

from .airfoil import Airfoil, AirfoilFormatError, read_dat
from .multielement import MultiElementAirfoil
from .transforms import FlapParameters, transform_flap
from .validation import ValidationIssue, ValidationReport, validate_elements
from .export import InvalidGeometryError
from .naca import naca4_points, naca4_surfaces

__all__ = ["Airfoil", "AirfoilFormatError", "read_dat", "MultiElementAirfoil",
           "FlapParameters", "transform_flap", "ValidationIssue", "ValidationReport",
           "validate_elements", "InvalidGeometryError", "naca4_points", "naca4_surfaces"]
