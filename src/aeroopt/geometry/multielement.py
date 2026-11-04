"""Immutable main + fixed single flap geometry configurations."""

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

from ._numeric import FloatPoints
from .airfoil import Airfoil
from .transforms import FlapParameters, transform_flap
from .validation import ValidationReport, validate_elements

if TYPE_CHECKING:
    from matplotlib.axes import Axes


@dataclass(frozen=True, eq=False)
class MultiElementAirfoil:
    """Two independent source contours plus one immutable parameter set."""

    main: Airfoil
    flap: Airfoil
    parameters: FlapParameters = FlapParameters()

    @classmethod
    def from_dat(cls, main: str | Path, flap: str | Path, *,
                 flap_chord_ratio: float = 0.30,
                 flap_reference_position: tuple[float, float] = (0.95, -0.02),
                 main_leading_edge_index: int | None = None,
                 flap_leading_edge_index: int | None = None) -> "MultiElementAirfoil":
        """Import and normalize independent sources (point counts may differ)."""
        return cls(
            Airfoil.from_dat(main, leading_edge_index=main_leading_edge_index),
            Airfoil.from_dat(flap, leading_edge_index=flap_leading_edge_index),
            FlapParameters(chord_ratio=flap_chord_ratio,
                           reference_position=flap_reference_position),
        )

    def configure(self, *, flap_deflection_deg: float = 0.0,
                  flap_dx: float = 0.0, flap_dy: float = 0.0) -> "MultiElementAirfoil":
        """Return a NEW configuration; supplied values are absolute, not increments.

        Omitted design variables reset to zero. Assign the returned object.
        """
        return replace(self, parameters=replace(
            self.parameters, deflection_deg=flap_deflection_deg, dx=flap_dx, dy=flap_dy
        ))

    @property
    def main_points(self) -> FloatPoints:
        """Main contour in the normalized global frame."""
        return self.main.normalized_points

    @property
    def flap_points(self) -> FloatPoints:
        """Generate placed flap afresh from its normalized reference contour."""
        return transform_flap(self.flap.normalized_points, self.parameters)

    def validate(self, *, min_clearance: float = 0.0) -> ValidationReport:
        """Return reasons and metrics, with all lengths in main chord units."""
        return validate_elements(self.main_points, self.flap_points, min_clearance=min_clearance)

    def export(self, directory: str | Path, *, configuration_id: str | None = None,
               min_clearance: float = 0.0) -> Path:
        """Validate again with the explicit export clearance, then export separately."""
        from .export import export_geometry

        return export_geometry(self, directory, configuration_id=configuration_id,
                               min_clearance=min_clearance)

    def plot(self, *, ax: "Axes | None" = None, show_reference: bool = True) -> "Axes":
        """Plot without opening a GUI; save/show through the returned axes."""
        from .plotting import plot_geometry

        return plot_geometry(self, ax=ax, show_reference=show_reference)

