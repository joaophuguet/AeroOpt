"""Matplotlib figures in the same closed-contour convention as validation."""

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from aeroopt import __version__

if TYPE_CHECKING:
    from .multielement import MultiElementAirfoil


def _closed(points: np.ndarray) -> np.ndarray:
    return points if np.array_equal(points[0], points[-1]) else np.vstack((points, points[0]))


def plot_geometry(geometry: "MultiElementAirfoil", *, ax: Axes | None = None,
                  show_reference: bool = True) -> Axes:
    """Draw both elements, caps included; return axes without calling show()."""
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4), layout="constrained")
    for points, label, color in ((geometry.main_points, "Principal", "#235789"),
                                 (geometry.flap_points, "Flap", "#d45d00")):
        contour = _closed(points)
        ax.plot(contour[:, 0], contour[:, 1], color=color, label=label, linewidth=1.7)
        ax.fill(contour[:, 0], contour[:, 1], color=color, alpha=0.12)
    if show_reference:
        ax.scatter(*geometry.parameters.global_reference, marker="x", color="#9b2226",
                   s=45, label="Referência LE do flap", zorder=5)
    p = geometry.parameters
    ax.set_title(f"δ = {p.deflection_deg:g}° | cf/c = {p.chord_ratio:g} | "
                 f"Δx/c = {p.dx:g}, Δy/c = {p.dy:g}\n"
                 f"LE nominal/c = ({p.reference_position[0]:g}, {p.reference_position[1]:g})",
                 fontsize=10)
    ax.set(xlabel="x/c principal", ylabel="y/c principal")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="best")
    ax.margins(0.08)
    return ax


def compare_configurations(configurations: Sequence["MultiElementAirfoil"], *,
                           labels: Sequence[str] | None = None,
                           min_clearance: float = 0.0) -> Figure:
    """Compare with shared limits and equal axes; mark geometric rejections."""
    if not configurations:
        raise ValueError("at least one configuration is required")
    if labels is not None and len(labels) != len(configurations):
        raise ValueError("labels must match configuration count")
    columns = min(2, len(configurations))
    rows = (len(configurations) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, figsize=(7 * columns, 3.8 * rows),
                             squeeze=False, sharex=True, sharey=True, layout="constrained")
    for index, (geometry, ax) in enumerate(zip(configurations, axes.flat)):
        plot_geometry(geometry, ax=ax)
        report = geometry.validate(min_clearance=min_clearance)
        label = labels[index] if labels is not None else f"Configuração {index + 1}"
        status = "VÁLIDA" if report.valid else "REJEITADA"
        ax.set_title(f"{label} — {status}\n{ax.get_title()}", fontsize=10)
    for ax in list(axes.flat)[len(configurations):]:
        ax.set_visible(False)
    fig.suptitle(f"AeroOpt v{__version__} — comparação geométrica (sem análise aerodinâmica)", fontsize=13)
    return fig


def save_figure(figure: Figure, stem: str | Path) -> tuple[Path, Path]:
    """Save the same figure as PNG (180 dpi) and SVG; caller closes the figure."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = (stem.with_suffix(".png"), stem.with_suffix(".svg"))
    for path in paths:
        figure.savefig(path, dpi=180)
    return paths

