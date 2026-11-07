"""Analytical constant source / CCW vortex panels, in main-chord units.

Private numerical kernel. No point reordering, resampling or geometry repair.
See docs/PANEL_METHOD.md for equations, signs and numerical acceptance policy.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from shapely.geometry import Polygon

from aeroopt.geometry._numeric import immutable_points

Array = NDArray[np.float64]
MIN_LENGTH_RATIO = 1e-10
MIN_CLEARANCE = 1e-6
MAX_TE_RATIO = .05
MAX_PANELS = 2000
CONDITION_LIMIT = 1e10
RESIDUAL_LIMIT = 1e-9


@dataclass(frozen=True)
class PanelMesh:
    starts: Array
    ends: Array
    centers: Array
    tangents: Array
    normals: Array
    lengths: Array
    orientation: int
    surface_count: int

    @property
    def size(self) -> int:
        return len(self.lengths)


def make_mesh(points: Array, leading_edge_index: int) -> PanelMesh:
    """Preserve traversal and close a finite TE with one explicit base panel.

    Both CCW and CW are supported using signed-area exterior normals. The TE
    declaration is inherited from Airfoil, with additional terminal-panel checks;
    these checks cannot infer the intended TE of an arbitrary polygon.
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or len(p) < 5 or not np.isfinite(p).all():
        raise ValueError("panel contour must contain at least five finite xy points")
    if not 2 <= leading_edge_index <= len(p) - 3:
        raise ValueError("TE-LE-TE ordering needs at least two surface panels on each side")
    polygon = Polygon(p)
    if not polygon.is_valid or polygon.area <= 0:
        raise ValueError("panel contour must be a simple nondegenerate polygon")
    te = p[0] + (p[-1] - p[0]) / 2
    chord = te - p[leading_edge_index]
    chord_length = np.linalg.norm(chord)
    if not np.isfinite(chord_length) or chord_length <= 0:
        raise ValueError("panel chord reference is degenerate")
    cap = np.linalg.norm(p[-1] - p[0])
    if cap / chord_length > MAX_TE_RATIO:
        raise ValueError("finite TE thickness exceeds supported 0.05 element chord")
    # Exact repeated endpoint denotes a sharp TE. A distinct tiny cap is never
    # merged with it: reject unresolved panels rather than modifying points.
    closed = np.vstack((p, p[0])) if cap != 0 else p
    starts, ends = closed[:-1], closed[1:]
    delta = ends - starts
    lengths = np.linalg.norm(delta, axis=1)
    if not np.isfinite(lengths).all() or np.any(lengths / chord_length <= MIN_LENGTH_RATIO):
        raise ValueError("zero-length or unresolved panel in element-chord scale")
    t = delta / lengths[:, None]
    axis = chord / chord_length
    if t[0] @ axis >= -.25 or t[len(p) - 2] @ axis <= .25:
        raise ValueError("invalid TE ordering: terminal surface panels must head upstream/downstream")
    if t[0] @ t[len(p) - 2] >= -.5:
        raise ValueError("unsupported TE: terminal surface tangents must be opposed")
    # Translation before the area calculation avoids cancellation far downstream.
    local = closed - closed[0]
    area2 = np.sum(local[:-1, 0] * local[1:, 1] - local[1:, 0] * local[:-1, 1])
    orientation = 1 if area2 > 0 else -1
    normals = orientation * np.column_stack((t[:, 1], -t[:, 0]))
    return PanelMesh(*(immutable_points(a) for a in
                       (starts, ends, starts + delta / 2, t, normals, lengths)),
                     orientation, len(p) - 1)


def source_velocity(targets: Array, starts: Array, ends: Array) -> Array:
    """Unit source-sheet velocity at off-sheet targets, shape (M,N,2).

    Local x follows each panel; y points left. The angular integral uses atan2,
    avoiding quadrant ambiguity. Endpoints and exact on-sheet targets are
    undefined here; collocation self terms are explicitly inserted elsewhere.
    """
    delta = np.asarray(ends) - np.asarray(starts)
    length = np.linalg.norm(delta, axis=1)
    if np.any(length <= 0) or not np.isfinite(length).all():
        raise ValueError("source panels must have finite positive length")
    t = delta / length[:, None]
    left = np.column_stack((-t[:, 1], t[:, 0]))
    rel = np.asarray(targets)[:, None, :] - starts
    x = np.einsum("ijk,jk->ij", rel, t)
    y = np.einsum("ijk,jk->ij", rel, left)
    r1, r2 = np.hypot(x, y), np.hypot(x - length, y)
    on_sheet = (y == 0) & (x >= 0) & (x <= length)
    if np.any(on_sheet) or np.any(r1 == 0) or np.any(r2 == 0):
        raise ValueError("source influence is undefined on a panel or its endpoints")
    u = (np.log(r1) - np.log(r2)) / (2 * np.pi)
    # Scale before the angle products to avoid overflow/underflow.
    scale = np.maximum(np.maximum(np.abs(x), np.abs(y)), length)
    xs, ys, ls = x / scale, y / scale, length / scale
    v = np.arctan2(ys * ls, xs * (xs - ls) + ys * ys) / (2 * np.pi)
    velocity = u[..., None] * t + v[..., None] * left
    if not np.isfinite(velocity).all():
        raise ValueError("nonfinite analytical panel influence")
    return velocity


def influences(meshes: tuple[PanelMesh, ...]) -> tuple[Array, Array]:
    """Source/vortex vector influences with exterior one-sided self terms.

    Positive vortex strength means counterclockwise point-vortex circulation.
    Its velocity is the source velocity rotated +90 degrees. Therefore source
    self = n/2 and vortex self = orientation*t/2, for either input traversal.
    """
    starts = np.vstack([m.starts for m in meshes])
    ends = np.vstack([m.ends for m in meshes])
    centers = np.vstack([m.centers for m in meshes])
    normals = np.vstack([m.normals for m in meshes])
    # Evaluate off-diagonal entries only: exact on-panel evaluation is undefined.
    size = len(starts)
    source = np.empty((size, size, 2))
    for j in range(size):
        others = np.arange(size) != j
        source[others, j] = source_velocity(centers[others], starts[j:j+1], ends[j:j+1])[:, 0]
        source[j, j] = .5 * normals[j]
    vortex = np.stack((-source[:, :, 1], source[:, :, 0]), axis=2)
    return source, vortex


@dataclass(frozen=True)
class PanelSystem:
    meshes: tuple[PanelMesh, ...]
    matrix: Array
    rhs: Array
    normal: Array
    tangent: Array
    freestream: Array
    kutta_pairs: tuple[tuple[int, int], ...]


def assemble(meshes: tuple[PanelMesh, ...], alpha_deg: float) -> PanelSystem:
    """N source strengths plus one constant vortex-sheet strength per body."""
    if not meshes or sum(m.size for m in meshes) > MAX_PANELS:
        raise ValueError("unsupported panel count (limit 2000)")
    source, vortex = influences(meshes)
    normal = np.vstack([m.normals for m in meshes])
    tangent = np.vstack([m.tangents for m in meshes])
    count, elements = len(normal), len(meshes)
    velocity = np.empty((count, count + elements, 2))
    velocity[:, :count] = source
    pairs = []
    offset = 0
    for e, mesh in enumerate(meshes):
        velocity[:, count + e] = np.sum(vortex[:, offset:offset + mesh.size], axis=1)
        pairs.append((offset, offset + mesh.surface_count - 1))
        offset += mesh.size
    vn = np.einsum("ijk,ik->ij", velocity, normal)
    vt = np.einsum("ijk,ik->ij", velocity, tangent)
    angle = np.deg2rad(alpha_deg)
    free = np.array([np.cos(angle), np.sin(angle)])  # V_inf = 1
    matrix = np.vstack([vn, *[vt[a] + vt[b] for a, b in pairs]])
    rhs = np.concatenate((-normal @ free, [-((tangent[a] + tangent[b]) @ free) for a, b in pairs]))
    if not np.isfinite(matrix).all() or not np.isfinite(rhs).all():
        raise ValueError("nonfinite coupled panel system")
    return PanelSystem(meshes, matrix, rhs, vn, vt, free, tuple(pairs))


def solve(system: PanelSystem, moment_x: float, moment_y: float, sense: str) -> dict:
    """Solve and gate diagnostics before computing any aerodynamic outputs."""
    a, b = system.matrix, system.rhs
    diagnostics = {"unknown_count": len(b), "condition_number_2": None}
    try:
        condition = float(np.linalg.cond(a))
        diagnostics["condition_number_2"] = condition if np.isfinite(condition) else None
        if not np.isfinite(condition) or condition > CONDITION_LIMIT:
            return {"accepted": False, "reason": "ill-conditioned or singular panel system", "diagnostics": diagnostics}
        strengths = np.linalg.solve(a, b)
    except np.linalg.LinAlgError as error:
        return {"accepted": False, "reason": f"linear algebra failure: {error}", "diagnostics": diagnostics}
    if not np.isfinite(strengths).all():
        return {"accepted": False, "reason": "nonfinite panel solution", "diagnostics": diagnostics}
    count = sum(m.size for m in system.meshes)
    n = np.vstack([m.normals for m in system.meshes])
    t = np.vstack([m.tangents for m in system.meshes])
    vn = system.normal @ strengths + n @ system.freestream
    vt = system.tangent @ strengths + t @ system.freestream
    kutta = np.array([vt[i] + vt[j] for i, j in system.kutta_pairs])
    linear = a @ strengths - b
    backward = np.max(np.abs(linear)) / (np.linalg.norm(a, np.inf) * np.linalg.norm(strengths, np.inf)
                                         + np.linalg.norm(b, np.inf))
    diagnostics.update(impermeability_max=float(np.max(np.abs(vn))),
                       kutta_max=float(np.max(np.abs(kutta))),
                       backward_error=float(backward), kutta_residuals=kutta.tolist())
    if (not np.isfinite(vt).all() or not np.isfinite(vn).all()
            or max(diagnostics[key] for key in ("impermeability_max", "kutta_max", "backward_error")) > RESIDUAL_LIMIT):
        return {"accepted": False, "reason": "numerical residual exceeds acceptance policy", "diagnostics": diagnostics}
    cp = 1 - vt**2
    lengths = np.concatenate([m.lengths for m in system.meshes])
    centers = np.vstack([m.centers for m in system.meshes])
    df = -cp[:, None] * n * lengths[:, None]  # force / (q_inf*c_main)
    force = np.sum(df, axis=0)
    lift_axis = np.array([-system.freestream[1], system.freestream[0]])
    cl = float(force @ lift_axis)
    arm = centers - [moment_x, moment_y]
    cm = float(np.sum(arm[:, 0] * df[:, 1] - arm[:, 1] * df[:, 0]))
    if sense == "clockwise":
        cm = -cm
    elements = []
    offset = 0
    for e, mesh in enumerate(system.meshes):
        sl = slice(offset, offset + mesh.size)
        gamma = float(strengths[count + e])
        elements.append({"panel_count": mesh.size, "surface_panel_count": mesh.surface_count,
                         "orientation": "CCW" if mesh.orientation == 1 else "CW",
                         "centers": mesh.centers.tolist(), "normals": mesh.normals.tolist(),
                         "lengths": mesh.lengths.tolist(), "cp": cp[sl].tolist(),
                         "tangential_velocity": vt[sl].tolist(), "normal_velocity": vn[sl].tolist(),
                         "source_strengths": strengths[sl].tolist(), "vortex_sheet_strength": gamma,
                         "circulation": gamma * float(np.sum(mesh.lengths)),
                         "source_flux": float(strengths[sl] @ mesh.lengths)})
        offset += mesh.size
    diagnostics["cl_from_circulation"] = -2 * sum(e["circulation"] for e in elements)
    diagnostics["pressure_circulation_cl_difference"] = cl - diagnostics["cl_from_circulation"]
    # Finite-panel pressure integration has a force component along freestream;
    # expose it only as a discretization diagnostic, never as physical drag.
    diagnostics["streamwise_pressure_force_diagnostic"] = float(force @ system.freestream)
    return {"accepted": True, "reason": "linear system and boundary residuals accepted",
            "cl": cl, "cm": cm, "cd": None, "elements": elements, "diagnostics": diagnostics}
