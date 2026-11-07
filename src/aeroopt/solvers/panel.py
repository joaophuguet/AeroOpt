"""Steady 2D incompressible inviscid Hess--Smith adapter for main + flap."""

import time
from dataclasses import dataclass

import numpy as np
from shapely.geometry import Polygon

from aeroopt.geometry import MultiElementAirfoil

from . import _panel_method as method
from .base import AerodynamicSolver, SolverCase, SolverIdentity
from .conditions import OperatingConditions, real_number
from .execution import AttemptContext, ExecutionEvidence
from .results import AerodynamicCoefficients, ConvergenceStatus, MomentReference, ParsedOutput
from .serialization import FrozenJSON, read_json

MOMENT_NORMALIZATION = "q_inf*c_main^2"
OUTPUT_SCHEMA = "aeroopt.panel-output.v1"


@dataclass(frozen=True)
class PanelSolver(AerodynamicSolver):
    """No physical defaults for moment reference; use the existing SolverRunner.

    Mach must be exactly zero (incompressible model). Reynolds is retained as
    metadata but is not used by this inviscid method. No transition extensions
    are supported. Dense NumPy execution is cooperative, not preemptible.

    Cp, circulation and residual arrays live in panel_output.json, inventoried
    with SHA-256 in SolverResult. No change to the solver-neutral result schema.
    """

    moment_reference: MomentReference

    def __post_init__(self) -> None:
        if not isinstance(self.moment_reference, MomentReference):
            raise TypeError("PanelSolver requires an explicit MomentReference")
        if self.moment_reference.normalization != MOMENT_NORMALIZATION:
            raise ValueError(f"PanelSolver moment normalization must be {MOMENT_NORMALIZATION!r}")

    @property
    def identity(self) -> SolverIdentity:
        return SolverIdentity("aeroopt.hess_smith", "0.3.0", False, True,
            ("steady", "2d", "incompressible", "inviscid", "main_flap_interaction",
             "surface_cp", "cl", "cm", "element_circulation", "boundary_residuals"),
            FrozenJSON({"numpy_version": np.__version__, "formulation": "constant_sources_constant_vortex_per_element",
                        "moment_reference": self.moment_reference.to_dict(), "freestream_speed": 1.0,
                        "force_normalization": "q_inf*c_main", "circulation_positive": "counterclockwise",
                        "te_treatment": "straight_cap_and_opposed_terminal_tangential_velocities",
                        "acceptance_scope": "linear_algebra_and_collocation_only",
                        "condition_limit_2": method.CONDITION_LIMIT, "residual_limit": method.RESIDUAL_LIMIT,
                        "minimum_clearance_c_main": method.MIN_CLEARANCE,
                        "minimum_panel_length_over_element_chord": method.MIN_LENGTH_RATIO,
                        "maximum_te_thickness_over_element_chord": method.MAX_TE_RATIO,
                        "maximum_panels": method.MAX_PANELS}))

    @staticmethod
    def _meshes(main, flap, geometry):
        return (method.make_mesh(main, geometry.main.leading_edge_index),
                method.make_mesh(flap, geometry.flap.leading_edge_index))

    def validate_input(self, geometry: MultiElementAirfoil,
                       conditions: OperatingConditions) -> tuple[str, ...]:
        reasons = []
        if conditions.mach != 0:
            reasons.append("PanelSolver requires mach=0 for the incompressible model")
        if conditions.additional:
            reasons.append("PanelSolver does not support additional/transition/turbulence parameters")
        try:
            main, flap = geometry.main_points, geometry.flap_points
            meshes = self._meshes(main, flap, geometry)
            if sum(m.size for m in meshes) > method.MAX_PANELS:
                reasons.append("PanelSolver supports at most 2000 panels")
            distance = Polygon(main).boundary.distance(Polygon(flap).boundary)
            if distance < method.MIN_CLEARANCE:
                reasons.append("panel numerical clearance is below 1e-6 main chord; no geometry repair performed")
        except ValueError as error:
            reasons.append(str(error))
        return tuple(reasons)

    def prepare(self, case: SolverCase, context: AttemptContext) -> None:
        meshes = self._meshes(case.main_points, case.flap_points, case.geometry)
        context.write_json("panel_request.json", {"case_id": case.case_id,
            "schema_version": OUTPUT_SCHEMA, "solver": self.identity.to_dict(),
            "elements": [{"name": name, "panel_count": mesh.size, "surface_panel_count": mesh.surface_count,
                          "starts": mesh.starts.tolist(), "ends": mesh.ends.tolist()}
                         for name, mesh in zip(("main", "flap"), meshes, strict=True)]})
        context.write_text("panel.log", "Hess-Smith: steady 2D incompressible inviscid.\n"
                           "Reynolds retained as metadata; no viscous physics. cd unavailable.\n"
                           "Numerical acceptance is not experimental aerodynamic validation.\n")

    def execute(self, case: SolverCase, context: AttemptContext) -> ExecutionEvidence:
        start = time.monotonic()

        def check_deadline():
            if context.deadline is not None and time.monotonic() >= context.deadline:
                raise TimeoutError("panel library execution exceeded cooperative deadline")

        check_deadline()
        meshes = self._meshes(case.main_points, case.flap_points, case.geometry)
        system = method.assemble(meshes, case.conditions.alpha_deg)
        check_deadline()
        ref = self.moment_reference
        output = method.solve(system, ref.x_over_c, ref.y_over_c, ref.positive_sense)
        distance = Polygon(case.main_points).boundary.distance(Polygon(case.flap_points).boundary)
        output["diagnostics"]["clearance_c_main"] = float(distance)
        output["diagnostics"]["clearance_over_largest_panel"] = float(distance / max(np.max(m.lengths) for m in meshes))
        context.write_json("panel_output.json", {"schema_version": OUTPUT_SCHEMA,
            "case_id": case.case_id, "moment_reference": ref.to_dict(), **output})
        check_deadline()  # A blocking BLAS call cannot be interrupted in process.
        context.event("execute", "panel numerical solution completed", accepted=output["accepted"],
                      reason=output["reason"])
        return ExecutionEvidence(None, time.monotonic() - start)

    def parse(self, case: SolverCase, context: AttemptContext,
              evidence: ExecutionEvidence) -> ParsedOutput:
        output = read_json(context.path("panel_output.json"))
        if (output["schema_version"] != OUTPUT_SCHEMA or output["case_id"] != case.case_id
                or output["moment_reference"] != self.moment_reference.to_dict()
                or type(output["accepted"]) is not bool):
            raise ValueError("panel output schema/case/reference/acceptance mismatch")
        if not output["accepted"]:
            return ParsedOutput(ConvergenceStatus.NOT_CONVERGED, messages=(output["reason"],))
        diagnostics = output["diagnostics"]
        condition = real_number(diagnostics["condition_number_2"], "condition_number_2")
        if not 1 <= condition <= method.CONDITION_LIMIT:
            raise ValueError("panel condition number failed acceptance")
        for key in ("impermeability_max", "kutta_max", "backward_error"):
            if not 0 <= real_number(diagnostics[key], key) <= method.RESIDUAL_LIMIT:
                raise ValueError(f"panel {key} failed acceptance")
        meshes = self._meshes(case.main_points, case.flap_points, case.geometry)
        if len(output["elements"]) != 2 or diagnostics["unknown_count"] != sum(m.size for m in meshes) + 2:
            raise ValueError("panel output element/system size mismatch")
        force = np.zeros(2)
        moment_ccw = 0.
        ref = self.moment_reference
        for mesh, element in zip(meshes, output["elements"], strict=True):
            if element["panel_count"] != mesh.size or element["surface_panel_count"] != mesh.surface_count:
                raise ValueError("panel output mesh size mismatch")
            for key in ("cp", "normal_velocity", "tangential_velocity", "source_strengths", "lengths"):
                values = np.asarray(element[key], dtype=float)
                if values.shape != (mesh.size,) or not np.isfinite(values).all():
                    raise ValueError(f"invalid panel {key} distribution")
            vt, vn = np.array(element["tangential_velocity"]), np.array(element["normal_velocity"])
            if (np.max(np.abs(vn)) > method.RESIDUAL_LIMIT
                    or abs(vt[0] + vt[mesh.surface_count - 1]) > method.RESIDUAL_LIMIT
                    or not np.allclose(element["cp"], 1 - vt**2, rtol=1e-13, atol=1e-13)):
                raise ValueError("panel surface data contradicts acceptance/Cp")
            if (not np.array_equal(element["lengths"], mesh.lengths)
                    or not np.array_equal(element["centers"], mesh.centers)
                    or not np.array_equal(element["normals"], mesh.normals)):
                raise ValueError("panel output geometry differs from input mesh")
            gamma = real_number(element["vortex_sheet_strength"], "vortex_sheet_strength")
            circulation = real_number(element["circulation"], "circulation")
            if not np.isclose(circulation, gamma * np.sum(mesh.lengths), rtol=1e-13, atol=1e-13):
                raise ValueError("inconsistent element circulation")
            df = -np.array(element["cp"])[:,None] * mesh.normals * mesh.lengths[:,None]
            force += np.sum(df, axis=0)
            arm = mesh.centers - [ref.x_over_c, ref.y_over_c]
            moment_ccw += float(np.sum(arm[:,0] * df[:,1] - arm[:,1] * df[:,0]))
        if output["cd"] is not None:
            raise ValueError("PanelSolver cannot publish physical drag")
        coefficients = AerodynamicCoefficients(output["cl"], None, output["cm"], self.moment_reference)
        angle = np.deg2rad(case.conditions.alpha_deg)
        cl = float(force @ [-np.sin(angle), np.cos(angle)])
        cm = moment_ccw if ref.positive_sense == "counterclockwise" else -moment_ccw
        if not np.allclose([coefficients.cl, coefficients.cm], [cl, cm], rtol=1e-12, atol=1e-12):
            raise ValueError("panel forces/moment disagree with surface pressure")
        return ParsedOutput(ConvergenceStatus.CONVERGED, coefficients,
            ("Inviscid solution: Re is metadata; cd unavailable; no separation or stall prediction.",
             "Cp, CCW circulation and numerical residuals: panel_output.json. Refine panels before interpretation."))
