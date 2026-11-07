"""Optional, commit-pinned Alula JavaScript adapter with supervised Node execution."""

from dataclasses import asdict, dataclass, field
import hashlib
import os
from pathlib import Path
import re
import shutil
from uuid import uuid4

import numpy as np

from .base import AerodynamicSolver, SolverIdentity
from .conditions import real_number
from .execution import AttemptContext
from .results import AerodynamicCoefficients, ConvergenceStatus, MomentReference, ParsedOutput
from .serialization import FrozenJSON, read_json

RESOURCES = Path(__file__).with_name("_alula")
LOCK = read_json(RESOURCES/"upstream-lock.json")
UPSTREAM_COMMIT = LOCK["commit"]
BRIDGE_VERSION = "0.4.0"


def file_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def node_path(path):
    """Node CLI/module URLs require ordinary absolute Windows/UNC spelling."""
    value = str(Path(path).absolute())
    if value.startswith("\\\\?\\UNC\\"):
        return "\\\\"+value[8:]
    return value[4:] if value.startswith("\\\\?\\") else value


def verify_installation(directory):
    """Check exact selected source bytes, including package type and notices."""
    root = Path(directory).resolve()
    for name,expected in LOCK["files"].items():
        path = root/name
        if not path.is_file():
            raise FileNotFoundError(f"Missing pinned Alula file: {path}; run scripts/install_alula.py explicitly")
        if file_sha(path) != expected:
            raise ValueError(f"Pinned Alula source hash mismatch: {name}")
    return root


@dataclass(frozen=True)
class AlulaViscousSettings:
    """Explicit transition physics; numerical defaults are upstream controls."""

    ncrit: float
    trips: tuple[float,float]
    max_iterations: int = 60
    wake_length: float = 2.
    wake_count: int = 24
    tolerance: float = 1e-8
    wake_tolerance: float = 1e-6
    max_wake_iterations: int = 40
    wake_relaxation: float = .8

    def __post_init__(self):
        for name in ("ncrit","wake_length","tolerance","wake_tolerance","wake_relaxation"):
            object.__setattr__(self,name,real_number(getattr(self,name),name))
        if not 0 <= self.ncrit <= 14:
            raise ValueError("Ncrit must be between 0 and 14")
        trips = tuple(self.trips)
        if len(trips) != 2:
            raise ValueError("Two explicit upper/lower trips required")
        trips = tuple(real_number(v,"trip") for v in trips)
        if not all(0 < v <= 1 for v in trips):
            raise ValueError("Trips must lie in (0,1] of each element chord")
        object.__setattr__(self,"trips",trips)
        for name,low,high in (("max_iterations",1,200),("wake_count",4,200),("max_wake_iterations",0,200)):
            value = getattr(self,name)
            if type(value) is not int or not low <= value <= high:
                raise ValueError(f"{name} must be an integer in [{low},{high}]")
        if min(self.wake_length,self.tolerance,self.wake_tolerance) <= 0 or not 0 < self.wake_relaxation <= 1:
            raise ValueError("Positive lengths/tolerances and wake relaxation in (0,1] required")


@dataclass(frozen=True)
class AlulaSolver(AerodynamicSolver):
    """Two immutable SolverCase elements; the standalone bridge supports 1--6.

    Dependency installation is always explicit. Construction performs only a
    supervised Node version probe and local source/hash checks. Paths and UUIDs
    are operational records, never part of the scientific identity.
    """

    mode: str
    moment_reference: MomentReference
    node_executable: str | Path | None = None
    alula_directory: str | Path | None = None
    viscous_settings: AlulaViscousSettings | None = None
    boundary_condition: str | None = None
    _identity: SolverIdentity = field(init=False,repr=False)
    _node_hash: str = field(init=False,repr=False)
    _bridge_hash: str = field(init=False,repr=False)

    def __post_init__(self):
        if self.mode not in ("inviscid","viscous_incompressible"):
            raise ValueError("Alula mode must be inviscid or viscous_incompressible")
        if not isinstance(self.moment_reference,MomentReference) or self.moment_reference.normalization != "q_inf*c_main^2":
            raise ValueError("An explicit MomentReference normalized by q_inf*c_main^2 is required")
        viscous = self.mode == "viscous_incompressible"
        if viscous != isinstance(self.viscous_settings,AlulaViscousSettings):
            raise ValueError("Viscous mode requires AlulaViscousSettings; inviscid mode must omit them")
        boundary = self.boundary_condition or ("streamfunction" if viscous else "normal-velocity")
        if boundary not in ("normal-velocity","streamfunction") or viscous and boundary != "streamfunction":
            raise ValueError("Unsupported outer boundary condition for selected mode")
        object.__setattr__(self,"boundary_condition",boundary)
        upstream = verify_installation(self.alula_directory or Path.cwd()/".deps/alula"/UPSTREAM_COMMIT)
        node = self.node_executable or shutil.which("node")
        if node is None:
            local = Path.cwd()/".deps"/f"node-{LOCK['tested_node']['version']}-win-x64/node.exe"
            node = local if local.is_file() else None
        if node is None or not Path(node).is_file():
            raise FileNotFoundError("Node.js missing: supply node_executable or install Node >=22 explicitly")
        node = Path(node).resolve()
        if os.name == "nt" and node.suffix.lower() != ".exe":
            raise ValueError("Windows Node must resolve to an absolute .exe path")
        probe = Path.cwd()/".tmp/alula-probes"/uuid4().hex
        probe.mkdir(parents=True,exist_ok=False)
        context = AttemptContext(probe,10)
        evidence = context.run_process([str(node),"--version"],label="version")
        version = (probe/"version.stdout.log").read_text(encoding="utf-8").strip()
        if evidence.return_code != 0 or evidence.timed_out or not re.fullmatch(r"v\d+\.\d+\.\d+",version) or int(version[1:].split('.')[0]) < 22:
            raise ValueError("Node version probe failed; Node >=22 is required")
        object.__setattr__(self,"node_executable",node)
        object.__setattr__(self,"alula_directory",upstream)
        object.__setattr__(self,"_node_hash",file_sha(node))
        object.__setattr__(self,"_bridge_hash",file_sha(RESOURCES/"bridge.mjs"))
        settings = {"upstream_commit":UPSTREAM_COMMIT,"upstream_lock_sha256":file_sha(RESOURCES/"upstream-lock.json"),
                    "bridge_version":BRIDGE_VERSION,"bridge_sha256":self._bridge_hash,
                    "node_version":version,"node_binary_sha256":self._node_hash,
                    "numerical_dependencies":LOCK["runtime_dependencies"],"mode":self.mode,
                    "boundary_condition":boundary,"moment_reference":self.moment_reference.to_dict(),
                    "reference_chord":1.,"alpha_convention":"global_freestream_cos_sin",
                    "upstream_cm_positive":"clockwise","eligibility_policy":"converged_complete_nonprovisional_v1",
                    "geometry_policy":"CCW_only_exact_points_explicit_straight_TE_cap",
                    "viscous":None if not viscous else asdict(self.viscous_settings),
                    "initialization":"auto_multielement_native_single","jacobian_mode":"analytic"}
        object.__setattr__(self,"_identity",SolverIdentity("aeroopt.alula",BRIDGE_VERSION,False,True,
            ("steady","2d","incompressible","main_flap_interaction","surface_cp","cl","cm",self.mode,
             *(('cd','integral_boundary_layer') if viscous else ())),FrozenJSON(settings)))

    @property
    def identity(self):
        return self._identity

    def validate_input(self,geometry,conditions):
        reasons = []
        viscous = self.mode == "viscous_incompressible"
        if conditions.mach != 0:
            reasons.append("Alula incompressible modes require mach=0; nonzero Mach is not reinterpreted")
        if abs(conditions.alpha_deg) > (20 if viscous else 90):
            reasons.append("Incidence outside selected Alula mode range")
        if viscous and not 5e4 <= conditions.reynolds <= 5e7:
            reasons.append("Viscous Reynolds must be in [50000,50000000], based on main chord")
        if conditions.additional:
            reasons.append("Additional conditions are unsupported; declare transition via AlulaViscousSettings")
        total = 0
        for name,points in (("main",geometry.main_points),("flap",geometry.flap_points)):
            closed = np.array_equal(points[0],points[-1])
            count = len(points)-1
            total += count
            p = points-points[0]
            area = np.sum(p[:,0]*np.roll(p[:,1],-1)-np.roll(p[:,0],-1)*p[:,1])
            if area <= 0:
                reasons.append(f"{name}: counterclockwise contour required; no automatic reordering")
            if len(points) < 9 or count+(not closed) > 600:
                reasons.append(f"{name}: requires >=9 supplied points and <=600 solid panels")
            if viscous and count < 40:
                reasons.append(f"{name}: viscous mode requires >=40 wetted panels")
            if not closed and self.boundary_condition != "streamfunction":
                reasons.append(f"{name}: finite base requires explicit streamfunction boundary_condition")
        if total > 700:
            reasons.append("Alula supports at most 700 total wetted panels")
        return tuple(reasons)

    def input_payload(self,case):
        """Exact snapshot translation, also rechecked after the external solve."""
        elements = []
        for role,points in (("main",case.main_points),("flap",case.flap_points)):
            elements.append({"name":role,"role":role,"points":[{"x":float(x),"y":float(y)} for x,y in points],
                             "closure":"sharp_repeated_endpoint" if np.array_equal(points[0],points[-1]) else "straight_TE_cap"})
        settings = self.identity.settings.to_dict()
        return {"schema_version":"aeroopt.alula-input.v1","case_id":case.case_id,
            "mode":self.mode,"elements":elements,
            "conditions":{k:getattr(case.conditions,k) for k in ("alpha_deg","reynolds","mach")},
            "moment_reference":self.moment_reference.to_dict(),
            "settings":{"boundary_condition":self.boundary_condition,"viscous":settings["viscous"]},
            "runtime":{"node_version":settings["node_version"],"upstream_commit":UPSTREAM_COMMIT}}

    def prepare(self,case,context):
        if os.name == "nt" and len(node_path(context.root).encode("utf-16-le"))//2 > 240:
            raise ValueError("Alula execution directory exceeds the supported Windows 240-character cwd budget; choose a shorter SolverRunner directory")
        verify_installation(self.alula_directory)
        if file_sha(self.node_executable) != self._node_hash or file_sha(RESOURCES/"bridge.mjs") != self._bridge_hash:
            raise ValueError("Node or bridge changed since identity was established")
        context.write_json("alula_input.json",self.input_payload(case))
        for name in ("bridge.mjs","upstream-lock.json","COPYING"):
            with context.path(name).open("xb") as stream:
                stream.write((RESOURCES/name).read_bytes())
        context.write_json("alula_runtime.json",{"node_executable_display_only":str(self.node_executable),
            "upstream_directory_display_only":str(self.alula_directory),"source_commit":UPSTREAM_COMMIT,
            "node_binary_sha256":self._node_hash})

    def execute(self,case,context):
        # Absolute ordinary spelling avoids Node main-module realpath failures
        # with a Win32 \\?\ cwd/entry. Python retains native paths for its IO.
        return context.run_process([str(self.node_executable),node_path(context.path("bridge.mjs")),
            "--upstream",node_path(self.alula_directory),"--input",node_path(context.path("alula_input.json")),
            "--output",node_path(context.path("alula_output.json"))],label="alula")

    def parse(self,case,context,evidence):
        output = read_json(context.path("alula_output.json"))
        requested = read_json(context.path("alula_input.json"))
        if requested != self.input_payload(case):
            raise ValueError("Alula input changed from immutable SolverCase snapshot")
        if (output["schema_version"] != "aeroopt.alula-output.v1" or output["case_id"] != case.case_id
                or output["mode"] != self.mode or output["moment_reference"] != self.moment_reference.to_dict()
                or output["conditions"] != requested["conditions"] or output["geometry"] != requested["elements"]
                or output["input_file_sha256"] != file_sha(context.path("alula_input.json"))):
            raise ValueError("Alula output schema/input/case/geometry/conventions mismatch")
        settings = self.identity.settings.to_dict()
        metadata = output["metadata"]
        for key,expected in (("upstream_commit",UPSTREAM_COMMIT),("bridge_version",BRIDGE_VERSION),
                             ("bridge_sha256",self._bridge_hash),("lock_sha256",settings["upstream_lock_sha256"]),
                             ("node_version",settings["node_version"]),("reference_chord",1)):
            if metadata[key] != expected:
                raise ValueError(f"Alula runtime/output identity mismatch: {key}")
        if metadata["physical_validation"] is not False or metadata["numerical_eligibility_only"] is not True:
            raise ValueError("Alula output omitted its scientific scope")
        expected_limits = ["RESEARCH_MODEL_NOT_PHYSICALLY_VALIDATED"]
        if any(e["closure"] == "straight_TE_cap" for e in requested["elements"]):
            expected_limits.append("ALULA_FINITE_BASE_RESEARCH")
        if output["limitations"] != expected_limits:
            raise ValueError("Alula output changed its mandatory limitation labels")
        if metadata["numerical_dependencies"] != LOCK["runtime_dependencies"]:
            raise ValueError("Alula numerical dependency metadata mismatch")
        status = output["status"]
        if type(status["eligible"]) is not bool:
            raise ValueError("Invalid Alula eligibility flag")
        coefficients = output["coefficients"]
        if status["upstream"] != output["upstream"]["status"]:
            raise ValueError("Contradictory upstream/bridge status")
        if not isinstance(output["warnings"],list) or not all(isinstance(w,str) for w in output["warnings"]):
            raise ValueError("Invalid Alula warnings")
        if set(coefficients) != {"cl","cd","cm"}:
            raise ValueError("Unexpected coefficient fields")
        for key,value in coefficients.items():
            if value is not None:
                real_number(value,key)
        messages = tuple(output["warnings"])+tuple(output["limitations"])+("Cp, BL and upstream diagnostics: alula_output.json; numerical eligibility is not physical validation.",)
        if not status["eligible"]:
            if any(value is not None for value in coefficients.values()):
                raise ValueError("Ineligible Alula coefficients must be withheld")
            return ParsedOutput(ConvergenceStatus.NOT_CONVERGED,messages=messages)
        if (status["upstream"] != "solved" or status["numerical_convergence"] != "converged"
                or status["coefficient_status"] != "available" or output["nonfinite_diagnostics"]
                or status.get("model_admissible") is not True
                or output["upstream"].get("coefficientStatus") not in (None,"available","converged")):
            raise ValueError("Alula eligibility contradicts convergence/availability")
        d = output["upstream"]["diagnostics"]
        gates = (("linearResidual",1e-9),("kuttaResidual",1e-9)) if self.mode == "inviscid" else (
            ("equationResidual",self.viscous_settings.tolerance),("wakeResidual",self.viscous_settings.wake_tolerance),
            ("maxRelativeDisplacement",.1))
        if self.mode == "inviscid":
            extras = ("normalVelocityResidual",) if self.boundary_condition == "normal-velocity" else (
                "nodalStreamfunctionResidual","interiorTEResidual")
            gates += tuple((key,1e-9) for key in extras)
        for key,limit in gates:
            if not 0 <= real_number(d[key],key) <= limit:
                raise ValueError(f"Alula diagnostic gate failed: {key}")
        required = ("cl","cm") if self.mode == "inviscid" else ("cl","cd","cm")
        if any(coefficients[k] is None for k in required) or self.mode == "inviscid" and coefficients["cd"] is not None:
            raise ValueError("Missing or physically unavailable Alula coefficients")
        raw = output["upstream"]["coefficients"]
        sign = -1 if self.moment_reference.positive_sense == "counterclockwise" else 1
        if coefficients["cl"] != raw["cl"] or coefficients["cm"] != sign*real_number(raw["cm_clockwise"],"upstream cm") or coefficients["cd"] != raw["cd"]:
            raise ValueError("Alula coefficient sign/translation mismatch")
        elements = output["upstream"]["elements"]
        if len(elements) != 2 or any(not e["cp"] for e in elements):
            raise ValueError("Missing Alula surface pressure distributions")
        for e,supplied in zip(elements,requested["elements"],strict=True):
            expected = supplied["points"]
            if supplied["closure"] == "straight_TE_cap":
                expected = expected+[expected[0]]
            if e["points"] != expected:
                raise ValueError("Alula returned solid geometry differs from supplied coordinates")
            for p in e["cp"]:
                for key in ("x","y","cp"):
                    real_number(p[key],f"surface {key}")
        return ParsedOutput(ConvergenceStatus.CONVERGED,
            AerodynamicCoefficients(**coefficients,moment_reference=self.moment_reference),messages)
