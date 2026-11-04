"""Attempt-local artifacts and controlled sequential external processes."""

import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from ._paths import native_path, reject_links, relative_path, safe_path
from .conditions import real_number
from .serialization import canonical_json, write_json


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class ExecutionEvidence:
    """Process return code is separate from parsed scientific convergence."""

    return_code: int | None
    duration_s: float
    timed_out: bool = False
    command: tuple[str, ...] = ()
    stdout_path: str | None = None
    stderr_path: str | None = None
    tree_cleanup: str = "not_applicable_python_library"
    supervisor_pid: int | None = None

    def __post_init__(self) -> None:
        if self.return_code is not None and type(self.return_code) is not int:
            raise TypeError("return_code must be int or absent")
        duration = real_number(self.duration_s, "execution duration")
        if duration < 0 or type(self.timed_out) is not bool:
            raise ValueError("invalid execution timing")
        object.__setattr__(self, "duration_s", duration)
        if isinstance(self.command, (str, bytes)):
            raise TypeError("command must be an argument sequence, not a shell string")
        command = tuple(self.command)
        if not all(isinstance(item, str) for item in command):
            raise TypeError("command must contain text arguments")
        object.__setattr__(self, "command", command)
        if self.supervisor_pid is not None and (type(self.supervisor_pid) is not int or self.supervisor_pid <= 0):
            raise ValueError("supervisor_pid must be a positive integer when present")
        for name in ("stdout_path", "stderr_path"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, relative_path(getattr(self, name)))

    def to_dict(self) -> dict:
        return {"return_code": self.return_code, "duration_s": self.duration_s,
                "timed_out": self.timed_out, "command": list(self.command),
                "stdout_path": self.stdout_path, "stderr_path": self.stderr_path,
                "tree_cleanup": self.tree_cleanup, "supervisor_pid": self.supervisor_pid}


@dataclass(frozen=True)
class AttemptContext:
    """Trusted adapter helpers. This path boundary is not an OS sandbox."""

    root: Path
    timeout_s: float | None
    deadline: float | None = None

    def __post_init__(self) -> None:
        root = native_path(self.root)
        reject_links(root)
        object.__setattr__(self, "root", root)
        if self.timeout_s is not None:
            value = real_number(self.timeout_s, "timeout_s")
            if value <= 0:
                raise ValueError("timeout_s must be positive")
            object.__setattr__(self, "timeout_s", value)
        if self.deadline is not None:
            object.__setattr__(self, "deadline", real_number(self.deadline, "deadline"))

    def path(self, relative: str | Path) -> Path:
        return safe_path(self.root, relative)

    def write_json(self, relative: str | Path, value: object) -> Path:
        path = self.path(relative)
        path.parent.mkdir(parents=True, exist_ok=True)
        write_json(path, value)
        return path

    def write_text(self, relative: str | Path, text: str) -> Path:
        path = self.path(relative)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(text)
        return path

    def event(self, stage: str, message: str, **data: object) -> None:
        record = {"at_utc": utc_now(), "stage": stage, "message": message, **data}
        with self.path("events.jsonl").open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(canonical_json(record) + "\n")

    def run_process(self, argv: list[str] | tuple[str, ...], *, label: str = "process",
                    timeout_s: float | None = None) -> ExecutionEvidence:
        """Run argv without a shell, capture files, terminate the ordinary tree.

        Windows: a gated Python worker joins a kill-on-close job before launching
        the executable. No fallback to an uncontained process if assignment fails.
        POSIX: a fresh process group. Descendants deliberately using OS escape
        mechanisms are outside the trusted-adapter contract.
        """
        if not isinstance(argv, (list, tuple)) or not argv:
            raise TypeError("argv must be a nonempty argument sequence")
        if not all(isinstance(item, str) and "\0" not in item for item in argv):
            raise ValueError("arguments must be NUL-free strings")
        executable = Path(argv[0])
        if not executable.is_absolute() or not executable.is_file():
            raise ValueError("executable must be an explicit existing absolute path")
        if os.name == "nt" and executable.suffix.lower() != ".exe":
            raise ValueError("Windows execution requires an .exe, not an implicit shell script")
        timeout = self.timeout_s if timeout_s is None else real_number(timeout_s, "timeout_s")
        if timeout is not None and timeout <= 0:
            raise ValueError("timeout_s must be positive")
        if self.deadline is not None:
            remaining = self.deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("execution deadline already expired")
            timeout = remaining if timeout is None else min(timeout, remaining)
        stdout = self.path(f"{label}.stdout.log")
        stderr = self.path(f"{label}.stderr.log")
        stdout.parent.mkdir(parents=True, exist_ok=True)
        stderr.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path("tmp")
        temporary.mkdir(exist_ok=True)
        environment = os.environ.copy()
        environment.update(TEMP=str(temporary), TMP=str(temporary), TMPDIR=str(temporary),
                           MPLCONFIGDIR=str(temporary / "matplotlib"))
        started = time.monotonic()
        self.event("execute", "external process requested", argv=list(argv), timeout_s=timeout)
        process = None
        job = None
        timed_out = False
        cleanup = "posix_process_group_terminated"
        with stdout.open("xb") as out, stderr.open("xb") as err:
            try:
                if os.name == "nt":
                    from ._windows_job import WindowsJob

                    job = WindowsJob()
                    worker = str(Path(__file__).with_name("_process_worker.py"))
                    process = subprocess.Popen([sys.executable, "-I", "-u", worker, *argv],
                        cwd=self.root, stdin=subprocess.PIPE, stdout=out, stderr=err,
                        env=environment, shell=False, creationflags=subprocess.CREATE_NO_WINDOW)
                    job.assign(process.pid)
                    process.stdin.write(b"G")
                    process.stdin.flush()
                    process.stdin.close()
                    cleanup = "windows_job_terminated_and_empty"
                else:
                    process = subprocess.Popen(argv, cwd=self.root, stdin=subprocess.DEVNULL,
                        stdout=out, stderr=err, env=environment, shell=False, start_new_session=True)
                try:
                    elapsed = time.monotonic() - started
                    process.wait(timeout=None if timeout is None else max(.001, timeout - elapsed))
                except subprocess.TimeoutExpired:
                    timed_out = True
            finally:
                try:
                    if job is not None:
                        job.terminate_and_wait()
                    elif process is not None:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                finally:
                    if job is not None:
                        job.close()
                    if process is not None:
                        if process.poll() is None:
                            process.kill()  # Also handles a gated worker whose job assignment failed.
                        process.wait(timeout=5.)
                        if process.stdin is not None and not process.stdin.closed:
                            process.stdin.close()
        evidence = ExecutionEvidence(process.returncode, time.monotonic() - started, timed_out,
                                     tuple(argv), stdout.relative_to(self.root).as_posix(),
                                     stderr.relative_to(self.root).as_posix(), cleanup, process.pid)
        self.write_json(f"{label}.execution.json", evidence.to_dict())
        self.event("execute", "external process completed", **evidence.to_dict())
        return evidence
