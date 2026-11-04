"""Fixtures for contract checks only: none implements aerodynamic physics."""

import ctypes
import os
import sys
from ctypes import wintypes
from pathlib import Path

from aeroopt.solvers import (AerodynamicSolver, ConvergenceStatus, ExecutionEvidence,
                            FrozenJSON, ParsedOutput, SolverIdentity)


class ProcessFixtureSolver(AerodynamicSolver):
    """Non-physical external program fixture, never a physical solver."""

    def __init__(self, mode: str = "success") -> None:
        self.mode = mode

    @property
    def identity(self):
        return SolverIdentity("non_physical.process_fixture", "1", True, False,
                              ("test_only",), FrozenJSON({"mode": self.mode}))

    def prepare(self, case, context):
        context.write_text("test_evidence.txt", "Non-physical subprocess fixture\n")

    def execute(self, case, context):
        fixture = Path(__file__).with_name("fixtures") / "process_fixture.py"
        return context.run_process([sys.executable, "-u", str(fixture.resolve()), self.mode])

    def parse(self, case, context, evidence):
        return ParsedOutput(ConvergenceStatus.NOT_APPLICABLE)


def process_running(pid: int) -> bool:
    if os.name != "nt":
        try:
            os.kill(pid, 0)
            return True
        except ProcessLookupError:
            return False
    api = ctypes.WinDLL("kernel32", use_last_error=True)
    api.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    api.OpenProcess.restype = wintypes.HANDLE
    api.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    api.WaitForSingleObject.restype = wintypes.DWORD
    api.CloseHandle.argtypes = [wintypes.HANDLE]
    api.CloseHandle.restype = wintypes.BOOL
    handle = api.OpenProcess(0x00100000 | 0x1000, False, pid)  # SYNCHRONIZE | QUERY_LIMITED_INFORMATION
    if not handle:
        if ctypes.get_last_error() == 87:  # no such PID
            return False
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        return api.WaitForSingleObject(handle, 0) == 258  # WAIT_TIMEOUT
    finally:
        api.CloseHandle(handle)
