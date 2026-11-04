"""Harmless local process tree for execution-infrastructure tests only."""

import os
import subprocess
import sys
import time
from pathlib import Path


def main():
    mode = sys.argv[1]
    print("NON-PHYSICAL PROCESS FIXTURE STDOUT", flush=True)
    print("NON-PHYSICAL PROCESS FIXTURE STDERR", file=sys.stderr, flush=True)
    if mode == "success":
        return
    if mode == "exit7":
        raise SystemExit(7)
    if mode == "child":
        Path("child.pid").write_text(str(os.getpid()), encoding="ascii")
        time.sleep(3600)
        return
    if mode in ("blocking_tree", "orphan_on_success"):
        Path("parent.pid").write_text(str(os.getpid()), encoding="ascii")
        subprocess.Popen([sys.executable, "-u", str(Path(__file__).resolve()), "child"],
                         stdin=subprocess.DEVNULL, shell=False,
                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        deadline = time.monotonic() + 5
        while not Path("child.pid").exists():
            if time.monotonic() > deadline:
                raise RuntimeError("fixture child did not start")
            time.sleep(.01)
        print("CHILD_READY", flush=True)
        if mode == "blocking_tree":
            time.sleep(3600)
        return
    if mode == "mutate_input":
        Path("geometry/main.dat").write_text("corrupted by non-physical test fixture", encoding="utf-8")
        return
    raise ValueError(mode)


if __name__ == "__main__":
    main()
