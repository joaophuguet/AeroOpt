"""Trusted Windows launch gate. No solver code runs until job assignment succeeds."""

import subprocess
import sys


def main() -> None:
    if sys.stdin.buffer.read(1) != b"G":
        raise SystemExit(125)
    completed = subprocess.run(sys.argv[1:], stdin=subprocess.DEVNULL, shell=False, check=False,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
