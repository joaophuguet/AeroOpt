"""Portable relative artifact paths; reject traversal, devices, ADS and links."""

import os
from pathlib import Path, PureWindowsPath


def native_path(path: str | Path) -> Path:
    """Use explicit Win32 long-path syntax, independent of machine policy."""
    result = Path(path).absolute()
    text = str(result)
    if os.name == "nt" and not text.startswith("\\\\?\\"):
        text = "\\\\?\\UNC\\" + text[2:] if text.startswith("\\\\") else "\\\\?\\" + text
        result = Path(text)
    return result


def relative_path(value: str | Path) -> str:
    text = str(value).replace("\\", "/")
    windows = PureWindowsPath(text)
    parts = text.split("/")
    if (windows.drive or text.startswith("/") or not text
            or any(part in ("", ".", "..") for part in parts)):
        raise ValueError("expected a relative path without traversal")
    reserved = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                *(f"LPT{i}" for i in range(1, 10))}
    for part in parts:
        if (any(ord(char) < 32 or char in ':<>"|?*' for char in part)
                or part.endswith((".", " ")) or part.split(".")[0].upper() in reserved):
            raise ValueError("unsafe or reserved Windows artifact name")
    return "/".join(parts)


def reject_links(path: Path) -> None:
    for item in (path, *path.parents):
        if item.is_symlink() or item.is_junction():
            raise ValueError(f"symlink/junction is not allowed in execution paths: {item}")


def safe_path(root: Path, value: str | Path) -> Path:
    root = native_path(root)
    relative = relative_path(value)
    target = root.joinpath(*relative.split("/"))
    reject_links(target)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError("artifact path escapes execution workspace")
    return target
