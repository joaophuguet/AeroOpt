"""Finite JSON, immutable metadata and deterministic float64 identity encoding."""

import hashlib
import json
import math
import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from uuid import uuid4

from ._paths import native_path

INPUT_SCHEMA = "aeroopt.solver-input.v1"
RESULT_SCHEMA = "aeroopt.solver-result.v1"


def _plain(value: object) -> object:
    if value is None or type(value) in (str, bool, int):
        return value
    if type(value) is float:
        if not math.isfinite(value):
            raise ValueError("JSON numbers must be finite")
        return value
    if isinstance(value, Mapping):
        if any(type(key) is not str for key in value):
            raise TypeError("JSON object keys must be strings")
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    raise TypeError(f"Unsupported JSON value: {type(value).__name__}")


def canonical_json(value: object) -> str:
    """Sorted UTF-8 JSON; float values are finite CPython round-trip decimals."""
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)


def scientific_hash(value: object) -> str:
    """Hash a tagged JSON tree; float nodes have exact hexadecimal encoding.

    Containers are tagged too, so a user dict cannot impersonate a float tag.
    Neither array order nor Unicode text is normalized. Dict order is ignored.
    """
    def tagged(item: object) -> object:
        if isinstance(item, dict):
            return ["object", [[key, tagged(val)] for key, val in sorted(item.items())]]
        if isinstance(item, list):
            return ["array", [tagged(val) for val in item]]
        if type(item) is float:
            return ["float64", (0.0 if item == 0 else item).hex()]
        if item is None:
            return ["null"]
        return [type(item).__name__, item]

    payload = canonical_json(tagged(_plain(value))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, init=False)
class FrozenJSON:
    """Deep immutable object snapshot; to_dict always returns a detached copy."""

    _text: str = field(repr=False)

    def __init__(self, value: Mapping | None = None) -> None:
        if value is not None and not isinstance(value, Mapping):
            raise TypeError("FrozenJSON requires a mapping")
        object.__setattr__(self, "_text", canonical_json({} if value is None else value))

    def to_dict(self) -> dict:
        return json.loads(self._text)


def read_json(path: str | Path) -> dict:
    """Reject nonstandard NaN/Infinity, duplicate keys and overflowed decimals."""
    def constant(value: str) -> None:
        raise ValueError(f"Non-finite JSON constant: {value}")

    def pairs(items: list[tuple]) -> dict:
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    value = json.loads(native_path(path).read_text(encoding="utf-8"), parse_constant=constant,
                       object_pairs_hook=pairs)
    _plain(value)
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object")
    return value


def write_json(path: str | Path, value: object) -> None:
    """Publish complete JSON exclusively via a same-directory hard link.

    Requires a hard-link-capable local filesystem (NTFS on the target Windows).
    On I/O failure the existing target is untouched; a pending file may remain.
    """
    path = native_path(path)
    payload = (canonical_json(value) + "\n").encode("utf-8")
    pending = path.with_name(f".{path.name}.{uuid4().hex}.pending")
    with pending.open("xb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        os.link(pending, path)
    finally:
        pending.unlink()
