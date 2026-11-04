"""Offline reproduction and byte-preserving opt-in download safeguards."""

import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

from examples import _reference_data

ROOT = Path(__file__).resolve().parents[1]


def _catalog(tmp_path, monkeypatch, raw):
    data = [{"file": "source.dat", "sha256": hashlib.sha256(raw).hexdigest(),
             "url": "https://example.invalid/source.dat", "name": "synthetic test"}]
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.setattr(_reference_data, "CATALOG_PATH", path)


def test_download_preserves_bytes_and_refuses_overwrite(tmp_path, monkeypatch):
    raw = b"Synthetic\r\n1 0\r\n.5 .1\r\n0 0\r\n.5 -.1\r\n1 0\r\n"
    _catalog(tmp_path, monkeypatch, raw)
    monkeypatch.setattr(_reference_data, "urlopen", lambda *args, **kwargs: io.BytesIO(raw))
    directory = tmp_path / "sources"
    records = _reference_data.published_sources(directory, fetch=True)
    assert (directory / "source.dat").read_bytes() == raw
    assert records[0]["obtained_on"]
    assert records[0]["redistribution"].startswith("local use only")
    with pytest.raises(FileExistsError):
        _reference_data.published_sources(directory, fetch=True)
    assert (directory / "source.dat").read_bytes() == raw


def test_changed_upstream_hash_requires_explicit_review(tmp_path, monkeypatch):
    _catalog(tmp_path, monkeypatch, b"original")
    monkeypatch.setattr(_reference_data, "urlopen", lambda *args, **kwargs: io.BytesIO(b"changed"))
    directory = tmp_path / "sources"
    with pytest.raises(ValueError, match="Source hash changed"):
        _reference_data.published_sources(directory, fetch=True)
    assert (directory / "source.dat").read_bytes() == b"changed"


def test_local_changed_bytes_fail_hash_check(tmp_path, monkeypatch):
    _catalog(tmp_path, monkeypatch, b"original")
    (tmp_path / "source.dat").write_bytes(b"changed")
    with pytest.raises(ValueError, match="Source hash changed"):
        _reference_data.published_sources(tmp_path)


def test_offline_diagnostics_cli_and_no_silent_overwrite(tmp_path):
    output = tmp_path / "diagnostics"
    command = [sys.executable, str(ROOT / "examples/validate_real_airfoils.py"), "--output", str(output)]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    assert summary["published_profile_count"] == 0
    assert summary["sensitivity_rows"] == 32
    assert summary["configurations"] == 10
    assert summary["roundtrips"] == 5
    assert summary["max_relative_area_error"] < 2e-12
    assert (output / "normalized_profiles.png").is_file()
    assert (output / "discretization_errors.svg").is_file()
    before = (output / "summary.json").read_bytes()
    repeated = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    assert repeated.returncode != 0
    assert "Refusing to overwrite" in repeated.stderr
    assert (output / "summary.json").read_bytes() == before
