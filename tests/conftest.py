"""Analytic fixtures; caches and headless plotting stay inside the workspace."""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil


def pytest_addoption(parser):
    parser.addoption("--reference-data", type=Path, default=None,
                     help="Directory with SHA-256-pinned UIUC originals (not redistributed)")
    parser.addoption("--williams-data", type=Path, default=None,
                     help="Local pinned Williams PDF and checked transcription (not redistributed)")


def pytest_configure(config):
    if config.option.basetemp is None:
        from uuid import uuid4

        config.option.basetemp = str(ROOT / ".tmp" / f"pytest_{uuid4().hex}")


@pytest.fixture(params=["n0012.dat", "naca4412.dat", "s1223.dat"])
def published_reference(request):
    import hashlib
    import json

    directory = request.config.getoption("--reference-data")
    if directory is None:
        pytest.skip("UIUC originals not redistributed; provide --reference-data explicitly")
    catalog = json.loads((ROOT / "examples/reference_catalog.json").read_text(encoding="utf-8"))
    entry = next(row for row in catalog if row["file"] == request.param)
    path = directory / entry["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    return path, entry


@pytest.fixture
def diamond_points():
    # Chord=1, perpendicular diagonals=1 and .2, area=.1.
    return np.array([[1., 0.], [.5, .1], [0., 0.], [.5, -.1], [1., 0.]])


@pytest.fixture
def diamond(diamond_points):
    return Airfoil.from_points(diamond_points, name="analytic diamond")


@pytest.fixture
def separated(diamond):
    return MultiElementAirfoil(diamond, diamond, FlapParameters(.3, (1.2, 0.)))

