"""Analytic fixtures; caches and headless plotting stay inside the workspace."""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pytest

from aeroopt.geometry import Airfoil, FlapParameters, MultiElementAirfoil


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

