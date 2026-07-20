"""Shared fixtures for the gliderflight test suite.

The test data set ``data/gliderflight_data.npy`` is a subset of a real glider
mission (the same data used in the JTECH paper and the ``examples`` scripts).
It is a 6 x N array holding, in order: ctd time, pressure, pitch, buoyancy
change, in-situ density and ADCP-measured incident velocity.
"""
from pathlib import Path

import numpy as np
import pytest

# The data lives at the repository root, one level above this tests/ directory.
DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "gliderflight_data.npy"

# Background density used throughout the reference data set.
RHO0 = 1004


@pytest.fixture(scope="session")
def raw_data():
    """The full (6, N) array loaded from disk."""
    assert DATA_FILE.is_file(), f"missing test data file: {DATA_FILE}"
    return np.load(DATA_FILE)


def _as_dict(array):
    time, pressure, pitch, buoyancy_change, density, u_adcp = array
    return dict(
        time=time,
        pressure=pressure,
        pitch=pitch,
        buoyancy_change=buoyancy_change,
        density=density,
        u_adcp=u_adcp,
    )


@pytest.fixture(scope="session")
def glider_data(raw_data):
    """The full reference data set as a dict of named arrays."""
    return _as_dict(raw_data)


@pytest.fixture(scope="session")
def glider_data_subset(raw_data):
    """A short (first 1500 samples) slice, for the slower dynamic tests."""
    return _as_dict(raw_data[:, :1500])
