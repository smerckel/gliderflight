"""End-to-end tests for the steady-state glider flight model.

These exercise the angle-of-attack solver (``SteadyStateGliderModel``), the code
path that regressed under NumPy 2 when a size-1 array is assigned to a scalar
slot -- so they double as a regression guard for that fix.
"""
import numpy as np
import pytest

from gliderflight import Modelresult, SteadyStateCalibrate, SteadyStateGliderModel


def _solve(data, **coeffs):
    gm = SteadyStateGliderModel(rho0=1004)
    gm.define(Vg=70 / 1004, ah=3.5, mg=70.100, Cd0=0.16)
    gm.define(**coeffs)
    gm.solve(dict(time=data["time"], pressure=data["pressure"], pitch=data["pitch"],
                  buoyancy_change=data["buoyancy_change"], density=data["density"]))
    return gm


def test_solve_returns_modelresult(glider_data):
    gm = _solve(glider_data)
    assert isinstance(gm.modelresult, Modelresult)
    assert gm.modelresult._fields == ("t", "u", "w", "U", "alpha", "pitch", "ww", "depth")


def test_solve_outputs_finite_and_correct_length(glider_data):
    gm = _solve(glider_data)
    n = len(glider_data["time"])
    for name in ("t", "U", "alpha", "wg", "w"):
        arr = getattr(gm, name)
        assert len(arr) == n, name
        assert np.isfinite(arr).all(), name


def test_incident_velocity_is_physically_sane(glider_data):
    gm = _solve(glider_data)
    # Speed through water is non-negative and stays well below 1 m/s for a glider.
    assert (gm.U >= 0).all()
    assert gm.U.max() < 1.0
    # The bulk of the dive has a small angle of attack (a few degrees).
    assert np.median(np.abs(gm.alpha)) < np.radians(10)


def test_calibrate_recovers_expected_parameters(glider_data):
    gm = SteadyStateCalibrate(rho0=1004)
    gm.define(ah=3.8, Cd1=10.5)
    gm.set_input_data(time=glider_data["time"], pressure=glider_data["pressure"],
                      pitch=glider_data["pitch"], buoyancy_change=glider_data["buoyancy_change"],
                      density=glider_data["density"])
    # Exclude the near-surface and deep transition points from the cost function.
    gm.OR(glider_data["pressure"] * 10 < 10)
    gm.OR(glider_data["pressure"] * 10 > 35)
    gm.define(Cd0=0.15, mg=70, Vg=70 / 1004)

    result = gm.calibrate("Cd0", "mg")

    assert set(result) == {"Cd0", "mg"}
    # Physically plausible parasite drag and a mass close to the ~70 kg input.
    assert 0.05 < result["Cd0"] < 0.5
    assert 60 < result["mg"] < 80
    # The model attributes are left at (approximately) the optimum in-place.
    assert gm.Cd0 == pytest.approx(result["Cd0"], rel=1e-2)
