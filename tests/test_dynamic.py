"""End-to-end tests for the dynamic (Runge-Kutta integrated) glider model.

The dynamic ``solve`` exercises the parallel time-interval computation that
regressed under NumPy 2 (``ndarray.ptp`` -> ``np.ptp``), so these act as a
regression guard for that fix as well.
"""
import numpy as np
import pytest

from gliderflight import DynamicCalibrate, DynamicGliderModel, Modelresult


def _solve(data, dt=2):
    dm = DynamicGliderModel(rho0=1004, k1=0.2, k2=0.98, dt=dt)
    dm.define(Vg=70 / 1004, ah=3.5, mg=70.100, Cd0=0.16)
    dm.solve(dict(time=data["time"], pressure=data["pressure"], pitch=data["pitch"],
                  buoyancy_change=data["buoyancy_change"], density=data["density"]))
    return dm


def test_dynamic_solve_finite_and_sane(glider_data_subset):
    dm = _solve(glider_data_subset)
    n = len(glider_data_subset["time"])
    assert isinstance(dm.modelresult, Modelresult)
    assert len(dm.U) == n
    assert np.isfinite(dm.U).all()
    assert (dm.U >= 0).all()
    assert dm.U.max() < 1.0


def test_dynamic_and_steady_are_comparable(glider_data_subset):
    """The dynamic and steady models should broadly agree on mean speed."""
    from gliderflight import SteadyStateGliderModel

    dm = _solve(glider_data_subset)
    sm = SteadyStateGliderModel(rho0=1004)
    sm.define(Vg=70 / 1004, ah=3.5, mg=70.100, Cd0=0.16)
    sm.solve(dict(time=glider_data_subset["time"], pressure=glider_data_subset["pressure"],
                  pitch=glider_data_subset["pitch"],
                  buoyancy_change=glider_data_subset["buoyancy_change"],
                  density=glider_data_subset["density"]))
    # Mean incident velocities agree to within 5 cm/s.
    assert abs(np.nanmean(dm.U) - np.nanmean(sm.U)) < 0.05


@pytest.mark.slow
def test_dynamic_calibrate_smoke(glider_data_subset):
    dm = DynamicCalibrate(rho0=1004, k1=0.2, k2=0.98, dt=2)
    dm.define(ah=3.8, Cd1=10.5)
    dm.set_input_data(time=glider_data_subset["time"], pressure=glider_data_subset["pressure"],
                      pitch=glider_data_subset["pitch"],
                      buoyancy_change=glider_data_subset["buoyancy_change"],
                      density=glider_data_subset["density"])
    dm.OR(glider_data_subset["pressure"] * 10 < 10)
    dm.OR(glider_data_subset["pressure"] * 10 > 35)
    dm.define(Cd0=0.15, mg=70, Vg=70 / 1004)

    result = dm.calibrate("Cd0")

    assert set(result) == {"Cd0"}
    assert 0.05 < result["Cd0"] < 0.5
    assert dm.Cd0 == pytest.approx(result["Cd0"], rel=1e-2)
