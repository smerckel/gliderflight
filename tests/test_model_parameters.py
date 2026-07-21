"""Unit tests for the parameter/masking machinery that needs no glider data."""
import numpy as np
import pytest

import gliderflight
from gliderflight import SteadyStateCalibrate, SteadyStateGliderModel


def test_define_and_get_settings():
    gm = SteadyStateGliderModel(rho0=1004)
    gm.define(mg=70.0, Vg=68e-3, Cd0=0.16)
    settings = gm.get_settings()
    assert settings["mg"] == 70.0
    assert settings["Vg"] == 68e-3
    assert settings["Cd0"] == 0.16


def test_undefined_parameters_tracks_required_coefficients():
    gm = SteadyStateGliderModel(rho0=1004)
    # A fresh model needs mass, volume and parasite drag to be supplied.
    assert set(gm.undefined_parameters()) == {"Cd0", "mg", "Vg"}
    gm.define(mg=70.0, Vg=68e-3)
    assert list(gm.undefined_parameters()) == ["Cd0"]
    # Once every required coefficient is set, nothing is reported as undefined.
    gm.define(Cd0=0.16)
    assert list(gm.undefined_parameters()) == []


def test_copy_settings_transfers_coefficients():
    src = SteadyStateGliderModel(rho0=1004)
    src.define(mg=71.0, Vg=69e-3, Cd0=0.17, ah=3.6, Cd1=9.5)
    dst = gliderflight.DynamicGliderModel(rho0=1004, k1=0.2, k2=0.98, dt=1.0)
    dst.copy_settings(src)
    for key, value in src.get_settings().items():
        assert dst.get_settings()[key] == value


def test_has_aoa_parameter_changed_tracks_lift_coefficients():
    gm = SteadyStateGliderModel(rho0=1004)
    # ah is an angle-of-attack parameter; changing it flips the flag.
    gm.define(ah=3.8)
    assert gm.has_aoa_parameter_changed() is True


def _make_calibrator(n=20):
    gm = SteadyStateCalibrate(rho0=1004)
    time = np.arange(n, dtype=float)
    pressure = np.linspace(0.0, 6.0, n)  # bar
    pitch = np.full(n, -0.4)
    buoyancy_change = np.zeros(n)
    density = np.full(n, 1004.0)
    gm.define(mg=70.0, Vg=68e-3, Cd0=0.16)
    gm.set_input_data(time=time, pressure=pressure, pitch=pitch,
                      buoyancy_change=buoyancy_change, density=density)
    return gm, pressure


def test_mask_defaults_to_all_false():
    gm, _ = _make_calibrator()
    # set_input_data initialises the mask to all-False (i.e. keep every point).
    assert gm.mask.dtype == bool
    assert not gm.mask.any()


def test_or_masks_matching_points():
    gm, pressure = _make_calibrator()
    condition = pressure * 10 < 10  # mask the shallow points
    gm.OR(condition)
    np.testing.assert_array_equal(gm.mask, condition)


def test_and_intersects_masks():
    gm, pressure = _make_calibrator()
    gm.OR(pressure * 10 < 40)          # start with a broad mask
    gm.AND(pressure * 10 > 10)         # intersect
    expected = (pressure * 10 < 40) & (pressure * 10 > 10)
    np.testing.assert_array_equal(gm.mask, expected)


def test_nand_clears_intersection():
    gm, pressure = _make_calibrator()
    gm.OR(pressure * 10 < 40)
    gm.NAND(pressure * 10 < 20)
    expected = (pressure * 10 < 40) & ~(pressure * 10 < 20)
    np.testing.assert_array_equal(gm.mask, expected)


def test_set_mask_bool_toggles_all():
    gm, _ = _make_calibrator()
    gm.set_mask(True)
    assert gm.mask.all()
    gm.set_mask(False)
    assert not gm.mask.any()
