import numpy as np
import pytest

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_DATA_PATH,
    SAMPL4_FIT_CONTROLS,
    SAMPL4_INITIAL_PARAMETERS,
)


@pytest.fixture
def clean_model():
    """Return a freshly initialised model so Fortran state is isolated."""

    return nlsl.nlsl()


def test_vary_stores_bounds_and_steps(clean_model):
    """Varying a scalar parameter should mirror limits and steps to parcom."""

    model = clean_model
    model["gxx"] = 2.005

    vary_mapping = model.fit_params.vary
    vary_mapping["gxx"] = {
        "minimum": 1.9,
        "maximum": 2.1,
        "scale": 0.5,
        "fdstep": 1.0e-4,
    }

    entry = vary_mapping["gxx"]

    assert np.isclose(entry["minimum"], 1.9)
    assert np.isclose(entry["maximum"], 2.1)
    assert np.isclose(entry["scale"], 0.5)

    # Verify the stored step mirrors the provided finite-difference value.
    expected_step_value = 1.0e-4
    assert np.isclose(entry["fdstep"], expected_step_value)
    assert entry["index"] == 0

    parameter_code = model.parameter_index("gxx")
    position = vary_mapping._entries(parameter_code)[0][1]
    parcom = model.fit_params._core.parcom

    assert int(parcom.nprm) == 1
    assert np.isclose(parcom.prmin[position], 1.9)
    assert np.isclose(parcom.prmax[position], 2.1)
    assert np.isclose(parcom.prscl[position], 0.5)
    assert np.isclose(parcom.xfdstp[position], expected_step_value)


def test_vary_records_fortran_slots(clean_model):
    """The Fortran parcom arrays should store limits and steps at each slot."""

    model = clean_model
    model["nsite"] = 2
    model["gxx"] = (2.01, 2.02)
    model["gzz"] = (2.15, 2.25)

    vary_mapping = model.fit_params.vary
    vary_mapping["gxx"] = {
        "minimum": [1.9, 2.05],
        "maximum": [2.2, 2.35],
        "scale": [0.6, 0.85],
        "fdstep": [1.0e-4, 2.5e-4],
    }
    vary_mapping["gzz"] = {
        "index": [2],
        "minimum": [2.1],
        "maximum": [2.6],
        "scale": [1.1],
        "fdstep": [3.0e-4],
    }

    parcom = model.fit_params._core.parcom
    gxx_code = model.parameter_index("gxx")
    gzz_code = model.parameter_index("gzz")

    # The bookkeeping structure assigns one slot per varied site.  Each slot
    # records the parameter code, site index, and bound/step arrays.
    assert int(parcom.nprm) == 3

    gxx_entries = vary_mapping._entries(gxx_code)
    ordered_gxx = sorted(gxx_entries, key=lambda item: item[0])
    assert len(ordered_gxx) == 2

    first_index, first_position = ordered_gxx[0]
    assert first_index == 1
    assert np.isclose(parcom.prmin[first_position], 1.9)
    assert np.isclose(parcom.prmax[first_position], 2.2)
    assert np.isclose(parcom.prscl[first_position], 0.6)
    assert np.isclose(parcom.xfdstp[first_position], 1.0e-4)
    assert int(parcom.ixpr[first_position]) == gxx_code
    assert int(parcom.ixst[first_position]) == first_index

    second_index, second_position = ordered_gxx[1]
    assert second_index == 2
    assert np.isclose(parcom.prmin[second_position], 2.05)
    assert np.isclose(parcom.prmax[second_position], 2.35)
    assert np.isclose(parcom.prscl[second_position], 0.85)
    assert np.isclose(parcom.xfdstp[second_position], 2.5e-4)
    assert int(parcom.ixpr[second_position]) == gxx_code
    assert int(parcom.ixst[second_position]) == second_index

    gzz_entries = vary_mapping._entries(gzz_code)
    assert len(gzz_entries) == 1
    gzz_index, gzz_position = gzz_entries[0]
    assert gzz_index == 2
    assert np.isclose(parcom.prmin[gzz_position], 2.1)
    assert np.isclose(parcom.prmax[gzz_position], 2.6)
    assert np.isclose(parcom.prscl[gzz_position], 1.1)
    assert np.isclose(parcom.xfdstp[gzz_position], 3.0e-4)
    assert int(parcom.ixpr[gzz_position]) == gzz_code
    assert int(parcom.ixst[gzz_position]) == gzz_index


def test_vary_toggle_removes_parameter(clean_model):
    """Assigning ``False`` should drop an existing vary entry."""

    model = clean_model
    model["gxx"] = 2.0

    vary_mapping = model.fit_params.vary
    vary_mapping["gxx"] = True

    assert "gxx" in vary_mapping
    assert int(model.fit_params._core.parcom.nprm) == 1

    vary_mapping["gxx"] = False

    assert "gxx" not in vary_mapping
    assert int(model.fit_params._core.parcom.nprm) == 0


def test_vary_rejects_out_of_range_indices(clean_model):
    """Index lists larger than the site count should raise immediately."""

    model = clean_model
    model["nsite"] = 2

    vary_mapping = model.fit_params.vary

    with pytest.raises(ValueError):
        vary_mapping["gxx"] = {"index": [1, 3]}

    assert int(model.fit_params._core.parcom.nprm) == 0


def test_vary_tracks_multiple_indices(clean_model):
    """Explicit index arrays should be preserved and ordered when queried."""

    model = clean_model
    model["nsite"] = 2
    model["gxx"] = (2.0, 2.2)

    vary_mapping = model.fit_params.vary
    vary_mapping["gxx"] = {
        "index": [2, 1],
        "minimum": [1.8, 1.9],
        "maximum": [2.4, 2.5],
        "scale": [0.5, 0.75],
        "fdstep": [2.0e-4, 1.0e-4],
    }

    entry = vary_mapping["gxx"]

    assert np.array_equal(entry["index"], np.array([1, 2], dtype=int))
    assert np.allclose(entry["minimum"], np.array([1.9, 1.8]))
    assert np.allclose(entry["maximum"], np.array([2.5, 2.4]))
    assert np.allclose(entry["scale"], np.array([0.75, 0.5]))

    expected_steps = np.array([1.0e-4, 2.0e-4])
    assert np.allclose(entry["fdstep"], expected_steps)

    parcom = model.fit_params._core.parcom
    assert int(parcom.nprm) == 2


def test_fixed_parameter_remains_constant_during_fit():
    """Fits should leave parameters unchanged when they are marked fixed."""

    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)

    model.load_data(
        SAMPL4_DATA_PATH,
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=False,
        derivative_mode=DERIVATIVE_MODE,
    )

    for key in SAMPL4_FIT_CONTROLS:
        model.fit_params[key] = SAMPL4_FIT_CONTROLS[key]

    vary_mapping = model.fit_params.vary
    vary_mapping["gib0"] = True
    vary_mapping["rx"] = {"index": [1]}

    starting_rx = np.array(model["rx"], dtype=float)
    starting_gib0 = np.array(model["gib0"], dtype=float)

    model.fit()

    final_rx = np.array(model["rx"], dtype=float)
    final_gib0 = np.array(model["gib0"], dtype=float)

    assert np.isclose(final_rx[1], starting_rx[1])
    assert not np.allclose(final_gib0, starting_gib0)
    assert not np.isclose(final_rx[0], starting_rx[0], rtol=1.0e-8, atol=1.0e-8)
