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


def test_parameter_bounds_and_steps_mirror_parcom(clean_model):
    """Assigning lmfit bounds should mirror into the Fortran parcom slots."""

    model = clean_model
    model["gxx"] = 2.005

    param = model.parameters["gxx_0"]
    param.min = 1.9
    param.max = 2.1
    param.user_data["scale"] = 0.5
    param.user_data["fdstep"] = 1.0e-4
    param.vary = True

    slot = param._position()
    parcom = nlsl.fortrancore.parcom

    assert slot is not None
    assert int(parcom.nprm) == 1
    assert np.isclose(parcom.prmin[slot], 1.9)
    assert np.isclose(parcom.prmax[slot], 2.1)
    assert np.isclose(parcom.prscl[slot], 0.5)
    assert np.isclose(parcom.xfdstp[slot], 1.0e-4)
    assert int(parcom.ixpr[slot]) == model.parameter_index("gxx")
    assert int(parcom.ixst[slot]) == 0


def test_multi_site_entries_follow_parameter_slots(clean_model):
    """Each site-specific parameter should allocate one parcom slot."""

    model = clean_model
    model["nsite"] = 2
    model["gxx"] = (2.01, 2.02)
    model["gzz"] = (2.15, 2.25)

    first = model.parameters["gxx_0"]
    second = model.parameters["gxx_1"]
    only_gzz = model.parameters["gzz_1"]

    first.min = 1.9
    first.max = 2.2
    first.user_data["scale"] = 0.6
    first.user_data["fdstep"] = 1.0e-4
    second.min = 2.05
    second.max = 2.35
    second.user_data["scale"] = 0.85
    second.user_data["fdstep"] = 2.5e-4
    only_gzz.min = 2.1
    only_gzz.max = 2.6
    only_gzz.user_data["scale"] = 1.1
    only_gzz.user_data["fdstep"] = 3.0e-4

    first.vary = True
    second.vary = True
    only_gzz.vary = True

    parcom = nlsl.fortrancore.parcom

    assert int(parcom.nprm) == 3

    first_slot = first._position()
    second_slot = second._position()
    gzz_slot = only_gzz._position()

    assert first_slot is not None and second_slot is not None and gzz_slot is not None

    assert np.isclose(parcom.prmin[first_slot], 1.9)
    assert np.isclose(parcom.prmax[first_slot], 2.2)
    assert np.isclose(parcom.prscl[first_slot], 0.6)
    assert np.isclose(parcom.xfdstp[first_slot], 1.0e-4)
    assert int(parcom.ixst[first_slot]) == 1

    assert np.isclose(parcom.prmin[second_slot], 2.05)
    assert np.isclose(parcom.prmax[second_slot], 2.35)
    assert np.isclose(parcom.prscl[second_slot], 0.85)
    assert np.isclose(parcom.xfdstp[second_slot], 2.5e-4)
    assert int(parcom.ixst[second_slot]) == 2

    assert np.isclose(parcom.prmin[gzz_slot], 2.1)
    assert np.isclose(parcom.prmax[gzz_slot], 2.6)
    assert np.isclose(parcom.prscl[gzz_slot], 1.1)
    assert np.isclose(parcom.xfdstp[gzz_slot], 3.0e-4)
    assert int(parcom.ixst[gzz_slot]) == 2


def test_vary_toggle_removes_parameter(clean_model):
    """Clearing ``vary`` should drop the slot from ``parcom``."""

    model = clean_model
    model["gxx"] = 2.0
    param = model.parameters["gxx_0"]

    param.vary = True

    assert int(nlsl.fortrancore.parcom.nprm) == 1

    param.vary = False

    assert int(nlsl.fortrancore.parcom.nprm) == 0


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

    model.fortran_lm_engine.update(SAMPL4_FIT_CONTROLS)

    model.parameters["gib0_0"].vary = True

    starting_rx = np.array(model["rx"], dtype=float)
    starting_gib0 = np.array(model["gib0"], dtype=float)

    model.fit()

    final_rx = np.array(model["rx"], dtype=float)
    final_gib0 = np.array(model["gib0"], dtype=float)

    assert np.isclose(final_rx[1], starting_rx[1])
    assert not np.allclose(final_gib0, starting_gib0)
    assert np.isclose(final_rx[0], starting_rx[0])
