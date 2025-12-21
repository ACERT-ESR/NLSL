import numpy as np
import pytest

import nlsl


def test_parameters_share_storage_and_sync():
    model = nlsl.nlsl()
    model.nsites = 2
    model["gxx"] = [2.0, 3.0]

    row = model._fepr_names.index("gxx")
    assert np.shares_memory(model["gxx"], model._fparm[row, :2])

    model.parameters["gxx_0"].value = 4.5

    assert model._fparm[row, 0] == pytest.approx(4.5)
    assert model["gxx"][0] == pytest.approx(4.5)

    model["gxx"] = [5.0, 6.0]
    assert model.parameters["gxx_1"].value == pytest.approx(6.0)

    model.parameters["gxx_1"].value = 7.25
    assert model._fparm[row, 1] == pytest.approx(7.25)
    assert model["gxx"][1] == pytest.approx(7.25)


def test_parameters_vary_updates_fortran_list():
    model = nlsl.nlsl()
    model.nsites = 1

    index_code = model.parameter_index("gxx")
    model.parameters["gxx_0"].vary = True

    parcom = model._core.parcom
    entries = [int(parcom.ixpr[idx]) for idx in range(int(parcom.nprm))]
    assert index_code in entries

    model.parameters["gxx_0"].vary = False
    entries = [int(parcom.ixpr[idx]) for idx in range(int(parcom.nprm))]
    assert index_code not in entries


def test_spectral_parameter_updates_expdat():
    model = nlsl.nlsl()
    model.nspec = 1

    model.parameters["b0"].value = 3400.0

    assert model._core.expdat.sb0[0] == pytest.approx(3400.0)
    assert model["b0"] == pytest.approx(3400.0)


def test_legacy_vary_path_raises():
    model = nlsl.nlsl()

    with pytest.raises(RuntimeError):
        _ = model.fit_params.vary

    with pytest.raises(RuntimeError):
        model.fit_params.vary = {}
