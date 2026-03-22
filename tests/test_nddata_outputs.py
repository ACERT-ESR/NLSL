import numpy as np
import pytest

import nlsl
from nlsl.data import process_spectrum
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    NSPLINE_POINTS,
    SAMPL4_DATA_PATH,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_INITIAL_PARAMETERS,
    SAMPL4_POINT_COUNT,
    SAMPL4_SPECTRAL_DATA,
)

pytest.importorskip(
    "pyspecdata", exc_type=ImportError, reason="pyspecdata is required"
)
from pyspecdata.core import nddata


def _load_processed_sampl4(model):
    processed = process_spectrum(
        SAMPL4_DATA_PATH,
        NSPLINE_POINTS,
        BASELINE_EDGE_POINTS,
        derivative_mode=model.derivative_mode,
        normalize=False,
    )
    idx = model.generate_coordinates(
        processed.start,
        processed.stop,
        processed.y.size,
    )
    model.data = processed.y
    model.name(str(SAMPL4_DATA_PATH.with_suffix("")), spectrum=idx)
    model.noise(processed.noise, spectrum=idx)


def _sampl4_dataset(offset=0.0, label="field", units="mT"):
    fields = (
        SAMPL4_FIELD_START
        + SAMPL4_FIELD_STEP * np.arange(SAMPL4_POINT_COUNT)
        + float(offset)
    )
    dataset = nddata(
        SAMPL4_SPECTRAL_DATA.copy(),
        [SAMPL4_POINT_COUNT],
        [label],
    )
    dataset.name("sampl4.dat")
    dataset.setaxis(label, fields)
    dataset.set_units(label, units)
    return dataset, fields


def test_load_raw_datafile_keeps_array_outputs():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    _load_processed_sampl4(model)

    experimental = model.data
    site_spectra = model.current_spectrum
    weights = model.weights

    assert model.return_nddata is False
    assert isinstance(experimental, np.ndarray)
    assert isinstance(site_spectra, np.ndarray)
    assert isinstance(weights, np.ndarray)


def test_multi_spectrum_raw_data_returns_list():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    model.procline("series psi 0 90")
    _load_processed_sampl4(model)
    _load_processed_sampl4(model)

    experimental = model.data

    assert isinstance(experimental, list)
    assert len(experimental) == 2
    assert all(isinstance(item, np.ndarray) for item in experimental)


def test_nddata_assignment_wraps_data_current_spectrum_and_weights():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    dataset, fields = _sampl4_dataset()

    model.data = dataset

    experimental = model.data
    site_spectra = model.current_spectrum
    weights = model.weights

    assert model.return_nddata is True
    assert isinstance(experimental, nddata)
    assert experimental.dimlabels == ["field"]
    assert np.allclose(experimental["field"], fields)
    assert experimental.get_units("field") == "mT"
    assert experimental.name() == "sampl4.dat"
    assert isinstance(site_spectra, nddata)
    assert site_spectra.dimlabels == ["sites", "field"]
    assert np.array_equal(site_spectra["sites"], np.arange(2))
    assert np.allclose(site_spectra["field"], fields)
    assert site_spectra.get_units("field") == "mT"

    assert isinstance(weights, nddata)
    assert weights.dimlabels == ["sites"]
    assert np.array_equal(weights["sites"], np.arange(2))

    total = weights @ site_spectra

    assert isinstance(total, nddata)
    assert total.dimlabels == ["field"]
    assert np.allclose(total["field"], fields)
    assert total.get_units("field") == "mT"
    manual_total = site_spectra["sites", 0] * weights["sites", 0].item()
    for site_index in range(1, len(weights.getaxis("sites"))):
        manual_total += (
            site_spectra["sites", site_index]
            * weights["sites", site_index].item()
        )
    manual_total_difference = total - manual_total
    manual_total_difference.run(np.abs)
    assert manual_total_difference.max() < 1.0e-12


def test_fit_returns_nddata_after_nddata_assignment():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    model.fit_params["maxitr"] = 0
    model.fit_params["maxfun"] = 1
    dataset, fields = _sampl4_dataset(units="G")

    model.data = dataset

    fit_site_spectra = model.fit()
    current_site_spectra = model.current_spectrum

    assert isinstance(fit_site_spectra, nddata)
    assert fit_site_spectra.dimlabels == ["sites", "field"]
    assert np.allclose(fit_site_spectra["field"], fields)
    assert fit_site_spectra.get_units("field") == "G"
    fit_difference = fit_site_spectra - current_site_spectra
    fit_difference.run(np.abs)
    assert fit_difference.max() < 1.0e-12


def test_multi_spectrum_nddata_data_returns_list():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    model.procline("series psi 0 90")
    dataset_a, fields_a = _sampl4_dataset(offset=0.0, units="mT")
    dataset_b, fields_b = _sampl4_dataset(
        offset=0.25 * SAMPL4_FIELD_STEP,
        units="G",
    )

    model.data = dataset_a
    model.data = dataset_b

    experimental = model.data

    assert isinstance(experimental, list)
    assert len(experimental) == 2
    assert all(isinstance(item, nddata) for item in experimental)
    assert np.allclose(experimental[0]["field"], fields_a)
    assert np.allclose(experimental[1]["field"], fields_b)
    assert experimental[0].get_units("field") == "mT"
    assert experimental[1].get_units("field") == "G"


def test_current_spectrum_requires_shared_multi_spectrum_field_axis():
    model = nlsl.nlsl()
    model.update(SAMPL4_INITIAL_PARAMETERS)
    model.shift = True
    model.procline("series psi 0 90")
    dataset_a, _ = _sampl4_dataset(offset=0.0)
    dataset_b, _ = _sampl4_dataset(offset=0.25 * SAMPL4_FIELD_STEP)

    model.data = dataset_a
    model.data = dataset_b

    with pytest.raises(ValueError, match="common field axis"):
        model.current_spectrum
