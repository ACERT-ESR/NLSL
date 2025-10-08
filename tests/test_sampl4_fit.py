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
    SAMPL4_INTENSITIES,
    SAMPL4_PARAMETERS_TO_VARY,
)


def run_pythonic_sampl4_fit():
    """Execute the two-stage ``sampl4`` fit and capture the resulting spectra."""

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

    for token in SAMPL4_PARAMETERS_TO_VARY:
        model.procline(f"vary {token}")

    for key in SAMPL4_FIT_CONTROLS:
        model.fit_params[key] = SAMPL4_FIT_CONTROLS[key]

    # The historical run issues ``fit`` twice; repeating it here mirrors the
    # published optimisation cycle and ensures the spectra are stored.
    model.fit()
    site_spectra = model.fit().copy()
    weights_value = np.array(model["weights"], copy=True)

    start_index = int(model.layout["ixsp"][0])
    point_count = int(model.layout["npts"][0])
    data_slice = slice(start_index, start_index + point_count)

    components = site_spectra[:, data_slice]
    weight_matrix = np.atleast_2d(weights_value)
    weighted_components = weight_matrix[:, :, np.newaxis] * components[np.newaxis, :, :]
    simulated_total = np.squeeze(weight_matrix.dot(components))

    experimental = SAMPL4_INTENSITIES[:point_count]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    return {
        "model": model,
        "site_spectra": site_spectra,
        "weights": weights_value,
        "components": components,
        "weighted_components": np.squeeze(weighted_components),
        "simulated_total": simulated_total,
        "experimental": experimental,
        "rel_rms": rel_rms,
        "data_slice": data_slice,
    }


@pytest.fixture(scope="module")
def sampl4_fit_result():
    """Share the cached ``sampl4`` fit across tests in this module."""

    return run_pythonic_sampl4_fit()


def test_pythonic_sampl4_fit_matches_data(sampl4_fit_result):
    """The refactored example should reach the published residual."""

    assert sampl4_fit_result["rel_rms"] < 0.0401


def test_current_spectrum_matches_fit_components(sampl4_fit_result):
    """A post-fit single-point evaluation must reproduce the fit spectra."""

    site_spectra_cs = sampl4_fit_result["model"].current_spectrum
    weights_matrix = np.array(sampl4_fit_result["model"]["weights"], copy=True)

    assert np.allclose(site_spectra_cs, sampl4_fit_result["site_spectra"], atol=5e-6)
    assert np.allclose(weights_matrix, sampl4_fit_result["weights"])

    recomputed = np.dot(
        np.atleast_2d(weights_matrix),
        site_spectra_cs[:, sampl4_fit_result["data_slice"]],
    )

    assert np.allclose(np.squeeze(recomputed), sampl4_fit_result["simulated_total"], atol=3e-6)
