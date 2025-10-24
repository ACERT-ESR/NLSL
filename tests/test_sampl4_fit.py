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

    components = site_spectra
    weights = np.array(model["weights"], copy=True)
    total_points = components.shape[1]
    weighted_components = weights[:, :, np.newaxis] * components[np.newaxis, :, :]
    simulated_total = np.squeeze(weights @ components)

    experimental_block = np.array(model.experimental_data, copy=True)
    experimental = np.squeeze(experimental_block)
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    return {
        "model": model,
        "site_spectra": site_spectra,
        "weights": weights,
        "components": components,
        "weighted_components": weighted_components.reshape(-1, total_points),
        "simulated_total": simulated_total,
        "experimental": experimental,
        "rel_rms": rel_rms,
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

    assert np.allclose(site_spectra_cs, sampl4_fit_result["site_spectra"], atol=5e-6)
    assert np.allclose(
        sampl4_fit_result["model"]["weights"],
        sampl4_fit_result["weights"],
    )

    recomputed = sampl4_fit_result["model"]["weights"] @ site_spectra_cs

    assert np.allclose(
        np.squeeze(recomputed),
        sampl4_fit_result["simulated_total"],
        atol=3e-6,
    )
