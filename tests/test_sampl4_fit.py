import numpy as np
import pytest
from pathlib import Path

import nlsl

NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

# These match the ``let`` statements executed by ``sampl4.run`` and serve as the
# starting guesses for the optimisation.
INITIAL_PARAMETERS = {
    "nsite": 2,
    "in2": 2,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "lemx": 12,
    "lomx": 10,
    "kmx": 7,
    "mmx": 7,
    "ipnmx": 2,
    "gib0": 0.5,
    "rx": np.array([np.log10(3.0e8), np.log10(1.0e7)]),
}

# The tolerance map reproduces the ``fit`` controls from the runfile.
FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
    "ftol": 1.0e-2,
    "gtol": 1.0e-6,
    "xtol": 1.0e-4,
}

# The classic interface uses these short commands to release parameters.
PARAMETERS_TO_VARY = ["gib0", "rbar(1)", "rbar(2)"]


def run_pythonic_sampl4_fit():
    """Execute the two-stage ``sampl4`` fit and capture the resulting spectra."""

    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    examples_dir = Path(__file__).resolve().parent.parent / "examples"
    data_path = examples_dir / "sampl4.dat"

    model.load_data(
        data_path,
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=False,
        derivative_mode=DERIVATIVE_MODE,
    )

    for token in PARAMETERS_TO_VARY:
        model.procline(f"vary {token}")

    for key in FIT_CONTROLS:
        model.fit_params[key] = FIT_CONTROLS[key]

    # The historical run issues ``fit`` twice; repeating it here mirrors the
    # published optimisation cycle and ensures the spectra are stored.
    model.fit()
    site_spectra, weights = model.fit()
    site_spectra = site_spectra.copy()
    weights = weights.copy()

    start_index = int(model.layout["ixsp"][0])
    point_count = int(model.layout["npts"][0])
    data_slice = slice(start_index, start_index + point_count)

    components = site_spectra[:, data_slice]
    weighted_components = weights[0][:, np.newaxis] * components
    simulated_total = np.sum(weighted_components, axis=0)

    experimental = nlsl.fortrancore.expdat.data[data_slice].copy()
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    return {
        "model": model,
        "site_spectra": site_spectra,
        "weights": weights,
        "components": components,
        "weighted_components": weighted_components,
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

    site_spectra_cs, weights_cs = sampl4_fit_result["model"].current_spectrum

    assert np.allclose(site_spectra_cs, sampl4_fit_result["site_spectra"], atol=5e-6)
    assert np.allclose(weights_cs, sampl4_fit_result["weights"])

    recomputed = weights_cs[0] @ site_spectra_cs[:, sampl4_fit_result["data_slice"]]

    assert np.allclose(recomputed, sampl4_fit_result["simulated_total"], atol=3e-6)
