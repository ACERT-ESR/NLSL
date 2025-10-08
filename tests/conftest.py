import numpy as np
import pytest
from pathlib import Path

import nlsl


NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

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

FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
    "ftol": 1.0e-2,
    "gtol": 1.0e-6,
    "xtol": 1.0e-4,
}

PARAMETERS_TO_VARY = ["gib0", "rbar(1)", "rbar(2)"]

DATA_PATH = Path(__file__).resolve().parent.parent / "examples" / "sampl4.dat"


@pytest.fixture(scope="module")
def sampl4_fit_result():
    """Run the Python reproduction of ``sampl4.run`` once for sharing."""

    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    model.load_data(
        DATA_PATH,
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=False,
        derivative_mode=DERIVATIVE_MODE,
    )

    for token in PARAMETERS_TO_VARY:
        model.procline(f"vary {token}")

    for key, value in FIT_CONTROLS.items():
        model.fit_params[key] = value

    model.fit()
    site_spectra, weights = model.fit()
    site_spectra = site_spectra.copy()
    weights = weights.copy()

    start_index = int(model.layout["ixsp"][0])
    point_count = int(model.layout["npts"][0])
    data_slice = slice(start_index, start_index + point_count)

    experimental = nlsl.fortrancore.expdat.data[data_slice].copy()
    components = site_spectra[:, data_slice]
    weighted_components = weights[0][:, np.newaxis] * components
    simulated_total = np.sum(weighted_components, axis=0)

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
