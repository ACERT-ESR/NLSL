import numpy as np
import pytest

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_FINAL_PARAMETERS,
    SAMPL4_FINAL_SB0,
    SAMPL4_FINAL_SRNG,
    SAMPL4_FINAL_ISHFT,
    SAMPL4_FINAL_SHFT,
    SAMPL4_FINAL_NRMLZ,
    SAMPL4_FINAL_WEIGHTS,
    SAMPL4_SPECTRAL_DATA,
    SAMPL4_POINT_COUNT,
)


def test_sampl4_best_parameters_match_data_without_fit():
    """Copy the converged parameters into a fresh model and verify the
    residual."""

    model = nlsl.nlsl()
    model["nsite"] = SAMPL4_FINAL_PARAMETERS["nsite"]

    model.generate_coordinates(
        SAMPL4_POINT_COUNT,
        start=SAMPL4_FIELD_START,
        step=SAMPL4_FIELD_STEP,
        derivative_mode=DERIVATIVE_MODE,
        baseline_points=BASELINE_EDGE_POINTS,
        normalize=False,
        nspline=NSPLINE_POINTS,
        shift=True,
        label="sampl4-known-parameters",
        reset=True,
    )

    # Assign the converged runfile-4 state without invoking the optimiser.
    model.update(SAMPL4_FINAL_PARAMETERS)
    model["sb0"] = SAMPL4_FINAL_SB0
    model["srng"] = SAMPL4_FINAL_SRNG
    model["ishft"] = SAMPL4_FINAL_ISHFT
    model["shft"] = SAMPL4_FINAL_SHFT
    model["nrmlz"] = SAMPL4_FINAL_NRMLZ
    model.weights = SAMPL4_FINAL_WEIGHTS

    site_spectra = model.current_spectrum
    simulated_total = np.squeeze(model.weights @ site_spectra)
    experimental = SAMPL4_SPECTRAL_DATA[: site_spectra.shape[1]]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
