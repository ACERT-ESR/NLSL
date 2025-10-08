import numpy as np
import pytest

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_FINAL_FPARM,
    SAMPL4_FINAL_IPARM,
    SAMPL4_FINAL_PARAMETERS,
    SAMPL4_FINAL_SB0,
    SAMPL4_FINAL_SHFT,
    SAMPL4_FINAL_SRNG,
    SAMPL4_FINAL_ISHFT,
    SAMPL4_FINAL_NRMLZ,
    SAMPL4_FINAL_WEIGHTS,
    SAMPL4_INTENSITIES,
    SAMPL4_POINT_COUNT,
)


def test_sampl4_best_parameters_match_data_without_fit():
    """Copy the converged parameters into a fresh model and verify the residual."""

    model = nlsl.nlsl()
    model['nsite'] = SAMPL4_FINAL_PARAMETERS['nsite']

    index, data_slice = model.generate_coordinates(
        SAMPL4_POINT_COUNT,
        start=SAMPL4_FIELD_START,
        step=SAMPL4_FIELD_STEP,
        derivative_mode=DERIVATIVE_MODE,
        baseline_points=BASELINE_EDGE_POINTS,
        normalize=False,
        nspline=NSPLINE_POINTS,
        shift=True,
        label='sampl4-known-parameters',
        reset=True,
    )

    count = data_slice.stop - data_slice.start
    model.set_data(data_slice, SAMPL4_INTENSITIES[:count])

    model.apply_parameter_state(SAMPL4_FINAL_FPARM, SAMPL4_FINAL_IPARM)
    model.set_spectral_state(
        sb0=SAMPL4_FINAL_SB0,
        srng=SAMPL4_FINAL_SRNG,
        ishift=SAMPL4_FINAL_ISHFT,
        shift=SAMPL4_FINAL_SHFT,
        normalize_flags=SAMPL4_FINAL_NRMLZ,
    )

    model.set_site_weights(index, SAMPL4_FINAL_WEIGHTS)

    site_spectra, weights = model.current_spectrum

    simulated_total = np.dot(weights, site_spectra[:, data_slice])
    experimental = SAMPL4_INTENSITIES[:count]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
