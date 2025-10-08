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
    SAMPL4_INTENSITIES,
    SAMPL4_POINT_COUNT,
    apply_sampl4_final_state,
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

    apply_sampl4_final_state(model)

    site_spectra = model.current_spectrum
    weights_matrix = np.array(model['weights'], copy=True)

    count = data_slice.stop - data_slice.start
    simulated_total = np.dot(np.atleast_2d(weights_matrix), site_spectra[:, data_slice])
    simulated_total = np.squeeze(simulated_total)
    experimental = SAMPL4_INTENSITIES[:count]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
