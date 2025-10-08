import numpy as np
import pytest

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_FINAL_COMMANDS,
    SAMPL4_FINAL_WEIGHTS,
    SAMPL4_INTENSITIES,
    SAMPL4_POINT_COUNT,
)


def test_sampl4_best_parameters_match_data_without_fit():
    """Copy the converged parameters into a fresh model and verify the residual."""

    model = nlsl.nlsl()
    model['nsite'] = 2

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

    nlsl.fortrancore.expdat.data[data_slice] = SAMPL4_INTENSITIES[: data_slice.stop - data_slice.start]

    for command in SAMPL4_FINAL_COMMANDS:
        model.procline(command)

    nlsl.fortrancore.mspctr.sfac[:, index] = 0.0
    nlsl.fortrancore.mspctr.sfac[: SAMPL4_FINAL_WEIGHTS.size, index] = SAMPL4_FINAL_WEIGHTS

    site_spectra, weights = model.current_spectrum

    simulated_total = np.dot(weights, site_spectra[:, data_slice])
    experimental = SAMPL4_INTENSITIES[: data_slice.stop - data_slice.start]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
