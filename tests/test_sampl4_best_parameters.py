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
    SAMPL4_FINAL_FPARM,
    SAMPL4_FINAL_IPARM,
    SAMPL4_FINAL_SB0,
    SAMPL4_FINAL_SRNG,
    SAMPL4_FINAL_ISHFT,
    SAMPL4_FINAL_SHFT,
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

    # Assign the converged runfile-4 state without invoking the optimiser.
    model.update(SAMPL4_FINAL_PARAMETERS)
    # Recreate the full set of fitted parameters captured during the runfile
    # execution so the synthetic spectrum matches the stored convergence point.
    rows = min(model._fparm.shape[0], SAMPL4_FINAL_FPARM.shape[0])
    cols = min(model._fparm.shape[1], SAMPL4_FINAL_FPARM.shape[1])
    model._fparm[:rows, :cols] = SAMPL4_FINAL_FPARM[:rows, :cols]
    rows = min(model._iparm.shape[0], SAMPL4_FINAL_IPARM.shape[0])
    cols = min(model._iparm.shape[1], SAMPL4_FINAL_IPARM.shape[1])
    model._iparm[:rows, :cols] = SAMPL4_FINAL_IPARM[:rows, :cols]
    model.set_spectral_state(
        sb0=SAMPL4_FINAL_SB0,
        srng=SAMPL4_FINAL_SRNG,
        ishift=SAMPL4_FINAL_ISHFT,
        shift=SAMPL4_FINAL_SHFT,
        normalize_flags=SAMPL4_FINAL_NRMLZ,
    )
    model['weights'] = SAMPL4_FINAL_WEIGHTS

    site_spectra = model.current_spectrum
    weights_matrix = np.array(model['weights'], copy=True)

    count = data_slice.stop - data_slice.start
    simulated_total = np.dot(np.atleast_2d(weights_matrix), site_spectra[:, data_slice])
    simulated_total = np.squeeze(simulated_total)
    experimental = SAMPL4_INTENSITIES[:count]
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
