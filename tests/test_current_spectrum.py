import numpy as np

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_FINAL_PARAMETERS,
    SAMPL4_POINT_COUNT,
    SAMPL4_INTENSITIES,
    apply_sampl4_final_state,
)


def test_generate_coordinates_enables_current_spectrum():
    model = nlsl.nlsl()
    model['nsite'] = 2
    model.update(SAMPL4_FINAL_PARAMETERS)

    # Generate the field grid used by the SAMPL4 data so the evaluation spans
    # the same axis as the published runfile.
    index, data_slice = model.generate_coordinates(
        SAMPL4_POINT_COUNT,
        start=SAMPL4_FIELD_START,
        step=SAMPL4_FIELD_STEP,
        derivative_mode=DERIVATIVE_MODE,
        baseline_points=BASELINE_EDGE_POINTS,
        normalize=False,
        nspline=NSPLINE_POINTS,
        shift=True,
        label='sampl4-single-eval',
        reset=True,
    )

    assert index == 0
    assert data_slice.start == 0 and data_slice.stop == SAMPL4_POINT_COUNT

    # Copy the processed intensities so the synthetic spectrum can be compared
    # directly against the reference trace.
    model.set_data(data_slice, SAMPL4_INTENSITIES[:SAMPL4_POINT_COUNT])

    # Store the final runfile-4 solution without touching the legacy front-end.
    apply_sampl4_final_state(model)

    site_spectra = model.current_spectrum
    weights_matrix = np.array(model['weights'], copy=True)

    assert site_spectra.shape == (2, SAMPL4_POINT_COUNT)
    assert np.all(np.isfinite(site_spectra))
    assert np.all(np.isfinite(weights_matrix))

    simulated_total = np.dot(np.atleast_2d(weights_matrix), site_spectra[:, data_slice])
    rel_rms = np.linalg.norm(np.squeeze(simulated_total) - SAMPL4_INTENSITIES[:SAMPL4_POINT_COUNT])
    rel_rms /= np.linalg.norm(SAMPL4_INTENSITIES[:SAMPL4_POINT_COUNT])

    assert rel_rms < 0.0401
