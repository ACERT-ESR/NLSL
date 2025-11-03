import numpy as np

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
    SAMPL4_POINT_COUNT,
    SAMPL4_SPECTRAL_DATA,
)


def test_generate_coordinates_enables_current_spectrum():
    model = nlsl.nlsl()
    model["nsite"] = 2
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
        label="sampl4-single-eval",
        reset=True,
    )

    assert index == 0
    assert data_slice.start == 0 and data_slice.stop == SAMPL4_POINT_COUNT

    # Copy the processed intensities so the synthetic spectrum can be compared
    # directly against the reference trace.
    model.set_data(data_slice, SAMPL4_SPECTRAL_DATA[:SAMPL4_POINT_COUNT])

    # Mirror the runfile-4 solution through the dictionary interface so the
    # synthetic spectrum is generated with the converged parameters.
    model.update(SAMPL4_FINAL_PARAMETERS)
    model["sb0"] = SAMPL4_FINAL_SB0
    model["srng"] = SAMPL4_FINAL_SRNG
    model["ishft"] = SAMPL4_FINAL_ISHFT
    model["shft"] = SAMPL4_FINAL_SHFT
    model["nrmlz"] = SAMPL4_FINAL_NRMLZ
    model.weights = SAMPL4_FINAL_WEIGHTS

    site_spectra = model.current_spectrum
    assert site_spectra.shape == (2, SAMPL4_POINT_COUNT)
    assert np.all(np.isfinite(site_spectra))
    assert np.all(np.isfinite(model.weights))

    simulated_total = np.squeeze(model.weights @ site_spectra)
    rel_rms = np.linalg.norm(
        np.squeeze(simulated_total) - SAMPL4_SPECTRAL_DATA[:SAMPL4_POINT_COUNT]
    )
    rel_rms /= np.linalg.norm(SAMPL4_SPECTRAL_DATA[:SAMPL4_POINT_COUNT])

    assert rel_rms < 0.0401
