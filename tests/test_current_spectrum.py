import numpy as np

import nlsl
from tests.test_sampl4_fit import run_pythonic_sampl4_fit
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_FIELD_START,
    SAMPL4_FIELD_STEP,
    SAMPL4_FINAL_PARAMETERS,
    SAMPL4_FINAL_WEIGHTS,
    SAMPL4_POINT_COUNT,
)


def test_generate_coordinates_enables_current_spectrum():
    # Capture the runfile-4 fit through the Python API so the reference weights
    # and experimental trace match the published example.
    sampl4_fit_result = run_pythonic_sampl4_fit()

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

    # Store the final site populations so the theoretical spectra reflect the
    # fitted composition from the runfile reproduction.
    nlsl.fortrancore.mspctr.sfac[:, index] = 0.0
    nlsl.fortrancore.mspctr.sfac[: SAMPL4_FINAL_WEIGHTS.size, index] = SAMPL4_FINAL_WEIGHTS

    site_spectra, weights = model.current_spectrum

    assert site_spectra.shape == (2, SAMPL4_POINT_COUNT)
    assert np.all(np.isfinite(site_spectra))
    assert weights.shape == (2,)
    assert np.all(np.isfinite(weights))
    assert np.allclose(weights, sampl4_fit_result["weights"])

    # Confirm that the component spectra stored during the Pythonic fit combine
    # to reproduce the experimental trace within the published tolerance.
    fit_total = np.dot(
        sampl4_fit_result["weights"],
        sampl4_fit_result["site_spectra"][:, sampl4_fit_result["data_slice"]],
    )
    residual = fit_total - sampl4_fit_result["experimental"]
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(sampl4_fit_result["experimental"])

    assert rel_rms < 0.0401
