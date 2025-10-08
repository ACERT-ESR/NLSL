import numpy as np

import nlsl
from tests.test_sampl4_fit import run_pythonic_sampl4_fit
from tests.sampl4_reference import SAMPL4_FINAL_PARAMETERS, SAMPL4_FINAL_WEIGHTS


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
        200,
        start=3350.0,
        step=0.5025125628139904,
        derivative_mode=1,
        baseline_points=20,
        normalize=False,
        nspline=200,
        shift=True,
        label='sampl4-single-eval',
        reset=True,
    )

    assert index == 0
    assert data_slice.start == 0 and data_slice.stop == 200

    # Store the final site populations so the theoretical spectra reflect the
    # fitted composition from the runfile reproduction.
    nlsl.fortrancore.mspctr.sfac[: SAMPL4_FINAL_WEIGHTS.size, index] = SAMPL4_FINAL_WEIGHTS

    site_spectra, weights = model.current_spectrum

    assert site_spectra.shape == (2, 200)
    assert np.all(np.isfinite(site_spectra))
    assert weights.shape == (1, 2)
    assert np.all(np.isfinite(weights))
    assert np.allclose(weights[0], sampl4_fit_result["weights"][0])

    # Confirm that the component spectra stored during the Pythonic fit combine
    # to reproduce the experimental trace within the published tolerance.
    fit_total = sampl4_fit_result["weights"][0] @ sampl4_fit_result["site_spectra"][:, sampl4_fit_result["data_slice"]]
    residual = fit_total - sampl4_fit_result["experimental"]
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(sampl4_fit_result["experimental"])

    assert rel_rms < 0.0401
