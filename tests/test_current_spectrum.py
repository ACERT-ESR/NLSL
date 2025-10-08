import numpy as np

import nlsl
from test_sampl4_fit import run_pythonic_sampl4_fit


def test_generate_coordinates_enables_current_spectrum():
    # Capture the runfile-4 fit through the Python API so the reference weights
    # and experimental trace match the published example.
    sampl4_fit_result = run_pythonic_sampl4_fit()

    final_params = {
        'phase': 0.0,
        'gib0': 1.9962757195220067,
        'gib2': 0.0,
        'wxx': 0.0,
        'wyy': 0.0,
        'wzz': 0.0,
        'gxx': 2.0089,
        'gyy': 2.0063,
        'gzz': 2.0021,
        'axx': 5.0,
        'ayy': 5.0,
        'azz': 33.0,
        'rx': np.array([7.14177897, 7.8396974]),
        'ry': 0.0,
        'rz': 0.0,
        'pml': 0.0,
        'pmxy': 0.0,
        'pmzz': 0.0,
        'djf': 0.0,
        'djfprp': 0.0,
        'oss': 0.0,
        'psi': 0.0,
        'alphad': 0.0,
        'betad': 0.0,
        'gammad': 0.0,
        'alpham': 0.0,
        'betam': 0.0,
        'gammam': 0.0,
        'c20': 0.0,
        'c22': 0.0,
        'c40': 0.0,
        'c42': 0.0,
        'c44': 0.0,
        'lb': 0.0,
        'dc20': 0.0,
        'b0': np.array([3400.50251256, 0.0]),
        'gamman': 0.0,
        'cgtol': 0.001,
        'shiftr': 0.001,
        'shifti': 0.0,
        'range': np.array([100.0, 0.0]),
        'in2': 2,
        'ipdf': 0,
        'ist': 0,
        'ml': 0,
        'mxy': 0,
        'mzz': 0,
        'lemx': 12,
        'lomx': 10,
        'kmn': 0,
        'kmx': 7,
        'mmn': 0,
        'mmx': 7,
        'ipnmx': 2,
        'nort': 0,
        'nstep': 0,
        'nfield': 200,
        'ideriv': 1,
    }

    model = nlsl.nlsl()
    model['nsite'] = 2
    model.update(final_params)

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
    nlsl.fortrancore.mspctr.sfac[: sampl4_fit_result["weights"].shape[1], index] = sampl4_fit_result["weights"][0]

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
