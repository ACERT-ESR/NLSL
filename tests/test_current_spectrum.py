import numpy as np

import nlsl


def test_generate_coordinates_enables_current_spectrum():
    # Use the published SAMPL4 run-4 fit parameters so that the spectral
    # evaluation exercise matches a known configuration with two sites.
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

    # Prepare a synthetic field grid that mirrors the processed data from the
    # sample 4 run so that the theoretical evaluation produces the site
    # spectra on that axis.
    index, data_slice = model.generate_coordinates(
        200,
        start=3350.0,
        step=0.5025125628139904,
        derivative_mode=1,
        baseline_points=20,
        normalize=True,
        nspline=200,
        shift=True,
        label='sampl4-single-eval',
        reset=True,
    )

    assert index == 0
    assert data_slice.start == 0 and data_slice.stop == 200

    # ``current_spectrum`` should now be able to call into the Fortran core and
    # fill both the component spectra and the associated weights with finite
    # floating point numbers.
    site_spectra, weights = model.current_spectrum

    assert site_spectra.shape == (2, 200)
    assert np.all(np.isfinite(site_spectra))
    assert weights.shape == (1, 2)
    assert np.all(np.isfinite(weights))
