import numpy as np


def test_pythonic_sampl4_fit_matches_data(sampl4_fit_result):
    """The refactored example should reach the published residual."""

    assert sampl4_fit_result["rel_rms"] < 0.0401


def test_current_spectrum_matches_fit_components(sampl4_fit_result):
    """A post-fit single-point evaluation must reproduce the fit spectra."""

    site_spectra_cs, weights_cs = sampl4_fit_result["model"].current_spectrum

    assert np.allclose(site_spectra_cs, sampl4_fit_result["site_spectra"], atol=5e-6)
    assert np.allclose(weights_cs, sampl4_fit_result["weights"])

    recomputed = np.sum(
        weights_cs[0][:, np.newaxis] * site_spectra_cs[:, sampl4_fit_result["data_slice"]],
        axis=0,
    )

    assert np.allclose(recomputed, sampl4_fit_result["simulated_total"], atol=3e-6)
