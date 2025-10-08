import numpy as np
import pytest

import nlsl
from tests.sampl4_reference import (
    BASELINE_EDGE_POINTS,
    DERIVATIVE_MODE,
    NSPLINE_POINTS,
    SAMPL4_DATA_PATH,
    SAMPL4_FINAL_COMMANDS,
    SAMPL4_FINAL_WEIGHTS,
)


def test_sampl4_best_parameters_match_data_without_fit():
    """Copy the converged parameters into a fresh model and verify the residual."""

    model = nlsl.nlsl()
    model['nsite'] = 2

    model.load_data(
        SAMPL4_DATA_PATH,
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=False,
        derivative_mode=DERIVATIVE_MODE,
    )

    for command in SAMPL4_FINAL_COMMANDS:
        model.procline(command)

    nlsl.fortrancore.mspctr.sfac[:, :] = 0.0
    nlsl.fortrancore.mspctr.sfac[: SAMPL4_FINAL_WEIGHTS.size, 0] = SAMPL4_FINAL_WEIGHTS

    site_spectra, weights = model.current_spectrum

    layout = model.layout
    start_index = int(layout["ixsp"][0])
    point_count = int(layout["npts"][0])
    data_slice = slice(start_index, start_index + point_count)

    simulated_total = weights[0] @ site_spectra[:, data_slice]
    experimental = nlsl.fortrancore.expdat.data[data_slice].copy()
    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)

    assert rel_rms == pytest.approx(0.040096, abs=1e-4)
