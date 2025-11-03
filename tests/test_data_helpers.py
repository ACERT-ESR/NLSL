import os

import numpy as np

import nlsl

from nlsl.data import fit_linear_baseline, natural_cubic_spline_resample


def test_fit_linear_baseline_recovers_line():
    x = np.linspace(0.0, 9.0, 10)
    y = 2.0 + 0.5 * x
    corrected, intercept, slope, noise = fit_linear_baseline(x, y, edge_points=4)

    assert np.allclose(corrected, np.zeros_like(x))
    assert np.isclose(intercept, 2.0)
    assert np.isclose(slope, 0.5)
    assert np.isclose(noise, 0.0)


def test_natural_cubic_spline_resample_matches_fortran(tmp_path):
    x = np.linspace(-2.0, 2.0, 9)
    y = x**3 - 2.0 * x + 1.0
    nspline = 17
    data_path = tmp_path / "sample.dat"
    with data_path.open("w", encoding="utf-8") as handle:
        for xi, yi in zip(x, y):
            handle.write(f"{xi: .8f} {yi: .8f}\n")

    fortran = nlsl.nlsl()
    nlsl.fortrancore.parcom.nser = 1
    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        fortran.procline(
            f"data sample ascii nspline {nspline} bc 0 deriv 1 noshift nonorm"
        )
    finally:
        os.chdir(cwd)
    expdat = nlsl.fortrancore.expdat
    ix0 = int(expdat.ixsp[0]) - 1
    count = int(expdat.npts[0])
    fortran_y = expdat.data[ix0 : ix0 + count].copy()

    _, python_y = natural_cubic_spline_resample(x, y, nspline)
    assert np.allclose(python_y, fortran_y)


def _capture_state():
    fc = nlsl.fortrancore
    expdat = fc.expdat
    mspctr = fc.mspctr
    lmcom = fc.lmcom

    state = {
        "data": expdat.data.copy(),
        "iform": expdat.iform.copy(),
        "ibase": expdat.ibase.copy(),
        "nft": expdat.nft.copy(),
        "rmsn": expdat.rmsn.copy(),
        "dataid": expdat.dataid.copy(),
        "wndoid": expdat.wndoid.copy(),
        "fvec": lmcom.fvec.copy(),
        "npts": expdat.npts.copy(),
        "ixsp": expdat.ixsp.copy(),
        "ishft": expdat.ishft.copy(),
        "idrv": expdat.idrv.copy(),
        "nrmlz": expdat.nrmlz.copy(),
        "sbi": expdat.sbi.copy(),
        "sdb": expdat.sdb.copy(),
        "srng": expdat.srng.copy(),
        "shft": expdat.shft.copy(),
        "tmpshft": expdat.tmpshft.copy(),
        "slb": expdat.slb.copy(),
        "sb0": expdat.sb0.copy(),
        "sphs": expdat.sphs.copy(),
        "spsi": expdat.spsi.copy(),
        "nspc": int(expdat.nspc),
        "ndatot": int(expdat.ndatot),
        "ishglb": int(expdat.ishglb),
        "normflg": int(expdat.normflg),
        "shftflg": int(expdat.shftflg),
        "bcmode": int(expdat.bcmode),
        "drmode": int(expdat.drmode),
        "nspline": int(expdat.nspline),
        "inform": int(expdat.inform),
    }
    for name in ("spectr", "wspec", "sfac"):
        if hasattr(mspctr, name):
            state[name] = getattr(mspctr, name).copy()
    return state
def test_load_data_matches_datac():
    samples = ["tests/sampl1", "tests/sampl3"]

    # Legacy path
    legacy = nlsl.nlsl()
    for sample in samples:
        cmd = f"data {sample} ascii nspline 200 bc 20 shift norm"
        legacy.procline(cmd)
    legacy_state = _capture_state()

    # Python path
    modern = nlsl.nlsl()
    for sample in samples:
        modern.load_data(sample, nspline=200, bc_points=20, shift=True, normalize=True)
    python_state = _capture_state()

    for key, legacy_value in legacy_state.items():
        python_value = python_state[key]
        if isinstance(legacy_value, np.ndarray):
            if np.issubdtype(legacy_value.dtype, np.floating):
                assert np.allclose(legacy_value, python_value)
            else:
                assert np.array_equal(legacy_value, python_value)
        else:
            assert legacy_value == python_value
