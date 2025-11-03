import numpy as np
from pathlib import Path

from nlsl.data import process_spectrum

# Reuse the runfile-4 setup across the regression tests.
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

SAMPL4_DATA_PATH = Path(__file__).resolve().parent / "sampl4.dat"

# Process the experimental trace with the same resampling parameters the
# classic runfile uses so the tests operate on the interpolated 200-point grid
# instead of the raw 256-point ASCII data.
_SAMPL4_PROCESSED = process_spectrum(
    SAMPL4_DATA_PATH,
    NSPLINE_POINTS,
    BASELINE_EDGE_POINTS,
    derivative_mode=DERIVATIVE_MODE,
    normalize=False,
)

SAMPL4_FIELD_START = float(_SAMPL4_PROCESSED.start)
SAMPL4_FIELD_STEP = float(_SAMPL4_PROCESSED.step)
SAMPL4_POINT_COUNT = _SAMPL4_PROCESSED.y.size
SAMPL4_SPECTRAL_DATA = _SAMPL4_PROCESSED.y.copy()

# Starting guesses copied from the ``let`` statements in ``sampl4.run``.
SAMPL4_INITIAL_PARAMETERS = {
    "nsite": 2,
    "in2": 2,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "lemx": 12,
    "lomx": 10,
    "kmx": 7,
    "mmx": 7,
    "ipnmx": 2,
    "gib0": 0.5,
    "rx": np.array([np.log10(3.0e8), np.log10(1.0e7)]),
}

# MINPACK control parameters recorded in ``sampl4.run``.
SAMPL4_FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
    "ftol": 1.0e-2,
    "gtol": 1.0e-6,
    "xtol": 1.0e-4,
}

# The procline tokens that release the parameters prior to fitting.
SAMPL4_PARAMETERS_TO_VARY = ["gib0", "rbar(1)", "rbar(2)"]

# Final parameters reported by the classic runfile after convergence.
SAMPL4_FINAL_PARAMETERS = {
    "nsite": 2,
    "phase": 0.0,
    "gib0": 1.9962757195220067,
    "gib2": 0.0,
    "wxx": 0.0,
    "wyy": 0.0,
    "wzz": 0.0,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "rx": np.array([7.8396974, 7.14177897]),
    "ry": 0.0,
    "rz": 0.0,
    "pml": 0.0,
    "pmxy": 0.0,
    "pmzz": 0.0,
    "djf": 0.0,
    "djfprp": 0.0,
    "oss": 0.0,
    "psi": 0.0,
    "alphad": 0.0,
    "betad": 0.0,
    "gammad": 0.0,
    "alpham": 0.0,
    "betam": 0.0,
    "gammam": 0.0,
    "c20": 0.0,
    "c22": 0.0,
    "c40": 0.0,
    "c42": 0.0,
    "c44": 0.0,
    "lb": 0.0,
    "dc20": 0.0,
    "b0": 3400.50251256,
    "fldi": 3350.5046000757857,
    "dfld": SAMPL4_FIELD_STEP,
    "gamman": 0.0,
    "cgtol": 0.001,
    "shiftr": 0.001,
    "shifti": 0.0,
    "range": 100.0,
    "in2": 2,
    "ipdf": 0,
    "ist": 0,
    "ml": 0,
    "mxy": 0,
    "mzz": 0,
    "lemx": 12,
    "lomx": 10,
    "kmn": 0,
    "kmx": 7,
    "mmn": 0,
    "mmx": 7,
    "ipnmx": 2,
    "nort": 0,
    "nstep": 0,
    "nfield": 200,
    "ideriv": 1,
    "iwflg": 0,
    "igflg": 0,
    "iaflg": 0,
    "jkmn": 0,
    "jmmn": 0,
    "irflg": 2,
    "ndim": 156,
}

# Site populations extracted from the ``sampl4`` runfile output.
SAMPL4_FINAL_WEIGHTS = np.array([0.2848810, 0.7155313])

# Spectral metadata required for reproducing the converged simulation.
SAMPL4_FINAL_SB0 = np.array([SAMPL4_FINAL_PARAMETERS["b0"]])
SAMPL4_FINAL_SRNG = np.array([SAMPL4_FINAL_PARAMETERS["range"]])
SAMPL4_FINAL_ISHFT = np.array([1], dtype=np.int32)
SAMPL4_FINAL_SHFT = np.array([0.0])
SAMPL4_FINAL_NRMLZ = np.array([0], dtype=np.int32)
