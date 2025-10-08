import numpy as np
from pathlib import Path

# Reuse the runfile-4 setup across the regression tests.
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

SAMPL4_DATA_PATH = Path(__file__).resolve().parent / "sampl4.dat"

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
    "rx": np.array([7.14177897, 7.8396974]),
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
    "b0": np.array([3400.50251256, 0.0]),
    "gamman": 0.0,
    "cgtol": 0.001,
    "shiftr": 0.001,
    "shifti": 0.0,
    "range": np.array([100.0, 0.0]),
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
}

# Site populations extracted from the ``sampl4`` runfile output.
SAMPL4_FINAL_WEIGHTS = np.array([0.7155313, 0.2848810])

# Parameters required to recreate the converged spectrum without invoking the
# optimiser again. These are sourced from the fitted model and assigned to a
# fresh instance before calling ``current_spectrum``.
SAMPL4_SPECTRAL_KEYS = [
    "phase",
    "gib0",
    "gib2",
    "wxx",
    "wyy",
    "wzz",
    "gxx",
    "gyy",
    "gzz",
    "axx",
    "ayy",
    "azz",
    "rx",
    "ry",
    "rz",
    "pml",
    "pmxy",
    "pmzz",
    "djf",
    "djfprp",
    "oss",
    "psi",
    "alphad",
    "betad",
    "gammad",
    "alpham",
    "betam",
    "gammam",
    "c20",
    "c22",
    "c40",
    "c42",
    "c44",
    "lb",
    "dc20",
    "b0",
    "gamman",
    "cgtol",
    "shiftr",
    "shifti",
    "range",
    "in2",
    "ipdf",
    "ist",
    "ml",
    "mxy",
    "mzz",
    "lemx",
    "lomx",
    "kmn",
    "kmx",
    "mmn",
    "mmx",
    "ipnmx",
    "nort",
    "nstep",
    "nfield",
    "ideriv",
]

# ``procline`` statements that recreate the SAMPL4 best-fit parameters without
# invoking ``fit`` again.
SAMPL4_FINAL_COMMANDS = [
    "let phase = 0.0",
    "let gxx,gyy,gzz = 2.0089,2.0063,2.0021",
    "let axx,ayy,azz = 5.0,5.0,33.0",
    "let gib0 = 1.9962762455456595",
    "let rbar(1) = 7.14177909",
    "let rbar(2) = 7.83969762",
    "let lemx,lomx,kmx,mmx,ipnmx=12,10,7,7,2",
    "let in2 = 2",
    "let b0 = 3400.50251256,0.0",
    "let range = 100.0,0.0",
    "let shiftr = 0.001",
    "let shifti = 0.0",
    "let wxx,wyy,wzz = 0.0,0.0,0.0",
]
