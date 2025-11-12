import math
from pathlib import Path

import nlsl

from .runfile_helpers import run_runfile


TESTS_DIR = Path(__file__).resolve().parent

RUNFILE_NAME = "sampl1.run"
DATA_NAME = "sampl1.dat"
FINAL_TOKENS_SAMPL1 = ("gib0", "rprp", "rpll")

RUNFILE_NAME_2 = "sampl2.run"
DATA_SERIES_2 = ("sampl200.dat", "sampl290.dat")
FINAL_TOKENS_SAMPL2 = ("gib0", "gib2", "rbar", "n", "betad", "c20", "c22")

RUNFILE_NAME_3 = "sampl3.run"
DATA_NAME_3 = "sampl3.dat"
FINAL_TOKENS_SAMPL3 = ("gib0", "gib2", "rbar", "c20")


def run_sampl1_pythonic():
    """Reproduce the ``sampl1`` runfile purely through the Python API."""

    model = nlsl.nlsl()

    # sampl1.run: call csl.par
    # csl.par: let gxx,gyy,gzz=2.0089,2.0021,2.0058
    # csl.par: let in2=2
    # csl.par: let axx,ayy,azz=5.6,33.8,5.3
    # csl.par: let betad=15
    # sampl1.run: let rpll, rprp = log(1.0e8), 8.0
    # sampl1.run: let gib0 = 1.5
    # sampl1.run: let lemx,lomx,kmx,mmx=6,5,4,2,2
    model.update(
        {
            "in2": 2,
            "gxx": 2.0089,
            "gyy": 2.0021,
            "gzz": 2.0058,
            "axx": 5.6,
            "ayy": 33.8,
            "azz": 5.3,
            "betad": 15,
            "lemx": 6,
            "lomx": 5,
            "kmx": 4,
            "mmx": (2, 2),
            "rpll": math.log10(1.0e8),
            "rprp": 8.0,
            "gib0": 1.5,
        }
    )

    # sampl1.run: data sampl1 ascii nspline 200 bc 20 shift
    data_path = TESTS_DIR / DATA_NAME
    model.load_data(
        data_path,
        nspline=200,
        bc_points=20,
        shift=True,
        normalize=False,
        derivative_mode=1,
    )

    # sampl1.run: vary rpll, rprp, gib0
    for token in ("rpll", "rprp", "gib0"):
        model.fit_params.vary[token] = True

    # sampl1.run: fit maxit 40 maxfun 1000 ftol 1e-3 xtol 1e-3
    for key, value in {"maxitr": 40, "maxfun": 1000, "ftol": 1e-3, "xtol": 1e-3}.items():
        model.fit_params[key] = value

    model.fit()
    return model


def run_legacy_runfile(runfile_name):
    """Feed a legacy runfile to ``procline`` and discard its log."""

    model = nlsl.nlsl()
    run_runfile(model, TESTS_DIR / runfile_name)
    log_path = TESTS_DIR / (runfile_name.replace(".run", "") + ".log")
    if log_path.exists():
        log_path.unlink()
    return model


def run_sampl1_procline():
    """Execute the original ``sampl1`` script through ``procline``."""

    return run_legacy_runfile(RUNFILE_NAME)


def run_sampl2_pythonic():
    """Mirror the ``sampl2`` workflow via the Python API."""

    model = nlsl.nlsl()

    # sampl2.run: let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2
    # sampl2.run: call 5pc.par
    model.update(
        {
            "lemx": 6,
            "lomx": 5,
            "kmx": 4,
            "mmx": 4,
            "ipnmx": 2,
            "in2": 2,
            "gxx": 2.0089,
            "gyy": 2.0058,
            "gzz": 2.0021,
            "axx": 4.9,
            "ayy": 4.9,
            "azz": 33.0,
        }
    )

    # sampl2.run: series psi = 0, 90
    core = nlsl.fortrancore
    core.parcom.nser = 2
    core.parcom.serval[:2] = (0.0, 90.0)

    # sampl2.run: data sampl200 ascii bcmode 20 nspline 200 / data sampl290
    for data_name in DATA_SERIES_2:
        model.load_data(
            TESTS_DIR / data_name,
            nspline=200,
            bc_points=20,
            shift=True,
            normalize=False,
            derivative_mode=1,
        )

    # sampl2.run: let b0 = 3400; let rbar, n = 7.5, 0.5; let c20 = 2; let gib0 = 1
    model.update(
        {
            "b0": 3400.0,
            "rbar": 7.5,
            "n": 0.5,
            "c20": 2.0,
            "gib0": 1.0,
            "gib2": 0.0,
            "c22": 0.0,
            "betad": 0.0,
        }
    )

    # sampl2.run: fit controls are left at their defaults; replicate the staged
    # variation blocks with the Python mapping interface.
    for token in ("gib0", "rbar", "c20"):
        model.fit_params.vary[token] = True
    model.fit()

    for token in ("gib2", "n", "c22", "betad"):
        model.fit_params.vary[token] = True
    model.fit()
    model.fit()

    return model


def run_sampl2_procline():
    """Execute the original ``sampl2`` script through ``procline``."""

    return run_legacy_runfile(RUNFILE_NAME_2)


def run_sampl3_pythonic():
    """Reproduce the ``sampl3`` MOMD fit via the Python bindings."""

    model = nlsl.nlsl()

    # sampl3.run: call 5pc.par; let nort=20; let rbar = 8.0; let c20 = 1.0; let gib0 = 1.5
    model.update(
        {
            "in2": 2,
            "gxx": 2.0089,
            "gyy": 2.0058,
            "gzz": 2.0021,
            "axx": 4.9,
            "ayy": 4.9,
            "azz": 33.0,
            "nort": 20,
            "rbar": 8.0,
            "c20": 1.0,
            "gib0": 1.5,
            "gib2": 0.0,
        }
    )

    # sampl3.run: data sampl3 ascii nspline 400 bc 20 shift (the example keeps
    # the baseline correction and shifting logic active).
    model.load_data(
        TESTS_DIR / DATA_NAME_3,
        nspline=400,
        bc_points=20,
        shift=True,
        normalize=True,
        derivative_mode=1,
    )

    for key, value in {"maxitr": 40, "maxfun": 400, "ftol": 1.0e-3, "xtol": 1.0e-3}.items():
        model.fit_params[key] = value

    # sampl3.run: vary rpll, rprp, c20, gib0
    for token in ("rpll", "rprp", "c20", "gib0"):
        model.fit_params.vary[token] = True
    model.fit()

    # sampl3.run: fix rpll, rprp; vary rbar (after converting the tensor form).
    for token in ("rpll", "rprp"):
        model.fit_params.vary[token] = False
    model.fit_params.vary["rbar"] = True
    model.fit()

    # sampl3.run: vary gib2 and perform the final refinement.
    model.fit_params.vary["gib2"] = True
    model.fit()

    return model


def run_sampl3_procline():
    """Execute the original ``sampl3`` script through ``procline``."""

    return run_legacy_runfile(RUNFILE_NAME_3)


def capture_final_parameters(model, tokens):
    """Collect the converged fit parameters for comparison."""

    return {token: float(model[token]) for token in tokens}


def test_sampl1_pythonic_matches_procline():
    """The pythonic and runfile workflows should converge to the same fit."""

    pythonic_model = run_sampl1_pythonic()
    procline_model = run_sampl1_procline()

    pythonic_params = capture_final_parameters(pythonic_model, FINAL_TOKENS_SAMPL1)
    procline_params = capture_final_parameters(procline_model, FINAL_TOKENS_SAMPL1)

    for token in FINAL_TOKENS_SAMPL1:
        assert math.isclose(
            pythonic_params[token],
            procline_params[token],
            rel_tol=1e-6,
            abs_tol=1e-8,
        )


def test_sampl2_pythonic_matches_procline():
    """Ensure the ``sampl2`` pythonic port reaches the legacy fit."""

    pythonic_model = run_sampl2_pythonic()
    procline_model = run_sampl2_procline()

    pythonic_params = capture_final_parameters(pythonic_model, FINAL_TOKENS_SAMPL2)
    procline_params = capture_final_parameters(procline_model, FINAL_TOKENS_SAMPL2)

    for token in FINAL_TOKENS_SAMPL2:
        assert math.isclose(
            pythonic_params[token],
            procline_params[token],
            rel_tol=1e-6,
            abs_tol=1e-8,
        )


def test_sampl3_pythonic_matches_procline():
    """Validate that the ``sampl3`` pythonic workflow matches ``procline``."""

    pythonic_model = run_sampl3_pythonic()
    procline_model = run_sampl3_procline()

    pythonic_params = capture_final_parameters(pythonic_model, FINAL_TOKENS_SAMPL3)
    procline_params = capture_final_parameters(procline_model, FINAL_TOKENS_SAMPL3)

    for token in FINAL_TOKENS_SAMPL3:
        assert math.isclose(
            pythonic_params[token],
            procline_params[token],
            rel_tol=1e-6,
            abs_tol=1e-8,
        )
