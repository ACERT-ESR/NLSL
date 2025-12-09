import math
from pathlib import Path

import nlsl

from .runfile_helpers import run_runfile


TESTS_DIR = Path(__file__).resolve().parent

# Each runfile lists parameters either directly or through referenced ``.par``
# files. All of these tokens must match between the pythonic reconstruction and
# the original ``procline`` workflow.
RUNFILE_TOKENS = {
    "sampl1": (
        "gxx",
        "gyy",
        "gzz",
        "in2",
        "axx",
        "ayy",
        "azz",
        "betad",
        "lemx",
        "lomx",
        "kmx",
        "mmx",
        "rpll",
        "rprp",
        "gib0",
    ),
    "sampl2": (
        "gxx",
        "gyy",
        "gzz",
        "in2",
        "axx",
        "ayy",
        "azz",
        "lemx",
        "lomx",
        "kmx",
        "mmx",
        "ipnmx",
        "b0",
        "rbar",
        "n",
        "c20",
        "gib0",
        "gib2",
        "c22",
        "betad",
    ),
    "sampl3": (
        "gxx",
        "gyy",
        "gzz",
        "in2",
        "axx",
        "ayy",
        "azz",
        "nort",
        "rbar",
        "c20",
        "gib0",
        "gib2",
        "rpll",
        "rprp",
    ),
}


def run_pythonic(runfile_name):
    """Recreate the requested runfile using the Python API."""

    match runfile_name:
        case "sampl1":
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
            model.load_data(
                TESTS_DIR / "sampl1.dat",
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

        case "sampl2":
            model = nlsl.nlsl()

            # sampl2.run: let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2
            model.update({"lemx": 6, "lomx": 5, "kmx": 4, "mmx": 4, "ipnmx": 2})

            # sampl2.run: series psi = 0, 90
            nlsl.fortrancore.parcom.nser = 2
            nlsl.fortrancore.parcom.serval[:2] = (0.0, 90.0)

            # sampl2.run: call 5pc.par
            # 5pc.par: let gxx,gyy,gzz=2.0089,2.0058,2.0021
            # 5pc.par: let in2=2
            # 5pc.par: let axx,ayy,azz=4.9,4.9,33.0
            model.update(
                {
                    "gxx": 2.0089,
                    "gyy": 2.0058,
                    "gzz": 2.0021,
                    "in2": 2,
                    "axx": 4.9,
                    "ayy": 4.9,
                    "azz": 33.0,
                }
            )

            # sampl2.run: data sampl200 ascii bcmode 20 nspline 200 / data sampl290
            for data_name in ("sampl200.dat", "sampl290.dat"):
                model.load_data(
                    TESTS_DIR / data_name,
                    nspline=200,
                    bc_points=20,
                    shift=True,
                    normalize=False,
                    derivative_mode=1,
                )

            # sampl2.run: let b0 = 3400; let rbar, n = 7.5, 0.5; let c20 = 2; let gib0 = 1
            model.update({"b0": 3400.0, "rbar": 7.5, "n": 0.5, "c20": 2.0, "gib0": 1.0})

            # sampl2.run: refine the remaining linewidth and ordering coefficients
            model.update({"gib2": 0.0, "c22": 0.0, "betad": 0.0})

            # sampl2.run: vary gib0, rbar, c20 / fit
            for token in ("gib0", "rbar", "c20"):
                model.fit_params.vary[token] = True
            model.fit()

            # sampl2.run: vary gib2, n, c22, betad / fit / fit
            for token in ("gib2", "n", "c22", "betad"):
                model.fit_params.vary[token] = True
            model.fit()
            model.fit()
            return model

        case "sampl3":
            model = nlsl.nlsl()

            # sampl3.run: call 5pc.par
            # 5pc.par: let gxx,gyy,gzz=2.0089,2.0058,2.0021
            # 5pc.par: let in2=2
            # 5pc.par: let axx,ayy,azz=4.9,4.9,33.0
            model.update(
                {
                    "gxx": 2.0089,
                    "gyy": 2.0058,
                    "gzz": 2.0021,
                    "in2": 2,
                    "axx": 4.9,
                    "ayy": 4.9,
                    "azz": 33.0,
                }
            )

            # sampl3.run: let nort=20; let rbar = 8.0; let c20 = 1.0; let gib0 = 1.5
            model.update({"nort": 20, "rbar": 8.0, "c20": 1.0, "gib0": 1.5, "gib2": 0.0})

            # sampl3.run: data sampl3 ascii nspline 400 bc 20 shift
            model.load_data(
                TESTS_DIR / "sampl3.dat",
                nspline=400,
                bc_points=20,
                shift=True,
                normalize=True,
                derivative_mode=1,
            )

            # sampl3.run: search rbar / search c20 / axial r
            model.procline("search rbar")
            model.procline("search c20")
            model.procline("axial r")

            # sampl3.run: vary rpll,rprp,c20,gib0
            for token in ("rpll", "rprp", "c20", "gib0"):
                model.fit_params.vary[token] = True

            # sampl3.run: fit maxit 40 maxfun 400 ftol 1e-3 xtol 1e-3
            for key, value in {"maxitr": 40, "maxfun": 400, "ftol": 1.0e-3, "xtol": 1.0e-3}.items():
                model.fit_params[key] = value
            model.fit()

            # sampl3.run: fix rpll / fix rprp / spherical r / vary rbar / fit
            for token in ("rpll", "rprp"):
                model.fit_params.vary[token] = False
            model.procline("spherical r")
            model.fit_params.vary["rbar"] = True
            model.fit()

            # sampl3.run: vary gib2 / fit
            model.fit_params.vary["gib2"] = True
            model.fit()
            return model

        case _:
            raise ValueError("Unsupported runfile requested")


def run_legacy_runfile(runfile_name):
    """Feed a legacy runfile to ``procline`` and discard its log."""

    model = nlsl.nlsl()
    run_runfile(model, TESTS_DIR / (runfile_name + ".run"))
    log_path = TESTS_DIR / (runfile_name + ".log")
    if log_path.exists():
        log_path.unlink()
    return model


def capture_final_parameters(model, tokens):
    """Collect the converged fit parameters for comparison."""

    values = {}
    for token in tokens:
        value = model[token]
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            values[token] = tuple(float(entry) for entry in value)
        else:
            values[token] = float(value)
    return values


def test_pythonic_runs_match_procline():
    """All runfile ports should match every reported parameter."""

    for runfile_name in ("sampl1", "sampl2", "sampl3"):
        pythonic_model = run_pythonic(runfile_name)
        procline_model = run_legacy_runfile(runfile_name)

        pythonic_params = capture_final_parameters(pythonic_model, RUNFILE_TOKENS[runfile_name])
        procline_params = capture_final_parameters(procline_model, RUNFILE_TOKENS[runfile_name])

        for token in RUNFILE_TOKENS[runfile_name]:
            if isinstance(pythonic_params[token], tuple):
                for index in range(len(pythonic_params[token])):
                    assert math.isclose(
                        pythonic_params[token][index],
                        procline_params[token][index],
                        rel_tol=1e-6,
                        abs_tol=1e-8,
                    )
            else:
                assert math.isclose(
                    pythonic_params[token],
                    procline_params[token],
                    rel_tol=1e-6,
                    abs_tol=1e-8,
                )
