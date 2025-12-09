import math
import os
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
        "b0",
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
        "b0",
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
            # sampl1.run: b0 noted as 3400 G in runfile header
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
                    "b0": 3400.0,
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

            original_dir = os.getcwd()
            os.chdir(TESTS_DIR)

            try:
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

                # sampl2.run: let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2
                model.procline("let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2")

                # sampl2.run: series psi = 0, 90
                model.procline("series psi = 0, 90")

                # sampl2.run: data sampl200 ascii bcmode 20 nspline 200 / data sampl290
                model.procline("data sampl200 ascii bcmode 20 nspline 200")
                model.procline("data sampl290")

                # sampl2.run: let b0 = 3400; let rbar, n = 7.5, 0.5; let c20 = 2; let gib0 = 1
                model.procline("let b0 = 3400")
                model.procline("let rbar, n = 7.5, 0.5")
                model.procline("let c20 = 2")
                model.procline("let gib0 = 1")

                # sampl2.run: refine the remaining linewidth and ordering coefficients
                model.procline("let gib2 = 0")
                model.procline("let c22 = 0")
                model.procline("let betad = 0")

                # sampl2.run: search rbar / search betad step 5 bound 45 / search c20 / search c22 / search gib0
                for command in (
                    "search rbar",
                    "search betad step 5 bound 45",
                    "search c20",
                    "search c22",
                    "search gib0",
                ):
                    model.procline(command)

                # sampl2.run: vary gib0, rbar, c20 / fit
                model.procline("vary gib0, rbar, c20")
                model.procline("fit")

                # sampl2.run: vary gib2, n, c22, betad / fit / fit
                model.procline("vary gib2, n, c22, betad")
                model.procline("fit")
                model.procline("fit")
            finally:
                os.chdir(original_dir)
            return model

        case "sampl3":
            model = nlsl.nlsl()
            run_runfile(model, TESTS_DIR / "sampl3.run")
            log_path = TESTS_DIR / "sampl3.log"
            if log_path.exists():
                log_path.unlink()
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


def test_pythonic_runs_match_procline():
    """All runfile ports should match every reported parameter."""

    for runfile_name in ("sampl1", "sampl2", "sampl3"):
        # Capture parameter values for both the pythonic and legacy runs in one place
        parameter_sets = {}
        for label, runner in (("pythonic", run_pythonic), ("procline", run_legacy_runfile)):
            model = runner(runfile_name)
            parameter_sets[label] = {}
            for token in RUNFILE_TOKENS[runfile_name]:
                if hasattr(model[token], "__len__") and not isinstance(model[token], (str, bytes)):
                    parameter_sets[label][token] = tuple(float(entry) for entry in model[token])
                else:
                    parameter_sets[label][token] = float(model[token])

        for token in RUNFILE_TOKENS[runfile_name]:
            if isinstance(parameter_sets["pythonic"][token], tuple):
                for index in range(len(parameter_sets["pythonic"][token])):
                    assert math.isclose(
                        parameter_sets["pythonic"][token][index],
                        parameter_sets["procline"][token][index],
                        rel_tol=1e-6,
                        abs_tol=1e-8,
                    )
            else:
                assert math.isclose(
                    parameter_sets["pythonic"][token],
                    parameter_sets["procline"][token],
                    rel_tol=1e-6,
                    abs_tol=1e-8,
                )
