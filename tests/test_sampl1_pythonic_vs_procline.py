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
    "sampl2a": (
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
    "sampl3a": (
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
    "sampl4": (
        "gxx",
        "gyy",
        "gzz",
        "in2",
        "axx",
        "ayy",
        "azz",
        "b0",
        "lemx",
        "lomx",
        "kmx",
        "mmx",
        "ipnmx",
        "rbar",
        "gib0",
    ),
    "sampl5": (
        "gxx",
        "gyy",
        "gzz",
        "in2",
        "axx",
        "ayy",
        "azz",
        "b0",
        "betad",
        "rbar",
        "n",
        "c20",
        "c22",
        "gib0",
        "gib2",
        "nort",
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

        case "sampl2" | "sampl2a":
            model = nlsl.nlsl()

            original_dir = os.getcwd()
            os.chdir(TESTS_DIR)

            try:
                # call 5pc.par: set magnetic parameters for the spin probe
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

                # let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2
                model.procline("let lemx,lomx,kmx,mmx,ipnmx = 6,5,4,4,2")

                # series psi = 0, 90
                model.procline("series psi = 0, 90")

                # pick the sampl200 loader that matches the requested runfile
                if runfile_name == "sampl2":
                    model.procline("data sampl200 ascii bcmode 20 nspline 200")
                else:
                    model.procline("data sampl200 ascii nspline 200")
                model.procline("data sampl290")

                # let b0 = 3400 / let rbar, n = 7.5, 0.5 / let c20 = 2 / let gib0 = 1
                model.procline("let b0 = 3400")
                model.procline("let rbar, n = 7.5, 0.5")
                model.procline("let c20 = 2")
                model.procline("let gib0 = 1")

                # refine linewidth and ordering coefficients before fitting
                model.procline("let gib2 = 0")
                model.procline("let c22 = 0")
                model.procline("let betad = 0")

                # search rbar / search betad step 5 bound 45 / search c20 / search c22 / search gib0
                for command in (
                    "search rbar",
                    "search betad step 5 bound 45",
                    "search c20",
                    "search c22",
                    "search gib0",
                ):
                    model.procline(command)

                # vary gib0, rbar, c20 / fit
                model.procline("vary gib0, rbar, c20")
                model.procline("fit")

                # vary gib2, n, c22, betad / fit / fit
                model.procline("vary gib2, n, c22, betad")
                model.procline("fit")
                model.procline("fit")
            finally:
                os.chdir(original_dir)
                log_path = TESTS_DIR / "sampl2.log"
                if log_path.exists():
                    log_path.unlink()
            return model

        case "sampl3" | "sampl3a":
            model = nlsl.nlsl()

            original_dir = os.getcwd()
            os.chdir(TESTS_DIR)

            try:
                # log sampl3
                model.procline("log sampl3")

                # call 5pc.par to set magnetic parameters
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

                # let nort=20 / let rbar = 8.0 / let c20 = 1.0 / let gib0 = 1.5
                model.procline("let nort=20")
                model.procline("let rbar = 8.0")
                model.procline("let c20 = 1.0")
                model.procline("let gib0 = 1.5")

                # basis sampl3
                model.procline("basis sampl3")

                # data sampl3 ascii ...
                if runfile_name == "sampl3":
                    model.procline("data sampl3 ascii nspline 400 bc 20 shift")
                else:
                    model.procline("data sampl3 ascii nspline 400 shift")

                # search rbar / search c20 / axial r
                model.procline("search rbar")
                model.procline("search c20")
                model.procline("axial r")

                # vary rpll,rprp,c20,gib0
                model.procline("vary rpll,rprp,c20,gib0")

                # fit maxit 40 maxfun 400 ftol 1e-3 xtol 1e-3
                model.procline("fit maxit 40 maxfun 400 ftol 1e-3 xtol 1e-3")

                # fix rpll / fix rprp / spherical r / vary rbar
                model.procline("fix rpll")
                model.procline("fix rprp")
                model.procline("spherical r")
                model.procline("vary rbar")

                # vary gib2 / fit
                model.procline("vary gib2")
                model.procline("fit")
            finally:
                os.chdir(original_dir)
                log_path = TESTS_DIR / "sampl3.log"
                if log_path.exists():
                    log_path.unlink()
            return model

        case "sampl4":
            model = nlsl.nlsl()

            original_dir = os.getcwd()
            os.chdir(TESTS_DIR)

            try:
                # sampl4.run: set a two-site model and magnetic parameters
                model.update(
                    {
                        "gxx": 2.0089,
                        "gyy": 2.0063,
                        "gzz": 2.0021,
                        "in2": 2,
                        "axx": 5.0,
                        "ayy": 5.0,
                        "azz": 33.0,
                    }
                )
                model.procline("sites 2")

                # sampl4.run: let lemx,lomx,kmx,mmx,ipnmx=12,10,7,7,2
                model.procline("let lemx,lomx,kmx,mmx,ipnmx=12,10,7,7,2")

                # sampl4.run: data sampl4 ascii nspline 200 bc 20 shift
                model.procline("data sampl4 ascii nspline 200 bc 20 shift")

                # sampl4.run: let rbar(1) = log(1.0e7) / let rbar(2) = log(3.0e8) / let gib0=0.5
                model.procline("let rbar(1) = log(1.0e7)")
                model.procline("let rbar(2) = log(3.0e8)")
                model.procline("let gib0=0.5")

                # sampl4.run: vary gib0, rbar(1), rbar(2) / fit maxit 40 maxfun 1000 ftol .01 / fit
                model.procline("vary gib0, rbar(1), rbar(2)")
                model.procline("fit maxit 40 maxfun 1000 ftol .01")
                model.procline("fit")
            finally:
                os.chdir(original_dir)
                log_path = TESTS_DIR / "sampl4.log"
                if log_path.exists():
                    log_path.unlink()
            return model

        case "sampl5":
            model = nlsl.nlsl()

            original_dir = os.getcwd()
            os.chdir(TESTS_DIR)

            try:
                # sampl5.run: mirror csl.par explicitly so each site starts with
                # non-zero g- and A-tensor elements.
                model.update(
                    {
                        "gxx": 2.0089,
                        "gyy": 2.0021,
                        "gzz": 2.0058,
                        "in2": 2,
                        "axx": 5.6,
                        "ayy": 33.8,
                        "azz": 5.3,
                        "betad": 15,
                    }
                )

                # sampl5.run: call csl.par for base magnetic values
                model.procline("call csl.par")

                # sampl5.run: define two sites and two tilt angles
                model.procline("sites 2")
                model.procline("series psi 0 90")

                # sampl5.run: initial estimates for motion, ordering, and linewidth
                model.procline("let b0=3405")
                model.procline("let rbar,n = 7.5,1")
                model.procline("let c20(1), c22(1) = 4, 1")
                model.procline("let c20(2), c22(2) = 0.2, 0")
                model.procline("let betad = 15")
                model.procline("let gib0=3")
                model.procline("let nort(2)=10")

                # sampl5.run: attempt to add a pruned basis set
                model.procline("basis sampl5")

                # sampl5.run: data sampl500 ascii nspline 200 bc 20 shift / data sampl590
                model.procline("data sampl500 ascii nspline 200 bc 20 shift")
                model.procline("data sampl590")

                # sampl5.run: vary gib0(*) rbar(*) / fit maxi 40 maxfun 1000
                model.procline("vary gib0(*) rbar(*)")
                model.procline("fit maxi 40 maxfun 1000")

                # sampl5.run: vary c20(*) / fit
                model.procline("vary c20(*)")
                model.procline("fit")

                # sampl5.run: vary n(*) c22(*) / fit
                model.procline("vary n(*) c22(*)")
                model.procline("fit")

                # sampl5.run: search gib2(1) / search gib2(2)
                model.procline("search gib2(1)")
                model.procline("search gib2(2)")

                # sampl5.run: search rbar(2) and refine the second site
                model.procline("search rbar(2)")
                model.procline("fix all")
                model.procline("vary rbar(2) n(2) c20(2) c22(2) gib0(2)")
                model.procline("fit")
            finally:
                os.chdir(original_dir)
                log_path = TESTS_DIR / "sampl5.log"
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
    # sampl2a.run and sampl3a.run reuse the base log names from sampl2.run and
    # sampl3.run respectively, so clean up the shared names as well.
    if runfile_name.endswith("a"):
        alternate_log = TESTS_DIR / (runfile_name[:-1] + ".log")
        if alternate_log.exists():
            alternate_log.unlink()
    return model


def test_pythonic_runs_match_procline():
    """All runfile ports should match every reported parameter."""

    for runfile_name in (
        "sampl1",
        "sampl2",
        "sampl2a",
        "sampl3",
        "sampl3a",
        "sampl4",
        "sampl5",
    ):
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
