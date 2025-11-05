import os
import math
import numpy as np
import nlsl


def run_sample1_manual():
    intro = (
        "*" * 70
        + """file SAMPL1.RUN:  sample NLSL script file

  Illustrates fitting of anisotropic rotation of CSL spin probe
  in an isotropic solvent at X-band.

  Test data in file SAMPL1.DAT calculated with the following parameters:
               {g}   = 2.0089, 2.0021, 2.0058
               {A}   = 5.6, 33.8, 5.3  (gauss)
               betad = 15 degrees
               Rpll  = 1e7
               Rperp = 1e8
               B0    = 3400 G
               GIB   = 2.0 G (p-p width of Gaussian inhomog. linewidth)
**********************************************************************

  --- Open file 'sampl1.log' to save a record of this session
"""
    )
    print(intro)

    examples_dir = os.path.join(
        os.path.dirname(__file__), os.pardir, "examples"
    )
    os.chdir(examples_dir)
    n = nlsl.nlsl()
    data_files_out = []

    def procline(cmd):
        if cmd.startswith("data "):
            n.procline(cmd)
            data_files_out.append(cmd[5:].strip().split(" ")[0])
        else:
            n.procline(cmd)

    procline("log sampl1")
    print()
    print("  --- Set magnetic parameters for CSL spin probe")
    print()
    n.update(
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
    print()
    print("  --- Specify spectrometer field and make initial estimates for")
    print('  --- fitting parameters using the "let" statement as shown.')
    print("  --- Note in particular that the rotational rate constants")
    print("  --- are fit in log space, so that the parameters RPLL and RPRP")
    print("  --- are log10 of the rate constants for rotation around the")
    print("  --- axes parallel and perpendicular to the long axis of")
    print("  --- the molecule, respectively.")
    print("  ---")
    print(
        "  --- Note also that the log function may be used in a let statement."
    )
    print("  ---")
    print("  --- GIB0 is the Gaussian inhomogeneous broadening.")
    n.update(
        {
            "rpll": math.log10(1.0e8),
            "rprp": 8.0,
            "gib0": 1.5,
            "lemx": 6,
            "lomx": 5,
            "kmx": 4,
            "mmx": (2, 2),
        }
    )
    print()
    print("  --- Specify basis set truncation parameters")
    print()
    print("   --- Read in ASCII datafile 'sampl1.dat':")
    print("   ---    (1) Spline interpolate the data to 200 points")
    print(
        "   ---    (2) baseline-correct by fitting a line to 20 points at"
        " each end"
    )
    print("   ---    (3) allow shifting of B0 to maximize overlap with data")
    procline("data sampl1 ascii nspline 200 bc 20 shift")
    print()
    print("   --- Specify parameters to be varied in fitting procedure")
    print()
    for token in ("rpll", "rprp", "gib0"):
        n.fit_params.vary[token] = True
    print()
    print("   --- Carry out nonlinear least-squares procedure:")
    print("   ---    (1) Stop after a maximum of 40 iterations")
    print("   ---    (2) Stop after a maximum of 600 spectral calculations")
    print("   ---    (3) Chi-squared convergence tolerance is 1 part in 10^3")
    n.fit_params["maxitr"] = 40
    n.fit_params["maxfun"] = 1000
    n.fit_params["ftol"] = 1e-3
    n.fit_params["xtol"] = 1e-3
    n.fit()
    n.write_spc()
    procline("log end")

    rel_rms_list = []
    for fname in data_files_out:
        data = np.loadtxt(fname + ".spc")
        exp_sq = np.sum(data[:, 1] ** 2)
        rms_sq = np.sum((data[:, 2] - data[:, 1]) ** 2)
        if exp_sq > 0:
            rel_rms_list.append(math.sqrt(rms_sq) / math.sqrt(exp_sq))

    final_params = {}
    import re

    with open("sampl1.log") as fp:
        lines = fp.readlines()
    start = None
    for i, line in enumerate(lines):
        if "Final Parameters" in line:
            start = i + 4  # skip header lines
            break
    if start is not None:
        for line in lines[start:]:
            if not line.strip() or line.lstrip().startswith("Confidence"):
                break
            m = re.match(r"\s*(\w+)\s*=\s*([-0-9.Ee+]+)", line)
            if m:
                final_params[m.group(1).lower()] = float(m.group(2))

    return rel_rms_list, final_params


def test_sample1_only_manual():
    rel_rms, params = run_sample1_manual()

    assert rel_rms and all(r < 0.0404 * 1.02 for r in rel_rms)

    expected = {
        "gib0": 2.088995,
        "rprp": 6.982988,
        "rpll": 8.017160,
        "spctrm1": 3.115812,
    }

    for key, val in expected.items():
        assert key in params
        assert math.isclose(params[key], val, rel_tol=2e-3)
