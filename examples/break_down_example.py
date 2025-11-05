import numpy as np
import nlsl

# create a fresh NLSL instance
n = nlsl.nlsl()
print("off")
print("**********************************************************************")
print("file SAMPL1.RUN:  sample NLSL script file\n")
print("  Illustrates fitting of anisotropic rotation of CSL spin probe ")
print("  in an isotropic solvent at X-band.\n")
print(
    "  Test data in file SAMPL1.DAT calculated with the following parameters:"
)
print("              {g}   = 2.0089, 2.0021, 2.0058")
print("              {A}   = 5.6, 33.8, 5.3  (gauss)")
print("              betad = 15 degrees")
print("              Rpll  = 1e7")
print("              Rperp = 1e8")
print("              B0    = 3400 G")
print("              GIB   = 2.0 G (p-p width of Gaussian inhomog. linewidth)")
print("*" * 70)
print('  --- Open file "sampl1.log" to save a record of this session\n')
n.procline("log sampl1")
print()
print("  --- Set magnetic parameters for CSL spin probe \n")
print("********************************************")
print("Magnetic parameters for CSL spin probe")
print("********************************************")
n.update({
    "gxx": 2.0089,
    "gyy": 2.0021,
    "gzz": 2.0058,
    "in2": 2,
    "axx": 5.6,
    "ayy": 33.8,
    "azz": 5.3,
    "betad": 15,
})
print("********************************************")
print("CSL spin probe parameters loaded")
print("********************************************\n")
print("  --- Specify spectrometer field and make initial estimates for")
print('  --- fitting parameters using the "let" statement as shown.')
print("  --- Note in particular that the rotational rate constants")
print("  --- are fit in log space, so that the parameters RPLL and RPRP")
print("  --- are log10 of the rate constants for rotation around the ")
print("  --- axes parallel and perpendicular to the long axis of")
print("  --- the molecule, respectively. ")
print("  ---")
print("  --- Note also that the log function may be used in a let statement.")
print("  ---")
print("  --- GIB0 is the Gaussian inhomogeneous broadening.")
n.update({
    "rpll": np.log(1.0e8),
    "rprp": 8.0,
    "gib0": 1.5,
})
print()
print("  --- Specify basis set truncation parameters\n")
n.update({"lemx": 6, "lomx": 5, "kmx": 4, "mmx": (2, 2)})
print()
print('   --- Read in ASCII datafile "sampl1.dat":')
print("   ---    (1) Spline interpolate the data to 200 points")
print(
    "   ---    (2) baseline-correct by fitting a line to 20 points at each end"
)
print("   ---    (3) allow shifting of B0 to maximize overlap with data\n")
n.procline("data sampl1 ascii nspline 200 bc 20 shift")
print()
print("   --- Specify parameters to be varied in fitting procedure\n")
for token in ("rpll", "rprp", "gib0"):
    n.fit_params.vary[token] = True
print()
print("   --- Carry out nonlinear least-squares procedure:")
print("   ---    (1) Stop after a maximum of 40 iterations")
print("   ---    (2) Stop after a maximum of 600 spectral calculations")
print("   ---    (3) Chi-squared convergence tolerance is 1 part in 10^3\n")
n.fit_params["maxitr"] = 40
n.fit_params["maxfun"] = 1000
n.fit_params["ftol"] = 1e-3
n.fit_params["xtol"] = 1e-3
n.fit()
print(dict(n.items()))
