# Fortran source layout

| Directory | Role | Representative files |
|-----------|------|----------------------|
| `cli/` | User-facing command interpreter, status reporters, and tokenizers that turn scripted input into actions. These routines dispatch commands such as `ASSIGN`, `FIT`, `SERIES`, and provide shared parsing utilities. | `nlsl.f90`, `assgnc.f90`, `fitc.f90`, `series.f90`, `strutl1.f90` |
| `least_squares/` | Nonlinear-Levenberg–Marquardt driver, trust-region controls, QR-based scaling/shift utilities, and statistical post-processing. Keeps MINPACK-derived code with the routines that prepare residuals and interpret fit quality. | `lmnls.f90`, `lmpar.f90`, `qrsolv.f90`, `sscale.f90`, `sshift.f90`, `stats.f90` |
| `matrix_model/` | Spectral model assembly: basis/pruning logic, Liouville matrix construction, parameter validation, and the LM function callback that marshals spectra. Houses the physics-heavy pieces shared by both fitting and standalone calculations. | `matrll.f90`, `pmatrl.f90`, `lbasix.f90`, `setmts.f90`, `lcheck.f90`, `momdls.f90`, `lfun.f90` |
| `krylov_solvers/` | Lanczos/CG machinery and related transforms that evaluate spectra from matrix data, including continued fractions, matrix–vector products, starting-vector builders, and integration kernels. | `eprls.f90`, `cscg.f90`, `scmvm.f90`, `cgltri.f90`, `ccrints.f90`, `stvect.f90`, `pstvec.f90`, `cfs.f90` |
| `io_preprocessing/` | Data ingestion, spectrum setup, and output routines that translate between files, experimental arrays, and solver-ready structures, plus FFT/correlation helpers used during shifting and convolution. | `getdat.f90`, `setnm.f90`, `setspc.f90`, `writec.f90`, `writr.f90`, `correl.f90`, `gconvl.f90` |
| `core/` | Shared state modules and constants that define global dimensions, parameter arrays, filenames, and error codes used across all subsystems. Keeping them together clarifies cross-cutting dependencies. | `nlsdim.f90`, `parcom.f90`, `expdat.f90`, `lmcom.f90`, `mspctr.f90`, `lpnam.f90`, `errmsg.f90`, `rnddbl.f90` |
| `numerics/` | Low-level linear-algebra, BLAS, and special-function helpers shared by least-squares and Krylov code (QR utilities, rotations, FFT kernels, associated Legendre/Wigner functions). | `qrfac.f90`, `qrutil.f90`, `daxpy.f90`, `dcopy.f90`, `enorm.f90`, `dchex.f90`, `ftfuns.f90`, `plgndr.f90`, `w3j.f90` |

The directory layout follows the execution pipeline from command parsing through spectral modelling and Krylov evaluation while isolating shared infrastructure and mathematical utilities.
