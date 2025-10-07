"""Plot the sampl4 multi-component model without running a fit.

This script mirrors the two-component ``sampl4.run`` example but skips the
least-squares optimisation entirely.  The optimal parameters reported for run 4
are copied into ``OPTIMAL_PARAMETERS`` and applied directly, so the Fortran core
only has to evaluate a single spectrum.

Experimental data are loaded purely for plotting via :mod:`nlsl.data` helper
functions; the measurements are *not* registered with the legacy Fortran
extension.  Instead we call :meth:`nlsl.nlsl.generate_coordinates` so that
:meth:`nlsl.nlsl.current_spectrum` produces the calculated components on the
same field axis.  This highlights how to generate the summed spectrum and the
individual site contributions using a fixed, known-good parameter set.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import nlsl
from nlsl.data import process_spectrum

# Final parameters taken from the converged ``sampl4`` fit reported in the run
# 4 log.  Values are copied verbatim from the ``items()`` dump so the simulated
# trace matches the published solution.
OPTIMAL_PARAMETERS = {
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

SAMPL4_SITE_WEIGHTS = np.array([0.7147334, 0.2852620])
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

def main() -> None:
    """Render the sampl4 components using the hard-coded optimal parameters."""

    examples_dir = Path(__file__).resolve().parent
    data_path = examples_dir / "sampl4.dat"

    processed = process_spectrum(
        data_path,
        NSPLINE_POINTS,
        BASELINE_EDGE_POINTS,
        derivative_mode=DERIVATIVE_MODE,
        normalize=True,
    )

    fields = processed.x
    experimental = processed.y

    model = nlsl.nlsl()
    model.update(OPTIMAL_PARAMETERS)

    model.generate_coordinates(
        fields.size,
        start=processed.start,
        step=processed.step,
        derivative_mode=DERIVATIVE_MODE,
        baseline_points=BASELINE_EDGE_POINTS,
        normalize=True,
        nspline=NSPLINE_POINTS,
        label="sampl4 synthetic",
        reset=True,
    )

    site_spectra, _ = model.current_spectrum

    if site_spectra.shape[0] != SAMPL4_SITE_WEIGHTS.size:
        raise RuntimeError("Unexpected number of site spectra returned")

    weighted_components = SAMPL4_SITE_WEIGHTS[:, np.newaxis] * site_spectra
    simulated_total = np.sum(weighted_components, axis=0)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(fields, experimental, color="black", linewidth=1.0, label="experimental (normalised)")
    ax.plot(fields, simulated_total, color="#d62728", linewidth=2.0, alpha=0.8, label="simulated sum")

    colours = ["#1f77b4", "#2ca02c"]
    for idx, component in enumerate(weighted_components):
        ax.plot(
            fields,
            component,
            color=colours[idx % len(colours)],
            linewidth=1.2,
            alpha=0.7,
            label=f"component {idx + 1}",
        )

    ax.set_xlabel("Magnetic field (G)")
    ax.set_ylabel("Normalised intensity")
    ax.set_title("sampl4 two-component simulation using fixed optimal parameters")
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
