"""Plot the sampl4 multi-component model without running a fit.

This script mirrors the two-component ``sampl4.run`` example but skips the
least-squares optimisation entirely.  The optimal parameters reported in the
reference log (``sampl4.log_ref``) are copied into ``OPTIMAL_PARAMETERS`` and
applied directly, so the Fortran core only has to evaluate a single spectrum.

Experimental data are loaded purely for plotting via :mod:`nlsl.data` helper
functions; the measurements are *not* registered with the legacy Fortran
extension.  Instead we call :meth:`nlsl.nlsl.generate_coordinates` so that
:meth:`nlsl.nlsl.current_spectrum` produces the calculated components on the
same field axis.  This highlights how to generate the summed spectrum and the
individual site contributions using a fixed, known-good parameter set.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np

import nlsl
from nlsl.data import process_spectrum

# Final parameters taken from the converged ``sampl4`` fit reported in
# ``examples/sampl4.log_ref``.  The rotational diffusion rates are expressed in
# the logarithmic units used by the runfile (base-10), and the site weights map
# directly to the 71.47% / 28.53% populations listed in the log output.
OPTIMAL_PARAMETERS: Dict[str, object] = {
    "nsite": 2,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "in2": 2,
    "b0": 3400.0,
    "gib0": 1.994103,
    "rbar": (7.141443, 7.838921),
    "lemx": 12,
    "lomx": 10,
    "kmx": 7,
    "mmx": (7, 7),
    "ipnmx": 2,
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

    site_spectra, _ = model.current_spectrum()

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
