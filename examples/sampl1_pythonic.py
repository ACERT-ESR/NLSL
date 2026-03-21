"""Reproduce the ``sampl1`` runfile through the Python interface."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import nlsl
from nlsl.data import process_spectrum

NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20

INITIAL_PARAMETERS = {
    "in2": 2,
    "gxx": 2.0089,
    "gyy": 2.0021,
    "gzz": 2.0058,
    "axx": 5.6,
    "ayy": 33.8,
    "azz": 5.3,
    "lemx": 6,
    "lomx": 5,
    "kmx": 4,
    "mmx": 2,
    "ipnmx": 2,
    "rpll": np.log10(1.0e8),
    "rprp": 8.0,
    "gib0": 1.5,
}

FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
    "ftol": 1.0e-3,
    "xtol": 1.0e-3,
}

PARAMETERS_TO_VARY = ("rpll", "rprp", "gib0")


def main():
    """Execute the ``sampl1`` optimisation and visualise the fit."""

    examples_dir = Path(__file__).resolve().parent
    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    model.shift = True
    # ``normalize=True`` asks ``process_spectrum`` to rescale the experimental
    # trace before it is copied into NLSL.  For first-derivative data this
    # still removes the constant offset and normalizes the double integral,
    # but the processed-data path does not preserve the old ``nrmlz``
    # bookkeeping flag.
    processed = process_spectrum(
        examples_dir / "sampl1.dat",
        NSPLINE_POINTS,
        BASELINE_EDGE_POINTS,
        derivative_mode=model.derivative_mode,
        normalize=True,
    )
    stop = float(processed.start) + float(processed.step) * max(
        int(processed.y.size) - 1,
        0,
    )
    idx = model.generate_coordinates(
        float(processed.start),
        stop,
        int(processed.y.size),
    )
    model.data = processed.y
    model.name(str(examples_dir / "sampl1"), spectrum=idx)
    model.noise(processed.noise, spectrum=idx)

    for token in PARAMETERS_TO_VARY:
        # ``fit_params.vary`` mirrors the Fortran vary list, so toggling each
        # entry exposes the same optimisation controls as the legacy runfile.
        model.fit_params.vary[token] = True

    for key in FIT_CONTROLS:
        model.fit_params[key] = FIT_CONTROLS[key]

    site_spectra = model.fit()
    weights = model.weights
    if weights.ndim == 1:
        simulated_total = weights @ site_spectra
        simulated_total = simulated_total[np.newaxis, :]
        component_curves = weights[:, np.newaxis] * site_spectra
        component_curves = component_curves[np.newaxis, :, :]
    else:
        simulated_total = weights @ site_spectra
        component_curves = (
            weights[:, :, np.newaxis] * site_spectra[np.newaxis, :, :]
        )

    experimental_block = model.experimental_data
    fields = model.field_axes
    windows = model.relative_windows
    experimental_series = tuple(
        experimental_block[idx, window] for idx, window in enumerate(windows)
    )
    simulated_series = tuple(
        simulated_total[idx, window] for idx, window in enumerate(windows)
    )
    component_series = tuple(
        component_curves[idx, :, window] for idx, window in enumerate(windows)
    )

    numerators = []
    denominators = []
    for idx, experimental in enumerate(experimental_series):
        residual = simulated_series[idx] - experimental
        numerator = float(np.linalg.norm(residual))
        denominator = float(np.linalg.norm(experimental))
        if denominator > 0.0:
            numerators.append(numerator)
            denominators.append(denominator)
            print(
                f"sampl1 spectrum {idx + 1}: relative rms ="
                f" {numerator / denominator:.6f}"
            )
    if denominators:
        combined_num = float(np.sum(np.square(numerators)))
        combined_den = float(np.sum(np.square(denominators)))
        if combined_den > 0.0:
            print(
                "sampl1: overall relative rms ="
                f" {np.sqrt(combined_num / combined_den):.6f}"
            )

    figure, axes = plt.subplots(len(fields), 1, figsize=(10, 5 * len(fields)))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    colours = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for idx in range(len(fields)):
        axis = axes[idx]
        axis.plot(
            fields[idx],
            experimental_series[idx],
            color="black",
            linewidth=1.0,
            label="experimental",
        )
        axis.plot(
            fields[idx],
            simulated_series[idx],
            color="#d62728",
            linewidth=2.0,
            alpha=0.8,
            label="simulated",
        )
        for comp_idx in range(component_series[idx].shape[0]):
            axis.plot(
                fields[idx],
                component_series[idx][comp_idx],
                color=colours[comp_idx % len(colours)],
                linewidth=1.2,
                alpha=0.7,
                label=f"component {comp_idx + 1}",
            )
        axis.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        axis.legend(loc="upper right")
        axis.set_ylabel("Intensity (arb. units)")
    axes[-1].set_xlabel("Magnetic field (G)")
    axes[0].set_title("sampl1 anisotropic rotation fit reproduced from Python")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
