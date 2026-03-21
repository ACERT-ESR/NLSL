"""Python equivalent of the ``sampl3a`` MOMD example."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import nlsl
from nlsl.data import process_spectrum

NSPLINE_POINTS = 400
BASELINE_EDGE_POINTS = 0

INITIAL_PARAMETERS = {
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

FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 400,
    "ftol": 1.0e-3,
    "xtol": 1.0e-3,
}

SETUP_COMMANDS = [
    "basis sampl3",
    "search rbar",
    "search c20",
    "axial r",
]

FIRST_VARY = ["vary rpll,rprp,c20,gib0"]
SECOND_PHASE = [
    "fix rpll",
    "fix rprp",
    "spherical r",
    "vary rbar",
]
FINAL_COMMANDS = ["vary gib2"]


def main():
    """Drive the ``sampl3a`` fit and display the converged spectra."""

    examples_dir = Path(__file__).resolve().parent
    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)
    model.shift = True

    for command in SETUP_COMMANDS:
        model.procline(command)

    # ``normalize=True`` preprocesses the experimental trace onto the
    # same normalized scale as the old loader, but the explicit
    # processed-data workflow does not preserve the legacy ``nrmlz``
    # bookkeeping flag.
    processed = process_spectrum(
        examples_dir / "sampl3.dat",
        NSPLINE_POINTS,
        BASELINE_EDGE_POINTS,
        derivative_mode=model.derivative_mode,
        normalize=True,
    )
    stop = float(processed.start) + float(processed.step) * max(
        int(processed.y.size) - 1,
        0,
    )
    # TODO ☐: explanatory comment needed for generate_coordinates (carry
    #         through to all examples).
    # TODO ☐: it should not be needed to convert to float here --
    #         generate_coordinates should do that as needed. (This should
    #         affect all examples)
    idx = model.generate_coordinates(
        float(processed.start),
        stop,
        int(processed.y.size),
    )
    model.data = processed.y
    model.name(str(examples_dir / "sampl3"), spectrum=idx)
    model.noise(processed.noise, spectrum=idx)

    for key in FIT_CONTROLS:
        model.fit_params[key] = FIT_CONTROLS[key]

    for command in FIRST_VARY:
        model.procline(command)
    model.fit()

    for command in SECOND_PHASE:
        model.procline(command)
    model.fit()

    for command in FINAL_COMMANDS:
        model.procline(command)
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
    experimental_series = []
    simulated_series = []
    component_series = []
    for idx in range(int(model.nspec)):
        experimental_series.append(
            experimental_block[idx, windows[idx]]
        )
        simulated_series.append(
            simulated_total[idx, windows[idx]]
        )
        component_series.append(
            component_curves[idx, :, windows[idx]]
        )

    combined_num = 0.0
    combined_den = 0.0
    for idx, experimental in enumerate(experimental_series):
        simulated = simulated_series[idx]
        residual = simulated - experimental
        numerator = float(np.linalg.norm(residual))
        denominator = float(np.linalg.norm(experimental))
        if denominator > 0.0:
            print(
                f"sampl3a spectrum {idx + 1}: relative rms ="
                f" {numerator / denominator:.6f}"
            )
        combined_num += numerator**2
        combined_den += denominator**2
    if combined_den > 0.0:
        print(
            "sampl3a: overall relative rms ="
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
    axes[0].set_title("sampl3a MOMD fit reproduced from Python")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
