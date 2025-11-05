"""Python-native fit mirroring the ``c16pc371e`` runfile."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import nlsl

NSPLINE_POINTS = 400
DERIVATIVE_MODE = 1

INITIAL_PARAMETERS = {
    "in2": 2,
    "gxx": 2.0084,
    "gyy": 2.0060,
    "gzz": 2.0020,
    "axx": 5.2,
    "ayy": 5.2,
    "azz": 33.2,
    "b0": 3320.0,
    "lemx": 8,
    "lomx": 7,
    "kmx": 4,
    "mmx": 4,
    "ipnmx": 2,
}

SITE_PARAMETER_VALUES = {
    "rpll": np.array([9.0, 9.5]),
    "rprp": np.array([8.21131, 9.0]),
    "c20": np.array([0.86071, 2.04448]),
    "c22": np.array([-0.67687, 1.00135]),
    "oss": np.array([0.0, 8.938]),
    "gib0": np.array([1.5, 2.0]),
    "gib2": np.array([-0.5, -0.5]),
}

GLOBAL_CONTROLS = {
    "nstep": 300,
    "cgtol": 1.0e-4,
    "shiftr": 1.0,
}

INITIAL_VARIATIONS = {
    "gib0": [1, 2],
    "rbar": [1, 2],
}

INITIAL_FIT = {
    "maxitr": 10,
    "maxfun": 250,
}

REFINEMENT_STEPS = (
    {
        "remove": (),
        "add": {
            "rprp": [1],
            "rpll": [1],
            "c20": [1],
            "c22": [1],
            "oss": [2],
        },
    },
    {
        "remove": ("rprp", "rpll", "c20", "c22", "oss"),
        "add": {
            "rprp": [2],
            "rpll": [2],
            "c20": [2],
            "c22": [2],
        },
    },
    {
        "remove": ("rprp", "rpll", "c20", "c22"),
        "add": {
            "gib0": [1, 2],
        },
    },
)


def apply_variation_changes(model, step):
    """Synchronise the vary list with the supplied stage configuration."""

    for token in step["remove"]:
        if token in model.fit_params.vary:
            del model.fit_params.vary[token]
    for token, indices in step["add"].items():
        model.fit_params.vary[token] = {"index": indices}


def main():
    """Fit the COS-7 spectrum without invoking the legacy runfile."""

    examples_dir = Path(__file__).resolve().parent
    model = nlsl.nlsl()
    model["nsites"] = 2
    model.update(INITIAL_PARAMETERS)
    model.update(SITE_PARAMETER_VALUES)
    model.update(GLOBAL_CONTROLS)

    for token, indices in INITIAL_VARIATIONS.items():
        model.fit_params.vary[token] = {"index": indices}

    model.load_data(
        examples_dir / "c16pc371e.dat",
        nspline=NSPLINE_POINTS,
        bc_points=0,
        shift=True,
        normalize=True,
        derivative_mode=DERIVATIVE_MODE,
    )

    for key, value in INITIAL_FIT.items():
        model.fit_params[key] = value

    site_spectra = model.fit()

    for step in REFINEMENT_STEPS:
        apply_variation_changes(model, step)
        site_spectra = model.fit()

    site_spectra = model.fit()

    model.write_spc()

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
    windows = model.layout["relative_windows"]
    experimental_series = tuple(
        experimental_block[idx, window] for idx, window in enumerate(windows)
    )
    simulated_series = tuple(
        simulated_total[idx, window] for idx, window in enumerate(windows)
    )
    component_series = tuple(
        component_curves[idx, :, window] for idx, window in enumerate(windows)
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
                f"c16pc371e spectrum {idx + 1}: relative rms ="
                f" {numerator / denominator:.6f}"
            )
        combined_num += numerator**2
        combined_den += denominator**2
    if combined_den > 0.0:
        print(
            "c16pc371e: overall relative rms ="
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
    axes[0].set_title("c16pc371e experimental fit reproduced from Python")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
