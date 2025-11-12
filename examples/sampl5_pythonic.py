"""Python conversion of the multi-spectrum ``sampl5`` workflow."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import nlsl

NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

INITIAL_PARAMETERS = {
    "nsite": 2,
    "in2": 2,
    "gxx": 2.0089,
    "gyy": 2.0021,
    "gzz": 2.0058,
    "axx": 5.6,
    "ayy": 33.8,
    "azz": 5.3,
    "b0": 3405.0,
    "rbar": 7.5,
    "n": 1.0,
    "gib0": 3.0,
    "gib2": 0.0,
    "betad": 15.0,
}

FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
}

def main():
    """Execute the ``sampl5`` MOMD refinement from Python."""

    examples_dir = Path(__file__).resolve().parent
    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    model.series("psi", (0.0, 90.0))
    model.load_basis("sampl5")
    model["c20"] = np.array([4.0, 0.2])
    model["c22"] = np.array([1.0, 0.0])
    model["nort"] = np.array([0.0, 10.0])

    model.load_data(
        examples_dir / "sampl500.dat",
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=True,
        derivative_mode=DERIVATIVE_MODE,
    )
    model.load_data(
        examples_dir / "sampl590.dat",
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=True,
        derivative_mode=DERIVATIVE_MODE,
    )

    model.weights = np.ones((2, 2))

    for key in FIT_CONTROLS:
        model.fit_params[key] = FIT_CONTROLS[key]

    all_sites = {"index": [1, 2]}
    model.fit_params.vary["gib0"] = all_sites
    model.fit_params.vary["rbar"] = all_sites
    model.fit()

    model.fit_params.vary["c20"] = all_sites
    model.fit()

    model.fit_params.vary["n"] = all_sites
    model.fit_params.vary["c22"] = all_sites
    model.fit()

    model.search("gib2", site=1)
    model.search("gib2", site=2)
    model.search("rbar", site=2)

    model.fit_params.vary.clear()
    for token in ("rbar", "n", "c20", "c22", "gib0"):
        model.fit_params.vary[token] = {"index": 2}
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
    fields = []
    experimental_series = []
    simulated_series = []
    component_series = []
    for idx in range(int(model.layout["nspc"])):
        fields.append(
            float(model.layout["sbi"][idx])
            + float(model.layout["sdb"][idx])
            * np.arange(int(model.layout["npts"][idx]))
        )
        experimental_series.append(
            experimental_block[idx, model.layout["relative_windows"][idx]]
        )
        simulated_series.append(
            simulated_total[idx, model.layout["relative_windows"][idx]]
        )
        component_series.append(
            component_curves[idx, :, model.layout["relative_windows"][idx]]
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
                f"sampl5 spectrum {idx + 1}: relative rms ="
                f" {numerator / denominator:.6f}"
            )
        combined_num += numerator**2
        combined_den += denominator**2
    if combined_den > 0.0:
        print(
            "sampl5: overall relative rms ="
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
    axes[0].set_title(
        "sampl5 two-orientation, two-site fit reproduced from Python"
    )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
