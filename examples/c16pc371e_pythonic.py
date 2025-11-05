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

# TODO: This is not pythonic.  Set all these
# through a dictionary-like format.  Check in
# general that NONE of the examples called
# "pythonic" use the procline function AT ALL.
SETUP_COMMANDS = [
    "sites 2",
    "let rpll(1) = 9",
    "let rprp(1) = 8.21131",
    "let rpll(2) = 9.5",
    "let rprp(2) = 9",
    "let c20(1) = 0.86071",
    "let c22(1) = -0.67687",
    "let c20(2) = 2.04448",
    "let c22(2) = 1.00135",
    "let oss(2) = 8.938",
    "let gib0(1) = 1.5",
    "let gib2(1) = -0.5",
    "let gib0(2) = 2.0",
    "let gib2(2) = -0.5",
    "let nstep,cgtol,shiftr = 300,1e-4,1.0",
]

INITIAL_FIT = {
    "maxitr": 10,
    "maxfun": 250,
}

# TODO: These should all be controlled
# pythonically through the
# FitParameterVaryMapping attribute of the model
# class.
# TODO: in keeping with the comments inside
# FitParameterVaryMapping, these need to be
# controlled in an array format now.
REFINEMENT_SEQUENCE = [
    ["vary rprp(1) rpll(1) c20(1) c22(1) oss(2)"],
    [
        "fix rprp(1) rpll(1) c20(1) c22(1) oss(2)",
        "vary rprp(2) rpll(2) c20(2) c22(2)",
    ],
    ["fix rprp(2) rpll(2) c20(2) c22(2)", "vary gib0(1)", "vary gib0(2)"],
]


def main():
    """Fit the COS-7 spectrum without invoking the legacy runfile."""

    examples_dir = Path(__file__).resolve().parent
    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    for command in SETUP_COMMANDS:
        model.procline(command)

    model.load_data(
        examples_dir / "c16pc371e.dat",
        nspline=NSPLINE_POINTS,
        bc_points=0,
        shift=True,
        normalize=True,
        derivative_mode=DERIVATIVE_MODE,
    )

    for key in INITIAL_FIT:
        model.fit_params[key] = INITIAL_FIT[key]

    site_spectra = model.fit()

    for commands in REFINEMENT_SEQUENCE:
        for command in commands:
            model.procline(command)
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
        component_curves = weights[:, :, np.newaxis] * site_spectra[np.newaxis, :, :]

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
                f"c16pc371e spectrum {idx + 1}: relative rms = {numerator / denominator:.6f}"
            )
        combined_num += numerator ** 2
        combined_den += denominator ** 2
    if combined_den > 0.0:
        print(
            f"c16pc371e: overall relative rms = {np.sqrt(combined_num / combined_den):.6f}"
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
