"""Reproduce the ``sampl4`` runfile entirely from Python."""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import nlsl

NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

# These entries mirror the manual ``let`` statements in ``sampl4.run`` and
# supply the same initial guesses for the fit.
INITIAL_PARAMETERS = {
    "nsite": 2,
    "in2": 2,
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "lemx": 12,
    "lomx": 10,
    "kmx": 7,
    "mmx": 7,
    "ipnmx": 2,
    "gib0": 0.5,
    "rx": np.array([np.log10(3.0e8), np.log10(1.0e7)]),
}

# Optimisation tolerances are applied through ``fit_params`` instead of the
# ``fit maxit ...`` directive from the runfile.
FIT_CONTROLS = {
    "maxitr": 40,
    "maxfun": 1000,
    "ftol": 1.0e-2,
    "gtol": 1.0e-6,
    "xtol": 1.0e-4,
}

# These parameters are refined during the optimisation.  The new
# ``fit_params.vary`` mapping mirrors the Fortran vary list so each entry
# below behaves the same way as the original ``vary`` commands in the runfile.
PARAMETERS_TO_VARY = {
    "gib0": [1, 2],
    "rbar": [1, 2],
}


def main():
    """Run the ``sampl4`` optimisation and plot the resulting spectra."""

    examples_dir = Path(__file__).resolve().parent
    data_path = examples_dir / "sampl4.dat"

    model = nlsl.nlsl()
    model.update(INITIAL_PARAMETERS)

    # ``data ... nspline 200 bc 20 shift`` from the runfile, executed through
    # the modern Python entry point so the processed intensities stay in
    # memory.
    model.load_data(
        data_path,
        nspline=NSPLINE_POINTS,
        bc_points=BASELINE_EDGE_POINTS,
        shift=True,
        normalize=False,
        derivative_mode=DERIVATIVE_MODE,
    )

    for token, indices in PARAMETERS_TO_VARY.items():
        model.fit_params.vary[token] = {"index": indices}

    for key, value in FIT_CONTROLS.items():
        model.fit_params[key] = value

    # Seed equal site populations once the data slot exists so the weight
    # matrix is initialised through the public mapping interface.
    model.weights = np.ones(2)

    # The original script issues two ``fit`` commands; calling ``fit()`` twice
    # reproduces the same refinement cycle and leaves the captured spectra
    # ready for plotting.
    model.fit()
    site_spectra = model.fit()

    fields = model.field_axes[0]

    experimental_block = model.experimental_data
    experimental = np.squeeze(experimental_block)

    weights = model.weights
    # ``weights`` behaves like a 1D vector when a single spectrum is active,
    # but the optimiser keeps a full ``(nspc, nsite)`` matrix in the
    # multi-spectrum case.  ``site_spectra`` always carries one row per site
    # with the sampled points arranged along the second axis.
    if weights.ndim == 1:
        weighted_components = weights[:, np.newaxis] * site_spectra
        simulated_total = weights @ site_spectra
        component_curves = weighted_components
    else:
        # For multiple spectra ``weights`` supplies one row per recorded trace.
        # Broadcasting the additional axis keeps each site's contribution tied
        # to the spectrum that owns it before we flatten the view for plotting.
        weighted_components = (
            weights[:, :, np.newaxis] * site_spectra[np.newaxis, :, :]
        )
        simulated_total = weights @ site_spectra
        component_curves = weighted_components.reshape(
            -1, site_spectra.shape[1]
        )

    residual = simulated_total - experimental
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(experimental)
    print(f"sampl4: relative rms = {rel_rms:.6f}")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(
        fields,
        experimental,
        color="black",
        linewidth=1.0,
        label="experimental",
    )
    ax.plot(
        fields,
        simulated_total,
        color="#d62728",
        linewidth=2.0,
        alpha=0.8,
        label="simulated sum",
    )

    colours = ["#1f77b4", "#2ca02c"]
    for idx, component in enumerate(component_curves):
        ax.plot(
            fields,
            component,
            color=colours[idx % len(colours)],
            linewidth=1.2,
            alpha=0.7,
            label=f"component {idx + 1}",
        )

    ax.set_xlabel("Magnetic field (G)")
    ax.set_ylabel("Intensity (arb. units)")
    ax.set_title("sampl4 two-component fit reproduced from Python")
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
