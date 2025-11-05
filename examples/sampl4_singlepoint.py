from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

import nlsl
from nlsl.data import process_spectrum

# --- Hard-coded SAMPL4 setup (no imports from tests or references) ---
NSPLINE_POINTS = 200
BASELINE_EDGE_POINTS = 20
DERIVATIVE_MODE = 1

# Final parameters (mirror of classic runfile-4 solution)
SAMPL4_FINAL_PARAMETERS = {
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
    "rx": np.array([7.8396974, 7.14177897]),
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
    "b0": 3400.50251256,
    "fldi": 3350.5046000757857,
    "dfld": None,  # set from processed data below
    "gamman": 0.0,
    "cgtol": 0.001,
    "shiftr": 0.001,
    "shifti": 0.0,
    "range": 100.0,
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
    "nfield": NSPLINE_POINTS,
    "ideriv": DERIVATIVE_MODE,
    "iwflg": 0,
    "igflg": 0,
    "iaflg": 0,
    "jkmn": 0,
    "jmmn": 0,
    "irflg": 2,
    "ndim": 156,
}

# Site populations and spectral metadata for reproducing the converged
# simulation
SAMPL4_FINAL_WEIGHTS = np.array([0.2848810, 0.7155313])
SAMPL4_FINAL_SB0 = None  # set from processed data below (uses b0)
SAMPL4_FINAL_SRNG = None  # set from processed data below (uses range)
SAMPL4_FINAL_ISHFT = np.array([1], dtype=np.int32)
SAMPL4_FINAL_SHFT = np.array([0.0])
SAMPL4_FINAL_NRMLZ = np.array([0], dtype=np.int32)


def main():
    # Load and process the experimental trace on the same 200-pt grid
    data_path = Path(__file__).resolve().parent / "sampl4.dat"
    proc = process_spectrum(
        data_path,
        NSPLINE_POINTS,
        BASELINE_EDGE_POINTS,
        derivative_mode=DERIVATIVE_MODE,
        normalize=False,
    )

    field_start = float(proc.start)
    field_step = float(proc.step)
    point_count = proc.y.size
    y_exp = proc.y.copy()

    # Fill in dfld and spectral metadata derived from processed data
    params = dict(SAMPL4_FINAL_PARAMETERS)
    params["dfld"] = field_step

    sb0 = np.array([params["b0"]])
    srng = np.array([params["range"]])

    # Build the model strictly from the converged (hard-coded) parameters
    model = nlsl.nlsl()
    model["nsite"] = 2

    # Generate the coordinate grid to match the processed data exactly
    index, sl = model.generate_coordinates(
        point_count,
        start=field_start,
        step=field_step,
        derivative_mode=DERIVATIVE_MODE,
        baseline_points=BASELINE_EDGE_POINTS,
        normalize=False,
        nspline=NSPLINE_POINTS,
        shift=True,
        label="sampl4-single-eval",
        reset=True,
    )

    # Set the experimental intensities on the same slice
    model.set_data(sl, y_exp[:point_count])

    # Mirror the runfile-4 solution through the dict interface (no fitting)
    model.update(params)
    model["sb0"] = sb0
    model["srng"] = srng
    model["ishft"] = SAMPL4_FINAL_ISHFT
    model["shft"] = SAMPL4_FINAL_SHFT
    model["nrmlz"] = SAMPL4_FINAL_NRMLZ
    model.weights = SAMPL4_FINAL_WEIGHTS

    # Pull the per-site spectra directly (no calls to fit)
    site_spectra = model.current_spectrum  # shape: (nsite, npts)

    # X-axis matching the processed field grid
    x = field_start + field_step * np.arange(point_count)

    # --- Plot using stylistic example conventions; components are WEIGHTED ---
    weights = model.weights
    site_spectra = model.current_spectrum  # (nsite, npts)

    if weights.ndim == 1:
        weighted_components = weights[:, np.newaxis] * site_spectra
        simulated_total = weights @ site_spectra
        component_curves = weighted_components
    else:
        weighted_components = (
            weights[:, :, np.newaxis] * site_spectra[np.newaxis, :, :]
        )
        simulated_total = weights @ site_spectra
        component_curves = weighted_components.reshape(
            -1, site_spectra.shape[1]
        )

    residual = simulated_total - y_exp
    rel_rms = np.linalg.norm(residual) / np.linalg.norm(y_exp)
    print(f"relative rms = {rel_rms:.6f}")

    fig, ax = plt.subplots(figsize=(10, 6))
    # experimental (black)
    ax.plot(x, y_exp, color="black", linewidth=1.0, label="experimental")
    # simulated total (red)
    ax.plot(
        x,
        simulated_total,
        color="#d62728",
        linewidth=2.0,
        alpha=0.8,
        label="simulated sum",
    )

    # components (blue, green, ...), weighted
    colours = ["#1f77b4", "#2ca02c"]
    for idx, component in enumerate(component_curves):
        ax.plot(
            x,
            component,
            color=colours[idx % len(colours)],
            linewidth=1.2,
            alpha=0.7,
            label=f"component {idx + 1}",
        )

    ax.set_xlabel("Magnetic field (G)")
    ax.set_ylabel("Intensity (arb. units)")
    ax.set_title("sampl4 components from model.current_spectrum (no fit)")
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
