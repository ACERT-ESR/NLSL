"""Evaluate the Sampl4 multi-component model at a single parameter point.

The script mirrors the published Sampl4 fit by loading the raw experimental
data from ``sampl4.dat``, applying the same spline resampling and baseline
correction requested by the original run file, and then calling
``nlsl.nlsl.current_spectrum`` to evaluate the two-component model without
invoking the legacy ``procline`` interpreter.  The simulator is configured
purely from the known fit parameters, so the experimental spectrum is only
used for plotting and never communicated back to :mod:`nlsl`.
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

import nlsl

DATA_FILE = Path(__file__).with_name("sampl4.dat")

# Magnetic, dynamic and basis parameters reported at the end of sampl4.log_ref.
FIT_PARAMETERS = {
    "gxx": 2.0089,
    "gyy": 2.0063,
    "gzz": 2.0021,
    "in2": 2,
    "axx": 5.0,
    "ayy": 5.0,
    "azz": 33.0,
    "lemx": 12,
    "lomx": 10,
    "kmx": 7,
    "mmx": (7, 7),
    "ipnmx": 2,
    "b0": 3400.5,
    "gib0": 1.994103,
    "rbar": [7.141443, 7.838921],
}

# Least-squares scale factors (site populations) from the published fit.
SITE_POPULATIONS = np.array([0.7147334, 0.2852620])

TARGET_POINTS = 200
BASELINE_POINTS = 20
FIELD_RANGE = 100.0
DERIVATIVE_ORDER = 1


def load_raw_spectrum() -> tuple[np.ndarray, np.ndarray]:
    raw = np.loadtxt(DATA_FILE)
    return raw[:, 0], raw[:, 1]


def resample_and_correct(
    target_fields: np.ndarray, raw_fields: np.ndarray, raw_signal: np.ndarray
) -> np.ndarray:
    """Spline resample the raw spectrum onto ``target_fields`` and remove a linear baseline."""

    resampled_signal = np.interp(target_fields, raw_fields, raw_signal)
    total_points = target_fields.size
    edge_idx = np.r_[0:BASELINE_POINTS, total_points - BASELINE_POINTS : total_points]
    baseline_coeffs = np.polyfit(target_fields[edge_idx], resampled_signal[edge_idx], 1)
    baseline = np.polyval(baseline_coeffs, target_fields)
    return resampled_signal - baseline


def main() -> None:
    raw_fields, raw_signal = load_raw_spectrum()

    sim = nlsl.nlsl()
    sim.nsites = SITE_POPULATIONS.size
    sim.update(FIT_PARAMETERS)
    sim["nfield"] = TARGET_POINTS
    sim["range"] = FIELD_RANGE
    sim["ideriv"] = DERIVATIVE_ORDER

    site_spectra, _ = sim.current_spectrum
    layout = sim.layout
    segment_points = int(layout["npts"][0])
    offsets = np.arange(segment_points, dtype=float)
    model_fields = layout["sbi"][0] + offsets * layout["sdb"][0]

    experimental = resample_and_correct(model_fields, raw_fields, raw_signal)

    weighted_components = site_spectra * SITE_POPULATIONS[:, None]
    combined_model = weighted_components.sum(axis=0)

    scale = (combined_model @ experimental) / (combined_model @ combined_model)
    scaled_components = weighted_components * scale
    scaled_model = combined_model * scale

    rel_rms = np.linalg.norm(scaled_model - experimental) / np.linalg.norm(experimental)
    print(f"Relative RMS deviation after linear rescaling: {rel_rms:.3f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(model_fields, experimental, label="Experimental", color="black", linewidth=1.0)
    ax.plot(model_fields, scaled_model, label="Model (sum)", color="tab:blue")
    for idx, component in enumerate(scaled_components, start=1):
        ax.plot(model_fields, component, linestyle="--", label=f"Component {idx}")

    ax.set_xlabel("Magnetic field / G")
    ax.set_ylabel("Signal (arb. units)")
    ax.set_title("Sampl4 multi-component single-point evaluation")
    ax.legend()
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
