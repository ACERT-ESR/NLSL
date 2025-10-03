"""Pure-Python helpers for ingesting experimental spectra.

The routines in this module mirror the algorithms used by the legacy
``datac`` Fortran command.  They provide ASCII spectrum loading,
baseline correction, natural cubic spline resampling and optional
normalization so that data prepared in Python is bitwise compatible with
the original code path.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np

_COMMENT_PREFIXES = ("c", "C", "!")


@dataclass
class ProcessedSpectrum:
    """Container describing a processed experimental spectrum."""

    x: np.ndarray
    y: np.ndarray
    start: float
    step: float
    baseline_intercept: float
    baseline_slope: float
    noise: float
    normalization: float
    derivative_mode: int


def read_ascii_spectrum(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Read a two-column ASCII spectrum from ``path``."""

    x_values: list[float] = []
    y_values: list[float] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line[0] in _COMMENT_PREFIXES:
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed data line: {raw_line!r}")
            try:
                x_val = float(parts[0])
                y_val = float(parts[1])
            except ValueError as exc:  # pragma: no cover - defensive
                raise ValueError(f"Could not parse data line: {raw_line!r}") from exc
            x_values.append(x_val)
            y_values.append(y_val)

    if not x_values:
        raise ValueError(f"Spectrum {path} contained no data points")

    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def fit_linear_baseline(
    x: np.ndarray, y: np.ndarray, edge_points: int
) -> Tuple[np.ndarray, float, float, float]:
    """Fit and optionally subtract a linear baseline from ``y``."""

    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")

    n = x.size
    if n < 2:
        raise ValueError("At least two points are required for baseline fitting")

    nend = min(max(int(edge_points), 0), n // 3)
    if nend <= 0:
        nend = max(10, n // 20)
        nend = max(1, min(nend, n // 2))

    sx = sy = sx2 = sxy = 0.0
    for idx in range(nend):
        k = n - idx - 1
        sx += x[idx] + x[k]
        sy += y[idx] + y[k]
        sx2 += x[idx] * x[idx] + x[k] * x[k]
        sxy += x[idx] * y[idx] + x[k] * y[k]

    sn = float(2 * nend)
    denom = sn * sx2 - sx * sx
    if denom == 0:
        raise ZeroDivisionError("Degenerate baseline fit")

    intercept = (sx2 * sy - sx * sxy) / denom
    slope = (sn * sxy - sx * sy) / denom

    residual = 0.0
    for idx in range(nend):
        k = n - idx - 1
        residual += (y[idx] + y[k] - 2.0 * intercept - slope * (x[idx] + x[k])) ** 2
    noise = np.sqrt(residual / float(2 * nend - 1)) if nend > 0 else 0.0

    corrected = y - (intercept + slope * x)
    if edge_points <= 0:
        corrected = y.copy()

    return corrected, intercept, slope, noise


def _natural_cubic_spline_second_derivatives(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = x.size
    if n < 3:
        raise ValueError("At least three points are required for spline interpolation")

    y2 = np.zeros_like(y)
    u = np.zeros_like(y)

    for i in range(1, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (
            6.0
            * (
                (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                - (y[i] - y[i - 1]) / (x[i] - x[i - 1])
            )
            / (x[i + 1] - x[i - 1])
            - sig * u[i - 1]
        ) / p

    y2[-1] = 0.0
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]

    return y2


def natural_cubic_spline_resample(
    x: np.ndarray, y: np.ndarray, n_points: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Resample ``y`` at ``n_points`` using a natural cubic spline."""

    if n_points <= 0:
        raise ValueError("Number of spline points must be positive")
    if n_points < 2:
        raise ValueError("At least two spline points are required")

    y2 = _natural_cubic_spline_second_derivatives(x, y)
    x0 = float(x[0])
    x1 = float(x[-1])
    new_x = np.linspace(x0, x1, n_points, dtype=float)
    new_y = np.empty_like(new_x)

    klo = 0
    khi = 1
    for idx, xv in enumerate(new_x):
        while xv > x[khi] and khi < x.size - 1:
            klo = khi
            khi += 1
        h = x[khi] - x[klo]
        if h == 0.0:  # pragma: no cover - defensive
            raise ZeroDivisionError("Duplicate X values in spline input")
        a = (x[khi] - xv) / h
        b = (xv - x[klo]) / h
        new_y[idx] = (
            a * y[klo]
            + b * y[khi]
            + ((a**3 - a) * y2[klo] + (b**3 - b) * y2[khi]) * (h * h) / 6.0
        )

    return new_x, new_y


def _sglint(values: np.ndarray) -> float:
    total = float(np.sum(values))
    return total - 0.5 * float(values[0] + values[-1])


def _dblint(values: np.ndarray) -> float:
    running = 0.0
    total = 0.0
    for val in values:
        running += val
        total += running
    return total - 0.5 * float(values[0] + running)


def normalize_spectrum(
    values: np.ndarray, spacing: float, derivative_mode: int
) -> Tuple[np.ndarray, float]:
    """Normalize ``values`` according to ``derivative_mode``."""

    norm_values = values.copy()
    if derivative_mode == 0:
        norm = spacing * _sglint(norm_values)
    else:
        integral = _sglint(norm_values)
        base = integral / float(norm_values.size)
        norm_values -= base
        norm = spacing * spacing * _dblint(norm_values)

    if norm != 0.0:
        norm_values /= norm

    return norm_values, norm


def process_spectrum(
    path: Path,
    nspline: int,
    bc_points: int,
    derivative_mode: int = 1,
    normalize: bool = True,
) -> ProcessedSpectrum:
    """Load and preprocess a spectrum from ``path``."""

    x_raw, y_raw = read_ascii_spectrum(path)
    corrected, intercept, slope, noise = fit_linear_baseline(x_raw, y_raw, bc_points)

    if nspline > 0:
        xs, ys = natural_cubic_spline_resample(x_raw, corrected, nspline)
    else:
        xs, ys = x_raw.copy(), corrected.copy()

    step = float(xs[1] - xs[0]) if xs.size > 1 else 0.0
    norm_factor = 0.0
    if normalize:
        ys, norm_factor = normalize_spectrum(ys, step, derivative_mode)
        if norm_factor != 0.0:
            noise /= norm_factor

    return ProcessedSpectrum(
        x=xs,
        y=ys,
        start=float(xs[0]),
        step=step,
        baseline_intercept=intercept,
        baseline_slope=slope,
        noise=noise,
        normalization=norm_factor,
        derivative_mode=derivative_mode,
    )


__all__ = [
    "ProcessedSpectrum",
    "fit_linear_baseline",
    "natural_cubic_spline_resample",
    "normalize_spectrum",
    "process_spectrum",
    "read_ascii_spectrum",
]
