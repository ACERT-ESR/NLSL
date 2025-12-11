import numpy as np
import pyspecdata as psd
import re
from pathlib import Path
from numpy import r_
from pyspecdata.datadir import pyspec_config
from scipy.optimize import differential_evolution
import nlsl

# Register the example directory with pyspecdata if it is not already present.
if not Path(psd.getDATADIR("nlsl_examples")).exists():
    pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(Path.cwd()))

# Load and prepare the experimental spectrum using the same workflow as the
# pyspecdata example.  All TODO markers from that script are intentionally left
# untouched here.
d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
d.set_units("$B_0$", None)
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]

n = nlsl.nlsl()

# ☐ TODO -- the max points should be an accessible property or attribute supplied by the
# class -- we should not need to do the following in our script
max_points = n._core.expdat.data.shape[0] // max(
    n._core.expdat.nft.shape[0], 1
)
# {{{ we use convolution to downsample the data
if d.data.shape[0] > max_points:
    divisor = d.shape["$B_0$"] // max_points + 1
    dB = np.diff(d["$B_0$"][r_[0, 1]]).item()
    d_orig_max = d.data.max()
    d.convolve("$B_0$", dB / 6 * divisor)
    d = d["$B_0$", 0::divisor]
    d *= d_orig_max / d.data.max()
# }}}

# Provide reasonable starting parameters
initial_params = {
    "in2": np.int32(2),
    "gxx": np.float64(2.0089),
    "gyy": np.float64(2.0058),
    "gzz": np.float64(2.0021),
    "axx": np.float64(4.9),
    "ayy": np.float64(4.9),
    "azz": np.float64(33.0),
    "lemx": np.int32(6),
    "lomx": np.int32(5),
    "kmx": np.int32(4),
    "mmx": np.int32(4),
    "ipnmx": np.int32(2),
    "rbar": np.float64(7.935044819885658),
    "n": np.float64(0.6932829362508512),
    "c20": np.float64(2.5061887785331742),
    "c22": np.float64(-1.2608085913642202),
    "betad": np.float64(29.942804557290664),
    "gib0": np.float64(2.0101124440329796),
    "gib2": np.float64(0.5050889997199121),
}
n.update(initial_params)
param_tokens = list(
    set(initial_params.keys())
    - set(["in2", "kmx", "mmx", "lemx", "lomx", "ipnmx", "b0"])
)
print(param_tokens)

d.data /= d.data.max() - d.data.min()

# Load nddata into the optimiser buffers
n.load_nddata(d)

# -------------------------
# Allow differential evolution to optimize ALL scalar parameters
# -------------------------

bounds = []
for k in param_tokens:
    v = n[k]
    if k in ["gxx", "gyy", "gzz"]:
        bounds.append((((v - 2) * 0.6) + 2, ((v - 2) * 1.4) + 2))
    elif k == "betad":
        bounds.append((0, 180))
    elif v < 0:
        bounds.append((v * 2.0, v * 0.1))
    else:
        bounds.append((v * 0.1, v * 2.0))


def objective(n):
    site_spectra = n.current_spectrum
    simulated_total = np.squeeze(n.weights @ site_spectra)
    simulated_total /= simulated_total.max() - simulated_total.min()
    return simulated_total


def residual_norm(candidate):
    for value, token in zip(candidate, param_tokens):
        n[token] = value
    return np.linalg.norm(d.data - objective(n))


# Target residual (L2 norm) for early stopping; tune as needed
target_residual = 5e-2


def de_callback(xk, convergence):
    """Progress callback for differential evolution.

    Prints iteration, current best residual, and convergence metric.
    Returns True if target_residual has been reached, causing early stop.
    """
    current_resid = residual_norm(xk)
    print(
        f"DE iter {de_callback.iteration}: residual={current_resid:.5g},"
        f" conv={convergence:.3g}"
    )
    de_callback.iteration += 1
    return current_resid < target_residual


de_callback.iteration = 0

print("params before", n.items())

result = differential_evolution(
    residual_norm,
    bounds,
    strategy="best1bin",
    maxiter=2000,
    popsize=200,
    tol=1e-4,
    mutation=(0.5, 1.0),
    recombination=0.7,
    polish=True,
    updating="deferred",
    workers=-1,
    callback=de_callback,
    disp=True,
)

# Apply optimized parameters
for value, token in zip(result.x, param_tokens):
    n[token] = value

site_spectra = n.current_spectrum
simulated_total = objective(n)
residual = d.data - simulated_total
print("final params", n.items())

with psd.figlist_var() as fl:
    fl.next("Differential evolution on pyspecdata trace")
    fl.plot(d, alpha=0.6, label="experimental")

    field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
    fl.plot(
        field_axis, simulated_total, alpha=0.9, label="differential evolution"
    )
    fl.plot(field_axis, residual, alpha=0.8, label="residual")
print("final params", n.items())
