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
# ☐ TODO -- while preserving this example, we ALSO should be able to use
# d.spline_lambda from pyspecdata to directly generate a spline, and then
# utilize that spline DIRECTLY, rather than either not making a spline,
# or relying on the internal spline mechanism.
# This definitely means we need to modify load_data so that it can accept
# a spline.
# This likely also means that the pyf file will need to be modified to
# allow us to directly access/supply the spline information from the
# output of spline_lambda
if d.data.shape[0] > max_points:
    divisor = d.shape["$B_0$"] // max_points + 1
    dB = np.diff(d["$B_0$"][r_[0, 1]]).item()
    d_orig_max = d.data.max()
    # (at 6σ, pretty much falls to zero between points
    d.convolve("$B_0$", dB / 6 * divisor)
    d = d["$B_0$", 0::divisor]
    # I'm a little confused b/c normalization should be preserved, and we end
    # up needing to do following
    d *= d_orig_max / d.data.max()
# }}}

# Provide reasonable starting parameters so the fit can run immediately.
# These match the values in the pyspecdata example and will be adjusted by the
# optimizer below.
n.update({
    "gxx": 2.0089,
    "gyy": 2.0021,
    "gzz": 2.0058,
    "in2": 2,
    "axx": 5.6,
    "ayy": 33.8,
    "azz": 5.3,
    "lemx": 6,
    "lomx": 5,
    "kmx": 4,
    "mmx": (2, 2),
    "rpll": np.log10(1.0e8),
    "rprp": 8.0,
    "gib0": 1.5,
})

# Load the nddata into the optimiser buffers without shifting the field.
n.load_nddata(d)

# Differential evolution will adjust a few relaxation-related parameters.  The
# parameter order here matches the bounds below and is reused when applying the
# optimized values back onto the model.
param_tokens = ["rpll", "rprp", "gib0"]

# Bounds are kept modest so the global search remains quick for the example.
bounds = [
    (np.log10(1.0e7), np.log10(5.0e9)),
    (4.0, 12.0),
    (0.5, 4.0),
]

# Evaluate the spectrum for a proposed parameter set and return the residual
# norm against the experimental trace for use by differential evolution.
def residual_norm(candidate):
    for value, token in zip(candidate, param_tokens):
        n[token] = value
    site_spectra = n.current_spectrum
    simulated_total = np.squeeze(n.weights @ site_spectra)
    return np.linalg.norm(d.data - simulated_total)

# Run a small number of iterations to demonstrate the search flow without
# making the example overly time-consuming.
result = differential_evolution(
    residual_norm,
    bounds,
    maxiter=8,
    popsize=8,
    polish=False,
    updating="deferred",
)

# Apply the optimized parameters back to the model to generate the final curves
# for plotting.
for value, token in zip(result.x, param_tokens):
    n[token] = value
site_spectra = n.current_spectrum
simulated_total = np.squeeze(n.weights @ site_spectra)
residual = d.data - simulated_total

with psd.figlist_var() as fl:
    fl.next("Differential evolution on pyspecdata trace")
    fl.plot(d, alpha=0.6, label="experimental")

    field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
    fl.plot(field_axis, simulated_total, alpha=0.9, label="differential evolution")
    fl.plot(field_axis, residual, alpha=0.8, label="residual")
    fl.show_legend()
