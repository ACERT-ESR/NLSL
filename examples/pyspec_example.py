import importlib.resources as resources
import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config
from numpy import r_

if not (Path(psd.getDATADIR("nlsl_examples")).exists()):
    # if we haven't registered the example directory then register it
    with resources.as_file(resources.files("nlsl").joinpath("examples")) as packaged_dir:
        target_dir = packaged_dir

    if not target_dir.exists():
        target_dir = Path(__file__).resolve().parent

    pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(target_dir))

d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
d.set_units(
    "$B_0$", None
)  # just for now, because I'm not prepared to deal with the weirdness, yet
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
n = nlsl.nlsl()

spline_func = d.spline_lambda()
field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
max_points = n.max_points
# {{{ we use convolution to downsample the data
# while preserving this example, we also demonstrate using a pyspecdata spline
# directly to keep the sampling under the solver's buffer limit.
if field_axis.size > max_points:
    divisor = d.shape["$B_0$"] // max_points + 1
    dB = np.diff(d["$B_0$"][r_[0, 1]]).item()
    d_orig_max = d.data.max()
    # (at 6σ, pretty much falls to zero between points
    d.convolve("$B_0$", dB / 6 * divisor)
    d = d["$B_0$", 0::divisor]
    # I'm a little confused b/c normalization should be preserved, and we end
    # up needing to do following
    d *= d_orig_max / d.data.max()
    spline_axis = np.linspace(field_axis[0], field_axis[-1], max_points)
else:
    spline_axis = field_axis
# }}}

# Provide reasonable starting parameters so the fit can run immediately.
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

for token in ("rpll", "rprp", "gib0"):
    n.fit_params.vary[token] = True

for key, value in {
    "maxitr": 20,
    "maxfun": 400,
    "ftol": 1.0e-3,
    "xtol": 1.0e-3,
}.items():
    n.fit_params[key] = value
with psd.figlist_var() as fl:
    fl.next("RS ESR figure")
    fl.plot(d)
    # Load the spline-evaluated nddata into the optimiser buffers without
    # shifting the field.
    n.load_spline(spline_func, spline_axis, bc_points=0, shift=False)

    # Run a quick fit using the single-site parameters above.
    site_spectra = n.fit()
    simulated_total = np.squeeze(n.weights @ site_spectra)

    # Overlay the simulated spectrum on the experimental trace.
    field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
    fl.plot(field_axis, simulated_total, alpha=0.8, label="NLSL fit")
