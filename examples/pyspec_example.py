import os
import importlib.resources as resources
import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config
from numpy import r_

if not Path(psd.getDATADIR("nlsl_examples")).exists():
    # Register the packaged examples directory.  We materialize a file that we
    # know is part of the wheel (``__init__.py``) to locate the installed
    # package root; this works with meson and other editable loaders that
    # refuse to materialize directories.  If the packaged DSC is unavailable,
    # fall back to the source tree copy.
    packaged_dir = None
    try:
        with resources.as_file(resources.files("nlsl").joinpath("__init__.py")) as init_file:
            example_root = Path(init_file).parent / "examples"
            if (example_root / "230621_w0_10.DSC").exists():
                packaged_dir = example_root
    except FileNotFoundError:
        packaged_dir = None

    if packaged_dir is None:
        packaged_dir = Path(__file__).resolve().parent
    pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(packaged_dir))

d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
d.set_units(
    "$B_0$", None
)  # just for now, because I'm not prepared to deal with the weirdness, yet
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
n = nlsl.nlsl()

field_axis = d[d.dimlabels[0]]
max_points = n.max_points
# {{{ we use convolution to downsample the data
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
    # Load the nddata into the optimiser buffers without shifting the field.
    n.load_nddata(d)

    # Run a quick fit using the single-site parameters above.
    site_spectra = n.fit()
    simulated_total = np.squeeze(n.weights @ site_spectra)

    # Overlay the simulated spectrum on the experimental trace.
    field_axis = d[d.dimlabels[0]]
    fl.plot(field_axis, simulated_total, alpha=0.8, label="NLSL fit")
