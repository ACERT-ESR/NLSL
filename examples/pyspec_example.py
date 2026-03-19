import os
import importlib.resources as resources
import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config
from matplotlib.pyplot import figure
from numpy import r_

if not Path(psd.getDATADIR("nlsl_examples")).exists():
    # Register the packaged examples directory.  We materialize a file that we
    # know is part of the wheel (``__init__.py``) to locate the installed
    # package root; this works with meson and other editable loaders that
    # refuse to materialize directories.  If the packaged DSC is unavailable,
    # fall back to the source tree copy.
    with resources.as_file(
        resources.files("nlsl").joinpath("__init__.py")
    ) as installed_loc:
        example_root = Path(installed_loc).parent.parent / "examples"
        if (example_root / "230621_w0_10.DSC").exists():
            packaged_dir = example_root
    pyspec_config.set_setting(
        "ExpTypes", "nlsl_examples", str(packaged_dir.absolute())
    )

d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
d.set_units(
    "$B_0$", None
)  # just for now, because I'm not prepared to deal with the weirdness, yet
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
d.name("230621_w0_10.DSC")
n = nlsl.nlsl()

field_axis = d[d.dimlabels[0]]
max_points = n.max_points
# {{{ we use convolution to downsample the data
if field_axis.size > max_points:
    divisor = d.shape["$B_0$"] // max_points + 1
    dB = np.diff(d["$B_0$"][r_[0, 1]]).item()
    d_orig_max = d.max()
    # (at 6σ, pretty much falls to zero between points
    d.convolve("$B_0$", dB / 6 * divisor)
    d = d["$B_0$", 0::divisor]
    # The convolution/downsampling step preserves area better than peak
    # height, so rescale the decimated trace back to the original maximum for
    # easier visual comparison.
    d *= d_orig_max / d.max()
# }}}

# Provide reasonable starting parameters so the fit can run immediately.
n["gxx"] = 2.0089
n["gyy"] = 2.0021
n["gzz"] = 2.0058
n["in2"] = 2
n["axx"] = 5.6
n["ayy"] = 33.8
n["azz"] = 5.3
n["lemx"] = 6
n["lomx"] = 5
n["kmx"] = 4
n["mmx"] = (2, 2)
n["rpll"] = np.log10(1.0e8)
n["rprp"] = 8.0
n["gib0"] = 1.5

for token in ("rpll", "rprp", "gib0"):
    n.fit_params.vary[token] = True

n.fit_params["maxitr"] = 20
n.fit_params["maxfun"] = 400
n.fit_params["ftol"] = 1.0e-3
n.fit_params["xtol"] = 1.0e-3

figure("RS ESR fit example")
psd.plot(d, alpha=0.45, label="experimental")
# ``normalize=False`` matches the old Fortran ``NONORM`` path: ``datac``
# would leave ``nrmlz=0`` for this spectrum, so the loader keeps the supplied
# nddata amplitudes and skips the integral-based rescaling step in ``getdat``.
n.load_nddata(d, normalize=False)

# Run a quick fit using the single-site parameters above so the
# least-squares weights are updated before plotting the nddata outputs.
n.fit()

# ``current_spectrum`` now carries the field axis and site labels directly.
psd.plot(n.current_spectrum, alpha=0.35, label="NLSL sites")
psd.plot(n.weights @ n.current_spectrum, alpha=0.8, label="NLSL fit")
