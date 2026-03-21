"""
pySpecData example
==================

pySpecData already directly loads Bruker datafiles, so we show how to directly
load a data file and fit it with NLSL.

As part of this, we show how to generate a somewhat reasonable initial guess,
and to show what it looks like.

(As elsewhere, comments make use of fold markers where {{{ and }}} indicate the
beginning and end of the section that the comment pertains to.)
"""

import importlib.resources as resources
import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config
import matplotlib.pyplot as plt
from numpy import r_

# {{{ This section seems complicated, but simply makes sure that pyspecdata
# knows where the examples that ship with NLSL are located.
# Thanks to R. Shathy for this data.
if not Path(psd.getDATADIR("nlsl_examples")).exists():
    with resources.as_file(
        resources.files("nlsl").joinpath("__init__.py")
    ) as installed_loc:
        example_root = Path(installed_loc).parent.parent / "examples"
        if (example_root / "230621_w0_10.DSC").exists():
            packaged_dir = example_root
    pyspec_config.set_setting(
        "ExpTypes", "nlsl_examples", str(packaged_dir.absolute())
    )
# }}}

# first, we load the data from the registered directory (exp_type is a weird
# name for registered directory)
d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
field_axis = d[d.dimlabels[0]]  # the field axis is called "$B_0$" so that it
#                                looks pretty, so we just store it out of
#                                laziness.
# next, this is a multi-harmonic file, so we pull the first harmonic
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
# and we give it a name
d.name("230621: RM with $w_0=10$")
n = nlsl.nlsl()

max_points = n.max_points
# {{{ using pySpecData allows us to do more advanced pre-processing → e.g.
#     here, we use convolution averaging to downsample our data, rather than
#     relying on splines as in other examples.
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

# {{{ Provide reasonable starting parameters so the fit can run immediately.
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
n["betad"] = 15.0
n.shift = True
# }}}
plt.figure("RS ESR fit example", figsize=(10,6))
# Show the raw data -- note that the field axis, units, etc, are preloaded as
# part of the pySpecData nddata object.
psd.plot(d, alpha=0.45, label="experimental: " + d.name())
n.data = d.real  # convolution introduced zero imag, so need to return
#                 real

# ``current_spectrum`` now carries the field axis and site labels directly.
# this shows the initial guess

# Before ``fit()``, ``current_spectrum`` is just the raw per-site calculation
# from Fortran ``single_point``.  It is not auto-rescaled to the data.
# During ``fit()``, Fortran ``lfun`` calls ``sscale`` (or ``sshift`` when
# shifting is enabled) to determine least-squares scale factors and stores
# them in ``sfac``, exposed in Python as ``weights``.
# Those coefficients are not constrained to sum to one, so for a single site
# they simply provide the best-fit overall amplitude, and for multiple sites
# they provide the best-fit linear combination of site spectra.
# They are only user-fixed if Fortran's ``iscal`` flag is turned off for a
# site; otherwise fitting recomputes them.  Python exposes that flag as
# ``n['iscal']``.  If you never touch it, Fortran startup leaves
# ``iscal(i)=1`` for every site, so autoscaling is on by default.

n.weights = r_[600]

# {{{ set parameters related to the fit
n.fit_params["maxitr"] = 20
n.fit_params["maxfun"] = 400
n.fit_params["ftol"] = 1.0e-3
n.fit_params["xtol"] = 1.0e-3
# }}}

# {{{ With no entries in ``fit_params.vary``, calling ``fit()`` still
#     lets NLSL run the least-squares scale/shift update once.  That
#     gives us the best weight and field shift for the current
#     hand-entered parameters before we let the dynamic parameters move.
initial_site_spectra = n.fit()
# }}}

# {{{ when plotting, we use matrix multiplication (@) even when there is
#     only one spectrum (the weights still scales the overall amplitude
#     to match the spectrum)
psd.plot(
    n.weights @ initial_site_spectra,
    alpha=0.35,
    label="initial guess (shift/scale optimized)",
)
# }}}

# {{{ now enable fit of dynamic parameters, and fit again
for token in ("rpll", "rprp", "gib0"):
    n.fit_params.vary[token] = True
n.fit_params.vary["betad"] = {"minimum": 0.0, "maximum": 90.0}
n.fit()
# }}}
# {{{ ``current_spectrum`` now carries the field axis and site labels directly.
psd.plot(n.weights @ n.current_spectrum, alpha=0.8, label="NLSL fit")
# }}}
# {{{ now enable fit of tensor parameters, and fit again
for token in (a+2*b for b in ["x","y","z"] for a in ["a"]):
    n.fit_params.vary[token] = True
n.fit()
# releasing all at once gives craziness, so do A first, which has greater
# effect
for token in (a+2*b for b in ["x","y","z"] for a in ["g"]):
    n.fit_params.vary[token] = True
n.fit()
# }}}
# {{{ ``current_spectrum`` now carries the field axis and site labels directly.
psd.plot(n.weights @ n.current_spectrum, alpha=0.8, label="NLSL fit, with tensors")
# }}}

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.tight_layout()
plt.show()
