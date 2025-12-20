import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config
from numpy import r_

if not (Path(psd.getDATADIR("nlsl_examples")).exists()):
    # if we haven't registered the example directory then register it
    # ☐ TODO -- we need to pull the example directory that's stored from the
    # packaging information, not just rely on the fact that it's the current
    # directory.  That way, we can also register this directory, e.g. from
    # within the tests (and we should check in the tests that we can remove
    # this directory from the pyspec config, and then register it)
    pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(Path.cwd()))

d = psd.find_file(re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples")
d.set_units(
    "$B_0$", None
)  # just for now, because I'm not prepared to deal with the weirdness, yet
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

# ☐ TODO: the st
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
    canonical = n.canonical_name(token)[0]
    key = f"{canonical}_0"
    n.parameters[key].vary = True

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
    field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
    fl.plot(field_axis, simulated_total, alpha=0.8, label="NLSL fit")
