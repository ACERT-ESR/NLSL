import numpy as np
import pyspecdata as psd
import re
import nlsl
from pathlib import Path
from pyspecdata.datadir import pyspec_config

pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(Path.cwd()))

d = psd.find_file(
    re.escape("230621_w0_10.DSC"), exp_type="nlsl_examples"
)
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
n = nlsl.nlsl()

# Trim the dataset to the workspace limit so it fits in the Fortran buffers.
max_points = n._core.expdat.data.shape[0] // max(n._core.expdat.nft.shape[0], 1)
if d.data.shape[0] > max_points:
    d = d[d.dimlabels[0], :max_points]

# Provide reasonable starting parameters so the fit can run immediately.
n.update(
    {
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
    }
)

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
    field_axis = np.asarray(d[d.dimlabels[0]], dtype=float)
    fl.plot(field_axis, simulated_total, alpha=0.8, label="NLSL fit")
