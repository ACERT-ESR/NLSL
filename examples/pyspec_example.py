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
with psd.figlist_var() as fl:
    fl.next("RS ESR figure")
    fl.plot(d)
    # ☐ TODO -- the desired behavior is that this should load, but it gives an
    # error.  Possibly, we need to reduce the number of datapoints, but
    # load_nddata should inform us of exactly this fact, and how many points
    # there are in the dimension, and what it must be reduced to.
    # ☐ TODO -- make shift a kwarg (default 0) of load_nddata, rather than a
    # required arg
    n.load_nddata(d, 0)
    # ☐ TODO -- now, we want to run a fit of the loaded data.  MOMD shouldn't
    # be required, so use any of the examples that doesn't use MOMD as an
    # example
