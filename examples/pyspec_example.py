import pyspecdata as psd
import re
import nlsl

d = psd.find_file(
    re.escape("230621_w0_10.DSC"), exp_type="francklab_esr/romana"
)
d = d.chunk_auto("harmonic")["harmonic", 0]["phase", 0]
n = nlsl.nlsl()
with psd.figlist_var() as fl:
    fl.next("RS ESR figure")
    fl.plot(d)
    n.load_nddata(d, 0)
