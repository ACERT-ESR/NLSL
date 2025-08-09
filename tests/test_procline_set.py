import pytest
import nlsl
import numpy as np

FE_PARAMS = [
    "b0",
    "phase",
    "gib0",
    "gib2",
    "wxx",
    "wyy",
    "wzz",
    "gxx",
    "gyy",
    "gzz",
    "axx",
    "ayy",
    "azz",
    "rx",
    "ry",
    "rz",
    "pml",
    "pmxy",
    "pmzz",
    "djf",
    "djfprp",
    "oss",
    "psi",
    "alphad",
    "betad",
    "gammad",
    "alpham",
    "betam",
    "gammam",
    "c20",
    "c22",
    "c40",
    "c42",
    "c44",
    "lb",
    "dc20",
    "gamman",
    "cgtol",
    "shiftr",
    "shifti",
    #"range",
]

IE_PARAMS = [
    "in2",
    "ipdf",
    "ist",
    "ml",
    "mxy",
    "mzz",
    "lemx",
    "lomx",
    "kmn",
    "kmx",
    "mmn",
    "mmx",
    "ipnmx",
    "nort",
    "nstep",
    "nfield",
    "ideriv",
]

ALL_PARAMS = [
    (n, 100.0) if n == "range" else (n, 1.234) for n in FE_PARAMS
] + [(n, 1) for n in IE_PARAMS]


@pytest.mark.parametrize("key,val", ALL_PARAMS)
def test_procline_sets_module(key, val):
    n = nlsl.nlsl()
    before_iparm = n._iparm.copy()
    before_fparm = n._fparm.copy()
    if isinstance(val, float):
        n.procline(f"let {key} = {val}")
    else:
        n.procline(f"let {key} = {int(val)}")
    which_f_changed = [
        n._fepr_names[j]
        for j, v in enumerate(abs(np.sum(n._fparm - before_fparm, axis=1)))
        if v > 0
    ]
    which_i_changed = [
        n._iepr_names[j]
        for j, v in enumerate(abs(np.sum(n._iparm - before_iparm, axis=1)))
        if v > 0
    ]
    change_formost = ["ideriv", "lb", "phase", ""]
    which_changed = (set(which_i_changed) | set(which_f_changed)) - set(
        change_formost
    )  # these seem to change
    if key not in change_formost:
        assert len(which_changed) == 1, f"for {key}, these changed: {which_changed}"
    assert pytest.approx(n[key]) == val, (
        f"Problem setting {key}, which sets to {n[key]} rather than {val}\n"
        f"fparm: {which_f_changed}\n"
        f"iparm: {which_i_changed}\n"
    )
