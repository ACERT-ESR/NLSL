import pytest
import nlsl

FE_PARAMS = [
    'phase', 'gib0', 'gib2', 'wxx', 'wyy', 'wzz', 'gxx', 'gyy', 'gzz',
    'axx', 'ayy', 'azz', 'rx', 'ry', 'rz', 'pml', 'pmxy', 'pmzz', 'djf',
    'djfprp', 'oss', 'psi', 'alphad', 'betad', 'gammad', 'alpham', 'betam',
    'gammam', 'c20', 'c22', 'c40', 'c42', 'c44', 'lb', 'dc20', 'b0',
    'gamman', 'cgtol', 'shiftr', 'shifti', 'range'
]

IE_PARAMS = [
    'in2', 'ipdf', 'ist', 'ml', 'mxy', 'mzz', 'lemx', 'lomx', 'kmn',
    'kmx', 'mmn', 'mmx', 'ipnmx', 'nort', 'nstep', 'nfield', 'ideriv'
]

ALL_PARAMS = [(n, 1.234) for n in FE_PARAMS] + [(n, 1) for n in IE_PARAMS]

@pytest.mark.parametrize("key,val", ALL_PARAMS)
def test_procline_sets_module(key, val):
    n = nlsl.nlsl()
    if isinstance(val, float):
        n.procline(f"let {key} = {val}")
    else:
        n.procline(f"let {key} = {int(val)}")
    assert pytest.approx(n[key]) == val

