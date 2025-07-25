import pytest
import nlsl
from nlsl import fortrancore as fc

@pytest.mark.parametrize("key,val", [
    ("gxx", 2.0123),
    ("lemx", 6),
])
def test_procline_sets_module(key, val):
    n = nlsl.nlsl()
    if isinstance(val, float):
        n.procline(f"let {key} = {val}")
    else:
        n.procline(f"let {key} = {int(val)}")
    assert pytest.approx(n[key]) == val

