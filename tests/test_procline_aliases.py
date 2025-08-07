import pytest
import nlsl
from nlsl import _ipfind_wrapper

# alias names taken from lpnam.f90 (both alias1 and alias2 arrays)
ALIASES = [
    'w1', 'w2', 'w3',
    'g1', 'g2', 'g3',
    'a1', 'a2', 'a3',
    'rbar', 'n', 'nxy',
    'wprp', 'wpll',
    'gprp', 'gpll',
    'aprp', 'apll',
    'rprp', 'rpll',
]


def canonical_name(n, alias):
    res = _ipfind_wrapper(alias.upper())
    if res > 0:
        if res > 100:
            idx = res - 101
            return n._iepr_names[idx]
        return n._fepr_names[res - 1]
    if res > -100:
        idx = -res - 1
    else:
        idx = -res - 101
    return n._fepr_names[idx]


@pytest.mark.parametrize("alias", ALIASES)
def test_procline_sets_alias(alias):
    n = nlsl.nlsl()
    n.procline(f"let {alias} = 1.234")
    canonical = canonical_name(n, alias)
    assert pytest.approx(n[canonical]) == 1.234
    assert pytest.approx(n[alias]) == 1.234

