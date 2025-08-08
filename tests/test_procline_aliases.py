import pytest
import nlsl

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

@pytest.mark.parametrize("alias", ALIASES)
def test_procline_sets_alias(alias):
    n = nlsl.nlsl()
    n.procline(f"let {alias} = 1.234")
    canonical = n.canonical_name(alias)
    assert pytest.approx(n[canonical]) == 1.234
    assert pytest.approx(n[alias]) == 1.234

