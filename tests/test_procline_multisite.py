import pytest
import nlsl


def test_procline_sets_individual_sites():
    n = nlsl.nlsl()
    n.procline("sites 2")
    n.procline("let gxx(1) = 1.1")
    n.procline("let gxx(2) = 2.2")
    n.procline("let wxx(1) = 3.3")
    n.procline("let wxx(2) = 4.4")

    gxx = n['gxx']
    wxx = n['wxx']

    assert pytest.approx(gxx[0]) == 1.1
    assert pytest.approx(gxx[1]) == 2.2
    assert pytest.approx(wxx[0]) == 3.3
    assert pytest.approx(wxx[1]) == 4.4
