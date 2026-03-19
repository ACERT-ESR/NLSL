import numpy as np
import pytest
import nlsl


def capture_vary_state(model):
    """Snapshot the Fortran vary tables for comparison."""
    parcom = model._core.parcom
    count = int(parcom.nprm)
    state = {
        "tag": [],
        "index": [],
        "parameter": [],
        "scale": [],
        "fdstep": [],
    }
    for position in range(count):
        label = parcom.tag[position]
        if isinstance(label, bytes):
            label = label.decode("ascii")
        state["tag"].append(label.strip().lower())
        state["index"].append(int(parcom.ixst[position]))
        state["parameter"].append(int(parcom.ixpr[position]))
        state["scale"].append(float(parcom.prscl[position]))
        state["fdstep"].append(float(parcom.xfdstp[position]))
    return state


def assert_two_site_g_tensor_expected_state(model, expected_scale, expected_fdstep):
    """Assert the vary table contains the expected two-site G-tensor state."""
    state = capture_vary_state(model)
    assert state["tag"] == ["gxx(1)", "gyy(2)"]
    assert state["index"] == [1, 2]
    assert state["parameter"] == [
        model.canonical_name("gxx")[2],
        model.canonical_name("gyy")[2],
    ]
    assert np.allclose(state["scale"], expected_scale)
    assert np.allclose(state["fdstep"], expected_fdstep)


@pytest.mark.parametrize("method", ["procline", "mapping"])
def test_fit_parameter_mapping_and_procline_reach_expected_vary_state(
    method, tmp_path, monkeypatch
):
    """Both API paths should produce the same expected vary-state mutation."""
    model = nlsl.nlsl()
    model.nsites = 2

    # Guard against cross-test contamination by proving the starting vary
    # table is empty before applying either parameterized update method.
    initial_state = capture_vary_state(model)
    assert initial_state["tag"] == []
    assert initial_state["index"] == []
    assert initial_state["parameter"] == []
    assert initial_state["scale"] == []
    assert initial_state["fdstep"] == []

    if method == "procline":
        model.procline(
            "vary gxx(1) minimum 1.0 maximum 3.0 scale 2.0 fdstep 1e-4"
        )
        first_state = capture_vary_state(model)
        assert len(first_state["tag"]) == 1
        assert first_state["tag"][0] == "gxx(1)"
        model.procline("vary gyy(2) scale 3.0")
        assert_two_site_g_tensor_expected_state(model, [2.0, 3.0], [1e-4, 1e-6])
    else:
        model.fit_params.vary["gxx"] = {
            "index": 1,
            "minimum": 1.0,
            "maximum": 3.0,
            "scale": 2.0,
            "fdstep": 1e-4,
        }
        first_state = capture_vary_state(model)
        assert len(first_state["tag"]) == 1
        assert first_state["tag"][0] == "gxx(1)"
        model.fit_params.vary["gyy"] = {"index": 2, "scale": 3.0}
        assert_two_site_g_tensor_expected_state(model, [2.0, 3.0], [1e-4, 1e-6])

    monkeypatch.chdir(tmp_path)
    model.procline("log vary_status")
    model.procline("status fit")
    model.procline("log end")
    output = (tmp_path / "vary_status.log").read_text()
    assert "GXX(1)" in output
    assert "GYY(2)" in output


@pytest.mark.parametrize("method", ["procline", "nlsl_parameters"])
def test_multisite_vary_updates_reach_expected_state(method):
    """Check multisite vary behavior for procline and NLSLParameters APIs."""
    model = nlsl.nlsl()
    model.nsites = 2

    # Guard against stale state from prior parameterized runs.
    initial_state = capture_vary_state(model)
    assert initial_state["tag"] == []
    assert initial_state["index"] == []
    assert initial_state["parameter"] == []
    assert initial_state["scale"] == []
    assert initial_state["fdstep"] == []

    if method == "procline":
        model.procline("vary gxx(1)")
        model.procline("vary gyy(2)")
    else:
        model.params["gxx_site0"].vary = True
        model.params["gyy_site1"].vary = True

    # Both approaches should apply the same default vary settings.
    assert_two_site_g_tensor_expected_state(model, [1.0, 1.0], [1e-6, 1e-6])
