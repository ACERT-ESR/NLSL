import numpy as np
import nlsl


# Ensure the lmfit-backed collection mirrors the live parameter values.
def test_parameter_values_track_nlsl_state():
    n = nlsl.nlsl()
    n["nsite"] = 2

    param = n.parameters["gxx_0"]
    code = n.parameter_index("gxx")
    original_other = nlsl.fortrancore.getprm(code, 2)

    new_value = param.value + 0.05
    param.value = new_value

    assert np.isclose(nlsl.fortrancore.getprm(code, 1), new_value)
    assert np.isclose(nlsl.fortrancore.getprm(code, 2), original_other)
    assert np.isclose(n.parameters["gxx_0"].value, new_value)


# Aliases should resolve to the same backing slot as canonical names.
def test_alias_and_canonical_share_storage():
    n = nlsl.nlsl()
    n["nsite"] = 1

    alias_key = "g1_0"
    canonical_key = "gxx_0"

    starting = n.parameters[canonical_key].value
    n.parameters[alias_key].value = starting + 0.02

    assert np.isclose(n.parameters[canonical_key].value, starting + 0.02)

    n.parameters[canonical_key].value = starting + 0.04
    assert np.isclose(n.parameters[alias_key].value, starting + 0.04)


# Adjusting the site count should grow and shrink the lmfit collection.
def test_parameter_entries_follow_site_count():
    n = nlsl.nlsl()

    n["nsite"] = 1
    assert "gxx_0" in n.parameters
    assert "gxx_1" not in n.parameters

    n["nsite"] = 2
    assert "gxx_1" in n.parameters

    n["nsite"] = 1
    assert "gxx_1" not in n.parameters


# Integer-only entries should not show up in the lmfit parameter collection.
def test_integer_parameters_are_skipped():
    n = nlsl.nlsl()
    n["nsite"] = 1

    assert "ndim_0" not in n.parameters
