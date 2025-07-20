import pytest

from runexample_first import run_example


def test_runexample_first():
    try:
        rel_rms = run_example()
    except ImportError as e:
        pytest.skip(f"required module missing: {e}")
    # run_example asserts internally, but also return value for clarity
    assert rel_rms is not None and rel_rms < 0.0404 * 1.01
