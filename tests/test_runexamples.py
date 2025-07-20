import pytest

from run_example import run_example

EXAMPLES = [
    (1, 0.0404),
    (2, 0.04585),
    (3, 0.06113),
    (4, 0.04010),
    (5, 0.11),
]

@pytest.mark.parametrize("example,allowed", EXAMPLES)
def test_runexample(example, allowed):
    try:
        rel_rms = run_example(example, allowed_rel_rms=allowed)
    except ImportError as e:
        pytest.skip(f"required module missing: {e}")
    assert rel_rms is not None and rel_rms < allowed * 1.01
