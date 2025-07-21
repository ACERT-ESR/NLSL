import pytest

from run_example import run_example

EXAMPLES = [
    (1, [0.0404]),
    (2, [0.0331, 0.0513]),
    (3, [0.06113]),
    (4, [0.04001]),
    (5, [0.0714, 0.1592]),
]

@pytest.mark.parametrize("example,allowed", EXAMPLES)
def test_runexample(example, allowed):
    try:
        rel_rms = run_example(example, allowed_rel_rms=allowed)
    except ImportError as e:
        pytest.skip(f"required module missing: {e}")
    assert rel_rms and all(r < a * 1.01 for r, a in zip(rel_rms, allowed))
