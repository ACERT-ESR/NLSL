import numpy as np
import pytest


def test_pyspecdata_import_and_nddata_probe():
    try:
        from pyspecdata import nddata
    except Exception as exc:
        pytest.fail(f"from pyspecdata import nddata failed: {exc!r}")

    try:
        dataset = nddata(np.r_[0:10], "t")
    except Exception as exc:
        pytest.fail(f"nddata(np.r_[0:10], 't') failed: {exc!r}")

    assert getattr(dataset, "dimlabels", None)
