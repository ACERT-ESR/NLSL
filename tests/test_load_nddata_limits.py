import numpy as np
import pytest
import nlsl

pytest.importorskip("pyspecdata", exc_type=ImportError, reason="pyspecdata is required")
from pyspecdata.core import nddata


def test_load_nddata_reports_capacity_limit():
    """Oversized nddata inputs should report the storage limit clearly."""

    model = nlsl.nlsl()
    max_points = model._core.expdat.data.shape[0] // max(
        model._core.expdat.nft.shape[0], 1
    )
    points = max_points + 1
    fields = np.arange(points, dtype=float)
    dataset = nddata(np.zeros(points), [points], ["field"])
    dataset.setaxis("field", fields)

    with pytest.raises(ValueError) as excinfo:
        model.load_nddata(dataset)

    message = str(excinfo.value)
    assert str(int(points)) in message
    assert str(int(max_points)) in message


def test_load_nddata_defaults_to_no_shift():
    """The shift flag should default to the unshifted behaviour."""

    model = nlsl.nlsl()
    points = 4
    fields = np.linspace(0.0, 3.0, points)
    dataset = nddata(np.linspace(-1.0, 1.0, points), [points], ["field"])
    dataset.setaxis("field", fields)

    model.load_nddata(dataset)

    assert int(model._core.expdat.ishft[0]) == 0
