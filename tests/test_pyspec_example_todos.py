import importlib.resources as resources
import os
import numpy as np
import nlsl
import pyspecdata as psd
from pathlib import Path
from pyspecdata import nddata
from pyspecdata.datadir import pyspec_config


def _clear_example_registration():
    """Remove the nlsl_examples entry from pyspecdata's config caches."""

    pyspec_config.get_setting("nlsl_examples", section="ExpTypes")
    if pyspec_config._config_parser is not None and pyspec_config._config_parser.has_option(
        "ExpTypes", "nlsl_examples"
    ):
        pyspec_config._config_parser.remove_option("ExpTypes", "nlsl_examples")
    if "ExpTypes" in pyspec_config.config_vars:
        pyspec_config.config_vars["ExpTypes"].pop("nlsl_examples", None)


def test_register_nlsl_examples_sets_packaged_path():
    original_path = pyspec_config.get_setting("nlsl_examples", section="ExpTypes")
    _clear_example_registration()

    try:
        example_root = resources.files("nlsl") / "examples"
        try:
            target_dir = Path(os.fspath(example_root))
        except TypeError:
            target_dir = None
        if target_dir is None or not target_dir.exists():
            target_dir = Path(__file__).resolve().parent.parent / "examples"

        if not Path(psd.getDATADIR("nlsl_examples")).exists():
            pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(target_dir))

        stored_path = pyspec_config.get_setting("nlsl_examples", section="ExpTypes")
        assert Path(stored_path) == target_dir
    finally:
        if original_path is None:
            _clear_example_registration()
        else:
            pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(original_path))


def test_max_points_property_matches_buffers():
    model = nlsl.nlsl()
    expected_points = model._core.expdat.data.shape[0] // max(
        model._core.expdat.nft.shape[0], 1
    )
    assert model.max_points == expected_points

def test_load_spline_accepts_pyspecdata_lambda():
    model = nlsl.nlsl()
    target_points = min(model.max_points, 10)
    field_axis = np.linspace(3360.0, 3370.0, target_points)

    dataset = nddata(np.linspace(0.0, 1.0, target_points), "$B_0$")
    dataset.setaxis("$B_0$", field_axis)
    spline = dataset.spline_lambda()

    model.load_spline(
        spline,
        field_axis,
        bc_points=0,
        shift=False,
        normalize=False,
    )

    stored_points = model._core.expdat.data[:target_points]
    assert np.allclose(stored_points, dataset.data)

    if target_points > 1:
        expected_step = field_axis[1] - field_axis[0]
    else:
        expected_step = 0.0
    assert np.isclose(model._core.expdat.sdb[0], expected_step)
