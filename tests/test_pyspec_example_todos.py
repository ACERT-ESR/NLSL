import importlib.resources as resources
import numpy as np
import nlsl
import warnings
from pathlib import Path
from pyspecdata.datadir import pyspec_config, getDATADIR


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
        packaged_dir = None
        with resources.as_file(resources.files("nlsl").joinpath("__init__.py")) as init_file:
            example_root = Path(init_file).parent / "examples"
            if (example_root / "230621_w0_10.DSC").exists():
                packaged_dir = example_root

        if packaged_dir is None:
            packaged_dir = Path(__file__).resolve().parent.parent / "examples"

        if not Path(getDATADIR("nlsl_examples")).exists():
            pyspec_config.set_setting("ExpTypes", "nlsl_examples", str(packaged_dir))

        stored_path = pyspec_config.get_setting("nlsl_examples", section="ExpTypes")
        assert Path(stored_path) == Path(str(packaged_dir))
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

def test_load_data_warns_for_spline_arguments():
    model = nlsl.nlsl()
    sample_path = Path(__file__).parent / "sampl1.dat"

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.load_data(
            sample_path,
            nspline=10,
            bc_points=0,
            shift=False,
            normalize=False,
        )
    assert any("backwards compatibility" in str(item.message) for item in caught)


def test_load_data_defaults_without_spline_arguments():
    model = nlsl.nlsl()
    sample_path = Path(__file__).parent / "sampl1.dat"

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("error")
        model.load_data(sample_path, bc_points=0, shift=False, normalize=False)

    data_span = model._core.expdat.data[: model.max_points]
    assert np.count_nonzero(data_span) > 0
