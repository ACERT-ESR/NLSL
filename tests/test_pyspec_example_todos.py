import importlib.resources as resources
import warnings
import numpy as np
import nlsl
import pytest
from pathlib import Path
from pyspecdata.datadir import pyspec_config, getDATADIR

from nlsl.data import process_spectrum


DEPRECATION_MATCH = "load_raw_datafile spline preprocessing is deprecated"


def _load_processed_with_public_api(model, processed, label):
    stop = float(processed.start) + float(processed.step) * max(
        int(processed.y.size) - 1,
        0,
    )
    idx = model.generate_coordinates(
        float(processed.start),
        stop,
        int(processed.y.size),
    )
    model.data = processed.y
    model.name(label, spectrum=idx)
    model.noise(processed.noise, spectrum=idx)
    return idx


def _clear_example_registration():
    """Remove the nlsl_examples entry from pyspecdata's config caches."""

    pyspec_config.get_setting("nlsl_examples", section="ExpTypes")
    if (
        pyspec_config._config_parser is not None
        and pyspec_config._config_parser.has_option(
            "ExpTypes", "nlsl_examples"
        )
    ):
        pyspec_config._config_parser.remove_option("ExpTypes", "nlsl_examples")
    if "ExpTypes" in pyspec_config.config_vars:
        pyspec_config.config_vars["ExpTypes"].pop("nlsl_examples", None)


def test_register_nlsl_examples_sets_packaged_path():
    original_path = pyspec_config.get_setting(
        "nlsl_examples", section="ExpTypes"
    )
    _clear_example_registration()

    try:
        packaged_dir = None
        with resources.as_file(
            resources.files("nlsl").joinpath("__init__.py")
        ) as init_file:
            example_root = Path(init_file).parent / "examples"
            if (example_root / "230621_w0_10.DSC").exists():
                packaged_dir = example_root

        if not Path(getDATADIR("nlsl_examples")).exists():
            pyspec_config.set_setting(
                "ExpTypes", "nlsl_examples", str(packaged_dir)
            )

        stored_path = pyspec_config.get_setting(
            "nlsl_examples", section="ExpTypes"
        )
        assert Path(stored_path) == Path(str(packaged_dir))
    finally:
        if original_path is None:
            _clear_example_registration()
        else:
            pyspec_config.set_setting(
                "ExpTypes", "nlsl_examples", str(original_path)
            )


def test_max_points_property_matches_buffers():
    model = nlsl.nlsl()
    expected_points = model._core.expdat.data.shape[0] // max(
        model._core.expdat.nft.shape[0], 1
    )
    assert model.max_points == expected_points


def test_load_raw_datafile_spline_arguments_warn_deprecated():
    model = nlsl.nlsl()
    sample_path = Path(__file__).parent / "sampl1.dat"

    model.shift = False
    model.derivative_mode = 1
    with pytest.warns(DeprecationWarning, match=DEPRECATION_MATCH):
        model.load_raw_datafile(
            sample_path,
            nspline=10,
            bc_points=0,
            normalize=False,
        )


def test_processed_data_sequence_runs_warning_free():
    model = nlsl.nlsl()
    sample_path = Path(__file__).parent / "sampl1.dat"

    model.shift = False
    model.derivative_mode = 1
    processed = process_spectrum(
        sample_path,
        10,
        0,
        derivative_mode=model.derivative_mode,
        normalize=False,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        index = _load_processed_with_public_api(
            model,
            processed,
            sample_path.stem,
        )

    assert caught == []
    assert np.count_nonzero(model.data) > 0
    assert model.name(spectrum=index) == sample_path.stem
    assert model["rmsn"] == pytest.approx(processed.noise)


def test_load_raw_datafile_defaults_without_spline_arguments():
    model = nlsl.nlsl()
    sample_path = Path(__file__).parent / "sampl1.dat"

    model.shift = False
    model.derivative_mode = 1
    model.load_raw_datafile(sample_path, bc_points=0, normalize=False)

    data_span = model._core.expdat.data[: model.max_points]
    assert np.count_nonzero(data_span) > 0


def test_name_defaults_to_most_recent_spectrum():
    model = nlsl.nlsl()
    model.series("psi", (0.0, 90.0))
    first = model.generate_coordinates(0.0, 1.0, 4, reset=True)
    model.data = np.arange(4, dtype=float)
    model.name("first", spectrum=first)

    second = model.generate_coordinates(1.0, 2.0, 4)
    model.data = np.arange(4, dtype=float) + 10.0
    model.name("second", spectrum=second)

    assert model.name() == "second"
    assert model.name(spectrum=0) == "first"
    assert model.name(spectrum=1) == "second"

    model.name("renamed-first", spectrum=0)
    assert model.name(spectrum=0) == "renamed-first"


def test_rmsn_mapping_tracks_noise_scale_per_spectrum():
    model = nlsl.nlsl()
    model.series("psi", (0.0, 90.0))
    first = model.generate_coordinates(0.0, 1.0, 4, reset=True)
    model.data = np.arange(4, dtype=float)
    second = model.generate_coordinates(1.0, 2.0, 4)
    model.data = np.arange(4, dtype=float) + 10.0

    model["rmsn"] = [0.25, 0.0]

    assert first == 0
    assert model["rmsn"][0] == pytest.approx(0.25)
    assert model["rmsn"][1] == pytest.approx(1.0)
    assert model.noise(spectrum=0) == pytest.approx(0.25)
    assert model.noise(spectrum=1) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "nsite,new_value,expected",
    [
        (1, [0], 0),
        (2, [1, 0], np.array([1, 0])),
        (3, [1, 0, 1], np.array([1, 0, 1])),
    ],
)
def test_iscal_exposed_through_mapping(nsite, new_value, expected):
    model = nlsl.nlsl()
    model["nsite"] = nsite

    assert int(model["iscal"]) == 1

    model["iscal"] = new_value

    if np.isscalar(expected):
        assert int(model["iscal"]) == int(expected)
    else:
        assert np.array_equal(model["iscal"], expected)
    assert int(model._core.mspctr.iscglb) == 0
