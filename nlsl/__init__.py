from . import fortrancore as _fortrancore
import os
from pathlib import Path
from collections.abc import Mapping
import numpy as np
from .data import process_spectrum

_SPECTRAL_PARAMETER_NAMES = {
    "phase",
    "psi",
    "b0",
    "lb",
    "range",
}


def _decode_lpnam_array(array, width):
    """Return uppercase tokens extracted from a Fortran character array."""

    matrix = np.asarray(array, order="F")
    if matrix.ndim == 1:
        if matrix.dtype.kind == "S" and matrix.dtype.itemsize == width:
            raw_entries = matrix.tolist()
        else:
            raw_entries = (
                matrix.reshape(-1, width)
                .view(dtype="|S%d" % width)[:, 0]
                .tolist()
            )
    else:
        raw_entries = (
            matrix.T.reshape(-1, width)
            .view(dtype="|S%d" % width)[:, 0]
            .tolist()
        )
    # Decode each entry and normalise the tokens in a single pass so the
    # surrounding logic does not need to maintain an accumulator.
    tokens = [
        (raw.decode("ascii") if isinstance(raw, bytes) else str(raw))
        .upper()
        .rstrip()
        for raw in raw_entries
    ]
    return tuple(tokens)


def _match_parameter_token(token, entries):
    """Return the first index whose entry begins with *token*."""

    for index, candidate in enumerate(entries):
        if candidate and candidate.startswith(token):
            return index
    return None


class FitParameterVaryMapping(object):
    """Expose the active set of variable parameters through a mapping."""

    def __init__(self, owner, model):
        self._owner = owner
        self._core = owner._core
        self._model = model

    def _count(self):
        return int(self._core.parcom.nprm)

    def _entries(self, parameter):
        """Return ``(index, position)`` pairs for the requested parameter."""

        records = []
        parcom = self._core.parcom
        for position in range(self._count()):
            if int(parcom.ixpr[position]) == parameter:
                # Preserve the recorded index so callers can manage all slots
                # associated with the parameter in question.
                records.append((int(parcom.ixst[position]), position))
        return records

    def _is_spectrum_parameter(self, parameter):
        spectral = []
        for attr in (
            "iphase",
            "ipsi",
            "ilb",
            "ib0",
            "ifldi",
            "idfld",
            "irange",
        ):
            if hasattr(self._core.eprprm, attr):
                spectral.append(int(getattr(self._core.eprprm, attr)))
        if parameter in spectral:
            return True
        integral = []
        for attr in ("infld", "iiderv"):
            if hasattr(self._core.eprprm, attr):
                integral.append(int(getattr(self._core.eprprm, attr)))
        return parameter in integral

    def _index_limit(self, spectral):
        if spectral:
            limit = int(getattr(self._core.parcom, "nser", 0))
            nspc = int(getattr(self._core.expdat, "nspc", 0))
            if nspc > limit:
                limit = nspc
        else:
            limit = int(self._core.parcom.nsite)
        if limit <= 0:
            limit = 1
        return limit

    def _parameter_value(self, parameter, index):
        site = index if index > 0 else 1
        return float(self._core.getprm(parameter, site))

    def __getitem__(self, key):
        if not isinstance(key, str):
            raise KeyError(key)
        token = key.strip().lower()
        if not token:
            raise KeyError(key)
        if "(" in token:
            raise KeyError(
                "parameter keys no longer accept explicit indices; use the"
                " array returned for multi-site parameters instead"
            )

        # Require an attached model so canonical resolution can proceed.
        if self._model is None:
            raise RuntimeError("parameter resolution requires an nlsl model")
        try:
            _, ix, parameter = self._model.canonical_name(token)
        except KeyError:
            raise KeyError(f"{token} doesn't seem to be a parameter")
        entries = self._entries(parameter)
        if not entries:
            raise KeyError(key)
        parcom = self._core.parcom
        ordered = sorted(entries, key=lambda item: item[0])
        minima = []
        maxima = []
        scales = []
        steps = []
        indices = []
        for index_value, position in ordered:
            minima.append(float(parcom.prmin[position]))
            maxima.append(float(parcom.prmax[position]))
            scales.append(float(parcom.prscl[position]))
            steps.append(float(parcom.xfdstp[position]))
            indices.append(int(index_value))
        if len(indices) == 1:
            return {
                "minimum": minima[0],
                "maximum": maxima[0],
                "scale": scales[0],
                "fdstep": steps[0],
                "index": indices[0],
            }
        return {
            "minimum": np.array(minima, dtype=float),
            "maximum": np.array(maxima, dtype=float),
            "scale": np.array(scales, dtype=float),
            "fdstep": np.array(steps, dtype=float),
            "index": np.array(indices, dtype=int),
        }

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise KeyError(key)
        token = key.strip().lower()
        if not token:
            raise KeyError(key)
        if "(" in token:
            raise KeyError(
                "parameter keys no longer accept explicit indices; provide"
                " arrays when multiple sites are present"
            )
        if isinstance(value, bool):
            if not value:
                self.__delitem__(key)
                return
            config = {}
        elif value is None:
            config = {}
        elif isinstance(value, Mapping):
            config = value
        else:
            raise TypeError("vary requires bool or mapping values")
        # Require an attached model so canonical resolution can proceed.
        if self._model is None:
            raise RuntimeError("parameter resolution requires an nlsl model")
        try:
            _, ix, parameter = self._model.canonical_name(token)
        except KeyError:
            raise KeyError(f"{token} doesn't seem to be a parameter")
        spectral = self._is_spectrum_parameter(parameter)
        limit = self._index_limit(spectral)
        indices = None
        if "index" in config:
            indices = config["index"]
        elif "site" in config:
            indices = config["site"]
        elif "spectrum" in config:
            indices = config["spectrum"]

        if indices is None:
            if limit > 1:
                indices = list(range(1, limit + 1))
            else:
                indices = [0]
        elif np.isscalar(indices):
            indices = [int(indices)]
        else:
            indices = [int(item) for item in indices]

        for index in indices:
            if index < 0:
                raise ValueError("negative indices are not supported")
            if index > limit:
                raise ValueError("index out of range")

        count = len(indices)
        boundary_flags = np.zeros(count, dtype=int)
        minima = np.zeros(count, dtype=float)
        maxima = np.zeros(count, dtype=float)
        scales = np.ones(count, dtype=float)
        steps = np.empty(count, dtype=float)
        step_mask = np.zeros(count, dtype=bool)

        if "minimum" in config:
            # Accept either a scalar or per-site sequence for each option; when
            # a scalar is supplied, broadcast it across the active indices.
            if np.isscalar(config["minimum"]):
                minima[:] = float(config["minimum"])
            else:
                array = np.asarray(config["minimum"], dtype=float)
                if array.size != count:
                    raise ValueError(
                        "minimum entries must match the index list"
                    )
                minima[:] = array
            boundary_flags += 1
        if "maximum" in config:
            if np.isscalar(config["maximum"]):
                maxima[:] = float(config["maximum"])
            else:
                array = np.asarray(config["maximum"], dtype=float)
                if array.size != count:
                    raise ValueError(
                        "maximum entries must match the index list"
                    )
                maxima[:] = array
            boundary_flags += 2
        if "scale" in config:
            if np.isscalar(config["scale"]):
                scales[:] = float(config["scale"])
            else:
                array = np.asarray(config["scale"], dtype=float)
                if array.size != count:
                    raise ValueError("scale entries must match the index list")
                scales[:] = array
        if "fdstep" in config:
            if np.isscalar(config["fdstep"]):
                steps[:] = float(config["fdstep"])
            else:
                array = np.asarray(config["fdstep"], dtype=float)
                if array.size != count:
                    raise ValueError(
                        "fdstep entries must match the index list"
                    )
                steps[:] = array
            step_mask[:] = True

        entries = self._entries(parameter)
        ident = token.upper().split("(", 1)[0][:9]
        for existing_index, _ in entries:
            # Clear any stale configuration so the updated records replace the
            # full vary list for the parameter.
            _fortrancore.rmvprm(ix, existing_index, ident.ljust(30))

        for idx, index in enumerate(indices):
            base_value = self._parameter_value(parameter, index)
            if not step_mask[idx]:
                default_step = 1.0e-6
                step_value = default_step * base_value
                if abs(step_value) < float(np.finfo(float).eps):
                    step_value = default_step
            else:
                step_value = steps[idx]
            step_factor = step_value
            if abs(base_value) >= float(np.finfo(float).eps):
                step_factor = step_value / base_value
            _fortrancore.addprm(
                ix,
                index,
                int(boundary_flags[idx]),
                float(minima[idx]),
                float(maxima[idx]),
                float(scales[idx]),
                float(step_factor),
                ident.ljust(9),
            )

    def __delitem__(self, key):
        if not isinstance(key, str):
            raise KeyError(key)
        token = key.strip().lower()
        if not token:
            raise KeyError(key)
        if "(" in token:
            raise KeyError(
                "parameter keys no longer accept explicit indices; remove"
                " entries by requesting the base name"
            )
        # Require an attached model so canonical resolution can proceed.
        if self._model is None:
            raise RuntimeError("parameter resolution requires an nlsl model")
        try:
            _, ix, parameter = self._model.canonical_name(token)
        except KeyError:
            raise KeyError(f"{token} doesn't seem to be a parameter")
        entries = self._entries(parameter)
        if not entries:
            raise KeyError(key)
        ident = token.upper().split("(", 1)[0][:9]
        for index, _ in entries:
            _fortrancore.rmvprm(ix, index, ident.ljust(30))

    def __contains__(self, key):
        try:
            self[key]
        except Exception:
            return False
        return True

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return self._count()

    def keys(self):
        result = []
        for pos in range(self._count()):
            label = self._core.parcom.tag[pos]
            if isinstance(label, bytes):
                label = label.decode("ascii")
            label = label.rstrip()
            if label:
                result.append(label.lower())
        return result

    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def values(self):
        return [self[key] for key in self.keys()]

    def clear(self):
        keys = self.keys()
        for key in reversed(keys):
            self.__delitem__(key)

    def update(self, other):
        if isinstance(other, dict):
            items = other.items()
        else:
            items = other
        for key, val in items:
            self[key] = val


class TensorSymmetryMapping(object):
    """Dictionary-style view of the tensor symmetry flags."""

    _TOKEN_MAP = {
        "r": ("R", "irflg"),
        "diffusion": ("R", "irflg"),
        "g": ("G", "igflg"),
        "magnetic": ("G", "igflg"),
        "a": ("A", "iaflg"),
        "hyperfine": ("A", "iaflg"),
        "w": ("W", "iwflg"),
        "linewidth": ("W", "iwflg"),
    }
    _MODE_CODES = {
        "cartesian": 1,
        "cart": 1,
        "spherical": 2,
        "spher": 2,
        "axial": 3,
    }
    _MODE_LABELS = {
        0: "unset",
        1: "cartesian",
        2: "spherical",
        3: "axial",
    }
    _BASE_PARAMETER = {
        "W": "WXX",
        "G": "GXX",
        "A": "AXX",
        "R": "RX",
    }

    def __init__(self, model):
        self._model = model

    def _resolve_token(self, key):
        if not isinstance(key, str):
            raise KeyError(key)
        token = key.strip().lower()
        if not token:
            raise KeyError(key)
        if token not in self._TOKEN_MAP:
            raise KeyError(key)
        return self._TOKEN_MAP[token]

    def __getitem__(self, key):
        _, flag_name = self._resolve_token(key)
        if flag_name not in self._model._iepr_names:
            raise KeyError(key)
        row = self._model._iepr_names.index(flag_name)
        active_sites = max(self._model.nsites, 1)
        values = self._model._iparm[row, :active_sites]
        labels = [self._MODE_LABELS.get(int(val), "unknown") for val in values]
        if len(labels) == 1:
            return labels[0]
        return np.array(labels, dtype=object)

    def __setitem__(self, key, value):
        canonical, flag_name = self._resolve_token(key)
        if isinstance(value, str):
            mode_label = value
            indices = None
        elif isinstance(value, Mapping):
            if "mode" not in value:
                raise ValueError("tensor symmetry assignments need a mode")
            mode_label = value["mode"]
            if "site" in value:
                indices = value["site"]
            elif "index" in value:
                indices = value["index"]
            elif "spectrum" in value:
                indices = value["spectrum"]
            else:
                indices = None
        else:
            raise TypeError(
                "tensor symmetry values must be strings or mappings"
            )
        if not isinstance(mode_label, str):
            raise ValueError("tensor symmetry mode must be text")
        mode = mode_label.strip().lower()
        if mode not in self._MODE_CODES:
            raise ValueError(f"unknown tensor symmetry '{mode_label}'")
        flag_row = self._model._iepr_names.index(flag_name)
        _, _, table_index = self._model.canonical_name(
            self._BASE_PARAMETER[canonical]
        )
        component_start = table_index - 1
        total_sites = self._model._iparm.shape[1]
        if indices is None:
            site_indices = range(total_sites)
        elif np.isscalar(indices):
            index = int(indices)
            if index < 0:
                raise ValueError("site indices must be non-negative")
            if index == 0:
                site_indices = range(total_sites)
            else:
                if index > total_sites:
                    raise IndexError("site index exceeds available entries")
                site_indices = [index - 1]
        else:
            site_indices = []
            for entry in indices:
                index = int(entry)
                if index <= 0:
                    raise ValueError(
                        "explicit site indices must be positive integers"
                    )
                if index > total_sites:
                    raise IndexError("site index exceeds available entries")
                site_indices.append(index - 1)
        # Look up the converter and bookkeeping tables through the model so the
        # mapping stays aligned with the instance it manages.
        core = self._model._core
        converter = None
        if mode in ("spherical", "spher"):
            converter = core.tensym.tosphr
        elif mode == "axial":
            converter = core.tensym.toaxil
        else:
            converter = core.tensym.tocart
        ixx_table = core.parcom.ixx
        tag_table = core.parcom.tag
        for site_index in site_indices:
            # Guard against converting tensors whose components are currently
            # varied so the optimiser's bookkeeping stays consistent.
            for offset in range(3):
                var_index = int(ixx_table[component_start + offset, site_index])
                if var_index != 0:
                    label = tag_table[var_index - 1].decode("ascii").strip()
                    raise RuntimeError(
                        f"tensor symmetry unchanged: {canonical} tensor "
                        f"component is being varied as {label}"
                    )
            current_mode = int(self._model._iparm[flag_row, site_index])
            if current_mode == 0:
                self._model._iparm[flag_row, site_index] = self._MODE_CODES[mode]
                continue
            if current_mode == self._MODE_CODES[mode]:
                continue
            # Apply the requested conversion in-place so the tensor components
            # match the new symmetry before updating the bookkeeping flag.
            converter(
                self._model._fparm[
                    component_start : component_start + 3, site_index
                ],
                current_mode,
            )
            self._model._iparm[flag_row, site_index] = self._MODE_CODES[mode]
        self._model._last_site_spectra = None

    def __contains__(self, key):
        try:
            self._resolve_token(key)
        except KeyError:
            return False
        return True

    def keys(self):
        return list(self._TOKEN_MAP.keys())

    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def values(self):
        return [self[key] for key in self.keys()]


class fit_params(dict):
    """Mapping-like interface for adjusting NLSL fit parameters.

    Keys correspond to the options listed in ``nlshlp.txt`` lines 20–38.
    The values are mirrored directly to the low level ``lmcom`` module so
    that no ``procline`` call is needed.
    """

    def __init__(self, model=None):
        super().__init__()
        self._core = _fortrancore
        self._model = model
        self._fl_names = [
            n.decode("ascii").strip().lower()
            for n in self._core.lmcom.flmprm_name.tolist()
        ]
        self._il_names = [
            n.decode("ascii").strip().lower()
            for n in self._core.lmcom.ilmprm_name.tolist()
        ]
        self._vary = FitParameterVaryMapping(self, model)

    def __setitem__(self, key, value):
        key = key.lower()
        if key in self._fl_names:
            idx = self._fl_names.index(key)
            self._core.lmcom.flmprm[idx] = value
        elif key in self._il_names:
            idx = self._il_names.index(key)
            self._core.lmcom.ilmprm[idx] = value
        else:
            raise KeyError(key)
        super().__setitem__(key, value)

    def __getitem__(self, key):
        key = key.lower()
        if key in self._fl_names:
            return self._core.lmcom.flmprm[self._fl_names.index(key)]
        elif key in self._il_names:
            return self._core.lmcom.ilmprm[self._il_names.index(key)]
        raise KeyError(key)

    def __contains__(self, key):
        key = key.lower()
        return key in self._fl_names or key in self._il_names

    def __iter__(self):
        return iter(self.keys())

    def keys(self):
        return list(self._fl_names) + list(self._il_names)

    def items(self):
        return [(k, self[k]) for k in self.keys() if len(k) > 0]

    def values(self):
        return [self[k] for k in self.keys()]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def update(self, other):
        if isinstance(other, dict):
            items = other.items()
        else:
            items = other
        for k, v in items:
            self[k] = v

    @property
    def vary(self):
        """Dictionary-like view of the optimiser's variable parameter list."""
        return self._vary


class nlsl(object):
    """Dictionary-like interface to the NLSL parameters."""

    def __init__(self):
        _fortrancore.nlsinit()

        # Keep a reference to the shared Fortran bindings so helper mappings can
        # reach the active runtime state through this instance.
        self._core = _fortrancore

        self._fepr_names = [
            name.decode("ascii").strip().lower()
            for name in _fortrancore.eprprm.fepr_name.reshape(-1).tolist()
        ]
        # The ``fepr_name`` table leaves the start/step descriptors blank even
        # though ``parcom.fparm`` exposes the associated slots.  ``ipfind``
        # still recognises the historic ``FLDI``/``DFLD`` mnemonics, so attach
        # those labels to the trailing empty entries and leave any populated
        # tokens untouched.
        missing_fepr = [
            i for i, token in enumerate(self._fepr_names) if len(token) == 0
        ]
        for idx, label in zip(missing_fepr, ["fldi", "dfld"]):
            self._fepr_names[idx] = label

        self._iepr_names = [
            name.decode("ascii").strip().lower()
            for name in _fortrancore.eprprm.iepr_name.reshape(-1).tolist()
        ]
        # ``iepr_name`` omits several control flags (``IWFLG``/``IGFLG``/etc.)
        # that runfiles manipulate directly.  Only the blank slots need
        # patching, so mirror their canonical mnemonics onto the empty entries
        # without disturbing the documented ones.
        missing_iepr = [
            i for i, token in enumerate(self._iepr_names) if len(token) == 0
        ]
        for idx, label in zip(
            missing_iepr,
            ["iwflg", "igflg", "iaflg", "irflg", "jkmn", "jmmn", "ndim"],
        ):
            self._iepr_names[idx] = label
        # Decode the Fortran ``lpnam`` tables once so later lookups can run
        # without touching the legacy ``ipfind`` resolver.  The string length
        # metadata exposed by ``lpnam`` keeps the NumPy reshaping logic in sync
        # with the fixed-width Fortran character arrays.
        self._lpnam_tables = {}
        lpnam_module = _fortrancore.lpnam
        for key, attr, length_attr in [
            ("parnam", "parnam", "parnam_strlen"),
            ("iprnam", "iprnam", "iprnam_strlen"),
            ("alias1", "alias1", "alias_strlen"),
            ("alias2", "alias2", "alias_strlen"),
        ]:
            # Decode each lpnam table once so future lookups avoid touching the
            # Fortran character arrays directly.
            self._lpnam_tables[key] = _decode_lpnam_array(
                getattr(lpnam_module, attr),
                int(getattr(lpnam_module, length_attr)),
            )
        self._lpnam_tables["iwxx"] = None
        for position, token in enumerate(self._lpnam_tables["parnam"]):
            if token == "WXX":
                self._lpnam_tables["iwxx"] = position + 1
                break
        self._fparm = _fortrancore.parcom.fparm
        self._iparm = _fortrancore.parcom.iparm
        self.fit_params = fit_params(self)
        self.tensor_symmetry = TensorSymmetryMapping(self)
        self._last_layout = None
        self._last_site_spectra = None
        self._last_experimental_data = None
        self._weight_shape = (0, 0)
        self._explicit_field_start = False
        self._explicit_field_step = False

    @property
    def nsites(self) -> int:
        """Number of active sites."""
        self._sync_weight_matrix()
        return int(_fortrancore.parcom.nsite)

    @nsites.setter
    def nsites(self, value: int) -> None:
        # Propagate the site count to the core and refresh ``sfac`` so newly
        # exposed rows default to unity populations.
        _fortrancore.parcom.nsite = int(value)
        self._sync_weight_matrix()

    @property
    def nspec(self):
        """Number of active spectra."""
        self._sync_weight_matrix()
        return int(_fortrancore.expdat.nspc)

    @nspec.setter
    def nspec(self, value):
        # Keep the spectrum count and the cached weight matrix in lock-step.
        _fortrancore.expdat.nspc = int(value)
        self._sync_weight_matrix()

    @property
    def weights(self):
        """Expose the active ``/mspctr/sfac/`` populations.

        ``sfac`` stores one scale factor per (site, spectrum) pair inside a
        static ``MXSITE × MXSPC`` workspace that the optimiser and the
        single-point evaluator share.  ``_sync_weight_matrix`` keeps that table
        aligned with the current ``nsite``/``nspc`` counters so previously
        fitted populations remain intact when callers change the active site or
        spectrum counts.  Returning a column view yields a 1D vector for the
        common single-spectrum case, while the multi-spectrum case exposes an
        ``(nspc, nsite)`` view via ``transpose``.  Both paths hand out live
        views, so any edits immediately update the Fortran state.
        """

        matrix = self._sync_weight_matrix()
        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        if nsite <= 0 or nspc <= 0:
            return np.empty(0, dtype=float)
        active = matrix[:nsite, :nspc]
        if nspc == 1:
            return active[:, 0]
        return active.T

    @weights.setter
    def weights(self, values):
        """Overwrite the active portion of ``sfac`` with ``values``.

        ``sfac`` is shared between the optimiser and any ad-hoc spectrum
        evaluations.  When a caller provides new weights we zero the visible
        block and rewrite it with the supplied populations so that any entries
        outside the active range remain at the default value of one.
        """

        matrix = self._sync_weight_matrix()
        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        if nsite <= 0 or nspc <= 0:
            raise RuntimeError("weights require positive nsite and nspc")
        array = np.asarray(values, dtype=float)
        if array.ndim == 1:
            if nspc != 1:
                raise ValueError("1D weight vector requires a single spectrum")
            if array.size < nsite:
                raise ValueError("insufficient weight values supplied")
            matrix[:, 0] = 0.0
            matrix[:nsite, 0] = array[:nsite]
            return
        if array.shape[0] < nspc or array.shape[1] < nsite:
            raise ValueError("weight matrix shape mismatch")
        matrix[:, :] = 0.0
        matrix[:nsite, :nspc] = array[:nspc, :nsite].T

    def procline(self, val):
        """Process a line of a traditional format text NLSL runfile."""
        _fortrancore.procline(val)

    def fit(self):
        """Run the nonlinear least-squares fit using current parameters."""
        _fortrancore.fitl()
        return self._capture_state()

    @property
    def current_spectrum(self):
        """Evaluate the current spectral model without running a full fit.

        The returned array contains one row per site; population weights remain
        available through ``model.weights``.
        """
        ndatot = int(_fortrancore.expdat.ndatot)
        nspc = int(_fortrancore.expdat.nspc)
        if ndatot <= 0 or nspc <= 0:
            raise RuntimeError("no spectra have been evaluated yet")
        _fortrancore.iterat.iter = 1
        _fortrancore.single_point(1)
        return self._capture_state()

    @property
    def experimental_data(self):
        """Return the trimmed experimental traces from the most recent capture.

        The matrix is shaped as ``(number of spectra, point span)`` so it
        aligns with ``model.weights @ model.site_spectra``.  Each row contains
        the measured intensities for the corresponding recorded spectrum,
        zeroing any samples that fall outside that spectrum's active window.
        """

        if self._last_experimental_data is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_experimental_data

    def write_spc(self):
        """Write the current spectra to ``.spc`` files."""
        _fortrancore.wrspc()

    def load_data(
        self,
        data_id: str | os.PathLike,
        *,
        nspline: int,
        bc_points: int,
        shift: bool,
        normalize: bool = True,
        derivative_mode: int | None = None,
    ) -> None:
        """Load experimental data and update the Fortran state.

        The workflow mirrors the legacy ``datac`` command but avoids the
        Fortran file I/O path so that tests can exercise the data
        preparation logic directly from Python.
        """

        path = Path(data_id)
        if not path.exists():
            if path.suffix:
                raise FileNotFoundError(path)
            candidate = path.with_suffix(".dat")
            if not candidate.exists():
                raise FileNotFoundError(candidate)
            path = candidate

        token = str(path)
        base_name = token[:-4] if token.lower().endswith(".dat") else token
        mxpt = _fortrancore.expdat.data.shape[0]
        mxspc = _fortrancore.expdat.nft.shape[0]
        mxspt = mxpt // max(mxspc, 1)

        requested_points = int(nspline)
        if requested_points > 0:
            requested_points = max(4, min(requested_points, mxspt))

        nser = max(0, int(getattr(_fortrancore.parcom, "nser", 0)))
        normalize_active = bool(normalize or (self.nsites > 1 and nser > 1))

        mode = int(derivative_mode) if derivative_mode is not None else 1
        spectrum = process_spectrum(
            path,
            requested_points,
            int(bc_points),
            derivative_mode=mode,
            normalize=normalize_active,
        )

        idx, data_slice = self.generate_coordinates(
            int(spectrum.y.size),
            start=spectrum.start,
            step=spectrum.step,
            derivative_mode=mode,
            baseline_points=int(bc_points),
            normalize=normalize_active,
            nspline=requested_points,
            shift=shift,
            label=base_name,
        )

        eps = float(np.finfo(float).eps)
        _fortrancore.expdat.rmsn[idx] = (
            spectrum.noise if spectrum.noise > eps else 1.0
        )

        _fortrancore.expdat.data[data_slice] = spectrum.y
        _fortrancore.lmcom.fvec[data_slice] = spectrum.y

        if shift:
            _fortrancore.expdat.ishglb = 1

    def load_basis(self, identifier, spectrum=None, site=None):
        """Load a basis index file and assign it to optional targets."""

        if not isinstance(identifier, str):
            raise ValueError("basis identifiers must be strings")
        token = identifier.strip()
        if not token:
            raise ValueError("basis identifier cannot be empty")
        parts = [token]
        if spectrum is not None:
            parts.append("spectrum")
            parts.append(str(int(spectrum)))
        if site is not None:
            parts.append("site")
            parts.append(str(int(site)))
        command = " ".join(parts)
        if len(command) > 80:
            raise ValueError("basis command exceeds the 80 character limit")
        _fortrancore.basisc(command)
        self._last_site_spectra = None

    def search(self, parameter, site=1, **options):
        """Perform a one-dimensional search for *parameter* at *site*."""

        if not isinstance(parameter, str):
            raise ValueError("search expects the parameter name as text")
        token = parameter.strip()
        if not token:
            raise ValueError("parameter name cannot be empty")
        parts = [token]
        if site is not None:
            parts.append(str(int(site)))
        # ``srchc`` recognises a small, fixed vocabulary of modifiers for the
        # golden-section driver.  Mirror that list so callers receive a
        # meaningful error whenever an unsupported keyword slips through.
        allowed = {
            "xtol",
            "ftol",
            "step",
            "bound",
            "maxfun",
            "srange",
        }
        for key, value in options.items():
            label = key.strip().lower()
            if label not in allowed:
                raise ValueError(f"unsupported search option '{key}'")
            parts.append(label)
            if label == "maxfun":
                parts.append(str(int(value)))
            else:
                parts.append(str(float(value)))
        command = " ".join(parts)
        if len(command) > 80:
            raise ValueError("search command exceeds the 80 character limit")
        _fortrancore.srchc(command)
        self._last_site_spectra = None

    def series(self, parameter, values):
        """Configure a series definition for *parameter* with *values*."""

        if not isinstance(parameter, str):
            raise ValueError("series expects the parameter name as text")
        token = parameter.strip()
        if not token:
            raise ValueError("series parameter cannot be empty")
        if np.isscalar(values):
            raise ValueError("series requires a sequence of at least two values")
        sequence = [float(item) for item in values]
        if len(sequence) < 2:
            raise ValueError("series requires at least two values")
        parts = [token]
        parts.extend(str(entry) for entry in sequence)
        command = " ".join(parts)
        if len(command) > 80:
            raise ValueError("series command exceeds the 80 character limit")
        _fortrancore.series(command)
        self._last_site_spectra = None

    # -- mapping protocol -------------------------------------------------

    def __getitem__(self, key):
        key = key.lower()
        if key in ("nsite", "nsites"):
            return self.nsites
        if key in ("nspc", "nspec", "nspectra"):
            return self.nspec
        if key in ("sb0", "b0"):
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                try:
                    idx = self.parameter_index("b0")
                    if idx > 0:
                        return float(_fortrancore.getprm(idx, 1))
                except Exception:
                    pass
                if "b0" in self._fepr_names:
                    row = self._fepr_names.index("b0")
                    columns = max(self.nsites, 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        values = self._fparm[row, :columns]
                        if np.allclose(values, values[0]):
                            return float(values[0])
                        return values.copy()
                return 0.0
            values = _fortrancore.expdat.sb0[:nspc].copy()
            if np.allclose(values, values[0]):
                return float(values[0])
            return values
        if key in ("srng", "range"):
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                try:
                    idx = self.parameter_index("range")
                    if idx > 0:
                        return float(_fortrancore.getprm(idx, 1))
                except Exception:
                    pass
                if "range" in self._fepr_names:
                    row = self._fepr_names.index("range")
                    columns = max(self.nsites, 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        values = self._fparm[row, :columns]
                        if np.allclose(values, values[0]):
                            return float(values[0])
                        return values.copy()
                return 0.0
            values = _fortrancore.expdat.srng[:nspc].copy()
            if np.allclose(values, values[0]):
                return float(values[0])
            return values
        if key == "fldi":
            # ``fldi`` mirrors the field origin stored in ``expdat.sbi`` so
            # callers can recover the absolute coordinates used for the most
            # recent spectrum.
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                if "fldi" in self._fepr_names:
                    row = self._fepr_names.index("fldi")
                    columns = max(self.nsites, 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        values = self._fparm[row, :columns]
                        if np.allclose(values, values[0]):
                            return float(values[0])
                        return values.copy()
                return 0.0
            values = _fortrancore.expdat.sbi[:nspc].copy()
            if np.allclose(values, values[0]):
                return float(values[0])
            return values
        if key == "dfld":
            # ``dfld`` exposes the constant spacing between consecutive field
            # points.  When no spectra have been registered yet we fall back to
            # the cached floating-parameter table populated by the runfile.
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                if "dfld" in self._fepr_names:
                    row = self._fepr_names.index("dfld")
                    columns = max(self.nsites, 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        values = self._fparm[row, :columns]
                        if np.allclose(values, values[0]):
                            return float(values[0])
                        return values.copy()
                return 0.0
            values = _fortrancore.expdat.sdb[:nspc].copy()
            if np.allclose(values, values[0]):
                return float(values[0])
            return values
        if key == "ishft":
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                return 0
            values = _fortrancore.expdat.ishft[:nspc].copy()
            if np.all(values == values[0]):
                return int(values[0])
            return values
        if key in ("shift", "shft"):
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                return 0.0
            values = _fortrancore.expdat.shft[:nspc].copy()
            if np.allclose(values, values[0]):
                return float(values[0])
            return values
        if key in ("normalize_flags", "nrmlz"):
            nspc = int(_fortrancore.expdat.nspc)
            if nspc <= 0:
                return 0
            values = _fortrancore.expdat.nrmlz[:nspc].copy()
            if np.all(values == values[0]):
                return int(values[0])
            return values
        if key in ("weights", "weight", "sfac"):
            return self.weights
        # Resolve the canonical mnemonic together with its index so the
        # fallback paths behave like the historic resolver.
        try:
            canonical_key, res, base = self.canonical_name(key)
        except KeyError:
            canonical_key, res, base = None, 0, 0
        if res == 0:
            if canonical_key in self._iepr_names:
                idx = (
                    base - 1 if base else self._iepr_names.index(canonical_key)
                )
                vals = self._iparm[idx, : self.nsites]
                if np.all(vals == vals[0]):
                    return int(vals[0])
                return vals.copy()
            if canonical_key in self._fepr_names:
                idx = (
                    base - 1 if base else self._fepr_names.index(canonical_key)
                )
                vals = self._fparm[idx, : self.nsites]
                if np.allclose(vals, vals[0]):
                    return float(vals[0])
                return vals.copy()
            raise KeyError(key)
        if res > 100:
            idx = base - 1 if base else self._iepr_names.index(canonical_key)
            vals = self._iparm[idx, : self.nsites]
        else:
            vals = np.array([
                _fortrancore.getprm(res, i) for i in range(1, self.nsites + 1)
            ])
        if np.allclose(vals, vals[0]):
            return vals[0]
        return vals

    def __setitem__(self, key, v):
        key = key.lower()
        if key in ("nsite", "nsites"):
            self.nsites = int(v)
            return
        elif key in ("nspc", "nspec", "nspectra"):
            self.nspec = int(v)
            return
        elif key in ("weights", "weight", "sfac"):
            self.weights = v
            return
        expdat = _fortrancore.expdat
        if key == "fldi":
            # ``fldi`` holds the absolute starting field for each spectrum.
            # Keep both the ``expdat`` cache and the floating-parameter table
            # in sync so future ``range`` adjustments can reuse the stored
            # origin.
            values = np.atleast_1d(np.asarray(v, dtype=float))
            if values.size == 0:
                raise ValueError("fldi requires at least one value")
            nspc = int(expdat.nspc)
            if nspc <= 0:
                nspc = 1
            fill_count = min(max(nspc, 1), expdat.sbi.shape[0])
            expanded = np.empty(fill_count, dtype=float)
            expanded[:] = float(values[0])
            limit = min(values.size, fill_count)
            expanded[:limit] = values[:limit]
            expdat.sbi[:fill_count] = expanded
            self._explicit_field_start = True
            if "fldi" in self._fepr_names:
                row = self._fepr_names.index("fldi")
                columns = max(int(_fortrancore.parcom.nsite), 1)
                columns = min(columns, self._fparm.shape[1])
                if columns > 0:
                    for col in range(columns):
                        if col < expanded.size:
                            self._fparm[row, col] = expanded[col]
                        else:
                            self._fparm[row, col] = expanded[0]
            self._last_site_spectra = None
            return
        if key == "dfld":
            # ``dfld`` records the field increment between points.  Preserve it
            # explicitly so synthetic spectra can reuse the converged sampling
            # without re-deriving it from the range and point count.
            values = np.atleast_1d(np.asarray(v, dtype=float))
            if values.size == 0:
                raise ValueError("dfld requires at least one value")
            nspc = int(expdat.nspc)
            if nspc <= 0:
                nspc = 1
            fill_count = min(max(nspc, 1), expdat.sdb.shape[0])
            expanded = np.empty(fill_count, dtype=float)
            expanded[:] = float(values[0])
            limit = min(values.size, fill_count)
            expanded[:limit] = values[:limit]
            expdat.sdb[:fill_count] = expanded
            self._explicit_field_step = True
            if "dfld" in self._fepr_names:
                row = self._fepr_names.index("dfld")
                columns = max(int(_fortrancore.parcom.nsite), 1)
                columns = min(columns, self._fparm.shape[1])
                if columns > 0:
                    for col in range(columns):
                        if col < expanded.size:
                            self._fparm[row, col] = expanded[col]
                        else:
                            self._fparm[row, col] = expanded[0]
            self._last_site_spectra = None
            return
        if key in ("b0", "sb0", "range", "srng"):
            values = np.atleast_1d(np.asarray(v, dtype=float))
            if values.size == 0:
                raise ValueError(f"{key} requires at least one value")

            nspc = int(expdat.nspc)
            if nspc < 0:
                nspc = 0

            canonical = "b0" if key in ("b0", "sb0") else "range"
            if canonical == "b0":
                fill_count = max(nspc, 1)
                if fill_count > expdat.sb0.shape[0]:
                    fill_count = expdat.sb0.shape[0]
                expanded = np.empty(fill_count, dtype=float)
                expanded[:] = float(values[0])
                limit = min(values.size, fill_count)
                expanded[:limit] = values[:limit]
                expdat.sb0[:fill_count] = expanded
                if "b0" in self._fepr_names:
                    row = self._fepr_names.index("b0")
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        if expanded.size >= columns:
                            self._fparm[row, :columns] = expanded[:columns]
                        else:
                            self._fparm[row, :columns] = expanded[0]
            else:
                fill_count = max(nspc, 1)
                if fill_count > expdat.srng.shape[0]:
                    fill_count = expdat.srng.shape[0]
                expanded = np.empty(fill_count, dtype=float)
                expanded[:] = float(values[0])
                limit = min(values.size, fill_count)
                expanded[:limit] = values[:limit]
                expdat.srng[:fill_count] = expanded
                if not self._explicit_field_start:
                    expdat.sbi[:fill_count] = (
                        expdat.sb0[:fill_count]
                        - 0.5 * expdat.srng[:fill_count]
                    )
                else:
                    self._explicit_field_start = True
                if not self._explicit_field_step:
                    steps = np.zeros(fill_count, dtype=float)
                    for spectrum in range(fill_count):
                        points = (
                            int(expdat.npts[spectrum])
                            if spectrum < expdat.npts.shape[0]
                            else 0
                        )
                        if points > 1:
                            steps[spectrum] = expdat.srng[spectrum] / float(
                                points - 1
                            )
                    expdat.sdb[:fill_count] = steps
                else:
                    self._explicit_field_step = True
                if "range" in self._fepr_names:
                    row = self._fepr_names.index("range")
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        if expanded.size >= columns:
                            self._fparm[row, :columns] = expanded[:columns]
                        else:
                            self._fparm[row, :columns] = expanded[0]

                if "fldi" in self._fepr_names:
                    row = self._fepr_names.index("fldi")
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        for col in range(columns):
                            if col < expdat.sbi.shape[0]:
                                self._fparm[row, col] = expdat.sbi[col]
                            else:
                                self._fparm[row, col] = expdat.sbi[0]
                if "dfld" in self._fepr_names:
                    row = self._fepr_names.index("dfld")
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        for col in range(columns):
                            if col < expdat.sdb.shape[0]:
                                self._fparm[row, col] = expdat.sdb[col]
                            else:
                                self._fparm[row, col] = expdat.sdb[0]

            update_geometry = False
            if "range" in self._fepr_names:
                update_geometry = True
            if update_geometry:
                start_row = self._fepr_names.index("range") + 1
                step_row = start_row + 1
                if start_row < self._fparm.shape[0]:
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        start_values = expdat.sbi[:columns]
                        if start_values.size >= columns:
                            self._fparm[start_row, :columns] = start_values[
                                :columns
                            ]
                        elif start_values.size > 0:
                            self._fparm[start_row, :columns] = start_values[0]
                        else:
                            self._fparm[start_row, :columns] = 0.0
                if step_row < self._fparm.shape[0]:
                    columns = max(int(_fortrancore.parcom.nsite), 1)
                    columns = min(columns, self._fparm.shape[1])
                    if columns > 0:
                        step_values = expdat.sdb[:columns]
                        if step_values.size >= columns:
                            self._fparm[step_row, :columns] = step_values[
                                :columns
                            ]
                        elif step_values.size > 0:
                            self._fparm[step_row, :columns] = step_values[0]
                        else:
                            self._fparm[step_row, :columns] = 0.0

            self._last_site_spectra = None
            return
        if key == "ishft":
            values = np.atleast_1d(np.asarray(v, dtype=int))
            if values.size == 0:
                raise ValueError("ishft requires at least one value")
            nspc = max(int(expdat.nspc), 1)
            filled = np.empty(nspc, dtype=np.int32)
            filled[:] = int(values[0])
            limit = min(values.size, nspc)
            filled[:limit] = values[:limit]
            expdat.ishft[:nspc] = filled
            self._last_site_spectra = None
            return
        if key in ("shift", "shft"):
            values = np.atleast_1d(np.asarray(v, dtype=float))
            if values.size == 0:
                raise ValueError("shift requires at least one value")
            nspc = max(int(expdat.nspc), 1)
            filled = np.empty(nspc, dtype=float)
            filled[:] = float(values[0])
            limit = min(values.size, nspc)
            filled[:limit] = values[:limit]
            expdat.shft[:nspc] = filled
            if np.any(filled != 0.0):
                expdat.ishglb = 1
            self._last_site_spectra = None
            return
        if key in ("normalize_flags", "nrmlz"):
            values = np.atleast_1d(np.asarray(v, dtype=int))
            if values.size == 0:
                raise ValueError("normalize flags require at least one value")
            nspc = max(int(expdat.nspc), 1)
            filled = np.empty(nspc, dtype=np.int32)
            filled[:] = int(values[0])
            limit = min(values.size, nspc)
            filled[:limit] = values[:limit]
            expdat.nrmlz[:nspc] = filled
            self._last_site_spectra = None
            return
        iterinput = isinstance(v, (list, tuple, np.ndarray))
        try:
            canonical_key, res, base = self.canonical_name(key)
        except KeyError:
            canonical_key, res, base = None, 0, 0
        if res == 0:
            if canonical_key in self._iepr_names:
                idx = (
                    base - 1 if base else self._iepr_names.index(canonical_key)
                )
                if iterinput:
                    limit = min(len(v), self.nsites)
                    self._iparm[idx, :limit] = np.asarray(v[:limit], dtype=int)
                else:
                    self._iparm[idx, : self.nsites] = int(v)
                return
            elif canonical_key in self._fepr_names:
                idx = (
                    base - 1 if base else self._fepr_names.index(canonical_key)
                )
                values = np.asarray(v, dtype=float)
                if values.ndim == 0:
                    self._fparm[idx, : self.nsites] = float(values)
                else:
                    limit = min(values.size, self.nsites)
                    self._fparm[idx, :limit] = values[:limit]
                return
            else:
                raise KeyError(key)
        else:
            is_spectral = canonical_key in _SPECTRAL_PARAMETER_NAMES
            if res > 100:
                if iterinput:
                    limit = len(v)
                    if not is_spectral:
                        limit = min(limit, self.nsites)
                    else:
                        limit = min(limit, int(_fortrancore.expdat.nspc))
                    for site_idx in range(limit):
                        _fortrancore.setipr(
                            res, site_idx + 1, int(v[site_idx])
                        )
                else:
                    _fortrancore.setipr(res, 0, int(v))
            else:
                if iterinput:
                    limit = len(v)
                    if not is_spectral:
                        limit = min(limit, self.nsites)
                    else:
                        limit = min(limit, int(_fortrancore.expdat.nspc))
                    for site_idx in range(limit):
                        _fortrancore.setprm(
                            res, site_idx + 1, float(v[site_idx])
                        )
                else:
                    _fortrancore.setprm(res, 0, float(v))

    def __contains__(self, key):
        key = key.lower()
        if key in ("nsite", "nsites"):
            return True
        if key in self._fepr_names or key in self._iepr_names:
            return True
        try:
            self.parameter_index(key)
        except KeyError:
            return False
        return True

    def canonical_name(self, name):
        """Return canonical metadata and index codes for *name*.

        The resolver mirrors the legacy ``ipfind`` routine so callers always
        receive the canonical mnemonic, the Fortran index code, and the 1-based
        table position used by the vary machinery.  ``KeyError`` is raised when
        the name cannot be resolved.
        """

        raw = name.strip()
        if not raw:
            raise KeyError(name)
        token = raw.upper()
        if token in ("NSITE", "NSITES"):
            return "nsite", 0, 0
        # Search the floating-parameter table first so canonical mnemonics and
        # their historic abbreviations map to the corresponding positive codes.
        token_idx = _match_parameter_token(token, self._lpnam_tables["parnam"])
        if token_idx is None:
            # ``alias1`` encodes the traditional ``W1``/``G1``/etc. shorthands.
            alias_index = _match_parameter_token(
                token, self._lpnam_tables["alias1"]
            )
            if alias_index is None:
                # ``alias2`` hosts the extended ``WPRP``/``GPLL``-style
                # mnemonics.
                alias_index = _match_parameter_token(
                    token, self._lpnam_tables["alias2"]
                )
                if alias_index is None:
                    # Integral parameters are offset by 100 so they remain
                    # distinct from the floating parameter table when mirrored
                    # into the Fortran layer.
                    token_idx = _match_parameter_token(
                        token, self._lpnam_tables["iprnam"]
                    )
                    if token_idx is None:
                        lowered = raw.lower()
                        if lowered in self._fepr_names:
                            idx = self._fepr_names.index(lowered)
                            return lowered, idx + 1, idx + 1
                        else:
                            if lowered in self._iepr_names:
                                idx = self._iepr_names.index(lowered)
                                return lowered, 100 + idx + 1, idx + 1
                            else:
                                raise KeyError(name)
                    else:
                        return (
                            self._iepr_names[token_idx],
                            100 + token_idx + 1,
                            token_idx + 1,
                        )
                else:
                    if self._lpnam_tables["iwxx"] is None:
                        raise KeyError(name)
                    else:
                        base_idx = self._lpnam_tables["iwxx"] + alias_index - 1
                        if base_idx < 0 or base_idx >= len(self._fepr_names):
                            raise KeyError(name)
                        else:
                            index_code = -(
                                99
                                + self._lpnam_tables["iwxx"]
                                + alias_index
                                + 1
                            )
                            return (
                                self._fepr_names[base_idx],
                                index_code,
                                base_idx + 1,
                            )
            else:
                if self._lpnam_tables["iwxx"] is None:
                    raise KeyError(name)
                else:
                    base_idx = self._lpnam_tables["iwxx"] + alias_index - 1
                    if base_idx < 0 or base_idx >= len(self._fepr_names):
                        raise KeyError(name)
                    else:
                        index_code = 1 - (
                            self._lpnam_tables["iwxx"] + alias_index + 1
                        )
                        return (
                            self._fepr_names[base_idx],
                            index_code,
                            base_idx + 1,
                        )
        else:
            return (
                self._fepr_names[token_idx],
                token_idx + 1,
                token_idx + 1,
            )

    def parameter_index(self, name):
        """Return the Fortran-style index code for *name*."""

        return self.canonical_name(name)[1]

    def __iter__(self):
        return iter(self.keys())

    @property
    def layout(self):
        """Metadata describing the most recent spectral evaluation."""
        if self._last_layout is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_layout

    @property
    def field_axes(self):
        """Return uniformly spaced field axes for each active spectrum."""

        layout = self.layout
        counts = np.asarray(layout["npts"], dtype=int)
        if counts.size == 0:
            return tuple()
        starts = np.asarray(layout["sbi"], dtype=float).reshape(-1, 1)
        steps = np.asarray(layout["sdb"], dtype=float).reshape(-1, 1)
        span = int(counts.max())
        base = np.arange(span, dtype=float)
        grid = starts + steps * base
        return tuple(grid[idx, :count] for idx, count in enumerate(counts))

    @property
    def site_spectra(self):
        """Return the most recently evaluated site spectra."""
        if self._last_site_spectra is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_site_spectra

    def generate_coordinates(
        self,
        points: int,
        *,
        start: float,
        step: float,
        derivative_mode: int,
        baseline_points: int,
        normalize: bool,
        nspline: int,
        shift: bool = False,
        label: str | None = None,
        reset: bool = False,
    ) -> tuple[int, slice]:
        """Initialise the Fortran buffers for a uniformly spaced spectrum.

        Parameters mirror the coordinate bookkeeping that
        :meth:`load_data` performs after processing an experimental trace.
        The method allocates a fresh spectrum slot, configures the shared
        ``expdat`` metadata, and clears the backing work arrays without
        copying any intensity values.  It returns the spectrum index together
        with the slice into the flattened intensity arrays so callers may
        populate them manually.

        The *reset* flag mirrors the behaviour of the legacy ``datac``
        command: when ``True`` the spectrum counter and accumulated point
        count are cleared before initialising the new slot.  This is useful
        when synthesising spectra without loading any measured data first.
        """

        if points <= 0:
            raise ValueError("points must be positive")

        core = _fortrancore

        mxpt = core.expdat.data.shape[0]
        mxspc = core.expdat.nft.shape[0]
        mxspt = mxpt // max(mxspc, 1)

        if points > mxspt:
            raise ValueError("insufficient storage for spectrum")

        nspline = int(nspline)
        if nspline > 0:
            nspline = max(4, min(nspline, mxspt))

        if reset:
            core.expdat.nspc = 0
            core.expdat.ndatot = 0

        nspc = int(core.expdat.nspc)

        if hasattr(core.parcom, "nser"):
            nser = max(0, int(core.parcom.nser))
        else:
            nser = 0
        if nspc >= nser:
            nspc = 0
            core.expdat.ndatot = 0

        normalize_active = bool(normalize or (self.nsites > 1 and nser > 1))

        idx = nspc
        ix0 = int(core.expdat.ndatot)

        if idx >= mxspc:
            raise ValueError("Maximum number of spectra exceeded")
        elif ix0 + points > mxpt:
            raise ValueError("insufficient storage for spectrum")

        core.expdat.nspc = idx + 1
        core.expdat.ixsp[idx] = ix0 + 1
        core.expdat.npts[idx] = points
        core.expdat.sbi[idx] = float(start)
        core.expdat.sdb[idx] = float(step)
        core.expdat.srng[idx] = float(step) * max(points - 1, 0)
        core.expdat.ishft[idx] = 1 if shift else 0
        core.expdat.idrv[idx] = int(derivative_mode)
        core.expdat.nrmlz[idx] = 1 if normalize_active else 0
        core.expdat.shft[idx] = 0.0
        core.expdat.tmpshft[idx] = 0.0
        core.expdat.slb[idx] = 0.0
        core.expdat.sb0[idx] = 0.0
        core.expdat.sphs[idx] = 0.0
        core.expdat.spsi[idx] = 0.0

        core.expdat.rmsn[idx] = 1.0

        core.expdat.iform[idx] = 0
        core.expdat.ibase[idx] = int(baseline_points)

        power = 1
        while power < points:
            power *= 2
        core.expdat.nft[idx] = power

        data_slice = slice(ix0, ix0 + points)

        # ``single_point`` only reads the coordinate metadata and the site
        # storage arrays, so clearing the data buffer is sufficient here.
        core.expdat.data[data_slice] = 0.0

        if hasattr(core.mspctr, "spectr"):
            spectr = core.mspctr.spectr
            row_stop = min(ix0 + points, spectr.shape[0])
            spectr[ix0:row_stop, :] = 0.0
        if hasattr(core.mspctr, "wspec"):
            wspec = core.mspctr.wspec
            row_stop = min(ix0 + points, wspec.shape[0])
            wspec[ix0:row_stop, :] = 0.0
        if hasattr(core.mspctr, "sfac"):
            sfac = core.mspctr.sfac
            if idx >= sfac.shape[1]:
                raise ValueError("Maximum number of spectra exceeded")
            sfac[:, idx] = 1.0

        core.expdat.shftflg = 1 if shift else 0
        core.expdat.normflg = 1 if normalize_active else 0
        core.expdat.bcmode = int(baseline_points)
        core.expdat.drmode = int(derivative_mode)
        core.expdat.nspline = nspline
        core.expdat.inform = 0

        if label is None:
            label = "synthetic"
        encoded = label.encode("ascii", "ignore")[:30]
        core.expdat.dataid[idx] = encoded.ljust(30, b" ")

        trimmed = label.strip()
        window_label = f"{idx + 1:2d}: {trimmed}"[:19] + "\0"
        core.expdat.wndoid[idx] = window_label.encode("ascii", "ignore").ljust(
            20, b" "
        )

        core.expdat.ndatot = ix0 + points

        self._explicit_field_start = False
        self._explicit_field_step = False
        self._sync_weight_matrix()

        return idx, data_slice

    def set_data(self, data_slice, values):
        """Copy processed intensity values into the flattened data buffer."""

        start = data_slice.start
        stop = data_slice.stop
        expected = stop - start
        flat = np.asarray(values, dtype=float).reshape(-1)
        if flat.size != expected:
            raise ValueError("intensity vector length mismatch")
        _fortrancore.expdat.data[data_slice] = flat
        _fortrancore.lmcom.fvec[data_slice] = flat

    def set_site_weights(self, spectrum_index, weights):
        """Update the population weights for a specific spectrum index."""

        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        if spectrum_index < 0 or spectrum_index >= nspc:
            raise IndexError("spectrum index out of range")
        if nsite <= 0:
            raise RuntimeError("no active sites to weight")
        vector = np.asarray(weights, dtype=float).reshape(-1)
        if vector.size < nsite:
            raise ValueError("insufficient weight values supplied")
        target = _fortrancore.mspctr.sfac
        target[:, spectrum_index] = 0.0
        target[:nsite, spectrum_index] = vector[:nsite]

    def _capture_state(self):
        nspc = int(_fortrancore.expdat.nspc)
        ndatot = int(_fortrancore.expdat.ndatot)
        nsite = int(_fortrancore.parcom.nsite)

        spectra_src = _fortrancore.mspctr.spectr

        nspc = min(
            nspc,
            _fortrancore.expdat.ixsp.shape[0],
            _fortrancore.expdat.npts.shape[0],
            _fortrancore.expdat.sbi.shape[0],
            _fortrancore.expdat.sdb.shape[0],
        )
        nsite = min(nsite, spectra_src.shape[1])
        ndatot = min(ndatot, spectra_src.shape[0])

        starts = _fortrancore.expdat.ixsp[:nspc] - 1
        counts = _fortrancore.expdat.npts[:nspc].copy()
        windows = tuple(
            slice(int(start), int(start + count))
            for start, count in zip(starts, counts)
        )

        if windows:
            min_start = min(window.start for window in windows)
            max_stop = max(window.stop for window in windows)
        else:
            min_start = 0
            max_stop = 0

        relative_windows = tuple(
            slice(window.start - min_start, window.stop - min_start)
            for window in windows
        )

        # ``windows`` preserve the absolute indices used by the Fortran work
        # arrays so callers can recover the recorded experimental data, while
        # ``relative_windows`` remap the same ranges onto the trimmed spectra
        # returned below.

        self._last_layout = {
            "ixsp": starts,
            "npts": counts,
            "sbi": _fortrancore.expdat.sbi[:nspc].copy(),
            "sdb": _fortrancore.expdat.sdb[:nspc].copy(),
            "ndatot": max_stop - min_start,
            "nsite": nsite,
            "nspc": nspc,
            "windows": windows,
            "relative_windows": relative_windows,
            "origin": min_start,
        }

        span = max_stop - min_start
        if span > 0 and nsite > 0:
            trimmed = spectra_src[min_start:max_stop, :nsite]
            site_spectra = trimmed.swapaxes(0, 1)
        else:
            site_spectra = np.empty((nsite, 0), dtype=float)

        if span > 0 and nspc > 0:
            trimmed_exp = _fortrancore.expdat.data[min_start:max_stop].copy()
            stacked = np.zeros((nspc, span), dtype=float)
            for idx, window in enumerate(relative_windows):
                stacked[idx, window] = trimmed_exp[window]
        else:
            stacked = np.empty((nspc, 0), dtype=float)

        self._last_site_spectra = site_spectra
        self._last_experimental_data = stacked
        return self._last_site_spectra

    def _sync_weight_matrix(self):
        """Keep ``/mspctr/sfac/`` aligned with the active site/spectrum counts.

        ``sfac`` is declared in ``mspctr.f90`` as a static ``MXSITE × MXSPC``
        array.  ``nlsinit`` fills every element with ``1.0`` even when no fit
        is running, and the Fortran code never reinitialises the table after
        the initial call.  Changing ``nsite`` or ``nspc`` therefore only
        updates the integer counters; without extra work the exposed
        populations would keep whatever values happened to be in memory.  This
        helper mirrors the housekeeping performed by the Fortran data loaders
        so Python callers always see a predictable set of populations.
        """

        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        new_shape = (nsite, nspc)

        weights = _fortrancore.mspctr.sfac
        if new_shape == self._weight_shape:
            return weights

        if nsite <= 0 or nspc <= 0:
            # No active spectra or sites: blank the entire ``sfac`` matrix so a
            # subsequent resize starts from zero populations.
            weights[:, :] = 0.0
            self._weight_shape = new_shape
            return weights

        # Preserve the overlapping block so previously fitted populations stay
        # in place when callers expand or shrink the grid of sites/spectra.
        row_stop = min(self._weight_shape[0], nsite)
        col_stop = min(self._weight_shape[1], nspc)
        if row_stop > 0 and col_stop > 0:
            preserved = weights[:row_stop, :col_stop].copy()
        else:
            preserved = None

        # The Fortran initialisation routines seed ``sfac`` with ones, so new
        # rows/columns must do the same.  Reset the full table before restoring
        # any surviving populations.
        weights[:, :] = 1.0

        if preserved is not None:
            weights[:row_stop, :col_stop] = preserved

        self._weight_shape = new_shape
        return weights

    def keys(self):
        return list(self._fepr_names) + list(self._iepr_names)

    def items(self):
        return [(k, self[k]) for k in self.keys() if len(k) > 0]

    def values(self):
        return [self[k] for k in self.keys()]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def update(self, other):
        """Update multiple parameters at once."""
        assert isinstance(other, dict)
        for k, v in other.items():
            self[k] = v


# expose the class for creating additional instances
NLSL = nlsl

__all__ = [x for x in dir() if x[0] != "_"]
