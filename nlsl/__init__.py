from . import fortrancore as _fortrancore
import os
from pathlib import Path
from collections.abc import Mapping
import numpy as np
import warnings
from .data import process_spectrum, read_ascii_spectrum

try:
    from pyspecdata import nddata

    _HAS_PYSPECDATA = True
except Exception:
    _HAS_PYSPECDATA = False

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
                var_index = int(
                    ixx_table[component_start + offset, site_index]
                )
                if var_index != 0:
                    label = tag_table[var_index - 1].decode("ascii").strip()
                    raise RuntimeError(
                        f"tensor symmetry unchanged: {canonical} tensor "
                        f"component is being varied as {label}"
                    )
            current_mode = int(self._model._iparm[flag_row, site_index])
            if current_mode == 0:
                self._model._iparm[flag_row, site_index] = self._MODE_CODES[
                    mode
                ]
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

        # Keep a reference to the shared Fortran bindings so helper mappings
        # can reach the active runtime state through this instance.
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
        self._last_site_spectra = None
        self._weight_shape = (0, 0)
        self._explicit_field_start = False
        self._explicit_field_step = False
        self.return_nddata = False
        # Global model-level controls for data loading and coordinate creation.
        self._shift = False
        self._derivative_mode = 1

    @property
    def max_points(self):
        """Maximum number of points available for a single spectrum."""

        mxpt = self._core.expdat.data.shape[0]
        mxspc = self._core.expdat.nft.shape[0]
        return int(mxpt // max(mxspc, 1))

    @property
    def shift(self):
        """Global shift flag used by data-loading and coordinate APIs."""

        return self._shift

    @shift.setter
    def shift(self, value):
        self._shift = bool(value)

    @property
    def derivative_mode(self):
        """Global derivative mode used by loading and coordinate APIs."""

        return self._derivative_mode

    @derivative_mode.setter
    def derivative_mode(self, value):
        if isinstance(value, bool):
            self._derivative_mode = 1 if value else 0
        else:
            self._derivative_mode = int(value)

    @property
    def nsites(self) -> int:
        """Number of active sites."""
        self._sync_weight_matrix()
        return int(_fortrancore.parcom.nsite)

    @nsites.setter
    def nsites(self, value: int) -> None:
        # Propagate the site count to the core and refresh ``sfac`` so newly
        # exposed rows default to unity scale factors.
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
        """Expose the active ``/mspctr/sfac/`` scale factors.

        ``sfac`` stores one least-squares coefficient per (site, spectrum)
        pair inside a static ``MXSITE × MXSPC`` workspace that the optimiser
        and the single-point evaluator share.  In Fortran, ``lfun`` updates
        those coefficients through ``sscale`` (or ``sshift`` when shifting is
        enabled) so the calculated spectra best match the experimental data.
        They therefore act as overall amplitudes or site contributions, and
        they are not constrained to sum to unity.  ``_sync_weight_matrix``
        keeps that table aligned with the current ``nsite``/``nspc`` counters
        so previously fitted scale factors remain intact when callers change
        the active site or spectrum counts.  Returning a column view yields a
        1D vector for the common single-spectrum case, while the
        multi-spectrum case exposes an ``(nspc, nsite)`` view via
        ``transpose``.  Both paths hand out live views, so any edits
        immediately update the Fortran state.
        """

        matrix = self._sync_weight_matrix()
        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        if nsite <= 0 or nspc <= 0:
            return np.empty(0, dtype=float)
        active = matrix[:nsite, :nspc]
        if nspc == 1:
            values = active[:, 0]
        else:
            values = active.T
        if not self.return_nddata:
            return values
        if not _HAS_PYSPECDATA:
            raise ImportError("pyspecdata is required for nddata outputs")
        array = np.asarray(values, dtype=float)
        if array.ndim == 1:
            wrapped = nddata(array.copy(), [array.shape[0]], ["sites"])
            wrapped.set_axis("sites", "#")
            return wrapped
        wrapped = nddata(
            array.copy(),
            [array.shape[0], array.shape[1]],
            ["spectrum", "sites"],
        )
        wrapped.set_axis("spectrum", "#")
        wrapped.set_axis("sites", "#")
        return wrapped

    @weights.setter
    def weights(self, values):
        """Overwrite the active portion of ``sfac`` with ``values``.

        ``sfac`` is shared between the optimiser and any ad-hoc spectrum
        evaluations.  When a caller provides new weights we zero the visible
        block and rewrite it with the supplied scale factors so that any
        entries outside the active range remain at the default value of one.
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
        """Run the nonlinear least-squares fit (``fitl``) using current
        parameters.

        Even when no model parameters are varying, Fortran still calls
        ``lfun`` once.  That means this method updates the least-squares
        scale factors in :attr:`weights` and, when shifting is enabled for a
        spectrum, updates the stored field shifts before returning the site
        spectra.
        """
        _fortrancore.fitl()
        return self._capture_state()

    @property
    def current_spectrum(self):
        """Evaluate the current spectral model without running ``fitl``.

        This performs the Fortran ``single_point`` calculation only.  It does
        not call ``lfun``, so it does not recompute least-squares scale
        factors or optimise any new data shift.  The returned array therefore
        contains the raw per-site calculated spectra at the currently stored
        shifts.  Use :meth:`fit` when you want NLSL to auto-rescale the model
        to the data before reading :attr:`weights`.

        The returned array contains one row per site unless ``return_nddata``
        is enabled, in which case the output is wrapped as labelled nddata.
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
        """Return experimental data on the shared trimmed fitting span.

        Returns
        -------
        ndarray
            A 2D array with shape ``(nspc, span)``.  Each row is padded with
            zeros outside the active region of that spectrum so the rows align
            with ``weights @ site_spectra`` in the legacy ndarray workflow.

        Notes
        -----
        This differs from :attr:`data`.  ``data`` returns each loaded
        spectrum by itself, either as a single trace or as a list of traces,
        preserving the original point counts and field axes.  By contrast,
        ``experimental_data`` repacks those same spectra onto one common
        trimmed span of the shared Fortran work arrays, which is convenient
        for matrix arithmetic against the similarly packed output of
        :attr:`current_spectrum` when ``return_nddata`` is ``False``.
        """

        nspc = int(_fortrancore.expdat.nspc)
        ndatot = int(_fortrancore.expdat.ndatot)
        if nspc <= 0 or ndatot <= 0:
            raise RuntimeError("no spectra have been evaluated yet")
        min_start = min(window.start for window in self.windows)
        max_stop = max(window.stop for window in self.windows)
        trimmed = _fortrancore.expdat.data[min_start:max_stop].copy()
        stacked = np.zeros((nspc, max_stop - min_start), dtype=float)
        for idx, window in enumerate(self.relative_windows):
            stacked[idx, window] = trimmed[window]
        return stacked

    def write_spc(self):
        """Write the current spectra to ``.spc`` files."""
        _fortrancore.wrspc()

    def load_raw_datafile(
        self,
        data_id,
        nspline=None,
        bc_points=None,
        normalize=True,
        preprocess=True,
    ):
        """Load a raw ASCII data file of experimental data and update the Fortran state.

        The workflow mirrors the legacy ``datac`` command, but performs the
        file read and optional preprocessing from Python before copying the
        prepared intensities into the active Fortran buffers.

        When ``preprocess`` is ``True`` this method calls
        :func:`nlsl.data.process_spectrum` before :meth:`generate_coordinates`
        and then copies intensities into :attr:`data`.

        To bypass spline/baseline/normalization preprocessing, set
        ``preprocess=False``.  In that mode the data is read directly from the
        ASCII file, must have a uniform field axis, and is loaded into
        :attr:`data` without calling :func:`nlsl.data.process_spectrum`.

        ``bc_points`` keeps the historical runfile naming: it mirrors the
        legacy ``bc`` argument used by ``data ... bc <points>``.
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
        mxspt = self.max_points

        spline_params_used = nspline is not None
        if spline_params_used:
            # NOTE TO AGENTS: leave this warning in place, in spite of
            # instructions about backwards compatibility.
            warnings.warn(
                "load_data spline arguments are retained only for backwards"
                " compatibility; prefer using pyspecdata or scipy to smooth"
                " the data before loading.",
                stacklevel=2,
            )

        requested_points = int(nspline) if nspline is not None else 0
        if requested_points > 0:
            requested_points = max(4, min(requested_points, mxspt))
        if bc_points is None:
            bc_points = 0

        normalize_active = bool(normalize)

        mode = int(self.derivative_mode)
        if preprocess:
            spectrum = process_spectrum(
                path,
                requested_points,
                int(bc_points),
                derivative_mode=mode,
                normalize=normalize_active,
            )
            spectrum_start = spectrum.start
            spectrum_stop = spectrum.start + spectrum.step * max(
                int(spectrum.y.size) - 1, 0
            )
            spectrum_y = spectrum.y
            spectrum_noise = spectrum.noise
            spectrum_norm = 1 if normalize_active else 0
            spectrum_nspline = requested_points
            spectrum_bcmode = int(bc_points)
        else:
            x_raw, y_raw = read_ascii_spectrum(path)
            if x_raw.size < 1:
                raise ValueError("spectrum contained no points")
            if x_raw.size > mxspt:
                raise ValueError("insufficient storage for spectrum")
            if x_raw.size > 1:
                steps = np.diff(x_raw)
                if not np.allclose(steps, steps[0]):
                    raise ValueError(
                        "preprocess=False requires a uniformly spaced field axis"
                    )
            spectrum_start = float(x_raw[0])
            spectrum_stop = float(x_raw[-1])
            spectrum_y = np.asarray(y_raw, dtype=float)
            spectrum_noise = float(np.std(spectrum_y))
            spectrum_norm = 0
            spectrum_nspline = 0
            spectrum_bcmode = 0

        # Keep Fortran metadata consistent with the loading options used for
        # this spectrum.
        _fortrancore.expdat.bcmode = spectrum_bcmode
        _fortrancore.expdat.nspline = spectrum_nspline
        _fortrancore.expdat.normflg = spectrum_norm
        _fortrancore.expdat.drmode = mode
        _fortrancore.expdat.shftflg = 1 if self.shift else 0

        idx = self.generate_coordinates(
            spectrum_start,
            spectrum_stop,
            int(spectrum_y.size),
            label=base_name,
        )

        eps = float(np.finfo(float).eps)
        _fortrancore.expdat.rmsn[idx] = (
            spectrum_noise if spectrum_noise > eps else 1.0
        )

        self.data = spectrum_y
        _fortrancore.expdat.nrmlz[idx] = spectrum_norm
        self.return_nddata = False

        if self.shift:
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
            raise ValueError(
                "series requires a sequence of at least two values"
            )
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
        if key == "iscal":
            nsite = int(_fortrancore.parcom.nsite)
            if nsite <= 0:
                return 1
            values = _fortrancore.mspctr.iscal[:nsite].copy()
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
        if key == "iscal":
            values = np.atleast_1d(np.asarray(v, dtype=int))
            if values.size == 0:
                raise ValueError("iscal requires at least one value")
            nsite = max(int(_fortrancore.parcom.nsite), 1)
            filled = np.empty(nsite, dtype=np.int32)
            filled[:] = int(values[0])
            limit = min(values.size, nsite)
            filled[:limit] = values[:limit]
            _fortrancore.mspctr.iscal[:nsite] = filled
            _fortrancore.mspctr.iscglb = 1 if np.all(filled != 0) else 0
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
    def field_axes(self):
        """Return uniformly spaced field axes for each active spectrum."""

        nspc = int(_fortrancore.expdat.nspc)
        counts = np.asarray(_fortrancore.expdat.npts[:nspc], dtype=int)
        if counts.size == 0:
            return tuple()
        starts = np.asarray(_fortrancore.expdat.sbi[:nspc], dtype=float).reshape(
            -1, 1
        )
        steps = np.asarray(_fortrancore.expdat.sdb[:nspc], dtype=float).reshape(
            -1, 1
        )
        span = int(counts.max())
        base = np.arange(span, dtype=float)
        grid = starts + steps * base
        return tuple(grid[idx, :count] for idx, count in enumerate(counts))

    @property
    def windows(self):
        """Return absolute slices for each spectrum in the shared buffers.

        Returns
        -------
        tuple of slice
            One Python ``slice`` per active spectrum.  Each slice is suitable
            for indexing the shared 1D Fortran work arrays such as
            ``expdat.data`` and ``lmcom.fvec``.

        Notes
        -----
        NLSL stores all experimental spectra in one flat workspace rather than
        in separate arrays.  Fortran records where each spectrum starts with
        ``ixsp(isp)`` and how many points it owns with ``npts(isp)``.
        ``windows`` is the direct 0-based Python translation of that layout,
        so ``self.windows[i]`` selects the actual segment of the shared memory
        that belongs to spectrum ``i``.
        """

        nspc = int(_fortrancore.expdat.nspc)
        counts = np.asarray(_fortrancore.expdat.npts[:nspc], dtype=int)
        if counts.size == 0:
            return tuple()
        starts = np.asarray(_fortrancore.expdat.ixsp[:nspc], dtype=int) - 1
        return tuple(
            slice(int(start), int(start + count))
            for start, count in zip(starts, counts)
        )

    @property
    def relative_windows(self):
        """Return per-spectrum slices into the shared trimmed span.

        Returns
        -------
        tuple of slice
            One Python ``slice`` per active spectrum, expressed relative to
            the smallest trimmed block that contains every active spectrum.

        Notes
        -----
        ``windows`` index the original shared Fortran arrays directly.
        ``relative_windows`` index a NumPy view of the trimmed block
        ``min_start:max_stop`` cut out of those arrays.  That makes them the
        convenient companion for code that first trims the shared storage and
        then wants per-spectrum views into the trimmed NumPy array without
        carrying the unused leading or trailing buffer space.  Here
        "trimmed" means that Python has rebased the shared Fortran storage
        onto the smallest contiguous NumPy block that still contains every
        active spectrum, so the first active point becomes index 0 in that
        temporary array.

        No fitting points are dropped by this remapping.  The actual Fortran
        fit uses the original shared-buffer layout described by ``ixsp`` and
        ``npts`` in ``expdat``, together with ``ndatot`` for the total stored
        point count.  ``relative_windows`` is only a Python indexing helper
        that re-expresses those same active regions relative to the trimmed
        NumPy block ``min_start:max_stop``.  Likewise, :attr:`current_spectrum`
        still contains the full active calculated span; these slices only tell
        Python where each individual spectrum lives within that packed result.
        ``windows`` are therefore the right slices when Python is indexing the
        original Fortran-backed storage directly, while ``relative_windows``
        are the right slices after Python has already created that trimmed,
        re-based NumPy block and now needs indices that are valid inside it.
        """

        if len(self.windows) == 0:
            return tuple()
        min_start = min(window.start for window in self.windows)
        return tuple(
            slice(window.start - min_start, window.stop - min_start)
            for window in self.windows
        )

    @property
    def site_spectra(self):
        """Return the most recently evaluated site spectra."""
        if self._last_site_spectra is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_site_spectra

    def generate_coordinates(
        self,
        start,
        stop,
        points,
        label=None,
        reset=False,
    ):
        """Initialise the Fortran buffers for a uniformly spaced spectrum.

        Use linspace-like ordering (start, stop, points) so this routine only
        describes the field axis and storage layout.

        """

        if points <= 0:
            raise ValueError("points must be positive")

        core = _fortrancore

        mxpt = core.expdat.data.shape[0]
        mxspc = core.expdat.nft.shape[0]
        mxspt = mxpt // max(mxspc, 1)

        if points > mxspt:
            raise ValueError("insufficient storage for spectrum")

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

        idx = nspc
        ix0 = int(core.expdat.ndatot)

        if idx >= mxspc:
            raise ValueError("Maximum number of spectra exceeded")
        if ix0 + points > mxpt:
            raise ValueError("insufficient storage for spectrum")

        start = float(start)
        stop = float(stop)
        if points == 1:
            step = 0.0
        else:
            step = (stop - start) / float(points - 1)

        core.expdat.nspc = idx + 1
        core.expdat.ixsp[idx] = ix0 + 1
        core.expdat.npts[idx] = points
        core.expdat.sbi[idx] = start
        core.expdat.sdb[idx] = step
        core.expdat.srng[idx] = step * max(points - 1, 0)
        core.expdat.ishft[idx] = 1 if self.shift else 0
        core.expdat.idrv[idx] = int(self.derivative_mode)
        core.expdat.drmode = int(self.derivative_mode)
        core.expdat.nrmlz[idx] = 0
        core.expdat.shft[idx] = 0.0
        core.expdat.tmpshft[idx] = 0.0
        core.expdat.slb[idx] = 0.0
        core.expdat.sb0[idx] = 0.0
        core.expdat.sphs[idx] = 0.0
        core.expdat.spsi[idx] = 0.0

        core.expdat.rmsn[idx] = 1.0
        core.expdat.iform[idx] = 0
        core.expdat.ibase[idx] = int(core.expdat.bcmode)

        power = 1
        while power < points:
            power *= 2
        core.expdat.nft[idx] = power

        data_slice = slice(ix0, ix0 + points)

        # Clear per-spectrum work buffers so each generated axis starts from
        # clean storage.
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

        core.expdat.shftflg = 1 if self.shift else 0
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

        return idx

    @property
    def data(self):
        """Return the loaded experimental traces spectrum-by-spectrum.

        A single loaded spectrum returns one array or one nddata.  Multiple
        spectra return a Python list so each trace keeps its own field axis
        and point count.

        Unlike :attr:`experimental_data`, this property does not repack
        spectra onto a shared trimmed span.  It simply returns the stored
        traces one spectrum at a time.
        """
        if len(self.windows) == 0:
            raise RuntimeError("no spectra have been allocated")
        if not self.return_nddata:
            traces = [
                _fortrancore.expdat.data[window].copy()
                for window in self.windows
            ]
            if len(traces) == 1:
                return traces[0]
            return traces
        if not _HAS_PYSPECDATA:
            raise ImportError("pyspecdata is required for nddata outputs")
        fields = self.field_axes
        traces = []
        for idx, window in enumerate(self.windows):
            encoded = _fortrancore.expdat.wndoid[idx]
            if isinstance(encoded, bytes):
                encoded = encoded.decode("ascii", "ignore")
            else:
                encoded = str(encoded)
            encoded = encoded.rstrip()
            field_label = "field"
            field_units = None
            if encoded.startswith("axis:"):
                field_label = encoded[len("axis:") :]
                if field_label.endswith("]") and "[" in field_label:
                    field_label, field_units = field_label[:-1].rsplit("[", 1)
                if len(field_label) == 0:
                    field_label = "field"
            array = _fortrancore.expdat.data[window].copy()
            wrapped = nddata(array, [array.size], [field_label])
            wrapped.set_axis(field_label, np.asarray(fields[idx], dtype=float))
            wrapped.set_units(field_label, field_units)
            data_id = _fortrancore.expdat.dataid[idx]
            if isinstance(data_id, bytes):
                data_id = data_id.decode("ascii", "ignore")
            else:
                data_id = str(data_id)
            data_id = data_id.rstrip()
            if len(data_id) > 0:
                wrapped.name(data_id)
            traces.append(wrapped)
        if len(traces) == 1:
            return traces[0]
        return traces

    @data.setter
    def data(self, values):
        """Set the experimental spectrum that we are trying to fit.

        Array-like input overwrites the most recently allocated spectrum
        window.  Passing a labelled 1D ``pyspecdata.nddata`` appends a new
        spectrum by taking its axis coordinates from the remaining dimension
        label and storing that label and its units so later outputs can be
        returned as nddata again.

        The overall amplitude match between simulation and experiment is not
        set here.  In the original Fortran flow, ``lfun`` calls ``sscale``
        (or ``sshift`` when shifting is enabled) to determine least-squares
        scale factors and stores them in ``sfac``.  Python exposes that same
        array as :attr:`weights`.  Those coefficients are not constrained to
        sum to unity, and if Fortran's ``iscal`` flag is zero for a site then
        the existing ``sfac`` value is treated as fixed instead of being
        re-optimised.  Python exposes that flag as ``model['iscal']``.  By
        default Fortran initialises ``iscal(i)=1`` for every site, so
        autoscaling is on unless the caller explicitly turns it off.  The
        data setter therefore only records the experimental trace; :meth:`fit`
        is what updates the autoscaled model.
        """

        if _HAS_PYSPECDATA and isinstance(values, nddata):
            if not values.dimlabels:
                raise ValueError(
                    "nddata objects must supply at least one dimension label"
                )
            values.squeeze()
            if len(values.dimlabels) != 1:
                raise ValueError(
                    "nddata input must have exactly one non-singleton"
                    " dimension"
                )
            field_label = values.dimlabels[0]
            field_coords = values.getaxis(field_label)
            if field_coords is None:
                raise ValueError("nddata field axis coordinates are missing")
            fields = np.asarray(field_coords, dtype=float).reshape(-1)
            intensities = np.array(
                [
                    values[field_label, j].item()
                    for j in range(len(fields))
                ],
                dtype=float,
            )
            if len(fields) == 0:
                raise ValueError("nddata field axis contained no points")
            if len(fields) != len(intensities):
                raise ValueError(
                    "field axis and intensity vector lengths do not match"
                )
            if len(fields) > 1:
                steps = np.diff(fields)
                if not np.allclose(steps, steps[0]):
                    raise ValueError(
                        "nddata field axis must be uniformly spaced"
                    )
            if len(fields) > self.max_points:
                raise ValueError(
                    "nddata field axis has "
                    + str(int(len(fields)))
                    + " points; reduce it to "
                    + str(int(self.max_points))
                    + " or fewer to fit the NLSL buffers"
                )
            idx = self.generate_coordinates(
                float(fields[0]),
                float(fields[-1]),
                int(len(fields)),
                label=str(field_label),
            )
            window = self.windows[idx]
            _fortrancore.expdat.data[window] = intensities
            _fortrancore.lmcom.fvec[window] = intensities
            noise = float(np.std(intensities))
            if noise <= float(np.finfo(float).eps):
                noise = 1.0
            _fortrancore.expdat.rmsn[idx] = noise
            data_id = values.name()
            if data_id is None or len(str(data_id).strip()) == 0:
                data_id = "nddata"
            _fortrancore.expdat.dataid[idx] = (
                str(data_id).encode("ascii", "ignore")[:30].ljust(30, b" ")
            )
            if values.get_units(field_label) is None:
                axis_id = "axis:" + str(field_label)
            else:
                axis_id = (
                    "axis:"
                    + f"{field_label}[{values.get_units(field_label)}]"
                )
            _fortrancore.expdat.wndoid[idx] = (
                axis_id.encode("ascii", "ignore")[:20].ljust(20, b" ")
            )
            _fortrancore.expdat.drmode = int(self.derivative_mode)
            _fortrancore.expdat.shftflg = 1 if self.shift else 0
            _fortrancore.expdat.normflg = 0
            _fortrancore.expdat.nrmlz[idx] = 0
            if self.shift:
                _fortrancore.expdat.ishglb = 1
            self.return_nddata = True
            return

        if len(self.windows) == 0:
            raise RuntimeError("no spectra have been allocated")
        flat = np.asarray(values, dtype=float).reshape(-1)
        expected = self.windows[-1].stop - self.windows[-1].start
        if flat.size != expected:
            raise ValueError("intensity vector length mismatch")
        _fortrancore.expdat.data[self.windows[-1]] = flat
        _fortrancore.lmcom.fvec[self.windows[-1]] = flat

    def set_site_weights(self, spectrum_index, weights):
        """Update the scale factors for a specific spectrum index."""

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

        if self.windows[:nspc]:
            min_start = min(window.start for window in self.windows[:nspc])
            max_stop = max(window.stop for window in self.windows[:nspc])
        else:
            min_start = 0
            max_stop = 0

        # Each spectrum occupies a contiguous slice of the shared Fortran work
        # arrays.  ``windows`` are those absolute slices in ``expdat.data`` and
        # ``mspctr.spectr``; ``relative_windows`` remap them onto the trimmed
        # block ``min_start:max_stop`` returned below so Python callers can
        # align multi-spectrum arrays without carrying the full unused storage.

        span = max_stop - min_start
        if span > 0 and nsite > 0:
            trimmed = spectra_src[min_start:max_stop, :nsite]
            site_spectra = trimmed.swapaxes(0, 1)
        else:
            site_spectra = np.empty((nsite, 0), dtype=float)

        if not self.return_nddata:
            site_payload = site_spectra
        else:
            if not _HAS_PYSPECDATA:
                raise ImportError("pyspecdata is required for nddata outputs")
            array = np.asarray(site_spectra, dtype=float)
            if nspc <= 0:
                raise RuntimeError("no field axis metadata is available")
            field_axes = self.field_axes
            field_label = "field"
            field_units = None
            field_coords = None
            for idx in range(nspc):
                coords = np.asarray(field_axes[idx], dtype=float).reshape(-1)
                encoded = _fortrancore.expdat.wndoid[idx]
                if isinstance(encoded, bytes):
                    encoded = encoded.decode("ascii", "ignore")
                else:
                    encoded = str(encoded)
                encoded = encoded.rstrip()
                label = "field"
                units = None
                if encoded.startswith("axis:"):
                    label = encoded[len("axis:") :]
                    if label.endswith("]") and "[" in label:
                        label, units = label[:-1].rsplit("[", 1)
                    if len(label) == 0:
                        label = "field"
                if field_coords is None:
                    field_label = label
                    field_units = units
                    field_coords = coords
                    continue
                if (
                    label != field_label
                    or units != field_units
                    or coords.shape != field_coords.shape
                    or not np.allclose(coords, field_coords)
                ):
                    raise ValueError(
                        "return_nddata requires all active spectra to share one"
                        " common field axis"
                    )

            if nspc <= 1:
                site_payload = nddata(
                    array.copy(),
                    [nsite, field_coords.size],
                    ["sites", field_label],
                )
                site_payload.set_axis("sites", "#")
                site_payload.set_axis(field_label, field_coords)
                site_payload.set_units(field_label, field_units)
            else:
                window_lengths = [
                    window.stop - window.start
                    for window in self.relative_windows[:nspc]
                ]
                if any(length != field_coords.size for length in window_lengths):
                    raise ValueError(
                        "return_nddata requires all active spectra to share one"
                        " common field axis"
                    )
                stacked_site_spectra = np.stack(
                    [
                        array[:, window].copy()
                        for window in self.relative_windows[:nspc]
                    ],
                    axis=0,
                )
                site_payload = nddata(
                    stacked_site_spectra,
                    [
                        stacked_site_spectra.shape[0],
                        stacked_site_spectra.shape[1],
                        stacked_site_spectra.shape[2],
                    ],
                    ["spectrum", "sites", field_label],
                )
                site_payload.set_axis("spectrum", "#")
                site_payload.set_axis("sites", "#")
                site_payload.set_axis(field_label, field_coords)
                site_payload.set_units(field_label, field_units)

        self._last_site_spectra = site_payload
        return self._last_site_spectra

    def _sync_weight_matrix(self):
        """Keep ``/mspctr/sfac/`` aligned with the active site/spectrum counts.

        ``sfac`` is declared in ``mspctr.f90`` as a static ``MXSITE × MXSPC``
        array.  ``nlsinit`` fills every element with ``1.0`` even when no fit
        is running, and the Fortran code never reinitialises the table after
        the initial call.  Changing ``nsite`` or ``nspc`` therefore only
        updates the integer counters; without extra work the exposed
        scale factors would keep whatever values happened to be in memory. This
        helper mirrors the housekeeping performed by the Fortran data loaders
        so Python callers always see a predictable set of coefficients.
        """

        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        new_shape = (nsite, nspc)

        weights = _fortrancore.mspctr.sfac
        if new_shape == self._weight_shape:
            return weights

        if nsite <= 0 or nspc <= 0:
            # No active spectra or sites: blank the entire ``sfac`` matrix so a
            # subsequent resize starts from zero scale factors.
            weights[:, :] = 0.0
            self._weight_shape = new_shape
            return weights

        # Preserve the overlapping block so previously fitted scale factors
        # stay in place when callers expand or shrink the grid of
        # sites/spectra.
        row_stop = min(self._weight_shape[0], nsite)
        col_stop = min(self._weight_shape[1], nspc)
        if row_stop > 0 and col_stop > 0:
            preserved = weights[:row_stop, :col_stop].copy()
        else:
            preserved = None

        # The Fortran initialisation routines seed ``sfac`` with ones, so new
        # rows/columns must do the same.  Reset the full table before restoring
        # any surviving scale factors.
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
