from . import fortrancore as _fortrancore
import os
from pathlib import Path
import numpy as np

_SPECTRAL_PARAMETER_NAMES = {
    "phase",
    "psi",
    "b0",
    "lb",
    "range",
}

from .data import process_spectrum


def _ipfind_wrapper(name: str) -> int:
    """Call the Fortran ``ipfind`` routine if available."""
    token = name.strip().upper()
    lth = len(token)
    if lth == 0:
        raise ValueError("zero-length token!")
    return int(_fortrancore.ipfind(token, lth))


class fit_params(dict):
    """Mapping-like interface for adjusting NLSL fit parameters.

    Keys correspond to the options listed in ``nlshlp.txt`` lines 20–38.
    The values are mirrored directly to the low level ``lmcom`` module so
    that no ``procline`` call is needed.
    """

    def __init__(self):
        super().__init__()
        self._core = _fortrancore
        self._fl_names = [
            n.decode("ascii").strip().lower()
            for n in self._core.lmcom.flmprm_name.tolist()
        ]
        self._il_names = [
            n.decode("ascii").strip().lower()
            for n in self._core.lmcom.ilmprm_name.tolist()
        ]

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


class nlsl(object):
    """Dictionary-like interface to the NLSL parameters."""

    def __init__(self):
        global _fortrancore
        _fortrancore.nlsinit()

        self._fepr_names = [
            name.decode("ascii").strip().lower()
            for name in _fortrancore.eprprm.fepr_name.reshape(-1).tolist()
        ]
        extra_fepr = iter(["fldi", "dfld"])
        for idx, name in enumerate(self._fepr_names):
            if name:
                continue
            try:
                self._fepr_names[idx] = next(extra_fepr)
            except StopIteration:
                break
        self._iepr_names = [
            name.decode("ascii").strip().lower()
            for name in _fortrancore.eprprm.iepr_name.reshape(-1).tolist()
        ]
        extra_iepr = iter(["iwflg", "igflg", "iaflg", "irflg", "jkmn", "jmmn", "ndim"])
        for idx, name in enumerate(self._iepr_names):
            if name:
                continue
            try:
                self._iepr_names[idx] = next(extra_iepr)
            except StopIteration:
                break
        self._fparm = _fortrancore.parcom.fparm
        self._iparm = _fortrancore.parcom.iparm
        self.fit_params = fit_params()
        self._last_layout = None
        self._last_site_spectra = None
        self._weight_shape = (0, 0)

    @property
    def nsites(self) -> int:
        """Number of active sites."""
        return int(_fortrancore.parcom.nsite)

    @nsites.setter
    def nsites(self, value: int) -> None:
        _fortrancore.parcom.nsite = int(value)
        self._resize_weight_matrix()

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
        available through ``model['weights']``.
        """
        ndatot = int(_fortrancore.expdat.ndatot)
        nspc = int(_fortrancore.expdat.nspc)
        if ndatot <= 0 or nspc <= 0:
            raise RuntimeError("no spectra have been evaluated yet")
        _fortrancore.iterat.iter = 1
        _fortrancore.single_point(1)
        return self._capture_state()

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

    # -- mapping protocol -------------------------------------------------

    def __getitem__(self, key):
        key = key.lower()
        if key in ("nsite", "nsites"):
            # ensure the weight matrix tracks external nsite adjustments
            self._resize_weight_matrix()
            return self.nsites
        if key in ("nspc", "nspec", "nspectra"):
            self._resize_weight_matrix()
            return int(_fortrancore.expdat.nspc)
        if key in ("weights", "weight", "sfac"):
            self._resize_weight_matrix()
            nsite = int(_fortrancore.parcom.nsite)
            nspc = int(_fortrancore.expdat.nspc)
            if nsite <= 0 or nspc <= 0:
                if nspc == 1:
                    return np.zeros(nsite, dtype=float)
                return np.zeros((nspc, nsite), dtype=float)
            matrix = _fortrancore.mspctr.sfac
            if nspc == 1:
                return matrix[:nsite, 0].copy()
            return matrix[:nsite, :nspc].swapaxes(0, 1).copy()
        res = _ipfind_wrapper(key)
        if res == 0:
            if key in self._iepr_names:
                idx = self._iepr_names.index(key)
                vals = self._iparm[idx, : self.nsites]
                if np.all(vals == vals[0]):
                    return int(vals[0])
                return vals.copy()
            if key in self._fepr_names:
                idx = self._fepr_names.index(key)
                vals = self._fparm[idx, : self.nsites]
                if np.allclose(vals, vals[0]):
                    return float(vals[0])
                return vals.copy()
            raise KeyError(key)
        if res > 100:
            idx = self._iepr_names.index(key)
            vals = self._iparm[idx, : self.nsites]
        else:
            vals = np.array(
                [_fortrancore.getprm(res, i) for i in range(1, self.nsites + 1)]
            )
        if np.allclose(vals, vals[0]):
            return vals[0]
        return vals

    def __setitem__(self, key, v):
        key = key.lower()
        if key in ("nsite", "nsites"):
            self.nsites = int(v)
            return
        if key in ("nspc", "nspec", "nspectra"):
            _fortrancore.expdat.nspc = int(v)
            self._resize_weight_matrix()
            return
        if key in ("weights", "weight", "sfac"):
            self._resize_weight_matrix()
            nsite = int(_fortrancore.parcom.nsite)
            nspc = int(_fortrancore.expdat.nspc)
            if nsite <= 0 or nspc <= 0:
                raise RuntimeError("weights require positive nsite and nspc")
            matrix = np.asarray(v, dtype=float)
            if matrix.ndim == 1:
                if nspc != 1:
                    raise ValueError("1D weight vector requires a single spectrum")
                if matrix.size < nsite:
                    raise ValueError("insufficient weight values supplied")
                _fortrancore.mspctr.sfac[:, 0] = 0.0
                _fortrancore.mspctr.sfac[:nsite, 0] = matrix[:nsite]
                return
            if matrix.shape[0] < nspc or matrix.shape[1] < nsite:
                raise ValueError("weight matrix shape mismatch")
            _fortrancore.mspctr.sfac[:, :] = 0.0
            _fortrancore.mspctr.sfac[:nsite, :nspc] = matrix[:nspc, :nsite].swapaxes(0, 1)
            return
        if key == "b0":
            self.set_spectral_state(sb0=v)
            return
        if key == "range":
            self.set_spectral_state(srng=v)
            return
        res = _ipfind_wrapper(key)
        iterinput = isinstance(v, (list, tuple, np.ndarray))
        if res == 0:
            if key in self._iepr_names:
                idx = self._iepr_names.index(key)
                if iterinput:
                    limit = min(len(v), self.nsites)
                    self._iparm[idx, :limit] = np.asarray(v[:limit], dtype=int)
                else:
                    self._iparm[idx, : self.nsites] = int(v)
                return
            if key in self._fepr_names:
                idx = self._fepr_names.index(key)
                values = np.asarray(v, dtype=float)
                if values.ndim == 0:
                    self._fparm[idx, : self.nsites] = float(values)
                else:
                    limit = min(values.size, self.nsites)
                    self._fparm[idx, :limit] = values[:limit]
                return
            raise KeyError(key)
        is_spectral = key in _SPECTRAL_PARAMETER_NAMES
        if res > 100:
            if iterinput:
                limit = len(v)
                if not is_spectral:
                    limit = min(limit, self.nsites)
                else:
                    limit = min(limit, int(_fortrancore.expdat.nspc))
                for site_idx in range(limit):
                    _fortrancore.setipr(res, site_idx + 1, int(v[site_idx]))
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
                    _fortrancore.setprm(res, site_idx + 1, float(v[site_idx]))
            else:
                _fortrancore.setprm(res, 0, float(v))

    def __contains__(self, key):
        key = key.lower()
        if key in ("nsite", "nsites"):
            return True
        if key in self._fepr_names or key in self._iepr_names:
            return True
        return _ipfind_wrapper(key) != 0

    def canonical_name(self, name: str) -> str:
        """Return the canonical parameter name for *name*.

        Uses the Fortran ``ipfind`` routine to resolve aliases.  If *name*
        is already canonical it is returned unchanged.  ``KeyError`` is raised
        when the name cannot be resolved.
        """
        key = name.lower()
        if key in ("nsite", "nsites"):
            return "nsite"
        if key in self._fepr_names or key in self._iepr_names:
            return key
        res = _ipfind_wrapper(key)
        if res == 0:
            raise KeyError(name)
        if res > 100:
            return self._iepr_names[res - 101]
        if res > 0:
            return self._fepr_names[res - 1]
        if res > -100:
            idx = -res - 1
        else:
            idx = -res - 101
        return self._fepr_names[idx]

    def __iter__(self):
        return iter(self.keys())

    @property
    def layout(self):
        """Metadata describing the most recent spectral evaluation."""
        if self._last_layout is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_layout

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
        if ix0 + points > mxpt:
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
        core.expdat.wndoid[idx] = (
            window_label.encode("ascii", "ignore").ljust(20, b" ")
        )

        core.expdat.ndatot = ix0 + points

        self._resize_weight_matrix()

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

    def set_spectral_state(
        self,
        *,
        sb0=None,
        srng=None,
        ishift=None,
        shift=None,
        normalize_flags=None,
    ):
        """Update the active spectrum metadata without reloading data files."""

        expdat = _fortrancore.expdat

        def _broadcast_to_spectra(row_index, values):
            count = int(_fortrancore.parcom.nsite)
            if count <= 0:
                count = 1
            limit = min(count, self._fparm.shape[1])
            if limit <= 0:
                return
            if values.size >= limit:
                self._fparm[row_index, :limit] = values[:limit]
            else:
                self._fparm[row_index, :limit] = values[0]
        if sb0 is not None:
            values = np.atleast_1d(np.asarray(sb0, dtype=float))
            expdat.sb0[: values.size] = values
            if "b0" in self._fepr_names:
                _broadcast_to_spectra(self._fepr_names.index("b0"), values)
        if srng is not None:
            values = np.atleast_1d(np.asarray(srng, dtype=float))
            expdat.srng[: values.size] = values
        if (sb0 is not None or srng is not None) and "range" in self._fepr_names:
            start_row = self._fepr_names.index("range") + 1
            step_row = start_row + 1
            if start_row < self._fparm.shape[0]:
                base = float(expdat.sb0[0]) if expdat.sb0.shape[0] > 0 else 0.0
                span = float(expdat.srng[0]) if expdat.srng.shape[0] > 0 else 0.0
                start_val = base - 0.5 * span
                _broadcast_to_spectra(start_row, np.array([start_val], dtype=float))
            if step_row < self._fparm.shape[0]:
                if expdat.npts.shape[0] > 0 and expdat.npts[0] > 1:
                    step_val = float(expdat.srng[0]) / float(expdat.npts[0] - 1)
                else:
                    step_val = 0.0
                _broadcast_to_spectra(step_row, np.array([step_val], dtype=float))
        if ishift is not None:
            values = np.atleast_1d(np.asarray(ishift, dtype=int))
            expdat.ishft[: values.size] = values
        if shift is not None:
            values = np.atleast_1d(np.asarray(shift, dtype=float))
            expdat.shft[: values.size] = values
        if normalize_flags is not None:
            values = np.atleast_1d(np.asarray(normalize_flags, dtype=int))
            expdat.nrmlz[: values.size] = values
        self._last_site_spectra = None

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

        self._last_layout = {
            "ixsp": _fortrancore.expdat.ixsp[:nspc] - 1,
            "npts": _fortrancore.expdat.npts[:nspc].copy(),
            "sbi": _fortrancore.expdat.sbi[:nspc].copy(),
            "sdb": _fortrancore.expdat.sdb[:nspc].copy(),
            "ndatot": ndatot,
            "nsite": nsite,
            "nspc": nspc,
        }

        if ndatot > 0 and nsite > 0:
            site_spectra = spectra_src[:ndatot, :nsite].swapaxes(0, 1)
        else:
            site_spectra = np.empty((nsite, ndatot), dtype=float)

        self._last_site_spectra = site_spectra
        return self._last_site_spectra

    def _resize_weight_matrix(self):
        # Keep the site-by-spectrum weight table sized to the active model.
        # The Levenberg–Marquardt core works with a rectangular ``sfac`` array
        # whose leading dimensions are controlled by ``nsite`` and ``nspc``.
        # Whenever either axis changes we zero-fill the new storage and copy
        # the overlapping block of any previously stored populations so fits
        # retain their converged weights after grid changes.
        nsite = int(_fortrancore.parcom.nsite)
        nspc = int(_fortrancore.expdat.nspc)
        new_shape = (nsite, nspc)
        if new_shape == self._weight_shape:
            return

        weights = _fortrancore.mspctr.sfac
        if nsite <= 0 or nspc <= 0:
            weights[:, :] = 0.0
            self._weight_shape = new_shape
            return

        preserved = None
        if self._weight_shape[0] > 0 and self._weight_shape[1] > 0:
            row_stop = min(self._weight_shape[0], nsite)
            col_stop = min(self._weight_shape[1], nspc)
            if row_stop > 0 and col_stop > 0:
                preserved = weights[:row_stop, :col_stop].copy()

        weights[:, :] = 1.0
        if preserved is not None:
            row_stop, col_stop = preserved.shape
            weights[:row_stop, :col_stop] = preserved

        self._weight_shape = new_shape

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
