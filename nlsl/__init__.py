from . import fortrancore as _fortrancore
import os
from pathlib import Path
import numpy as np

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
        self._iepr_names = [
            name.decode("ascii").strip().lower()
            for name in _fortrancore.eprprm.iepr_name.reshape(-1).tolist()
        ]
        self._fparm = _fortrancore.parcom.fparm
        self._iparm = _fortrancore.parcom.iparm
        self.fit_params = fit_params()
        self._last_layout = None
        self._last_site_spectra = None
        self._last_weights = None

    @property
    def nsites(self) -> int:
        """Number of active sites."""
        return int(_fortrancore.parcom.nsite)

    @nsites.setter
    def nsites(self, value: int) -> None:
        _fortrancore.parcom.nsite = int(value)

    def procline(self, val):
        """Process a line of a traditional format text NLSL runfile."""
        _fortrancore.procline(val)

    def fit(self):
        """Run the nonlinear least-squares fit using current parameters."""
        _fortrancore.fitl()
        return self._capture_state()

    @property
    def current_spectrum(self):
        """Evaluate the current spectral model without running a full fit."""
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
            return self.nsites
        res = _ipfind_wrapper(key)
        if res == 0:
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
        res = _ipfind_wrapper(key)
        iterinput = isinstance(v, (list, tuple, np.ndarray))
        if res == 0:
            raise KeyError(key)
        if res > 100:
            if iterinput:
                for site_idx in range(len(v)):
                    _fortrancore.setipr(res, site_idx + 1, int(v[site_idx]))
            else:
                for site_idx in range(self.nsites):
                    _fortrancore.setipr(res, site_idx + 1, int(v))
        else:
            if iterinput:
                for site_idx in range(len(v)):
                    _fortrancore.setprm(res, site_idx, float(v[site_idx]))
            else:
                for site_idx in range(self.nsites):
                    _fortrancore.setprm(res, site_idx, float(v))

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

    @property
    def weights(self):
        """Return the most recently evaluated site weights."""
        if self._last_weights is None:
            raise RuntimeError("no spectra have been evaluated yet")
        return self._last_weights

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

        return idx, data_slice

    def _capture_state(self):
        nspc = int(_fortrancore.expdat.nspc)
        ndatot = int(_fortrancore.expdat.ndatot)
        nsite = int(_fortrancore.parcom.nsite)

        spectra_src = _fortrancore.mspctr.spectr
        weights_src = _fortrancore.mspctr.sfac

        nspc = min(
            nspc,
            _fortrancore.expdat.ixsp.shape[0],
            _fortrancore.expdat.npts.shape[0],
            _fortrancore.expdat.sbi.shape[0],
            _fortrancore.expdat.sdb.shape[0],
            weights_src.shape[1],
        )
        nsite = min(nsite, spectra_src.shape[1], weights_src.shape[0])
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

        if nspc > 0 and nsite > 0:
            weight_matrix = weights_src[:nsite, :nspc].swapaxes(0, 1)
        else:
            weight_matrix = np.empty((nspc, nsite), dtype=float)

        self._last_site_spectra = site_spectra
        self._last_weights = weight_matrix

        return self._last_site_spectra, self._last_weights

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
