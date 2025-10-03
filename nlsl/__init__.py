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
        return [(k, self[k]) for k in self.keys() if len(k)>0]

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

        nspc = int(_fortrancore.expdat.nspc)
        if nspc >= nser:
            nspc = 0
            _fortrancore.expdat.ndatot = 0

        nspc += 1
        if nspc > mxspc:
            raise ValueError("Maximum number of spectra exceeded")
        _fortrancore.expdat.nspc = nspc
        idx = nspc - 1

        ix0 = int(_fortrancore.expdat.ndatot)
        n_points = int(spectrum.y.size)
        if ix0 + n_points > mxpt:
            raise ValueError("Insufficient storage for spectrum")

        _fortrancore.expdat.ishft[idx] = 1 if shift else 0
        _fortrancore.expdat.idrv[idx] = mode
        _fortrancore.expdat.nrmlz[idx] = 1 if normalize_active else 0
        _fortrancore.expdat.npts[idx] = n_points
        _fortrancore.expdat.ixsp[idx] = ix0 + 1
        _fortrancore.expdat.sbi[idx] = spectrum.start
        _fortrancore.expdat.sdb[idx] = spectrum.step
        _fortrancore.expdat.srng[idx] = spectrum.step * max(n_points - 1, 0)
        _fortrancore.expdat.shft[idx] = 0.0
        _fortrancore.expdat.tmpshft[idx] = 0.0
        _fortrancore.expdat.slb[idx] = 0.0
        _fortrancore.expdat.sb0[idx] = 0.0
        _fortrancore.expdat.sphs[idx] = 0.0
        _fortrancore.expdat.spsi[idx] = 0.0

        eps = float(np.finfo(float).eps)
        _fortrancore.expdat.rmsn[idx] = (
            spectrum.noise if spectrum.noise > eps else 1.0
        )

        _fortrancore.expdat.iform[idx] = 0
        _fortrancore.expdat.ibase[idx] = int(bc_points)
        power = 1
        while power < n_points:
            power *= 2
        _fortrancore.expdat.nft[idx] = power

        _fortrancore.expdat.data[ix0 : ix0 + n_points] = spectrum.y
        _fortrancore.expdat.ndatot = ix0 + n_points

        _fortrancore.lmcom.fvec[ix0 : ix0 + n_points] = spectrum.y
        if hasattr(_fortrancore.mspctr, "spectr"):
            _fortrancore.mspctr.spectr[ix0 : ix0 + n_points, :] = 0.0
        if hasattr(_fortrancore.mspctr, "wspec"):
            _fortrancore.mspctr.wspec[ix0 : ix0 + n_points, :] = 0.0
        if hasattr(_fortrancore.mspctr, "sfac"):
            _fortrancore.mspctr.sfac[:, idx] = 1.0

        label = base_name.encode("ascii", "ignore")[:30]
        _fortrancore.expdat.dataid[idx] = label.ljust(30, b" ")

        trimmed_name = base_name.strip()
        window_label = f"{nspc:2d}: {trimmed_name}"[:19] + "\0"
        _fortrancore.expdat.wndoid[idx] = (
            window_label.encode("ascii", "ignore").ljust(20, b" ")
        )

        if shift:
            _fortrancore.expdat.ishglb = 1
        _fortrancore.expdat.shftflg = 1 if shift else 0
        _fortrancore.expdat.normflg = 1 if normalize_active else 0
        _fortrancore.expdat.bcmode = int(bc_points)
        _fortrancore.expdat.drmode = mode
        _fortrancore.expdat.nspline = requested_points
        _fortrancore.expdat.inform = 0

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
            vals = self._iparm[idx, :self.nsites]
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

    def keys(self):
        return list(self._fepr_names) + list(self._iepr_names)

    def items(self):
        return [(k, self[k]) for k in self.keys() if len(k)>0]

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
