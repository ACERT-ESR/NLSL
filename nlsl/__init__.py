from . import fortrancore as _fortrancore
import importlib
import numpy as np


def _ipfind_wrapper(name: str) -> int:
    """Call the Fortran ``ipfind`` routine if available."""
    upper = name.upper()
    return int(_fortrancore.ipfind(upper, len(upper)))


class fit_params(dict):
    """Mapping-like interface for adjusting NLSL fit parameters.

    Keys correspond to the options listed in ``nlshlp.txt`` lines 20–38.
    The values are mirrored directly to the low level ``lmcom`` module so
    that no ``procline`` call is needed.
    """

    def __init__(self):
        super().__init__()
        self._core = _fortrancore
        self._fl_names = [n.decode('ascii').strip().lower() for n in self._core.lmcom.flmprm_name.tolist()]
        self._il_names = [n.decode('ascii').strip().lower() for n in self._core.lmcom.ilmprm_name.tolist()]

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
        return [(k, self[k]) for k in self.keys()]

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
        # initialize the Fortran core so the parameter name arrays are set
        _fortrancore.nlsinit()

        self._fepr_names = [
            name.decode('ascii').strip().lower()
            for name in _fortrancore.eprprm.fepr_name.reshape(-1).tolist()
        ]
        self._iepr_names = [
            name.decode('ascii').strip().lower()
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

    # -- mapping protocol -------------------------------------------------

    def __getitem__(self, key):
        key = key.lower()
        if key in ("nsite", "nsites"):
            return self.nsites
        res = _ipfind_wrapper(key)
        if res == 0:
            raise KeyError(key)
        if res > 100:
            ixp = abs(res) % 100
            if ixp == _fortrancore.eprprm.INFLD:
                vals = _fortrancore.expdat.npts[: self.nsites]
            elif ixp == _fortrancore.eprprm.IIDERV:
                vals = _fortrancore.expdat.idrv[: self.nsites]
            else:
                vals = self._iparm[ixp - 1, : self.nsites]
            if np.all(vals == vals[0]):
                return int(vals[0])
            return np.array(vals, dtype=int)
        vals = np.array([
            _fortrancore.getprm(res, i) for i in range(1, self.nsites + 1)
        ])
        if np.allclose(vals, vals[0]):
            return vals[0]
        return vals

    def __setitem__(self, key, value):
        key = key.lower()
        if key in ("nsite", "nsites"):
            self.nsites = int(value)
            return
        res = _ipfind_wrapper(key)
        if res == 0:
            raise KeyError(key)
        if isinstance(value, (list, tuple, np.ndarray)):
            for i, v in enumerate(value, start=1):
                if i > self.nsites:
                    break
                if res > 100:
                    _fortrancore.setipr(res, i, int(v))
                else:
                    _fortrancore.setprm(res, i, float(v))
        else:
            if res > 100:
                _fortrancore.setipr(res, 0, int(value))
            else:
                _fortrancore.setprm(res, 0, float(value))

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
        return [(k, self[k]) for k in self.keys()]

    def values(self):
        return [self[k] for k in self.keys()]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def update(self, other):
        """Update multiple parameters at once."""
        if isinstance(other, dict):
            items = other.items()
        else:
            items = other
        for k, v in items:
            self[k] = v


# expose the class for creating additional instances
NLSL = nlsl

__all__ = [x for x in dir() if x[0] != "_"]
