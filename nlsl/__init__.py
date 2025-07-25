from . import fortrancore as _fortrancore
import re
from pathlib import Path


def _parse_eprprm_names():
    """Fallback parser for fepr_name and iepr_name from the Fortran source."""
    fe_names = []
    ie_names = []
    src = Path(__file__).resolve().parent.parent / 'src' / 'fortran' / 'eprprm.f90'
    if not src.exists():
        return fe_names, ie_names

    dimsrc = Path(__file__).resolve().parent.parent / 'src' / 'fortran' / 'nlsdim.f90'
    dims = {}
    if dimsrc.exists():
        dims = {name: int(val) for name, val in re.findall(r"(NFPRM|NIPRM)=([0-9]+)", dimsrc.read_text())}

    text = src.read_text()
    # map constant names to indices
    consts = {
        name: int(val)
        for name, val in re.findall(r"(\w+)=([0-9]+)", text)
    }

    fe_names = [''] * dims.get('NFPRM', 0)
    ie_names = [''] * dims.get('NIPRM', 0)

    for m in re.finditer(r"fepr_name\((\w+)\)\s*=\s*\"(\w+)\"", text):
        idx = consts.get(m.group(1))
        if idx:
            fe_names[idx - 1] = m.group(2).lower()

    for m in re.finditer(r"iepr_name\((\w+)\)\s*=\s*\"(\w+)\"", text):
        idx = consts.get(m.group(1))
        if idx:
            ie_names[idx - 1] = m.group(2).lower()

    # remove empty placeholders
    fe_names = [n for n in fe_names if n]
    ie_names = [n for n in ie_names if n]
    return fe_names, ie_names


class fit_params(dict):
    """Mapping-like interface for adjusting NLSL fit parameters.

    Keys correspond to the options listed in ``nlshlp.txt`` lines 20–38.
    The values are mirrored directly to the low level ``lmcom`` module so
    that no ``procline`` call is needed.
    """

    def __init__(self):
        super().__init__()
        self._core = _fortrancore
        self._fl_names = [n.strip().lower() for n in self._core.lmcom.flmprm_name.tolist()]
        self._il_names = [n.strip().lower() for n in self._core.lmcom.ilmprm_name.tolist()]

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

        try:
            self._fepr_names = [
                name.strip().lower()
                for name in _fortrancore.eprprm.fepr_name.reshape(-1).tolist()
            ]
            self._iepr_names = [
                name.strip().lower()
                for name in _fortrancore.eprprm.iepr_name.reshape(-1).tolist()
            ]
        except AttributeError:
            # allow import of the high level API even if low level
            # modules were not compiled successfully
            self._fepr_names, self._iepr_names = _parse_eprprm_names()
            self._fepr = _fortrancore.parcom.fparm[:, 0]
            self._iepr = _fortrancore.parcom.iparm[:, 0]
        else:
            self._fepr = _fortrancore.eprprm.fepr
            self._iepr = _fortrancore.eprprm.iepr

        self.fit_params = fit_params()

    def procline(self, val):
        """Process a line of a traditional format text NLSL runfile."""
        _fortrancore.procline(val)

    def fit(self):
        """Run the nonlinear least-squares fit using current parameters."""
        _fortrancore.fitl()

    # -- mapping protocol -------------------------------------------------
    def __getitem__(self, key):
        key = key.lower()
        if key in self._fepr_names:
            return self._fepr[self._fepr_names.index(key)]
        if key in self._iepr_names:
            return self._iepr[self._iepr_names.index(key)]
        raise KeyError(key)

    def __setitem__(self, key, value):
        key = key.lower()
        if key in self._fepr_names:
            idx = self._fepr_names.index(key)
            if isinstance(value, (list, tuple)):
                for i, v in enumerate(value):
                    self._fepr[idx + i] = v
            else:
                self._fepr[idx] = value
        elif key in self._iepr_names:
            idx = self._iepr_names.index(key)
            if isinstance(value, (list, tuple)):
                for i, v in enumerate(value):
                    self._iepr[idx + i] = v
            else:
                self._iepr[idx] = value
        else:
            raise KeyError(key)

    def __contains__(self, key):
        key = key.lower()
        return key in self._fepr_names or key in self._iepr_names

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

