from . import fortrancore as _fortrancore


class fit_params(dict):
    """Mapping-like interface for adjusting NLSL fit parameters.

    Keys correspond to the options listed in ``nlshlp.txt`` lines 20–38.
    Setting an item issues the appropriate ``fit`` command via
    :func:`fortrancore.procline`.
    """

    def __init__(self):
        super().__init__()
        self._core = _fortrancore

    def __setitem__(self, key, value):
        super().__setitem__(key.lower(), value)

    def __getitem__(self, key):
        return super().get(key.lower())



class nlsl(object):
    """Dictionary-like interface to the NLSL parameters."""

    def __init__(self):
        self._params = {}
        # automatically initialize parameters when the object is created
        _fortrancore.nlsinit()
        self.fit_params = fit_params()

    def procline(self, val):
        """Process a line of a traditional format text NLSL runfile."""
        _fortrancore.procline(val)

    def fit(self):
        """Run the nonlinear least-squares fit using current parameters."""
        opts = []
        for k, v in self.fit_params.items():
            opts.append(f"{k} {v}")
        if opts:
            _fortrancore.procline("fit " + " ".join(opts))
        else:
            _fortrancore.fitl()

    # -- utility -----------------------------------------------------------
    def _refresh(self):
        """Refresh the cached parameter dictionary from the Fortran core."""
        try:
            names = (
                _fortrancore.eprprm.fepr_name.T.reshape(
                    -1, _fortrancore.eprprm.nfprm
                )
                .view(dtype=f"|S{_fortrancore.eprprm.nfprm}")[:, 0]
                .tolist()
            )
            names = [x[: x.find(" ")] for x in names]
            values = _fortrancore.eprprm.fepr
            self._params = dict(zip(names, values))
        except AttributeError:
            # if low level arrays are unavailable, return cached values
            pass

    # -- mapping protocol -------------------------------------------------
    def __getitem__(self, key):
        self._refresh()
        return self._params[key]

    def __setitem__(self, key, value):
        if isinstance(value, (list, tuple)):
            val_str = ", ".join(str(v) for v in value)
        else:
            val_str = str(value)
        _fortrancore.procline(f"let {key} = {val_str}")
        self._params[key] = value

    def __contains__(self, key):
        self._refresh()
        return key in self._params

    def __iter__(self):
        self._refresh()
        return iter(self._params)

    def keys(self):
        self._refresh()
        return self._params.keys()

    def items(self):
        self._refresh()
        return self._params.items()

    def values(self):
        self._refresh()
        return self._params.values()

    def get(self, key, default=None):
        self._refresh()
        return self._params.get(key, default)

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

