from . import fortrancore as _fortrancore

# now I want to use "decorated names" to reassign things

def nlsinit():
    "Initializes all NLSL parameters, and readies a run"
    _fortrancore.nlsinit()

def procline(val):
    "process a line of a traditional format text NLSL runfile"
    _fortrancore.procline(val)

class _parameter_class(object):
    def __init__(self):
        self._asdict = {}

    @property
    def asdict(self):  # this is the getter
        """A dictionary containing the various floating-point ESR parameters"""
        try:
            names = _fortrancore.eprprm.fepr_name.T.reshape(-1, _fortrancore.eprprm.nfprm).view(dtype=f"|S{_fortrancore.eprprm.nfprm}")[:, 0].tolist()
            names = [x[: x.find(" ")] for x in names]
            values = _fortrancore.eprprm.fepr
            self._asdict = dict(zip(names, values))
        except AttributeError:
            # if low level arrays are unavailable, return cached values
            pass
        return self._asdict
    @asdict.setter
    def asdict(self,value):
        """Update parameters using procline commands."""
        if not isinstance(value, dict):
            raise TypeError("value must be a dictionary")
        for key, val in value.items():
            if isinstance(val, (list, tuple)):
                val_str = ', '.join(str(v) for v in val)
            else:
                val_str = str(val)
            _fortrancore.procline(f"let {key} = {val_str}")
            self._asdict[key] = val

parameters = _parameter_class()

__all__ = [x for x in dir() if x[0] != '_']
#print "all is",__all__# on import, this generates a reasonable result, proving there is no craziness, ad also that the init function above is run
