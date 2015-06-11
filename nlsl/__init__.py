from . import fortrancore as _fortrancore

# now I want to use "decorated names" to reassign things

def nlsinit():
    "Initializes all NLSL parameters, and readies a run"
    _fortrancore.nlsinit()

def procline(val):
    "process a line of a traditional format text NLSL runfile"
    _fortrancore.procline(val)

class _parameter_class (object):
    def __init__(self):
        _x = 0
    @property
    def asdict(self):#this is the getter
        "A dictionary containing the various floating-point ESR parameters"
        self._listofparms = _fortrancore.eprprm.fepr_name.T.reshape(-1,10).view(dtype="|S10")[:,0].tolist()
        self._listofparms = [x[0:x.find(' ')] for x in self._listofparms]
        self._asdict = dict(zip(self._listofparms,_fortrancore.eprprm.fepr))
        # this is not the full array yet, needs to be adjusted
        return self._asdict
    @asdict.setter
    def asdict(self,value):
        print "not supported to set from here yet"

parameters = _parameter_class()

__all__ = filter(lambda x: x[0] != '_',dir())
#print "all is",__all__# on import, this generates a reasonable result, proving there is no craziness, ad also that the init function above is run
