from . import fortrancore as _fortrancore

# now I want to use "decorated names" to reassign things

def nlsinit():
    "Initializes all NLSL parameters, and readies a run"
    _fortrancore.nlsinit()

def procline(val):
    "process a line of a traditional format text NLSL runfile"
    _fortrancore.procline(val)

class parameter_class (object):
    def __init__(self):
        _x = 0
    @property
    def param1(self):
        pass
        return
    @param1.setter
    def param1(self,value):
        print "you tried to set param1 to",value
    @param1.getter
    def param1(self):
        "This is param 1 it does stuff"
        print "you tried to get param1"

parameters = parameter_class()

__all__ = filter(lambda x: x[0] != '_',dir())
#print "all is",__all__# on import, this generates a reasonable result, proving there is no craziness, ad also that the init function above is run
