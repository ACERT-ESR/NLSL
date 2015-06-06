from .fortrancore import nlsinit
from .fortrancore import procline

# now I want to use "decorated names" to reassign things

__all__ = filter(lambda x: x[0] != '_',dir())
#print "all is",__all__# on import, this generates a reasonable result, proving there is no craziness, ad also that the init function above is run
