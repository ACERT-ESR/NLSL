#!/usr/bin/python
from pylab import *
import os
def read_column_data(filename):
    fp = open(filename,'r')
    data = []
    for j in fp.readlines():
        data.append(j.split())
    data = array(data,dtype = double)
    return data
os.system('make') # check that everything is up to date
#os.system('./nlsl < c16pc371e.run') # actually run nlsl
data = read_column_data('c16pc371e.spc')
fields = data[:,0]
experimental = data[:,1]
fit = data[:,2]
normalization = sum(cumsum(experimental))
normalization = 1
fig = figure(figsize = (9,6))
fig.add_axes([0.1,0.1,0.6,0.8]) # l b w h
if data.shape[1] > 3:
    components = data[:,3:]
plot(fields,experimental/normalization,'k',linewidth = 1,label = 'experimental')
plot(fields,fit/normalization,'k',alpha = 0.5,linewidth = 2,label = 'fit')
plot(fields,components/normalization,alpha = 0.3,linewidth = 1,label = 'component')
ax = gca()
ylims = ax.get_ylim()
plot(fields,cumsum(experimental)/normalization,'k:',alpha = 0.5,linewidth = 1,label = '$\int \int$')
legend(bbox_to_anchor=(1.05,0,0.5,1), # bounding box l b w h
        loc = 2, # upper left (of the bounding box)
        borderaxespad=0.)
ax.set_ylim(ylims)
rms = mean((fit/normalization-experimental/normalization)**2)
ax.text(0.75, 0.75, 'rms = %0.2g'%rms,
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes)
show()
