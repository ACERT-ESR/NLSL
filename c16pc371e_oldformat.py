#!/usr/bin/python
from pylab import *
from subprocess import Popen, PIPE, STDOUT
import os
def read_column_data(filename):
    fp = open(filename,'r')
    data = []
    for j in fp.readlines():
        data.append(j.split())
    data = array(data,dtype = double)
    fp.close()
    return data
os.system('make') # check that everything is up to date
#cmd = './nlsl < c16pc371e.run'
#print "about to run '"+cmd+"'"
#os.system(cmd) # actually run nlsl
print("about to run nlsl")
#proc = Popen(['nlsl'],stdout = PIPE, stdin = PIPE, stderr = STDOUT)
if os.name == 'posix':
    proc = Popen(['./nlsl'],stdin = PIPE, stderr = STDOUT)
else:
    proc = Popen(['nlsl'],stdin = PIPE, stderr = STDOUT)
fp = open('c16pc371e.run')
output = proc.communicate(input = fp.read())
fp.close()
print("output was:",output)
data = read_column_data('c16pc371e.spc')
fields = data[:,0]
experimental = data[:,1]
fit = data[:,2]
integral_of_spectrum = cumsum(experimental)
normalization = abs(sum(integral_of_spectrum))
fig = figure(figsize = (9,6))
fig.add_axes([0.1,0.1,0.6,0.8]) # l b w h
if data.shape[1] > 3:
    components = data[:,3:]
else:
    components = None
plot(fields,experimental/normalization,'k',linewidth = 1,label = 'experimental')
plot(fields,fit/normalization,'k',alpha = 0.5,linewidth = 2,label = 'fit')
max_of_fit = max(fit)/normalization
if components is not None:
    plot(fields,components/normalization,alpha = 0.3,linewidth = 1,label = 'component')
ax = gca()
ylims = ax.get_ylim()
scale_integral = max_of_fit/max(integral_of_spectrum)
plot(fields,integral_of_spectrum * scale_integral,'k:',alpha = 0.5,linewidth = 1,label = '$\int \int$ (scaled by %0.2g)'%scale_integral)
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
