#!/usr/bin/python
# "nicely" plot the output of image_extractval.py

import sys, os
import numpy as np
import matplotlib.pyplot as plt

datafile = sys.argv[1]

def f(x, B0, B1):
    return B0*x + B1
# extimate errors and accept errors on ydata
def linear_fit(x, y, yerr=None):
    from scipy.optimize import curve_fit
    if yerr == None: yerr = np.ones(len(y))
    for i,e in enumerate(yerr):
        if e == 0: yerr[i] = 1
    out = curve_fit(f, x, y, [-1. ,0.], yerr)
    # return B0, B1, errB0, errB1 (err are in std dev)
    return (out[0][0], out[0][1], np.sqrt(out[1][0][0]), np.sqrt(out[1][1][1]))

data = np.loadtxt(datafile, comments='#', dtype=np.dtype({'names':['freq','flux','rms'], 'formats':[float,float,float]}))

data_nolog = {'flux':np.copy(data['flux']),'freq':np.copy(data['freq']),'rms':np.copy(data['rms'])}
data['rms'] = 0.434*data['rms']/data['flux']
data['flux'] = np.log10(data['flux'])
data['freq'] = np.log10(data['freq'])
B = linear_fit(data['freq'], data['flux'], yerr=data['rms'])
print "Regression:", B

# in linear plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.set_xlabel(r'Log Freq [Hz]')
ax.set_ylabel(r'Log Flux density [Jy]')
ax.set_xlim(xmin=min(data['freq'])-min(data['freq'])/10.,\
            xmax=max(data['freq'])+max(data['freq'])/10.)
ax.set_ylim(ymin=min(data['flux'])-min(data['flux'])/10.,\
            ymax=max(data['flux'])+max(data['flux'])/10.)
ax.errorbar(data['freq'], data['flux'], yerr=data['rms'], fmt='ko')
ax.errorbar(data['freq'], data['flux'], fmt='k-')
ax.plot(np.arange(6,10,0.1), [f(freq, B[0], B[1]) for freq in np.arange(6,10,0.1)], \
        label=r'y={:.1f}$\pm${:.1f}*x+{:2.0f}$\pm${:2.0f}'.format(B[0],B[2],B[1],B[3]))
ax.legend(loc=1)
#fig.savefig(datafile.replace('.dat','-lin.pdf'))
del fig

# in log plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'Freq [Hz]')
ax.set_ylabel(r'Flux density [Jy]')
ax.set_xlim(xmin=min(data_nolog['freq'])-min(data_nolog['freq'])/10.,\
            xmax=max(data_nolog['freq'])+max(data_nolog['freq'])/10.)
ax.set_ylim(ymin=min(data_nolog['flux'])-min(data_nolog['flux'])/10.,\
            ymax=max(data_nolog['flux'])+max(data_nolog['flux'])/10.)
# workaround for too big errors inlog plot
ymaxerr = data_nolog['rms']
yminerr = data_nolog['rms']
yminerr[ data_nolog['rms'] >= data_nolog['flux'] ] = \
        data_nolog['flux'][ data_nolog['rms'] >= data_nolog['flux'] ]*.9999
ax.errorbar(data_nolog['freq'], data_nolog['flux'], yerr=[ymaxerr,yminerr], fmt='k-')
ax.errorbar(data_nolog['freq'], data_nolog['flux'], fmt='ko')
freqs = np.logspace(6, 10, num=100)
ax.plot(freqs, [10**f(np.log10(freq), B[0], B[1]) for freq in freqs], \
        label=r'y={:.1f}$\pm${:.1f}*x+{:2.0f}$\pm${:2.0f}'.format(B[0],B[2],B[1],B[3]))
ax.legend(loc=1)
fig.savefig(datafile.replace('.dat','.pdf'))
del fig


