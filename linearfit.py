#!/usr/bin/python
# "nicely" plot the output of image_extractval.py

import sys, os
import numpy as np
import matplotlib.pyplot as plt

datafile = sys.argv[1]

def f(B, x):
      return B[0]*x + B[1]

# the errors are not estimatable...
def linear_fit_fast(x, y, yerr=None):
    from scipy.optimize import leastsq
    if yerr == None: yerr = np.ones(len(y))
    for i,e in enumerate(yerr):
        if e == 0: yerr[i] = 1
    def errf(B, x, y, yerr):
        return f(B,x)-y/yerr
    out = leastsq(errf, [-1. ,0.], args=(x, y, yerr), full_output=1)
    return (out[0][0], out[0][1], 0, 0)

data = np.loadtxt(datafile, comments='#', dtype=np.dtype({'names':['freq','flux','rms'], 'formats':[float,float,float]}))

data['rms'] = 5*0.434*data['rms']/data['flux']
B = linear_fit_fast(np.log10(data['freq']), np.log10(data['flux'], np.log10(data['rms'])))
print "Regression:", B

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.set_xlim(xmin=100e6,xmax=2e9)
ax.set_ylim(ymin=0.01,ymax=3)
ax.set_xlabel(r'Freq [Hz]')
ax.set_ylabel(r'Flux density [Jy]')
ax.set_yscale('log')
ax.set_xscale('log')
ax.errorbar(data['freq'], data['flux'], yerr=data['rms'], fmt='ko')
freqs = np.logspace(6, 10, num=100)
ax.plot(freqs, [10**f(B, np.log10(freq)) for freq in freqs], label='y={:.2}*x+{:2.2}'.format(B[0],B[1]))
ax.legend(loc=1)
fig.savefig(datafile.replace('.dat','.png'))
