#!/usr/bin/python
# "nicely" plot the output of image_extractval.py

import sys, os
import numpy as np
import matplotlib.pyplot as plt

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

def plotlinax(data, plotname):
    """Plot spectra using linear axes
    data are a dict: {flux:[],freq:[],rms:[]}
    """
    #reorder following freq
    srtidx = np.argsort(data['freq'])
    data = {'flux':data['flux'][srtidx], 'freq':data['freq'][srtidx], 'rms':data['rms'][srtidx]}
    thisdata = {'flux': np.log10(data['flux']), 'freq': np.log10(data['freq']), 'rms': 0.434*data['rms']/data['flux']}
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'Log Freq [Hz]')
    ax.set_ylabel(r'Log Flux density [Jy]')
    ax.set_xlim(xmin=min(thisdata['freq'])-abs(min(thisdata['freq'])/100.),\
            xmax=max(thisdata['freq'])+abs(max(thisdata['freq'])/100.))
    ax.set_ylim(ymin=min(thisdata['flux'])-abs(min(thisdata['flux'])/100.),\
            ymax=max(thisdata['flux'])+abs(max(thisdata['flux'])/100.))
    ax.errorbar(thisdata['freq'], thisdata['flux'], yerr=thisdata['rms'], fmt='ko')
    ax.errorbar(thisdata['freq'], thisdata['flux'], fmt='k-')
    B = linear_fit(thisdata['freq'], thisdata['flux'], yerr=thisdata['rms'])
    print "Regression:", B
    ax.plot(np.arange(6,10,0.1), [f(freq, B[0], B[1]) for freq in np.arange(6,10,0.1)], \
        label=r'y={:.1f}$\pm${:.1f}*x+{:2.0f}$\pm${:2.0f}'.format(B[0],B[2],B[1],B[3]))
    ax.legend(loc=1)
    fig.savefig(plotname)
    del fig

def plotlogax(data, plotname):
    """Plot spectra using log axes
    data are a dict: {flux:[],freq:[],rms:[]}
    """
    #reorder following freq
    srtidx = np.argsort(data['freq'])
    data = {'flux':data['flux'][srtidx], 'freq':data['freq'][srtidx], 'rms':data['rms'][srtidx]}
    print data
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'Freq [Hz]')
    ax.set_ylabel(r'Flux density [Jy]')
    ax.set_xlim(xmin=min(data['freq'])-min(data['freq'])/10.,\
            xmax=max(data['freq'])+max(data['freq'])/10.)
    ax.set_ylim(ymin=min(data['flux'])-min(data['flux'])/10.,\
            ymax=max(data['flux'])+max(data['flux'])/10.)
    # workaround for too big errors inlog plot
    ymaxerr = data['rms']
    yminerr = data['rms']
    yminerr[ data['rms'] >= data['flux'] ] = \
        data['flux'][ data['rms'] >= data['flux'] ]*.9999
    ax.errorbar(data['freq'], data['flux'], yerr=[ymaxerr,yminerr], fmt='k-')
    ax.errorbar(data['freq'], data['flux'], fmt='ko')
    freqs = np.logspace(6, 10, num=100)
    B = linear_fit(np.log10(data['freq']), np.log10(data['flux']),\
        yerr = 0.434*data['rms']/data['flux'])
    print "Regression:", B
    ax.plot(freqs, [10**f(np.log10(freq), B[0], B[1]) for freq in freqs], \
        label=r'y={:.1f}$\pm${:.1f}*x+{:2.0f}$\pm${:2.0f}'.format(B[0],B[2],B[1],B[3]))
    ax.legend(loc=1)
    fig.savefig(plotname)
    del fig

if __name__ == "__main__":
    import optparse
    opt = optparse.OptionParser(usage="%prog images", version="%prog 0.1")
    opt.add_option('-d', '--datafile', help='Input data file with freq, flux and rms', default=None)
    opt.add_option('-o', '--output', help='Name of the output plot [default = datafile.pdf]', default=None)
    opt.add_option('-l', help='Output plot shows the log10 of the values', action="store_true", dest="log")
    options, _null = opt.parse_args()
    datafile = options.datafile
    if datafile == None: sys.exit('missing data file')
    print "Data file = "+datafile
    output = options.output
    if output == None: output = datafile.replace('.dat','.pdf')
    print "Output file = "+output
    log = options.log

    data = np.loadtxt(datafile, comments='#', dtype=np.dtype({'names':['freq','flux','rms'], 'formats':[float,float,float]}))

    if log:
        plotlinax(data, output)
    else:
        plotlogax(data, output)
