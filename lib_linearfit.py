#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# "nicely" plot the output of image_extractval.py

import sys, os
import numpy as np
import matplotlib.pyplot as plt


def f(x, B0, B1):
    return B0*x + B1


def twopoint_spidx_bootstrap(freq, flux, flux_err, niter=10000):
    """
    Quick bootstrap for spectral index calulcation
    freq: 2 array
    flux: 2 or 2xN array
    flux_err: 2 or 2xN array
    N is the number of sources
    """
    # calculate spidx assuming [iter,source,freq_point] shapes
    def spidx(freq, flux):
        return np.log10(flux[:,:,0]/flux[:,:,1])/np.log10(freq[:,:,0]/freq[:,:,1])

    freq = np.array(freq).astype(float)
    flux = np.array(flux).astype(float)
    flux_err = np.array(flux_err).astype(float)
    # if only 1 source, add degenerate axis
    if flux.shape == (2,): flux = np.expand_dims(flux, axis=1)
    if flux_err.shape == (2,): flux_err = np.expand_dims(flux_err, axis=1)
    flux = flux.T
    flux_err = flux_err.T
    nsource = flux.shape[0]

    results = np.zeros(shape=(niter,nsource))
    random_flux = np.resize(flux, (niter, nsource, 2)) + np.resize(flux_err, (niter, nsource, 2)) * np.random.randn(niter, nsource, 2)
    random_flux[random_flux <= 0] = np.nan # remove negative, this create a bias
    freq = np.resize(freq, (niter, nsource, 2))
    results = spidx(freq, random_flux)

    mean = np.nanmean(results,axis=0)
    err = np.nanstd(results,axis=0)
    return mean, err


def linear_fit_bootstrap(x, y, yerr, niter=1000, tolog=False):
    # An issue arises with scipy.curve_fit when errors in the y data points
    # are given.  Only the relative errors are used as weights, so the fit
    # parameter errors, determined from the covariance do not depended on the
    # magnitude of the errors in the individual data points.  This is clearly wrong. 
    # 
    # To circumvent this problem I have implemented a simple bootstraping 
    # routine that uses some Monte-Carlo to determine the errors in the fit
    # parameters.  This routines generates random datay points starting from
    # the given datay plus a random variation. 
    #
    # The random variation is determined from average standard deviation of y
    # points in the case where no errors in the y data points are avaiable.
    #
    # If errors in the y data points are available, then the random variation 
    # in each point is determined from its given error. 
    # 
    # A large number of random data sets are produced, each one of the is fitted
    # an in the end the variance of the large number of fit results is used as 
    # the error for the fit parameters. 

    # Estimate the confidence interval of the fitted parameter using
    # the bootstrap Monte-Carlo method
    # http://phe.rockefeller.edu/LogletLab/whitepaper/node17.html
    #
    # tolog : convert in log space x, y, and yerr before doing linear regression
    # use: (a, b, sa, sb) = linear_fit_bootstrap(x, y, yerr)

    from scipy import optimize
    errfunc = lambda B, x, y: f(x, B[0], B[1]) - y

    x = np.array(x)
    y = np.array(y)
    if yerr is not None: yerr = np.array(yerr)

    if tolog:
        if yerr is not None: yerr = 0.434*yerr/y
        x=np.log10(x)
        y=np.log10(y)

    pfit, pcov, infodict, errmsg, success = optimize.leastsq( errfunc, [-1, 0], args=(x, y), full_output=1)

    # 2 vals without error, cannot estimate sigmas
    if len(y) == 2 and yerr is None: return (pfit[0], pfit[1], 0, 0)

    residuals = errfunc( pfit, x, y )
    s_res = np.std(residuals)
    ps = []
    # TODO: remove cycle and use only array shapes
    # n random data sets are generated and fitted
    for i in range(niter):
      if yerr is None:
          randomDelta = np.random.normal(0., s_res, len(y))
          randomdataY = y + randomDelta
      else:
          randomDelta = np.array( [ np.random.normal(0., derr, 1)[0] for derr in yerr ] ) 
          randomdataY = y + randomDelta
      randomfit, randomcov = optimize.leastsq( errfunc, [-1, 0], args=(x, randomdataY), full_output=0)
      ps.append( randomfit ) 

    ps = np.array(ps)
    mean_pfit = np.mean(ps,0)
    Nsigma = 1. # 1sigma gets approximately the same as methods above
                # 1sigma corresponds to 68.3% confidence interval
                # 2sigma corresponds to 95.44% confidence interval
    err_pfit = Nsigma * np.std(ps,0) 

    return (mean_pfit[0], mean_pfit[1], err_pfit[0], err_pfit[1])


# extimate errors and accept errors on ydata
def linear_fit(x, y, yerr=None, tolog=False):
    # Using OLS (X|Y)" # for more algo read: Isobe et al 1990
    # tolog : convert in log space x, y, and yerr before doing linear regression
    from scipy.optimize import curve_fit
    if tolog:
        if not yerr is None: yerr = 0.434*yerr/y
        x=np.log10(x)
        y=np.log10(y)
    #if yerr is None: yerr = np.ones(len(y))
    #for i,e in enumerate(yerr):
    #    if e == 0: yerr[i] = 1
    out = curve_fit(f, x, y, [-1. ,0.], yerr)
    # return B0, B1, errB0, errB1 (err are in std dev)
    if type(out[1]) is np.ndarray:
        return (out[0][0], out[0][1], np.sqrt(out[1][0][0]), np.sqrt(out[1][1][1]))
    else:
        return (out[0][0], out[0][1], 0, 0)


# extimate errors and accept errors on x and y-data
def linear_fit_odr(x, y, xerr=None, yerr=None):
    from scipy import odr
    def f(B, x):
        return B[0]*x + B[1]
    linear = odr.Model(f)
    if xerr == None: xerr = np.ones(len(x))
    if yerr == None: yerr = np.ones(len(y))
    for i,e in enumerate(yerr):
       if e == 0: yerr[i] = 1
    mydata = odr.RealData(x, y, sx=xerr, sy=yerr)
    myodr = odr.ODR(mydata, linear, beta0=[-1., 0.])
    myoutput = myodr.run()
    return(myoutput.beta[0],myoutput.beta[1],myoutput.sd_beta[0],myoutput.sd_beta[1])


def armonizeXY(dataX, dataY, errY):
    """
    Return xmin,xmax,ymin,ymax in order to have the two axis
    covering almost the same amount of orders of magnitudes
    input must be the log10 of data!!!
    """

    minY = min(dataY) - abs(errY[np.where(dataY == min(dataY))])[0]
    maxY = max(dataY) + abs(errY[np.where(dataY == max(dataY))])[0]
    minX = min(dataX)
    maxX = max(dataX)

    diffX = maxX - minX # X separation
    diffY = maxY - minY # Y separation (including errors)
    maxdiff = max(diffX, diffY)*1.1
    xmin = np.floor(((minX+diffX/2.) - maxdiff/2.)*10.)/10.
    xmax = np.ceil(((minX+diffX/2.) + maxdiff/2.)*10.)/10.
    ymin = np.floor(((minY+diffY/2.) - maxdiff/2.)*10.)/10.
    ymax = np.ceil(((minY+diffY/2.) + maxdiff/2.)*10.)/10.
    return xmin, xmax, ymin, ymax
    

def plotlinax(data, plotname):
    """Plot spectra using linear axes
    data are a dict: {flux:[],freq:[],rms:[]}
    """
    #reorder following freq
    srtidx = np.argsort(data['freq'])
    data = {'flux':data['flux'][srtidx], 'freq':data['freq'][srtidx], 'rms':data['rms'][srtidx]}
    
    # take the log10
    thisdata = {'flux': np.log10(data['flux']), 'freq': np.log10(data['freq']), 'rms': 0.434*data['rms']/data['flux']}

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(wspace=0)
    ax = fig.add_subplot(111)
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.label_outer()
    ax.set_xlabel(r'Log Frequency [Hz]')
    ax.set_ylabel(r'Log Flux density [Jy]')
    xmin, xmax, ymin, ymax = armonizeXY(thisdata['freq'], thisdata['flux'], thisdata['rms'])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.errorbar(thisdata['freq'], thisdata['flux'], yerr=thisdata['rms'], fmt='ko')
#    ax.errorbar(thisdata['freq'], thisdata['flux'], fmt='k-')
    B = linear_fit(thisdata['freq'], thisdata['flux'], yerr=thisdata['rms'])
#    B = linear_fit_odr(thisdata['freq'], thisdata['flux'], yerr=thisdata['rms'])
    print("Regression:", B)
    ax.plot(freqs, [f(freq, B[0], B[1]) for freq in freqs], \
        label=r'$\alpha$={:.2f}$\pm${:.2f}'.format(B[0],B[2]))
    ax.legend(loc=1)
    print("Writing "+plotname)
    fig.savefig(plotname, bbox_inches='tight')
    del fig

def plotlogax(data, plotname):
    """Plot spectra using log axes
    data are a dict: {flux:[],freq:[],rms:[]}
    """
    #reorder following freq
    srtidx = np.argsort(data['freq'])
    data = {'flux':data['flux'][srtidx], 'freq':data['freq'][srtidx], 'rms':data['rms'][srtidx]}

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(wspace=0)
    ax = fig.add_subplot(111)
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.label_outer()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'Frequency [Hz]')
    ax.set_ylabel(r'Flux density [Jy]')
    xmin, xmax, ymin, ymax = armonizeXY(np.log10(data['freq']), np.log10(data['flux']), 0.434*data['rms']/data['flux'])
    ax.set_xlim(10**xmin, 10**xmax)
    ax.set_ylim(10**ymin, 10**ymax)
    # workaround for too big errors inlog plot
    ymaxerr = data['rms']
    yminerr = data['rms']
    # i.e. if it would fall in the negative part of the plot
    yminerr[ data['rms'] >= data['flux'] ] = \
        data['flux'][ data['rms'] >= data['flux'] ]*.9999 # let it be just ~0
    ax.errorbar(data['freq'], data['flux'], yerr=[ymaxerr,yminerr], fmt='ko')
#    ax.errorbar(data['freq'], data['flux'], fmt='k-')
    freqs = np.logspace(np.log10(min(data['freq'])), np.log10(max(data['freq'])), num=100)
    B = linear_fit(np.log10(data['freq']), np.log10(data['flux']),\
#    B = linear_fit_odr(np.log10(data['freq']), np.log10(data['flux']),\
        yerr = 0.434*data['rms']/data['flux'])
    print("Regression:", B)
    ax.plot(freqs, [10**f(np.log10(freq), B[0], B[1]) for freq in freqs], \
        label=r'$\alpha$={:.2f}$\pm${:.2f}'.format(B[0],B[2]))
    ax.legend(loc=1)

    # minor tics
    #ax.xaxis.set_minor_formatter(plt.LogFormatter(base=10.0, labelOnlyBase=False))
    #ax.xaxis.set_major_formatter(plt.LogFormatter(base=10.0, labelOnlyBase=False))
    #ax.yaxis.set_minor_formatter(plt.FormatStrFormatter('%.2f'))
    #ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    count = 0
    for i in ax.xaxis.get_minorticklabels():
        if (count%4 == 0):
            i.set_fontsize(12)
        else:
            i.set_visible(False)
        count+=1
    for i in ax.yaxis.get_minorticklabels():
        if (count%4 == 0):
            i.set_fontsize(12)
        else:
            i.set_visible(False)
        count+=1

    print("Writing "+plotname)
    fig.savefig(plotname, bbox_inches='tight')
    del fig

if __name__ == "__main__":
    import optparse
    opt = optparse.OptionParser(usage="%prog -d datafile", version="%prog 0.1")
    opt.add_option('-d', '--datafile', help='Input data file with freq, flux and rms', default=None)
    opt.add_option('-o', '--output', help='Name of the output plot [default = datafile.pdf]', default=None)
    opt.add_option('-l', help='Output plot shows the log10 of the values', action="store_true", dest="log")
    options, _null = opt.parse_args()
    datafile = options.datafile
    if datafile == None: sys.exit('missing data file')
    print("Data file = "+datafile)
    output = options.output
    if output == None: output = datafile+'.pdf'
    print("Output file = "+output)
    log = options.log

    data = np.loadtxt(datafile, comments='#', dtype=np.dtype({'names':['freq','flux','rms'], 'formats':[float,float,float]}))

    if log:
        plotlinax(data, output)
    else:
        plotlogax(data, output)



