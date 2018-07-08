#!/usr/bin/env python
""" Provides acces to the weights of an MS file. """

import math
import os
import sys

import numpy as np

from casacore.tables import taql
from matplotlib import cm
from matplotlib.pyplot import figure, show
from matplotlib import ticker

__author__ = 'Frits Sweijen'
__credits__= 'Francesco de Gasperin'
__version__ = '1.0.0'
__maintainer__ = 'Frits Sweijen'
__email__ = 'sweijen <at> strw.leidenuniv.nl'
__status__ = 'Development'

def normalize(x, nmin, nmax):
    normed = (x - nmin) / (nmax - nmin)
    return normed

def plot_weight_channel_time_average(msfile, pol=0, delta=16, per_antenna=False, threshold=1e5):
    ''' Plot the weights in the WEIGHT_SPECTRUM column of an MS file as function of frequency.

    Args:
        msfile (str): MS file to use.
        pol (int): integer value of 0, 1, 2 or 3 (as applicable) denoting which polarization to use.
        delta (int): the number of channels to use when calculating the variance.
        threshold (int, float): upper limit for plotting weights to deal with outliers. Selects only weights below this value. This does not affect calculations.
    Returns:
        None
    '''
    # Polarization indices are 0, 1, 2, 3 = XX, YY, XY, YX, respectively.
    print 'Plotting weights vs. channels for %s' % (msfile,)
    # Select only rows where not all data is flagged.
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM, FLAG FROM $msfile WHERE ALL(FLAG)==False && ALL(ISNAN(DATA))==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT ANTENNA, MEANS(GAGGR(DATA), 0) AS DATA_REAL, MEANS(GAGGR(WEIGHT_SPECTRUM),0), FLAG AS WEIGHT FROM $t1 GROUPBY TIME')
    #t = taql('SELECT TIME, ANTENNA, BOXEDMEAN(GAGGR(REAL(DATA)), 10, 1, 1) AS DATAR FROM $t1 WHERE !ANY(ISNAN(REAL(DATA))) GROUPBY TIME')
    # Get the polarization setup; X/Y or R/L.
    temp = taql('SELECT CORR_TYPE from '+msfile+'/POLARIZATION')
    if temp.getcol('CORR_TYPE')[0] in np.asarray([5, 6, 7, 8]):
        # Circular polarization.
        polarization = ['RR', 'LL', 'RL', 'LR']
    elif temp.getcol('CORR_TYPE')[0] in np.asarray([9, 10, 11, 12]):
        polarization = ['XX', 'YY', 'XY', 'YX']
    # Average the data over all time stamps.
    flags = t.getcol('FLAG')
    datar = t.getcol('DATA_REAL')
    datar = np.ma.MaskedArray(data=datar, mask=flags)
    datar = datar[:,:,pol]
    data_taverage = datar.mean(axis=0)
    print datar.shape
    weights = t.getcol('WEIGHT')
    weights = np.ma.MaskedArray(data=weights, mask=flags)
    weights = weights[:, :, pol]
    # Obtain channel frequencies in Hz.
    chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    # Select the first table, column CHAN_FREQ and convert to MHz.
    freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    #print 'Frequencies [MHz]: '
    #print freq
    # Calculate the variance in the visibilities over channels.
    print 'Calculating visibility variance.'
    #variance = np.ones(shape=(weights.shape[0], weights.shape[1]//delta, weights.shape[2]))
    variance = np.zeros(shape=(data_taverage.shape[0]//delta))
    # Subtract adjacent channels to eliminate physical signal.
    data_shifted = np.roll(data_taverage, -1, axis=0)
    datas = data_taverage - data_shifted
    print datas.shape
    datar = datas.real
    datai = datas.imag
    for i in xrange(datas.shape[0]//delta):
        # Take a frequency bin of delta channels.
        vr = np.nanvar(datar[delta*i: delta*i+delta], axis=0)
        vi = np.nanvar(datai[delta*i: delta*i+delta], axis=0)
        v = (vr + vi) / 2.
        print np.any(np.isfinite(v))
        if not np.any(np.isfinite(v)):
            variance[i] = np.nan
        else:
            variance[i] = np.where(np.isfinite(1. / v), 1. / v, 0)
    print variance
    print variance.shape
    # Plot the results.
    print 'Plotting weights...'
    if not per_antenna:
        weights = np.mean(weights, axis=0, keepdims=True)
    imgname ='weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+str(delta)+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    weights = normalize(weights, np.min(weights), np.max(weights))

    for w in weights:
        if len(freq) > weights.shape[1]:
            ax.scatter(freq[:-1], w, color='k', marker='.')
        elif len(freq) == weights.shape[1]:
            ax.scatter(freq, w, color='k', marker='.')
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.min(weights), np.max(weights))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Normalized Weights')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    
    print 'Plotting variance weights...'
    imgname ='var_chan_tavg_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+str(delta)+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    #variance = variance[:, :, pol]
    #variance = normalize(variance, np.nanmin(variance), np.nanmax(variance))
    #f = freq[indices]
    f = np.zeros(data_taverage.shape[0]//delta)
    for i in xrange(data_taverage.shape[0]//delta):
        f[i] = np.mean(freq[delta*i:delta*i + delta])

    ax.plot(f, variance, 'd--', color='k')
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.nanmin(variance), np.nanmax(variance))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Variance Normalized w.r.t. '+polarization[0])
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

def plot_weight_channel(msfile, pol=0, delta=16, per_antenna=False, threshold=1e5):
    ''' Plot the weights in the WEIGHT_SPECTRUM column of an MS file as function of frequency.

    Args:
        msfile (str): MS file to use.
        pol (int): integer value of 0, 1, 2 or 3 (as applicable) denoting which polarization to use.
        delta (int): the number of channels to use when calculating the variance.
        threshold (int, float): upper limit for plotting weights to deal with outliers. Selects only weights below this value. This does not affect calculations.
    Returns:
        None
    '''
    # Polarization indices are 0, 1, 2, 3 = XX, YY, XY, YX, respectively.
    print 'Plotting weights vs. channels for %s' % (msfile,)
    # Select only rows where not all data is flagged.
    t1 = taql('SELECT TIME, ANTENNA1 AS ANTENNA, DATA, WEIGHT_SPECTRUM, FLAG FROM $msfile WHERE ALL(FLAG)==False && ALL(ISNAN(DATA))==False')
    # Average over all baselines (axis 0) and group the resulting data by antenna.
    t = taql('SELECT ANTENNA, MEANS(GAGGR(DATA), 0) AS DATA_REAL, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, FLAG FROM $t1 GROUPBY ANTENNA')
    #t = taql('SELECT TIME, ANTENNA, BOXEDMEAN(GAGGR(REAL(DATA)), 10, 1, 1) AS DATAR FROM $t1 WHERE !ANY(ISNAN(REAL(DATA))) GROUPBY TIME')
    # Get the polarization setup; X/Y or R/L.
    temp = taql('SELECT CORR_TYPE from '+msfile+'/POLARIZATION')
    if temp.getcol('CORR_TYPE')[0] in np.asarray([5, 6, 7, 8]):
        # Circular polarization.
        polarization = ['RR', 'LL', 'RL', 'LR']
    elif temp.getcol('CORR_TYPE')[0] in np.asarray([9, 10, 11, 12]):
        polarization = ['XX', 'YY', 'XY', 'YX']
    flags = t.getcol('FLAG')
    datar = t.getcol('DATA_REAL')
    datar = np.ma.MaskedArray(data=datar, mask=flags)
    datar = datar[:,:,pol]
    print datar.shape
    weights = t.getcol('WEIGHT')
    weights = np.ma.MaskedArray(data=weights, mask=flags)
    weights = weights[:, :, pol]

    antennas = t.getcol('ANTENNA')
    print len(antennas), ' antennas'
    antenna_names = taql('SELECT NAME FROM '+msfile+'/ANTENNA')
    antenna_names = antenna_names.getcol('NAME')
    print antenna_names
    # Obtain channel frequencies in Hz.
    chan_freq = taql('SELECT CHAN_FREQ FROM '+msfile+'/SPECTRAL_WINDOW')
    # Select the first table, column CHAN_FREQ and convert to MHz.
    freq = chan_freq[0]['CHAN_FREQ'] * 1e-6
    #print 'Frequencies [MHz]: '
    #print freq
    # Calculate the variance in the visibilities over channels.
    print 'Calculating visibility variance.'
    #variance = np.ones(shape=(weights.shape[0], weights.shape[1]//delta, weights.shape[2]))
    variance = np.ones(shape=(weights.shape[0], weights.shape[1]//delta))
    # Subtract adjacent channels to eliminate physical signal.
    datar_shifted = np.roll(datar, -1, axis=1)
    datars = datar - datar_shifted
    datar = datars.real
    datai = datars.imag
    if not per_antenna:
        datar = np.mean(datar, axis=0, keepdims=True)
        datai = np.mean(datai, axis=0, keepdims=True)
        antenna_names = [{'NAME':''}]
    for i in xrange(datar.shape[1]//delta):
        # Take a frequency bin of delta channels.
        vr = np.nanvar(datar[:,delta*i: delta*i+delta], axis=1)
        vi = np.nanvar(datai[:,delta*i: delta*i+delta], axis=1)
        v = (vr + vi) / 2.
        print np.any(np.isfinite(v))
        if not np.any(np.isfinite(v)):
            variance[:,i] = np.nan
        else:
            variance[:,i] = np.where(np.isfinite(1. / v), 1. / v, 0)
    
    # Plot the results.
    print 'Plotting weights...'
    if not per_antenna:
        weights = np.mean(weights, axis=0, keepdims=True)
    imgname ='weight_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+str(delta)+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    weights = normalize(weights, np.min(weights), np.max(weights))

    for w, c, a in zip(weights, colors, antennas):
        if len(freq) > weights.shape[1]:
            ax.scatter(freq[:-1], w, color=c, marker='.', label=antenna_names[a]['NAME'])
        elif len(freq) == weights.shape[1]:
            ax.scatter(freq, w, color=c, marker='.', label=antenna_names[a]['NAME'])
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.min(weights), np.max(weights))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Normalized Weights')
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    
    print 'Plotting variance weights...'
    imgname ='var_chan_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+str(delta)+'.png'
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    #variance = variance[:, :, pol]
    variance = normalize(variance, np.nanmin(variance), np.nanmax(variance))

    #f = freq[indices]
    f = freq[::delta]
    f = np.zeros(shape=(weights.shape[1]//delta))
    for i in xrange(weights.shape[1]//delta):
        f[i] = np.mean(freq[delta*i:delta*i + delta])
    for v, c, a in zip(variance, colors, antennas):
        print antenna_names[a]['NAME']
        ax.plot(f, v, '--.', color=c, label=antenna_names[a]['NAME'])        
        '''
        if len(f) > variance.shape[1]:
            #ax.scatter(f[:-1], v, color=c, marker='.', label=antenna_names[a]['NAME'])
            ax.plot(f[:-1], v, '--.', color=c, label=antenna_names[a]['NAME'])
        elif len(f) == variance.shape[1]:
            #ax.scatter(f, v, color=c, marker='.', label=antenna_names[a]['NAME'])
            ax.plot(f, v, '--.', color=c, label=antenna_names[a]['NAME'])
        '''
    ax.set_xlim(min(freq), max(freq))
    ax.set_ylim(np.nanmin(variance), np.nanmax(variance))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    major_formatter = ticker.FuncFormatter(lambda x, pos: '%.2f'%(x,))
    ax.xaxis.set_major_formatter(major_formatter)
    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Variance Normalized w.r.t. '+polarization[0])
    leg = ax.legend(bbox_to_anchor=(1.05, 1.0), ncol=3, borderaxespad=0.0)
    
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    return 

def plot_weight_time(msfile, delta=10, plot_time_unit='h'):
    ''' Plot the weights in the WEIGHT_SPECTRUM column of an MS file as function of time for each polarization.

    Args:
        msfile (str): MS file to use.
        delta (int): bin width in timestamps to use for calculating the variance in the visibilities.
        plot_time_unit (str): either 'h' (implemented) or 's' (not implemented) to use hours or seconds on the time axis.
    Returns:
        None
    '''
    print 'Plotting weights vs. time for %s' % (msfile,)
    imgname ='weight_time_' + msfile[msfile.find('SB'):msfile.find('SB')+5]+'_'+str(delta)+'.png'
    # Select the time, weights and elevation of ANTENNA1 averaging over baselines/antennas.
    # Select only unflagged data.
    t1 = taql('select ANTENNA1, ANTENNA2, DATA, WEIGHT_SPECTRUM, TIME, FIELD_ID, FLAG from $msfile where ALL(FLAG)==False')
    # Get the polarization setup; X/Y or R/L.
    temp = taql('SELECT CORR_TYPE from '+msfile+'/POLARIZATION')
    if temp.getcol('CORR_TYPE')[0] in np.asarray([5, 6, 7, 8]):
        # Circular polarization.
        polarization = ['RR', 'LL', 'RL', 'LR']
    elif temp.getcol('CORR_TYPE')[0] in np.asarray([9, 10, 11, 12]):
        polarization = ['XX', 'YY', 'XY', 'YX']
    # Select time, weights and elevation.
    # This gets the average elevation w.r.t. to antenna 1, where the average is taken over all baselines.
    t = taql('SELECT TIME, MEANS(GAGGR(WEIGHT_SPECTRUM),0) AS WEIGHT, MEANS(GAGGR(MSCAL.AZEL1()[1]), 0) AS ELEV FROM $t1 GROUPBY TIME')
    t2 = taql('SELECT TIME, MEANS(GAGGR(DATA),0) AS DATA_REAL FROM $t1 GROUPBY TIME')
    time = t.getcol('TIME')
    elevation = t.getcol('ELEV')

    flags = t2.getcol('FLAG')
    datar = t2.getcol('DATA_REAL')
    datar = np.ma.MaskedArray(data=datar, mask=flags)
    weights = t.getcol('WEIGHT')
    weights = np.ma.MaskedArray(data=weights, mask=flags)

    # Calculate variance over 10 timestamp intervals. The last interval may be shorter.
    # datar shape is (timestamps, channels, polarizations)
    variance = np.ones(shape=(len(time)//delta, weights.shape[1], weights.shape[2]))
    datar_shifted = np.roll(datar, -1, axis=1)
    datars = datar - datar_shifted
    datar = datars.real
    datai = datars.imag
    for i in xrange(len(time)//delta):
        vr = np.nanvar(datar[delta*i: delta*i+delta,:,:], axis=0)
        vi = np.nanvar(datai[delta*i: delta*i+delta,:,:], axis=0)
        v = (vr + vi) / 2.
        if np.any(np.isfinite(v)):
            variance[i] = np.where(np.isfinite(1. / v), 1. / v, 0)
        else:
            variance[i] = np.nan

    # Select the weights for all timestamps one channel and one polarization (in that order).
    weights = t.getcol('WEIGHT')
    # Plot weights for elevation and save the image.
    print 'Plotting...'
    colors = iter(cm.rainbow(np.linspace(0, 1, len(weights))))
    fig = figure()
    fig.suptitle(msfile, fontweight='bold')
    ax = fig.add_subplot(111)

    if plot_time_unit.lower() == 'h':
        time_h = time / 3600
        time_h = (time_h - math.floor(time_h[0]/24) * 24)
        
        # Deal with outlier weights if necessary.        
        wexponents = np.floor(np.log10(weights)).astype(int)
        print np.unique(wexponents)
        weights = np.where(wexponents < 0, weights, 1e-10)
        
        # Normalize the weights w.r.t. the XX/RR polarization.
        nmin = np.nanmin(weights[:, :, 0])
        nmax = np.nanmax(weights[:, :, 0])
        
        ax.scatter(time_h, normalize(weights[:, 5, 0], nmin, nmax), marker='.', color='C1', alpha=0.25, label=polarization[0]+' Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 1], nmin, nmax), marker='.', color='C2', alpha=0.25, label=polarization[1]+' Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 2], nmin, nmax), marker='.', color='C3', alpha=0.25, label=polarization[2]+' Weights')
        ax.scatter(time_h, normalize(weights[:, 5, 3], nmin, nmax), marker='.', color='C4', alpha=0.25, label=polarization[3]+' Weights')
        del nmin, nmax
        
        # Normalize the statistic w.r.t. the XX/RR polarization.
        print variance[:, :, 0]
        nmin = np.nanmin(variance[:, :, 0])
        nmax = np.nanmax(variance[:, :, 0])
        print nmin, nmax
        indices = ((np.asarray(range(0, len(variance[:, 5, 0]))) + 0.5) * delta).astype(int)
        ax.plot(time_h[indices], normalize(variance[:, 5, 0], nmin, nmax), '--d', color='C1', label=polarization[0]+' Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 1], nmin, nmax), '--d', color='C2', label=polarization[1]+' Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 2], nmin, nmax), '--d', color='C3', label=polarization[2]+' Boxed variance $\\Delta=%d$'%(delta,))
        ax.plot(time_h[indices], normalize(variance[:, 5, 3], nmin, nmax), '--d', color='C4', label=polarization[3]+' Boxed variance $\\Delta=%d$'%(delta,))
        # Plot the elevation as a function of time.
        ax_elev = ax.twinx()
        ax_elev.plot(time_h, elevation * 180/np.pi, 'k', linewidth=2, label='Elevation')

        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%.2d:%.2d:%.2d' % (int(x), (x%1)*60, (((x%1)*60)%1 * 60))))
        ax.set_xlim(min(time_h), max(time_h))
        ax.set_ylim(0, 1.5)
        ax.set_xlabel('Time [h]')
    elif plot_time_unit.lower() == 's':
        # To do.
        pass
    ax.set_ylabel('Weights Normalized w.r.t. '+polarization[0])
    ax_elev.set_ylabel('Elevation [deg]')

    # Deal with the legend.
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax_elev.get_legend_handles_labels()
    leg = ax.legend(handles+handles2, labels+labels2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, borderaxespad=0.0)
    print 'Saving plot as %s' % (imgname,)
    fig.savefig(imgname, bbox_inches='tight', additional_artists=leg, dpi=250)
    #show()
    return

if __name__ == '__main__':
    # Get the MS filename.
    msfile = sys.argv[1]
    #plot_weight_channel_time_average(msfile, delta=32)
    plot_weight_channel(msfile, delta=32)
    #plot_weight_time(msfile, delta=100)
