#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 - Francesco de Gasperin, Henrik Edler
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

# Usage: BLsmooth.py vis.MS [--options]
# Load a MS, smooth visibilities according to the baseline length,
# i.e. shorter BLs are averaged more, and write a new column to the MS

import os, sys
import optparse
import logging
import numpy as np
from scipy.ndimage import gaussian_filter1d as smooth_real1d

import casacore.tables as pt

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')
logging.info('BL-based smoother - Francesco de Gasperin, Henrik Edler')


def smooth_imag1d(values, std, axis):
    """ smooth complex values along one axis. """
    valuesR = smooth_real1d(values.real, std, axis=axis, truncate=3)
    valuesI = smooth_real1d(values.imag, std, axis=axis, truncate=3)
    return valuesR + 1j * valuesI


def addcol(ms, incol, outcol, setvals=True):
    """ Add a new column to a MS. """
    if outcol not in ms.colnames():
        logging.info('Adding column: '+outcol)
        coldmi = ms.getdminfo(incol)
        coldmi['NAME'] = outcol
        ms.addcols(pt.makecoldesc(outcol, ms.getcoldesc(incol)), coldmi)
    if setvals and (outcol != incol):
        # copy columns val
        logging.info('Set '+outcol+'='+incol)
        pt.taql("UPDATE $ms SET "+outcol+"="+incol)


opt = optparse.OptionParser(usage="%prog [options] MS", version="%prog 3.0")
opt.add_option('-f', '--ionfactor', help='Gives an indication on how strong is the ionosphere [default: 0.01]', type='float', default=0.01)
opt.add_option('-s', '--bscalefactor', help='Gives an indication on how the smoothing varies with BL-lenght [default: 0.8]', type='float', default=0.8)
opt.add_option('-i', '--incol', help='Column name to smooth [default: DATA]', type='string', default='DATA')
opt.add_option('-o', '--outcol', help='Output column [default: SMOOTHED_DATA]', type="string", default='SMOOTHED_DATA')
opt.add_option('-w', '--weight', help='Save the newly computed WEIGHT_SPECTRUM, this action permanently modify the MS! [default: False]', action="store_true", default=False)
opt.add_option('-r', '--restore', help='If WEIGHT_SPECTRUM_ORIG exists then restore it before smoothing [default: False]', action="store_true", default=False)
opt.add_option('-b', '--nobackup', help='Do not backup the old WEIGHT_SPECTRUM in WEIGHT_SPECTRUM_ORIG [default: do backup if -w]', action="store_true", default=False)
opt.add_option('-a', '--onlyamp', help='Smooth only amplitudes [default: smooth real/imag]', action="store_true", default=False)
opt.add_option('-t', '--notime', help='Do not do smoothing in time [default: False]', action="store_true", default=False)
opt.add_option('-q', '--nofreq', help='Do not do smoothing in frequency [default: False]', action="store_true", default=False)
opt.add_option('-c', '--chunks', help='Split the I/O in n chunks. If you run out of memory, set this to a value > 2.', default=2, type='int')
(options, msfile) = opt.parse_args()

if msfile == []:
    opt.print_help()
    sys.exit(0)
msfile = msfile[0]
if not os.path.exists(msfile):
    logging.error("Cannot find MS file.")
    sys.exit(1)
# open input/output MS
ms = pt.table(msfile, readonly=False, ack=False)

with pt.table(msfile + '::SPECTRAL_WINDOW', ack=False) as freqtab:
    freq = freqtab.getcol('REF_FREQUENCY')[0]
    freqpersample = np.mean(freqtab.getcol('RESOLUTION'))
    timepersample = ms.getcell('INTERVAL',0)

# get info on all baselines
with pt.taql("SELECT ANTENNA1,ANTENNA2,sqrt(sumsqr(UVW)),GCOUNT() FROM $ms GROUPBY ANTENNA1,ANTENNA2") as BL:
    ants1, ants2 = BL.getcol('ANTENNA1'), BL.getcol('ANTENNA2')
    dists = BL.getcol('Col_3')/1e3 # baseleline length in km
    n_t = BL.getcol('Col_4')[0] # number of timesteps
    n_bl = len(ants1)

# check if ms is time-ordered
times = ms.getcol('TIME_CENTROID')
if not all(np.diff(times) >= 0):
    logging.critical('This code cannot handle MS that are not time-sorted.')
    sys.exit(1)
del times

# create column to smooth
addcol(ms, options.incol, options.outcol, setvals=False) # setvals=F reduces overhead
# restore WEIGHT_SPECTRUM
if 'WEIGHT_SPECTRUM_ORIG' in ms.colnames() and options.restore:
    addcol(ms, 'WEIGHT_SPECTRUM_ORIG', 'WEIGHT_SPECTRUM')
# backup WEIGHT_SPECTRUM
elif options.weight and not options.nobackup:
    addcol(ms, 'WEIGHT_SPECTRUM', 'WEIGHT_SPECTRUM_ORIG')

# Iterate over chunks of baselines
for c, idx in enumerate(np.array_split(np.arange(n_bl), options.chunks)):
    
    logging.debug('### Fetching chunk {}/{}'.format(c,options.chunks))

    # get input data for this chunk
    ants1_chunk, ants2_chunk = ants1[idx], ants2[idx]
    chunk = pt.taql("SELECT FROM $ms WHERE any(ANTENNA1== $ants1_chunk && ANTENNA2==$ants2_chunk)")
    data_chunk = chunk.getcol(options.incol)
    weights_chunk = chunk.getcol('WEIGHT_SPECTRUM')
    # flag NaNs and set weights to zero
    flags = chunk.getcol('FLAG')
    flags[np.isnan(data_chunk)] = True
    weights_chunk[flags] = 0
    del flags
    # prepare output cols
    smoothed_data = data_chunk.copy()
    if options.weight:
        new_weights = np.zeros_like(weights_chunk)

    # Iterate on each baseline in this chunk
    for i_chunk, (ant1, ant2, dist) in enumerate(zip(ants1_chunk, ants2_chunk, dists[idx])):
        if ant1 == ant2: continue # skip autocorrelations
        elif np.isnan(dist): continue # fix for missing antennas
        logging.debug('Working on baseline: {} - {} (dist = {:.2f}km)'.format(ant1, ant2, dist))

        in_bl = slice(i_chunk,  -1, len(ants1_chunk)) # All times for 1 BL
        data = data_chunk[in_bl]
        weights = weights_chunk[in_bl]

        stddev_t = options.ionfactor * (25.e3 / dist) ** options.bscalefactor * (freq / 60.e6)  # in sec
        stddev_t = stddev_t / timepersample  # in samples
        # TODO: for freq this is hardcoded, it should be thought better
        # However, the limitation is probably smearing here
        stddev_f = 1e6 / dist  # Hz
        stddev_f = stddev_f / freqpersample  # in samples
        logging.debug("-Time: sig={:.1f} samples ({:.1f}s) -Freq: sig={:.1f} samples ({:.2f}MHz)".format(
            stddev_t, timepersample * stddev_t, stddev_f, freqpersample * stddev_f/1e6))

        if stddev_t == 0: continue  # fix for flagged antennas
        if stddev_t < 0.5: continue  # avoid very small smoothing

        # Multiply every element of the data by the weights, convolve both the
        # scaled data and the weights, and then divide the convolved data by the
        # convolved weights (translating flagged data into weight=0).
        # That's basically the equivalent of a running weighted average with a
        # Gaussian window function.

        # set bad data to 0 so nans do not propagate
        data = np.nan_to_num(data * weights)

        if np.isnan(data).all(): continue # flagged ants
        # smear weighted data and weights
        if options.onlyamp: # smooth only amplitudes
            dataAMP, dataPH = np.abs(data), np.angle(data)
            if not options.notime:
                dataAMP = smooth_real1d(dataAMP, stddev_t, axis=0, truncate=3)
            if not options.nofreq:
                dataAMP = smooth_real1d(dataAMP, stddev_f, axis=1, truncate=3)
            data = dataAMP * (np.cos(dataPH) + 1j * np.sin(dataPH)) # recreate data
        else:
            if not options.notime and not options.nofreq: # smooth in t and f
                data = smooth_imag1d(data, stddev_t, 0)
                data = smooth_imag1d(data, stddev_f, 1)
            elif not options.notime:
                data = smooth_imag1d(data, stddev_t, axis=0)
            elif not options.nofreq:
                data = smooth_imag1d(data, stddev_f, axis=1)
        if not options.notime:
            weights = smooth_real1d(weights, stddev_t, axis=0, truncate=3)
        if not options.nofreq:
            weights = smooth_real1d(weights, stddev_f, axis=1, truncate=3)
        data[(weights != 0)] /= weights[(weights != 0)]  # avoid divbyzero

        # print( np.count_nonzero(data[~flags[in_bl]]), np.count_nonzero(data[flags[in_bl]]), 100*np.count_nonzero(data[flags[in_bl]])/np.count_nonzero(data))
        # print( "NANs in flagged data: ", np.count_nonzero(np.isnan(data[flags[in_bl]])))
        # print( "NANs in unflagged data: ", np.count_nonzero(np.isnan(data[~flags[in_bl]])))
        # print( "NANs in weights: ", np.count_nonzero(np.isnan(weights)))
        smoothed_data[in_bl] = data
        if options.weight:
            new_weights[in_bl] = weights
    logging.info('Writing %s column.' % options.outcol)
    chunk.putcol(options.outcol, smoothed_data)
    if options.weight:
        logging.warning('Writing WEIGHT_SPECTRUM column.')
        chunk.putcol('WEIGHT_SPECTRUM', new_weights)
    chunk.close()
ms.close()
logging.info("Done.")