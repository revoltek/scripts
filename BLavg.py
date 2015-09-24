#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2015 - Francesco de Gasperin
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

# Usage: BLavg.py vis.MS
# Load a MS, average visibilities according to the baseline leght,
# i.e. shorter BLs are averaged more, and write a new MS

import os, sys
import optparse, itertools
import logging
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import pyrap.tables as pt
logging.basicConfig(level=logging.DEBUG)

opt = optparse.OptionParser(usage="%prog [options] MS", version="%prog 0.1")
opt.add_option('-i', '--ionfactor', help='Gives an indication on how strong is the ionosphere [default: 0.2]', type='float', default=0.2)
opt.add_option('-o', '--overwrite', help='If active, overwrite the input MS or create a new one [default: False, output is "inMS"-BLavg.MS]', action="store_true", default=False)
opt.add_option('-c', '--clobber', help='If active, delete the output file if different from input and exists [default: False]', action="store_true", default=False)
opt.add_option('-m', '--memory', help='Work with all the data in memory, much faster [default: False]', action="store_true", default=False)
opt.add_option('-l', '--column', help='Column name to average, output will always be the same column [default: DATA]', type='string', default='DATA')
(options, msfile) = opt.parse_args()
ionfactor = options.ionfactor

if msfile == []:
    opt.print_help()
    sys.exit(0)
msfile = msfile[0]

if not os.path.exists(msfile):
    logging.error("Cannot find MS file.")
    sys.exit(1)

# prepare new ms
if not options.overwrite:
    msfile_new = msfile.replace('.MS','-BLavg.MS')
    if os.path.exists(msfile_new):
        if not options.clobber:
            logging.error("Output file exists and clobber=False")
            sys.exit(1)
        os.system('rm -r '+msfile_new)
    logging.info("Copying MS, this may take a while.")
    os.system('cp -r '+msfile+' '+msfile_new)
    msfile = msfile_new

# open input/output MS
ms = pt.table(msfile, readonly=False, ack=False)
        
freqtab = pt.table(msfile + '/SPECTRAL_WINDOW', ack=False)
freq = freqtab.getcol('REF_FREQUENCY')
freqtab.close()
wav = 299792458. / freq
timepersample = ms.getcell('INTERVAL',0)
all_time = ms.getcol('TIME_CENTROID')

# check if ms is time-ordered
if not all(all_time[i] <= all_time[i+1] for i in xrange(len(all_time)-1)):
    logging.critical('This code cannot handle MS that are not time-sorted.')
    sys.exit(1)

# check if loading all data in memory
if options.memory:
    ant1 = ms.getcol('ANTENNA1')
    ant2 = ms.getcol('ANTENNA2')
    all_data = ms.getcol(options.column)
    all_weights = ms.getcol('WEIGHT_SPECTRUM')
    all_flags = ms.getcol('FLAG')
    all_uvw = ms.getcol('UVW')

    # iteration on baseline combination
    for ant in itertools.product(set(ant1), set(ant2)):

        if ant[0] >= ant[1]: continue
        sel1 = np.where(ant1 == ant[0])[0]
        sel2 = np.where(ant2 == ant[1])[0]
        sel = sorted(list(frozenset(sel1).intersection(sel2)))
    
        # compute the FWHM
        uvw = all_uvw[sel,:]
        uvw_dist = np.sqrt(uvw[:, 0]**2 + uvw[:, 1]**2 + uvw[:, 2]**2)
        dist = np.mean(uvw_dist) / 1.e3
        stddev = options.ionfactor * np.sqrt((25.e3 / dist)) * (freq / 60.e6) # in sec
        stddev = stddev/timepersample # in samples
        logging.debug("For BL %i - %i (dist = %.1f km): sigma=%.2f samples." % (ant[0], ant[1], dist, stddev))
    
#    Multiply every element of the data by the weights, convolve both the scaled data and the weights, and then
#    divide the convolved data by the convolved weights (translating flagged data into weight=0). That's basically the equivalent of a
#    running weighted average with a Gaussian window function.
    
        # get weights
        flags = all_flags[sel,:,:]
        weights = all_weights[sel,:,:]*~flags # set flagged data weight to 0
        # get data
        data = all_data[sel,:,:]*weights
    
        # smear weighted data and weights
        dataR = gfilter(np.real(data), stddev, axis=0)#, truncate=4.)
        dataI = gfilter(np.imag(data), stddev, axis=0)#, truncate=4.)
        weights = gfilter(weights, stddev, axis=0)#, truncate=4.)
    
        # re-create data
        all_data[sel,:,:] = (dataR + 1j * dataI)/weights # can I do it?
        all_weights[sel,:,:] = weights

    ms.putcol('DATA', all_data)
    ms.putcol('WEIGHT_SPECTRUM', all_weights)

# load BL per BL, cleaner but much slower
else:

    # iteration on baselines
    for t in ms.iter(["ANTENNA1", "ANTENNA2"]):
        ant1 = t.getcell('ANTENNA1', 0)
        ant2 = t.getcell('ANTENNA2', 0)
        if ant1 >= ant2: continue
        
        # compute the FWHM
        uvw = t.getcol('UVW')
        uvw_dist = np.sqrt(uvw[:, 0]**2 + uvw[:, 1]**2 + uvw[:, 2]**2)
        dist = np.mean(uvw_dist) / 1.e3
        stddev = options.ionfactor * np.sqrt((25.e3 / dist)) * (freq / 60.e6) # in sec
        stddev = stddev/timepersample # in samples
        logging.debug("For BL %i - %i (dist = %.1f km): sigma=%.2f samples." % (ant1, ant2, dist, stddev))
    
#    Multiply every element of the data by the weights, convolve both the scaled data and the weights, and then
#    divide the convolved data by the convolved weights (translating flagged data into weight=0). That's basically the equivalent of a
#    running weighted average with a Gaussian window function.
    
        # get weights
        flags = t.getcol('FLAG')
        weights = t.getcol('WEIGHT_SPECTRUM')*~flags # set flagged data weight to 0
        # get data
        data = t.getcol(options.column)*weights
    
        # smear weighted data and weights
        dataR = gfilter(np.real(data), stddev, axis=0)#, truncate=4.)
        dataI = gfilter(np.imag(data), stddev, axis=0)#, truncate=4.)
        weights = gfilter(weights, stddev, axis=0)#, truncate=4.)
    
        # re-create data
        data = (dataR + 1j * dataI)/weights # can I do it?
    
        # write the BL
        t.putcol(options.column, data)
        t.putcol('WEIGHT_SPECTRUM', weights)

ms.close()
logging.info("Done.")
