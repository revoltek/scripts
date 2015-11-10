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

def addcol(ms, incol, outcol):
    if outcol not in ms.colnames():
        logging.info('Adding column: '+outcol)
        coldmi = ms.getdminfo(incol)
        coldmi['NAME'] = outcol
        datatype = ms.col(incol).datatype()
        print datatype
        ms.addcols(pt.maketabdesc(pt.makearrcoldesc(outcol, 0., valuetype=datatype, shape=np.array(ms.getcell(incol,0)).shape)), coldmi)
    # copy columns val
    logging.info('Set '+outcol+'='+incol)
    data = ms.getcol(incol)
    ms.putcol(outcol, data)

opt = optparse.OptionParser(usage="%prog [options] MS", version="%prog 0.1")
opt.add_option('-f', '--ionfactor', help='Gives an indication on how strong is the ionosphere [default: 0.2]', type='float', default=0.2)
opt.add_option('-i', '--incol', help='Column name to smooth [default: DATA]', type='string', default='DATA')
opt.add_option('-o', '--outcol', help='Output column [default: DATA_SMOOTH]', type="string", default='DATA_SMOOTH')
opt.add_option('-w', '--weight', help='Save the newly computed WEIGHT_SPECTRUM, this action permanently modify the MS! [default: False]', action="store_true", default=False)
opt.add_option('-r', '--restore', help='If WEIGHT_SPECTRUM_ORIG exists then restore it before smoothing [default: False]', action="store_true", default=False)
opt.add_option('-b', '--nobackup', help='Do not backup the old WEIGHT_SPECTRUM in WEIGHT_SPECTRUM_ORIG [default: do backup if -w]', action="store_true", default=False)
(options, msfile) = opt.parse_args()
ionfactor = options.ionfactor

if msfile == []:
    opt.print_help()
    sys.exit(0)
msfile = msfile[0]

if not os.path.exists(msfile):
    logging.error("Cannot find MS file.")
    sys.exit(1)

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

# create column to smooth
addcol(ms, options.incol, options.outcol)

# retore WEIGHT_SPECTRUM
if 'WEIGHT_SPECTRUM_ORIG' in ms.colnames() and options.restore:
    addcol(ms,'WEIGHT_SPECTRUM_ORIG','WEIGHT_SPECTRUM')
# backup WEIGHT_SPECTRUM
elif options.weight and not options.nobackup:
    addcol(ms, 'WEIGHT_SPECTRUM', 'WEIGHT_SPECTRUM_ORIG')

ant1 = ms.getcol('ANTENNA1')
ant2 = ms.getcol('ANTENNA2')
all_data = ms.getcol(options.outcol)
all_weights = ms.getcol('WEIGHT_SPECTRUM')
all_flags = ms.getcol('FLAG')
all_uvw = ms.getcol('UVW')

all_flags[ np.isnan(all_data) ] = True # flag NaNs
all_weights = all_weights * ~all_flags # set weight of flagged data to 0
    
# Check that all NaNs are flagged
if np.count_nonzero(np.isnan(all_data[~all_flags])) > 0:
    logging.error('NaNs in unflagged data!')

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

    #Multiply every element of the data by the weights, convolve both the scaled data and the weights, and then
    #divide the convolved data by the convolved weights (translating flagged data into weight=0). That's basically the equivalent of a
    #running weighted average with a Gaussian window function.

    # get cycle values
    weights = all_weights[sel,:,:]
    data = all_data[sel,:,:,]

    # get data, and set bad data to 0 so nans do not propagate
    data = np.nan_to_num(data*weights)

    # smear weighted data and weights
    dataR = gfilter(np.real(data), stddev, axis=0)#, truncate=4.)
    dataI = gfilter(np.imag(data), stddev, axis=0)#, truncate=4.)
    weights = gfilter(weights, stddev, axis=0)#, truncate=4.)

    # re-create data
    data = (dataR + 1j * dataI)
    data[(weights != 0)] /= weights[(weights != 0)] # avoid divbyzero
    all_data[sel,:,:] = data
    all_weights[sel,:,:] = weights

#    print np.count_nonzero(data[~flags]), np.count_nonzero(data[flags]), 100*np.count_nonzero(data[flags])/np.count_nonzero(data)
#    print "NANs in flagged data: ", np.count_nonzero(np.isnan(data[flags]))
#    print "NANs in unflagged data: ", np.count_nonzero(np.isnan(data[~flags]))
#    print "NANs in weights: ", np.count_nonzero(np.isnan(weights))

ms.putcol(options.outcol, all_data)
ms.putcol('FLAG', all_flags) # this saves flags of nans, which is always good
if options.weight:
    logging.warning('Changing WEIGHT_SPECTRUM column.')
    ms.putcol('WEIGHT_SPECTRUM', all_weights)

ms.close()
logging.info("Done.")
