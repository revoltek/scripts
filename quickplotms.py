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

# Usage: quickplotms vis.MS
# Load an MS data and create a quick plot

import os, sys, optparse
import logging
import numpy as np
import pyrap.tables as pt
import matplotlib.pyplot as plt
logging.basicConfig(level=logging.INFO)

opt = optparse.OptionParser(usage="%prog [options] MS", version="%prog 0.1")
opt.add_option('-s', '--save', help='Instead of showing the image it saves it in this filename (png) [default: '']', type='string', default='')
opt.add_option('-c', '--col', help='Column name to plot [default: DATA]', type='string', default='DATA')
opt.add_option('-a', '--ant', help='Baseline to plot [default: 0&1]', type="string", default='0&1')
opt.add_option('-p', '--pol', help='Pol to plot, accept: 0,1,2,3 [default: 0]', type="int", default=0)
opt.add_option('-n', '--chan', help='Which chan to plot [default: 0]', type="int", default=0)
opt.add_option('-f', '--flag', help='Plot flags? [default: False]', action="store_true", default=False)
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
ant1 = int(options.ant.split('&')[0])
ant2 = int(options.ant.split('&')[1])
msb = pt.taql('select from $ms where ANTENNA1 = $ant1 and ANTENNA2 = $ant2')
data = np.absolute(msb.getcol(options.col))
flags = msb.getcol('FLAG')
flags[ np.isnan(data) ] = True # flag NaNs
time = msb.getcol('TIME')
    
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(wspace=0)
ax = fig.add_subplot(111)
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.set_xlabel(r'Time [s]')
ax.set_ylabel(r'Amplitude [Jy]')
ax.label_outer()
print('Shape:', data.shape)
ax.plot( time[~flags[:,options.chan,options.pol]], data[:,options.chan,options.pol][~flags[:,options.chan,options.pol]], 'k,', ls='' )
if options.flag: ax.plot( time[flags[:,options.chan,options.pol]], data[:,options.chan,options.pol][flags[:,options.chan,options.pol]], 'r,', ls='' ) # flags

if options.save != '':
    logging.info('Save file: '+options.save)
    fig.savefig(options.save, bbox_inches='tight')
    fig.clf()
else:
    plt.show()

logging.info("Done.")
