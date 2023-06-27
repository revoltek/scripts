#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2022 - Francesco de Gasperin
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

# Usage: ./elev_from_coord.py ra [deg] dec [deg]
# Calculate the elevation at each LST [assuming LOFAR]

import sys, os, re
import numpy as np
from casacore import tables
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u

import matplotlib as mpl
import matplotlib.pyplot as plt

lofar = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m)

sources = []
for source in sys.argv[1:]:
    print(source)
    if source.count(',') == 2:
        name, ra, dec = source.split(',')
    elif source.count(',') == 1:
        ra, dec = source.split(',')
        name = 'No name'
    else:
        print('Use ./elev_from_coord.py [name,]ra,dec')
        sys.exit()

    sources.append([name,float(ra),float(dec)])

fig = plt.figure(figsize=(8, 12))
fig.subplots_adjust(wspace=0)
ax = fig.add_subplot(111)

for source in sources:
    name = source[0]
    coord_ra = source[1]
    coord_dec = source[2]

    coord = SkyCoord(coord_ra*u.deg, coord_dec*u.deg)
    
    midnight = Time('2012-7-13 00:00:00')
    delta_midnight = np.linspace(-12, 12, 240)*u.hour
    elev = coord.transform_to(AltAz(obstime=midnight+delta_midnight,location=lofar)).alt.deg
    lst = (midnight+delta_midnight).sidereal_time('mean', lofar.lon).deg
    for l,e in zip(lst,elev):
        if e < 40:
            print('LST: %03.1f hr -> %02d' % (l/15.,e))
        else:
            print('LST: %03.1f hr -> %02d *** GOOD' % (l/15.,e))

    ax.plot(lst/15, elev, marker='.', ls='', label=name)
    
ax.tick_params('both', length=10, width=2, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.set_xlabel(r'LST [hr]')
ax.set_ylabel(r'Elevation [deg]')
ax.set_ylim(ymin=-10,ymax=90)
ax.set_xlim(xmin=0,xmax=24)
ax.hlines(40, xmin=0, xmax=25, ls=':', color='red')
ax.hlines(0, xmin=0, xmax=25, ls='-', color='red')
ax.label_outer()
ax.legend(loc='upper right', frameon=False)
plt.show(block=True)
