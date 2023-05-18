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

coord_ra = float(sys.argv[1])
coord_dec = float(sys.argv[2])

lofar = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m)
coord = SkyCoord(coord_ra*u.deg, coord_dec*u.deg)

midnight = Time('2012-7-13 00:00:00')
delta_midnight = np.linspace(-12, 12, 24)*u.hour
elev = coord.transform_to(AltAz(obstime=midnight+delta_midnight,location=lofar)).alt.deg
lst = (midnight+delta_midnight).sidereal_time('mean', lofar.lon).deg
for l,e in zip(lst,elev):
    if e < 40:
        print('LST: %03.1f hr -> %02d' % (l/15.,e))
    else:
        print('LST: %03.1f hr -> %02d *** GOOD' % (l/15.,e))
