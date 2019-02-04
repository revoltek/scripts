#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Virginia Cuciti
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

# Usage: ./sun_dist.py MSname
# Calculate the elevation of the sun at each hour during the observation and the minimum distance
# between the sun and the observed target

import sys, os, re
import numpy as np
from casacore import tables
from astropy.coordinates import get_sun, SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u

lofar = EarthLocation(lat=52.90889*u.deg, lon=6.86889*u.deg, height=0*u.m)

tb=tables.table(sys.argv[1], ack=False)
times = tb.getcol('TIME')
time=Time(np.arange(times[0], times[-1], 3600)/86400,format='mjd')
sun=get_sun(time)

sun_elevation=sun.transform_to(AltAz(obstime=time,location=lofar)).alt
T = Time(time, format='iso', scale='utc')

for i, t in enumerate(T):
        print "Sun elevation at time %s is: %s" %(t, str(sun_elevation[i]))

tb=tables.table(sys.argv[1]+'/FIELD', ack=False)
phase_dir=tb.getcol('PHASE_DIR')[0,0]

coord=SkyCoord(float(phase_dir[0]*180./np.pi)*u.degree, float(phase_dir[1]*180./np.pi)*u.degree)
dist=coord.separation(sun)
print "Minimum distance from the Sun:", np.min(dist)

