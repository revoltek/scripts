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

# flux2lum.py redshift flux(in Jy) [error]

import sys, os
import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.27)
print "Using Flat LmbdaCDM H0=0.73 Om0=0.27"

z = float(sys.argv[1]) # redshift
flux = float(sys.argv[2]) # in Jy
if len(sys.argv) > 2: flux_err = float(sys.argv[3])
else: flux_err = 0.

dist = cosmo.luminosity_distance(z).value # in Mpc
print "Distance: ", dist, "Mpc"

# 1 pc = 3.08567758e16 m
print flux * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "±", flux_err * 1e-26 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "W/Hz"
print flux * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "±", flux_err * 1e-26 * 1e7 * ( 4*np.pi * (dist * 1e6 * 3.08567758e16)**2),
print "erg/s/Hz"
print "(For a distance of "+str(dist)+" Mpc)"
