#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016 - Francesco de Gasperin
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

# Usage: jybeam2jyarcsec.py flux bmaj bmin
# Convert flux from jy/beam to jy/arcsec

import os, sys
import optparse, itertools
import logging
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d as gfilter
import pyrap.tables as pt
logging.basicConfig(level=logging.DEBUG)

def jybeam2jyarcsec(flux, beam):
    """
    flux -- Flux in jansky
    beam -- [min, maj] in arcsec
    """
    g=2*np.pi/(2.0*np.sqrt(2.0*np.log(2.0)))**2
    return flux/(g*beam[0]*beam[1])

if len(sys.argv) < 3:
    print "Use", sys.argv[0], " flux [jy], bmaj [arcsec], bmin [arcsec]"
    sys.exit(0)

flux = float(sys.argv[1])
bmaj = float(sys.argv[2])
bmin = float(sys.argv[3])

print "Flux in Jy/beam: %f (beam: %f arcsec x %f arcsec)" % (flux, bmaj, bmin)
print "Flux in Jy/arcsec: %f" % jybeam2jyarcsec(flux, [bmaj,bmin])
