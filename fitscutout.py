#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 - Francesco de Gasperin
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

#./fitscutout.py fitsfile [set position and size below]
position = (1178, 1164) # pixel
size = (1852, 1882) # pixel

import os, sys
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from lib_fits import flatten

filename = sys.argv[1]
header, data = flatten(filename)

# Load the image and the WCS
wcs = WCS(header)

# Make the cutout, including the WCS
cutout = Cutout2D(data, position=position, size=size, wcs=wcs)

# Update the FITS header with the cutout WCS
header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = filename.replace('.fits','-cut.fits')
fits.writeto(cutout_filename, cutout.data, header, overwrite=False)
