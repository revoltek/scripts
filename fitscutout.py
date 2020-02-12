#!/usr/bin/env python
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
position = (500, 500) # pixel
size = (400, 400) # pixel

import os, sys
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

filename = sys.argv[1]

# Load the image and the WCS
data, hdr = fits.open(filename)[0]
wcs = WCS(hdu.header)

# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

# Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = filename.replace('.fits','-cut.fits')
hdu.writeto(cutout_filename, overwrite=False)
