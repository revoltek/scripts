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

# Usage: flipfitsaxis.py -a x|y fitsfile

import sys, optparse
import numpy as np
import pyfits

opt = optparse.OptionParser(usage="%prog [-a #] fitsfile", version="%prog 0.1")
opt.add_option('-a', '--axis', help='select which image axis to flip (give an integer)', default=None)
(options, img) = opt.parse_args()
axis = options.axis

try:
	hdulist = pyfits.open(img[0], mode='update')
except:
	print("ERROR: problems opening file "+img[0])
	sys.exit(1)

if ( axis is None ):
	print("ERROR: no axis selected, use '-a x' or '-a y'.")
	sys.exit(1)

try:
    axis = str(int(axis))
except:
	print("ERROR: axis option must be an integer")
	sys.exit(1)

# flip data
data = hdulist[0].data
if len(data.shape) > 2:
	print("ERROR: this program works only with 2D data, here we have: "+str(data.shape))
	sys.exit(1)

print("Image shape:", data.shape)

if axis == '1':
    hdulist[0].data = np.copy(data[:,::-1])
elif axis == '2':
    hdulist[0].data = np.copy(data[::-1])
else:
    print("ERROR: axis must be 1 or 2")
    sys.exit(1)

prihdr = hdulist[0].header
print("Flipping axis "+axis+" ("+prihdr.get('CTYPE'+axis).split('-')[0]+")")

# reset reference pixel
refpix = prihdr['CRPIX'+axis]
axlen = prihdr['NAXIS'+axis]
print("updating "+'CRPIX'+axis)
prihdr['CRPIX'+axis] = axlen-refpix # mirror reference pixel

# invert axis increase value
print("updating: "+'CD'+axis+'_'+axis)
prihdr['CD'+axis+'_'+axis] = -1 * prihdr['CD'+axis+'_'+axis]

hdulist.flush()
print("Done!")
