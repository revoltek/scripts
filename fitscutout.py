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

import os, sys
import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from lib_fits import flatten, find_freq
import argparse

parser = argparse.ArgumentParser(description='Make a cutout of a fits-image')
parser.add_argument('filename', type=str)
parser.add_argument('-s', '--size', type=float, nargs='*', help='size in pixels')
parser.add_argument('-p', '--position', type=int, nargs=2, help='position in pixels')
parser.add_argument('--size_sky', type=float, nargs='*', help='size in arcmin')
parser.add_argument('--position_sky', type=float, nargs='*', help='position in degrees')
parser.add_argument('--skip', action='store_true', help='Skip existing cutouts?')
parser.add_argument('--freq', type=float, help='Manually set frequency')
parser.add_argument('-o', '--output', help='prefix and title of output image (default: input-cut.fits)')
args = parser.parse_args()

filename = args.filename
if args.output:
    cutout_filename = args.output
else:
    cutout_filename = filename.replace('.fits','-cut.fits')
if os.path.isfile(cutout_filename) and args.skip:
    print(f'Skipping {cutout_filename}')
    sys.exit(0)

header, data = flatten(filename)
# Load the image and the WCS
wcs = WCS(header)

if args.size and args.position:
    position = args.position
    size = args.size
    if len(size) == 1:
        size = [size[0], size[0]]
    elif len(size) == 2:
        pass
    else:
        raise ValueError('-s/--size requires one or two arguments.')
elif (args.size_sky and args.position_sky):
    position = args.position_sky
    size = np.array(args.size_sky)
    if len(size) == 1:
        size = [size[0], size[0]]
    elif len(size) == 2:
        pass
    else:
        raise ValueError('---size_sky requires one or two arguments.')
    position = wcs.wcs_world2pix([position[0]], [position[1]], 0)
    cdelt = np.abs(np.array(wcs.wcs.cdelt[0:2]) * 60)
    size = size / cdelt

    print(size)
else:
    raise ValueError('Provide either -s/--size and -p/--position or --size_sky and --position_sky.')

# Make the cutout, including the WCS
cutout = Cutout2D(data, position=position, size=size, wcs=wcs)

# Update the FITS header with the cutout WCS
hdr_new = cutout.wcs.to_header()
if args.freq:
    hdr_new["FREQ"] = args.freq
else:
    hdr_new["FREQ"] = find_freq(header)
header.update(hdr_new)
if os.path.isfile(cutout_filename):
    os.system(f'rm {cutout_filename}')
fits.writeto(cutout_filename, cutout.data, header, overwrite=False)
