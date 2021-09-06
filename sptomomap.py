#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021 - Francesco de Gasperin
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

import os, sys, argparse, logging, copy
import numpy as np
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
from lib_linearfit import linear_fit, linear_fit_bootstrap
from lib_fits import AllImages
# https://github.com/astrofrog/reproject
from reproject import reproject_interp, reproject_exact
reproj = reproject_exact
logging.root.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser(description='Make spectral tomography maps, e.g. sptomomap.py --region ds9.reg --noise --sigma 5 --save *fits')
parser.add_argument('images', nargs=2, help='Two images to use for tomography')
parser.add_argument('--beam', dest='beam', nargs=3, type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
parser.add_argument('--alpha', dest='alpha', nargs=3, type=float, help='alpha min, max and num. Num is equal to the number of slices to make.')
parser.add_argument('--region', dest='region', type=str, help='Ds9 region to restrict analysis')
parser.add_argument('--noiseregion', dest='noiseregion', type=str, help='Ds9 region to calculate rms noise (default: do not use)')
parser.add_argument('--noisesigma', dest='noisesigma', default=5, type=int, help='Sigma used in the calc_noise function when no region is specified (default: 5)')
parser.add_argument('--size', dest='size', nargs=2, type=float, help='Size (ra and dec) of final image in degree (example: 3.5 4.0)')
parser.add_argument('--radec', dest='radec', nargs=2, type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx (default: false)')
parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image, necessary for the error map (default: false)')
parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results (default: false)')
parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
parser.add_argument('--circbeam', dest='circbeam', action='store_true', help='Force final beam to be circular (default: False, use minimum common beam area)')
parser.add_argument('--output', dest='output', default='sptomo.fits', type=str, help='Name of output mosaic (default: sptomo.fits)')

args = parser.parse_args()

alpha_min = args.alpha[0]
alpha_max = args.alpha[1]
alpha_num = int(args.alpha[2])

########################################################
# prepare images and convolve+regrid if needed
if np.all([os.path.exists(name.replace('.fits', '-conv-regrid.fits')) for name in args.images]):
    logging.info('Found convolved+regridded image... restoring.')
    all_images = AllImages([name.replace('.fits', '-conv-regrid.fits') for name in args.images])
    regrid_hdr = copy.copy(all_images.images[0].img_hdr)
    del regrid_hdr['FREQ']
    del regrid_hdr['RESTFREQ']
else:
    all_images = AllImages(args.images)
    all_images.convolve_to(beam=args.beam, circbeam=args.circbeam)
    if args.save: all_images.write('conv')
    regrid_hdr = all_images.regrid_common(size=args.size, region=args.region, radec=args.radec, action='regrid_header')
    if args.save: all_images.write('conv-regrid')

#####################################################
# find+apply shift w.r.t. lowest noise image
if args.shift:
    all_images.align_catalogue()
 
#########################################################
# only after regrid+convolve apply mask and find noise
for image in all_images:
    if args.noise:
        if args.sigma is not None:
            # usually the sigma used for the blanking is rather low, better to increase it for the calc_noise
            image.calc_noise(sigma=args.noisesigma, bg_reg=args.noiseregion, force_recalc=True) # after convolution
            image.blank_noisy(args.sigma)
        else:
            image.calc_noise(force_recalc=True) # after mask?/convolution

    if args.region is not None:
        image.apply_region(args.region, invert=True) # after convolution to minimise bad pixels

#########################################################
# do spdix and write output
xsize = regrid_hdr['NAXIS1']
ysize = regrid_hdr['NAXIS2']
frequencies = [ image.get_freq() for image in all_images ]

if args.noise: yerr = [ image.noise for image in all_images ]
else: yerr = None
tomo_data = np.empty(shape=(alpha_num,ysize,xsize))
tomo_data[:] = np.nan
tomo_err_data = np.empty(shape=(alpha_num,ysize,xsize))
tomo_err_data[:] = np.nan

logging.info("Alphas: "+str(np.linspace(alpha_min, alpha_max, num=alpha_num)))
for a, alpha in enumerate(np.linspace(alpha_min, alpha_max, num=alpha_num)):
    if frequencies[0] < frequencies[1]:
        tomo = all_images.images[0].img_data - ((frequencies[0]/frequencies[1])**alpha) * all_images.images[1].img_data
        tomo_err = 0
    else:
        tomo = all_images.images[1].img_data - ((frequencies[1]/frequencies[0])**alpha) * all_images.images[0].img_data 
        tomo_err = 0
    tomo_data[a] = tomo
    tomo_err_data[a] = tomo_err

if 'FREQ' in regrid_hdr.keys():
    del regrid_hdr['FREQ']
if 'RESTFREQ' in regrid_hdr.keys():
    del regrid_hdr['RESTFREQ']
regrid_hdr['BTYPE'] = 'SPTOMO'

# add headers for alpha axes
regrid_hdr['NAXIS'] = 3
regrid_hdr['NAXIS3'] = alpha_num
regrid_hdr['CRPIX3'] = 1.
regrid_hdr['CRVAL3'] = 1.#alpha_min
regrid_hdr['CDELT3'] = 1.#(alpha_max-alpha_min)/alpha_num
regrid_hdr['CUNIT3'] = ' '
del regrid_hdr['WCSAXES']

if args.output[-5:] == '.fits':
    filename_out = args.output
else:
    filename_out = args.output+'.fits'
logging.info('Save %s (and errors)' % filename_out)
regrid_hdr['CTYPE3'] = 'ALPHA'
pyfits.writeto(filename_out, tomo_data, regrid_hdr, overwrite=True, output_verify='fix')
regrid_hdr['CTYPE3'] = 'ALPHAERR'
pyfits.writeto(filename_out.replace('.fits','-err.fits'), tomo_err_data, regrid_hdr, overwrite=True, output_verify='fix')
