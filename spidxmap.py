#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 - Francesco de Gasperin
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

import os, sys, argparse, logging
from lib_linearfit import linear_fit_bootstrap as linearfit
from lib_fits import flatten
from lib_beamdeconv import *
from astropy.io import fits as pyfits
from astropy import convolution
import pyregion
# https://github.com/astrofrog/reproject
from reproject import reproject_interp, reproject_exact
reproj = reproject_exact
logging.root.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser(description='Mosaic ddf-pipeline directories')
parser.add_argument('images', nargs='+', help='List of images to use for spidx')
parser.add_argument('--beam', dest='beam', nargs='+', type=float, help='3 parameters final beam to convolve all images in degrees (BMAJ, BMIN, BPA)')
parser.add_argument('--region', dest='region', help='Ds9 region to restrict analysis')
parser.add_argument('--size', dest='size', type=float, help='Size of final image in degree')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx')
parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image')
parser.add_argument('--output', dest='output', default='spidx.fits', help='Name of output mosaic (default: mosaic.fits)')

args = parser.parse_args()

# check input
if len(args.images) < 2:
    logging.error('Requires at lest 2 images.')
    sys.exit(1)

if len(args.beam) != 3:
    logging.error('Beam must be in the form of "BMAJ BMIN BPA" (3 floats).')
    sys.exit(1)

if len(args.radec) != 2:
    logging.error('--radec must be in the form of "RA DEC" (2 floats).')
    sys.exit(1)

from lib_fits import Image
class ImageSpidx(Image):

    def __init__(self, imagefile):
        Image.__init__(self, imagefile)

    def concolve(self, target_beam):
        """
        Convolve *to* this rsolution
        beam = [bmaj, bmin, bpa]
        """
        # first find beam to convolve with
        convolve_beam = deconvolve_ell(target_beam[0], target_beam[1], target_beam[2], self.beam[0], self.beam[1], self.beam[2])
        logging.debug('%s - convolve_beam: %f %f %f' % (self.fistfile, convolve_beam[0], convolve_beam[1], convolve_beam[2]))
        # do convolution on data
        bmaj, bmin, bpa = convolve_beam
        gauss_kern = EllipticalGaussian2DKernel(bmaj, bmin, bpa)
        self.data = convolution.convolve(self.data, gauss_kern, boundary=None)

    def regrid(self, regrid_hdr):
        logging.debug('%s - regridding' % (self.fistfile))
        self.data_regrid, __footprint = reproj((d.img_data, d.img_hdr), regrid_hdr, parallel=True)
        
    def calc_shift(self, ref_cat):
        """
        Find a shift cross-matching source extracted from the image and a given catalog
        """
        # TODO: relative shift between images
        from lofar import bdsm
        from astropy.coordinates import match_coordinates_sky
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        img_cat = self.imagefile+'.cat'
        bdsm_img = bdsm.process_image(self.imagefile, rms_box=(100,30), \
            thresh_pix=5, thresh_isl=3, atrous_do=False, \
            adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10), quiet=True)
        bdsm_img.write_catalog(outfile=img_cat, catalog_type='srl', format='fits', clobber=True)

        # read catlogue
        ref_t = Table.read(ref_cat)
        img_t = Table.read(img_cat)

        # cross match
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(ref_t['RA'], ref_t['DEC']),\
                                                  SkyCoord(img_t['RA'], img_t['DEC']))
        idx_match = idx_match[sep<15*u.arcsec]

        # find & apply shift
        if len(idx_match) == 0:
            logging.warning('No match found in TGSS.')
            return
        dra = np.mean(ref_t['RA'][idx_match] - img_t['RA'][idx_match])
        ddec = np.mean(ref_t['DEC'][idx_match] - img_t['DEC'][idx_match])
        logging.debug('Shift for %s: %f %f (arcsec)' % (self.imagefile, dra*3600, ddec*3600))
        ra = self.img_hdr['CRVAL1']
        dec = self.img_hdr['CRVAL2']
        self.img_hdr['CRVAL1'] -= dra/(3600.*np.cos(np.pi*dec/180.))
        self.img_hdr['CRVAL2'] -= ddec/3600.

        # clean up
        os.system('rm '+img_cat)


########################################################
# prepare images and shift
all_images = []
all_beams = []
for imagefile in args.images:
    image = Image(imagefile)
    if args.shift:
        image.calc_shift()
    all_beams.append(image.beam)
    all_images.append(image)

#####################################################
# find the smallest common beam
if args.beam is not None:
    target_beam = lib_findCommonBeam(all_beams)
else:
    target_beam = args.beam

logging.info('Final beam: maj=%f min=%f pa=%f' % (target_beam[0]*3600., target_beam[1]*3600., target_beam[2]*3600.))

######################################################
# Generate regrid headers
rwcs = pywcs(naxis=2)
rwcs.wcs.ctype = all_images[0].get_wcs().wcs.ctype
cdelt = np.rint(target_beam[1]/5.) # 1/5 of minor axes (deg)
logging.info('Pixel scale: %i"' % cdelt/3600.)
rwcs.wcs.cdelt = cdelt
if args.radec is not None:
    mra = radec[0]*np.pi/180
    mdec = radec[1]*np.pi/180
else:
    mra = all_images[0].img_hdr['CVAL0']
    mdec = all_images[0].img_hdr['CVAL1']
rwcs.wcs.crval = [mra,mdec]

# if size is not give is taken from the mask
if args.size is None and args.region is not None:
    r = pyregion.open(args.region)
    mask = r.get_mask(header=all_images[0].img_hdr, shape=all_images[0].img_data.shape)
    y, x = mask.nonzero()
    ra_max, dec_max = w.all_pix2world(np.max(x), np.max(y), 0, ra_dec_order=True)
    ra_min, dec_min = w.all_pix2world(np.min(x), np.min(y), 0, ra_dec_order=True)
    args.size = np.max( [ np.abs(ra_max-ra_min), np.abs(dec_max-dec_min) ] )
else:
    logging.warning('No size or region provided, use entire size of first image.')
    args.size = 1 # TODO
    sys.exit('not implemented')
logging.info('Image size: %f deg' % args.size)
xsize = args.size/cdelt
ysize = args.size/cdelt
rwcs.wcs.crpix = [xsize/2,ysize/2]

regrid_hdr = rwcs.to_header()
regrid_hdr['NAXIS'] = 2
regrid_hdr['NAXIS1'] = xsize
regrid_hdr['NAXIS2'] = ysize

#########################################################
# convolve and only after apply mask and find noise
for image in all_images:
    image.convolve(target_beam)

    if args.region is not None:
        image.applyRegion(args.region, invert=True) # after convolution
    if args.noise:
        image.calc_noise() # after mask/convolution

    image.regrid(regrid_hdr) # finally regrid
 
#########################################################
# do spdix and write output
frequencies = [ image.freq for image in all_images ]
if args.noise: yerr = [ image.noise for image in all_images ]
else: yerr = None
spidx_data = np.nans(shape=(xsize,ysize))
spidx_err_data = np.nans(shape=(xsize,ysize))

for i in xrange(xsize):
    for j in xrange(ysize):
        val4reg = [ image.img_data[i,j] for image in all_images ]
        if np.isnan(val4reg).any(): continue
        (a, b, sa, sb) = linearfit(x=np.log10(frequencies), y=np.log10(val4reg), yerr=yerr)
        print '.',

spidx = fits.PrimaryHDU(spidx_data, regrid_hdr)
spidx_err = fits.PrimaryHDU(spidx_err_data, regrid_hdr)
spidx.writeto(args.output)
spidx_err.writeto(args.output.replace('.fits','-err.fits'))
