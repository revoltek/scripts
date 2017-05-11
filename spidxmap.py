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
from astropy.wcs import WCS as pywcs
from astropy import convolution
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
# https://github.com/astrofrog/reproject
from reproject import reproject_interp, reproject_exact
reproj = reproject_exact
logging.root.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser(description='Mosaic ddf-pipeline directories')
parser.add_argument('images', nargs='+', help='List of images to use for spidx')
parser.add_argument('--beam', dest='beam', nargs='+', type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
parser.add_argument('--region', dest='region', help='Ds9 region to restrict analysis')
parser.add_argument('--size', dest='size', type=float, help='Size (horizontal and vertical) of final image in degree')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx')
parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image')
parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results')
parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
parser.add_argument('--output', dest='output', default='spidx.fits', help='Name of output mosaic (default: mosaic.fits)')

args = parser.parse_args()

# check input
if len(args.images) < 2:
    logging.error('Requires at lest 2 images.')
    sys.exit(1)

if args.beam is not None and len(args.beam) != 3:
    logging.error('Beam must be in the form of "BMAJ BMIN BPA" (3 floats).')
    sys.exit(1)

if args.radec is not None and len(args.radec) != 2:
    logging.error('--radec must be in the form of "RA DEC" (2 floats).')
    sys.exit(1)

from lib_fits import Image
class ImageSpidx(Image):

    def __init__(self, imagefile):
        Image.__init__(self, imagefile)

    def convolve(self, target_beam):
        """
        Convolve *to* this rsolution
        beam = [bmaj, bmin, bpa]
        """
        # if difference between beam is negligible, skip - it mostly happens when beams are exactly the same
        beam = self.get_beam()
        if ((target_beam[0] - beam[0]) < 1e-4) and ((target_beam[1] - beam[1]) < 1e-4) and ((target_beam[2] - beam[2]) < 1):
            return
        # first find beam to convolve with
        convolve_beam = deconvolve_ell(target_beam[0], target_beam[1], target_beam[2], beam[0], beam[1], beam[2])
        if convolve_beam[0] is None: 
            logging.error('Cannot deconvolve this beam.')
            sys.exit(1)
        logging.debug('%s: Convolve beam: %.3f" %.3f" (pa %.1f deg)' \
                % (self.imagefile, convolve_beam[0]*3600, convolve_beam[1]*3600, convolve_beam[2]))
        # do convolution on data
        bmaj, bmin, bpa = convolve_beam
        assert abs(self.img_hdr['CDELT1']) == abs(self.img_hdr['CDELT2'])
        pixsize = abs(self.img_hdr['CDELT1'])
        fwhm2sigma = 1./np.sqrt(8.*np.log(2.))
        gauss_kern = EllipticalGaussian2DKernel((bmaj*fwhm2sigma)/pixsize, (bmin*fwhm2sigma)/pixsize, (90+bpa)*np.pi/180.) # bmaj and bmin are in pixels
        self.img_data = convolution.convolve(self.img_data, gauss_kern, boundary=None)
        self.img_data *= (target_beam[0]*target_beam[1])/(beam[0]*beam[1]) # since we are in Jt/b we need to renormalise
        self.set_beam(target_beam) # update beam

    def regrid(self, regrid_hdr):
        logging.debug('%s: regridding' % (self.imagefile))
        self.img_data, __footprint = reproj((self.img_data, self.img_hdr), regrid_hdr, parallel=True)
        beam = self.get_beam()
        self.img_hdr = regrid_hdr
        self.set_beam(beam) # retain beam info if not present in regrd_hdr

    def blank_noisy(self, nsigma):
        """
        Set to nan pixels below nsigma*noise
        """
        nans_before = np.sum(np.isnan(self.img_data))
        self.img_data[ np.isnan(self.img_data) ] = 0 # temporary set nans to 0 to prevent error in "<"
        self.img_data[ np.where(self.img_data <= nsigma*image.noise) ] = np.nan
        nans_after = np.sum(np.isnan(self.img_data))
        logging.debug('%s: Blanked pixels %i -> %i' % (self.imagefile, nans_before, nans_after))
        
    def make_catalogue(self):
        """
        Create catalogue for this image
        """
        import bdsf
        from astropy.table import Table

        img_cat = self.imagefile+'.cat'
        if not os.path.exists(img_cat):
            bdsf_img = bdsf.process_image(self.imagefile, rms_box=(100,30), \
                thresh_pix=5, thresh_isl=3, atrous_do=False, \
                adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10), quiet=True)
            bdsf_img.write_catalog(outfile=img_cat, catalog_type='srl', format='fits', clobber=True)
        else:
            logging.warning('%s already exists, using it.' % img_cat)

        self.cat = Table.read(img_cat)
        logging.debug('%s: Number of sources detected: %i' % (self.imagefile, len(self.cat)) )

    def apply_shift(self, dra, ddec):
        """
        Shift header by dra/ddec
        """
        logging.debug('%s: Shift %.2f %.2f (arcsec)' % (self.imagefile, dra*3600, ddec*3600))
        dec = self.img_hdr['CRVAL2']
        self.img_hdr['CRVAL1'] += dra/(np.cos(np.pi*dec/180.))
        self.img_hdr['CRVAL2'] += ddec


########################################################
# prepare images and make catalogues if necessary
all_images = []
all_beams = []
for imagefile in args.images:
    image = ImageSpidx(imagefile)
    all_beams.append(image.get_beam())
    all_images.append(image)
    if args.shift:
        image.make_catalogue()

#####################################################
# find the smallest common beam
if args.beam is None:
    target_beam = findCommonBeam(all_beams)
else:
    target_beam = [args.beam[0]/3600., args.beam[1]/3600., args.beam[2]]

logging.info('Final beam: %.1f" %.1f" (pa %.1f deg)' \
    % (target_beam[0]*3600., target_beam[1]*3600., target_beam[2]))

#####################################################
# find+apply shift w.r.t. first image
if args.shift:
    ref_cat = all_images[0].cat
    for image in all_images[1:]:
        # cross match
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(ref_cat['RA'], ref_cat['DEC']),\
                                             SkyCoord(image.cat['RA'], image.cat['DEC']))
        idx_matched_ref = np.arange(0,len(ref_cat))[sep<target_beam[0]*u.degree]
        idx_matched_img = idx_match[sep<target_beam[0]*u.degree]

        # find & apply shift
        if len(idx_match) == 0:
            logging.warning('%s: No match found, assume no shift.' % image.imagefile)
            continue
            
        dra = np.mean(ref_cat['RA'][idx_matched_ref] - image.cat['RA'][idx_matched_img])
        ddec = np.mean(ref_cat['DEC'][idx_matched_ref] - image.cat['DEC'][idx_matched_img])
        image.apply_shift(dra, ddec)

    # clean up
    #for image in all_images:
    #    os.system(rm ...)

 
######################################################
# Generate regrid headers
rwcs = pywcs(naxis=2)
rwcs.wcs.ctype = all_images[0].get_wcs().wcs.ctype
cdelt = target_beam[1]/5. # 1/5 of minor axes (deg)
logging.info('Pixel scale: %f"' % (cdelt*3600.))
rwcs.wcs.cdelt = [-cdelt, cdelt]
if args.radec is not None:
    mra = radec[0]*np.pi/180
    mdec = radec[1]*np.pi/180
else:
    mra = all_images[0].img_hdr['CRVAL1']
    mdec = all_images[0].img_hdr['CRVAL2']
rwcs.wcs.crval = [mra,mdec]

# if size is not give is taken from the mask
if args.size is None:
    if args.region is not None:
        r = pyregion.open(args.region)
        mask = r.get_mask(header=all_images[0].img_hdr, shape=all_images[0].img_data.shape)
        w = all_images[0].get_wcs()
        y, x = mask.nonzero()
        ra_max, dec_max = w.all_pix2world(np.max(x), np.max(y), 0, ra_dec_order=True)
        ra_min, dec_min = w.all_pix2world(np.min(x), np.min(y), 0, ra_dec_order=True)
        args.size = np.max( [ np.abs(ra_max-ra_min), np.abs(dec_max-dec_min) ] )
    else:
        logging.warning('No size or region provided, use entire size of first image.')
        sys.exit('not implemented')

xsize = int(np.rint(args.size/cdelt))
ysize = int(np.rint(args.size/cdelt))
if xsize % 2 != 0: xsize += 1
if ysize % 2 != 0: ysize += 1
rwcs.wcs.crpix = [xsize/2,ysize/2]

regrid_hdr = rwcs.to_header()
regrid_hdr['NAXIS'] = 2
regrid_hdr['NAXIS1'] = xsize
regrid_hdr['NAXIS2'] = ysize
logging.info('Image size: %f deg (%i %i pixels)' % (args.size,xsize,ysize))

#########################################################
# regrid, convolve and only after apply mask and find noise
for image in all_images:
    image.regrid(regrid_hdr) # regrid initially to speed up convolution
 
    image.convolve(target_beam)
    if args.save:
        intermediate = pyfits.PrimaryHDU(image.img_data, image.img_hdr)
        intermediate.writeto(image.imagefile+'-regrid-conv.fits', overwrite=True)

    if args.region is not None:
        image.apply_region(args.region, invert=True) # after convolution to minimise bad pixels
    if args.noise:
        image.calc_noise() # after mask/convolution
        if args.sigma is not None:
            image.blank_noisy(args.sigma)


#########################################################
# do spdix and write output
frequencies = [ image.freq for image in all_images ]
if args.noise: yerr = [ image.noise for image in all_images ]
else: yerr = None
spidx_data = np.empty(shape=(xsize, ysize))
spidx_data[:] = np.nan
spidx_err_data = np.empty(shape=(xsize, ysize))
spidx_err_data[:] = np.nan

for i in xrange(xsize):
    print '.',
    sys.stdout.flush()
    for j in xrange(ysize):
        val4reg = [ image.img_data[i,j] for image in all_images ]
        if np.isnan(val4reg).any(): continue
        (a, b, sa, sb) = linearfit(x=np.log10(frequencies), y=np.log10(val4reg), yerr=yerr)
        spidx_data[i,j] = a
        spidx_err_data[i,j] = sa

spidx = pyfits.PrimaryHDU(spidx_data, regrid_hdr)
spidx_err = pyfits.PrimaryHDU(spidx_err_data, regrid_hdr)
spidx.writeto(args.output, overwrite=True)
spidx_err.writeto(args.output.replace('.fits','-err.fits'), overwrite=True)
