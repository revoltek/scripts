#!/usr/bin/env python3
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
import numpy as np
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
from lib_linearfit import linear_fit, linear_fit_bootstrap
from lib_fits import flatten
# https://github.com/astrofrog/reproject
from reproject import reproject_interp, reproject_exact
reproj = reproject_exact
logging.root.setLevel(logging.DEBUG)

parser = argparse.ArgumentParser(description='Make spectral index maps, e.g. spidxmap.py --region ds9.reg --noise --sigma 5 --save *fits')
parser.add_argument('images', nargs='+', help='List of images to use for spidx')
parser.add_argument('--beam', dest='beam', nargs='+', type=float, help='3 parameters final beam to convolve all images (BMAJ (arcsec), BMIN (arcsec), BPA (deg))')
parser.add_argument('--region', dest='region', help='Ds9 region to restrict analysis')
parser.add_argument('--noiseregion', dest='noiseregion', help='Ds9 region to calculate rms noise')
parser.add_argument('--noisesigma', dest='noisesigma', default=5, help='Sigma used in the calc_noise function (default: 5)')
parser.add_argument('--size', dest='size', type=float, help='Size (horizontal and vertical) of final image in degree')
parser.add_argument('--radec', dest='radec', nargs='+', type=float, help='RA/DEC where to center final image in deg (if not given, center on first image)')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before calculating spidx')
parser.add_argument('--noise', dest='noise', action='store_true', help='Calculate noise of each image')
parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results')
parser.add_argument('--sigma', dest='sigma', type=float, help='Restrict to pixels above this sigma in all images')
parser.add_argument('--circbeam', dest='circbeam', action='store_true', help='Force final beam to be circular (default: False, use minimum common beam area)')
parser.add_argument('--bootstrap', dest='bootstrap', action='store_true', help='Use bootstrap to estimate errors (default: use normal X|Y with errors)')
parser.add_argument('--output', dest='output', default='spidx.fits', help='Name of output mosaic (default: spidx.fits)')

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

from lib_fits import Image, find_freq
class ImageSpidx(Image):

    def __init__(self, imagefile):
        Image.__init__(self, imagefile)

    def regrid(self, regrid_hdr):
        logging.debug('%s: regridding' % (self.imagefile))
        self.img_data, __footprint = reproj((self.img_data, self.img_hdr), regrid_hdr, parallel=True)
        beam = self.get_beam()
        freq = find_freq(self.img_hdr)
        self.img_hdr = regrid_hdr
        self.img_hdr['FREQ'] = freq
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
    
    if all_beams.count(all_beams[0]) == len(all_beams):
        # all beams are already exactly the same
        target_beam = all_beams[0]

    elif args.circbeam:
        maxmaj = np.max([b[0] for b in all_beams])
        target_beam = [maxmaj*1.01, maxmaj*1.01, 0.] # add 1% to prevent crash in convolution
    else:
        from radio_beam import Beams
        my_beams = Beams([b[0] for b in all_beams] * u.deg, [b[1] for b in all_beams] * u.deg, [b[2] for b in all_beams] * u.deg)
        common_beam = my_beams.common_beam()
        target_beam = [common_beam.major.value, common_beam.minor.value, common_beam.pa.value]
else:
    target_beam = [args.beam[0]/3600., args.beam[1]/3600., args.beam[2]]

logging.info('Final beam: %.1f" %.1f" (pa %.1f deg)' \
    % (target_beam[0]*3600., target_beam[1]*3600., target_beam[2]))

#####################################################
# find+apply shift w.r.t. first image
if args.shift:
    ref_cat = all_images[0].cat
    # keep only point sources
    print('Reference catalogue:', ref_cat)
    for image in all_images[1:]:
        # cross match
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(ref_cat['RA'], ref_cat['DEC']),\
                                             SkyCoord(image.cat['RA'], image.cat['DEC']))
        idx_matched_ref = np.arange(0,len(ref_cat))[sep<target_beam[0]*u.degree]
        idx_matched_img = idx_match[sep<target_beam[0]*u.degree]

        # find & apply shift
        if len(idx_match) < 3:
            logging.warning('%s: Not enough matches found, assume no shift.' % image.imagefile)
            continue
            
        dra = ref_cat['RA'][idx_matched_ref] - image.cat['RA'][idx_matched_img]
        dra[ dra>180 ] -= 360
        dra[ dra<-180 ] += 360
        ddec = ref_cat['DEC'][idx_matched_ref] - image.cat['DEC'][idx_matched_img]
        flux = ref_cat['Peak_flux'][idx_matched_ref]
        image.apply_shift(np.average(dra, weights=flux), np.average(ddec, weights=flux))

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
        intermediate = pyfits.PrimaryHDU(mask.astype(float), all_images[0].img_hdr)
        intermediate.writeto('mask.fits', overwrite=True)
        w = all_images[0].get_wcs()
        y, x = mask.nonzero()
        ra_max, dec_max = w.all_pix2world(np.max(x), np.max(y), 0, ra_dec_order=True)
        ra_min, dec_min = w.all_pix2world(np.min(x), np.min(y), 0, ra_dec_order=True)
        args.size = 1.2*2.*np.max( [ np.max([np.abs(ra_max-mra),np.abs(ra_min-mra)]), np.max([np.abs(dec_max-mdec),np.abs(dec_min-mdec)]) ] )
        #print(ra_min,ra_max,dec_min,dec_max)
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
regrid_hdr['EQUINOX'] = 2000.0
regrid_hdr['RADESYSa'] = 'J2000 '
logging.info('Image size: %f deg (%i %i pixels)' % (args.size,xsize,ysize))

#########################################################
# regrid, convolve and only after apply mask and find noise
for image in all_images:

    if os.path.exists(image.imagefile+'-conv.fits'):
        logging.warning('Load: '+image.imagefile+'-conv.fits')
        data, hdr = pyfits.getdata(image.imagefile+'-conv.fits', 0, header=True)
        image.img_data = data
        image.img_hdr = hdr
        image.set_beam([hdr['BMAJ'], hdr['BMIN'], hdr['BPA']])
    else:
        image.convolve(target_beam)
        if args.save:
            intermediate = pyfits.PrimaryHDU(image.img_data, image.img_hdr)
            intermediate.writeto(image.imagefile+'-conv.fits', overwrite=True)

    if os.path.exists(image.imagefile+'-regrid-conv.fits'):
        logging.warning('Load: '+image.imagefile+'-regrid-conv.fits')
        data, hdr = pyfits.getdata(image.imagefile+'-regrid-conv.fits', 0, header=True)
        image.img_data = data
        image.img_hdr = hdr
        image.set_beam([hdr['BMAJ'], hdr['BMIN'], hdr['BPA']])
    else:
        image.regrid(regrid_hdr)
        if args.save:
            intermediate = pyfits.PrimaryHDU(image.img_data, image.img_hdr)
            intermediate.writeto(image.imagefile+'-regrid-conv.fits', overwrite=True)
 
    if args.noise:
        if args.sigma is not None:
            # usually the sigma used for the blanking is rather low, better to increase it for the calc_noise
            image.calc_noise(sigma=args.noisesigma, bg_reg=args.noiseregion) # after mask?/convolution
            image.blank_noisy(args.sigma)
        else:
            image.calc_noise() # after mask?/convolution

    if args.region is not None:
        image.apply_region(args.region, invert=True) # after convolution to minimise bad pixels


#########################################################
# do spdix and write output
frequencies = [ image.freq for image in all_images ]
if args.noise: yerr = [ image.noise for image in all_images ]
else: yerr = None
spidx_data = np.empty(shape=(xsize, ysize))
spidx_data[:] = np.nan
spidx_err_data = np.empty(shape=(xsize, ysize))
spidx_err_data[:] = np.nan

for i in range(xsize):
    print('%i/%i' % (i,xsize), end='\r')
    sys.stdout.flush()
    for j in range(ysize):
        val4reg = [ image.img_data[i,j] for image in all_images ]
        if np.isnan(val4reg).any() or (np.array(val4reg) < 0).any(): continue
        if args.bootstrap:
            (a, b, sa, sb) = linear_fit_bootstrap(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
        else:
            (a, b, sa, sb) = linear_fit(x=frequencies, y=val4reg, yerr=yerr, tolog=True)
        spidx_data[i,j] = a
        spidx_err_data[i,j] = sa

spidx = pyfits.PrimaryHDU(spidx_data, regrid_hdr)
spidx_err = pyfits.PrimaryHDU(spidx_err_data, regrid_hdr)
logging.info('Save %s (and errors)' % args.output)
if args.output[-5:] == '.fits':
    spidx.writeto(args.output, overwrite=True)
    spidx_err.writeto(args.output.replace('.fits','-err.fits'), overwrite=True)
else:
    spidx.writeto(args.output+'.fits', overwrite=True)
    spidx_err.writeto(args.output+'-err.fits', overwrite=True)
