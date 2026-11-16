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

# Mosaic images

import os, sys, argparse, logging
import numpy as np
import astropy.units as u
from lib_fits import flatten, Image, correct_beam_header
from astropy.io import fits as pyfits
from astropy.wcs import WCS as pywcs
from astropy.table import Table
# https://github.com/astrofrog/reproject
from reproject import reproject_interp
reproj = reproject_interp

# cleanup logs
from astropy.wcs import FITSFixedWarning
import warnings
warnings.filterwarnings('ignore', category=FITSFixedWarning)

ref_catalog = '/homes/fdg/scripts/FIRST_14dec17.fits.gz'

parser = argparse.ArgumentParser(description='Mosaic for LiLF dd-pipeline.')
parser.add_argument('--images', dest='images', nargs='+', help='List of images to combine')
parser.add_argument('--regions', dest='regions', nargs='+', help='List of regions to blank images')
parser.add_argument('--mask', dest='mask', help='One mask with a number per direction, numbers must be in the same order of those given in the "images" parameter.')
parser.add_argument('--beams', dest='beams', nargs='+', help='List of beams')
parser.add_argument('--beamcut', dest='beamcut', type=float, default=0.3, help='Beam level to cut at (default: 0.3, use 0.0 to deactivate)')
parser.add_argument('--beamcorr', dest='beamcorr', action='store_true', help='Pre-correct for beam before combining (default: do not apply)')
parser.add_argument('--beamarm', dest='beamarm', action='store_true', help='Convolve all images to minimum common beam (default: False)')
parser.add_argument('--beamcirc', dest='beamcirc', action='store_true', help='If beamarm is set, then forces the beam to be circular (default: False)')
parser.add_argument('--header', dest='header', help='An image/header to use for the output mosaic coordinates')
parser.add_argument('--noises', dest='noises', type=float, nargs='+', help='UNSCALED Central noise level for weighting: must match numbers of maps')
parser.add_argument('--scales', dest='scales', type=float, nargs='+', help='Scale factors by which maps should be multiplied: must match numbers of maps')
parser.add_argument('--shift', dest='shift', action='store_true', help='Shift images before mosaicing')
parser.add_argument('--find_noise', dest='find_noise', action='store_true', help='Find noise from image (default: assume equal weights, ignored if noises are given)')
parser.add_argument('--use_channel', dest='use_channel', type=int, default=0, help='Channel to be used in a cube image (default: 0)')
parser.add_argument('--use_stokes', dest='use_stokes', type=int, default=0, help='Stokes to be used in a cube image (default: 0)')
parser.add_argument('--save', dest='save', action='store_true', help='Save intermediate results (default: False)')
parser.add_argument('--output', dest='output', default='mosaic.fits', help='Name of output mosaic (default: mosaic.fits)')

args = parser.parse_args()
logging.root.setLevel(logging.DEBUG)

#######################################################
# input check

if args.images is None or len(args.images) < 1:
    logging.error('Requires at least 1 image.')
    sys.exit(1)

if args.scales is not None:
    if len(args.scales) != len(args.images):
        logging.error('Scales provided must match images.')
        sys.exit(1)

if args.noises is not None:
    if len(args.noises) != len(args.images):
        logging.error('Noises provided must match images.')
        sys.exit(1)

if args.beamcirc and not args.beamarm:
    logging.error('--beamcirc requires --beamarm.')
    sys.exit(1)

if args.mask is not None:
    logging.debug('Reading mask: %s.' % args.mask)
    with pyfits.open(args.mask) as hdul:
        mask_n = hdul[0].copy()

if args.shift and not args.beamcorr:
    logging.warning('Attempting shift calculation without beam correction, this is not the best.')

#############################################################

def reproject_fast(data, wcs_in, wcs_out, shape_out):
    """
    Reproject a 2D image, or a stack of images on the same grid, in a single call.
    The coordinate round trip of reproject (a check that takes half of the time) is skipped.
    """
    try:
        out, _ = reproj((data, wcs_in), wcs_out, shape_out=shape_out, roundtrip_coords=False)
    except TypeError: # old reproject without roundtrip_coords
        out, _ = reproj((data, wcs_in), wcs_out, shape_out=shape_out)
    return out


def same_grid(hdr1, hdr2):
    """True if the two headers describe the same pixel grid."""
    for key in ('NAXIS1', 'NAXIS2', 'CTYPE1', 'CTYPE2'):
        if hdr1.get(key) != hdr2.get(key): return False
    for key in ('CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2'):
        if hdr1.get(key) is None or hdr2.get(key) is None or not np.isclose(hdr1[key], hdr2[key], rtol=0, atol=1e-9): return False
    return True


class Direction(Image):

    def __init__(self, imagefile, channel=0, stokes=0):
        logging.debug('Create direction for %s' % imagefile)
        Image.__init__(self, imagefile, channel, stokes)
        # the mosaic never goes back to the original data: do not keep a second copy in memory
        self.img_data_orig = None
        self.img_hdu = None
        self.img_data = np.asarray(self.img_data, dtype=np.float32)
        self.scale = 1.
        self.shift = 0.
        self.beamfile = None
        self.noise = 1.
        self.imagefile = imagefile

    def set_beam_file(self, beamfile):
        if not os.path.exists(beamfile):
            logging.error('Beam file %s not found.' % beamfile)
            sys.exit(1)
        self.beamfile = beamfile
        self.beam_hdr, beam_data = flatten(self.beamfile)
        logging.debug('%s: set beam file %s' % (self.imagefile, beamfile))
        if beam_data.shape != self.img_data.shape:
            beamfile = self.imagefile+'__beam.fits'
            if os.path.exists(beamfile):
                logging.debug('Beam and image shape are different, load regridded beam %s...' % beamfile)
                self.beam_hdr, beam_data = flatten(beamfile)
            else:
                logging.warning('Beam and image shape are different, regrid beam and save to %s...' % beamfile)
                beam_data = reproject_fast(beam_data, pywcs(self.beam_hdr), pywcs(self.img_hdr), self.img_data.shape)
                # other mosaics running at the same time use the same file: write it
                # under a temporary name so that nobody can read it half-written
                beamfile_tmp = '%s.tmp%i' % (beamfile, os.getpid())
                pyfits.writeto(beamfile_tmp, header=self.img_hdr, data=beam_data.astype(np.float32), overwrite=True)
                os.replace(beamfile_tmp, beamfile)
                self.beam_hdr = self.img_hdr
            self.beamfile = beamfile
        self.beam_data = np.asarray(beam_data, dtype=np.float32)

    def apply_beam_cut(self, beamcut=0.3):
        if self.beamfile is None: return
        # blank with nan (not 0): the interpolation of the regridding then drops the pixels
        # on the edge of the cut instead of mixing them with zeros. "not >=" blanks also nan beams
        self.img_data[~(self.beam_data >= beamcut)] = np.nan

    def apply_beam_corr(self):
        if self.beamfile is None: return
        with np.errstate(divide='ignore', invalid='ignore'):
            self.img_data /= self.beam_data

    def get_valid(self):
        """Pixels with data: blanked pixels are nan, 0 marks regions without data"""
        return np.isfinite(self.img_data) & (self.img_data != 0)

    def calc_weight(self):
        self.weight_data = np.ones_like(self.img_data)
        if self.beamfile is not None:
            self.weight_data *= self.beam_data
        # at this point this is the beam factor: we want 1/sigma**2.0, so divide by central noise and square
        self.weight_data /= self.noise * self.scale
        # https://en.wikipedia.org/wiki/Inverse-variance_weighting
        self.weight_data = self.weight_data**2.0
        # no data, no weight: nan in both so that image and weight are dropped together when regridding
        invalid = ~self.get_valid() | ~np.isfinite(self.weight_data)
        self.weight_data[invalid] = np.nan
        self.img_data[invalid] = np.nan
        # the beam is not needed anymore
        if self.beamfile is not None:
            self.beam_data = None

    def reproject(self, regrid_hdr):
        """
        Regrid image and weights on the output grid, only in the part of it covered by valid data.
        Return image, weights (0 where there are no data, never nan) and the slices
        of the output grid they refer to; None if the image does not contribute.
        """
        valid = np.isfinite(self.weight_data) & (self.weight_data > 0)
        if not valid.any():
            return None
        ysize, xsize = regrid_hdr['NAXIS2'], regrid_hdr['NAXIS1']
        wcs_in = self.get_wcs().celestial
        wcs_out = pywcs(regrid_hdr).celestial

        # box of the input with valid data
        ys = np.where(valid.any(axis=1))[0]
        xs = np.where(valid.any(axis=0))[0]
        y0, y1, x0, x1 = ys[0], ys[-1]+1, xs[0], xs[-1]+1
        del valid

        # corresponding box of the output: follow the edge of the input box
        n = 200
        ex = np.concatenate([np.linspace(x0-1, x1, n), np.full(n, x1), np.linspace(x0-1, x1, n), np.full(n, x0-1)])
        ey = np.concatenate([np.full(n, y0-1), np.linspace(y0-1, y1, n), np.full(n, y1), np.linspace(y0-1, y1, n)])
        ra, dec = wcs_in.wcs_pix2world(ex, ey, 0)
        ox, oy = wcs_out.wcs_world2pix(ra, dec, 0)
        if np.all(np.isfinite(ox)) and np.all(np.isfinite(oy)):
            ox0, ox1 = max(0, int(np.floor(ox.min()))-2), min(xsize, int(np.ceil(ox.max()))+3)
            oy0, oy1 = max(0, int(np.floor(oy.min()))-2), min(ysize, int(np.ceil(oy.max()))+3)
        else: # part of the image cannot be projected on the output: do not try to be smart
            ox0, ox1, oy0, oy1 = 0, xsize, 0, ysize
        if ox1 <= ox0 or oy1 <= oy0:
            return None

        # image and weights share the coordinate computation
        stack = np.stack([self.img_data[y0:y1, x0:x1], self.weight_data[y0:y1, x0:x1]])
        try:
            r, w = reproject_fast(stack, wcs_in[y0:y1, x0:x1], wcs_out[oy0:oy1, ox0:ox1], (oy1-oy0, ox1-ox0))
        except ValueError: # old reproject that does not broadcast over the extra axis
            r = reproject_fast(stack[0], wcs_in[y0:y1, x0:x1], wcs_out[oy0:oy1, ox0:ox1], (oy1-oy0, ox1-ox0))
            w = reproject_fast(stack[1], wcs_in[y0:y1, x0:x1], wcs_out[oy0:oy1, ox0:ox1], (oy1-oy0, ox1-ox0))
        del stack

        # a pixel counts only if both image and weight are defined: a nan in the image
        # must never enter the mosaic as a 0 with a weight
        bad = ~np.isfinite(r) | ~np.isfinite(w) | (w <= 0)
        r[bad] = 0
        w[bad] = 0
        return r, w, (slice(oy0, oy1), slice(ox0, ox1))

    def calc_shift(self, ref_cat):
        """
        Find a shift cross-matching source extracted from the image and a given catalog.
        """
        import bdsf
        from astropy.coordinates import match_coordinates_sky
        from astropy.coordinates import SkyCoord
        from astropy.stats import median_absolute_deviation

        # if there are no data in the image prevent the crash of bdsf
        if np.isnan(self.img_data).all(): return
        img_cat = self.imagefile+'.cat'
        if not os.path.exists(img_cat):
            bdsf_img = bdsf.process_image(self.imagefile, rms_box=(100,30), \
                thresh_pix=5, thresh_isl=3, atrous_do=False, \
                adaptive_rms_box=True, adaptive_thresh=100, rms_box_bright=(30,10))
            bdsf_img.write_catalog(outfile=img_cat, catalog_type='srl', format='fits', clobber=True)

        # read catlogue
        ref_t = Table.read(ref_cat)
        img_t = Table.read(img_cat)
        logging.debug('SHIFT: Initial len: %i (ref:%i)' % (len(img_t),len(ref_t)))
 
        # reduce to isolated sources LOFAR
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(img_t['RA'], img_t['DEC']),\
                                                  SkyCoord(img_t['RA'], img_t['DEC']), nthneighbor=2)
    
        idx_match_img = np.arange(0,len(img_t))[sep>3*self.get_beam()[0]*u.arcsec]
        img_t = img_t[idx_match_img]
        # reduce to isolated sources REF
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(ref_t['RA'], ref_t['DEC']),\
                                                  SkyCoord(ref_t['RA'], ref_t['DEC']), nthneighbor=2)
    
        idx_match_ref = np.arange(0,len(ref_t))[sep>self.get_beam()[0]*u.arcsec]
        ref_t = ref_t[idx_match_ref]
        logging.debug('SHIFT: After isaolated sources len: %i (ref:%i)' % (len(img_t),len(ref_t)))
    
        # reduce to compact sources
        img_t = img_t[ (img_t['S_Code'] == 'S') ]
        img_t = img_t[ (img_t['Total_flux']/img_t['Peak_flux']) < 2 ]
        ref_t = ref_t[ (ref_t['FINT']/ref_t['FPEAK']) < 1.2 ]
        logging.debug('SHIFT: After compact source len: %i (ref:%i)' % (len(img_t),len(ref_t)))
    
        # cross match
        idx_match, sep, _ = match_coordinates_sky(SkyCoord(img_t['RA'], img_t['DEC']),\
                                                  SkyCoord(ref_t['RA'], ref_t['DEC']))

        if len(idx_match) == 0:
            logging.warning('No match found in the reference catalogue.')
            return

        sep_mad = median_absolute_deviation(sep[np.where(sep < (3*self.get_beam()[0])*u.deg)])
        sep_med = np.median(sep[np.where(sep < (3*self.get_beam()[0])*u.deg)])
        logging.debug('SHIFT: Sep init Med: %f" - MAD: %f"' % (sep_med.arcsec, sep_mad.arcsec))
        sep_mad_old = 0. * u.deg
        i = 0
        while sep_mad != sep_mad_old and not i > 100:
            sep_mad_old = sep_mad
            idx = np.where(sep-sep_med < 7 * sep_mad)
            sep_mad = median_absolute_deviation(sep[idx])
            sep_med = np.median(sep[idx])
            logging.debug('SHIFT: Sep Med: %.2f" - MAD: %.2f" (n sources:%i)' % \
                    (sep_med.arcsec, sep_mad.arcsec, len(sep[idx])))
            if np.isnan(sep_mad):
                logging.error('MAD diverged')
                sys.exit(1)
            i+=1
    
        idx_match_ref = idx_match[sep-sep_med < 3*sep_mad]
        idx_match_img = np.arange(0,len(img_t))[sep-sep_med < 3*sep_mad]
        img_t = img_t[idx_match_img]
    
        logging.debug('SHIFT: After match source len: %i' % len(img_t))

        # find & apply shift
        ddec = ref_t['DEC'][idx_match_ref] - img_t['DEC']
        dra = ref_t['RA'][idx_match_ref] - img_t['RA']
        dra[ dra>180 ] -= 360
        dra[ dra<-180 ] += 360

        self.apply_shift(np.mean(dra), np.mean(ddec)) # ra is in deg on the sphere, no dec correction

        # debug
        #self.write(self.imagefile.replace('fits','shift.fits'))
        #img_t['RA'] += np.mean(dra)
        #img_t['DEC'] += np.mean(ddec)
        #img_t.write(img_cat.replace('cat','shiftedcat'), format='fits', overwrite=True)
        #ref_t[idx_match_ref].write(img_cat.replace('cat','refcat'), format='fits', overwrite=True)

        # clean up
        #if not args.save:
        #    os.system('rm '+img_cat)

def get_image_beam(imagefile):
    """Beam of an image from its header, without loading the data"""
    with pyfits.open(imagefile) as hdul:
        header = correct_beam_header(hdul[0].header)
        return [header['BMAJ'], header['BMIN'], header['BPA']]


if args.beamarm:

    logging.info('Reading beams...')
    beams = [get_image_beam(image) for image in args.images]
    if args.beamcirc:
        maxmaj = np.max([b[0] for b in beams])
        common_beam = [maxmaj*1.01, maxmaj*1.01, 0.] # add 1% to prevent crash in convolution
    elif beams.count(beams[0]) == len(beams):
        # all beams are already exactly the same
        common_beam = list(beams[0])
    else:
        from radio_beam import Beams
        my_beams = Beams([b[0] for b in beams] * u.deg, [b[1] for b in beams] * u.deg, [b[2] for b in beams] * u.deg)
        common_beam = my_beams.common_beam()
        common_beam = [common_beam.major.value, common_beam.minor.value, common_beam.pa.value]

    logging.debug('Minimum common beam: %.1f" %.1f" (pa %.1f deg)' % \
             (common_beam[0]*3600., common_beam[1]*3600., common_beam[2]))


def prepare_direction(i):
    """Load image i and get it ready to be regridded: beam, cuts, noise, weights"""
    d = Direction(args.images[i], channel=args.use_channel, stokes=args.use_stokes)

    if args.beamarm:
        d.convolve(common_beam)

    if args.beams is not None:
        d.set_beam_file(args.beams[i])
        d.apply_beam_cut(beamcut = args.beamcut)

    if args.regions is not None:
        d.apply_region(args.regions[i], blankvalue=np.nan, invert=True)

    if args.noises is not None: d.noise = args.noises[i]
    elif args.find_noise:
        d.calc_noise(force_recalc=True) # after beam cut/mask
        if not d.noise: # None or 0: no data left in the image
            logging.warning('%s: no valid data to estimate the noise, image ignored.' % d.imagefile)
            d.noise = np.nan # all weights become nan

    if args.scales is not None: d.scale = args.scales[i]

    if args.beamcorr: d.apply_beam_corr() # after noise calculation

    d.calc_weight() # after setting: beam, noise, scale

    if args.shift:
        d.calc_shift(ref_catalog)

    return d


logging.info("Working on %i images..." % len(args.images))

# prepare header for final gridding
directions = None
if args.header is None:
    logging.warning('Calculate output headers...')
    # the extent of all images is needed before starting: keep them in memory
    directions = [prepare_direction(i) for i in range(len(args.images))]
    mra = np.mean( np.array([d.get_wcs().wcs.crval[0]%360 for d in directions]) )
    mdec = np.mean( np.array([d.get_wcs().wcs.crval[1] for d in directions]) )

    logging.info('Will make mosaic at %f %f' % (mra,mdec))
    
    # we make a reference WCS and use it to find the extent in pixels
    # needed for the combined image

    rwcs = pywcs(naxis=2)
    rwcs.wcs.ctype = directions[0].get_wcs().wcs.ctype
    rwcs.wcs.cdelt = directions[0].get_wcs().wcs.cdelt
    rwcs.wcs.crval = [mra,mdec]
    rwcs.wcs.crpix = [1,1]

    xmin=0
    xmax=0
    ymin=0
    ymax=0
    for d in directions:
        w = d.get_wcs()
        ys, xs = np.where(d.get_valid())
        if xs.size == 0:
            continue
        axmin = xs.min()
        aymin = ys.min()
        axmax = xs.max()
        aymax = ys.max()
        del(xs)
        del(ys)
        for x,y in ((axmin,aymin),(axmax,aymin),(axmin,aymax),(axmax,aymax)):
            ra, dec = [float(f) for f in w.wcs_pix2world(x,y,0)]
            #print ra,dec
            nx, ny = [float (f) for f in rwcs.wcs_world2pix(ra,dec,0)]
            #print nx,ny
            if nx < xmin: xmin=nx
            if nx > xmax: xmax=nx
            if ny < ymin: ymin=ny
            if ny > ymax: ymax=ny

    #print 'co-ord range:', xmin, xmax, ymin, ymax

    xsize = int(xmax-xmin)
    ysize = int(ymax-ymin)

    rwcs.wcs.crpix = [-int(xmin)+1,-int(ymin)+1]
    #print 'checking:', rwcs.wcs_world2pix(mra,mdec,0)

    regrid_hdr = rwcs.to_header()
    regrid_hdr['NAXIS'] = 2
    regrid_hdr['NAXIS1'] = xsize
    regrid_hdr['NAXIS2'] = ysize

else:
    try:
        logging.info("Using %s header for final gridding." % args.header)
        with pyfits.open(args.header) as hdul:
            regrid_hdr = hdul[0].header.copy()
        xsize = regrid_hdr['NAXIS1']
        ysize = regrid_hdr['NAXIS2']
    except Exception:
        logging.error("--header must be a fits file.")
        sys.exit(1)

logging.info('Making mosaic...')
isum = np.zeros([ysize,xsize])
wsum = np.zeros_like(isum)
if args.mask is not None:
    logging.debug('Reprojecting mask...')
    outname = args.mask.replace('.fits','-reproj.fits')
    mask_data = None
    if os.path.exists(outname):
        with pyfits.open(outname) as hdul:
            if same_grid(hdul[0].header, regrid_hdr):
                logging.debug('Loading %s...' % outname)
                mask_data = hdul[0].data.copy()
            else:
                logging.warning('%s is on a different grid, ignored.' % outname)
    if mask_data is None:
        mask_data, footprint = reproj((mask_n.data, mask_n.header), regrid_hdr, order='bilinear')#, parallel=True)
        if args.save:
            pyfits.writeto(outname, header=regrid_hdr, data=mask_data, overwrite=True)
    del mask_n

    # get numbers into mask in increasing order
    mask_numbers = sorted(np.unique(mask_data))

def load_saved(outname):
    """Load a regridded image saved by a previous run with --save, if it is on the output grid"""
    if not os.path.exists(outname):
        return None
    with pyfits.open(outname) as hdul:
        if not same_grid(hdul[0].header, regrid_hdr):
            logging.warning('%s is on a different grid, ignored.' % outname)
            return None
        logging.debug('Loading %s...' % outname)
        return np.nan_to_num(hdul[0].data)

for i in range(len(args.images)):
    logging.info('Working on: %s' % args.images[i])

    # regridded image and weights saved by a previous run (full output grid)
    outname_r = args.images[i].replace('.fits','-reproj.fits')
    outname_w = args.images[i].replace('.fits','-reprojW.fits')
    r, w = load_saved(outname_r), load_saved(outname_w)
    if r is not None and w is not None:
        isum += r*w
        wsum += w
        del r, w
        continue

    # one image at a time in memory, unless they are all already there
    d = directions[i] if directions is not None else prepare_direction(i)
    logging.debug('Reprojecting data and weights...')
    result = d.reproject(regrid_hdr)
    if directions is not None: directions[i] = None
    del d
    if result is None:
        logging.warning('%s: no valid data on the output grid, image ignored.' % args.images[i])
        continue
    r, w, box = result

    if args.mask is not None:
        w[mask_data[box] != mask_numbers[i]] = 0
    if args.save:
        for outname, data in ((outname_r, r), (outname_w, w)):
            full = np.zeros([ysize,xsize], dtype=np.float32)
            full[box] = data
            pyfits.writeto(outname, header=regrid_hdr, data=full, overwrite=True)
        del full
    logging.debug('Add to mosaic...')
    isum[box] += r*w
    wsum[box] += w
    del r, w, result

logging.debug('Write mosaic: %s...' % args.output)
# no weight means that no image has valid data there
valid = wsum != 0
isum[valid] /= wsum[valid]
isum[~valid] = np.nan
del wsum, valid

#set beam
try:
    regrid_hdr['BMAJ'] = common_beam[0]
    regrid_hdr['BMIN'] = common_beam[1]
    regrid_hdr['BPA'] = common_beam[2]
except (NameError, KeyError):
    logging.warning('Setting beam in the headers equal to: %s' % args.images[0])
    with pyfits.open(args.images[0]) as hdul:
        for ch in ('BMAJ', 'BMIN', 'BPA'):
            regrid_hdr[ch] = hdul[0].header[ch]

regrid_hdr['ORIGIN'] = 'LiLF-pipeline-mosaic'
regrid_hdr['UNITS'] = 'Jy/beam'

# single precision is plenty for an image and halves the file size
pyfits.writeto(args.output, header=regrid_hdr, data=isum.astype(np.float32), overwrite=True)

logging.debug('Done.')
